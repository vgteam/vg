/**
 * \file minimizer_mapper_from_chains.cpp
 * Defines the code for the long-read code path for the
 * minimizer-and-GBWT-based mapper (long read Giraffe).
 */

#include "minimizer_mapper.hpp"

#include "annotation.hpp"
#include "crash.hpp"
#include "path_subgraph.hpp"
#include "multipath_alignment.hpp"
#include "split_strand_graph.hpp"
#include "subgraph.hpp"
#include "statistics.hpp"
#include "algorithms/count_covered.hpp"
#include "algorithms/intersect_path_offsets.hpp"
#include "algorithms/extract_containing_graph.hpp"
#include "algorithms/extract_connecting_graph.hpp"
#include "algorithms/extract_extending_graph.hpp"
#include "algorithms/chain_items.hpp"

#include <bdsg/overlays/strand_split_overlay.hpp>
#include <gbwtgraph/algorithms.h>
#include <gbwtgraph/cached_gbwtgraph.h>

#include <iostream>
#include <algorithm>
#include <cmath>
#include <cfloat>

// Turn on debugging prints
//#define debug
// Turn on printing of minimizer fact tables
//#define print_minimizer_table
// Dump local graphs that we align against 
//#define debug_dump_graph
// Dump fragment length distribution information
//#define debug_fragment_distr
//Do a brute force check that clusters are correct
//#define debug_validate_clusters

namespace vg {

using namespace std;

/// Class for an error representing that chaining has backed us into some kind
/// of corner and we can't actually produce an alignment. We can throw this to
/// leave the read unmapped, complain, and try the next read.
class ChainAlignmentFailedError : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

static void set_coverage_flags(std::vector<bool>& flags, size_t start, size_t end) {
    for (size_t i = start; i < end; i++) {
        flags[i] = true;
    }
}

static double get_fraction_covered(const std::vector<bool>& flags) {
    size_t covered_bases = 0;
    for (bool flag : flags) {
        if (flag) {
            covered_bases++;
        }
    }
    return (double) covered_bases / flags.size();
}

/// Get the forward-relative-to-the-read version of a seed's position. Will
/// have the correct orientation, but won't necessarily be to any particular
/// (i.e. first or last) base of the seed.
static pos_t forward_pos(const MinimizerMapper::Seed& seed, const VectorView<MinimizerMapper::Minimizer>& minimizers, const HandleGraph& graph) {
    pos_t position = seed.pos;
    if (minimizers[seed.source].value.is_reverse) {
        // Need to flip the position, for which we need to fetch the node length.
        position = reverse_base_pos(position, graph.get_length(graph.get_handle(id(position), is_rev(position))));
    }
    return position;
}

void MinimizerMapper::dump_debug_dotplot(const std::string& name, const std::string& marker, const VectorView<Minimizer>& minimizers, const std::vector<Seed>& seeds, const std::vector<size_t>& included_seeds, const std::vector<size_t>& highlighted_seeds, const PathPositionHandleGraph* path_graph) {
    if (!path_graph) {
        // We don't have a path positional graph for this
        return;
    }

    // Log the best bucket's seed positions in read and linear reference
    TSVExplainer exp(true, name + "-dotplot");

    // We need to know which seeds to highlight
    std::unordered_set<size_t> highlight_set;
    for (auto& seed_num : highlighted_seeds) {
        highlight_set.insert(seed_num);
    }

    for (auto& seed_num : included_seeds) {
        // For each seed in the best bucket
        auto& seed = seeds.at(seed_num);
        
        // Get its effective path positions again
        auto offsets = algorithms::nearest_offsets_in_paths(path_graph, seed.pos, 100);

        for (auto& handle_and_positions : offsets) {
            std::string path_name = path_graph->get_path_name(handle_and_positions.first);
            for (auto& position : handle_and_positions.second) {
                // For each position on a ref path that this seed is at, log a line
                exp.line();
                if (highlight_set.count(seed_num)) {
                    // Contig and a marker
                    exp.field(path_name + "-" + marker);
                } else {
                    // Contig
                    exp.field(path_name);
                }
                // Offset on contig
                exp.field(position.first);
                // Offset in read
                exp.field(minimizers[seed.source].forward_offset());
            }
        }

    }
}

std::vector<MinimizerMapper::Seed> MinimizerMapper::reseed_between(
    size_t read_region_start,
    size_t read_region_end,
    pos_t left_graph_pos,
    pos_t right_graph_pos,
    const HandleGraph& graph,
    const VectorView<Minimizer>& minimizers,
    const std::function<void(const Minimizer&, const std::vector<nid_t>&, const std::function<void(const pos_t&)>&)>& for_each_pos_for_source_in_subgraph
) const {
    
    // We are going to make up some seeds
    std::vector<MinimizerMapper::Seed> forged_items;                                    
    
    
    std::vector<pos_t> seed_positions;
    seed_positions.reserve(2);
    std::vector<size_t> position_forward_max_dist;
    position_forward_max_dist.reserve(2);
    std::vector<size_t> position_backward_max_dist;
    position_backward_max_dist.reserve(2);
    
    if (!is_empty(left_graph_pos)) {
        // We have a left endpoint
        seed_positions.emplace_back(left_graph_pos);
        position_forward_max_dist.emplace_back(this->reseed_search_distance);
        position_backward_max_dist.emplace_back(0);
    }
    
    if (!is_empty(right_graph_pos)) {
        // We have a left endpoint
        seed_positions.emplace_back(right_graph_pos);
        position_forward_max_dist.emplace_back(0);
        position_backward_max_dist.emplace_back(this->reseed_search_distance);
    }
    
    std::vector<nid_t> sorted_ids;
    {
        bdsg::HashGraph subgraph;
        // TODO: can we use connecting graph again?
        // TODO: Should we be using more seeds from the cluster?
        algorithms::extract_containing_graph(&graph, &subgraph, seed_positions, this->reseed_search_distance);
        sorted_ids.reserve(subgraph.get_node_count());
        subgraph.for_each_handle([&](const handle_t& h) {
            sorted_ids.push_back(subgraph.get_id(h));
        });
    }
    std::sort(sorted_ids.begin(), sorted_ids.end());
    
    if (this->show_work) {
        #pragma omp critical (cerr)
        {
            std::cerr << log_name() << "Reseeding against nodes ";
            // Dump the nodes as consecutive ranges 
            nid_t prev_node;
            nid_t printed_node;
            for (size_t i = 0; i < sorted_ids.size(); i++) {
                if (i == 0 || prev_node + 1 != sorted_ids[i]) {
                    if (i > 0) {
                        std::cerr << "-" << prev_node << ", ";
                    }
                    std::cerr << sorted_ids[i];
                    printed_node = sorted_ids[i];
                }
                prev_node = sorted_ids[i];
            }
            if (!sorted_ids.empty() && printed_node != sorted_ids.back()) {
                std::cerr << "-" << sorted_ids.back();
            }
            std::cerr << endl;
        }
    }
    
    for (size_t i = 0; i < minimizers.size(); i++) {
        auto& m = minimizers[i];
        
        if (m.forward_offset() < read_region_start || m.forward_offset() + m.length > read_region_end) {
            // Minimizer is not in the range we care about.
            // TODO: Find a faster way to find the relevant minimizers that doesn't require a scan! Sort them by start position or something.
            continue;
        }
        
        if (this->show_work) {
            #pragma omp critical (cerr)
            {
                std::cerr << log_name() << "Query minimizer #" << i << " at " << m.forward_offset() << " which overall has " << m.hits << " hits" << std::endl;
            }
        }
        
        // We may see duplicates, so we want to do our own deduplication.
        unordered_set<pos_t> seen;
        
        size_t hit_count = 0;
        
        // Find all its hits in the part of the graph between the bounds
        for_each_pos_for_source_in_subgraph(m, sorted_ids, [&](const pos_t& pos) {
            // So now we know pos corresponds to read base
            // m.value.offset, in the read's forward orientation.
            
            // Forge an item.
            forged_items.emplace_back();
            forged_items.back().pos = pos;
            forged_items.back().source = i;
            
            // Record the hit
            hit_count++;
        });
        
        if (this->show_work) {
            #pragma omp critical (cerr)
            {
                std::cerr << log_name() << "\tFound " << hit_count << "/" << m.hits << " hits" << std::endl;
            }
        }
    }
    
    // TODO: sort and deduplicate the new seeds
    
    return forged_items;
                                        
}

MinimizerMapper::chain_set_t MinimizerMapper::chain_clusters(const Alignment& aln, const VectorView<Minimizer>& minimizers, const std::vector<Seed>& seeds, const ZipCodeTree& zip_code_tree, const std::vector<Cluster>& clusters, const chain_config_t& cfg, size_t old_seed_count, size_t new_seed_start, Funnel& funnel, size_t seed_stage_offset, size_t reseed_stage_offset, LazyRNG& rng) const {

    // Convert the seeds into chainable anchors in the same order
    vector<algorithms::Anchor> seed_anchors = this->to_anchors(aln, minimizers, seeds);
    
    // We need to remember which order we did the chains in, independent of the provenance funnel.
    // TODO: Drop this when we are done with fragment statistics!
    vector<size_t> cluster_nums;
    cluster_nums.reserve(clusters.size());
    
    // These are the collections of chains for all the clusters, as score and sequence of visited seeds.
    vector<vector<pair<int, vector<size_t>>>> cluster_chains;
    cluster_chains.reserve(clusters.size());
    
    // To compute the windows for explored minimizers, we need to get
    // all the minimizers that are explored.
    SmallBitset minimizer_explored(minimizers.size());
    //How many hits of each minimizer ended up in each cluster we kept?
    vector<vector<size_t>> minimizer_kept_cluster_count; 

    size_t kept_cluster_count = 0;
    
    // What cluster seeds define the space for clusters' chosen chains?
    vector<vector<size_t>> cluster_chain_seeds;
    cluster_chain_seeds.reserve(clusters.size());
    
    //Process clusters sorted by both score and read coverage
    process_until_threshold_c<double>(clusters.size(), [&](size_t i) -> double {
            return clusters[i].coverage;
        }, [&](size_t a, size_t b) -> bool {
            return ((clusters[a].coverage > clusters[b].coverage) ||
                    (clusters[a].coverage == clusters[b].coverage && clusters[a].score > clusters[b].score));
        }, cfg.cluster_coverage_threshold, cfg.min_clusters_to_chain, cfg.max_clusters_to_chain, rng, [&](size_t cluster_num) -> bool {
            // Handle sufficiently good clusters in descending coverage order
            
            const Cluster& cluster = clusters[cluster_num];
            if (track_provenance) {
                funnel.pass("cluster-coverage", cluster_num, cluster.coverage);
                funnel.pass("max-clusters-to-chain", cluster_num);
            }
            
            // Collect some cluster statistics in the graph
            size_t cluster_node_count = 0;
            nid_t cluster_min_node = std::numeric_limits<nid_t>::max();
            nid_t cluster_max_node = 0;
            {
                // Count the distinct node IDs in the cluster (as seed starts)
                // to get an idea of its size in the reference
                std::unordered_set<nid_t> id_set;
                for (auto seed_index : cluster.seeds) {
                    auto& seed = seeds[seed_index];
                    nid_t node_id = id(seed.pos);
                    cluster_min_node = std::min(cluster_min_node, node_id);
                    cluster_max_node = std::max(cluster_max_node, node_id);
                    id_set.insert(node_id);
                }
                cluster_node_count = id_set.size();
            }
            
            // First check against the additional score filter
            if (cfg.cluster_score_cutoff_enabled && cluster.score < cfg.cluster_score_cutoff 
                && kept_cluster_count >= cfg.min_clusters_to_chain) {
                //If the score isn't good enough and we already kept at least cfg.min_clusters_to_chain clusters,
                //ignore this cluster
                if (track_provenance) {
                    funnel.fail("cluster-score", cluster_num, cluster.score);
                }
                if (show_work) {
                    #pragma omp critical (cerr)
                    {
                        cerr << log_name() << "Cluster " << cluster_num << " fails cluster score cutoff" <<  endl;
                        cerr << log_name() << "Covers " << clusters[cluster_num].coverage << "/best-" << cfg.cluster_coverage_threshold << " of read" << endl;
                        cerr << log_name() << "Involves " << cluster_node_count << " nodes in " << cluster_min_node << "-" << cluster_max_node << endl;
                        cerr << log_name() << "Scores " << clusters[cluster_num].score << "/" << cfg.cluster_score_cutoff << endl;
                    }
                }
                return false;
            }
            
            if (track_provenance) {
                funnel.pass("cluster-score", cluster_num, cluster.score);
            }
            

            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Cluster " << cluster_num << endl;
                    cerr << log_name() << "Covers " << cluster.coverage << "/best-" << cfg.cluster_coverage_threshold << " of read" << endl;
                    cerr << log_name() << "Involves " << cluster_node_count << " nodes in " << cluster_min_node << "-" << cluster_max_node << endl;
                    cerr << log_name() << "Scores " << cluster.score << "/" << cfg.cluster_score_cutoff << endl;
                }
            }
            
            if (track_provenance) {
                // Say we're working on this cluster
                funnel.processing_input(cluster_num);
            }
           
            // Count how many of each minimizer is in each cluster that we kept.
            // TODO: deduplicate with extend_cluster
            minimizer_kept_cluster_count.emplace_back(minimizers.size(), 0);
            for (auto seed_index : cluster.seeds) {
                auto& seed = seeds[seed_index];
                minimizer_kept_cluster_count.back()[seed.source]++;
            }
            ++kept_cluster_count;
            
            if (show_work) {
                dump_debug_seeds(minimizers, seeds, cluster.seeds);
            }
           
            // Sort all the seeds used in the cluster by start position, so we can chain them.
            std::vector<size_t> cluster_seeds_sorted = cluster.seeds;
            
            // Sort seeds by read start of seeded region, and remove indexes for seeds that are redundant
            algorithms::sort_and_shadow(seed_anchors, cluster_seeds_sorted);
            
            if (track_provenance) {
                funnel.substage("find_chain");
            }
            
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Computing chain over " << cluster_seeds_sorted.size() << " seeds for cluster " << cluster_num << endl;
                }
            }
            
            if (show_work) {
                // Log the chaining problem so we can try it again elsewhere.
                this->dump_chaining_problem(seed_anchors, cluster_seeds_sorted, gbwt_graph);
            }
            
            // Compute the best chain
            cluster_nums.push_back(cluster_num);
            cluster_chains.emplace_back();
            cluster_chain_seeds.emplace_back();
                
            // Find chains from this cluster
            algorithms::transition_iterator for_each_transition = algorithms::zip_tree_transition_iterator(
                seeds,
                zip_code_tree,
                cfg.max_lookback_bases
            ); 
            VectorView<algorithms::Anchor> cluster_view {seed_anchors, cluster_seeds_sorted};
            std::vector<std::pair<int, std::vector<size_t>>> chains = algorithms::find_best_chains(
                cluster_view,
                *distance_index,
                gbwt_graph,
                get_regular_aligner()->gap_open,
                get_regular_aligner()->gap_extension,
                cfg.max_chains_per_cluster,
                for_each_transition,
                cfg.item_bonus,
                cfg.item_scale,
                cfg.max_indel_bases
            );
            if (show_work) {
                #pragma omp critical (cerr)
                cerr << log_name() << "Asked for " << cfg.max_chains_per_cluster << " and found " << chains.size() << " chains in cluster " << cluster_num << std::endl;
                for (auto& scored_chain : chains) {
                    if (!scored_chain.second.empty()) {
                        #pragma omp critical (cerr)
                        {
                            
                            cerr << log_name() << "Cluster " << cluster_num << " running " << seed_anchors[cluster_seeds_sorted.front()] << " to " << seed_anchors[cluster_seeds_sorted.back()]
                                << " has chain with score " << scored_chain.first
                                << " and length " << scored_chain.second.size()
                                << " running R" << cluster_view[scored_chain.second.front()].read_start()
                                << " to R" << cluster_view[scored_chain.second.back()].read_end() << std::endl;
                        }
                    }
                }
            }
            
            cluster_chains.back() = std::move(chains);
            cluster_chain_seeds.back() = std::move(cluster_seeds_sorted);
            
            if (track_provenance) {
                funnel.substage_stop();
            }
            
            if (track_provenance) {
                for (auto& chain : cluster_chains.back()) {
                    // Record with the funnel that there is now a chain that comes
                    // from all the seeds that participate in the chain.
                    funnel.introduce();
                    funnel.score(funnel.latest(), chain.first);
                    // Accumulate the old and new seed funnel numbers to connect to.
                    // TODO: should we just call into the funnel every time instead of allocating?
                    std::vector<size_t> old_seed_ancestors;
                    std::vector<size_t> new_seed_ancestors;
                    for (auto& sorted_seed_number : chain.second) {
                        // Map each seed back to its canonical seed order
                        size_t seed_number = cluster_chain_seeds.back().at(sorted_seed_number);
                        if (seed_number < old_seed_count) {
                            // Seed is original, from "seed" stage
                            old_seed_ancestors.push_back(seed_number);
                        } else {
                            // Seed is new, from "reseed" stage. Came
                            // after all the fragments which also live in the reseed stage.
                            new_seed_ancestors.push_back(seed_number - old_seed_count + new_seed_start);
                        }
                    }
                    
                    if (!old_seed_ancestors.empty()) {
                        // We came from all the original seeds
                        funnel.also_merge_group(seed_stage_offset, old_seed_ancestors.begin(), old_seed_ancestors.end());
                    }
                    
                    if (!new_seed_ancestors.empty()) {
                        // We came from all the new seeds
                        funnel.also_merge_group(reseed_stage_offset, new_seed_ancestors.begin(), new_seed_ancestors.end());
                    }
                    
                    // We're also related to the source cluster from the
                    // immediately preceeding stage.
                    funnel.also_relevant(1, cluster_num);
                }
                
                // Say we finished with this cluster, for now.
                funnel.processed_input();
            }
            
            return true;
            
        }, [&](size_t cluster_num) -> void {
            // There are too many sufficiently good clusters
            const Cluster& cluster = clusters[cluster_num];
            if (track_provenance) {
                funnel.pass("cluster-coverage", cluster_num, cluster.coverage);
                funnel.fail("max-clusters-to-chain", cluster_num);
            }
            
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    
                    cerr << log_name() << "Cluster " << cluster_num << " passes cluster cutoffs but we have too many" <<  endl;
                    cerr << log_name() << "Covers " << cluster.coverage << "/best-" << cfg.cluster_coverage_threshold << " of read" << endl;
                    cerr << log_name() << "Scores " << cluster.score << "/" << cfg.cluster_score_cutoff << endl;
                }
            }
            
        }, [&](size_t cluster_num) -> void {
            // This cluster is not sufficiently good.
            if (track_provenance) {
                funnel.fail("cluster-coverage", cluster_num, clusters[cluster_num].coverage);
            }
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Cluster " << cluster_num << " fails cluster coverage cutoffs" <<  endl;
                    cerr << log_name() << "Covers " << clusters[cluster_num].coverage << "/best-" << cfg.cluster_coverage_threshold << " of read" << endl;
                    cerr << log_name() << "Scores " << clusters[cluster_num].score << "/" << cfg.cluster_score_cutoff << endl;
                }
            }
        });
   
   // Now give back the chains and the context needed to interpret them.
   return {cluster_nums, cluster_chains, cluster_chain_seeds, seed_anchors, minimizer_explored, minimizer_kept_cluster_count, kept_cluster_count}; 
    
}
    

vector<Alignment> MinimizerMapper::map_from_chains(Alignment& aln) {
    
    if (show_work) {
        #pragma omp critical (cerr)
        dump_debug_query(aln);
    }
    
    // Make a new funnel instrumenter to watch us map this read.
    Funnel funnel;
    funnel.start(aln.name());
    
    // Prepare the RNG for shuffling ties, if needed
    LazyRNG rng([&]() {
        return aln.sequence();
    });


    // Minimizers sorted by position
    std::vector<Minimizer> minimizers_in_read = this->find_minimizers(aln.sequence(), funnel);
    // Indexes of minimizers, sorted into score order, best score first
    std::vector<size_t> minimizer_score_order = sort_minimizers_by_score(minimizers_in_read);
    // Minimizers sorted by best score first
    VectorView<Minimizer> minimizers{minimizers_in_read, minimizer_score_order};

    vector<ZipCodeDecoder> decoders;
    
    // Find the seeds and mark the minimizers that were located.
    vector<Seed> seeds = this->find_seeds(minimizers_in_read, minimizers, aln, funnel);
    // Make them into a zip code tree
    ZipCodeTree zip_code_tree;
    crash_unless(distance_index);
    zip_code_tree.fill_in_tree(seeds, *distance_index);
    
    if (show_work) {
        #pragma omp critical (cerr)
        {
            std::cerr << log_name() << "Zip code tree:";
            zip_code_tree.print_self();
        }
    }

    // Pre-cluster just the seeds we have. Get sets of input seed indexes that go together.
    if (track_provenance) {
        funnel.stage("bucket");
        funnel.substage("compute-buckets");
    }

    // Bucket the hits coarsely into sets that might be able to interact.
    std::vector<Cluster> buckets = clusterer.cluster_seeds(seeds, aln.sequence().size() * bucket_scale);
    //std::vector<Cluster> buckets = zip_clusterer.coarse_cluster_seeds(seeds, aln.sequence().size() * bucket_scale);
    
    // Score all the buckets
    if (track_provenance) {
        funnel.substage("score-buckets");
    }
    double best_bucket_score = 0;
    double second_best_bucket_score = 0;
    for (size_t i = 0; i < buckets.size(); i++) {
        Cluster& bucket = buckets[i];
        
        if (this->track_provenance) {
            // Say we're making it
            funnel.producing_output(i);
        } 
        this->score_cluster(bucket, i, minimizers, seeds, aln.sequence().size());
        if (bucket.score > best_bucket_score) {
            second_best_bucket_score = best_bucket_score;
            best_bucket_score = bucket.score;
        } else if (bucket.score > second_best_bucket_score) {
            second_best_bucket_score = bucket.score;
        }
        if (this->track_provenance) {
            // Record the cluster in the funnel as a group of the size of the number of items.
            funnel.merge_group(bucket.seeds.begin(), bucket.seeds.end());
            funnel.score(funnel.latest(), bucket.score);

            if (show_work) {
                auto bucket_positions = funnel.get_positions(funnel.latest());
                #pragma omp critical (cerr)
                {
                    std::cerr << log_name() << "Positions for bucket " << i << ":" << std::endl;
                    for (auto& handle_and_range : bucket_positions) {
                        // Log each range on a path associated with the bucket.
                        std::cerr << log_name() << "\t"
                            << this->path_graph->get_path_name(handle_and_range.first)
                            << ":" << handle_and_range.second.first
                            << "-" << handle_and_range.second.second << std::endl;
                    }
                }
            }

            // Say we made it.
            funnel.produced_output();
        }
    }
    
    // Now we need to chain into fragments.
    // Each fragment needs to end up with a seeds array of seed numbers, and a
    // coverage float on the read, just like a cluster, for downstream
    // processing.
    if (track_provenance) {
        funnel.stage("fragment");
        funnel.substage("fragment");
    }
    
    chain_config_t fragment_cfg;
    
    // Make fragments be compact
    fragment_cfg.max_lookback_bases = this->fragment_max_lookback_bases;
    fragment_cfg.min_lookback_items = this->fragment_min_lookback_items;
    fragment_cfg.lookback_item_hard_cap = this->fragment_lookback_item_hard_cap;
    fragment_cfg.initial_lookback_threshold = this->initial_lookback_threshold;
    fragment_cfg.lookback_scale_factor = this->lookback_scale_factor;
    fragment_cfg.min_good_transition_score_per_base = this->min_good_transition_score_per_base;
    
    fragment_cfg.item_bonus = this->item_bonus;
    fragment_cfg.item_scale = this->item_scale;
    fragment_cfg.max_indel_bases = this->fragment_max_indel_bases;
    
    // Do all the ones that are 75% as good as the best, or down to 50% as good
    // as the best if that is what it takes to get the second best
    double bucket_score_cutoff = best_bucket_score / 0.75;
    if (bucket_score_cutoff - (bucket_score_cutoff / 0.25) < second_best_bucket_score) {
        bucket_score_cutoff = std::min(bucket_score_cutoff, second_best_bucket_score);
    }
    fragment_cfg.cluster_score_cutoff = bucket_score_cutoff;
    fragment_cfg.cluster_score_cutoff_enabled = true;
    fragment_cfg.cluster_coverage_threshold = 1.0;
    fragment_cfg.min_clusters_to_chain = this->min_buckets_to_fragment;
    fragment_cfg.max_clusters_to_chain = this->max_buckets_to_fragment;
    
    fragment_cfg.max_chains_per_cluster = this->max_fragments_per_bucket;
    
    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "=====Creating fragments=====" << endl;
        }
    }

    // Go get fragments from the buckets. Note that this doesn't process all buckets! It will really only do the best ones!
    auto fragment_results = this->chain_clusters(aln, minimizers, seeds, zip_code_tree, buckets, fragment_cfg, seeds.size(), seeds.size(), funnel, 2, std::numeric_limits<size_t>::max(), rng);
    
    if (track_provenance) {
        funnel.substage("translate-fragments");
    }
   
    // Turn fragments into several corresponding lists.
    // What seeds are visited in what order in the fragment?
    std::vector<std::vector<size_t>> fragments;
    // What score does each fragment have?
    std::vector<double> fragment_scores;
    // Which bucket did each fragment come from (for stats)
    std::vector<size_t> fragment_source_bucket;
    // How many of each minimizer ought to be considered explored by each fragment?
    std::vector<std::vector<size_t>> minimizer_kept_fragment_count;
   
    for (size_t i = 0; i < fragment_results.cluster_chains.size(); i++) {
        // For each source bucket (in exploration order)
        for (auto& chain : fragment_results.cluster_chains[i]) {
            // For each fragment found in the bucket
        
            // Convert format
            fragments.emplace_back();
        
            if (this->track_provenance) {
                // Say we're making it
                funnel.producing_output(fragments.size());
            } 
            // Copy all the seeds in the chain over
            fragments.back().reserve(chain.second.size());
            for (auto& chain_visited_index : chain.second) {
                // Make sure to translate to real seed space
                fragments.back().push_back(fragment_results.cluster_chain_seeds[i].at(chain_visited_index));
            }
            
            // Record score
            fragment_scores.push_back(chain.first);
            
            // Work out the source bucket (in bucket order) that the fragment came from
            size_t source_bucket = fragment_results.cluster_nums.at(i);
            if (this->track_provenance) {
                // Record the fragment in the funnel as coming from the bucket 
                funnel.project(source_bucket);
                funnel.score(funnel.latest(), chain.first);

                // Say we made it.
                funnel.produced_output();
            }
            
            // Remember outside the funnel what bucket it came from, for statistics
            fragment_source_bucket.push_back(source_bucket);
            
            // Remember how many of each minimizer's hits were in the bucket for each fragment. These are ordered by visited bucket, so index with i.
            // TODO: Is there a faster way to do this? Do we even care about this for MAPQ anymore? 
            minimizer_kept_fragment_count.push_back(fragment_results.minimizer_kept_cluster_count.at(i));
        }
    }
    
    // Now glom the fragments together into chains 
    if (track_provenance) {
        funnel.stage("chain");
        funnel.substage("fragment-stats");
    }
    
    // Select the "best" bucket.
    // Bucket with the best fragment score
    size_t best_bucket = 0;
    // That fragment
    size_t best_fragment = 0;
    // That score
    double best_bucket_fragment_score = 0;
    for (size_t i = 0; i < fragment_scores.size(); i++) {
        if (fragment_scores[i] >= best_bucket_fragment_score) {
            best_bucket_fragment_score = fragment_scores[i];
            best_fragment = i;
            best_bucket = fragment_source_bucket[i];
        }
    }
    if (show_work) {
        #pragma omp critical (cerr)
        std::cerr << log_name() << "Bucket " << best_bucket << " is best with fragment " << best_fragment << " with score " << best_bucket_fragment_score << std::endl;
    }
    size_t best_bucket_seed_count = buckets.at(best_bucket).seeds.size();

    // Count up all the minimizers in the best bucket
    size_t best_bucket_minimizer_count;
    {
        std::unordered_set<size_t> best_bucket_minimizers;
        for (auto& seed : buckets.at(best_bucket).seeds) {
            best_bucket_minimizers.insert(seeds.at(seed).source);
        }
        best_bucket_minimizer_count = best_bucket_minimizers.size();
    }

    if (show_work) {
        // Dump the best bucket's best fragment
        dump_debug_dotplot("best-fragment", "fragment", minimizers, seeds, buckets.at(best_bucket).seeds, fragments.at(best_fragment), this->path_graph);
    }
    
    // Find the fragments that are in the best bucket
    std::vector<size_t> best_bucket_fragments;
    for (size_t i = 0; i < fragments.size(); i++) {
        if (show_work) {
            #pragma omp critical (cerr)
            std::cerr << log_name() << "Fragment " << i << " with score " << fragment_scores.at(i) << " came from bucket " << fragment_source_bucket.at(i) << std::endl;
        }
        if (fragment_source_bucket.at(i) == best_bucket) {
            // Get all the fragment indexes that are from the best bucket
            best_bucket_fragments.push_back(i);
        }
    }
    
    // Sort fragments in best bucket by score, descending
    std::sort(best_bucket_fragments.begin(), best_bucket_fragments.end(), [&](const size_t& a, const size_t& b) {
        // Return true if a has a larger score and should come before b.
        // Make sure to use chaining scores and not scores as clusters.
        return fragment_scores.at(a) > fragment_scores.at(b);
        
    });
    
    // Work out of read with top k fragments by score, in best bucket
    const size_t TOP_FRAGMENTS = 4;
    std::vector<double> best_bucket_fragment_coverage_at_top(TOP_FRAGMENTS + 1, 0.0);
    for (size_t fragment_count = 0; fragment_count <= TOP_FRAGMENTS && fragment_count < fragments.size(); fragment_count++) {
        // Do O(n^2) easy way to compute coverage in top k fragments up to this many.
        std::vector<size_t> top_fragments;
        top_fragments.reserve(fragment_count);
        for (size_t i = 0; i < fragment_count && i < best_bucket_fragments.size(); i++) {
            top_fragments.push_back(best_bucket_fragments.at(i));
        }
        best_bucket_fragment_coverage_at_top[fragment_count] = get_read_coverage(aln, {fragments, top_fragments}, seeds, minimizers);
    }
    
    if (track_provenance) {
        funnel.substage("chain");
    }
    
    // For each chain, we need:
    // The chain itself, pointing into seeds
    std::vector<std::vector<size_t>> chains;
    // The bucket it came from
    std::vector<size_t> chain_source_buckets;
    // An estimated alignment score
    std::vector<int> chain_score_estimates;
    // A count, for each minimizer, of how many hits of it could have been in the chain, or were considered when making the chain.
    std::vector<std::vector<size_t>> minimizer_kept_chain_count;
    
    // We also need a set of anchors for all the seeds. We will extend this if we reseed more seeds.
    std::vector<algorithms::Anchor>& seed_anchors = fragment_results.seed_anchors;
    
    // Make a list of anchors where we have each fragment as itself an anchor
    std::vector<algorithms::Anchor> fragment_anchors;
    fragment_anchors.reserve(fragments.size());
    for (size_t i = 0; i < fragments.size(); i++) {
        auto& fragment = fragments.at(i);
        auto& score = fragment_scores.at(i);
        fragment_anchors.push_back(algorithms::Anchor(seed_anchors.at(fragment.front()), seed_anchors.at(fragment.back()), score));
    }
    
    // Get all the fragment numbers for each bucket we actually used, so we can chain each bucket independently again.
    // TODO: Stop reswizzling so much.
    std::unordered_map<size_t, std::vector<size_t>> bucket_fragment_nums;
    for (size_t i = 0; i < fragment_source_bucket.size(); i++) {
        bucket_fragment_nums[fragment_source_bucket[i]].push_back(i);
    }
    
    // Get the score of the top-scoring fragment per bucket.
    std::unordered_map<size_t, double> bucket_best_fragment_score;
    for (auto& kv : bucket_fragment_nums) {
        for (auto& fragment_num : kv.second) {
            // Max in the score of each fragmrnt in the bucket
            bucket_best_fragment_score[kv.first] = std::max(bucket_best_fragment_score[kv.first], fragment_scores.at(fragment_num));
        }
    }
    
    // Filter down to just the good ones, sorted by read start
    std::unordered_map<size_t, std::vector<size_t>> bucket_good_fragment_nums;
    for (auto& kv : bucket_fragment_nums) {
        // Decide on how good fragments have to be to keep.
        double fragment_score_threshold = bucket_best_fragment_score.at(kv.first) * fragment_score_fraction;
    
        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "Keeping, of the " << kv.second.size() << " fragments in bucket " << kv.first << ", those with score of at least "  << fragment_score_threshold << endl;
            }
        }
    
        // Keep the fragments that have good scores.
        for (auto& fragment_num : kv.second) {
            // For each fragment in the bucket
            if (fragment_scores.at(fragment_num) >= fragment_score_threshold) {
                // If its score is high enough, keep it.
                // TODO: Tell the funnel.
                bucket_good_fragment_nums[kv.first].push_back(fragment_num);
            }
        }
        
        // Now sort anchors by read start. Don't bother with shadowing.
        algorithms::sort_anchor_indexes(fragment_anchors, bucket_good_fragment_nums[kv.first]);

        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "\tKept " << bucket_good_fragment_nums[kv.first].size() << " fragments." << endl;
            }
        }
    }
    
    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "=====Creating chains=====" << endl;
        }
    }

    for (auto& kv : bucket_good_fragment_nums) {
        auto& bucket_num = kv.first;
        // Get a view of all the good fragments in the bucket.
        // TODO: Should we just not make a global fragment anchor list?
        VectorView<algorithms::Anchor> bucket_fragment_view {fragment_anchors, kv.second};

        if (bucket_fragment_view.empty()) {
            // Nothing to chain!
            if (show_work) {
                #pragma omp critical (cerr)
                std::cerr << log_name() << "Bucket " << bucket_num << " has no good fragments to chain!" << std::endl;
            } 
            continue;
        }
        
        if (show_work) {
            #pragma omp critical (cerr)
            std::cerr << log_name() << "Chaining bucket " << bucket_num << std::endl;
        } 

        // Chain up the fragments
        algorithms::transition_iterator for_each_transition = algorithms::zip_tree_transition_iterator(
            seeds,
            zip_code_tree,
            this->max_lookback_bases
        ); 
        std::vector<std::pair<int, std::vector<size_t>>> chain_results = algorithms::find_best_chains(
            bucket_fragment_view,
            *distance_index,
            gbwt_graph,
            get_regular_aligner()->gap_open,
            get_regular_aligner()->gap_extension,
            2,
            for_each_transition,
            this->item_bonus,
            this->item_scale,
            this->max_indel_bases
        );
        
        for (auto& chain_result: chain_results) {
            // Each chain of fragments becomes a chain of seeds
            chains.emplace_back();
            auto& chain = chains.back();
            // With a bucket
            chain_source_buckets.push_back(bucket_num);
            // With a score
            chain_score_estimates.emplace_back(0);
            int& score = chain_score_estimates.back();
            // And counts of each minimizer kept
            minimizer_kept_chain_count.emplace_back();
            auto& minimizer_kept = minimizer_kept_chain_count.back();
            
            // We record the fragments that merge into each chain for reporting.
            std::vector<size_t> chain_fragment_nums_overall;
            chain_fragment_nums_overall.reserve(chain_result.second.size());
            
            for (const size_t& fragment_in_bucket: chain_result.second) {
                // For each fragment in the chain
                           
                // Get its fragment number out of all fragments
                size_t fragment_num_overall = kv.second.at(fragment_in_bucket);
                
                // Save it
                chain_fragment_nums_overall.push_back(fragment_num_overall);
                
                // Go get that fragment
                auto& fragment = fragments.at(fragment_num_overall);
                    
                // And append all the seed numbers to the chain
                std::copy(fragment.begin(), fragment.end(), std::back_inserter(chain));
                
                // And count the score
                score += fragment_scores.at(fragment_num_overall);
                
                // And count the kept minimizers
                auto& fragment_minimizer_kept = minimizer_kept_fragment_count.at(fragment_num_overall);
                if (minimizer_kept.size() < fragment_minimizer_kept.size()) {
                    minimizer_kept.resize(fragment_minimizer_kept.size());
                }
                for (size_t i = 0; i < fragment_minimizer_kept.size(); i++) {
                    minimizer_kept[i] += fragment_minimizer_kept[i];
                }
            }
            if (track_provenance) {
                // Say all those fragments became a chain
                funnel.merge_group(chain_fragment_nums_overall.begin(), chain_fragment_nums_overall.end());
                // With the total score
                funnel.score(funnel.latest(), score);
            }
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    std::cerr << log_name() << "Chain " << (chains.size() - 1) << " with score " << score << " is composed from fragments:";
                    for (auto& f : chain_fragment_nums_overall) {
                        std::cerr << " " << f;
                    } 
                    std::cerr << std::endl;
                }
            } 
        }
    }
    
    // Find the best chain
    size_t best_chain = std::numeric_limits<size_t>::max();
    int best_chain_score = 0;
    for (size_t i = 0; i < chains.size(); i++) {
        if (best_chain == std::numeric_limits<size_t>::max() || chain_score_estimates.at(i) > best_chain_score) {
            // Friendship ended with old chain
            best_chain = i;
            best_chain_score = chain_score_estimates[i];
        }
    }
    bool best_chain_correct = false;
    if (track_correctness && best_chain != std::numeric_limits<size_t>::max()) {
        // We want to explicitly check if the best chain was correct, for looking at stats about it later.
        if (funnel.is_correct(best_chain)) {
            best_chain_correct = true;
        }
    }

    if (show_work) {
        // Dump the best chain
        dump_debug_dotplot("best-chain", "chain", minimizers, seeds, buckets.at(chain_source_buckets.at(best_chain)).seeds, chains.at(best_chain), this->path_graph);
    }
    
    // Find its coverage
    double best_chain_coverage = get_read_coverage(aln, std::vector<std::vector<size_t>> {chains.at(best_chain)}, seeds, minimizers);
    
    // Find out how gappy it is. We can get the longest and the average distance maybe.
    size_t best_chain_longest_jump = 0;
    size_t best_chain_total_jump = 0;
    for (size_t i = 1; i < chains.at(best_chain).size(); i++) {
        // Find the pair of anchors we go between
        auto& left_anchor = seed_anchors.at(chains.at(best_chain).at(i - 1));
        auto& right_anchor = seed_anchors.at(chains.at(best_chain).at(i));
        // And get the distance between them in the read
        size_t jump = right_anchor.read_start() - left_anchor.read_end();
        // Max and add it in
        best_chain_longest_jump = std::max(best_chain_longest_jump, jump);
        best_chain_total_jump += jump;
    }
    double best_chain_average_jump = chains.at(best_chain).size() > 1 ? best_chain_total_jump / (chains.at(best_chain).size() - 1) : 0.0;

    // Also count anchors in the chain
    size_t best_chain_anchors = chains.at(best_chain).size();

    // And total length of anchors in the chain
    size_t best_chain_anchor_length = 0;
    for (auto& item : chains.at(best_chain)) {
        best_chain_anchor_length += seed_anchors.at(item).length();
    }
    
    // Now do reseeding inside chains. Not really properly a funnel stage; it elaborates the chains
    if (track_provenance) {
        funnel.substage("reseed");
    }
    
    // Remember how many seeds we had before reseeding
    size_t old_seed_count = seeds.size();
    
    // We are going to need a widget for finding minimizer hit
    // positions in a subgraph, in the right orientation.
    auto find_minimizer_hit_positions = [&](const Minimizer& m, const vector<id_t>& sorted_ids, const std::function<void(const pos_t)>& iteratee) -> void {
        gbwtgraph::hits_in_subgraph(m.hits, m.occs, sorted_ids, [&](pos_t pos, gbwtgraph::Payload) {
            if (m.value.is_reverse) {
                // Convert to face along forward strand of read.
                size_t node_length = this->gbwt_graph.get_length(this->gbwt_graph.get_handle(id(pos)));
                pos = reverse_base_pos(pos, node_length);
            }
            // Show the properly stranded position to the iteratee.
            iteratee(pos);
        });
    };
    
    // We are going to need our existing seeds in the form of something we can deduplicate.
    // TODO: Also remove overlap?
    std::unordered_set<std::pair<size_t, pos_t>> seen_seeds;
    for (auto& seed : seeds) {
        seen_seeds.emplace(minimizers[seed.source].forward_offset(), seed.pos);
    }

    // TODO: Do any reseeding. For now we do none.
    // TODO: Rescore the reseeded chains.
    
    if (track_provenance) {
        funnel.stage("align");
    }

    //How many of each minimizer ends up in a chain that actually gets turned into an alignment?
    vector<size_t> minimizer_kept_count(minimizers.size(), 0);
    
    // Now start the alignment step. Everything has to become an alignment.

    // We will fill this with all computed alignments in estimated score order.
    vector<Alignment> alignments;
    alignments.reserve(chain_score_estimates.size());
    // This maps from alignment index back to chain index, for
    // tracing back to minimizers for MAPQ. Can hold
    // numeric_limits<size_t>::max() for an unaligned alignment.
    vector<size_t> alignments_to_source;
    alignments_to_source.reserve(chain_score_estimates.size());

    // Create a new alignment object to get rid of old annotations.
    {
      Alignment temp;
      temp.set_sequence(aln.sequence());
      temp.set_name(aln.name());
      temp.set_quality(aln.quality());
      aln = std::move(temp);
    }

    // Annotate the read with metadata
    if (!sample_name.empty()) {
        aln.set_sample_name(sample_name);
    }
    if (!read_group.empty()) {
        aln.set_read_group(read_group);
    }
    
    // We need to be able to discard a chain because its score isn't good enough.
    // We have more components to the score filter than process_until_threshold_b supports.
    auto discard_chain_by_score = [&](size_t processed_num) -> void {
        // This chain is not good enough.
        if (track_provenance) {
            funnel.fail("chain-score", processed_num, chain_score_estimates[processed_num]);
        }
        
        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "chain " << processed_num << " failed because its score was not good enough (score=" << chain_score_estimates[processed_num] << ")" << endl;
                if (track_correctness && funnel.was_correct(processed_num)) {
                    cerr << log_name() << "\tCORRECT!" << endl;
                }
            }
        }
    };
    
    // Track if minimizers were explored by alignments
    SmallBitset minimizer_explored(minimizers.size());
    
    // Go through the chains in estimated-score order.
    process_until_threshold_b<int>(chain_score_estimates,
        chain_score_threshold, min_chains, max_alignments, rng, [&](size_t processed_num) -> bool {
            // This chain is good enough.
            // Called in descending score order.
            
            if (chain_score_estimates[processed_num] < chain_min_score) {
                // Actually discard by score
                discard_chain_by_score(processed_num);
                return false;
            }
            
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "chain " << processed_num << " is good enough (score=" << chain_score_estimates[processed_num] << ")" << endl;
                    if (track_correctness && funnel.was_correct(processed_num)) {
                        cerr << log_name() << "\tCORRECT!" << endl;
                    }
                }
            }
            if (track_provenance) {
                funnel.pass("chain-score", processed_num, chain_score_estimates[processed_num]);
                funnel.pass("max-alignments", processed_num);
                funnel.processing_input(processed_num);
            }

            // Collect the top alignments. Make sure we have at least one always, starting with unaligned.
            vector<Alignment> best_alignments(1, aln);

            // Align from the chained-up seeds
            if (do_dp) {
                // We need to do base-level alignment.
            
                if (track_provenance) {
                    funnel.substage("align");
                }
                
                // We currently just have the one best score and chain per cluster
                vector<size_t>& chain = chains[processed_num];
                
                try {
                    // Do the DP between the items in the chain. 
                    best_alignments[0] = find_chain_alignment(aln, seed_anchors, chain);
                } catch (ChainAlignmentFailedError& e) {
                    // We can't actually make an alignment from this chain
                    #pragma omp critical (cerr)
                    cerr << log_name() << "Error creating alignment from chain for " << aln.name() << ": " << e.what() << endl;
                    // Leave the read unmapped.
                }
                    
                // TODO: Come up with a good secondary for the cluster somehow.
            } else {
                // We would do base-level alignment but it is disabled.
                // Leave best_alignment unaligned
            }
           
            // Have a function to process the best alignments we obtained
            auto observe_alignment = [&](Alignment& aln) {
                alignments.emplace_back(std::move(aln));
                alignments_to_source.push_back(processed_num);

                if (track_provenance) {
    
                    funnel.project(processed_num);
                    funnel.score(alignments.size() - 1, alignments.back().score());
                }
                if (show_work) {
                    #pragma omp critical (cerr)
                    {
                        cerr << log_name() << "Produced alignment from chain " << processed_num
                            << " with score " << alignments.back().score() << ": " << log_alignment(alignments.back()) << endl;
                    }
                }
            };
            
            if (!best_alignments.empty() && best_alignments[0].score() <= 0) {
                if (show_work) {
                    // Alignment won't be observed but log it anyway.
                    #pragma omp critical (cerr)
                    {
                        cerr << log_name() << "Produced terrible best alignment from chain " << processed_num << ": " << log_alignment(best_alignments[0]) << endl;
                    }
                }
            }
            for(auto aln_it = best_alignments.begin() ; aln_it != best_alignments.end() && aln_it->score() != 0 && aln_it->score() >= best_alignments[0].score() * 0.8; ++aln_it) {
                //For each additional alignment with score at least 0.8 of the best score
                observe_alignment(*aln_it);
            }
           
            if (track_provenance) {
                // We're done with this input item
                funnel.processed_input();
            }

            for (size_t i = 0 ; i < minimizer_kept_chain_count[processed_num].size() ; i++) {
                minimizer_kept_count[i] += minimizer_kept_chain_count[processed_num][i];
                if (minimizer_kept_chain_count[processed_num][i] > 0) {
                    // This minimizer is in a cluster that gave rise
                    // to at least one alignment, so it is explored.
                    minimizer_explored.insert(i);
                }
            }
            
            return true;
        }, [&](size_t processed_num) -> void {
            // There are too many sufficiently good chains
            if (track_provenance) {
                funnel.pass("chain-score", processed_num, chain_score_estimates[processed_num]);
                funnel.fail("max-alignments", processed_num);
            }
            
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "chain " << processed_num << " failed because there were too many good chains (score=" << chain_score_estimates[processed_num] << ")" << endl;
                    if (track_correctness && funnel.was_correct(processed_num)) {
                        cerr << log_name() << "\tCORRECT!" << endl;
                    }
                }
            }
        }, discard_chain_by_score);
    
    if (alignments.size() == 0) {
        // Produce an unaligned Alignment
        alignments.emplace_back(aln);
        alignments_to_source.push_back(numeric_limits<size_t>::max());
        
        if (track_provenance) {
            // Say it came from nowhere
            funnel.introduce();
        }
    }
    
    if (track_provenance) {
        // Now say we are finding the winner(s)
        funnel.stage("winner");
    }
    
    // Fill this in with the alignments we will output as mappings
    vector<Alignment> mappings;
    mappings.reserve(min(alignments.size(), max_multimaps));
    
    // Grab all the scores in order for MAPQ computation.
    vector<double> scores;
    scores.reserve(alignments.size());
    
    process_until_threshold_a(alignments.size(), (std::function<double(size_t)>) [&](size_t i) -> double {
        return alignments.at(i).score();
    }, 0, 1, max_multimaps, rng, [&](size_t alignment_num) {
        // This alignment makes it
        // Called in score order
        
        // Remember the score at its rank
        scores.emplace_back(alignments[alignment_num].score());
        
        // Remember the output alignment
        mappings.emplace_back(std::move(alignments[alignment_num]));
        
        if (track_provenance) {
            // Tell the funnel
            funnel.pass("max-multimaps", alignment_num);
            funnel.project(alignment_num);
            funnel.score(funnel.latest(), scores.back());
        }
        
        return true;
    }, [&](size_t alignment_num) {
        // We already have enough alignments, although this one has a good score
        
        // Remember the score at its rank anyway
        scores.emplace_back(alignments[alignment_num].score());
        
        if (track_provenance) {
            funnel.fail("max-multimaps", alignment_num);
        }
    }, [&](size_t alignment_num) {
        // This alignment does not have a sufficiently good score
        // Score threshold is 0; this should never happen
        crash_unless(false);
    });
    
    if (track_provenance) {
        funnel.substage("mapq");
    }

    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "Picked best alignment " << log_alignment(mappings[0]) << endl;
            cerr << log_name() << "For scores";
            for (auto& score : scores) cerr << " " << score << ":" << endl;
        }
    }

    crash_unless(!mappings.empty());
    // Compute MAPQ if not unmapped. Otherwise use 0 instead of the 50% this would give us.
    // Use exact mapping quality 
    double mapq = (mappings.front().path().mapping_size() == 0) ? 0 : 
        get_regular_aligner()->compute_max_mapping_quality(scores, false) ;

#ifdef print_minimizer_table
    double uncapped_mapq = mapq;
#endif
    
    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "uncapped MAPQ is " << mapq << endl;
        }
    }
    
    // TODO: give SmallBitset iterators so we can use it instead of an index vector.
    vector<size_t> explored_minimizers;
    for (size_t i = 0; i < minimizers.size(); i++) {
        if (minimizer_explored.contains(i)) {
            explored_minimizers.push_back(i);
        }
    }
    // Compute caps on MAPQ. TODO: avoid needing to pass as much stuff along.
    double escape_bonus = mapq < std::numeric_limits<int32_t>::max() ? 1.0 : 2.0;
    double mapq_explored_cap = escape_bonus * faster_cap(minimizers, explored_minimizers, aln.sequence(), aln.quality());

    // Remember the uncapped MAPQ and the caps
    set_annotation(mappings.front(),"secondary_scores", scores);
    set_annotation(mappings.front(), "mapq_uncapped", mapq);
    set_annotation(mappings.front(), "mapq_explored_cap", mapq_explored_cap);

    // Apply the caps and transformations
    mapq = round(min(mapq_explored_cap, min(mapq, 60.0)));

    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "Explored cap is " << mapq_explored_cap << endl;
            cerr << log_name() << "MAPQ is " << mapq << endl;
        }
    }
        
    // Make sure to clamp 0-60.
    mappings.front().set_mapping_quality(max(min(mapq, 60.0), 0.0));
   
    
    if (track_provenance) {
        funnel.substage_stop();
    }
    
    for (size_t i = 0; i < mappings.size(); i++) {
        // For each output alignment in score order
        auto& out = mappings[i];
        
        // Assign primary and secondary status
        out.set_is_secondary(i > 0);
    }
    
    // Stop this alignment
    funnel.stop();
    
    // Annotate with whatever's in the funnel
    funnel.annotate_mapped_alignment(mappings[0], track_correctness);
    
    if (track_provenance) {
        if (track_correctness) {
            annotate_with_minimizer_statistics(mappings[0], minimizers, seeds, old_seed_count, fragments.size(), funnel);
        }
        // Annotate with parameters used for the filters and algorithms.
        
        set_annotation(mappings[0], "param_hit-cap", (double) hit_cap);
        set_annotation(mappings[0], "param_hard-hit-cap", (double) hard_hit_cap);
        set_annotation(mappings[0], "param_score-fraction", (double) minimizer_score_fraction);
        set_annotation(mappings[0], "param_max-unique-min", (double) max_unique_min);
        set_annotation(mappings[0], "param_num-bp-per-min", (double) num_bp_per_min);
        set_annotation(mappings[0], "param_exclude-overlapping-min", exclude_overlapping_min);
        set_annotation(mappings[0], "param_align-from-chains", align_from_chains);
        set_annotation(mappings[0], "param_chaining-cluster-distance", (double) chaining_cluster_distance);
        set_annotation(mappings[0], "param_fragment-connection-coverage-threshold", fragment_connection_coverage_threshold);
        set_annotation(mappings[0], "param_min-fragment-connections", (double) min_fragment_connections);
        set_annotation(mappings[0], "param_max-fragment-connections", (double) max_fragment_connections);
        set_annotation(mappings[0], "param_min-buckets-to-fragment", (double) min_buckets_to_fragment);
        set_annotation(mappings[0], "param_max-buckets-to-fragment", (double) max_buckets_to_fragment);
        set_annotation(mappings[0], "param_reseed-search-distance", (double) reseed_search_distance);
        
        // Chaining algorithm parameters
        set_annotation(mappings[0], "param_max-lookback-bases", (double) max_lookback_bases);
        set_annotation(mappings[0], "param_initial-lookback-threshold", (double) initial_lookback_threshold);
        set_annotation(mappings[0], "param_lookback-scale-factor", lookback_scale_factor);
        set_annotation(mappings[0], "param_min-good-transition-score-per-base", min_good_transition_score_per_base);
        set_annotation(mappings[0], "param_item-bonus", (double) item_bonus);
        set_annotation(mappings[0], "param_item-scale", (double) item_scale);
        set_annotation(mappings[0], "param_max-indel-bases", (double) max_indel_bases);
        
        set_annotation(mappings[0], "param_max-chain-connection", (double) max_chain_connection);
        set_annotation(mappings[0], "param_max-tail-length", (double) max_tail_length);
        set_annotation(mappings[0], "param_max-alignments", (double) max_alignments);
        set_annotation(mappings[0], "param_cluster-score", (double) cluster_score_threshold);
        set_annotation(mappings[0], "param_cluster-coverage", (double) cluster_coverage_threshold);
        set_annotation(mappings[0], "param_chain-score", (double) chain_score_threshold);
        set_annotation(mappings[0], "param_chain-min-score", (double) chain_min_score);
        set_annotation(mappings[0], "param_min-chains", (double) min_chains);
        
    }
    
    // Special fragment and chain statistics
    set_annotation(mappings[0], "fragment_scores", fragment_scores);
    set_annotation(mappings[0], "best_bucket_fragment_coverage_at_top", best_bucket_fragment_coverage_at_top);
    set_annotation(mappings[0], "best_bucket_seed_count", (double)best_bucket_seed_count);
    set_annotation(mappings[0], "best_bucket_minimizer_count", (double)best_bucket_minimizer_count);
    if (track_correctness) {
        set_annotation(mappings[0], "best_chain_correct", best_chain_correct);
    }
    set_annotation(mappings[0], "best_chain_coverage", best_chain_coverage);
    set_annotation(mappings[0], "best_chain_longest_jump", (double) best_chain_longest_jump);
    set_annotation(mappings[0], "best_chain_average_jump", best_chain_average_jump);
    set_annotation(mappings[0], "best_chain_anchors", (double) best_chain_anchors);
    set_annotation(mappings[0], "best_chain_anchor_length", (double) best_chain_anchor_length);
    
#ifdef print_minimizer_table
    cerr << aln.sequence() << "\t";
    for (char c : aln.quality()) {
        cerr << (char)(c+33);
    }
    cerr << "\t" << clusters.size();
    for (size_t i = 0 ; i < minimizers.size() ; i++) {
        auto& minimizer = minimizers[i];
        cerr << "\t"
             << minimizer.value.key.decode(minimizer.length) << "\t"
             << minimizer.forward_offset() << "\t"
             << minimizer.agglomeration_start << "\t"
             << minimizer.agglomeration_length << "\t"
             << minimizer.hits << "\t"
             << minimizer_kept_count[i];
         if (minimizer_kept_count[i]>0) {
             assert(minimizer.hits<=hard_hit_cap) ;
         }
    }
    cerr << "\t" << uncapped_mapq << "\t" << mapq_explored_cap << "\t"  << mappings.front().mapping_quality() << "\t";
    cerr << "\t";
    for (auto& score : scores) {
        cerr << score << ",";
    }
    if (track_correctness) {
        cerr << "\t" << funnel.last_correct_stage() << endl;
    } else {
        cerr << "\t" << "?" << endl;
    }
#endif

    if (track_provenance) {
        if (show_work && aln.sequence().size() < LONG_LIMIT) {
            // Dump the funnel info graph to standard error
            #pragma omp critical (cerr)
            {
                funnel.to_dot(cerr);
            }
        }
        
        // Otherwise/also, if we are dumping explanations, dump it to a file
        DotDumpExplainer<Funnel> explainer(true, funnel);
    }

    return mappings;
}

double MinimizerMapper::get_read_coverage(
    const Alignment& aln,
    const VectorView<std::vector<size_t>>& seed_sets,
    const std::vector<Seed>& seeds,
    const VectorView<Minimizer>& minimizers) const {
    
    std::vector<bool> covered(aln.sequence().size(), false);
    
    for (auto& list : seed_sets) {
        // We will fill in the range it occupies in the read
        std::pair<size_t, size_t> read_range {std::numeric_limits<size_t>::max(), 0};
        
        for (auto& seed_index : list) {
            // Which means we look at the minimizer for each seed
            auto& seed = seeds.at(seed_index);
            crash_unless(seed.source < minimizers.size());
            auto& minimizer = minimizers[seed.source];
            
            if (minimizer.forward_offset() < read_range.first) {
                // Min all their starts to get the start
                read_range.first = minimizer.forward_offset();
            }
            
            if (minimizer.forward_offset() + minimizer.length > read_range.second) {
                // Max all their past-ends to get the past-end
                read_range.second = minimizer.forward_offset() + minimizer.length;
            }
        }
        
        // Then mark its coverage
        set_coverage_flags(covered, read_range.first, read_range.second);
    }
    
    // And return the fraction covered.
    return get_fraction_covered(covered);
}

#define debug_chaining
Alignment MinimizerMapper::find_chain_alignment(
    const Alignment& aln,
    const VectorView<algorithms::Anchor>& to_chain,
    const std::vector<size_t>& chain) const {
    
    if (chain.empty()) {
        throw ChainAlignmentFailedError("Cannot find an alignment for an empty chain!");
    }
    
    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "Align chain of";
            if (chain.size() < MANY_LIMIT) {
                cerr << ": ";
                for (auto item_number : chain) {
                    cerr << " " << item_number;
                }
            } else {
                cerr << " " << chain.size() << " items";
            }
            cerr << " in " << to_chain.size() << " items" << endl;
        }
    }
    
    // We need an Aligner for scoring.
    const Aligner& aligner = *get_regular_aligner();
    
    // We need a WFAExtender to do tail and intervening alignments.
    // Note that the extender expects anchoring matches!!!
    WFAExtender extender(gbwt_graph, aligner); 
    
    // Keep a couple cursors in the chain: extension before and after the linking up we need to do.
    auto here_it = chain.begin();
    auto next_it = here_it;
    ++next_it;
    
    const algorithms::Anchor* here = &to_chain[*here_it];
    
#ifdef debug_chaining
    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "First item " << *here_it
                << " with overall index " << to_chain.backing_index(*here_it)
                << " aligns " << (*here).read_start() << "-" << (*here).read_end()
                << " with " << (*here).graph_start() << "-" << (*here).graph_end()
                << endl;
        }
    }
#endif
    
    // We compose into a Path, since sometimes we may have to drop back to
    // aligners that aren't the WFAAligner and don't make WFAAlignments.
    Path composed_path;
    // We also track the total score of all the pieces.
    int composed_score = 0;
    
    // Do the left tail, if any.
    size_t left_tail_length = (*here).read_start();
    if (left_tail_length > 0) {
        // We need to do a left tail.
        // Anchor position will not be covered. 
        string left_tail = aln.sequence().substr(0, left_tail_length);
        WFAAlignment left_alignment;
        pos_t right_anchor = (*here).graph_start();
        if (left_tail.size() <= max_tail_length) {
            // Tail is short so keep to the GBWT.
            // We align the left tail with prefix(), which creates a prefix of the alignment.
            left_alignment = extender.prefix(left_tail, right_anchor);
            if (left_alignment && left_alignment.seq_offset != 0) {
                // We didn't get all the way to the left end of the read without
                // running out of score.
                // Prepend a softclip.
                // TODO: Can we let the aligner know it can softclip for free?
                WFAAlignment prepend = WFAAlignment::make_unlocalized_insertion(0, left_alignment.seq_offset, 0);
                prepend.join(left_alignment);
                left_alignment = std::move(prepend);
            }
            if (left_alignment.length != (*here).read_start()) {
                // We didn't get the alignment we expected.
                stringstream ss;
                ss << "Aligning left tail " << left_tail << " from " << (*here).graph_start() << " produced wrong-length alignment ";
                left_alignment.print(ss);
                throw ChainAlignmentFailedError(ss.str());
            }
        }
        if (left_alignment) {
            // We got an alignment, so make it a path
            left_alignment.check_lengths(gbwt_graph);
            
#ifdef debug_chaining
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Start with left tail of " << left_alignment.length << " with score of " << left_alignment.score << endl;
                }
            }
#endif
            
            composed_path = left_alignment.to_path(this->gbwt_graph, aln.sequence());
            composed_score = left_alignment.score;
        } else {
            // We need to fall back on alignment against the graph
            
            if (left_tail_length > MAX_DP_LENGTH) {
                // Left tail is too long to align.
                
                #pragma omp critical (cerr)
                {
                    cerr << "warning[MinimizerMapper::find_chain_alignment]: Refusing to align " << left_tail_length << " bp left tail against " << right_anchor << " in " << aln.name() << " to avoid overflow" << endl;
                }
                
                // Make a softclip for it.
                left_alignment = WFAAlignment::make_unlocalized_insertion(0, left_tail.size(), 0);
                composed_path = left_alignment.to_path(this->gbwt_graph, aln.sequence());
                composed_score = left_alignment.score;
            } else {
            
#ifdef debug_chaining
                if (show_work) {
                    #pragma omp critical (cerr)
                    {
                        cerr << log_name() << "Start with long left tail fallback alignment" << endl;
                    }
                }
#endif
                
                #pragma omp critical (cerr)
                {
                    cerr << "warning[MinimizerMapper::find_chain_alignment]: Falling back to non-GBWT alignment of " << left_tail_length << " bp left tail against " << right_anchor << " in " << aln.name() << endl;
                }
                
                Alignment tail_aln;
                tail_aln.set_sequence(left_tail);
                if (!aln.quality().empty()) {
                    tail_aln.set_quality(aln.quality().substr(0, left_tail_length));
                }
                
                // Work out how far the tail can see
                size_t graph_horizon = left_tail_length + this->get_regular_aligner()->longest_detectable_gap(aln, aln.sequence().begin());
                // Align the left tail, anchoring the right end.
                align_sequence_between(empty_pos_t(), right_anchor, graph_horizon, &this->gbwt_graph, this->get_regular_aligner(), tail_aln, this->max_dp_cells);
                // Since it's the left tail we can just clobber the path
                composed_path = tail_aln.path();
                composed_score = tail_aln.score();
            }
        }
    }
        
    size_t longest_attempted_connection = 0;
    while(next_it != chain.end()) {
        // Do each region between successive gapless extensions
        
        // We have to find the next item we can actually connect to
        const algorithms::Anchor* next;
        // And the actual connecting alignment to it
        WFAAlignment link_alignment;
        
        while (next_it != chain.end()) {
            next = &to_chain[*next_it];
            // Try and find a next thing to connect to
            
            if (algorithms::get_read_distance(*here, *next) == std::numeric_limits<size_t>::max()) {
                // There's overlap between these items. Keep here and skip next.
#ifdef debug_chaining
                if (show_work) {
                    #pragma omp critical (cerr)
                    {
                        cerr << log_name() << "Don't try and connect " << *here_it << " to " << *next_it << " because they overlap" << endl;
                    }
                }
#endif
            
                ++next_it;
            } else {
                // No overlap, so try it.
                break;
            }
        }
        
        if (next_it == chain.end()) {
            // We couldn't find anything to connect to
            break;
        }
            
#ifdef debug_chaining
        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "Add current item " << *here_it << " of length " << (*here).length() << " with score of " << (*here).score() << endl;
            }
        }
#endif
        
        // Make an alignment for the bases used in this item, and
        // concatenate it in.
        WFAAlignment here_alignment = this->to_wfa_alignment(*here);
        append_path(composed_path, here_alignment.to_path(this->gbwt_graph, aln.sequence()));
        composed_score += here_alignment.score;
        
#ifdef debug_chaining
        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "Next connectable item " << *next_it
                    << " with overall index " << to_chain.backing_index(*next_it)
                    << " aligns " << (*next).read_start() << "-" << (*next).read_end()
                    << " with " << (*next).graph_start() << "-" << (*next).graph_end()
                    << endl;
            }
        }
#endif
        
        // Pull out the intervening string to the next, if any.
        size_t link_start = (*here).read_end();
        size_t link_length = (*next).read_start() - link_start;
        string linking_bases = aln.sequence().substr(link_start, link_length);
        size_t graph_length = algorithms::get_graph_distance(*here, *next, *distance_index, gbwt_graph);
        
#ifdef debug_chaining
        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "Need to align graph from " << (*here).graph_end() << " to " << (*next).graph_start()
                    << " separated by " << graph_length << " bp and sequence \"" << linking_bases << "\"" << endl;
            }
        }
#endif
        
        if (link_length == 0 && graph_length == 0) {
            // These items abut in the read and the graph, so we assume we can just connect them.
            // WFAExtender::connect() can't handle an empty read sequence, and
            // our fallback method to align just against the graph can't handle
            // an empty graph region.
            // TODO: We can be leaving the GBWT's space here!
            
#ifdef debug_chaining
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Treat as empty link" << endl;
                }
            }
#endif
            
            link_alignment = WFAAlignment::make_empty();
        } else if (link_length > 0 && link_length <= max_chain_connection) {
            // If it's not empty and is a reasonable size, align it.
            // Make sure to walk back the left anchor so it is outside of the region to be aligned.
            pos_t left_anchor = (*here).graph_end();
            get_offset(left_anchor)--;
            
            link_alignment = extender.connect(linking_bases, left_anchor, (*next).graph_start());
            
            longest_attempted_connection = std::max(longest_attempted_connection, linking_bases.size());
            
            if (!link_alignment) {
                // We couldn't align.
                if (graph_length == 0) {
                    // We had read sequence but no graph sequence.
                    // Try falling back to a pure insertion.
                    // TODO: We can be leaving the GBWT's space here!
                    // TODO: What if this is forcing an insertion that could also be in the graph already?
#ifdef debug_chaining
                    if (show_work) {
                        #pragma omp critical (cerr)
                        {
                            cerr << log_name() << "connect() failed; treat as insertion" << endl;
                        }
                    }
#endif
                    link_alignment = WFAAlignment::make_unlocalized_insertion((*here).read_end(), link_length, aligner.score_gap(link_length));
                }
            } else if (link_alignment.length != linking_bases.size()) {
                // We could align, but we didn't get the alignment we expected. This shouldn't happen for a middle piece that can't softclip.
                stringstream ss;
                ss << "Aligning anchored link " << linking_bases << " (" << linking_bases.size() << " bp) from " << left_anchor << " - " << (*next).graph_start() << " against graph distance " << graph_length << " produced wrong-length alignment ";
                link_alignment.print(ss);
                throw ChainAlignmentFailedError(ss.str());
            } else {
                // We got the right alignment.
                // Put the alignment back into full read space
                link_alignment.seq_offset += (*here).read_end();
            }
        }
        
        if (link_alignment) {
            // We found a link alignment
            
#ifdef debug_chaining
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Add link of length " << link_alignment.length << " with score of " << link_alignment.score << endl;
                }
            }
#endif
        
            link_alignment.check_lengths(gbwt_graph);
            
            // Then the link (possibly empty)
            append_path(composed_path, link_alignment.to_path(this->gbwt_graph, aln.sequence()));
            composed_score += link_alignment.score;
        } else {
            // The sequence to the next thing is too long, or we couldn't reach it doing connect().
            // Fall back to another alignment method
            
            if (linking_bases.size() > MAX_DP_LENGTH) {
                // This would be too long for GSSW to handle and might overflow 16-bit scores in its matrix.
                #pragma omp critical (cerr)
                {
                    cerr << "warning[MinimizerMapper::find_chain_alignment]: Refusing to align " << link_length << " bp connection between chain items " << to_chain.backing_index(*here_it) << " and " << to_chain.backing_index(*next_it) << " which are " << graph_length << " apart at " << (*here).graph_end() << " and " << (*next).graph_start() << " in " << aln.name() << " to avoid overflow" << endl;
                }
                // Just jump to right tail
                break;
            }
            
            // We can't actually do this alignment, we'd have to align too
            // long of a sequence to find a connecting path.
            #pragma omp critical (cerr)
            {
                cerr << "warning[MinimizerMapper::find_chain_alignment]: Falling back to non-GBWT alignment of " << link_length << " bp connection between chain items " << to_chain.backing_index(*here_it) << " and " << to_chain.backing_index(*next_it) << " which are " << graph_length << " apart at " << (*here).graph_end() << " and " << (*next).graph_start() << " in " << aln.name() << endl;
            }
            
            Alignment link_aln;
            link_aln.set_sequence(linking_bases);
            if (!aln.quality().empty()) {
                link_aln.set_quality(aln.quality().substr(link_start, link_length));
            }
            // Guess how long of a graph path we ought to allow in the alignment.
            size_t path_length = std::max(graph_length, link_length) + this->get_regular_aligner()->longest_detectable_gap(aln, aln.sequence().begin() + link_start);
            MinimizerMapper::align_sequence_between((*here).graph_end(), (*next).graph_start(), path_length, &this->gbwt_graph, this->get_regular_aligner(), link_aln, this->max_dp_cells);
            
#ifdef debug_chaining
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Add link of length " << path_to_length(link_aln.path()) << " with score of " << link_aln.score() << endl;
                }
            }
#endif
            
            // Then tack that path and score on
            append_path(composed_path, link_aln.path());
            composed_score += link_aln.score();
        }
        
        // Advance here to next and start considering the next after it
        here_it = next_it;
        ++next_it;
        here = next;
    }
    
#ifdef debug_chaining
    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "Add last extension " << *here_it << " of length " << (*here).length() << " with score of " << (*here).score() << endl;
        }
    }
#endif
    
    WFAAlignment here_alignment = this->to_wfa_alignment(*here);
    
    here_alignment.check_lengths(gbwt_graph);
    
    // Do the final GaplessExtension itself (may be the first)
    append_path(composed_path, here_alignment.to_path(this->gbwt_graph, aln.sequence()));
    composed_score += here_alignment.score;
   
    // Do the right tail, if any. Do as much of it as we can afford to do.
    size_t right_tail_length = aln.sequence().size() - (*here).read_end();
    if (right_tail_length > 0) {
        // We need to do a right tail
        string right_tail = aln.sequence().substr((*here).read_end(), right_tail_length);
        WFAAlignment right_alignment;
        pos_t left_anchor = (*here).graph_end();
        get_offset(left_anchor)--;
        if (right_tail_length <= max_tail_length) {
            // We align the right tail with suffix(), which creates a suffix of the alignment.
            // Make sure to walk back the anchor so it is outside of the region to be aligned.
            right_alignment = extender.suffix(right_tail, left_anchor);
        }
        
        if (right_alignment) {
            // Right tail did align. Put the alignment back into full read space.
            right_alignment.seq_offset += (*here).read_end();
            if (right_alignment.seq_offset + right_alignment.length != aln.sequence().size()) {
                // We didn't get all the way to the right end of the read without
                // running out of score.
                // Append a softclip.
                // TODO: Can we let the aligner know it can softclip for free?
                size_t right_end = right_alignment.seq_offset + right_alignment.length;
                size_t remaining = aln.sequence().size() - right_end;
                right_alignment.join(WFAAlignment::make_unlocalized_insertion(right_end, remaining, 0));
            }
            if (right_alignment.length != right_tail_length) {
                // We didn't get the alignment we expected.
                stringstream ss;
                ss << "Aligning right tail " << right_tail << " from " << left_anchor << " produced wrong-length alignment ";
                right_alignment.print(ss);
                throw ChainAlignmentFailedError(ss.str());
            }
#ifdef debug_chaining
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Add right tail of " << right_tail.size() << " with score of " << right_alignment.score << endl;
                }
            }
#endif
            
            right_alignment.check_lengths(gbwt_graph);
            
            append_path(composed_path, right_alignment.to_path(this->gbwt_graph, aln.sequence()));
            composed_score += right_alignment.score;
        } else {
            // We need to fall back on alignment against the graph
            
#ifdef debug_chaining
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "End with long right tail fallback alignment" << endl;
                }
            }
#endif

            if (right_tail.size() > MAX_DP_LENGTH) {
                // Right tail is too long to align.
                
                #pragma omp critical (cerr)
                {
                    cerr << "warning[MinimizerMapper::find_chain_alignment]: Refusing to align " << right_tail.size() << " bp right tail against " << left_anchor << " in " << aln.name() << " to avoid overflow" << endl;
                }
                
                // Make a softclip for it.
                right_alignment = WFAAlignment::make_unlocalized_insertion((*here).read_end(), aln.sequence().size() - (*here).read_end(), 0);
                append_path(composed_path, right_alignment.to_path(this->gbwt_graph, aln.sequence()));
                composed_score += right_alignment.score;
            } else {

                #pragma omp critical (cerr)
                {
                    cerr << "warning[MinimizerMapper::find_chain_alignment]: Falling back to non-GBWT alignment of " << right_tail_length << " bp right tail against " << left_anchor << " in " << aln.name() << endl;
                }
                
                Alignment tail_aln;
                tail_aln.set_sequence(right_tail);
                if (!aln.quality().empty()) {
                    tail_aln.set_quality(aln.quality().substr((*here).read_end(), right_tail_length));
                }
                
                // Work out how far the tail can see
                size_t graph_horizon = right_tail_length + this->get_regular_aligner()->longest_detectable_gap(aln, aln.sequence().begin() + (*here).read_end());
                // Align the right tail, anchoring the left end.
                align_sequence_between(left_anchor, empty_pos_t(), graph_horizon, &this->gbwt_graph, this->get_regular_aligner(), tail_aln, this->max_dp_cells);
                // Since it's the right tail we have to add it on
                append_path(composed_path, tail_aln.path());
                composed_score += tail_aln.score();
            }
        } 
    }
    
    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "Composed alignment is length " << path_to_length(composed_path) << " with score of " << composed_score << endl;
            if (composed_path.mapping_size() > 0) {
                cerr << log_name() << "Composed alignment starts with: " << pb2json(composed_path.mapping(0)) << endl;
                cerr << log_name() << "Composed alignment ends with: " << pb2json(composed_path.mapping(composed_path.mapping_size() - 1)) << endl;
            }
        }
    }
    
    // Convert to a vg Alignment.
    Alignment result(aln);
    *result.mutable_path() = std::move(simplify(composed_path));
    result.set_score(composed_score);
    if (!result.sequence().empty()) {
        result.set_identity(identity(result.path()));
    }
    
    set_annotation(result, "left_tail_length", (double) left_tail_length);
    set_annotation(result, "longest_attempted_connection", (double) longest_attempted_connection); 
    set_annotation(result, "right_tail_length", (double) right_tail_length); 
    
    return result;
}
#undef debug_chaining

void MinimizerMapper::wfa_alignment_to_alignment(const WFAAlignment& wfa_alignment, Alignment& alignment) const {
    *(alignment.mutable_path()) = wfa_alignment.to_path(this->gbwt_graph, alignment.sequence());
    alignment.set_score(wfa_alignment.score);
    if (!alignment.sequence().empty()) {
        alignment.set_identity(identity(alignment.path()));
    }
}

void MinimizerMapper::with_dagified_local_graph(const pos_t& left_anchor, const pos_t& right_anchor, size_t max_path_length, const HandleGraph& graph, const std::function<void(DeletableHandleGraph&, const std::function<std::pair<nid_t, bool>(const handle_t&)>&)>& callback) {
    
    if (is_empty(left_anchor) && is_empty(right_anchor)) {
        throw ChainAlignmentFailedError("Cannot align sequence between two unset positions");
    }
    
    // We need to get the graph to align to.
    bdsg::HashGraph local_graph;
    unordered_map<id_t, id_t> local_to_base;
    if (!is_empty(left_anchor) && !is_empty(right_anchor)) {
        // We want a graph actually between two positions
        local_to_base = algorithms::extract_connecting_graph(
            &graph,
            &local_graph,
            max_path_length,
            left_anchor, right_anchor
        );
    } else if (!is_empty(left_anchor)) {
        // We only have the left anchor
        local_to_base = algorithms::extract_extending_graph(
            &graph,
            &local_graph,
            max_path_length,
            left_anchor,
            false,
            false
        );
    } else {
        // We only have the right anchor
        local_to_base = algorithms::extract_extending_graph(
            &graph,
            &local_graph,
            max_path_length,
            right_anchor,
            true,
            false
        );
    }
    
    // To find the anchoring nodes in the extracted graph, we need to scan local_to_base.
    nid_t local_left_anchor_id = 0;
    nid_t local_right_anchor_id = 0;
    for (auto& kv : local_to_base) {
        if (kv.second == id(left_anchor)) {
            local_left_anchor_id = kv.first;
        }
        if (kv.second == id(right_anchor)) {
            local_right_anchor_id = kv.first;
        }
        // TODO: Stop early when we found them all.
    }

    if (!is_empty(left_anchor) && local_left_anchor_id == 0) {
        // Somehow the left anchor didn't come through. Complain.
        std::stringstream ss;
        ss << "Extracted graph from " << left_anchor;
        if (!is_empty(right_anchor)) {
            ss << " to " << right_anchor;
        }
        ss << " with max path length of " << max_path_length;
        ss << " but from node was not present in the resulting translation";
        local_graph.serialize("crashdump.vg");
        throw ChainAlignmentFailedError(ss.str());
    }

    if (!is_empty(right_anchor) && local_right_anchor_id == 0) {
        // Somehow the right anchor didn't come through. Complain.
        std::stringstream ss;
        ss << "Extracted graph";
        if (!is_empty(left_anchor)) {
            ss << " from " << left_anchor;
        }
        ss << " to " << right_anchor;
        ss << " with max path length of " << max_path_length;
        ss << " but to node was not present in the resulting translation";
        local_graph.serialize("crashdump.vg");
        throw ChainAlignmentFailedError(ss.str());
    }
    
    // And split by strand since we can only align to one strand
    StrandSplitGraph split_graph(&local_graph);
    
    // And make sure it's a DAG of the stuff reachable from our anchors
    bdsg::HashGraph dagified_graph;
    // For which we need the handles that anchor the graph, facing inwards
    std::vector<handle_t> bounding_handles;
    if (!is_empty(left_anchor)) {
        // Dagify from the forward version of the left anchor
        
        // Grab the left anchor in the local graph
        if (!local_graph.has_node(local_left_anchor_id)) {
            std::stringstream ss;
            ss << "Extracted graph from " << left_anchor;
            if (!is_empty(right_anchor)) {
                ss << " to " << right_anchor;
            }
            ss << " with max path length of " << max_path_length;
            ss << " but from node local ID " << local_left_anchor_id << " was not present in the resulting graph";
            local_graph.serialize("crashdump.vg");
            throw ChainAlignmentFailedError(ss.str());
        }
        handle_t local_handle = local_graph.get_handle(local_left_anchor_id, is_rev(left_anchor));
        
        // And get the node that that orientation of it is in the strand-split graph
        handle_t overlay_handle = split_graph.get_overlay_handle(local_handle);
        
        // And use that
        bounding_handles.push_back(overlay_handle);
    }
    if (!is_empty(right_anchor)) {
        // Dagify from the reverse version of the node for the forward version of the right anchor
        
        // Grab the right anchor from the local graph
        if (!local_graph.has_node(local_right_anchor_id)) {
            std::stringstream ss;
            ss << "Extracted graph";
            if (!is_empty(left_anchor)) {
                ss << " from " << left_anchor;
            }
            ss << " to " << right_anchor;
            ss << " with max path length of " << max_path_length;
            ss << " but to node local ID " << local_right_anchor_id << " was not present in the resulting graph";
            local_graph.serialize("crashdump.vg");
            throw ChainAlignmentFailedError(ss.str());
        }
        handle_t local_handle = local_graph.get_handle(local_right_anchor_id, is_rev(right_anchor));
        
        // And get the node that that orientation of it is in the strand-split graph
        // But flip it because we want to dagify going inwards from the right
        handle_t overlay_handle = split_graph.flip(split_graph.get_overlay_handle(local_handle));
        
        // And use that
        bounding_handles.push_back(overlay_handle);
    }
    
    auto dagified_to_split = handlegraph::algorithms::dagify_from(&split_graph, bounding_handles, &dagified_graph, max_path_length);
    
#ifdef debug
    std::cerr << "Dagified from " << bounding_handles.size() << " bounding handles in " << split_graph.get_node_count() << " node strand-split graph to " << dagified_graph.get_node_count() << " node DAG" << std::endl;
#endif

    // Make an accessor for getting back to the base graph space
    auto dagified_handle_to_base = [&](const handle_t& h) -> pair<nid_t, bool> {
        nid_t dagified_id = dagified_graph.get_id(h);
        bool dagified_is_reverse = dagified_graph.get_is_reverse(h);
        auto found_in_split = dagified_to_split.find(dagified_id);
        if (found_in_split == dagified_to_split.end()) {
            throw ChainAlignmentFailedError("ID " + std::to_string(dagified_id) + " from dagified graph not found in strand-split graph");
        }
        nid_t split_id = found_in_split->second;
        handle_t split_handle = split_graph.get_handle(split_id, dagified_is_reverse);
        // We rely on get_underlying_handle understanding reversed handles in the split graph
        handle_t local_handle = split_graph.get_underlying_handle(split_handle);
        nid_t local_id = local_graph.get_id(local_handle);
        bool local_is_reverse = local_graph.get_is_reverse(local_handle);
        auto found_in_base = local_to_base.find(local_id);
        if (found_in_base == local_to_base.end()) {
            throw ChainAlignmentFailedError("ID " + std::to_string(local_id) + " from local graph not found in full base graph");
        }
        nid_t base_id = found_in_base->second;
        return std::make_pair(base_id, local_is_reverse);
    };
    
    // Show the graph we made and the translation function
    callback(dagified_graph, dagified_handle_to_base);
}

void MinimizerMapper::align_sequence_between(const pos_t& left_anchor, const pos_t& right_anchor, size_t max_path_length, const HandleGraph* graph, const GSSWAligner* aligner, Alignment& alignment, size_t max_dp_cells) {
    
    // Get the dagified local graph, and the back translation
    MinimizerMapper::with_dagified_local_graph(left_anchor, right_anchor, max_path_length, *graph,
        [&](DeletableHandleGraph& dagified_graph, const std::function<std::pair<nid_t, bool>(const handle_t&)>& dagified_handle_to_base) {
    
        // Then trim off the tips that are either in the wrong orientation relative
        // to whether we want them to be a source or a sink, or extraneous
        
        std::vector<handle_t> tip_handles = handlegraph::algorithms::find_tips(&dagified_graph);
        bool trimmed;
        size_t trim_count = 0;
        do {
            trimmed = false;
            // We need to make sure to remove only one orientation of each handle
            // we remove.
            std::unordered_set<nid_t> to_remove_ids;
            std::vector<handle_t> to_remove_handles;
            for (auto& h : tip_handles) {
                auto base_coords = dagified_handle_to_base(h);
                if (!dagified_graph.get_is_reverse(h) && (is_empty(left_anchor) || base_coords.first == id(left_anchor))) {
                    // Tip is inward forward, so it's a source.
                    // This is a head in the subgraph, and either matches a left
                    // anchoring node or we don't have any, so keep it.
#ifdef debug
                    std::cerr << "Dagified graph node " << dagified_graph.get_id(h) << " " << dagified_graph.get_is_reverse(h) << " is an acceptable source (" << base_coords.first << " " << base_coords.second << ")" << std::endl;
#endif
                } else if (dagified_graph.get_is_reverse(h) && (is_empty(right_anchor) || base_coords.first == id(right_anchor))) {
                    // Tip is inward reverse, so it's a sink.
                    // This is a tail in the subgraph, and either matches a right
                    // anchoring node or we don't have any, so keep it.
#ifdef debug
                    std::cerr << "Dagified graph node " << dagified_graph.get_id(h) << " " << dagified_graph.get_is_reverse(h) << " is an acceptable sink (" << base_coords.first << " " << base_coords.second << ")" << std::endl;
#endif
                } else {
                    // This is a wrong orientation of an anchoring node, or some other tip.
                    // We don't want to keep this handle
#ifdef debug
                    std::cerr << "Dagified graph node " << dagified_graph.get_id(h) << " " << dagified_graph.get_is_reverse(h) << " is an unacceptable tip (" << base_coords.first << " " << base_coords.second << ")" << std::endl;
#endif
                    nid_t dagified_id = dagified_graph.get_id(h);
                    if (!to_remove_ids.count(dagified_id)) {
                        to_remove_ids.insert(dagified_id);
                        to_remove_handles.push_back(h);
                    }
                }
            }
            for (auto& h : to_remove_handles) {
                dagified_graph.destroy_handle(h);
                trimmed = true;
            }
            if (trimmed) {
                // TODO: This is going to be O(slow) if we actually have to
                // prune back a dangling run. We should look at what is
                // connected to the tip and the tip only, and make that the new
                // tip. Or keep some kind of online tip info. Or use an
                // algorithm function that we make actually good.
                tip_handles = handlegraph::algorithms::find_tips(&dagified_graph);
                trim_count++;
            }
        } while (trimmed);
        if (trim_count > 0) {
            #pragma omp critical (cerr)
            std::cerr << "warning[MinimizerMapper::align_sequence_between]: Trimmed back tips " << trim_count << " times on graph between " << left_anchor << " and " << right_anchor << " leaving " <<  dagified_graph.get_node_count() << " nodes and " << tip_handles.size() << " tips" << std::endl;
        }
        
        if (!is_empty(left_anchor) && !is_empty(right_anchor)) {
            // Then align the linking bases, with global alignment so they have
            // to go from a source to a sink. Banded alignment means we can safely do big problems.
            aligner->align_global_banded(alignment, dagified_graph);
        } else {
            // Do pinned alignment off the anchor we actually have.
            // Don't use X-Drop because Dozeu is known to just overwrite the
            // stack with garbage whenever alignments are "too big", and these
            // alignments are probably often too big.
            // But if we don't use Dozeu this uses GSSW and that can *also* be too big.
            // So work out how big it will be
            size_t cell_count = dagified_graph.get_total_length() * alignment.sequence().size();
            if (cell_count > max_dp_cells) {
                #pragma omp critical (cerr)
                std::cerr << "warning[MinimizerMapper::align_sequence_between]: Refusing to fill " << cell_count << " DP cells in tail with GSSW" << std::endl;
                // Fake a softclip right in input graph space
                alignment.clear_path();
                Mapping* m = alignment.mutable_path()->add_mapping();
                // TODO: Is this fake position OK regardless of anchoring side?
                m->mutable_position()->set_node_id(is_empty(left_anchor) ? id(right_anchor) : id(left_anchor));
                m->mutable_position()->set_is_reverse(is_empty(left_anchor) ? is_rev(right_anchor) : is_rev(left_anchor));
                m->mutable_position()->set_offset(is_empty(left_anchor) ? offset(right_anchor) : offset(left_anchor));
                Edit* e = m->add_edit();
                e->set_to_length(alignment.sequence().size());
                e->set_sequence(alignment.sequence());
                return;
            } else {
#ifdef debug_chaining
                #pragma omp critical (cerr)
                std::cerr << "debug[MinimizerMapper::align_sequence_between]: Fill " << cell_count << " DP cells in tail with GSSW" << std::endl;
#endif
                aligner->align_pinned(alignment, dagified_graph, !is_empty(left_anchor), false);
            }
        }
        
        // And translate back into original graph space
        for (size_t i = 0; i < alignment.path().mapping_size(); i++) {
            // Translate each mapping's ID and orientation down to the base graph
            Mapping* m = alignment.mutable_path()->mutable_mapping(i);
            
            handle_t dagified_handle = dagified_graph.get_handle(m->position().node_id(), m->position().is_reverse());
            auto base_coords = dagified_handle_to_base(dagified_handle); 
            
            m->mutable_position()->set_node_id(base_coords.first);
            m->mutable_position()->set_is_reverse(base_coords.second);
        }
        if (!is_empty(left_anchor) && alignment.path().mapping_size() > 0 && offset(left_anchor) != 0) {
            // Get the positions of the leftmost mapping
            Position* left_pos = alignment.mutable_path()->mutable_mapping(0)->mutable_position();
            // Add on the offset for the missing piece of the left anchor node
            left_pos->set_offset(left_pos->offset() + offset(left_anchor));
        }
        if (alignment.path().mapping_size() > 0) {
            // Make sure we don't have an empty mapping on the end
            auto* last_mapping = alignment.mutable_path()->mutable_mapping(alignment.path().mapping_size() - 1);
            if (last_mapping->edit_size() > 0) {
                // Make sure we don't have an empty edit on the end
                auto& last_edit = last_mapping->edit(last_mapping->edit_size() - 1);
                if (last_edit.from_length() == 0 && last_edit.to_length() == 0 && last_edit.sequence().empty()) {
                    // Last edit is empty so drop from the mapping
                    last_mapping->mutable_edit()->RemoveLast();
                }
            }
            if (last_mapping->edit_size() == 0) {
                // Last mapping is empty, so drop it.
                alignment.mutable_path()->mutable_mapping()->RemoveLast();
            }
        }
    
        // Now the alignment is filled in!
    });
}

std::vector<algorithms::Anchor> MinimizerMapper::to_anchors(const Alignment& aln, const VectorView<Minimizer>& minimizers, const std::vector<Seed>& seeds) const {
    std::vector<algorithms::Anchor> to_return;
    to_return.reserve(seeds.size());
    for (size_t i = 0; i < seeds.size(); i++) {
        to_return.push_back(this->to_anchor(aln, minimizers, seeds, i));
    }
    return to_return;
}

algorithms::Anchor MinimizerMapper::to_anchor(const Alignment& aln, const VectorView<Minimizer>& minimizers, const std::vector<Seed>& seeds, size_t seed_number) const {
    // Turn each seed into the part of its match on the node where the
    // anchoring end (start for forward-strand minimizers, end for
    // reverse-strand minimizers) falls.
    auto& seed = seeds[seed_number];
    auto& source = minimizers[seed.source];
    size_t length;
    pos_t graph_start;
    size_t read_start;
    size_t hint_start;
    if (source.value.is_reverse) {
        // Seed stores the final base of the match in the graph.
        // So get the past-end position.
        pos_t graph_end = make_pos_t(id(seed.pos), is_rev(seed.pos), offset(seed.pos) + 1);
        
        // Work out how much of the node it could use before there.
        length = std::min((size_t) source.length, offset(graph_end));
        // And derive the graph start
        graph_start = make_pos_t(id(graph_end), is_rev(graph_end), offset(graph_end) - length);
        // And the read start
        read_start = source.value.offset + 1 - length;
        // The seed is actually the last 1bp interval
        hint_start = length - 1;
    } else {
        // Seed stores the first base of the match in the graph
        graph_start = seed.pos;
        
        // Get the handle to the node it's on.
        handle_t start_handle = gbwt_graph.get_handle(id(graph_start), is_rev(graph_start));
        // Work out how much of the node it could use before there.
        length = std::min((size_t) source.length, gbwt_graph.get_length(start_handle) - offset(graph_start));
        
        // And we store the read start position already in the item
        read_start = source.value.offset;
        // The seed is actually at the start
        hint_start = 0;
    }
    // Work out how many points the anchor is
    // TODO: Always make sequence and quality available for scoring!
    int score = get_regular_aligner()->score_exact_match(aln, read_start, length);
    return algorithms::Anchor(read_start, graph_start, length, score, seed_number, seed.zipcode_decoder.get(), hint_start); 
}

WFAAlignment MinimizerMapper::to_wfa_alignment(const algorithms::Anchor& anchor) const {
    return {
        {gbwt_graph.get_handle(id(anchor.graph_start()), is_rev(anchor.graph_start()))},
        {{WFAAlignment::match, (uint32_t)anchor.length()}},
        (uint32_t)offset(anchor.graph_start()),
        (uint32_t)anchor.read_start(),
        (uint32_t)anchor.length(),
        anchor.score(),
        true
    };
}

}
