/**
 * \file minimizer_mapper_from_chains.cpp
 * Defines the code for the long-read code path for the
 * minimizer-and-GBWT-based mapper (long read Giraffe).
 */

#include "minimizer_mapper.hpp"

#include "annotation.hpp"
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

void MinimizerMapper::score_merged_cluster(Cluster& cluster, 
                                           size_t i,
                                           const VectorView<Minimizer>& minimizers,
                                           const std::vector<Seed>& seeds,
                                           size_t first_new_seed,
                                           const std::vector<size_t>& seed_to_bucket,
                                           const std::vector<Cluster>& buckets,
                                           size_t seq_length,
                                           Funnel& funnel) const {
    

    if (this->track_provenance) {
        // Say we're making it
        funnel.producing_output(i);
    }

    // Initialize the values.
    cluster.score = 0.0;
    cluster.coverage = 0.0;
    cluster.present = SmallBitset(minimizers.size()); // TODO: This is probably usually too big to really be "small" now.
    
    // Collect the old clusters and new seeds we are coming from
    // TODO: Skip if not tracking provenance?
    std::vector<size_t> to_combine;
    // Deduplicate old clusters with a bit set
    SmallBitset buckets_seen(buckets.size());
    

    // Determine the minimizers that are present in the cluster.
    for (auto hit_index : cluster.seeds) {
        // We have this seed's minimizer
        cluster.present.insert(seeds[hit_index].source);
        
        if (hit_index < first_new_seed) {
            // An old seed.
            // We can also pick up an old cluster.
            size_t old_cluster = seed_to_bucket.at(hit_index);
            if (old_cluster != std::numeric_limits<size_t>::max()) {
                // This seed came form an old cluster, so we must have eaten it
                if (!buckets_seen.contains(old_cluster)) {
                    // Remember we used this old cluster
                    to_combine.push_back(old_cluster);
                    buckets_seen.insert(old_cluster);
                }
            }
        } else {
            // Make sure we tell the funnel we took in this new seed.
            // Translate from a space that is old seeds and then new seeds to a
            // space that is old *clusters* and then new seeds
            to_combine.push_back(hit_index - first_new_seed + buckets.size());
        }
    }
    if (show_work) {
        #pragma omp critical (cerr)
        dump_debug_clustering(cluster, i, minimizers, seeds);
    }

    // Compute the score and cluster coverage.
    sdsl::bit_vector covered(seq_length, 0);
    for (size_t j = 0; j < minimizers.size(); j++) {
        if (cluster.present.contains(j)) {
            const Minimizer& minimizer = minimizers[j];
            cluster.score += minimizer.score;

            // The offset of a reverse minimizer is the endpoint of the kmer
            size_t start_offset = minimizer.forward_offset();
            size_t k = minimizer.length;

            // Set the k bits starting at start_offset.
            covered.set_int(start_offset, sdsl::bits::lo_set[k], k);
        }
    }
    // Count up the covered positions and turn it into a fraction.
    cluster.coverage = sdsl::util::cnt_one_bits(covered) / static_cast<double>(seq_length);

    if (this->track_provenance) {
        // Record the cluster in the funnel as a group combining the previous groups.
        funnel.merge_groups(to_combine.begin(), to_combine.end());
        funnel.score(funnel.latest(), cluster.score);

        // Say we made it.
        funnel.produced_output();
    }

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

MinimizerMapper::chain_set_t MinimizerMapper::chain_clusters(const Alignment& aln, const VectorView<Minimizer>& minimizers, const std::vector<Seed>& seeds, const std::vector<Cluster>& clusters, const chain_config_t& cfg, size_t old_seed_count, size_t new_seed_start, Funnel& funnel, size_t seed_stage_offset, size_t reseed_stage_offset, LazyRNG& rng) const {

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
                    cerr << log_name() << "Computing chain over " << cluster_seeds_sorted.size() << " seeds" << endl;
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
            VectorView<algorithms::Anchor> cluster_view {seed_anchors, cluster_seeds_sorted};
            std::vector<std::pair<int, std::vector<size_t>>> chains = algorithms::find_best_chains(
                cluster_view,
                *distance_index,
                gbwt_graph,
                get_regular_aligner()->gap_open,
                get_regular_aligner()->gap_extension,
                cfg.max_chains_per_cluster,
                cfg.max_lookback_bases,
                cfg.min_lookback_items,
                cfg.lookback_item_hard_cap,
                cfg.initial_lookback_threshold,
                cfg.lookback_scale_factor,
                cfg.min_good_transition_score_per_base,
                cfg.item_bonus,
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
    // We may or may not need to invert this view, but if we do we will want to
    // keep the result. So have a place to lazily keep an inverse.
    std::unique_ptr<VectorViewInverse> minimizer_score_sort_inverse;
    
    // Find the seeds and mark the minimizers that were located.
    vector<Seed> seeds = this->find_seeds(minimizers, aln, funnel);
    
    // Pre-cluster just the seeds we have. Get sets of input seed indexes that go together.
    if (track_provenance) {
        funnel.stage("bucket");
        funnel.substage("compute-buckets");
    }

    // Bucket the hits coarsely into sets that might be able to interact.
    std::vector<Cluster> buckets = clusterer.cluster_seeds(seeds, aln.sequence().size() * bucket_scale);
    
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
    fragment_cfg.max_lookback_bases = 200;
    fragment_cfg.min_lookback_items = 0;
    fragment_cfg.lookback_item_hard_cap = 3;
    fragment_cfg.initial_lookback_threshold = this->initial_lookback_threshold;
    fragment_cfg.lookback_scale_factor = this->lookback_scale_factor;
    fragment_cfg.min_good_transition_score_per_base = this->min_good_transition_score_per_base;
    
    fragment_cfg.item_bonus = this->item_bonus;
    fragment_cfg.max_indel_bases = 50;
    
    // Do all the ones that are 75% as good as the best, or down to 50% as good
    // as the best if that is what it takes to get the second best
    double bucket_score_cutoff = best_bucket_score / 0.75;
    if (bucket_score_cutoff - (bucket_score_cutoff / 0.25) < second_best_bucket_score) {
        bucket_score_cutoff = std::min(bucket_score_cutoff, second_best_bucket_score);
    }
    fragment_cfg.cluster_score_cutoff = bucket_score_cutoff;
    fragment_cfg.cluster_score_cutoff_enabled = true;
    fragment_cfg.cluster_coverage_threshold = 1.0;
    fragment_cfg.min_clusters_to_chain = std::numeric_limits<size_t>::max();
    fragment_cfg.max_clusters_to_chain = std::numeric_limits<size_t>::max();
    
    fragment_cfg.max_chains_per_cluster = this->max_fragments_per_bucket;
    
    auto fragment_results = this->chain_clusters(aln, minimizers, seeds, buckets, fragment_cfg, seeds.size(), seeds.size(), funnel, 2, std::numeric_limits<size_t>::max(), rng);
    
    if (track_provenance) {
        funnel.substage("translate-fragments");
    }
    
    // Translate fragment chains into faked clusters, which downstream code expects. They need a seeds[] and a coverage.
    std::vector<Cluster> fragments;
    for (size_t i = 0; i < fragment_results.cluster_chains.size(); i++) {
        // For each source bucket
        for (auto& chain : fragment_results.cluster_chains[i]) {
            // For each fragment found in the bucket
        
            // Convert format
            fragments.emplace_back();
        
            if (this->track_provenance) {
                // Say we're making it
                funnel.producing_output(fragments.size());
            } 
            // Copy all the seeds in the chain over
            fragments.back().seeds.reserve(chain.second.size());
            for (auto& chain_visited_index : chain.second) {
                // Make sure to translate to real seed space
                fragments.back().seeds.push_back(fragment_results.cluster_chain_seeds[i].at(chain_visited_index));
            }
            // Rescore as a cluster
            this->score_cluster(fragments.back(), fragments.size() - 1, minimizers, seeds, aln.sequence().size());
            if (this->track_provenance) {
                // Record the fragment in the funnel as coming from the bucket 
                funnel.project(i);
                funnel.score(funnel.latest(), fragments.back().score);

                // Say we made it.
                funnel.produced_output();
            }
        }
    }
    
    // Find pairs of "adjacent" fragments
    if (track_provenance) {
        funnel.stage("reseed");
        funnel.substage("pair-fragments");
    }
    
    // To do that, we need start end end positions for each fragment, in the read
    std::vector<std::pair<size_t, size_t>> fragment_read_ranges(fragments.size(), {std::numeric_limits<size_t>::max(), 0});
    // And the lowest-numbered seeds in the fragment from those minimizers.
    std::vector<std::pair<size_t, size_t>> fragment_bounding_seeds(fragments.size(), {std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()});
    for (size_t i = 0; i < fragments.size(); i++) {
        // For each fragment
        auto& fragment = fragments[i];
        // We will fill in the range it occupies in the read
        auto& read_range = fragment_read_ranges[i];
        auto& graph_seeds = fragment_bounding_seeds[i];
        for (auto& seed_index : fragment.seeds) {
            // Which means we look at the minimizer for each seed
            auto& minimizer = minimizers[seeds[seed_index].source];
            
            if (minimizer.forward_offset() < read_range.first) {
                // Min all their starts to get the fragment start
                read_range.first = minimizer.forward_offset();
                if (seed_index < graph_seeds.first) {
                    // And keep a seed hit
                    graph_seeds.first = seed_index;
                }
            }
            
            if (minimizer.forward_offset() + minimizer.length > read_range.second) {
                // Max all their past-ends to get the fragment past-end
                read_range.second = minimizer.forward_offset() + minimizer.length;
                if (seed_index < graph_seeds.second) {
                    // And keep a seed hit
                    graph_seeds.second = seed_index;
                }
            }
        }
    }
    
    // Record fragment statistics
    // Chaining score (and implicitly fragment count)
    std::vector<double> fragment_scores;
    // Chain length
    std::vector<double> fragment_item_counts;
    // Best fragment score in each bucket
    std::vector<double> bucket_best_fragment_scores;
    // Score of each bucket
    std::vector<double> bucket_scores;
    // Coverage of each bucket
    std::vector<double> bucket_coverages;
    for (size_t bucket_num = 0; bucket_num < fragment_results.cluster_chains.size(); bucket_num++) {
        auto& bucket = fragment_results.cluster_chains[bucket_num];
        double best_fragment_score = 0;
        for (auto& fragment : bucket) {
            fragment_scores.push_back(fragment.first);
            fragment_item_counts.push_back(fragment.second.size());
            best_fragment_score = std::max(best_fragment_score, (double) fragment.first);
        }
        bucket_best_fragment_scores.push_back(best_fragment_score);
    }
    // Bucket with the best fragment score
    size_t best_bucket = 0;
    // That score
    double best_bucket_fragment_score = 0;
    for (size_t i = 0; i < fragment_scores.size(); i++) {
        if (fragment_scores[i] >= best_bucket_fragment_score) {
            best_bucket_fragment_score = fragment_scores[i];
            best_bucket = fragment_results.cluster_nums[i];
        }
    }
    for (auto& bucket_num : fragment_results.cluster_nums) {
        // Record the info about the buckets that the fragments came from
        bucket_scores.push_back(buckets.at(bucket_num).score);
        bucket_coverages.push_back(buckets.at(bucket_num).coverage);
    }
    
    // Coverage of read by each fragment, using outer bounds 
    std::vector<double> fragment_bound_coverages;
    for (size_t i = 0; i < fragments.size(); i++) {
        auto& fragment = fragments[i];
        fragment_bound_coverages.push_back((double) (fragment_read_ranges[i].second - fragment_read_ranges[i].first) / aln.sequence().size());
    }
    // Overall coverage of read with fragments of item count k or greater, in best bucket
    // Remember: best bucket was the one that had the fragment with the best score.
    std::vector<double> best_bucket_fragment_coverage_at_length(21, 0.0);
    std::vector<bool> fragment_covered(aln.sequence().size(), false);
    for (int threshold = best_bucket_fragment_coverage_at_length.size() - 1; threshold >= 0; threshold--) {
        for (size_t i = 0; i < fragments.size(); i++) {
            if (fragment_results.cluster_nums[i] != best_bucket) {
                // Only look at the best bucket's fragments here.
                continue;
            }
            if (threshold == (best_bucket_fragment_coverage_at_length.size() - 1) && fragments[i].seeds.size() > threshold || fragments[i].seeds.size() == threshold) {
                // Need to mark this fragment at this step.
                auto& range = fragment_read_ranges.at(i);
                for (size_t i = range.first; i < range.second; i++) {
                    fragment_covered[i] = true;
                }
            }
        }
        size_t covered_bases = 0;
        for (bool flag : fragment_covered) {
            if (flag) {
                covered_bases++;
            }
        }
        double fragment_overall_coverage = (double) covered_bases / aln.sequence().size();
        best_bucket_fragment_coverage_at_length[threshold] = fragment_overall_coverage;
    }
    // Overall coverage of read with top k fragments by score, in best bucket
    std::vector<double> best_bucket_fragment_coverage_at_top(6, 0.0);
    fragment_covered = std::vector<bool>(aln.sequence().size(), false);
    std::vector<size_t> best_bucket_fragments;
    for (size_t i = 0; i < fragments.size(); i++) {
        if (fragment_results.cluster_nums[i] == best_bucket) {
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
    for (size_t i = 0; i < best_bucket_fragment_coverage_at_top.size() - 2; i++) {
        if (i < best_bucket_fragments.size()) {
            // Add coverage from the fragment at this rank, if any
            auto& range = fragment_read_ranges.at(best_bucket_fragments.at(i));
            for (size_t j = range.first; j < range.second; j++) {
                fragment_covered[j] = true;
            }
        }
    
        // Compute coverage
        size_t covered_bases = 0;
        for (bool flag : fragment_covered) {
            if (flag) {
                covered_bases++;
            }
        }
        double fragment_overall_coverage = (double) covered_bases / aln.sequence().size();
        best_bucket_fragment_coverage_at_top[i + 1] = fragment_overall_coverage;
    }
    
    // Fraction of minimizers with seeds used in fragments of k or more items
    std::vector<size_t> minimizer_fragment_max_items(minimizers.size(), 0);
    std::vector<bool> minimizer_has_seeds(minimizers.size(), false);
    for (auto& seed : seeds) {
        minimizer_has_seeds[seed.source] = true;
    }
    for (auto& fragment : fragments) {
        for (auto& seed_index : fragment.seeds) {
            auto& slot = minimizer_fragment_max_items[seeds[seed_index].source];
            slot = std::max(slot, fragment.seeds.size());
        }
    }
    std::vector<double> seeded_minimizer_fraction_used_in_fragment_of_items;
    seeded_minimizer_fraction_used_in_fragment_of_items.reserve(10);
    for (size_t cutoff = 0; cutoff <= 10; cutoff++) {
        size_t minimizers_eligible = 0;
        size_t fragment_minimizers_used = 0;
        for (size_t i = 0; i < minimizers.size(); i++) {
            if (minimizer_has_seeds[i]) {
                minimizers_eligible++;
                if (minimizer_fragment_max_items[i] >= cutoff) {
                    fragment_minimizers_used++;
                }
            }
        }
        double fraction_used = minimizers_eligible == 0 ? 0.0 : (double) fragment_minimizers_used / minimizers_eligible;
        seeded_minimizer_fraction_used_in_fragment_of_items.push_back(fraction_used);
    }
    
    
    
    // Now we want to find, for each interval, the next interval that starts after it ends
    // So we put all the intervals in an ordered map by start position.
    std::map<size_t, size_t> fragments_by_start;
    // We're also going to need to know which seeds went into which fragments.
    // TODO: We could get away with one seed per fragment here probably.
    // TODO: Can we skip building this if not tracking provenance?
    std::vector<size_t> seed_to_fragment(seeds.size(), std::numeric_limits<size_t>::max());
    for (size_t i = 0; i < fragments.size(); i++) {
        auto found = fragments_by_start.find(fragment_read_ranges[i].first);
        if (found == fragments_by_start.end()) {
            // First thing we've found starting here
            fragments_by_start.emplace_hint(found, fragment_read_ranges[i].first, i);
        } else {
            // When multiple fragments start at a position, we always pick the one with the most seeds.
            // TODO: score the fragments and use the scores?
            if (fragments[found->second].seeds.size() < fragments[i].seeds.size()) {
                // If the one in the map has fewer seeds, replace it.
                found->second = i;
            }
        }
        for (auto& seed : fragments[i].seeds) {
            // Record which fragment this seed went into.
            seed_to_fragment.at(seed) = i;
        }
    }
    // And we need to know the unconnected-to fragments with nothing to their
    // left, which also won the contest for most seeds at their start position
    // (and so could have been connected to)
    std::unordered_set<size_t> unconnected_fragments;
    for (auto& kv : fragments_by_start) {
        unconnected_fragments.insert(kv.second);
    }
    // And then we do bound lookups for each cluster to find the next one
    // And we put those pairs here.
    using fragment_connection_t = std::pair<size_t, size_t>;
    std::vector<fragment_connection_t> fragment_connections;
    for (size_t i = 0; i < fragments.size(); i++) {
        size_t past_end = fragment_read_ranges[i].second;
        // Find the cluster with the most seeds that starts the soonest after the last base in this cluster.
        auto found = fragments_by_start.lower_bound(past_end);
        if (found != fragments_by_start.end()) {
            // We found one. Can we connect them?
            fragment_connections.emplace_back(i, found->second);
            // Something might connect to them
            unconnected_fragments.erase(found->second);
        } else {
            // There's nothing after us, so connect to nowhere.
            fragment_connections.emplace_back(i, std::numeric_limits<size_t>::max());
            if (show_work) {
                #pragma omp critical (cerr)
                std::cerr << log_name() << "Fragment at {R:" << fragment_read_ranges[i].first << "-" << fragment_read_ranges[i].second << "} has nowhere to reseed to" << std::endl;
            }
        }
    }
    for (auto& unconnected : unconnected_fragments) {
        // These fragments could have been connected to but weren't, so look left off of them.
        fragment_connections.emplace_back(std::numeric_limits<size_t>::max(), unconnected);
    }
    
    if (track_provenance) {
        funnel.substage("reseed");
    }
    
    if (track_provenance) {
        // We project all fragments into the funnel
        for (size_t i = 0; i < fragments.size(); i++) {
            funnel.project_group(i, fragments[i].seeds.size());
        }
    }
    
    // Remember how many seeds we had before reseeding
    size_t old_seed_count = seeds.size();
    
    // We are going to need a widget for finding minimizer hit
    // positions in a subgraph, in the right orientation.
    auto find_minimizer_hit_positions = [&](const Minimizer& m, const vector<id_t>& sorted_ids, const std::function<void(const pos_t)>& iteratee) -> void {
        gbwtgraph::hits_in_subgraph(m.hits, m.occs, sorted_ids, [&](pos_t pos, gbwtgraph::payload_type) {
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
    
    // Connections don't appear in the funnel so we track them ourselves.
    size_t fragment_connection_explored_count = 0;
    
    process_until_threshold_a(fragment_connections.size(), (std::function<double(size_t)>) [&](size_t i) -> double {
        // Best pairs to connect are those with the highest average coverage
        if (fragment_connections[i].first == std::numeric_limits<size_t>::max()) {
            return fragments[fragment_connections[i].second].coverage;
        } else if (fragment_connections[i].second == std::numeric_limits<size_t>::max()) {
            return fragments[fragment_connections[i].first].coverage;
        } else {
            return (fragments[fragment_connections[i].first].coverage + fragments[fragment_connections[i].second].coverage) / 2;
        }
    },
    fragment_connection_coverage_threshold,
    min_fragment_connections,
    max_fragment_connections,
    rng,
    [&](size_t connection_num) -> bool {
        // This connection is good enough
        
        // TODO: Add provenance tracking/stage for connections?
    
        // Reseed between each pair of fragments and dump into seeds
        auto& connected = fragment_connections[connection_num];
        
        // Where should we start in the read
        size_t left_read;
        // And in the graph
        pos_t left_pos;
        if (connected.first == std::numeric_limits<size_t>::max()) {
            // Nothing is on the left side of this connection
            left_read = 0;
            left_pos = empty_pos_t();
        } else {
            // Get the information from the fragment on the left side of this connection.
            left_read = fragment_read_ranges[connected.first].second;
            // Make sure graph position points forward along the read.
            left_pos = forward_pos(seeds.at(fragment_bounding_seeds[connected.first].second), minimizers, this->gbwt_graph);
        }
        
        // Where should we end in the read
        size_t right_read;
        // And in the graph
        pos_t right_pos;
        if (connected.second == std::numeric_limits<size_t>::max()) {
            // Nothing is on the right side of this connection
            right_read = aln.sequence().size();
            right_pos = empty_pos_t();
        } else {
            // Get the information from the fragment on the right side of this connection.
            right_read = fragment_read_ranges[connected.second].first;
            // Make sure graph position points forward along the read.
            right_pos = forward_pos(seeds.at(fragment_bounding_seeds[connected.second].first), minimizers, this->gbwt_graph);
        }
        
        if (show_work) {
            if (connected.first == std::numeric_limits<size_t>::max()) {
                #pragma omp critical (cerr)
                {
                    std::cerr << log_name() << "Reseeding before fragment " << connected.second << " at {R:" << right_read << "-" << fragment_read_ranges[connected.second].second << " = G:" << right_pos
                        << "}" << std::endl;
                }
            } else if (connected.second == std::numeric_limits<size_t>::max()) {
                #pragma omp critical (cerr)
                {
                    std::cerr << log_name() << "Reseeding after fragment " << connected.first << " at {R:" << fragment_read_ranges[connected.first].first << "-" << left_read << " = G:" << left_pos
                        << "}" << std::endl;
                }
            } else {
                #pragma omp critical (cerr)
                {
                    std::cerr << log_name() << "Reseeding between fragments " << connected.first << " at {R:" << fragment_read_ranges[connected.first].first << "-" << left_read << " = G:" << left_pos
                        << "} and " << connected.second << " at {R:" << right_read << "-" << fragment_read_ranges[connected.second].second << " = G:" << right_pos
                        << "}" << std::endl;
                }
            }
                    
            // Dump the minimizers in the region
            this->dump_debug_minimizers(minimizers, aln.sequence(), nullptr, left_read, right_read - left_read);
        }
        
        // Do the reseed
        std::vector<Seed> new_seeds = reseed_between(left_read, right_read, left_pos, right_pos, this->gbwt_graph, minimizers, find_minimizer_hit_positions);
        
        // Concatenate and deduplicate with existing seeds
        size_t seeds_before = seeds.size();
        seeds.reserve(seeds_before + new_seeds.size());
        for (auto& seed : new_seeds) {
            // Check if we have seen it before
            std::pair<size_t, pos_t> key {minimizers[seed.source].forward_offset(), seed.pos};
            auto found = seen_seeds.find(key);
            if (found == seen_seeds.end()) {
                // Keep this new seed
                seeds.emplace_back(std::move(seed));
                seen_seeds.emplace_hint(found, std::move(key));
                
                if (this->track_provenance) {
                    funnel.introduce();
                    // Tell the funnel we came from these fragments together
                    if (connected.first != std::numeric_limits<size_t>::max()) {
                        funnel.also_relevant(1, connected.first);
                    }
                    if (connected.second != std::numeric_limits<size_t>::max()) {
                        funnel.also_relevant(1, connected.second);
                    }
                    // TODO: Tie these back to the minimizers, several stages ago.
                }
            }
        }
        
        if (show_work) {
            #pragma omp critical (cerr)
            {
                std::cerr << log_name() << "Found " << new_seeds.size() << " seeds, of which " << (seeds.size() - seeds_before) << " are new" << std::endl;
                std::vector<size_t> new_seeds;
                for (size_t i = seeds_before; i < seeds.size(); i++) {
                    new_seeds.push_back(i);
                } 
                this->dump_debug_seeds(minimizers, seeds, new_seeds);
            }
        }
        
        fragment_connection_explored_count++;
        
        return true;
    }, [&](size_t connection_num) -> void {
        // There are too many sufficiently good connections
        // TODO: Add provenance tracking
    }, [&](size_t connection_num) -> void {
        // This connection is not sufficiently good.
        // TODO: Add provenance tracking
    });
    
    if (this->track_provenance) {
        // Make items in the funnel for all the new seeds, basically as one-seed fragments.
        if (this->track_correctness) {
            // Tag newly introduced seed items with correctness 
            funnel.substage("correct");
        } else {
            // We're just tagging them with read positions
            funnel.substage("placed");
        }
        this->tag_seeds(aln, seeds.cbegin() + old_seed_count, seeds.cend(), minimizers, fragments.size(), funnel);
    }
    
    // Make the main clusters that include the recovered seeds
    if (track_provenance) {
        funnel.stage("cluster");
    }
    
    std::vector<Cluster> clusters = clusterer.cluster_seeds(seeds, chaining_cluster_distance);
    
    // Determine the scores and read coverages for each cluster.
    // Also find the best and second-best cluster scores.
    if (this->track_provenance) {
        funnel.substage("score");
    }
    double best_cluster_score = 0.0, second_best_cluster_score = 0.0;
    for (size_t i = 0; i < clusters.size(); i++) {
        Cluster& cluster = clusters[i];
        
        if (this->track_provenance) {
            // Say we're making it
            funnel.producing_output(i);
        }
        // Since buckets/chains don't straightforwardly merge into clusters we need to completely re-score.
        this->score_cluster(cluster, i, minimizers, seeds, aln.sequence().size());
        // Tell the funnel about where the cluster came from.
        if (this->track_provenance) {
            // Record the cluster in the funnel.
            funnel.introduce();
            funnel.score(funnel.latest(), cluster.score);
            
            // TODO: add source links

            // Say we made it.
            funnel.produced_output();
        }
        if (cluster.score > best_cluster_score) {
            second_best_cluster_score = best_cluster_score;
            best_cluster_score = cluster.score;
        } else if (cluster.score > second_best_cluster_score) {
            second_best_cluster_score = cluster.score;
        }
    }
    
    // Throw out some scratch
    seed_to_fragment.clear();
    seen_seeds.clear();

    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "Found " << clusters.size() << " clusters" << endl;
        }
    }
    
    // We will set a score cutoff based on the best, but move it down to the
    // second best if it does not include the second best and the second best
    // is within pad_cluster_score_threshold of where the cutoff would
    // otherwise be. This ensures that we won't throw away all but one cluster
    // based on score alone, unless it is really bad.
    double cluster_score_cutoff = best_cluster_score - cluster_score_threshold;
    if (cluster_score_cutoff - pad_cluster_score_threshold < second_best_cluster_score) {
        cluster_score_cutoff = std::min(cluster_score_cutoff, second_best_cluster_score);
    }

    if (track_provenance) {
        // Now we go from clusters to chains
        funnel.stage("chain");
    }
    
    chain_config_t chain_cfg;
    
    chain_cfg.max_lookback_bases = this->max_lookback_bases;
    chain_cfg.min_lookback_items = this->min_lookback_items;
    chain_cfg.lookback_item_hard_cap = this->lookback_item_hard_cap;
    chain_cfg.initial_lookback_threshold = this->initial_lookback_threshold;
    chain_cfg.lookback_scale_factor = this->lookback_scale_factor;
    chain_cfg.min_good_transition_score_per_base = this->min_good_transition_score_per_base;
    
    chain_cfg.item_bonus = this->item_bonus;
    chain_cfg.max_indel_bases = this->max_indel_bases;
    
    chain_cfg.cluster_score_cutoff = cluster_score_cutoff;
    chain_cfg.cluster_score_cutoff_enabled = (cluster_score_threshold != 0);
    chain_cfg.cluster_coverage_threshold = this->cluster_coverage_threshold;
    chain_cfg.min_clusters_to_chain = this->min_clusters_to_chain;
    chain_cfg.max_clusters_to_chain = this->max_clusters_to_chain;
    
    chain_cfg.max_chains_per_cluster = 1;
    
    auto chain_results = this->chain_clusters(aln, minimizers, seeds, clusters, chain_cfg, old_seed_count, fragments.size(), funnel, 5, 2, rng);
    // Throw out all but the best chain. There should be one chain per cluster, like we asked.
    vector<pair<int, vector<size_t>>> cluster_chains;
    cluster_chains.reserve(chain_results.cluster_chains.size());
    for (auto& all_chains : chain_results.cluster_chains) {
        cluster_chains.emplace_back(std::move(all_chains.front()));
    }
    auto& cluster_chain_seeds = chain_results.cluster_chain_seeds;
    auto& seed_anchors = chain_results.seed_anchors;
    auto& minimizer_explored = chain_results.minimizer_explored;
    auto& minimizer_kept_cluster_count = chain_results.minimizer_kept_cluster_count;
    auto& kept_cluster_count = chain_results.kept_cluster_count;
    
        
    // We now estimate the best possible alignment score for each cluster.
    std::vector<int> cluster_alignment_score_estimates;
    // Copy cluster chain scores over
    cluster_alignment_score_estimates.resize(cluster_chains.size());
    for (size_t i = 0; i < cluster_chains.size(); i++) {
        cluster_alignment_score_estimates[i] = cluster_chains[i].first;
    }
    
    if (track_provenance) {
        funnel.stage("align");
    }

    //How many of each minimizer ends up in a cluster that actually gets turned into an alignment?
    vector<size_t> minimizer_kept_count(minimizers.size(), 0);
    
    // Now start the alignment step. Everything has to become an alignment.

    // We will fill this with all computed alignments in estimated score order.
    vector<Alignment> alignments;
    alignments.reserve(cluster_alignment_score_estimates.size());
    // This maps from alignment index back to chain index, for
    // tracing back to minimizers for MAPQ. Can hold
    // numeric_limits<size_t>::max() for an unaligned alignment.
    vector<size_t> alignments_to_source;
    alignments_to_source.reserve(cluster_alignment_score_estimates.size());

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
    
    // We need to be able to discard a processed cluster because its score isn't good enough.
    // We have more components to the score filter than process_until_threshold_b supports.
    auto discard_processed_cluster_by_score = [&](size_t processed_num) -> void {
        // This chain is not good enough.
        if (track_provenance) {
            funnel.fail("chain-score", processed_num, cluster_alignment_score_estimates[processed_num]);
        }
        
        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "processed cluster " << processed_num << " failed because its score was not good enough (score=" << cluster_alignment_score_estimates[processed_num] << ")" << endl;
                if (track_correctness && funnel.was_correct(processed_num)) {
                    cerr << log_name() << "\tCORRECT!" << endl;
                }
            }
        }
    };
    
    // Go through the processed clusters in estimated-score order.
    process_until_threshold_b<int>(cluster_alignment_score_estimates,
        chain_score_threshold, min_chains, max_alignments, rng, [&](size_t processed_num) -> bool {
            // This processed cluster is good enough.
            // Called in descending score order.
            
            if (cluster_alignment_score_estimates[processed_num] < chain_min_score) {
                // Actually discard by score
                discard_processed_cluster_by_score(processed_num);
                return false;
            }
            
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "processed cluster " << processed_num << " is good enough (score=" << cluster_alignment_score_estimates[processed_num] << ")" << endl;
                    if (track_correctness && funnel.was_correct(processed_num)) {
                        cerr << log_name() << "\tCORRECT!" << endl;
                    }
                }
            }
            if (track_provenance) {
                funnel.pass("chain-score", processed_num, cluster_alignment_score_estimates[processed_num]);
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
                auto& eligible_seeds = cluster_chain_seeds[processed_num];
                auto& score_and_chain = cluster_chains[processed_num]; 
                vector<size_t>& chain = score_and_chain.second;
                
                // Do the DP between the items in the cluster as specified by the chain we got for it. 
                best_alignments[0] = find_chain_alignment(aln, {seed_anchors, eligible_seeds}, chain);
                    
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
                        cerr << log_name() << "Produced alignment from processed cluster " << processed_num
                            << " with score " << alignments.back().score() << ": " << log_alignment(alignments.back()) << endl;
                    }
                }
            };
            
            for(auto aln_it = best_alignments.begin() ; aln_it != best_alignments.end() && aln_it->score() != 0 && aln_it->score() >= best_alignments[0].score() * 0.8; ++aln_it) {
                //For each additional alignment with score at least 0.8 of the best score
                observe_alignment(*aln_it);
            }

           
            if (track_provenance) {
                // We're done with this input item
                funnel.processed_input();
            }

            for (size_t i = 0 ; i < minimizer_kept_cluster_count[processed_num].size() ; i++) {
                minimizer_kept_count[i] += minimizer_kept_cluster_count[processed_num][i];
                if (minimizer_kept_cluster_count[processed_num][i] > 0) {
                    // This minimizer is in a cluster that gave rise
                    // to at least one alignment, so it is explored.
                    minimizer_explored.insert(i);
                }
            }
            
            return true;
        }, [&](size_t processed_num) -> void {
            // There are too many sufficiently good processed clusters
            if (track_provenance) {
                funnel.pass("chain-score", processed_num, cluster_alignment_score_estimates[processed_num]);
                funnel.fail("max-alignments", processed_num);
            }
            
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "processed cluster " << processed_num << " failed because there were too many good processed clusters (score=" << cluster_alignment_score_estimates[processed_num] << ")" << endl;
                    if (track_correctness && funnel.was_correct(processed_num)) {
                        cerr << log_name() << "\tCORRECT!" << endl;
                    }
                }
            }
        }, discard_processed_cluster_by_score);
    
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
        assert(false);
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

    assert(!mappings.empty());
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
        set_annotation(mappings[0], "param_min-clusters-to-chain", (double) min_clusters_to_chain);
        set_annotation(mappings[0], "param_max-clusters-to-chain", (double) max_clusters_to_chain);
        set_annotation(mappings[0], "param_reseed-search-distance", (double) reseed_search_distance);
        
        // Chaining algorithm parameters
        set_annotation(mappings[0], "param_max-lookback-bases", (double) max_lookback_bases);
        set_annotation(mappings[0], "param_initial-lookback-threshold", (double) initial_lookback_threshold);
        set_annotation(mappings[0], "param_lookback-scale-factor", lookback_scale_factor);
        set_annotation(mappings[0], "param_min-good-transition-score-per-base", min_good_transition_score_per_base);
        set_annotation(mappings[0], "param_item-bonus", (double) item_bonus);
        set_annotation(mappings[0], "param_max-indel-bases", (double) max_indel_bases);
        
        set_annotation(mappings[0], "param_max-chain-connection", (double) max_chain_connection);
        set_annotation(mappings[0], "param_max-tail-length", (double) max_tail_length);
        set_annotation(mappings[0], "param_max-alignments", (double) max_alignments);
        set_annotation(mappings[0], "param_cluster-score", (double) cluster_score_threshold);
        set_annotation(mappings[0], "param_cluster-coverage", (double) cluster_coverage_threshold);
        set_annotation(mappings[0], "param_chain-score", (double) chain_score_threshold);
        set_annotation(mappings[0], "param_chain-min-score", (double) chain_min_score);
        set_annotation(mappings[0], "param_min-chains", (double) min_chains);
        
        set_annotation(mappings[0], "fragment_connections_explored", (double)fragment_connection_explored_count);
        set_annotation(mappings[0], "fragment_connections_total", (double)fragment_connections.size());
    }
    
    // Special fragment statistics
    set_annotation(mappings[0], "fragment_scores", fragment_scores);
    set_annotation(mappings[0], "fragment_item_counts", fragment_item_counts);
    set_annotation(mappings[0], "fragment_bound_coverages", fragment_bound_coverages);
    set_annotation(mappings[0], "best_bucket_fragment_coverage_at_length", best_bucket_fragment_coverage_at_length);
    set_annotation(mappings[0], "best_bucket_fragment_coverage_at_top", best_bucket_fragment_coverage_at_top);
    set_annotation(mappings[0], "bucket_best_fragment_scores", bucket_best_fragment_scores);
    set_annotation(mappings[0], "bucket_scores", bucket_scores);
    set_annotation(mappings[0], "bucket_coverages", bucket_coverages);
    set_annotation(mappings[0], "seeded_minimizer_fraction_used_in_fragment_of_items", seeded_minimizer_fraction_used_in_fragment_of_items);
    
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
        DotDumpExplainer<Funnel> explainer(funnel);
    }

    return mappings;
}

Alignment MinimizerMapper::find_chain_alignment(
    const Alignment& aln,
    const VectorView<algorithms::Anchor>& to_chain,
    const std::vector<size_t>& chain) const {
    
    if (chain.empty()) {
        throw std::logic_error("Cannot find an alignment for an empty chain!");
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
                << " aligns source " << here->source
                << " at " << (*here).read_start() << "-" << (*here).read_end()
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
                throw std::runtime_error(ss.str());
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
                    << " aligns source " << next->source
                    << " at " << (*next).read_start() << "-" << (*next).read_end()
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
                throw std::runtime_error(ss.str());
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
                    cerr << "warning[MinimizerMapper::find_chain_alignment]: Refusing to align " << link_length << " bp connection between chain items " << graph_length << " apart at " << (*here).graph_end() << " and " << (*next).graph_start() << " in " << aln.name() << " to avoid overflow" << endl;
                }
                // Just jump to right tail
                break;
            }
            
            // We can't actually do this alignment, we'd have to align too
            // long of a sequence to find a connecting path.
            #pragma omp critical (cerr)
            {
                cerr << "warning[MinimizerMapper::find_chain_alignment]: Falling back to non-GBWT alignment of " << link_length << " bp connection between chain items " << graph_length << " apart at " << (*here).graph_end() << " and " << (*next).graph_start() << " in " << aln.name() << endl;
            }
            
            Alignment link_aln;
            link_aln.set_sequence(linking_bases);
            if (!aln.quality().empty()) {
                link_aln.set_quality(aln.quality().substr(link_start, link_length));
            }
            assert(graph_length != 0); // TODO: Can't handle abutting graph positions yet
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
                throw std::runtime_error(ss.str());
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

void MinimizerMapper::wfa_alignment_to_alignment(const WFAAlignment& wfa_alignment, Alignment& alignment) const {
    *(alignment.mutable_path()) = wfa_alignment.to_path(this->gbwt_graph, alignment.sequence());
    alignment.set_score(wfa_alignment.score);
    if (!alignment.sequence().empty()) {
        alignment.set_identity(identity(alignment.path()));
    }
}

void MinimizerMapper::with_dagified_local_graph(const pos_t& left_anchor, const pos_t& right_anchor, size_t max_path_length, const HandleGraph& graph, const std::function<void(DeletableHandleGraph&, const std::function<std::pair<nid_t, bool>(const handle_t&)>&)>& callback) {
    
    if (is_empty(left_anchor) && is_empty(right_anchor)) {
        throw std::runtime_error("Cannot align sequence between two unset positions");
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
    
    // And split by strand since we can only align to one strand
    StrandSplitGraph split_graph(&local_graph);
    
    // And make sure it's a DAG of the stuff reachable from our anchors
    bdsg::HashGraph dagified_graph;
    // For which we need the handles that anchor the graph, facing inwards
    std::vector<handle_t> bounding_handles;
    if (!is_empty(left_anchor)) {
        // Dagify from the forward version of the left anchor
        
        // Grab the left anchor in the local graph
        assert(local_graph.has_node(local_left_anchor_id));
        handle_t local_handle = local_graph.get_handle(local_left_anchor_id, is_rev(left_anchor));
        
        // And get the node that that orientation of it is in the strand-split graph
        handle_t overlay_handle = split_graph.get_overlay_handle(local_handle);
        
        // And use that
        bounding_handles.push_back(overlay_handle);
    }
    if (!is_empty(right_anchor)) {
        // Dagify from the reverse version of the node for the forward version of the right anchor
        
        // Grab the right anchor from the local graph
        assert(local_graph.has_node(local_right_anchor_id));
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
            throw std::runtime_error("ID " + std::to_string(dagified_id) + " from dagified graph not found in strand-split graph");
        }
        nid_t split_id = found_in_split->second;
        handle_t split_handle = split_graph.get_handle(split_id, dagified_is_reverse);
        // We rely on get_underlying_handle understanding reversed handles in the split graph
        handle_t local_handle = split_graph.get_underlying_handle(split_handle);
        nid_t local_id = local_graph.get_id(local_handle);
        bool local_is_reverse = local_graph.get_is_reverse(local_handle);
        auto found_in_base = local_to_base.find(local_id);
        if (found_in_base == local_to_base.end()) {
            throw std::runtime_error("ID " + std::to_string(local_id) + " from local graph not found in full base graph");
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
    
        // Now the alignment is filled in!
    });
}

std::vector<algorithms::Anchor> MinimizerMapper::to_anchors(const Alignment& aln, const VectorView<Minimizer>& minimizers, const std::vector<Seed>& seeds) const {
    std::vector<algorithms::Anchor> to_return;
    to_return.reserve(seeds.size());
    for (auto& seed : seeds) {
        to_return.push_back(this->to_anchor(aln, minimizers, seed));
    }
    return to_return;
}

algorithms::Anchor MinimizerMapper::to_anchor(const Alignment& aln, const VectorView<Minimizer>& minimizers, const Seed& seed) const {
    // Turn each seed into the part of its match on the node where the
    // anchoring end (start for forward-strand minimizers, ane for
    // reverse-strand minimizers) falls.
    auto& source = minimizers[seed.source];
    size_t length;
    pos_t graph_start;
    size_t read_start;
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
    } else {
        // Seed stores the first base of the match in the graph
        graph_start = seed.pos;
        
        // Get the handle to the node it's on.
        handle_t start_handle = gbwt_graph.get_handle(id(graph_start), is_rev(graph_start));
        // Work out how much of the node it could use before there.
        length = std::min((size_t) source.length, gbwt_graph.get_length(start_handle) - offset(graph_start));
        
        // And we store the read start position already in the item
        read_start = source.value.offset;
    }
    // Work out how many points the anchor is
    // TODO: Always make sequence and quality available for scoring!
    int score = get_regular_aligner()->score_exact_match(aln, read_start, length);
    return algorithms::Anchor(read_start, graph_start, length, score); 
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
