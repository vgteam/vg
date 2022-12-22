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
                                           const std::vector<size_t>& seed_to_precluster,
                                           const std::vector<Cluster>& preclusters,
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
    SmallBitset preclusters_seen(preclusters.size());
    

    // Determine the minimizers that are present in the cluster.
    for (auto hit_index : cluster.seeds) {
        // We have this seed's minimizer
        cluster.present.insert(seeds[hit_index].source);
        
        if (hit_index < first_new_seed) {
            // An old seed.
            // We can also pick up an old cluster.
            size_t old_cluster = seed_to_precluster.at(hit_index);
            if (old_cluster != std::numeric_limits<size_t>::max()) {
                // This seed came form an old cluster, so we must have eaten it
                if (!preclusters_seen.contains(old_cluster)) {
                    // Remember we used this old cluster
                    to_combine.push_back(old_cluster);
                    preclusters_seen.insert(old_cluster);
                }
            }
        } else {
            // Make sure we tell the funnel we took in this new seed.
            // Translate from a space that is old seeds and then new seeds to a
            // space that is old *clusters* and then new seeds
            to_combine.push_back(hit_index - first_new_seed + preclusters.size());
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
        funnel.stage("precluster");
        funnel.substage("compute-preclusters");
    }

    // Find the clusters up to a flat distance limit
    std::vector<Cluster> preclusters = clusterer.cluster_seeds(seeds, chaining_cluster_distance);
    
    if (track_provenance) {
        funnel.substage("score-preclusters");
    }
    for (size_t i = 0; i < preclusters.size(); i++) {
        Cluster& precluster = preclusters[i];
        this->score_cluster(precluster, i, minimizers, seeds, aln.sequence().length(), funnel);
    }
    
    // Find pairs of "adjacent" preclusters
    if (track_provenance) {
        funnel.substage("pair-preclusters");
    }
    
    // To do that, we need start end end positions for each precluster, in the read
    std::vector<std::pair<size_t, size_t>> precluster_read_ranges(preclusters.size(), {std::numeric_limits<size_t>::max(), 0});
    // And the lowest-numbered seeds in the precluster from those minimizers.
    std::vector<std::pair<size_t, size_t>> precluster_bounding_seeds(preclusters.size(), {std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()});
    for (size_t i = 0; i < preclusters.size(); i++) {
        // For each precluster
        auto& precluster = preclusters[i];
        // We will fill in the range it ocvcupies in the read
        auto& read_range = precluster_read_ranges[i];
        auto& graph_seeds = precluster_bounding_seeds[i];
        for (auto& seed_index : precluster.seeds) {
            // Which means we look at the minimizer for each seed
            auto& minimizer = minimizers[seeds[seed_index].source];
            
            if (minimizer.forward_offset() < read_range.first) {
                // Min all their starts to get the precluster start
                read_range.first = minimizer.forward_offset();
                if (seed_index < graph_seeds.first) {
                    // And keep a seed hit
                    graph_seeds.first = seed_index;
                }
            }
            
            if (minimizer.forward_offset() + minimizer.length > read_range.second) {
                // Max all their past-ends to get the precluster past-end
                read_range.second = minimizer.forward_offset() + minimizer.length;
                if (seed_index < graph_seeds.second) {
                    // And keep a seed hit
                    graph_seeds.second = seed_index;
                }
            }
        }
    }
    
    // Now we want to find, for each interval, the next interval that starts after it ends
    // So we put all the intervals in an ordered map by start position.
    std::map<size_t, size_t> preclusters_by_start;
    // We're also going to need to know which seeds went into which preclusters.
    // TODO: We could get away with one seed per precluster here probably.
    // TODO: Can we skip building this if not tracking provenance?
    std::vector<size_t> seed_to_precluster(seeds.size(), std::numeric_limits<size_t>::max());
    for (size_t i = 0; i < preclusters.size(); i++) {
        auto found = preclusters_by_start.find(precluster_read_ranges[i].first);
        if (found == preclusters_by_start.end()) {
            // First thing we've found starting here
            preclusters_by_start.emplace_hint(found, precluster_read_ranges[i].first, i);
        } else {
            // When multiple preclusters start at a position, we always pick the one with the most seeds.
            // TODO: score the preclusters and use the scores?
            if (preclusters[found->second].seeds.size() < preclusters[i].seeds.size()) {
                // If the one in the map has fewer seeds, replace it.
                found->second = i;
            }
        }
        for (auto& seed : preclusters[i].seeds) {
            // Record which precluster this seed went into.
            seed_to_precluster.at(seed) = i;
        }
    }
    // And we need to know the unconnected-to preclusters with nothing to their
    // left, which also won the contest for most seeds at their start position
    // (and so could have been connected to)
    std::unordered_set<size_t> unconnected_preclusters;
    for (auto& kv : preclusters_by_start) {
        unconnected_preclusters.insert(kv.second);
    }
    // And then we do bound lookups for each cluster to find the next one
    // And we put those pairs here.
    using precluster_connection_t = std::pair<size_t, size_t>;
    std::vector<precluster_connection_t> precluster_connections;
    for (size_t i = 0; i < preclusters.size(); i++) {
        size_t past_end = precluster_read_ranges[i].second;
        // Find the cluster with the most seeds that starts the soonest after the last base in this cluster.
        auto found = preclusters_by_start.lower_bound(past_end);
        if (found != preclusters_by_start.end()) {
            // We found one. Can we connect them?
            precluster_connections.emplace_back(i, found->second);
            // Something might connect to them
            unconnected_preclusters.erase(found->second);
        } else {
            // There's nothing after us, so connect to nowhere.
            precluster_connections.emplace_back(i, std::numeric_limits<size_t>::max());
            if (show_work) {
                #pragma omp critical (cerr)
                std::cerr << log_name() << "Precluster at {R:" << precluster_read_ranges[i].first << "-" << precluster_read_ranges[i].second << "} has nowhere to reseed to" << std::endl;
            }
        }
    }
    for (auto& unconnected : unconnected_preclusters) {
        // These preclusters could have been connected to but weren't, so look left off of them.
        precluster_connections.emplace_back(std::numeric_limits<size_t>::max(), unconnected);
    }
    
    if (track_provenance) {
        funnel.stage("reseed");
    }
    
    if (track_provenance) {
        // We project all preclusters into the funnel
        for (size_t i = 0; i < preclusters.size(); i++) {
            funnel.project_group(i, preclusters[i].seeds.size());
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
    
    process_until_threshold_a(precluster_connections.size(), (std::function<double(size_t)>) [&](size_t i) -> double {
        // Best pairs to connect are those with the highest average coverage
        if (precluster_connections[i].first == std::numeric_limits<size_t>::max()) {
            return preclusters[precluster_connections[i].second].coverage;
        } else if (precluster_connections[i].second == std::numeric_limits<size_t>::max()) {
            return preclusters[precluster_connections[i].first].coverage;
        } else {
            return (preclusters[precluster_connections[i].first].coverage + preclusters[precluster_connections[i].second].coverage) / 2;
        }
    },
    precluster_connection_coverage_threshold,
    min_precluster_connections,
    max_precluster_connections,
    rng,
    [&](size_t precluster_num) -> bool {
        // This connection is good enough
        
        // TODO: Add provenance tracking/stage for connections?
    
        // Reseed between each pair of preclusters and dump into seeds
        auto& connected = precluster_connections[precluster_num];
        
        // Where should we start in the read
        size_t left_read;
        // And in the graph
        pos_t left_pos;
        if (connected.first == std::numeric_limits<size_t>::max()) {
            // Nothing is on the left side of this connection
            left_read = 0;
            left_pos = empty_pos_t();
        } else {
            // Get the information from the precluster on the left side of this connection.
            left_read = precluster_read_ranges[connected.first].second;
            // Make sure graph position points forward along the read.
            left_pos = forward_pos(seeds.at(precluster_bounding_seeds[connected.first].second), minimizers, this->gbwt_graph);
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
            // Get the information from the precluster on the right side of this connection.
            right_read = precluster_read_ranges[connected.second].first;
            // Make sure graph position points forward along the read.
            right_pos = forward_pos(seeds.at(precluster_bounding_seeds[connected.second].first), minimizers, this->gbwt_graph);
        }
        
        if (show_work) {
            if (connected.first == std::numeric_limits<size_t>::max()) {
                #pragma omp critical (cerr)
                {
                    std::cerr << log_name() << "Reseeding before precluster " << connected.second << " at {R:" << right_read << "-" << precluster_read_ranges[connected.second].second << " = G:" << right_pos
                        << "}" << std::endl;
                }
            } else if (connected.second == std::numeric_limits<size_t>::max()) {
                #pragma omp critical (cerr)
                {
                    std::cerr << log_name() << "Reseeding after precluster " << connected.first << " at {R:" << precluster_read_ranges[connected.first].first << "-" << left_read << " = G:" << left_pos
                        << "}" << std::endl;
                }
            } else {
                #pragma omp critical (cerr)
                {
                    std::cerr << log_name() << "Reseeding between preclusters " << connected.first << " at {R:" << precluster_read_ranges[connected.first].first << "-" << left_read << " = G:" << left_pos
                        << "} and " << connected.second << " at {R:" << right_read << "-" << precluster_read_ranges[connected.second].second << " = G:" << right_pos
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
                    // Tell the funnel we came from these preclusters together
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
        
        return true;
    }, [&](size_t connection_num) -> void {
        // There are too many sufficiently good connections
        // TODO: Add provenance tracking
    }, [&](size_t connection_num) -> void {
        // This connection is not sufficiently good.
        // TODO: Add provenance tracking
    });
    
    if (this->track_provenance) {
        // Make items in the funnel for all the new seeds, basically as one-seed preclusters.
        if (this->track_correctness) {
            // Tag newly introduced seed items with correctness 
            funnel.substage("correct");
        } else {
            // We're just tagging them with read positions
            funnel.substage("placed");
        }
        this->tag_seeds(aln, seeds.cbegin() + old_seed_count, seeds.cend(), minimizers, preclusters.size(), funnel);
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
        this->score_merged_cluster(cluster,
                                   i,
                                   minimizers,
                                   seeds,
                                   old_seed_count,
                                   seed_to_precluster,
                                   preclusters,
                                   aln.sequence().length(),
                                   funnel);
        if (cluster.score > best_cluster_score) {
            second_best_cluster_score = best_cluster_score;
            best_cluster_score = cluster.score;
        } else if (cluster.score > second_best_cluster_score) {
            second_best_cluster_score = cluster.score;
        }
    }
    
    // Throw out some scratch
    seed_to_precluster.clear();
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
    
    // Convert the seeds into chainable anchors in the same order
    vector<algorithms::Anchor> seed_anchors = this->to_anchors(aln, minimizers, seeds);
    
    // These are the chains for all the clusters, as score and sequence of visited seeds.
    vector<pair<int, vector<size_t>>> cluster_chains;
    cluster_chains.reserve(clusters.size());
    
    // To compute the windows for explored minimizers, we need to get
    // all the minimizers that are explored.
    SmallBitset minimizer_explored(minimizers.size());
    //How many hits of each minimizer ended up in each cluster we kept?
    vector<vector<size_t>> minimizer_kept_cluster_count; 

    size_t kept_cluster_count = 0;
    
    // What cluster seeds define the space for clusters' chosen chains?
    vector<vector<size_t>> cluster_chain_seeds;
    
    //Process clusters sorted by both score and read coverage
    process_until_threshold_c<double>(clusters.size(), [&](size_t i) -> double {
            return clusters[i].coverage;
        }, [&](size_t a, size_t b) -> bool {
            return ((clusters[a].coverage > clusters[b].coverage) ||
                    (clusters[a].coverage == clusters[b].coverage && clusters[a].score > clusters[b].score));
        }, cluster_coverage_threshold, min_clusters_to_chain, max_clusters_to_chain, rng, [&](size_t cluster_num) -> bool {
            // Handle sufficiently good clusters in descending coverage order
            
            Cluster& cluster = clusters[cluster_num];
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
            if (cluster_score_threshold != 0 && cluster.score < cluster_score_cutoff 
                && kept_cluster_count >= min_clusters_to_chain) {
                //If the score isn't good enough and we already kept at least min_clusters_to_chain clusters,
                //ignore this cluster
                if (track_provenance) {
                    funnel.fail("cluster-score", cluster_num, cluster.score);
                }
                if (show_work) {
                    #pragma omp critical (cerr)
                    {
                        cerr << log_name() << "Cluster " << cluster_num << " fails cluster score cutoff" <<  endl;
                        cerr << log_name() << "Covers " << clusters[cluster_num].coverage << "/best-" << cluster_coverage_threshold << " of read" << endl;
                        cerr << log_name() << "Involves " << cluster_node_count << " nodes in " << cluster_min_node << "-" << cluster_max_node << endl;
                        cerr << log_name() << "Scores " << clusters[cluster_num].score << "/" << cluster_score_cutoff << endl;
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
                    cerr << log_name() << "Covers " << cluster.coverage << "/best-" << cluster_coverage_threshold << " of read" << endl;
                    cerr << log_name() << "Involves " << cluster_node_count << " nodes in " << cluster_min_node << "-" << cluster_max_node << endl;
                    cerr << log_name() << "Scores " << cluster.score << "/" << cluster_score_cutoff << endl;
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
            cluster_chains.emplace_back();
            cluster_chains.back().first = std::numeric_limits<int>::min();
            cluster_chain_seeds.emplace_back();
                
            // Find a chain from this cluster
            VectorView<algorithms::Anchor> cluster_view {seed_anchors, cluster_seeds_sorted};
            auto candidate_chain = algorithms::find_best_chain(cluster_view,
                                                               *distance_index,
                                                               gbwt_graph,
                                                               get_regular_aligner()->gap_open,
                                                               get_regular_aligner()->gap_extension,
                                                               max_lookback_bases,
                                                               min_lookback_items,
                                                               initial_lookback_threshold,
                                                               lookback_scale_factor,
                                                               min_good_transition_score_per_base,
                                                               item_bonus,
                                                               max_indel_bases);
            if (show_work && !candidate_chain.second.empty()) {
                #pragma omp critical (cerr)
                {
                    
                    cerr << log_name() << "Cluster " << cluster_num << " running " << seed_anchors[cluster_seeds_sorted.front()] << " to " << seed_anchors[cluster_seeds_sorted.back()]
                        << " has chain with score " << candidate_chain.first
                        << " and length " << candidate_chain.second.size()
                        << " running R" << cluster_view[candidate_chain.second.front()].read_start()
                        << " to R" << cluster_view[candidate_chain.second.back()].read_end() << std::endl;
                }
            }
            if (candidate_chain.first > cluster_chains.back().first) {
                // Keep it if it is better
                cluster_chains.back() = std::move(candidate_chain);
                cluster_chain_seeds.back() = cluster_seeds_sorted;
            }
            
            if (track_provenance) {
                funnel.substage_stop();
            }
            
            if (track_provenance) {
                // Record with the funnel that there is now a chain that comes
                // from all the seeds that participate in the chain.
                funnel.introduce();
                funnel.score(funnel.latest(), cluster_chains.back().first);
                // Accumulate the old and new seed funnel numbers to connect to.
                // TODO: should we just call into the funnel every time instead of allocating?
                std::vector<size_t> old_seed_ancestors;
                std::vector<size_t> new_seed_ancestors;
                for (auto& sorted_seed_number : cluster_chains.back().second) {
                    // Map each seed back to its canonical seed order
                    size_t seed_number = cluster_chain_seeds.back().at(sorted_seed_number);
                    if (seed_number < old_seed_count) {
                        // Seed is original, from "seed" stage 4 stages ago
                        old_seed_ancestors.push_back(seed_number);
                    } else {
                        // Seed is new, from "reseed" stage 2 stages ago. Came
                        // after all the preclusters which also live in the reseed stage.
                        new_seed_ancestors.push_back(seed_number - old_seed_count + preclusters.size());
                    }
                }
                // We came from all the original seeds, 4 stages ago
                funnel.also_merge_group(4, old_seed_ancestors.begin(), old_seed_ancestors.end());
                // We came from all the new seeds, 2 stages ago
                funnel.also_merge_group(2, new_seed_ancestors.begin(), new_seed_ancestors.end());
                // We're also related to the source cluster from the
                // immediately preceeding stage.
                funnel.also_relevant(1, cluster_num);
                
                // Say we finished with this cluster, for now.
                funnel.processed_input();
            }
            
            return true;
            
        }, [&](size_t cluster_num) -> void {
            // There are too many sufficiently good clusters
            Cluster& cluster = clusters[cluster_num];
            if (track_provenance) {
                funnel.pass("cluster-coverage", cluster_num, cluster.coverage);
                funnel.fail("max-clusters-to-chain", cluster_num);
            }
            
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    
                    cerr << log_name() << "Cluster " << cluster_num << " passes cluster cutoffs but we have too many" <<  endl;
                    cerr << log_name() << "Covers " << cluster.coverage << "/best-" << cluster_coverage_threshold << " of read" << endl;
                    cerr << log_name() << "Scores " << cluster.score << "/" << cluster_score_cutoff << endl;
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
                    cerr << log_name() << "Covers " << clusters[cluster_num].coverage << "/best-" << cluster_coverage_threshold << " of read" << endl;
                    cerr << log_name() << "Scores " << clusters[cluster_num].score << "/" << cluster_score_cutoff << endl;
                }
            }
        });
        
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
        extension_set_score_threshold, min_extension_sets, max_alignments, rng, [&](size_t processed_num) -> bool {
            // This processed cluster is good enough.
            // Called in descending score order.
            
            if (cluster_alignment_score_estimates[processed_num] < extension_set_min_score) {
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
                // Traceback over the remaining extensions?
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
            annotate_with_minimizer_statistics(mappings[0], minimizers, seeds, funnel);
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
        set_annotation(mappings[0], "param_precluster-connection-coverage-threshold", precluster_connection_coverage_threshold);
        set_annotation(mappings[0], "param_min-precluster-connections", (double) min_precluster_connections);
        set_annotation(mappings[0], "param_max-precluster-connections", (double) max_precluster_connections);
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
        set_annotation(mappings[0], "param_max-multimaps", (double) max_multimaps);
    }
    
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
                cerr << "warning[MinimizerMapper::find_chain_alignment]: Falling back to non-GBWT alignment of " << left_tail_length << " bp left tail in " << aln.name() << endl;
            }
            
            Alignment tail_aln;
            tail_aln.set_sequence(left_tail);
            if (!aln.quality().empty()) {
                tail_aln.set_quality(aln.quality().substr(0, left_tail_length));
            }
            
            // Work out how far the tail can see
            size_t graph_horizon = left_tail_length + this->get_regular_aligner()->longest_detectable_gap(aln, aln.sequence().begin());
            // Align the left tail, anchoring the right end.
            align_sequence_between(empty_pos_t(), right_anchor, graph_horizon, &this->gbwt_graph, this->get_regular_aligner(), tail_aln);
            // Since it's the left tail we can just clobber the path
            composed_path = tail_aln.path();
            composed_score = tail_aln.score();
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
            
            // We can't actually do this alignment, we'd have to align too
            // long of a sequence to find a connecting path.
            #pragma omp critical (cerr)
            {
                cerr << "warning[MinimizerMapper::find_chain_alignment]: Falling back to non-GBWT alignment of " << link_length << " bp connection between chain items " << graph_length << " apart in " << aln.name() << endl;
            }
            
            Alignment link_aln;
            link_aln.set_sequence(linking_bases);
            if (!aln.quality().empty()) {
                link_aln.set_quality(aln.quality().substr(link_start, link_length));
            }
            assert(graph_length != 0); // TODO: Can't handle abutting graph positions yet
            // Guess how long of a graph path we ought to allow in the alignment.
            size_t path_length = std::max(graph_length, link_length) + this->get_regular_aligner()->longest_detectable_gap(aln, aln.sequence().begin() + link_start);
            MinimizerMapper::align_sequence_between((*here).graph_end(), (*next).graph_start(), path_length, &this->gbwt_graph, this->get_regular_aligner(), link_aln);
            
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

            #pragma omp critical (cerr)
            {
                cerr << "warning[MinimizerMapper::find_chain_alignment]: Falling back to non-GBWT alignment of " << left_tail_length << " bp right tail in " << aln.name() << endl;
            }
            
            Alignment tail_aln;
            tail_aln.set_sequence(right_tail);
            if (!aln.quality().empty()) {
                tail_aln.set_quality(aln.quality().substr((*here).read_end(), right_tail_length));
            }
            
            // Work out how far the tail can see
            size_t graph_horizon = right_tail_length + this->get_regular_aligner()->longest_detectable_gap(aln, aln.sequence().begin() + (*here).read_end());
            // Align the right tail, anchoring the left end.
            align_sequence_between(left_anchor, empty_pos_t(), graph_horizon, &this->gbwt_graph, this->get_regular_aligner(), tail_aln);
            // Since it's the right tail we have to add it on
            append_path(composed_path, tail_aln.path());
            composed_score += tail_aln.score();
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

void MinimizerMapper::align_sequence_between(const pos_t& left_anchor, const pos_t& right_anchor, size_t max_path_length, const HandleGraph* graph, const GSSWAligner* aligner, Alignment& alignment) {
    
    if (is_empty(left_anchor) && is_empty(right_anchor)) {
        throw std::runtime_error("Cannot align sequence between two unset positions");
    }
    
    // We need to get the graph to align to.
    bdsg::HashGraph local_graph;
    unordered_map<id_t, id_t> local_to_base;
    if (!is_empty(left_anchor) && !is_empty(right_anchor)) {
        // We want a graph actually between two positions
        local_to_base = algorithms::extract_connecting_graph(
            graph,
            &local_graph,
            max_path_length,
            left_anchor, right_anchor
        );
    } else if (!is_empty(left_anchor)) {
        // We only have the left anchor
        local_to_base = algorithms::extract_extending_graph(
            graph,
            &local_graph,
            max_path_length,
            left_anchor,
            false,
            false
        );
    } else {
        // We only have the right anchor
        local_to_base = algorithms::extract_extending_graph(
            graph,
            &local_graph,
            max_path_length,
            right_anchor,
            true,
            false
        );
    }
    
    // And split by strand since we can only align to one strand
    StrandSplitGraph split_graph(&local_graph);
    
    // And make sure it's a DAG of the stuff reachable from our anchors
    bdsg::HashGraph dagified_graph;
    // For which we need the handles that anchor the graph, facing inwards
    std::vector<handle_t> bounding_handles;
    if (!is_empty(left_anchor)) {
        // Dagify from the forward version of the left anchor
        bounding_handles.push_back(split_graph.get_overlay_handle(local_graph.get_handle(id(left_anchor), is_rev(left_anchor))));
    }
    if (!is_empty(right_anchor)) {
        // Dagify from the reverse version of the node for the forward version of the right anchor
        bounding_handles.push_back(split_graph.flip(split_graph.get_overlay_handle(local_graph.get_handle(id(right_anchor), is_rev(right_anchor)))));
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
        // to go from a source to a sink.
        aligner->align_global_banded(alignment, dagified_graph);
    } else {
        // Do pinned alignment off the anchor we actually have, with x-drop
        aligner->align_pinned(alignment, dagified_graph, !is_empty(left_anchor), true);
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
