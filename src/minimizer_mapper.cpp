/**
 * \file minimizer_mapper.cpp
 * Defines the code for the minimizer-and-GBWT-based mapper.
 */

#include "minimizer_mapper.hpp"
#include "annotation.hpp"
#include "path_subgraph.hpp"
#include "multipath_alignment.hpp"
#include "funnel.hpp"

#include "algorithms/dijkstra.hpp"

#include <iostream>
#include <algorithm>
#include <cmath>

namespace vg {

using namespace std;

MinimizerMapper::MinimizerMapper(const gbwtgraph::GBWTGraph& graph,
    const std::vector<std::unique_ptr<gbwtgraph::DefaultMinimizerIndex>>& minimizer_indexes,
    MinimumDistanceIndex& distance_index, const PathPositionHandleGraph* path_graph) :
    path_graph(path_graph), minimizer_indexes(minimizer_indexes),
    distance_index(distance_index), gbwt_graph(graph),
    extender(gbwt_graph, *(get_regular_aligner())), clusterer(distance_index) {
    
    // Nothing to do!
}

#define debug

double MinimizerMapper::cheapest_cover_cost(const vector<Minimizer>& minimizers,
    vector<size_t>& broken, const string& sequence, const string& quality_bytes) const {
    
    if (broken.empty() || quality_bytes.empty()) {
        // If we have no agglomerations or no qualities, bail
        return numeric_limits<double>::infinity();
    }
    
    assert(sequence.size() == quality_bytes.size());

    // Sort the agglomerations by start position
    std::sort(broken.begin(), broken.end(), [&](const size_t& a, const size_t& b) {
        return minimizers[a].agglomeration_start < return minimizers[b].agglomeration_start;
    });

    // A window in flight is a pair of start position, inclusive end position
    struct window_t {
        size_t first;
        size_t last;
    };

    // Have a priority queue of window starts and ends, prioritized earliest-ending best, and then latest-starting best.
    // This is less, and greater priority is at the front of the queue (better).
    auto priority = [&](const window_t& a, const window_t& b) {
        // Returns true if a is worse (ends later, or ends at the same place and starts earlier).
        return (a.last > b.last) || (a.last == b.last && a.first < b.first);
    };
    priority_queue<window_t, vector<window_t>, decltype(priority)> active_windows(priority);
    
    // Have a cursor for which agglomeration should come in next.
    auto next_agglomeration = broken.begin();
    
    // Have a DP table with the cost of the cheapest solution to the problem up to here, including a hit at this base.
    // Or numeric_limits<double>::infinity() if base cannot be hit.
    vector<double> costs;
    costs.reserve(sequence.size());
    
    // Keep track of the latest-starting window ending before here. If none, this will be two numeric_limits<size_t>::max() values.
    window_t latest_starting_ending_before = { numeric_limits<size_t>::max(), numeric_limits<size_t>::max() };
    
    for (size_t i = 0; i < sequence.size(); i++) {
        // For each base in the read
        
#ifdef debug
        cerr << "At base " << i << endl;
#endif
        
        // We have the start and end of the latest-starting window ending before here (may be none)
        
        if (isATGC(sequence[i]) &&
            (!active_windows.empty() || 
            (next_agglomeration != broken.end() && minimizers[*next_agglomeration].agglomeration_start == i))) {
            
            // This base is not N, and it is either covered by an agglomeration
            // that hasn't ended yet, or a new agglomeration starts here.
            
#ifdef debug
            cerr << "\tBase is acceptable (" << sequence[i] << ", " << active_windows.size() << " active windows, "
                << ((next_agglomeration != broken.end() && minimizers[*next_agglomeration].agglomeration_start == i) ? "new starting" : "") 
                << ")" << endl;
#endif
            
            // Look at the start of that latest-starting window ending before here.
            
            if (latest_starting_ending_before.first == numeric_limits<size_t>::max()) {
                // If there is no such window, this is the first base hit, so
                // record the cost of hitting it.
                
                costs.push_back(quality_bytes[i]);
                
#ifdef debug
                cerr << "\tFirst base hit, costs " << costs.back() << endl;
#endif
            } else {
                // Else, scan from that window's start to its end in the DP
                // table, and find the min cost.
                auto min_prev_cost_at = std::min_element(costs.begin() + latest_starting_ending_before.first,
                    costs.begin() + latest_starting_ending_before.last + 1);
                double min_prev_cost = *min_prev_cost_at;
                    
                // Add the cost of hitting this base
                costs.push_back(min_prev_cost + quality_bytes[i]);
                
#ifdef debug
                cerr << "\tComes from prev base at " << (min_prev_cost_at - costs.begin()) << ", costs " << costs.back() << endl;
#endif
            }
            
        } else {
            // This base is N, or not covered by an agglomeration.
            // Say we can't hit it.
            costs.push_back(numeric_limits<double>::infinity());
        }
        
        // Now we compute the start of the latest-starting window ending here or before (and update the active windows).
        
        while (next_agglomeration != broken.end() && minimizers[*next_agglomeration].agglomeration_start == i) {
            // While the next agglomeration starts here
            
            // Look it up
            auto& minimizer = minimizers[*next_agglomeration];
            
            // Determine its window size from its index.
            auto& source_index = *minimizer_indexes[minimizer.source];
            size_t window_size = source_index.k() + source_index.w() - 1;
            
            
            for (size_t start = minimizer.agglomeration_start; start + window_size - 1 < minimizer.agglomeration_length; start++) {
                // Add all the agglomeration's windows to the queue
                active_windows.emplace(start, start + window_size - 1);
            }
            
#ifdef debug
            cerr << "\tBegin agglomeration of " << (minimizer.agglomeration_length - window_size + 1)
                << " windows of " << window_size << " bp each" << endl;
#endif
            
            // And advance the cursor
            ++next_agglomeration;
        }

        while (active_windows.top().last == i) {
            // The look at the queue to see if a window ends here. This is second so that we can handle 1-base windows.
            
#ifdef debug
            cerr << "\tEnd window " << active_windows.top().first << " - " << active_windows.top().last << endl;
#endif
            
            if (latest_starting_ending_before.first == numeric_limits<size_t>::max() ||
                active_windows.top().first > latest_starting_ending_before.first) {
                
#ifdef debug
                cerr << "\t\tNew latest-starting-before-here!" << endl;
#endif
                
                // If so, use the latest-starting of all such windows as our latest starting window ending here or before result.
                latest_starting_ending_before = active_windows.top();
            }
            
            
            
            // And pop them all off.
            active_windows.pop();
        }
        // If not, use the latest-starting window ending at the previous base or before (i.e. do nothing).
        
        // Loop around; we will have the latest-starting window ending before the next here.
    }
    
    // When we get here, all the agglomerations should have been handled
    assert(next_agglomeration == broken.end());
    // And all the windows should be off the queue.
    assert(active_windows.empty());
    
    // When we get to the end, we have the latest-starting window overall. It must exist.
    assert(latest_starting_ending_before.first != numeric_limits<size_t>::max());
    
    // Scan it for the best final base to hit and return the cost there.
    auto min_cost_at = std::min_element(costs.begin() + latest_starting_ending_before.first,
        costs.begin() + latest_starting_ending_before.last + 1);
#ifdef debug
    cerr << "Overall min cost: " << *min_cost_at << " at base " << (min_cost_at - costs.begin()) << endl;
#endif
    return *min_cost_at;
}

double MinimizerMapper::window_breaking_quality(const vector<Minimizer>& minimizers, vector<size_t>& broken,
    const string& sequence, const string& quality_bytes) const {
    
    return cheapest_cover_cost(minimizers, broken, sequence, quality_bytes);
    
   
    // Work out how big a window that gives rise to a minimizer actually is.
    // Use the kmer size and the number of 1-bp-slid kmers in a window
    size_t window_size = minimizer_indexes[0]->k() + minimizer_indexes[0]->w() - 1;
   
    // For now, sum up a totla windows figure for all the minimizers form index 0.
    // TODO: Implement the real algorithm.
    size_t total_windows = 0;
    
    for (auto& minimizer_num : broken) {
        // For each minimizer we have to break/create
        auto& m = minimizers[minimizer_num];
        if (m.origin != 0) {
            // Skip minimizers from other indexes
            continue;
        }
        
        // Compute a real sequence start position
        size_t start_pos = m.value.offset;
        if (m.value.is_reverse) {
            // We have the last base and we need the first base.
            start_pos -= (minimizer_indexes[m.origin]->k() - 1);
        }
        
#ifdef debug
        cerr << "Minimizer " << minimizer_num << " (" << m.value.key.decode(minimizer_indexes[m.origin]->k())
            << ") @ " << m.value.offset  << " " << (m.value.is_reverse ? '-' : '+') << " has agglomeration at "
            << m.agglomeration_start << " running " << m.agglomeration_length
            << " bp and contributes " << (m.agglomeration_length - window_size + 1)
            << " windows of size " << window_size << endl;
#endif
        
        // Count the number of slid windows in the agglomeration.
        // TODO: we assume windows aren't double-credited, which admittedly is rare.
        total_windows += m.agglomeration_length - window_size + 1;
    }
   
#ifdef debug
    cerr << "Cap for " << total_windows << ": ";
#endif
    if (total_windows == 0 || total_windows == numeric_limits<size_t>::max()) {
        // No limit can be imposed (because no windows need breaking, or
        // because we got a sentinel max value which means basically the same
        // thing, or because we are using multiple minimizer indexes and don't
        // keep windows straight between them).
#ifdef debug
        cerr << "uncapped" << endl;
#endif
        return numeric_limits<double>::infinity();
    } else {
        // Some windows were located, so we can use this algorithm safely
        
        // How many bases need to be modified to have been modified to have
        // created all of these windows?
        size_t minimum_disruptive_bases = 0;
        
        // Compute how the windows could most overlap, since we forgot where
        // they actually were in the read. If they all overlap, they cover the
        // window size, plus one for each window after the first.
        size_t overlapped_window_length = window_size + total_windows - 1;
        // Since each base can disrupt window length - 1 windows on either side, we
        // tile the overlapped windows with windows to work out how many bases
        // would need to change to create them all. This is at least 1.
        minimum_disruptive_bases = overlapped_window_length / window_size;
        
        // Then we find that many of the lowest base qualities that aren't for Ns (which don't participate in windows)
        // TODO: Replace with C++ 20 https://en.cppreference.com/w/cpp/ranges/filter_view and C++17 std::partial_sort_copy
        // This will have at least one entry since we got at least 1 window.
        vector<uint8_t> worst_qualities;
        worst_qualities.reserve(quality_bytes.size());
        for (size_t i = 0; i < min(sequence.size(), quality_bytes.size()); i++) {
            if (isATGC(sequence[i])) {
                // The base qwuality is not for an N
                worst_qualities.push_back(quality_bytes[i]);
            }
        }
       
        // Make sure the first minimum_disruptive_bases elements are the smallest.
        std::partial_sort(worst_qualities.begin(),
            worst_qualities.begin() + min(worst_qualities.size() - 1, minimum_disruptive_bases - 1), worst_qualities.end());
            
        // Sum them
        size_t total_errored_quality = std::accumulate(worst_qualities.begin(),
            worst_qualities.begin() + min(worst_qualities.size(), minimum_disruptive_bases), (size_t) 0);
        
#ifdef debug
        cerr << "Sum of " << minimum_disruptive_bases << " worst qualities: " << total_errored_quality << endl;
#endif
        
        // We can't be more confident in our mapping than we are that all our
        // located minimizers were not from windows that were created by errors
        // in the read. This can't be <0.
        return total_errored_quality;
    }
}

void MinimizerMapper::map(Alignment& aln, AlignmentEmitter& alignment_emitter) {
    // For each input alignment
    
#ifdef debug
    cerr << "Read " << aln.name() << ": " << aln.sequence() << endl;
#endif

    // Make a new funnel instrumenter to watch us map this read.
    Funnel funnel;
    // Start this alignment 
    funnel.start(aln.name());
   
    if (track_provenance) {
        // Start the minimizer finding stage
        funnel.stage("minimizer");
    }
    
    // This is where we will store the minimizers we find.
    std::vector<Minimizer> minimizers;

    // Find the minimizers and score them.
    double base_target_score = 0.0;
    double base_score = 1.0 + std::log(hard_hit_cap);
    for (size_t i = 0; i < minimizer_indexes.size(); i++) {
        
        // Get minimizers and their window agglomeration starts and lengths
        vector<tuple<gbwtgraph::DefaultMinimizerIndex::minimizer_type, size_t, size_t>> current_minimizers = 
            minimizer_indexes[i]->minimizer_regions(aln.sequence());

        for (auto& m : current_minimizers) {
            double score = 0.0;
            auto hits = minimizer_indexes[i]->count_and_find(get<0>(m));
            if (hits.first > 0) {
                if (hits.first <= hard_hit_cap) {
                    score = base_score - std::log(hits.first);
                } else {
                    score = 1.0;
                }
            }
            minimizers.push_back({ get<0>(m), get<1>(m), get<2>(m), hits.first, hits.second, i, score });
            base_target_score += score;
        }
    }
    double target_score = (base_target_score * minimizer_score_fraction) + 0.000001;
    // Sort (and renumber in best-to-worst order) all the minimizers
    std::sort(minimizers.begin(), minimizers.end());

    if (track_provenance) {
        // Record how many we found, as new lines.
        funnel.introduce(minimizers.size());
        
        // Start the minimizer locating stage
        funnel.stage("seed");
    }

    // Store the seeds and their source minimizers in separate vectors.
    std::vector<pos_t> seeds;
    std::vector<size_t> seed_to_source;

    // Select the minimizers we use for seeds.
    size_t rejected_count = 0;
    double selected_score = 0.0;
    // And count how many hits we had available and we actually located
    size_t total_hits = 0;
    size_t located_hits = 0;
    // And log the hit counts for minimizers we keep and reject
    vector<size_t> used_minimizer_hit_counts;
    used_minimizer_hit_counts.reserve(minimizers.size());
    vector<size_t> unused_minimizer_hit_counts;
    unused_minimizer_hit_counts.reserve(minimizers.size());
    // And flag whether each minimizer in the read was located or not.
    // TODO: should we count minimizers with no hits? We do right now because
    // we need them to be created in the read when coming from a cluster that
    // doesn't have them.
    vector<bool> minimizer_located(minimizers.size(), false);
    // In order to consistently take either all or none of the minimizers in
    // the read with a particular sequence, we track whether we took the
    // previous one.
    bool took_last = false;
    for (size_t i = 0; i < minimizers.size(); i++) {
        if (track_provenance) {
            // Say we're working on it
            funnel.processing_input(i);
        }

        // Select the minimizer if it is informative enough or if the total score
        // of the selected minimizers is not high enough.
        const Minimizer& minimizer = minimizers[i];
        
#ifdef debug
        cerr << "Minimizer " << i << " = " << minimizer.value.key.decode(minimizer_indexes[minimizer.origin]->k())
             << " has " << minimizer.hits << " hits" << endl;
#endif

        // Record that we had this many hits of a minimizer
        total_hits += minimizer.hits;
        
        if (minimizer.hits == 0) {
            // A minimizer with no hits can't go on.
            took_last = false;
            // But we should treat it as located, because we know it isn't anywhere.
            // TODO: should we also include it as needing to be covered by errors for MAPQ capping?
            minimizer_located[i] = true;
            if (track_provenance) {
                funnel.fail("any-hits", i);
            }
        } else if (seeds.size() == 1 ||
            minimizer.hits <= hit_cap || 
            (minimizer.hits <= hard_hit_cap && selected_score + minimizer.score <= target_score) ||
            (took_last && i > 0 && minimizer.value.key == minimizers[i - 1].value.key)) {
            // We should keep this minimizer instance because we only have one
            // seed location already, or it is sufficiently rare, or we want it
            // to make target_score, or it is the same sequence as the previous
            // minimizer which we also took.
        
            // Locate the hits.
            for (size_t j = 0; j < minimizer.hits; j++) {
                pos_t hit = gbwtgraph::DefaultMinimizerIndex::decode(minimizer.occs[j]);
                // Reverse the hits for a reverse minimizer
                if (minimizer.value.is_reverse) {
                    size_t node_length = gbwt_graph.get_length(gbwt_graph.get_handle(id(hit)));
                    hit = reverse_base_pos(hit, node_length);
                }
                // For each position, remember it and what minimizer it came from.
                // Minimizers that occur multiple times will have multiple
                // copies of the same seed locations in the reference, paired
                // with different places in the read.
                seeds.push_back(hit);
                seed_to_source.push_back(i);
            }
            
            if (!(took_last && i > 0 && minimizer.value.key == minimizers[i - 1].value.key)) {
                // We did not also take a previous identical-sequence minimizer, so count this one towards the score.
                selected_score += minimizer.score;
            }
            
            // Remember that we located these hits
            located_hits += minimizer.hits;
            // And that we took this minimizer
            took_last = true;
            
            if (track_provenance) {
                // Record in the funnel that this minimizer gave rise to these seeds.
                funnel.pass("any-hits", i);
                funnel.pass("hard-hit-cap", i);
                funnel.pass("hit-cap||score-fraction", i, selected_score  / base_target_score);
                funnel.expand(i, minimizer.hits);
                
                used_minimizer_hit_counts.push_back(minimizer.hits);
            }
        } else if (minimizer.hits <= hard_hit_cap) {
            // Passed hard hit cap but failed score fraction/normal hit cap
            took_last = false;
            rejected_count++;
            if (track_provenance) {
                funnel.pass("any-hits", i);
                funnel.pass("hard-hit-cap", i);
                funnel.fail("hit-cap||score-fraction", i, (selected_score + minimizer.score) / base_target_score);
                
                unused_minimizer_hit_counts.push_back(minimizer.hits);
            }
        } else {
            // Failed hard hit cap
            took_last = false;  
            rejected_count++;
            if (track_provenance) {
                funnel.pass("any-hits", i);
                funnel.fail("hard-hit-cap", i);
                
                unused_minimizer_hit_counts.push_back(minimizer.hits);
            }
        }
        if (track_provenance) {
            // Say we're done with this input item
            funnel.processed_input();
        }
    }

    if (track_provenance && track_correctness) {
        // Tag seeds with correctness based on proximity along paths to the input read's refpos
        funnel.substage("correct");
      
        if (path_graph == nullptr) {
            cerr << "error[vg::MinimizerMapper] Cannot use track_correctness with no XG index" << endl;
            exit(1);
        }
        
        if (aln.refpos_size() != 0) {
            // Take the first refpos as the true position.
            auto& true_pos = aln.refpos(0);
            
            for (size_t i = 0; i < seeds.size(); i++) {
                // Find every seed's reference positions. This maps from path name to pairs of offset and orientation.
                auto offsets = algorithms::nearest_offsets_in_paths(path_graph, seeds[i], 100);
                for (auto& hit_pos : offsets[path_graph->get_path_handle(true_pos.name())]) {
                    // Look at all the ones on the path the read's true position is on.
                    if (abs((int64_t)hit_pos.first - (int64_t) true_pos.offset()) < 200) {
                        // Call this seed hit close enough to be correct
                        funnel.tag_correct(i);
                    }
                }
            }
        }
    }
        
#ifdef debug
    cerr << "Found " << seeds.size() << " seeds from " << (minimizers.size() - rejected_count) << " minimizers, rejected " << rejected_count << endl;
#endif

    if (track_provenance) {
        // Begin the clustering stage
        funnel.stage("cluster");
    }
        
    // Cluster the seeds. Get sets of input seed indexes that go together.
    vector<vector<size_t>> clusters = clusterer.cluster_seeds(seeds, distance_limit);
    
    if (track_provenance) {
        funnel.substage("score");
    }

    // Cluster score is the sum of minimizer scores.
    std::vector<double> cluster_score(clusters.size(), 0.0);
    vector<double> read_coverage_by_cluster;
    read_coverage_by_cluster.reserve(clusters.size());
    
    // Also track which minimizers from the read participate in each cluster.
    // We need this for scoring the cluster, and also for bounding MAPQ by
    // considering creating subsets of the minimizers.
    vector<vector<bool>> present_in_cluster(clusters.size(), vector<bool>(minimizers.size(), false));

    for (size_t i = 0; i < clusters.size(); i++) {
        // For each cluster
        auto& cluster = clusters[i];
        
        if (track_provenance) {
            // Say we're making it
            funnel.producing_output(i);
        }

        // Which minimizers are present in the cluster.
        for (auto hit_index : cluster) {
            present_in_cluster[i][seed_to_source[hit_index]] = true;
#ifdef debug
            if (present_in_cluster[i][seed_to_source[hit_index]]) {
                cerr << "Minimizer " << seed_to_source[hit_index] << " is present in cluster " << i << endl;
            }
#endif
        }

        // Compute the score.
        for (size_t j = 0; j < minimizers.size(); j++) {
            if (present_in_cluster[i][j]) {
                cluster_score[i] += minimizers[j].score;
            }
        }
        
        if (track_provenance) {
            // Record the cluster in the funnel as a group of the size of the number of items.
            funnel.merge_group(cluster.begin(), cluster.end());
            funnel.score(funnel.latest(), cluster_score[i]);
            
            // Say we made it.
            funnel.produced_output();
        }

        // Get the cluster coverage
        // We set bits in here to true when query anchors cover them
        sdsl::bit_vector covered(aln.sequence().size(), 0);
        for (auto hit_index : cluster) {
            // For each hit in the cluster, work out what anchor sequence it is from.
            const Minimizer& minimizer = minimizers[seed_to_source[hit_index]];

            // The offset of a reverse minimizer is the endpoint of the kmer
            size_t start_offset = minimizer.value.offset;
            size_t k = minimizer_indexes[minimizer.origin]->k();
            if (minimizer.value.is_reverse) {
                start_offset = start_offset + 1 - k;
            }

            // Set the k bits starting at start_offset.
            covered.set_int(start_offset, sdsl::bits::lo_set[k], k);
        }

        // Count up the covered positions
        size_t covered_count = sdsl::util::cnt_one_bits(covered);

        // Turn that into a fraction
        read_coverage_by_cluster.push_back(covered_count / (double) covered.size());
    }

#ifdef debug
    cerr << "Found " << clusters.size() << " clusters" << endl;
#endif

    // Find the best and second-best cluster scores.
    double best_cluster_score = 0;
    double second_best_cluster_score = 0;
    for (auto& score : cluster_score) {
        // For each cluster's score
        if (score > best_cluster_score) {
            // If it is the new best bump down the old best
            second_best_cluster_score = best_cluster_score;
            best_cluster_score = score;
        } else if (score > second_best_cluster_score) {
            // If it only beats the second best, replace that instead.
            second_best_cluster_score = score;
        }
    }
    
    // We will set a score cutoff based on the best, but move it down to the
    // second best if it does not include the second best and the second best
    // is within pad_cluster_score_threshold of where the cutoff would
    // otherwise be. This ensures that we won't throw away all but one cluster
    // based on score alone, unless it is really bad.
                                    
    // Retain clusters only if their score is better than this, in addition to the coverage cutoff
    double cluster_score_cutoff = best_cluster_score - cluster_score_threshold;
    
    if (cluster_score_cutoff - pad_cluster_score_threshold < second_best_cluster_score) {
        // The second best cluster score is high enough that we might want to snap down to it as the cutoff instead.
        cluster_score_cutoff = std::min(cluster_score_cutoff, second_best_cluster_score);
    }
    
    if (track_provenance) {
        // Now we go from clusters to gapless extensions
        funnel.stage("extend");
    }
    
    // These are the GaplessExtensions for all the clusters.
    vector<vector<GaplessExtension>> cluster_extensions;
    cluster_extensions.reserve(clusters.size());
    // These are the clusters the extensions came from, which we need to trace
    // back from alignments to minimizers for annotation.
    // TODO: Do we need those annotations?
    vector<size_t> cluster_extensions_to_source;
    cluster_extensions_to_source.reserve(clusters.size());
    // To compute the windows present in any extended cluster, we need to get
    // all the minimizers in any extended cluster.
    vector<bool> present_in_any_extended_cluster(minimizers.size(), false);
    //For each cluster, what fraction of "equivalent" clusters did we keep?
    vector<double> probability_cluster_lost;
    //What is the score and coverage we are considering and how many reads
    size_t curr_coverage = 0;
    size_t curr_score = 0;
    size_t curr_kept = 0;
    size_t curr_count = 0;
    
    // We track unextended clusters.
    vector<size_t> unextended_clusters;
    unextended_clusters.reserve(clusters.size());
    
    // When a cluster is rejected for extension, call this to record it for MAPQ capping.
    auto cluster_not_extended = [&](size_t cluster_num) {
        // Since we didn't take this cluster, we need to cap the MAPQ with
        // the probability that the read came from it.
        
        // Remember that it is not extended.
        unextended_clusters.push_back(cluster_num);
    };
    
    //Process clusters sorted by both score and read coverage
    process_until_threshold(clusters, read_coverage_by_cluster,
        [&](size_t a, size_t b) {
            if (read_coverage_by_cluster[a] == read_coverage_by_cluster[b]){
                return cluster_score[a] > cluster_score[b];
            } else {
                return read_coverage_by_cluster[a] > read_coverage_by_cluster[b];
            }
        },
        cluster_coverage_threshold, 1, max_extensions,
        [&](size_t cluster_num) {
            // Handle sufficiently good clusters in descending coverage order
            
            if (track_provenance) {
                funnel.pass("cluster-coverage", cluster_num, read_coverage_by_cluster[cluster_num]);
                funnel.pass("max-extensions", cluster_num);
            }
            
            // First check against the additional score filter
            if (cluster_score_threshold != 0 && cluster_score[cluster_num] < cluster_score_cutoff) {
                //If the score isn't good enough, ignore this cluster
                if (track_provenance) {
                    funnel.fail("cluster-score", cluster_num, cluster_score[cluster_num]);
                }
                
                // Record MAPQ implications of not extending this cluster.
                cluster_not_extended(cluster_num);
                return false;
            }
            
            if (track_provenance) {
                funnel.pass("cluster-score", cluster_num, cluster_score[cluster_num]);
                funnel.processing_input(cluster_num);
            }
            if (read_coverage_by_cluster[cluster_num] == curr_coverage &&
                cluster_score[cluster_num] == curr_score &&
                curr_kept < max_extensions * 0.75) {
                curr_kept++;
                curr_count++;
            } else if (read_coverage_by_cluster[cluster_num] != curr_coverage ||
                    cluster_score[cluster_num] != curr_score) {
                    //If this is a cluster that has scores different than the previous one
                    for (size_t i = 0 ; i < curr_kept ; i++ ) {
                        probability_cluster_lost.push_back(1.0 - (double(curr_kept) / double(curr_count)));
                    }
                    curr_coverage = read_coverage_by_cluster[cluster_num];
                    curr_score = cluster_score[cluster_num];
                    curr_kept = 1;
                    curr_count = 1;
            } else {
                //If this cluster is equivalent to the previous one and we already took enough
                //equivalent clusters
                curr_count ++;
                
                // TODO: shouldn't we fail something for the funnel here?
                
                // Record MAPQ implications of not extending this cluster.
                cluster_not_extended(cluster_num);
                return false;
            }
            

            //Only keep this cluster if we have few enough equivalent clusters

            vector<size_t>& cluster = clusters[cluster_num];

#ifdef debug
            cerr << "Cluster " << cluster_num << endl;
            cerr << "Covers " << read_coverage_by_cluster[cluster_num] << "/best-" << cluster_coverage_threshold << " of read" << endl;
            cerr << "Scores " << cluster_score[cluster_num] << "/" << cluster_score_cutoff << endl;
#endif
             
            // Pack the seeds for GaplessExtender.
            GaplessExtender::cluster_type seed_matchings;
            for (auto& seed_index : cluster) {
                // Insert the (graph position, read offset) pair.
                seed_matchings.insert(GaplessExtender::to_seed(seeds[seed_index], minimizers[seed_to_source[seed_index]].value.offset));
#ifdef debug
                const Minimizer& minimizer = minimizers[seed_to_source[seed_index]];
                cerr << "Seed read:" << minimizer.value.offset << " = " << seeds[seed_index]
                    << " from minimizer " << seed_to_source[seed_index] << "(" << minimizer.hits << ")" << endl;
#endif
            }
            
            // Extend seed hits in the cluster into one or more gapless extensions
            cluster_extensions.emplace_back(std::move(extender.extend(seed_matchings, aln.sequence())));
            cluster_extensions_to_source.push_back(cluster_num);
            
            for (size_t i = 0; i < minimizers.size(); i++) {
                // Since the cluster was extended, OR in its minimizers with
                // those in all the other extended clusters
                present_in_any_extended_cluster[i] = (present_in_any_extended_cluster[i] || present_in_cluster[cluster_num][i]);
#ifdef debug
                if (present_in_cluster[cluster_num][i]) {
                    cerr << "Minimizer " << i << " is present in extended cluster " << cluster_num << endl;
                }
#endif
            }
            
            if (track_provenance) {
                // Record with the funnel that the previous group became a group of this size.
                // Don't bother recording the seed to extension matching...
                funnel.project_group(cluster_num, cluster_extensions.back().size());
                
                // Say we finished with this cluster, for now.
                funnel.processed_input();
            }
            return true;
            
        }, [&](size_t cluster_num) {
            // There are too many sufficiently good clusters
            if (track_provenance) {
                funnel.pass("cluster-coverage", cluster_num, read_coverage_by_cluster[cluster_num]);
                funnel.fail("max-extensions", cluster_num);
            }
            if (read_coverage_by_cluster[cluster_num] == curr_coverage &&
                cluster_score[cluster_num] == curr_score) {
                curr_count ++;
            } else {
    
                for (size_t i = 0 ; i < curr_kept ; i++ ) {
                    probability_cluster_lost.push_back(1.0 - (double(curr_kept) / double(curr_count)));
                }
                curr_score = 0;
                curr_coverage = 0;
                curr_kept = 0;
                curr_count = 0;
            }
            
            // Record MAPQ implications of not extending this cluster.
            cluster_not_extended(cluster_num);
            
        }, [&](size_t cluster_num) {
            // This cluster is not sufficiently good.
            if (track_provenance) {
                funnel.fail("cluster-coverage", cluster_num, read_coverage_by_cluster[cluster_num]);
            }
            for (size_t i = 0 ; i < curr_kept ; i++ ) {
                probability_cluster_lost.push_back(1.0 - (double(curr_kept) / double(curr_count)));
            }
            curr_kept = 0;
            curr_count = 0;
            curr_score = 0;
            curr_coverage = 0;
            
            // Record MAPQ implications of not extending this cluster.
            cluster_not_extended(cluster_num);
        });
        
    for (size_t i = 0 ; i < curr_kept ; i++ ) {
        probability_cluster_lost.push_back(1.0 - (double(curr_kept) / double(curr_count)));
    }
    
    if (track_provenance) {
        funnel.substage("score");
    }

    // We now estimate the best possible alignment score for each cluster.
    vector<int> cluster_extension_scores;
    cluster_extension_scores.reserve(cluster_extensions.size());
    for (size_t i = 0; i < cluster_extensions.size(); i++) {
        // For each group of GaplessExtensions
        
        if (track_provenance) {
            funnel.producing_output(i);
        }
        
        auto& extensions = cluster_extensions[i];
        // Work out the score of the best chain of extensions, accounting for gaps/overlaps in the read only
        cluster_extension_scores.push_back(score_extension_group(aln, extensions,
            get_regular_aligner()->gap_open, get_regular_aligner()->gap_extension));
        
        if (track_provenance) {
            // Record the score with the funnel
            funnel.score(i, cluster_extension_scores.back());
            funnel.produced_output();
        }
    }
    
    if (track_provenance) {
        funnel.stage("align");
    }
    
    // Now start the alignment step. Everything has to become an alignment.

    // We will fill this with all computed alignments in estimated score order.
    vector<Alignment> alignments;
    alignments.reserve(cluster_extensions.size());
    // This maps from alignment index back to cluster extension index, for
    // tracing back to minimizers for MAPQ. Can hold
    // numeric_limits<size_t>::max() for an unaligned alignment.
    vector<size_t> alignments_to_source;
    alignments_to_source.reserve(cluster_extensions.size());
    //probability_cluster_lost but ordered by alignment
    vector<double> probability_alignment_lost;
    probability_alignment_lost.reserve(cluster_extensions.size());

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
    
    // Go through the gapless extension groups in score order.
    process_until_threshold(cluster_extensions, cluster_extension_scores,
        extension_set_score_threshold, 2, max_alignments,
        [&](size_t extension_num) {
            // This extension set is good enough.
            // Called in descending score order.
            
            if (track_provenance) {
                funnel.pass("extension-set", extension_num, cluster_extension_scores[extension_num]);
                funnel.pass("max-alignments", extension_num);
                funnel.processing_input(extension_num);
            }
            
            auto& extensions = cluster_extensions[extension_num];
            
            // Get an Alignments best_ and second_best_alignment of it somehow, and throw it in.
            Alignment best_alignment = aln;
            Alignment second_best_alignment = aln;
            
            if (extensions[0].full()) {
                // We got full-length extensions, so directly convert to an Alignment.
                
                if (track_provenance) {
                    funnel.substage("direct");
                }

                //Fill in the best alignments from the extension
                
                *best_alignment.mutable_path() = extensions.front().to_path(gbwt_graph, best_alignment.sequence());
                size_t mismatch_count = extensions.front().mismatches();
                double identity = best_alignment.sequence().size() == 0 ? 0.0 : (best_alignment.sequence().size() - mismatch_count) / (double) best_alignment.sequence().size();
                
                // Fill in the score and identity
                best_alignment.set_score(extensions.front().score);
                best_alignment.set_identity(identity);

                if (extensions.size() > 1) {
                    //Do the same thing for the second extension, if one exists
                    *second_best_alignment.mutable_path() = extensions.back().to_path(gbwt_graph, second_best_alignment.sequence());
                    size_t mismatch_count = extensions.back().mismatches();
                    double identity = second_best_alignment.sequence().size() == 0 ? 0.0 : (second_best_alignment.sequence().size() - mismatch_count) / (double) second_best_alignment.sequence().size();
                    
                    // Fill in the score and identity
                    second_best_alignment.set_score(extensions.back().score);
                    second_best_alignment.set_identity(identity);
                }

                if (track_provenance) {
                    // Stop the current substage
                    funnel.substage_stop();
                }
            } else if (do_dp) {
                // We need to do chaining.
                
                if (track_provenance) {
                    funnel.substage("chain");
                }
                
                // Do the DP and compute alignment into best_alignment and
                // second_best_extension, if there is a second best
                find_optimal_tail_alignments(aln, extensions, best_alignment, second_best_alignment);

                
                if (track_provenance) {
                    // We're done chaining. Next alignment may not go through this substage.
                    funnel.substage_stop();
                }
            } else {
                // We would do chaining but it is disabled.
                // Leave best_alignment unaligned
            }
            
            
            
            if (second_best_alignment.score() != 0 && 
                second_best_alignment.score() > best_alignment.score() * 0.8) {
                //If there is a second extension and its score is at least half of the best score, bring it along
                alignments.emplace_back(std::move(second_best_alignment));
                alignments_to_source.push_back(extension_num);
                probability_alignment_lost.push_back(probability_cluster_lost[extension_num]);

                if (track_provenance) {
    
                    funnel.project(extension_num);
                    funnel.score(alignments.size() - 1, alignments.back().score());
                    // We're done with this input item
                    funnel.processed_input();
                }
            }

            alignments.push_back(std::move(best_alignment));
            alignments_to_source.push_back(extension_num);
            probability_alignment_lost.push_back(probability_cluster_lost[extension_num]);

            if (track_provenance) {

                funnel.project(extension_num);
                funnel.score(alignments.size() - 1, alignments.back().score());
                
                // We're done with this input item
                funnel.processed_input();
            }
            
            return true;
        }, [&](size_t extension_num) {
            // There are too many sufficiently good extensions
            if (track_provenance) {
                funnel.pass("extension-set", extension_num, cluster_extension_scores[extension_num]);
                funnel.fail("max-alignments", extension_num);
            }
        }, [&](size_t extension_num) {
            // This extension is not good enough.
            if (track_provenance) {
                funnel.fail("extension-set", extension_num, cluster_extension_scores[extension_num]);
            }
        });
    
    if (alignments.size() == 0) {
        // Produce an unaligned Alignment
        alignments.emplace_back(aln);
        alignments_to_source.push_back(numeric_limits<size_t>::max());
        probability_alignment_lost.push_back(0);
        
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
    // Track which Alignments they are
    vector<size_t> mappings_to_source;
    mappings_to_source.reserve(min(alignments.size(), max_multimaps));
    
    // Grab all the scores in order for MAPQ computation.
    vector<double> scores;
    scores.reserve(alignments.size());
    
    vector<double> probability_mapping_lost;
    process_until_threshold(alignments, (std::function<double(size_t)>) [&](size_t i) -> double {
        return alignments.at(i).score();
    }, 0, 1, max_multimaps, [&](size_t alignment_num) {
        // This alignment makes it
        // Called in score order
        
        // Remember the score at its rank
        scores.emplace_back(alignments[alignment_num].score());
        
        // Remember the output alignment
        mappings.emplace_back(std::move(alignments[alignment_num]));
        mappings_to_source.push_back(alignment_num);
        probability_mapping_lost.push_back(probability_alignment_lost[alignment_num]);
        
        if (track_provenance) {
            // Tell the funnel
            funnel.pass("max-multimaps", alignment_num);
            funnel.project(alignment_num);
            funnel.score(alignment_num, scores.back());
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
    
#ifdef debug
    cerr << "For scores ";
    for (auto& score : scores) cerr << score << " ";
#endif

    size_t winning_index;
    // Compute MAPQ if not unmapped. Otherwise use 0 instead of the 50% this would give us.
    double mapq = (mappings.empty() || mappings.front().path().mapping_size() == 0) ? 0 : 
        get_regular_aligner()->maximum_mapping_quality_exact(scores, &winning_index) / 2;
        
#ifdef debug
    cerr << "uncapped MAPQ is " << mapq << endl;
#endif

    if (probability_mapping_lost.front() > 0) {
        mapq = min(mapq,round(prob_to_phred(probability_mapping_lost.front())));
    }
    
    // We want to have a MAPQ cap based on minimum error bases it would take to
    // have created all our windows we located in the read we have.
#ifdef debug
    cerr << "Cap based on located minimizers all being faked by errors..." << endl;
#endif

    // Convert our flag vector to a list of the minimizers actually located
    vector<size_t> located_minimizers;
    located_minimizers.reserve(minimizers.size());
    for (size_t i = 0; i < minimizers.size(); i++) {
        if (minimizer_located[i]) {
            located_minimizers.push_back(i);
        }
    }
    double mapq_locate_cap = window_breaking_quality(minimizers, located_minimizers, aln.sequence(), aln.quality());
    
    // We need to cap MAPQ based on the likelihood of generating all the windows in the extended clusters by chance, too.
#ifdef debug
    cerr << "Cap based on extended clusters' minimizers all being faked by errors..." << endl;
#endif

    // Convert our flag vector to a list of the minimizers actually in extended clusters
    vector<size_t> extended_cluster_minimizers;
    extended_cluster_minimizers.reserve(minimizers.size());
    for (size_t i = 0; i < minimizers.size(); i++) {
        if (present_in_any_extended_cluster[i]) {
            extended_cluster_minimizers.push_back(i);
        }
    }
    double mapq_extended_cap = window_breaking_quality(minimizers, extended_cluster_minimizers, aln.sequence(), aln.quality());
    
    // And we also need to cap based on the probability of creating the windows
    // in the read that distinguish it from the most plausible (minimum created
    // windows needed) non-extended cluster.
#ifdef debug
    cerr << "Cap based on read's minimizers not in non-extended clusters all being wrong (and the read actually having come from the non-extended clusters)..." << endl;
#endif

    double mapq_non_extended_cap = numeric_limits<double>::infinity();
    for (auto& cluster_num : unextended_clusters) {
        // For each unextended cluster that might have created the read
        
        // Collect the minimizers not in it but in the read
        vector<size_t> synthesized_minimizers;
        synthesized_minimizers.reserve(minimizers.size());
        for (size_t i = 0; i < minimizers.size(); i++) {
            if (!present_in_cluster[cluster_num][i] && minimizer_located[i]) {
                synthesized_minimizers.push_back(i);
            }
        }
        
        // Cap MAPQ with MAPQ for creating all these minimizers.
        // TODO: only run the capping once if we can figure out the easiest-to-create set of minimizers in advance.
        // TODO: Find the set of minimizers *only* in unextended clusters and use that as the set we have to create.
        mapq_non_extended_cap = std::min(mapq_non_extended_cap,
            window_breaking_quality(minimizers, synthesized_minimizers, aln.sequence(), aln.quality()));
    }
    
    // Remember the uncapped MAPQ and the caps
    set_annotation(mappings[0], "mapq_uncapped", mapq);
    set_annotation(mappings[0], "mapq_locate_cap", mapq_locate_cap);
    set_annotation(mappings[0], "mapq_extended_cap", mapq_extended_cap);
    set_annotation(mappings[0], "mapq_non_extended_cap", mapq_non_extended_cap);
    
    // Apply the caps
    mapq = min(min(mapq, mapq_locate_cap), min(mapq_extended_cap, mapq_non_extended_cap));
    
#ifdef debug
    cerr << "MAPQ is " << mapq << endl;
#endif

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
    
    if (track_provenance) {
    
        // Annotate with statistics we may want to use to compute MAPQ
        
        // Some of them are medians, so we define median of a *sorted* vector.
        auto median = [](const std::vector<size_t>& nums) -> size_t {
            return nums.size() ? nums[nums.size() / 2] : 0;
        };
        
        set_annotation(mappings[0], "mapq_index_hits_total", (double)total_hits);
        set_annotation(mappings[0], "mapq_index_hits_located", (double)located_hits);
        set_annotation(mappings[0], "mapq_index_hits_fraction", (double)located_hits/max(1.0, (double)total_hits));
        set_annotation(mappings[0], "mapq_cluster_equivalent_extended_fraction", 1.0 - probability_mapping_lost.front());
        
        // Get all the minimizer hit counts for minimizers in the cluster that
        // contributed to the winning mapping, if any.
        vector<size_t> cluster_minimizer_hit_counts;
        // Look up the extension number. TODO: query the funnel?
        size_t alignment_num = mappings_to_source.at(0);
        size_t extension_num = alignments_to_source.at(alignment_num);
        if (extension_num != numeric_limits<size_t>::max()) {
            // The mapping came from an extension and is not unmapped. So find the cluster.
            size_t cluster_num = cluster_extensions_to_source.at(extension_num);
            const vector<size_t>& cluster = clusters.at(cluster_num);
            
            // Which minimizers are present in the cluster?
            // Deduplicate with a bit set.
            // TODO: don't recompute this
            vector<bool> present(minimizers.size(), false);
            for (auto hit_index : cluster) {
                present[seed_to_source[hit_index]] = true;
            }

            for (size_t minimizer_num = 0; minimizer_num < present.size(); minimizer_num++) {
                if (present[minimizer_num]) {
                    // For each minimizer in the cluster, record its hit count
                    cluster_minimizer_hit_counts.push_back(minimizers.at(minimizer_num).hits);
                }
            }
            
            std::sort(cluster_minimizer_hit_counts.begin(), cluster_minimizer_hit_counts.end());
        }
        // Annotate each read with the minimizer counts from the winning cluster
        set_annotation(mappings[0], "mapq_winning_cluster_minimizer_counts", cluster_minimizer_hit_counts);
        // And the median count, or 0 if not present.
        set_annotation(mappings[0], "mapq_winning_cluster_median_minimizer_count", median(cluster_minimizer_hit_counts));
       
        // Also annotate with MAPQ-relevant facts about the read and not just the cluster
        vector<size_t> minimizer_hit_counts;
        minimizer_hit_counts.reserve(minimizers.size());
        for (auto& minimizer : minimizers) {
            minimizer_hit_counts.push_back(minimizer.hits);
        }
        std::sort(minimizer_hit_counts.begin(), minimizer_hit_counts.end());
        set_annotation(mappings[0], "mapq_minimizer_counts", minimizer_hit_counts);
        set_annotation(mappings[0], "mapq_median_minimizer_count", median(minimizer_hit_counts));
        std::sort(used_minimizer_hit_counts.begin(), used_minimizer_hit_counts.end());
        set_annotation(mappings[0], "mapq_used_minimizer_counts", used_minimizer_hit_counts);
        set_annotation(mappings[0], "mapq_median_used_minimizer_count", median(used_minimizer_hit_counts));
        std::sort(unused_minimizer_hit_counts.begin(), unused_minimizer_hit_counts.end());
        set_annotation(mappings[0], "mapq_unused_minimizer_counts", unused_minimizer_hit_counts);
        set_annotation(mappings[0], "mapq_median_unused_minimizer_count", median(unused_minimizer_hit_counts));
        
        // Annotate with the number of results in play at each stage
        funnel.for_each_stage([&](const string& stage, const vector<size_t>& result_sizes) {
            // Save the number of items
            set_annotation(mappings[0], "stage_" + stage + "_results", (double)result_sizes.size());
        });
        
        if (track_correctness) {
            // And with the last stage at which we had any descendants of the correct seed hit locations
            set_annotation(mappings[0], "last_correct_stage", funnel.last_correct_stage());
        }
        
        // Annotate with the performances of all the filters
        // We need to track filter number
        size_t filter_num = 0;
        funnel.for_each_filter([&](const string& stage, const string& filter,
            const Funnel::FilterPerformance& by_count, const Funnel::FilterPerformance& by_size,
            const vector<double>& filter_statistics_correct, const vector<double>& filter_statistics_non_correct) {
            
            string filter_id = to_string(filter_num) + "_" + filter + "_" + stage;
            
            // Save the stats
            set_annotation(mappings[0], "filter_" + filter_id + "_passed_count_total", (double) by_count.passing);
            set_annotation(mappings[0], "filter_" + filter_id + "_failed_count_total", (double) by_count.failing);
            
            set_annotation(mappings[0], "filter_" + filter_id + "_passed_size_total", (double) by_size.passing);
            set_annotation(mappings[0], "filter_" + filter_id + "_failed_size_total", (double) by_size.failing);
            
            if (track_correctness) {
                set_annotation(mappings[0], "filter_" + filter_id + "_passed_count_correct", (double) by_count.passing_correct);
                set_annotation(mappings[0], "filter_" + filter_id + "_failed_count_correct", (double) by_count.failing_correct);
                
                set_annotation(mappings[0], "filter_" + filter_id + "_passed_size_correct", (double) by_size.passing_correct);
                set_annotation(mappings[0], "filter_" + filter_id + "_failed_size_correct", (double) by_size.failing_correct);
            }
            
            // Save the correct and non-correct filter statistics, even if
            // everything is non-correct because correctness isn't computed
            set_annotation(mappings[0], "filterstats_" + filter_id + "_correct", filter_statistics_correct);
            set_annotation(mappings[0], "filterstats_" + filter_id + "_noncorrect", filter_statistics_non_correct);
            
            filter_num++;
        });
        
        // Annotate with parameters used for the filters.
        set_annotation(mappings[0], "param_hit-cap", (double) hit_cap);
        set_annotation(mappings[0], "param_hard-hit-cap", (double) hard_hit_cap);
        set_annotation(mappings[0], "param_score-fraction", (double) minimizer_score_fraction);
        set_annotation(mappings[0], "param_max-extensions", (double) max_extensions);
        set_annotation(mappings[0], "param_max-alignments", (double) max_alignments);
        set_annotation(mappings[0], "param_cluster-score", (double) cluster_score_threshold);
        set_annotation(mappings[0], "param_pad-cluster-score", (double) pad_cluster_score_threshold);
        set_annotation(mappings[0], "param_cluster-coverage", (double) cluster_coverage_threshold);
        set_annotation(mappings[0], "param_extension-set", (double) extension_set_score_threshold);
        set_annotation(mappings[0], "param_max-multimaps", (double) max_multimaps);
    }
    
    // Ship out all the aligned alignments
    alignment_emitter.emit_mapped_single(std::move(mappings));

#undef debug

#ifdef debug
    // Dump the funnel info graph.
    funnel.to_dot(cerr);
#endif
}

int MinimizerMapper::score_extension_group(const Alignment& aln, const vector<GaplessExtension>& extended_seeds,
    int gap_open_penalty, int gap_extend_penalty) {
        
    if (extended_seeds.empty()) {
        // TODO: We should never see an empty group of extensions
        return 0;
    } else if (extended_seeds.front().full()) {
        // These are length matches. We already have the score.
        int best_score = 0;
        for (auto& extension : extended_seeds) {
            best_score = max(best_score, extension.score);
        }
        return best_score;
    } else {
        // This is a collection of one or more non-full-length extended seeds.
        
        if (aln.sequence().size() == 0) {
            // No score here
            return 0;
        }
       
        // We use a sweep line algorithm to find relevant points along the read: extension starts or ends.
        // This records the last base to be covered by the current sweep line.
        int64_t sweep_line = 0;
        // This records the first base not covered by the last sweep line.
        int64_t last_sweep_line = 0;
        
        // And we track the next unentered gapless extension
        size_t unentered = 0;
        
        // Extensions we are in are in this min-heap of past-end position and gapless extension number.
        vector<pair<size_t, size_t>> end_heap;
        // The heap uses this comparator
        auto min_heap_on_first = [](const pair<size_t, size_t>& a, const pair<size_t, size_t>& b) {
            // Return true if a must come later in the heap than b
            return a.first > b.first;
        };
        
        // We track the best score for a chain reaching the position before this one and ending in a gap.
        // We never let it go below 0.
        // Will be 0 when there's no gap that can be open
        int best_gap_score = 0;
        
        // We track the score for the best chain ending with each gapless extension
        vector<int> best_chain_score(extended_seeds.size(), 0);
        
        // And we're after the best score overall that we can reach when an extension ends
        int best_past_ending_score_ever = 0;
        
        // Overlaps are more complicated.
        // We need a heap of all the extensions for which we have seen the
        // start and that we can thus overlap.
        // We filter things at the top of the heap if their past-end positions
        // have occurred.
        // So we store pairs of score we get backtracking to the current
        // position, and past-end position for the thing we are backtracking
        // from.
        vector<pair<int, size_t>> overlap_heap;
        // We can just use the standard max-heap comparator
        
        // We encode the score relative to a counter that we increase by the
        // gap extend every base we go through, so we don't need to update and
        // re-sort the heap.
        int overlap_score_offset = 0;
        
        while(last_sweep_line <= aln.sequence().size()) {
            // We are processed through the position before last_sweep_line.
            
            // Find a place for sweep_line to go
            
            // Find the next seed start
            int64_t next_seed_start = numeric_limits<int64_t>::max();
            if (unentered < extended_seeds.size()) {
                next_seed_start = extended_seeds[unentered].read_interval.first;
            }
            
            // Find the next seed end
            int64_t next_seed_end = numeric_limits<int64_t>::max();
            if (!end_heap.empty()) {
                next_seed_end = end_heap.front().first;
            }
            
            // Whichever is closer between those points and the end, do that.
            sweep_line = min(min(next_seed_end, next_seed_start), (int64_t) aln.sequence().size());
            
            // So now we're only interested in things that happen at sweep_line.
            
            // Compute the distance from the previous sweep line position
            // Make sure to account for last_sweep_line's semantics as the next unswept base.
            int sweep_distance = sweep_line - last_sweep_line + 1;
            
            // We need to track the score of the best thing that past-ended here
            int best_past_ending_score_here = 0;
            
            while(!end_heap.empty() && end_heap.front().first == sweep_line) {
                // Find anything that past-ends here
                size_t past_ending = end_heap.front().second;
                
                // Mix it into the score
                best_past_ending_score_here = std::max(best_past_ending_score_here, best_chain_score[past_ending]);
                
                // Remove it from the end-tracking heap
                std::pop_heap(end_heap.begin(), end_heap.end(), min_heap_on_first);
                end_heap.pop_back();
            }
            

            // Mix that into the best score overall
            best_past_ending_score_ever = std::max(best_past_ending_score_ever, best_past_ending_score_here);
            
            if (sweep_line == aln.sequence().size()) {
                // We don't need to think about gaps or backtracking anymore since everything has ended
                break;
            }
            
            // Update the overlap score offset by removing some gap extends from it.
            overlap_score_offset += sweep_distance * gap_extend_penalty;
            
            // The best way to backtrack to here is whatever is on top of the heap, if anything, that doesn't past-end here.
            int best_overlap_score = 0;
            while (!overlap_heap.empty()) {
                // While there is stuff on the heap
                if (overlap_heap.front().second <= sweep_line) {
                    // We are already past this thing, so drop it
                    std::pop_heap(overlap_heap.begin(), overlap_heap.end());
                    overlap_heap.pop_back();
                } else {
                    // This is at the top of the heap and we aren't past it
                    // Decode and use its score offset if we only backtrack to here.
                    best_overlap_score = overlap_heap.front().first + overlap_score_offset;
                    // Stop looking in the heap
                    break;
                }
            }
            
            // The best way to end 1 before here in a gap is either:
            
            if (best_gap_score != 0) {
                // Best way to end 1 before our last sweep line position with a gap, plus distance times gap extend penalty
                best_gap_score -= sweep_distance * gap_extend_penalty;
            }
            
            // Best way to end 1 before here with an actual extension, plus the gap open part of the gap open penalty.
            // (Will never be taken over an actual adjacency)
            best_gap_score = std::max(0, std::max(best_gap_score, best_past_ending_score_here - (gap_open_penalty - gap_extend_penalty)));
            
            while (unentered < extended_seeds.size() && extended_seeds[unentered].read_interval.first == sweep_line) {
                // For each thing that starts here
                
                // Compute its chain score
                best_chain_score[unentered] = std::max(best_overlap_score,
                    std::max(best_gap_score, best_past_ending_score_here)) + extended_seeds[unentered].score;
                
                // Compute its backtrack-to-here score and add it to the backtracking heap
                // We want how far we would have had to have backtracked to be
                // able to preceed the base we are at now, where this thing
                // starts.
                size_t extension_length = extended_seeds[unentered].read_interval.second - extended_seeds[unentered].read_interval.first;
                int raw_overlap_score = best_chain_score[unentered] - gap_open_penalty - gap_extend_penalty * extension_length;
                int encoded_overlap_score = raw_overlap_score - overlap_score_offset;
                
                // Stick it in the heap
                overlap_heap.emplace_back(encoded_overlap_score, extended_seeds[unentered].read_interval.second);
                std::push_heap(overlap_heap.begin(), overlap_heap.end());
                
                // Add it to the end finding heap
                end_heap.emplace_back(extended_seeds[unentered].read_interval.second, unentered);
                std::push_heap(end_heap.begin(), end_heap.end(), min_heap_on_first);
                
                // Advance and check the next thing to start
                unentered++;
            }
            
            // Move last_sweep_line to sweep_line.
            // We need to add 1 since last_sweep_line is the next *un*included base
            last_sweep_line = sweep_line + 1;
        }
        

        // When we get here, we've seen the end of every extension and so we
        // have the best score at the end of any of them.
        return best_past_ending_score_ever;
    }


}

void MinimizerMapper::find_optimal_tail_alignments(const Alignment& aln, const vector<GaplessExtension>& extended_seeds, Alignment& best, Alignment& second_best) const {

#ifdef debug
    cerr << "Trying to find tail alignments for " << extended_seeds.size() << " extended seeds" << endl;
#endif

    // Make paths for all the extensions
    vector<Path> extension_paths;
    vector<double> extension_path_scores;
    extension_paths.reserve(extended_seeds.size());
    extension_path_scores.reserve(extended_seeds.size());
    for (auto& extended_seed : extended_seeds) {
        // Compute the path for each extension
        extension_paths.push_back(extended_seed.to_path(gbwt_graph, aln.sequence()));
        // And the extension's score
        extension_path_scores.push_back(get_regular_aligner()->score_partial_alignment(aln, gbwt_graph, extension_paths.back(),
            aln.sequence().begin() + extended_seed.read_interval.first));
    }
    
    // We will keep the winning alignment here, in pieces
    Path winning_left;
    Path winning_middle;
    Path winning_right;
    size_t winning_score = 0;

    Path second_left;
    Path second_middle;
    Path second_right;
    size_t second_score = 0;
    
    // Handle each extension in the set
    process_until_threshold(extended_seeds, extension_path_scores,
        extension_score_threshold, 1, max_local_extensions,
        (function<double(size_t)>) [&](size_t extended_seed_num) {
       
            // This extended seed looks good enough.
            
            // TODO: We don't track this filter with the funnel because it
            // operates within a single "item" (i.e. cluster/extension set).
            // We track provenance at the item level, so throwing out wrong
            // local alignments in a correct cluster would look like throwing
            // out correct things.
            // TODO: Revise how we track correctness and provenance to follow
            // sub-cluster things.
       
            // We start with the path in extension_paths[extended_seed_num],
            // scored in extension_path_scores[extended_seed_num]
            
            // We also have a left tail path and score
            pair<Path, int64_t> left_tail_result {{}, 0};
            // And a right tail path and score
            pair<Path, int64_t> right_tail_result {{}, 0};
           
            if (extended_seeds[extended_seed_num].read_interval.first != 0) {
                // There is a left tail
    
                // Get the forest of all left tail placements
                auto forest = get_tail_forest(extended_seeds[extended_seed_num], aln.sequence().size(), true);
           
                // Grab the part of the read sequence that comes before the extension
                string before_sequence = aln.sequence().substr(0, extended_seeds[extended_seed_num].read_interval.first);
                
                // Do right-pinned alignment
                left_tail_result = std::move(get_best_alignment_against_any_tree(forest, before_sequence,
                    extended_seeds[extended_seed_num].starting_position(gbwt_graph), false));
            }
            
            if (extended_seeds[extended_seed_num].read_interval.second != aln.sequence().size()) {
                // There is a right tail
                
                // Get the forest of all right tail placements
                auto forest = get_tail_forest(extended_seeds[extended_seed_num], aln.sequence().size(), false);
            
                // Find the sequence
                string trailing_sequence = aln.sequence().substr(extended_seeds[extended_seed_num].read_interval.second);
        
                // Do left-pinned alignment
                right_tail_result = std::move(get_best_alignment_against_any_tree(forest, trailing_sequence,
                    extended_seeds[extended_seed_num].tail_position(gbwt_graph), true));
            }

            // Compute total score
            size_t total_score = extension_path_scores[extended_seed_num] + left_tail_result.second + right_tail_result.second;

            //Get the node ids of the beginning and end of each alignment

            id_t winning_start = winning_score == 0 ? 0 : (winning_left.mapping_size() == 0
                                          ? winning_middle.mapping(0).position().node_id()
                                          : winning_left.mapping(0).position().node_id());
            id_t current_start = left_tail_result.first.mapping_size() == 0
                                     ? extension_paths[extended_seed_num].mapping(0).position().node_id()
                                     : left_tail_result.first.mapping(0).position().node_id();
            id_t winning_end = winning_score == 0 ? 0 : (winning_right.mapping_size() == 0
                                  ? winning_middle.mapping(winning_middle.mapping_size() - 1).position().node_id()
                                  : winning_right.mapping(winning_right.mapping_size()-1).position().node_id());
            id_t current_end = right_tail_result.first.mapping_size() == 0
                                ? extension_paths[extended_seed_num].mapping(extension_paths[extended_seed_num].mapping_size() - 1).position().node_id()
                                : right_tail_result.first.mapping(right_tail_result.first.mapping_size()-1).position().node_id();
            //Is this left tail different from the currently winning left tail?
            bool different_left = winning_start != current_start;
            bool different_right = winning_end != current_end;


            if (total_score > winning_score || winning_score == 0) {
                // This is the new best alignment seen so far.

                
                if (winning_score != 0 && different_left && different_right) {
                //The previous best scoring alignment replaces the second best
                    second_score = winning_score;
                    second_left = std::move(winning_left);
                    second_middle = std::move(winning_middle);
                    second_right = std::move(winning_right);
                }

                // Save the score
                winning_score = total_score;
                // And the path parts
                winning_left = std::move(left_tail_result.first);
                winning_middle = std::move(extension_paths[extended_seed_num]);
                winning_right = std::move(right_tail_result.first);

            } else if ((total_score > second_score || second_score == 0) && different_left && different_right) {
                // This is the new second best alignment seen so far and it is 
                // different from the best alignment.
                
                // Save the score
                second_score = total_score;
                // And the path parts
                second_left = std::move(left_tail_result.first);
                second_middle = std::move(extension_paths[extended_seed_num]);
                second_right = std::move(right_tail_result.first);
            }

            return true;
        }, [&](size_t extended_seed_num) {
            // This extended seed is good enough by its own score, but we have too many.
            // Do nothing
        }, [&](size_t extended_seed_num) {
            // This extended seed isn't good enough by its own score.
            // Do nothing
        });
        
    // Now we know the winning path and score. Move them over to out
    best.set_score(winning_score);
    second_best.set_score(second_score);

    // Concatenate the paths. We know there must be at least an edit boundary
    // between each part, because the maximal extension doesn't end in a
    // mismatch or indel and eats all matches.
    // We also don't need to worry about jumps that skip intervening sequence.
    *best.mutable_path() = std::move(winning_left);

    for (auto* to_append : {&winning_middle, &winning_right}) {
        // For each path to append
        for (auto& mapping : *to_append->mutable_mapping()) {
            // For each mapping to append
            
            if (mapping.position().offset() != 0 && best.path().mapping_size() > 0) {
                // If we have a nonzero offset in our mapping, and we follow
                // something, we must be continuing on from a previous mapping to
                // the node.
                assert(mapping.position().node_id() == best.path().mapping(best.path().mapping_size() - 1).position().node_id());

                // Find that previous mapping
                auto* prev_mapping = best.mutable_path()->mutable_mapping(best.path().mapping_size() - 1);
                for (auto& edit : *mapping.mutable_edit()) {
                    // Move over all the edits in this mapping onto the end of that one.
                    *prev_mapping->add_edit() = std::move(edit);
                }
            } else {
                // If we start at offset 0 or there's nothing before us, we need to just move the whole mapping
                *best.mutable_path()->add_mapping() = std::move(mapping);
            }
        }
    }
    best.set_identity(identity(best.path()));
    //Do the same for the second best
    *second_best.mutable_path() = std::move(second_left);

    for (auto* to_append : {&second_middle, &second_right}) {
        // For each path to append
        for (auto& mapping : *to_append->mutable_mapping()) {
            // For each mapping to append
            
            if (mapping.position().offset() != 0 && second_best.path().mapping_size() > 0) {
                // If we have a nonzero offset in our mapping, and we follow
                // something, we must be continuing on from a previous mapping to
                // the node.
                assert(mapping.position().node_id() == second_best.path().mapping(second_best.path().mapping_size() - 1).position().node_id());

                // Find that previous mapping
                auto* prev_mapping = second_best.mutable_path()->mutable_mapping(second_best.path().mapping_size() - 1);
                for (auto& edit : *mapping.mutable_edit()) {
                    // Move over all the edits in this mapping onto the end of that one.
                    *prev_mapping->add_edit() = std::move(edit);
                }
            } else {
                // If we start at offset 0 or there's nothing before us, we need to just move the whole mapping
                *second_best.mutable_path()->add_mapping() = std::move(mapping);
            }
        }
    }

    // Compute the identity from the path.
    second_best.set_identity(identity(second_best.path()));
}

pair<Path, size_t> MinimizerMapper::get_best_alignment_against_any_tree(const vector<TreeSubgraph>& trees,
    const string& sequence, const Position& default_position, bool pin_left) const {
   
    // We want the best alignment, to the base graph, done against any target path
    Path best_path;
    // And its score
    int64_t best_score = 0;
    
    if (!sequence.empty()) {
        // We start out with the best alignment being a pure softclip.
        // If we don't have any trees, or all trees are empty, or there's nothing beter, this is what we return.
        Mapping* m = best_path.add_mapping();
        Edit* e = m->add_edit();
        e->set_from_length(0);
        e->set_to_length(sequence.size());
        e->set_sequence(sequence);
        // Since the softclip consumes no graph, we place it on the node we are going to.
        *m->mutable_position() = default_position;
        
#ifdef debug
        cerr << "First best alignment: " << pb2json(best_path) << " score " << best_score << endl;
#endif
    }
    
    // We can align it once per target tree
    for (auto& subgraph : trees) {
        // For each tree we can map against, map pinning the correct edge of the sequence to the root.
        
        if (subgraph.get_node_count() != 0) {
            // This path has bases in it and could potentially be better than
            // the default full-length softclip

            // Do alignment to the subgraph with GSSWAligner.
            Alignment current_alignment;
            // If pinning right, we need to reverse the sequence, since we are
            // always pinning left to the left edge of the tree subgraph.
            current_alignment.set_sequence(pin_left ? sequence : reverse_complement(sequence));
#ifdef debug
            cerr << "Align " << pb2json(current_alignment) << " pinned left";

#ifdef debug_dump_graph
            cerr << " vs graph:" << endl;
            subgraph.for_each_handle([&](const handle_t& here) {
                cerr << subgraph.get_id(here) << " (" << subgraph.get_sequence(here) << "): " << endl;
                subgraph.follow_edges(here, true, [&](const handle_t& there) {
                    cerr << "\t" << subgraph.get_id(there) << " (" << subgraph.get_sequence(there) << ") ->" << endl;
                });
                subgraph.follow_edges(here, false, [&](const handle_t& there) {
                    cerr << "\t-> " << subgraph.get_id(there) << " (" << subgraph.get_sequence(there) << ")" << endl;
                });
            });
#else
            cerr << endl;
#endif
#endif
            
            // Align, accounting for full length bonus.
            // We *always* do left-pinned alignment internally, since that's the shape of trees we get.
            get_regular_aligner()->get_xdrop()->align_pinned(current_alignment, subgraph, subgraph.get_topological_order(), true);
            
#ifdef debug
            cerr << "\tScore: " << current_alignment.score() << endl;
#endif
            
            if (current_alignment.score() > best_score) {
                // This is a new best alignment.
                best_path = current_alignment.path();
                
                if (!pin_left) {
                    // Un-reverse it if we were pinning right
                    best_path = reverse_complement_path(best_path, [&](id_t node) { 
                        return subgraph.get_length(subgraph.get_handle(node, false));
                    });
                }
                
                // Translate from subgraph into base graph and keep it.
                best_path = subgraph.translate_down(best_path);
                best_score = current_alignment.score();
                
#ifdef debug
                cerr << "New best alignment is "
                    << pb2json(best_path) << " score " << best_score << endl;
#endif
            }
        }
    }

    return make_pair(best_path, best_score);
}

vector<TreeSubgraph> MinimizerMapper::get_tail_forest(const GaplessExtension& extended_seed,
    size_t read_length, bool left_tails) const {

    // We will fill this in with all the trees we return
    vector<TreeSubgraph> to_return;

    // Now for this extension, walk the GBWT in the appropriate direction
    
#ifdef debug
    cerr << "Look for " << (left_tails ? "left" : "right") << " tails from extension" << endl;
#endif

    // TODO: Come up with a better way to do this with more accessors on the extension and less get_handle
    // Get the Position reading out of the extension on the appropriate tail
    Position from;
    // And the length of that tail
    size_t tail_length;
    // And the GBWT search state we want to start with
    const gbwt::SearchState* base_state = nullptr;
    if (left_tails) {
        // Look right from start 
        from = extended_seed.starting_position(gbwt_graph);
        // And then flip to look the other way at the prev base
        from = reverse(from, gbwt_graph.get_length(gbwt_graph.get_handle(from.node_id(), false)));
       
        // Use the search state going backward
        base_state = &extended_seed.state.backward;
       
        tail_length = extended_seed.read_interval.first;
    } else {
        // Look right from end
        from = extended_seed.tail_position(gbwt_graph);
        
        // Use the search state going forward
        base_state = &extended_seed.state.forward;
        
        tail_length = read_length - extended_seed.read_interval.second;
    }

    if (tail_length == 0) {
        // Don't go looking for places to put no tail.
        return to_return;
    }

    // This is one tree that we are filling in
    vector<pair<int64_t, handle_t>> tree;
    
    // This is a stack of indexes at which we put parents in the tree
    list<int64_t> parent_stack;
    
    // Get the handle we are starting from
    // TODO: is it cheaper to get this out of base_state? 
    handle_t start_handle = gbwt_graph.get_handle(from.node_id(), from.is_reverse());
    
    // Decide if the start node will end up included in the tree, or if we cut it all off with the offset.
    bool start_included = (from.offset() < gbwt_graph.get_length(start_handle));
    
    // How long should we search? It should be the longest detectable gap plus the remaining sequence.
    size_t search_limit = get_regular_aligner()->longest_detectable_gap(tail_length, read_length) + tail_length;
    
    // Do a DFS over the haplotypes in the GBWT out to that distance.
    dfs_gbwt(*base_state, from.offset(), search_limit, [&](const handle_t& entered) {
        // Enter a new handle.
        
        if (parent_stack.empty()) {
            // This is the root of a new tree in the forrest
            
            if (!tree.empty()) {
                // Save the old tree and start a new one.
                // We need to cut off from.offset() from the root, unless we would cut off the whole root.
                // In that case, the GBWT DFS will have skipped the empty root entirely, so we cut off nothing.
                to_return.emplace_back(&gbwt_graph, std::move(tree), start_included ? from.offset() : 0);
                tree.clear();
            }
            
            // Add this to the tree with no parent
            tree.emplace_back(-1, entered);
        } else {
            // Just say this is visitable from our parent.
            tree.emplace_back(parent_stack.back(), entered);
        }
        
        // Record the parent index
        parent_stack.push_back(tree.size() - 1);
    }, [&]() {
        // Exit the last visited handle. Pop off the stack.
        parent_stack.pop_back();
    });
    
    if (!tree.empty()) {
        // Now save the last tree
        to_return.emplace_back(&gbwt_graph, std::move(tree), start_included ? from.offset() : 0);
        tree.clear();
    }
    
#ifdef debug
    cerr << "Found " << to_return.size() << " trees" << endl;
#endif
    
    // Now we have all the trees!
    return to_return;
}

size_t MinimizerMapper::immutable_path_from_length(const ImmutablePath& path) {
    size_t to_return = 0;
    for (auto& m : path) {
        // Sum up the from lengths of all the component Mappings
        to_return += mapping_from_length(m);
    }
    return to_return;
}

Path MinimizerMapper::to_path(const ImmutablePath& path) {
    Path to_return;
    for (auto& m : path) {
        // Copy all the Mappings into the Path.
        *to_return.add_mapping() = m;
    }
    
    // Flip the order around to actual path order.
    std::reverse(to_return.mutable_mapping()->begin(), to_return.mutable_mapping()->end());
    
    // Return the completed path
    return to_return;
}

void MinimizerMapper::dfs_gbwt(const Position& from, size_t walk_distance,
    const function<void(const handle_t&)>& enter_handle, const function<void(void)> exit_handle) const {
   
    // Get a handle to the node the from position is on, in the position's forward orientation
    handle_t start_handle = gbwt_graph.get_handle(from.node_id(), from.is_reverse());
    
    // Delegate to the handle-based version
    dfs_gbwt(start_handle, from.offset(), walk_distance, enter_handle, exit_handle);
    
}

void MinimizerMapper::dfs_gbwt(handle_t from_handle, size_t from_offset, size_t walk_distance,
    const function<void(const handle_t&)>& enter_handle, const function<void(void)> exit_handle) const {
    
    // Turn from_handle into a SearchState for everything on it.
    gbwt::SearchState start_state = gbwt_graph.get_state(from_handle);
    
    // Delegate to the state-based version
    dfs_gbwt(start_state, from_offset, walk_distance, enter_handle, exit_handle);
}
    
void MinimizerMapper::dfs_gbwt(const gbwt::SearchState& start_state, size_t from_offset, size_t walk_distance,
    const function<void(const handle_t&)>& enter_handle, const function<void(void)> exit_handle) const {
    
    // Holds the gbwt::SearchState we are at, and the distance we have consumed
    using traversal_state_t = pair<gbwt::SearchState, size_t>;
    
    if (start_state.empty()) {
        // No haplotypes even visit the first node. Stop.
        return;
    }
    
    // Get the handle we are starting on
    handle_t from_handle = gbwt_graph.node_to_handle(start_state.node);

    // The search state represents searching through the end of the node, so we have to consume that much search limit.

    // Tack on how much search limit distance we consume by going to the end of
    // the node. Our start position is a cut *between* bases, and we take everything after it.
    // If the cut is at the offset of the whole length of the node, we take 0 bases.
    // If it is at 0, we take all the bases in the node.
    size_t distance_to_node_end = gbwt_graph.get_length(from_handle) - from_offset;
    
#ifdef debug
    cerr << "DFS starting at offset " << from_offset << " on node of length "
        << gbwt_graph.get_length(from_handle) << " leaving " << distance_to_node_end << " bp" << endl;
#endif


    // Have a recursive function that does the DFS. We fire the enter and exit
    // callbacks, and the user can keep their own stack.
    function<void(const gbwt::SearchState&, size_t, bool)> recursive_dfs = [&](const gbwt::SearchState& here_state,
        size_t used_distance, bool hide_root) {
        
        handle_t here_handle = gbwt_graph.node_to_handle(here_state.node);
        
        if (!hide_root) {
            // Enter this handle if there are any bases on it to visit
            
#ifdef debug
            cerr << "Enter handle " << gbwt_graph.get_id(here_handle) << " " << gbwt_graph.get_is_reverse(here_handle) << endl;
#endif
            
            enter_handle(here_handle);
        }
        
        // Up the used distance with our length
        used_distance += gbwt_graph.get_length(here_handle);
        
        if (used_distance < walk_distance) {
            // If we haven't used up all our distance yet
            
            gbwt_graph.follow_paths(here_state, [&](const gbwt::SearchState& there_state) -> bool {
                // For each next state
                
                // Otherwise, do it with the new distance value.
                // Don't hide the root on any child subtrees; only the top root can need hiding.
                recursive_dfs(there_state, used_distance, false);
                
                return true;
            });
        }
            
        if (!hide_root) {
            // Exit this handle if we entered it
            
#ifdef debug
            cerr << "Exit handle " << gbwt_graph.get_id(here_handle) << " " << gbwt_graph.get_is_reverse(here_handle) << endl;
#endif
            
            exit_handle();
        }
    };
    
    // Start the DFS with our stating node, consuming the distance from our
    // offset to its end. Don't show the root state to the user if we don't
    // actually visit any bases on that node.
    recursive_dfs(start_state, distance_to_node_end, distance_to_node_end == 0);

}

}


