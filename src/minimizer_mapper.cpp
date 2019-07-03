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

// Set this to track provenance of intermediate results
//#define TRACK_PROVENANCE
// With TRACK_PROVENANCE on, set this to track correctness, which requires some expensive XG queries
//#define TRACK_CORRECTNESS

namespace vg {

using namespace std;

MinimizerMapper::MinimizerMapper(const XG* xg_index, const gbwt::GBWT* gbwt_index, const MinimizerIndex* minimizer_index,
     MinimumDistanceIndex* distance_index) :
    xg_index(xg_index), gbwt_index(gbwt_index), minimizer_index(minimizer_index),
    distance_index(distance_index), gbwt_graph(*gbwt_index, *xg_index),
    extender(gbwt_graph), clusterer(*distance_index) {
    
    // Nothing to do!
}

void MinimizerMapper::map(Alignment& aln, AlignmentEmitter& alignment_emitter) {
    // For each input alignment

    // Make a new funnel instrumenter to watch us map this read.
    Funnel funnel;
    // Start this alignment 
    funnel.start(aln.name());
    
    // Annotate the original read with metadata
    if (!sample_name.empty()) {
        aln.set_sample_name(sample_name);
    }
    if (!read_group.empty()) {
        aln.set_read_group(read_group);
    }
   
#ifdef TRACK_PROVENANCE
    // Start the minimizer finding stage
    funnel.stage("minimizer");
#endif
    
    // We will find all the seed hits
    vector<pos_t> seeds;
    
    // This will hold all the minimizers in the query
    vector<MinimizerIndex::minimizer_type> minimizers;
    // And either way this will map from seed to minimizer that generated it
    vector<size_t> seed_to_source;
    
    // Find minimizers in the query
    minimizers = minimizer_index->minimizers(aln.sequence());
    
#ifdef TRACK_PROVENANCE
    // Record how many we found, as new lines.
    funnel.introduce(minimizers.size());
    
    // Start the minimizer locating stage
    funnel.stage("seed");
#endif

    // Compute minimizer scores for all minimizers as 1 + ln(hard_hit_cap) - ln(hits).
    std::vector<double> minimizer_score(minimizers.size(), 0.0);
    double target_score = 0.0;
    for (size_t i = 0; i < minimizers.size(); i++) {
        size_t hits = minimizer_index->count(minimizers[i]);
        if (hits > 0) {
            if (hits <= hard_hit_cap) {
                minimizer_score[i] = 1.0 + std::log(hard_hit_cap) - std::log(hits);
            } else {
                minimizer_score[i] = 1.0;
            }
        }
        target_score += minimizer_score[i];
    }
    target_score *= minimizer_score_fraction;

    // Sort the minimizers by score.
    std::vector<size_t> minimizers_in_order(minimizers.size());
    for (size_t i = 0; i < minimizers_in_order.size(); i++) {
        minimizers_in_order[i] = i;
    }
    std::sort(minimizers_in_order.begin(), minimizers_in_order.end(), [&minimizer_score](const size_t a, const size_t b) {
        return (minimizer_score[a] > minimizer_score[b]);
    });

    // Select the minimizers we use for seeds.
    size_t rejected_count = 0;
    double selected_score = 0.0;
    for (size_t i = 0; i < minimizers.size(); i++) {
        size_t minimizer_num = minimizers_in_order[i];

#ifdef TRACK_PROVENANCE
        // Say we're working on it
        funnel.processing_input(minimizer_num);
#endif

        // Select the minimizer if it is informative enough or if the total score
        // of the selected minimizers is not high enough.
        size_t hits = minimizer_index->count(minimizers[minimizer_num]);
        if (hits <= hit_cap || (hits <= hard_hit_cap && selected_score + minimizer_score[minimizer_num] <= target_score)) {

            // Locate the hits.
            for (auto& hit : minimizer_index->find(minimizers[minimizer_num])) {
                // Reverse the hits for a reverse minimizer
                if (minimizers[minimizer_num].is_reverse) {
                    size_t node_length = gbwt_graph.get_length(gbwt_graph.get_handle(id(hit)));
                    hit = reverse_base_pos(hit, node_length);
                }
                // For each position, remember it and what minimizer it came from
                seeds.push_back(hit);
                seed_to_source.push_back(minimizer_num);
            }
            selected_score += minimizer_score[minimizer_num];
            
#ifdef TRACK_PROVENANCE
            // Record in the funnel that this minimizer gave rise to these seeds.
            funnel.expand(minimizer_num, hits);
#endif
        } else {
            rejected_count++;
            
#ifdef TRACK_PROVENANCE
            // Record in the funnel thast we rejected it
            funnel.kill(minimizer_num);
#endif
        }
        
#ifdef TRACK_PROVENANCE
        // Say we're done with this input item
        funnel.processed_input();
#endif
    }


#ifdef TRACK_PROVENANCE
#ifdef TRACK_CORRECTNESS
    // Tag seeds with correctness based on proximity along paths to the input read's refpos
    funnel.substage("correct");
    
    if (aln.refpos_size() != 0) {
        // Take the first refpos as the true position.
        auto& true_pos = aln.refpos(0);
        
        for (size_t i = 0; i < seeds.size(); i++) {
            // Find every seed's reference positions. This maps from path name to pairs of offset and orientation.
            auto offsets = xg_index->nearest_offsets_in_paths(seeds[i], 100);
            for (auto& hit_pos : offsets[true_pos.name()]) {
                // Look at all the ones on the path the read's true position is on.
                if (abs((int64_t)hit_pos.first - (int64_t) true_pos.offset()) < 200) {
                    // Call this seed hit close enough to be correct
                    funnel.tag_correct(i);
                }
            }
        }
    }
#endif
#endif
        
#ifdef debug
    cerr << "Read " << aln.name() << ": " << aln.sequence() << endl;
    cerr << "Found " << seeds.size() << " seeds from " << (minimizers.size() - rejected_count) << " minimizers, rejected " << rejected_count << endl;
#endif

#ifdef TRACK_PROVENANCE
    // Begin the clustering stage
    funnel.stage("cluster");
#endif
        
    // Cluster the seeds. Get sets of input seed indexes that go together.
    vector<vector<size_t>> clusters = clusterer.cluster_seeds(seeds, distance_limit);
    
#ifdef TRACK_PROVENANCE
    funnel.substage("score");
#endif

    // Cluster score is the sum of minimizer scores.
    std::vector<double> cluster_score(clusters.size(), 0.0);
    for (size_t i = 0; i < clusters.size(); i++) {
        // For each cluster
        auto& cluster = clusters[i];
        
#ifdef TRACK_PROVENANCE
        // Say we're making it
        funnel.producing_output(i);
#endif

        // Which minimizers are present in the cluster.
        vector<bool> present(minimizers.size(), false);
        for (auto hit_index : cluster) {
            present[seed_to_source[hit_index]] = true;
        }

        // Compute the score.
        for (size_t j = 0; j < minimizers.size(); j++) {
            if (present[j]) {
                cluster_score[i] += minimizer_score[j];
            }
        }
        
#ifdef TRACK_PROVENANCE
        // Record the cluster in the funnel as a group of the size of the number of items.
        funnel.merge_group(cluster.begin(), cluster.end());
        funnel.score(funnel.latest(), cluster_score[i]);
        
        // Say we made it.
        funnel.produced_output();
#endif
    }

#ifdef debug
    cerr << "Found " << clusters.size() << " clusters" << endl;
#endif
    
    // Make a vector of cluster indexes to sort
    vector<size_t> cluster_indexes_in_order;
    cluster_indexes_in_order.reserve(clusters.size());
    for (size_t i = 0; i < clusters.size(); i++) {
        cluster_indexes_in_order.push_back(i);
    }

    // Put the most covering cluster's index first
    std::sort(cluster_indexes_in_order.begin(), cluster_indexes_in_order.end(), [&](const size_t& a, const size_t& b) -> bool {
        // Return true if a must come before b, and false otherwise
        return (cluster_score[a] > cluster_score[b]);
    });
    
#ifdef TRACK_PROVENANCE
    // Now we go from clusters to gapless extensions
    funnel.stage("extend");
#endif
    
    // These are the GaplessExtensions for all the clusters, in cluster_indexes_in_order order.
    vector<vector<GaplessExtension>> cluster_extensions;
    cluster_extensions.reserve(cluster_indexes_in_order.size());
    
    for (size_t i = 0; i < clusters.size() && i < max_extensions; i++) {
        // For each cluster, in sorted order
        size_t& cluster_num = cluster_indexes_in_order[i];
        
#ifdef TRACK_PROVENANCE
        funnel.processing_input(cluster_num);
#endif

        vector<size_t>& cluster = clusters[cluster_num];

#ifdef debug
        cerr << "Cluster " << cluster_num << " rank " << i << ": " << endl;
#endif
         
        // Pack the seeds into (read position, graph position) pairs.
        vector<pair<size_t, pos_t>> seed_matchings;
        seed_matchings.reserve(cluster.size());
        for (auto& seed_index : cluster) {
            // For each seed in the cluster, generate its matching pair
            seed_matchings.emplace_back(minimizers[seed_to_source[seed_index]].offset, seeds[seed_index]);
#ifdef debug
            cerr << "Seed read:" << minimizers[seed_to_source[seed_index]].offset << " = " << seeds[seed_index]
                << " from minimizer " << seed_to_source[seed_index] << "(" << minimizer_index->count(minimizers[seed_to_source[seed_index]]) << ")" << endl;
#endif
        }
        
        // Extend seed hits in the cluster into one or more gapless extensions
        cluster_extensions.emplace_back(extender.extend(seed_matchings, aln.sequence()));
        
#ifdef TRACK_PROVENANCE
        // Record with the funnel that the previous group became a group of this size.
        // Don't bother recording the seed to extension matching...
        funnel.project_group(cluster_num, cluster_extensions.back().size());
        
        // Say we finished with this cluster, for now.
        funnel.processed_input();
#endif
    }
    
#ifdef TRACK_PROVENANCE
    funnel.substage("score");
#endif
    
    // Now score all the gapless extensions by max match count accounted for.
    vector<int> cluster_extension_scores;
    cluster_extension_scores.reserve(cluster_extensions.size());
    for (size_t i = 0; i < cluster_extensions.size(); i++) {
        // For each group of GaplessExtensions
        
#ifdef TRACK_PROVENANCE
        funnel.producing_output(i);
#endif
        
        auto& extensions = cluster_extensions[i];
        // Count the matches suggested by the group and use that as a score.
        cluster_extension_scores.push_back(estimate_extension_group_score(aln, extensions));
        
#ifdef TRACK_PROVENANCE
        // Record the score with the funnel
        funnel.score(i, cluster_extension_scores.back());
        funnel.produced_output();
#endif
    }
    
    // Now sort them by score
    vector<size_t> extension_indexes_in_order;
    extension_indexes_in_order.reserve(cluster_extension_scores.size());
    for (size_t i = 0; i < cluster_extension_scores.size(); i++) {
        extension_indexes_in_order.push_back(i);
    }
    
    // Put the most matching group of extensions first
    std::sort(extension_indexes_in_order.begin(), extension_indexes_in_order.end(), [&](const int& a, const int& b) -> bool {
        // Return true if a must come before b, and false otherwise
        return cluster_extension_scores.at(a) > cluster_extension_scores.at(b);
    });
    
#ifdef TRACK_PROVENANCE
    funnel.stage("align");
#endif
    
    // Now start the alignment step. Everything has to become an alignment.
    
    // We will fill this with all computed alignments in estimated score order.
    vector<Alignment> alignments;
    alignments.reserve(extension_indexes_in_order.size());
    
    // Clear any old refpos annotation and path
    aln.clear_refpos();
    aln.clear_path();
    aln.set_score(0);
    aln.set_identity(0);
    aln.set_mapping_quality(0);
    
    // Go through the gapless extension groups in score order.
    // Keep track of best and second best scores.
    int best_score = 0;
    int second_best_score = 0;
    for (size_t i = 0; i < extension_indexes_in_order.size() && i < max_alignments; i++) {
        // Find the extension group we are talking about
        size_t& extension_num = extension_indexes_in_order[i];
        
#ifdef TRACK_PROVENANCE
        funnel.processing_input(extension_num);
#endif

        auto& extensions = cluster_extensions[extension_num];
        
        if (i < 2 || score_is_significant(cluster_extension_scores[extension_num], best_score, second_best_score)) {
            // Always take the first and second.
            // For later ones, check if this score is significant relative to the running best and second best scores.
            
            // If so, get an Alignment out of it somehow, and throw it in.
            alignments.emplace_back(aln);
            Alignment& out = alignments.back();
            
            if (extensions.size() == 1 && extensions[0].full()) {
                // We got a full-length extension, so directly convert to an Alignment.
                
#ifdef TRACK_PROVENANCE
                funnel.substage("direct");
#endif

                *out.mutable_path() = extensions.front().to_path(gbwt_graph, out.sequence());
                
                // The score estimate is exact.
                int alignment_score = cluster_extension_scores[extension_num];
                
                // Compute identity from mismatch count.
                size_t mismatch_count = extensions[0].mismatches();
                double identity = out.sequence().size() == 0 ? 0.0 : (out.sequence().size() - mismatch_count) / (double) out.sequence().size();
                
                // Fill in the score and identity
                out.set_score(alignment_score);
                out.set_identity(identity);
                
#ifdef TRACK_PROVENANCE
                // Stop the current substage
                funnel.substage_stop();
#endif
            } else if (do_chaining) {
                // We need to do chaining.
                
#ifdef TRACK_PROVENANCE
                funnel.substage("chain");
#endif
                
                // Do the chaining and compute an alignment into out.
                bool chain_success = chain_extended_seeds(aln, extensions, out);
                
#ifdef TRACK_PROVENANCE
                // We're done chaining. Next alignment may not go through this substage.
                funnel.substage_stop();
#endif

                if (!chain_success) {
                    // We thought chaining would be too hard. Fall back on context extraction
                    
#ifdef TRACK_PROVENANCE
                    funnel.substage("context");
#endif

                    align_to_local_haplotypes(aln, extensions, out);
                 
#ifdef TRACK_PROVENANCE
                    funnel.substage_stop();
#endif
                
                }
            } else {
                // We would do chaining but it is disabled.
                // Leave out unaligned
            }
            
            // Update the running best and second best scores.
            if (out.score() > best_score) {
                second_best_score = best_score;
                best_score = out.score();
            } else if (out.score() > second_best_score) {
                second_best_score = out.score();
            }
            
#ifdef TRACK_PROVENANCE
            // Record the Alignment and its score with the funnel
            funnel.project(extension_num);
            funnel.score(i, out.score());
            
            // We're done with this input item
            funnel.processed_input();
#endif
        } else {
            // If this score is insignificant, nothing past here is significant.
            // Don't do any more.
            
#ifdef TRACK_PROVENANCE
            funnel.kill_all(extension_indexes_in_order.begin() + i, extension_indexes_in_order.end());
            funnel.processed_input();
#endif
            
            break;
        }
    }
    
    if (alignments.size() == 0) {
        // Produce an unaligned Alignment
        alignments.emplace_back(aln);
        
#ifdef TRACK_PROVENANCE
        // Say it came from nowhere
        funnel.introduce();
#endif
    }
    
    // Order the Alignments by score
    vector<size_t> alignments_in_order;
    alignments_in_order.reserve(alignments.size());
    for (size_t i = 0; i < alignments.size(); i++) {
        alignments_in_order.push_back(i);
    }
    
    // Sort again by actual score instead of extennsion score
    std::sort(alignments_in_order.begin(), alignments_in_order.end(), [&](const size_t& a, const size_t& b) -> bool {
        // Return true if a must come before b (i.e. it has a larger score)
        return alignments[a].score() > alignments[b].score();
    });
    
#ifdef TRACK_PROVENANCE
    // Now say we are finding the winner(s)
    funnel.stage("winner");
#endif
    
    vector<Alignment> mappings;
    mappings.reserve(min(alignments_in_order.size(), max_multimaps));
    for (size_t i = 0; i < alignments_in_order.size() && i < max_multimaps; i++) {
        // For each output slot, fill it with the alignment at that rank if available.
        size_t& alignment_num = alignments_in_order[i];
        mappings.emplace_back(std::move(alignments[alignment_num]));
        
#ifdef TRACK_PROVENANCE
        // Tell the funnel
        funnel.project(alignment_num);
        funnel.score(alignment_num, mappings.back().score());
#endif
    }

#ifdef TRACK_PROVENANCE
    if (max_multimaps < alignments_in_order.size()) {
        // Some things stop here
        funnel.kill_all(alignments_in_order.begin() + max_multimaps, alignments_in_order.end());
    }
    
    funnel.substage("mapq");
#endif
    
    // Grab all the scores for MAPQ computation.
    vector<double> scores;
    scores.reserve(alignments.size());
    for (size_t i = 0; i < mappings.size(); i++) {
        // Grab the scores of the alignments we are outputting
        scores.push_back(mappings[i].score());
    }
    for (size_t i = mappings.size(); i < alignments_in_order.size(); i++) {
        // And of the alignments we aren't
        scores.push_back(alignments[alignments_in_order[i]].score());
    }
        
#ifdef debug
    cerr << "For scores ";
    for (auto& score : scores) cerr << score << " ";
#endif

    size_t winning_index;
    double mapq = get_regular_aligner()->maximum_mapping_quality_exact(scores, &winning_index);
    
#ifdef debug
    cerr << "MAPQ is " << mapq << endl;
#endif
        
    // Make sure to clamp 0-60.
    mappings.front().set_mapping_quality(max(min(mapq, 60.0), 0.0));
    
#ifdef TRACK_PROVENANCE
    funnel.substage_stop();
#endif
    
    for (size_t i = 0; i < mappings.size(); i++) {
        // For each output alignment in score order
        auto& out = mappings[i];
        
        // Assign primary and secondary status
        out.set_is_secondary(i > 0);
    }
    
    // Stop this alignment
    funnel.stop();
    
#ifdef TRACK_PROVENANCE
    
    // And with the number of results in play at each stage
    funnel.for_each_stage([&](const string& stage, size_t result_count) {
        set_annotation(mappings[0], "stage_" + stage + "_results", (double)result_count);
    });

#ifdef TRACK_CORRECTNESS
    // And with the last stage at which we had any descendants of the correct seed hit locations
    set_annotation(mappings[0], "last_correct_stage", funnel.last_correct_stage());
#endif
#endif
    
    // Ship out all the aligned alignments
    alignment_emitter.emit_mapped_single(std::move(mappings));

#ifdef debug
    // Dump the funnel info graph.
    funnel.to_dot(cerr);
#endif
}

int MinimizerMapper::estimate_extension_group_score(const Alignment& aln, vector<GaplessExtension>& extended_seeds) const {
    if (extended_seeds.empty()) {
        // TODO: We should never see an empty group of extensions
        return 0;
    } else if (extended_seeds.size() == 1 && extended_seeds.front().full()) {
        // This is a full length match. Compute exact score.
        // TODO: Should we use the aligner instead of computing the score here?

        const Aligner* aligner = get_regular_aligner();
        return (aln.sequence().length() - extended_seeds.front().mismatches()) * aligner->match -
               extended_seeds.front().mismatches() * aligner->mismatch +
               2 * aligner->full_length_bonus;
    } else {
        // This is a collection of one or more non-full-length extended seeds.
        
        if (aln.sequence().size() == 0) {
            // No score here
            return 0;
        }
        
        // Now we compute an estimate of the score: match count for all the
        // flank bases that aren't universal mismatches, mismatch count for
        // those that are.
        int score_estimate = 0;
        
        // We use a sweep line algorithm.
        // This records the last base to be covered by the current sweep line.
        int64_t sweep_line = 0;
        // This records the first base not covered by the last sweep line.
        int64_t last_sweep_line = 0;
        
        // And we track the next unentered gapless extension
        size_t unentered = 0;
        
        // Extensions we are in are in this min-heap of past-end position and gapless extension number.
        vector<pair<size_t, size_t>> end_heap;
        // The heap uses this comparator
        auto compare = [](const pair<size_t, size_t>& a, const pair<size_t, size_t>& b) {
            // Return true if a must come later in the heap than b
            return a.first > b.first;
        };
        
        while(last_sweep_line < aln.sequence().size()) {
            // We are processed through the position before last_sweep_line.
            
            // Find a place for sweep_line to go
            
            // Find the next seed start
            int64_t next_seed_start = numeric_limits<int64_t>::max();
            if (unentered < extended_seeds.size()) {
                next_seed_start = extended_seeds[unentered].flanked_interval.first;
            }
            
            // Find the next mismatch
            int64_t next_mismatch = numeric_limits<int64_t>::max();
            for (auto& overlapping : end_heap) {
                // For each gapless extension we overlap, find its sorted mismatches
                auto& mismatches = extended_seeds[overlapping.second].mismatch_positions;
                for (auto& mismatch : mismatches) {
                    if (mismatch < last_sweep_line) {
                        // Already accounted for
                        continue;
                    }
                    if (mismatch < next_mismatch) {
                        // We found a new one
                        next_mismatch = mismatch;
                    }
                    // We only care about the first one not too early.
                    break;
                }
            }
            
            // Find the next seed end
            int64_t next_seed_end = numeric_limits<int64_t>::max();
            if (!end_heap.empty()) {
                next_seed_end = end_heap.front().first;
            }
            
            // Whichever is closer between those points and the end, do that.
            sweep_line = min(min(min(next_seed_end, next_mismatch), next_seed_start), (int64_t) aln.sequence().size() - 1);
            
            // So now we're only interested in things that happen at sweep_line.
            
            if (!end_heap.empty()) {
                // If we were covering anything, count matches between last_sweep_line and here, not including at last_sweep_line.
                score_estimate += get_regular_aligner()->score_exact_match(aln, last_sweep_line, sweep_line - last_sweep_line);
            }
            
            while(!end_heap.empty() && end_heap.front().first == sweep_line) {
                // Take out anything that past-ends here
                std::pop_heap(end_heap.begin(), end_heap.end());
                end_heap.pop_back();
            }
            
            while (unentered < extended_seeds.size() && extended_seeds[unentered].flanked_interval.first == sweep_line) {
                // Bring in anything that starts here
                end_heap.emplace_back(extended_seeds[unentered].flanked_interval.second, unentered);
                std::push_heap(end_heap.begin(), end_heap.end());
                unentered++;
            }
            
            if (!end_heap.empty()) {
                // We overlap some seeds
            
                // Count up mismatches that are here and extended seeds that overlap here
                size_t mismatching_count = 0;
                for (auto& overlapping : end_heap) {
                    // For each gapless extension we overlap, find its sorted mismatches
                    auto& mismatches = extended_seeds[overlapping.second].mismatch_positions;
                    for (auto& mismatch : mismatches) {
                        if (mismatch < last_sweep_line) {
                            // Already accounted for
                            continue;
                        }
                        if (mismatch == sweep_line) {
                            // We found a new one here
                            mismatching_count++;
                        }
                        // We only care about the first one not too early.
                        break;
                    }
                }
                
                if (mismatching_count == end_heap.size()) {
                    // This is a universal mismatch
                    // Add a mismatch to the score
                    score_estimate += get_regular_aligner()->score_mismatch(1);
                } else {
                    // Add a 1-base match to the score
                    score_estimate += get_regular_aligner()->score_exact_match(aln, sweep_line, 1);
                }
            }
            
            // If we don't overlap any seeds here, we won't score any matches or mismatches.
            
            // Move last_sweep_line to sweep_line.
            // We need to add 1 since last_sweep_line is the next *un*included base
            last_sweep_line = sweep_line + 1;
        }
        
        // TODO: should we apply full length bonuses?
        
        // When we get here, the score estimate is finished.
        return score_estimate;
    }
    
}

bool MinimizerMapper::score_is_significant(int score_estimate, int best_score, int second_best_score) const {
    // mpmap uses a heuristic of if the read coverage of the cluster is less than half the best cluster's read coverage, stop.
    // We do something similar, but with scores. And we make sure to get at least one second best score.
    // If it's not more than half the best score, it doesn't matter if it beats the second best score; both secondaries are sufficiently bad.
    // TODO: real scores from full-length gapless extensions aren't quite directly comparable with estimates.
    if (score_estimate * 2 >= best_score || second_best_score < 1) {
        return true;
    }
    return false;
}

bool MinimizerMapper::chain_extended_seeds(const Alignment& aln, const vector<GaplessExtension>& extended_seeds, Alignment& out) const {

#ifdef debug
    cerr << "Trying again to chain " << extended_seeds.size() << " extended seeds" << endl;
#endif

    // Find the paths between pairs of extended seeds that agree with haplotypes.
    // We don't actually need the read sequence for this, just the read length for longest gap computation.
    // The paths in the seeds know the hit length.
    // We assume all overlapping hits are exclusive.
    unordered_map<size_t, unordered_map<size_t, vector<Path>>> paths_between_seeds = find_connecting_paths(extended_seeds,
        aln.sequence().size());
        
        
    // Now we need to identify the sources and sinks in the reachability graph (again)
    // TODO: keep from find_connecting_paths
    unordered_set<size_t> source_extensions;
    unordered_set<size_t> sink_extensions;
    
    for (size_t i = 0; i < extended_seeds.size(); i++) {
        // start out assuming everything is a source and a sink
        source_extensions.insert(i);
        sink_extensions.insert(i);
    }
    
    for (auto& from_and_dests : paths_between_seeds) {
        // For each reachability edge from extension
        if (from_and_dests.first == numeric_limits<size_t>::max()) {
            // Skip edges from nowhere
            continue;
        }
        
        if ((from_and_dests.second.size() == 1 && !from_and_dests.second.count(numeric_limits<size_t>::max())) ||
            (from_and_dests.second.size() > 1)) {
            // Mark as not a sink if it goes anywhere other than out of the cluster 
            sink_extensions.erase(from_and_dests.first);
        }
        
        for (auto& to_and_paths : from_and_dests.second) {
            // For everywhere we can get from here
            if (to_and_paths.first == numeric_limits<size_t>::max()) {
                // Discount going nowhere
                continue;
            }
            
            if (!to_and_paths.second.empty()) {
                // If there are any actual paths, mark the edge to extension as
                // reachable from somewhere else.
                source_extensions.erase(to_and_paths.first);
            }
        }
    }
    
    assert(!source_extensions.empty());
    assert(!sink_extensions.empty());
   
    if (source_extensions.size() + sink_extensions.size() > max_tails) {
        return false;
    }
   
    // We're going to record source and sink path count distributions, for debugging
    vector<double> tail_path_counts;
    
    // We're also going to record read sequence lengths for tails
    vector<double> tail_lengths;
    
    // Make a MultipathAlignment and feed in all the extended seeds as subpaths
    MultipathAlignment mp;
    // Pull over all the non-alignment data (to get copied back out when linearizing)
    transfer_read_metadata(aln, mp);
    for (auto& extended_seed : extended_seeds) {
        Subpath* s = mp.add_subpath();
        // Copy in the path.
        *s->mutable_path() = extended_seed.to_path(gbwt_graph, out.sequence());
        // Score it
        s->set_score(get_regular_aligner()->score_partial_alignment(aln, gbwt_graph, s->path(),
            aln.sequence().begin() + extended_seed.core_interval.first));
        // The position in the read it occurs at will be handled by the multipath topology.
        if (extended_seed.core_interval.first == 0) {
            // But if it occurs at the very start of the read we need to mark that now.
            mp.add_start(mp.subpath_size() - 1);
        }
    }
    
    // Handle the leading/left tails
    if (linear_tails) {
        // Handle left tails as several parallel strings
        
#ifdef debug
        cerr << "Handle " << paths_between_seeds[numeric_limits<size_t>::max()].size() << " left tails linearly" << endl;
#endif

        for (auto& kv : paths_between_seeds[numeric_limits<size_t>::max()]) {
            // For each extended seed that can come from something outside the cluster
            const size_t& source = kv.first;
            
            // Grab the part of the read sequence that comes before it
            string before_sequence = aln.sequence().substr(0, extended_seeds[source].core_interval.first);
            
            // Record that a source has this many incoming haplotypes to process.
            tail_path_counts.push_back(kv.second.size());
            // Against a sequence this long
            tail_lengths.push_back(before_sequence.size());
            
            // Do right-pinned alignment
            pair<Path, int64_t> result = get_best_alignment_against_any_path(kv.second, before_sequence,
                extended_seeds[source].starting_position(gbwt_graph), true, false);

            // Grab the best path in backing graph space (which may be empty)
            Path& best_path = result.first;
            // And its score
            int64_t& best_score = result.second;
            
            // Put it in the MultipathAlignment
            Subpath* s = mp.add_subpath();
            *s->mutable_path() = std::move(best_path);
            s->set_score(best_score);
            
            // And make the edge from it to the correct source
            s->add_next(source);
            
#ifdef debug
            cerr << "Add leading tail " << (mp.subpath_size() - 1) << " -> " << source << endl;
#endif
            
            // And mark it as a start subpath
            mp.add_start(mp.subpath_size() - 1);
        }
    } else {
        // Handle left tails as a forest of trees
       
        // Get the forests of all left tails by extension they belong to
        auto tails_by_extension = get_tail_forests(extended_seeds, aln.sequence().size(), paths_between_seeds, true);
        
#ifdef debug
        cerr << "Handle " << tails_by_extension.size() << " left tail forests" << endl;
#endif
        
        for (auto& kv : tails_by_extension) {
            // For each source extension
            const size_t& source = kv.first;
            
            // Grab the part of the read sequence that comes before it
            string before_sequence = aln.sequence().substr(0, extended_seeds[source].core_interval.first);
            
            // Do right-pinned alignment
            pair<Path, int64_t> result = get_best_alignment_against_any_tree(kv.second, before_sequence,
                extended_seeds[source].starting_position(gbwt_graph), false);

            // Grab the best path in backing graph space (which may be empty)
            Path& best_path = result.first;
            // And its score
            int64_t& best_score = result.second;
            
            // Put it in the MultipathAlignment
            Subpath* s = mp.add_subpath();
            *s->mutable_path() = std::move(best_path);
            s->set_score(best_score);
            
            // And make the edge from it to the correct source
            s->add_next(source);
            
#ifdef debug
            cerr << "Add leading tail " << (mp.subpath_size() - 1) << " -> " << source << endl;
#endif
            
            // And mark it as a start subpath
            mp.add_start(mp.subpath_size() - 1);
        }
    }
    
    // We must have somewhere to start.
    assert(mp.start_size() > 0);

    for (auto& from_and_edges : paths_between_seeds) {
        const size_t& from = from_and_edges.first;
        if (from == numeric_limits<size_t>::max()) {
            continue;
        }
        // Then for all the other from extended seeds

        // Work out where the extended seed ends in the read
        size_t from_end = extended_seeds[from].core_interval.second;
        
        for (auto& to_and_paths : from_and_edges.second) {
            const size_t& to = to_and_paths.first;
            // For all the edges to other extended seeds
            
#ifdef debug
            cerr << "Consider " << to_and_paths.second.size() << " paths between extended seeds " << to << " and " << from << endl;
#endif

            if (to == numeric_limits<size_t>::max()) {
                // We can go to something outside the cluster
                
                // Do a bunch of left pinned alignments for the tails.
                
                // Find the sequence
                string trailing_sequence = aln.sequence().substr(from_end);
                
                if (!trailing_sequence.empty()) {
                    // There is actual trailing sequence to align on this escape path
                    
                    // Record that a sink has this many outgoing haplotypes to process.
                    tail_path_counts.push_back(to_and_paths.second.size());
                    // Against a sequence this size
                    tail_lengths.push_back(trailing_sequence.size());

                    // Do left-pinned alignment
                    pair<Path, int64_t> result = get_best_alignment_against_any_path(to_and_paths.second, trailing_sequence,
                        extended_seeds[from].tail_position(gbwt_graph), true, true);

                    // Grab the best path in backing graph space (which may be empty)
                    Path& best_path = result.first;
                    // And its score
                    int64_t& best_score = result.second;

                    // Put it in the MultipathAlignment
                    Subpath* s = mp.add_subpath();
                    *s->mutable_path() = std::move(best_path);
                    s->set_score(best_score);
                    
                    // And make the edge to hook it up
                    mp.mutable_subpath(from)->add_next(mp.subpath_size() - 1);
                    
#ifdef debug
                    cerr << "Add trailing tail " << from << " -> " << (mp.subpath_size() - 1) << endl;
#endif
                }
                
                // If there's no sequence to align on the path going off to nowhere, don't do anything.
            } else {
                // Do alignments between from and to

                // Find the sequence
                assert(extended_seeds[to].core_interval.first >= from_end);
                string intervening_sequence = aln.sequence().substr(from_end, extended_seeds[to].core_interval.first - from_end);
                
#ifdef debug
                cerr << "Connect " << pb2json(extended_seeds[from].tail_position(gbwt_graph))
                    << " and " << pb2json(extended_seeds[to].starting_position(gbwt_graph)) << endl;
#endif

                // Do un-pinned alignment
                pair<Path, int64_t> result = get_best_alignment_against_any_path(to_and_paths.second, intervening_sequence,
                    extended_seeds[from].tail_position(gbwt_graph), false, false);

                // Grab the best path in backing graph space (which may be empty)
                Path& best_path = result.first;
                // And its score
                int64_t& best_score = result.second;

                // We may have an empty path. That's fine.

                if (best_path.mapping_size() == 0 && intervening_sequence.empty()) {
                    // We just need an edge from from to to
#ifdef debug
                    cerr << "Add direct edge " << from << " -> " << to << endl;
#endif
                    mp.mutable_subpath(from)->add_next(to);
                } else {
                    // We need to connect from and to with a Subpath with this path

                    // We really should have gotten something
                    assert(best_path.mapping_size() != 0);

                    // Put it in the MultipathAlignment
                    Subpath* s = mp.add_subpath();
                    *s->mutable_path() = std::move(best_path);
                    s->set_score(best_score);
                    
                    // And make the edges to hook it up
                    s->add_next(to);
                    mp.mutable_subpath(from)->add_next(mp.subpath_size() - 1);
                    
#ifdef debug
                    cerr << "Add material along edge " << from << " -> " << (mp.subpath_size() - 1) << " -> " << to << endl;
#endif
                    
                }

            }
            
        }

    }
    
    if (!linear_tails) {
        // Handle right tails as a forest of trees
   
        // Get the forests of all right tails by extension they belong to
        auto tails_by_extension = get_tail_forests(extended_seeds, aln.sequence().size(), paths_between_seeds, false);
        
        for (auto& kv : tails_by_extension) {
            // For each source extension
            const size_t& from = kv.first;
            
#ifdef debug
            cerr << "Consider right tails for extension " << from << " with interval "
                << extended_seeds[from].core_interval.first << " - " << extended_seeds[from].core_interval.second << endl;
#endif
            
            // Find the sequence
            string trailing_sequence = aln.sequence().substr(extended_seeds[from].core_interval.second);
            
            // There should be actual trailing sequence to align on this escape path
            assert(!trailing_sequence.empty());
            
            // Do left-pinned alignment
            pair<Path, int64_t> result = get_best_alignment_against_any_tree(kv.second, trailing_sequence,
                extended_seeds[from].tail_position(gbwt_graph), true);

            // Grab the best path in backing graph space (which may be empty)
            Path& best_path = result.first;
            // And its score
            int64_t& best_score = result.second;
            
            // Put it in the MultipathAlignment
            Subpath* s = mp.add_subpath();
            *s->mutable_path() = std::move(best_path);
            s->set_score(best_score);
            
            // And make the edge to hook it up
            mp.mutable_subpath(from)->add_next(mp.subpath_size() - 1);
            
#ifdef debug
            cerr << "Add trailing tail " << from << " -> " << (mp.subpath_size() - 1) << endl;
#endif
                
        }
    }

        
    // Then we take the best linearization of the full MultipathAlignment.
    // Make sure to force source to sink
    topologically_order_subpaths(mp);

    if (!validate_multipath_alignment(mp, gbwt_graph)) {
        // If we generated an invalid multipath alignment, we did something wrong and need to stop
        cerr << "error[vg::MinimizerMapper]: invalid MultipathAlignment generated: " << pb2json(mp) << endl;
        exit(1);
    }
    
    // Linearize into the out alignment, copying path, score, and also sequence and other read metadata
    optimal_alignment(mp, out, true);
    // Compute the identity from the path.
    out.set_identity(identity(out.path()));
    
    // Save all the tail alignment debugging statistics
    set_annotation(out, "tail_path_counts", tail_path_counts);
    set_annotation(out, "tail_lengths", tail_lengths);
    
    return true;
}

pair<Path, size_t> MinimizerMapper::get_best_alignment_against_any_path(const vector<Path>& paths,
    const string& sequence, const Position& default_position, bool pinned, bool pin_left) const {
    
    // We want the best alignment, to the base graph, done against any target path
    Path best_path;
    // And its score
    int64_t best_score = numeric_limits<int64_t>::min();
    
    // We must have some target paths
    assert(!paths.empty());
    
    // We can align it once per target path
    for (auto& path : paths) {
        // For each path we can take to get to the source
        
#ifdef debug
        cerr << "Consider " << sequence.size() << " bp against path of " << path_from_length(path) << " bp" << endl;
#endif
        
        if (path.mapping_size() == 0) {
            // There's no graph bases here
            if (pinned) {
        
                // We might have extra read outside the graph. Handle leading insertions.
                // We consider a pure softclip.
                // We don't consider an empty sequence because if that were the case
                // we would not have any paths_between_seeds entries for the dangling-left-sequence sentinel.
                if (best_score < 0) {
                
                    best_score = 0;
                    best_path.clear_mapping();
                    Mapping* m = best_path.add_mapping();
                    Edit* e = m->add_edit();
                    e->set_from_length(0);
                    e->set_to_length(sequence.size());
                    e->set_sequence(sequence);
                    // Since the softclip consumes no graph, we place it on the node we are going to.
                    *m->mutable_position() = default_position;
                    
#ifdef debug
                    cerr << "New best alignment: " << pb2json(best_path) << " score " << best_score << endl;
#endif
                }
            } else {
                // We're aligning against nothing globally
                if (sequence.empty()) {
                    // Consider the nothing to nothing alignment, score 0
                    if (best_score < 0) {
                        best_score = 0;
                        best_path.clear_mapping();
#ifdef debug
                        cerr << "New best alignment: " << pb2json(best_path) << " score " << best_score << endl;
#endif
                    }
                } else {
                    // Consider the something to nothing alignment.
                    // We can't use the normal code path because the BandedGlobalAligner 
                    // wouldn't be able to generate a position form an empty graph.
                    
                    // We know the extended seeds we are between won't start/end with gaps, so we own the gap open.
                    int64_t score = get_regular_aligner()->score_gap(sequence.size());
                    if (score > best_score) {
                        best_score = score;
                        best_path.clear_mapping();
                        Mapping* m = best_path.add_mapping();
                        Edit* e = m->add_edit();
                        e->set_from_length(0);
                        e->set_to_length(sequence.size());
                        e->set_sequence(sequence);
                        // We can copy the position of where we are going to, since we consume no graph.
                        *m->mutable_position() = default_position;
                    
#ifdef debug
                        cerr << "New best alignment: " << pb2json(best_path) << " score " << best_score << endl;
#endif
                    
                    }
                }
            }
        } else {
            // This path has bases in it

            // Make a subgraph.
            // TODO: don't copy the path
            PathSubgraph subgraph(&gbwt_graph, path);
            
            // Do alignment to the path subgraph with GSSWAligner.
            Alignment current_alignment;
            current_alignment.set_sequence(sequence);
#ifdef debug
            cerr << "Align " << pb2json(current_alignment) << (pinned ? (pin_left ? " pinned left" : " pinned right") : " global");

#ifdef debug_dump_graph
            cerr << " vs:" << endl;
            subgraph.for_each_handle([&](const handle_t& here) {
                cerr << subgraph.get_id(here) << " (" << subgraph.get_sequence(here) << "): " << endl;
                subgraph.follow_edges(here, true, [&](const handle_t& there) {
                    cerr << "\t" << subgraph.get_id(there) << " (" << subgraph.get_sequence(there) << ") ->" << endl;
                });
                subgraph.follow_edges(here, false, [&](const handle_t& there) {
                    cerr << "\t-> " << subgraph.get_id(there) << " (" << subgraph.get_sequence(there) << ")" << endl;
                });
            });
            cerr << "Path: " << pb2json(path) << endl;
#else
            cerr << endl;
#endif
#endif
            
            // Align, accounting for full length bonus
            
            if (pinned) {
            
                if (use_xdrop_for_tails) {
#ifdef debug
                    Alignment clone = current_alignment;
                    get_regular_aligner()->align_pinned(clone, subgraph, subgraph.get_topological_order(), pin_left);
#endif
                    get_regular_aligner()->get_xdrop()->align_pinned(current_alignment, subgraph, subgraph.get_topological_order(), pin_left);
#ifdef debug
                    cerr << "Xdrop: " << pb2json(current_alignment) << endl;
                    cerr << "Normal: " << pb2json(clone) << endl;
#endif
                } else {
                    get_regular_aligner()->align_pinned(current_alignment, subgraph, subgraph.get_topological_order(), pin_left);
                }
            } else {
                get_regular_aligner()->align_global_banded(current_alignment, subgraph, 5, true);
            }
            
#ifdef debug
            cerr << "\tScore: " << current_alignment.score() << endl;
#endif
            
            if (current_alignment.score() > best_score) {
                // This is a new best alignment. Translate from subgraph into base graph and keep it
                best_path = subgraph.translate_down(current_alignment.path());
                best_score = current_alignment.score();
                
#ifdef debug
                cerr << "New best alignment against: " << pb2json(path) << " is "
                    << pb2json(best_path) << " score " << best_score << endl;
#endif
            }
        }
    }

    // We really should have gotten something to replace the placeholder score
    assert(best_score != numeric_limits<int64_t>::min());
    
    return make_pair(best_path, best_score);
}

pair<Path, size_t> MinimizerMapper::get_best_alignment_against_any_tree(const vector<TreeSubgraph>& trees,
    const string& sequence, const Position& default_position, bool pin_left) const {
    
    // We want the best alignment, to the base graph, done against any target path
    Path best_path;
    // And its score
    int64_t best_score = numeric_limits<int64_t>::min();
    
    // We can align it once per target tree
    for (auto& subgraph : trees) {
        // For each tree we can map against, map pinning the correct edge of the sequence to the root.
        
        if (subgraph.get_node_count() == 0) {
            // There's no graph bases here
            // We might have extra read outside the graph. Handle leading insertions.
            // We consider a pure softclip, since all alignment here is pinning.
            // We don't consider an empty sequence since that would produce no trees.
            if (best_score < 0) {
            
                best_score = 0;
                best_path.clear_mapping();
                Mapping* m = best_path.add_mapping();
                Edit* e = m->add_edit();
                e->set_from_length(0);
                e->set_to_length(sequence.size());
                e->set_sequence(sequence);
                // Since the softclip consumes no graph, we place it on the node we are going to.
                *m->mutable_position() = default_position;
                
#ifdef debug
                cerr << "New best alignment: " << pb2json(best_path) << " score " << best_score << endl;
#endif
            }
        } else {
            // This path has bases in it

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
            
            
            if (use_xdrop_for_tails) {
#ifdef debug
                Alignment clone = current_alignment;
                get_regular_aligner()->align_pinned(clone, subgraph, subgraph.get_topological_order(), true);
#endif
                get_regular_aligner()->get_xdrop()->align_pinned(current_alignment, subgraph, subgraph.get_topological_order(), true);
#ifdef debug
                cerr << "Xdrop: " << pb2json(current_alignment) << endl;
                cerr << "Normal: " << pb2json(clone) << endl;
#endif
            } else {
                get_regular_aligner()->align_pinned(current_alignment, subgraph, subgraph.get_topological_order(), true);
            }
            
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

    // We really should have gotten something
    assert(best_path.mapping_size() != 0);
    
    return make_pair(best_path, best_score);
}

void MinimizerMapper::align_to_local_haplotypes(const Alignment& aln, const vector<GaplessExtension>& extended_seeds, Alignment& out) const {
    
#ifdef debug
    cerr << "Aligning to haplotypes" << endl;
#endif
    
    // Find all the handles/nodes that are relevant.
    // This holds handles on the strand visited by each relevant extension.
    unordered_set<handle_t> start_points;
    for (auto& extension : extended_seeds) {
        for (auto& core_handle : extension.path) {
            // Collect together all the handles that extensions touch
            start_points.insert(core_handle);
        }
    }
    
    // Walk a stranded Dijkstra traversal the appropriate distance from all the starting points, both upstream and downstream. 
    size_t distance_limit = aln.sequence().size() * 2;
    
#ifdef debug
    cerr << "Search out " << distance_limit << " bp from " << start_points.size() << " handles" << endl;
    for (auto& start : start_points) {
        cerr << "\t" << gbwt_graph.get_id(start) << " " << gbwt_graph.get_is_reverse(start) << endl;
    }
#endif
    
    // We track the handles we get from the search.
    unordered_set<handle_t> context;
    
    // When we reach something going left or right, mark it in the
    // single-direction context and its ID in the total context.
    auto reached_callback = [&](const handle_t& reached, size_t distance) -> bool {
        if (distance > distance_limit) {
            // We hit the limit, so stop the search.
            return false;
        }
#ifdef debug
        cerr << "\tAdd " << gbwt_graph.get_id(reached) << " " << gbwt_graph.get_is_reverse(reached) << " to context" << endl;
#endif
        context.insert(reached);
        return true;
    };
    
    // Do the right pass
    algorithms::dijkstra(&gbwt_graph, start_points, reached_callback, false);
    
#ifdef debug
    cerr << "Got right context of " << context.size() << endl;
#endif
    
    // Pull out all the boundary nodes of the right-walking context graph that
    // have contact with the outside world, or that have no edges.
    // When looking for haplotypes, we have to walk left from them.
    unordered_set<handle_t> right_boundaries;
    for (auto& handle : context) {
        // See if each handle has anything to its right that is not the context.
        // If so, we stop early, and the handle is a boundary in the context graph.
        bool is_enterable = false;
        bool is_tip = true;
        gbwt_graph.follow_edges(handle, false, [&](const handle_t& neighbor) {
            is_tip = false;
            if (!context.count(neighbor)) {
                is_enterable = true;
                return false;
            }
            return true;
        });
        
        if (is_tip || is_enterable) {
            right_boundaries.insert(handle);
        }
    }
    context.clear();
    
#ifdef debug
    cerr << "Got " << right_boundaries.size() << " right boundaries" << endl;
    for (auto& boundary : right_boundaries) {
        cerr << "\t" << gbwt_graph.get_id(boundary) << " " << gbwt_graph.get_is_reverse(boundary) << endl;
    }
#endif
    
    // And the left pass
    algorithms::dijkstra(&gbwt_graph, start_points, reached_callback, true);

#ifdef debug
    cerr << "Got left context of " << context.size() << endl;
#endif
    
    // Pull out all the boundary nodes of the left-walking context graph.
    // When looking for haplotypes, we have to walk right from them.
    unordered_set<handle_t> left_boundaries;
    for (auto& handle : context) {
        bool is_enterable = false;
        bool is_tip = true;
        gbwt_graph.follow_edges(handle, true, [&](const handle_t& neighbor) {
            is_tip = false;
            if (!context.count(neighbor)) {
                is_enterable = true;
                return false;
            }
            return true;
        });
        
        if (is_tip || is_enterable) {
            left_boundaries.insert(handle);
        }
    }
    context.clear();
    
#ifdef debug
    cerr << "Got " << left_boundaries.size() << " left boundaries" << endl;
    for (auto& boundary : left_boundaries) {
        cerr << "\t" << gbwt_graph.get_id(boundary) << " " << gbwt_graph.get_is_reverse(boundary) << endl;
    }
#endif
    
    // Now we have the left and right boundaries, so we can find the haplotypes.
    vector<Path> haplotypes;
    
    for (auto& left_boundary : left_boundaries) {
        // For each left boundary
        
        // Explore the GBWT forward until we hit a right boundary.
        // We want to go until we hit the opposing boundary or run out of paths.
        Position boundary_start = make_position(gbwt_graph.get_id(left_boundary), gbwt_graph.get_is_reverse(left_boundary), 0);
        explore_gbwt(boundary_start, numeric_limits<size_t>::max(), [&](const ImmutablePath& path_to, const handle_t& here) -> bool {
            // When we actually touch something
            
#ifdef debug
            cerr << "Visit " << gbwt_graph.get_id(here) << " " << gbwt_graph.get_is_reverse(here) << endl;
#endif
            
            if (right_boundaries.count(here)) {
                // This is a right boundary.
                
#ifdef debug
                cerr << "\tReached right boundary" << endl;
#endif
                
                // Complete the path
                Mapping m;
                m.mutable_position()->set_node_id(gbwt_graph.get_id(here));
                m.mutable_position()->set_is_reverse(gbwt_graph.get_is_reverse(here));
                Edit* e = m.add_edit();
                e->set_from_length(gbwt_graph.get_length(here));
                e->set_to_length(gbwt_graph.get_length(here));
                
                ImmutablePath path_through = path_to.push_front(m);
                
                // Convert to a Path and save
                haplotypes.emplace_back(to_path(path_through));
                
                // Don't go past here
                return false;
                
            } else {
                // Keep going
                return true;
            }
        }, [&](const ImmutablePath&) {
            // When we hit the length limit or a dead end, do nothing.
            
#ifdef debug
            cerr << "Ran out of distance or hit a dead end" << endl;
#endif
            
        });
        
        // If we actually reach a boundary, that's a Path, so keep it
    }
    
#ifdef debug
    cerr << "Got " << haplotypes.size() << " haplotypes" << endl;
#endif
        
    // Then if there are any remaining right boundaries, we don't actually care, because they aren't in haplotypes with any left boundaries.
    // Nor do we care about left boundaries that couldn't get anywhere.
    // Because they are boundaries, boundaries must be on some haplotype that exits the local graph region through them.
    
    // Align to each haplotype path we found
    
    // Track how much DP is done when operating on haplotypes
    vector<double> haplotype_dp_areas;
    
    // Now look for the best alignment to any path.
    // We set thid to ~-inf because we even want results with 0 score.
    int64_t best_score = numeric_limits<int64_t>::min();
    // we will store the winner in out's path, so set the rest of out's parameters.

    for (auto& path : haplotypes) {
        // For each path, make sure it isn't empty.
        assert(path.mapping_size() > 0);

        // Make a subgraph.
        // TODO: don't copy the path
        PathSubgraph subgraph(&gbwt_graph, path);
        
        // Do global alignment to the path subgraph
        Alignment path_alignment;
        path_alignment.set_sequence(aln.sequence());

#ifdef debug
        cerr << "Align " << pb2json(path_alignment) << " local";

#ifdef debug_dump_graph
        cerr << " vs:" << endl;
        cerr << "Defining path: " << pb2json(path) << endl;
        subgraph.for_each_handle([&](const handle_t& here) {
            cerr << subgraph.get_id(here) << " len " << subgraph.get_length(here)
                << " (" << subgraph.get_sequence(here) << "): " << endl;
            subgraph.follow_edges(here, true, [&](const handle_t& there) {
                cerr << "\t" << subgraph.get_id(there) << " len " << subgraph.get_length(there)
                    << " (" << subgraph.get_sequence(there) << ") ->" << endl;
            });
            subgraph.follow_edges(here, false, [&](const handle_t& there) {
                cerr << "\t-> " << subgraph.get_id(there) << " len " << subgraph.get_length(there)
                    << " (" << subgraph.get_sequence(there) << ")" << endl;
            });
        });
#else
        cerr << endl;
#endif
#endif

        // Do a local alignment with traceback but no score matrix printing.
        get_regular_aligner()->align(path_alignment, subgraph, true, false);
        
#ifdef debug
        cerr << "\tScore: " << path_alignment.score() << endl;
#endif

        if (path_alignment.score() > best_score) {
            // This is a new best alignment. Translate from subgraph into base graph and keep it
            *out.mutable_path() = subgraph.translate_down(path_alignment.path());
            
            // Preserve the identity and score too.
            out.set_identity(path_alignment.identity());
            out.set_score(path_alignment.score());
#ifdef debug
            cerr << "\tNew best: " << pb2json(path_alignment) << endl;
#endif
            
            best_score = path_alignment.score();
        }
        
        // Record how much DP we did.
        haplotype_dp_areas.push_back(path_alignment.sequence().size() * path_from_length(path));
    }
    
    // Now the best alignment from any path is in out.
    
    // Annotate it with how many paths we had to align it against.
    set_annotation(out, "haplotype_alignment_paths", (double)haplotypes.size());
    // And how much area we did
    set_annotation(out, "haplotype_dp_areas", haplotype_dp_areas);
}

unordered_map<size_t, unordered_map<size_t, vector<Path>>>
MinimizerMapper::find_connecting_paths(const vector<GaplessExtension>& extended_seeds, size_t read_length) const {

    // Now this will hold, for each extended seed, for each other
    // reachable extended seed, the graph Paths that the
    // intervening sequence needs to be aligned against in the
    // graph.
    unordered_map<size_t, unordered_map<size_t, vector<Path>>> to_return;

    // All the extended seeds are forward in the read. So we index them by start.
    // Maps from handle in the GBWT graph to offset on that orientation that an extension starts at and index of the extension.
    unordered_map<handle_t, vector<pair<size_t, size_t>>> extensions_by_handle;

    // We track which extended seeds are sources (nothing else can reach them)
    unordered_set<size_t> sources;

    for (size_t i = 0; i < extended_seeds.size(); i++) {
        // For each extension

        // Get the handle it is on
        handle_t handle = extended_seeds[i].path.front();

        // Record that this extension starts at this offset along that handle
        extensions_by_handle[handle].emplace_back(extended_seeds[i].offset, i);

        // Assume it is a source
        sources.insert(i);

#ifdef debug
        auto pos = extended_seeds[i].starting_position(gbwt_graph);
        cerr << "Extended seed " << i << " starts on node " << pos.node_id() << " " << pos.is_reverse()
            << " at offset " << pos.offset() << " corresponding to read "
            << extended_seeds[i].core_interval.first << " - " << extended_seeds[i].core_interval.second << endl;
#endif
    }

    for (auto& kv : extensions_by_handle) {
        // Sort everything on the same handle
        std::sort(kv.second.begin(), kv.second.end(), [&](const pair<size_t, size_t>& a, const pair<size_t, size_t>& b) -> bool {
            return a.first < b.first;
        });
    }

    // For each seed in read order, walk out right in the haplotypes by the max length and see what other seeds we encounter.
    // Remember the read bounds and graph Path we found, for later alignment.
    for (size_t i = 0; i < extended_seeds.size(); i++) {
        // For each starting seed

        // Where do we cut the graph just after its end?
        handle_t start_handle = extended_seeds[i].path.back();
        Position cut_pos_graph = extended_seeds[i].tail_position(gbwt_graph);
        // And the read?
        size_t cut_pos_read = extended_seeds[i].core_interval.second;

#ifdef debug
        cerr << "Extended seed " << i << ": cut after on node " << cut_pos_graph.node_id() << " " << cut_pos_graph.is_reverse()
            << " at point " << cut_pos_graph.offset() << " corresponding to read base " << cut_pos_read << endl;
#endif

        assert(cut_pos_graph.offset() >= 0);

        // Decide if we need to actually do GBWT search, or if we can find a destination on the same node we ended on
        bool do_gbwt_search = true;

        // Look on the same graph node.
        // See if we hit any other extensions on this node.
        auto same_node_found = extensions_by_handle.find(start_handle);
        if (same_node_found != extensions_by_handle.end()) {
            // If we have extended seeds on this node
            
            for (auto& next_offset_and_index : same_node_found->second) {
                // Scan them in order.
                // TODO: Skip to after ourselves.
                
                if (extended_seeds[next_offset_and_index.second].core_interval.first >= cut_pos_read &&
                    next_offset_and_index.first >= cut_pos_graph.offset()) { 
                    
                    // As soon as we find one that starts after we end in both the read and the node

                    // Emit a connecting Path 
                    Path connecting;

                    if (next_offset_and_index.first > cut_pos_graph.offset()) {
                        // There actually is intervening graph material
                        Mapping* m = connecting.add_mapping();
                        *m->mutable_position() = cut_pos_graph;
                        Edit* e = m->add_edit();
                        e->set_from_length(next_offset_and_index.first - m->position().offset());
                        e->set_to_length(next_offset_and_index.first - m->position().offset());
                    }

                    // Emit that connection
                    to_return[i][next_offset_and_index.second].emplace_back(std::move(connecting));

                    // Record that the destination is not a source
                    sources.erase(next_offset_and_index.second);

                    // Don't look at any more destinations
                    do_gbwt_search = false;
                    break;
                }
            }
        }

        if (!do_gbwt_search) {
            // Skip the GBWT search from this extended hit and try the next one
            continue;
        }

        // How long should we search? It should be the longest detectable gap plus the remaining sequence.
        size_t search_limit = get_regular_aligner()->longest_detectable_gap(cut_pos_read, read_length) + (read_length - cut_pos_read);

        // Have we found a way to get to any other extended seeds yet?
        bool reachable_extended_seeds = false;

        // Search everything in the GBWT graph right from the end of the start extended seed, up to the limit.
        explore_gbwt(cut_pos_graph, search_limit, [&](const ImmutablePath& here_path, const handle_t& there_handle) -> bool {
            // When we encounter a new handle visited by haplotypes extending off of the last node in a Path

            // See if we hit any other extensions on this next node
            auto found = extensions_by_handle.find(there_handle);
            if (found != extensions_by_handle.end()) {
                // If we do
                
                for (auto& next_offset_and_index : found->second) {
                    // Look at them in order along the node
                    
                    if (extended_seeds[next_offset_and_index.second].core_interval.first >= cut_pos_read) { 
                        // As soon as we find one that starts in the read after our start extended seed ended

                        // Extend the Path to connect to it.
                        // TODO: Make these shared tail lists for better algorithmics
                        ImmutablePath extended = here_path;
                        
                        if (next_offset_and_index.first > 0) {
                            // There is actual material on this new node before the extended seed we have to hit.
                            Mapping m; 
                            m.mutable_position()->set_node_id(gbwt_graph.get_id(there_handle));
                            m.mutable_position()->set_is_reverse(gbwt_graph.get_is_reverse(there_handle));
                            Edit* e = m.add_edit();
                            // Make sure it runs through the last base *before* the extended seed we are going for
                            e->set_from_length(next_offset_and_index.first);
                            e->set_to_length(next_offset_and_index.first);
                            
                            extended = extended.push_front(m);
                        }

                        // And emit that connection
                        to_return[i][next_offset_and_index.second].emplace_back(to_path(extended));
                        reachable_extended_seeds = true;

                        // Record that the destination is not a source
                        sources.erase(next_offset_and_index.second);

                        // Don't look at any more destinations on this node, or
                        // any extensions past this node of this search state.
                        return false;
                    }
                }
            }

            // Otherwise we didn't hit anything we can stop at. Keep extending.
            return true;
        }, [&](const ImmutablePath& limit_path) {
            // When we blow past the walk distance limit or hit a dead end
            
            if (cut_pos_read < read_length && !reachable_extended_seeds && linear_tails) {
                // We have sequence to align and a way to escape and align it,
                // and nowhere else we know of (yet) to go with it, and we want
                // tails.
                
                // Save that as a tail path.
                to_return[i][numeric_limits<size_t>::max()].emplace_back(to_path(limit_path));
                // If we end up with paths anywhere else after all, we will destroy it, so we
                // will only keep it for sinks.
            }
        });
        
        if (reachable_extended_seeds && linear_tails) {
            // Make sure that if we can go anywhere else we *don't* consider wandering off to nowhere.
            auto found = to_return[i].find(numeric_limits<size_t>::max());
            if (to_return[i].size() > 1 && found != to_return[i].end()) {
                // We have a going-off-to-nothing path and also a path to somewhere else.
                // Don't go off to nothing.
                to_return[i].erase(found);
            }
        }
    }

#ifdef debug
    cerr << "After rightward extensions, have " << sources.size() << " sources" << endl;
#endif

    if (!linear_tails) {
        // We have all the non-tail paths, which is what we need
        return to_return;
    }

    // Otherwise we want tails.

    // We need the paths *from* numeric_limits<size_t>::max() to sources.
    // Luckily we know the sources.
    for (const size_t& i : sources) {
        // For each source
        
#ifdef debug
        cerr << "Extended seed " << i << " is a source" << endl;
#endif
        
        if (extended_seeds[i].core_interval.first > 0) {
#ifdef debug
            cerr << "\tIt is not at the start of the read, so there is a left tail" << endl;
#endif

            // Find its start
            Position start = extended_seeds[i].starting_position(gbwt_graph);
            
#ifdef debug
            cerr << "\tPosition read-forward to search left before: " << pb2json(start) << endl;
#endif

            // Flip it around to face left
            start = reverse(start, gbwt_graph.get_length(gbwt_graph.get_handle(start.node_id())));
            
#ifdef debug
            cerr << "\tPosition read-reverse to search right after: " << pb2json(start) << endl;
#endif

            // Now the search limit is all the read *before* the seed, plus the detectable gap
            size_t search_limit = get_regular_aligner()->longest_detectable_gap(read_length, extended_seeds[i].core_interval.first) +
                extended_seeds[i].core_interval.first;

            // Start another search, but going left.
            explore_gbwt(start, search_limit, [&](const ImmutablePath& here_path, const handle_t& there_handle) -> bool {
                // If we weren't reachable from anyone, nobody should be reachable from us going the other way.
                // So always keep going.
                return true;
            }, [&](const ImmutablePath& limit_path) {
                // We have this path going right from start and hitting the walk limit or the edge of the graph.

                for (auto& mapping : limit_path) {
                    // Make sure nothing has a negative offset before flipping.
                    assert(mapping.position().offset() >= 0);
                }

                // Convert to Path and flip around
                Path flipped = reverse_complement_path(to_path(limit_path), [&](id_t id) -> size_t {
                    return gbwt_graph.get_length(gbwt_graph.get_handle(id));
                });
                
                for (auto& mapping : flipped.mapping()) {
                    // Make sure nothing has a negative offset after flipping.
                    assert(mapping.position().offset() >= 0);
                }

                // Record that as a path from numeric_limits<size_t>::max() to i.
                to_return[numeric_limits<size_t>::max()][i].emplace_back(std::move(flipped));
            });
            
            assert(to_return[numeric_limits<size_t>::max()][i].size() > 0);
        }
    }
    
    // Now this should be filled in with all the connectivity, so return.
    return to_return;
    
}

unordered_map<size_t, vector<TreeSubgraph>> MinimizerMapper::get_tail_forests(const vector<GaplessExtension>& extended_seeds,
    size_t read_length, const unordered_map<size_t, unordered_map<size_t, vector<Path>>>& connecting_paths, bool left_tails) const {

    // We will fill this in with all the trees we return, by parent extension.
    unordered_map<size_t, vector<TreeSubgraph>> to_return;

    // First, find all the source/sink extensions we actually want to do.
    unordered_set<size_t> tail_havers;
    if (left_tails) {
        for (size_t i = 0; i < extended_seeds.size(); i++) {
            // We have a left tail if nothing in connecting_paths comes to us.
        
            // So everything has left tails to start with
            if (extended_seeds[i].core_interval.first != 0) {
                // As long as it has some read before it
                tail_havers.insert(i);
#ifdef debug
                cerr << "Extension " << i << " running " << extended_seeds[i].core_interval.first
                    << " - " << extended_seeds[i].core_interval.second << " may have a left tail" << endl;
#endif
            }
        }
        
        for (auto& kv : connecting_paths) {
            for (auto& dest_and_path : kv.second) {
                // And then we remove things when they are visited.
                // Note that we assume paths to/from
                // numeric_limits<size_t>::max() (i.e. tail paths) don't
                // appear.
                tail_havers.erase(dest_and_path.first);
            }
        }
    } else {
        for (size_t i = 0; i < extended_seeds.size(); i++) {
            // We might have a right tail if we go nowhere in connecting_paths
            
            auto found = connecting_paths.find(i);
            if ((found == connecting_paths.end() || found->second.empty()) && extended_seeds[i].core_interval.second < read_length) {
                // So if we go nowhere and actually have bases after us, add us
                tail_havers.insert(i);
                
#ifdef debug
                cerr << "Extension " << i << " running " << extended_seeds[i].core_interval.first
                    << " - " << extended_seeds[i].core_interval.second << " may have a right tail in read of length " << read_length << endl;
#endif
                
            }
        }
    }
    
    for (auto& extension_number : tail_havers) {
        // Now for each extension that can have tails, walk the GBWT in the appropriate direction
        
        // TODO: Come up with a better way to do this with more accessors on the extension and less get_handle
        // Get the Position reading out of the extnsion on the appropriate tail
        Position from;
        // And the length of that tail
        size_t tail_length;
        if (left_tails) {
            // Look right from start 
            from = extended_seeds[extension_number].starting_position(gbwt_graph);
            // And then flip to look the other way at the prev base
            from = reverse(from, gbwt_graph.get_length(gbwt_graph.get_handle(from.node_id(), false)));
            
            tail_length = extended_seeds[extension_number].core_interval.first;
        } else {
            // Look right from end
            from = extended_seeds[extension_number].tail_position(gbwt_graph);
            
            tail_length = read_length - extended_seeds[extension_number].core_interval.second;
        }
        
        // This is one tree that we are filling in
        vector<pair<int64_t, handle_t>> tree;
        
        // This is a stack of indexes at which we put parents in the tree
        list<int64_t> parent_stack;
        
        // Get the handle we are starting from
        handle_t start_handle = gbwt_graph.get_handle(from.node_id(), from.is_reverse());
        
        // Decide if the start node will end up included in the tree, or if we cut it all off with the offset.
        bool start_included = (from.offset() < gbwt_graph.get_length(start_handle));
        
        // How long should we search? It should be the longest detectable gap plus the remaining sequence.
        size_t search_limit = get_regular_aligner()->longest_detectable_gap(tail_length, read_length) + tail_length;
        
        // Do a DFS over the haplotypes in the GBWT out to that distance.
        dfs_gbwt(start_handle, from.offset(), search_limit, [&](const handle_t& entered) {
            // Enter a new handle.
            
            if (parent_stack.empty()) {
                // This is the root of a new tree in the forrest
                
                if (!tree.empty()) {
                    // Save the old tree and start a new one.
                    // We need to cut off from.offset() from the root, unless we would cut off the whole root.
                    // In that case, the GBWT DFS will have skipped the empty root entirely, so we cut off nothing.
                    to_return[extension_number].emplace_back(&gbwt_graph, std::move(tree), start_included ? from.offset() : 0);
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
            to_return[extension_number].emplace_back(&gbwt_graph, std::move(tree), start_included ? from.offset() : 0);
            tree.clear();
        }
    }
    
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

void MinimizerMapper::explore_gbwt(const Position& from, size_t walk_distance,
    const function<bool(const ImmutablePath&, const handle_t&)>& visit_callback,
    const function<void(const ImmutablePath&)>& limit_callback) const {
   
    // Get a handle to the node the from position is on, in the position's forward orientation
    handle_t start_handle = gbwt_graph.get_handle(from.node_id(), from.is_reverse());
    
    // Delegate to the handle-based version
    explore_gbwt(start_handle, from.offset(), walk_distance, visit_callback, limit_callback);
    
}

void MinimizerMapper::explore_gbwt(handle_t from_handle, size_t from_offset, size_t walk_distance,
    const function<bool(const ImmutablePath&, const handle_t&)>& visit_callback,
    const function<void(const ImmutablePath&)>& limit_callback) const {
    
#ifdef debug
    cerr << "Exploring GBWT out from " << gbwt_graph.get_id(from_handle) << " " << gbwt_graph.get_is_reverse(from_handle)
        << " + " << from_offset << " to distance " << walk_distance << endl;
#endif
    
    // Holds the gbwt::SearchState we are at, and the ImmutablePath (backward)
    // from the end of the starting seed up through the end of the node we just
    // searched. The from_length of the path tracks our consumption of distance
    // limit.
    using traversal_state_t = pair<gbwt::SearchState, ImmutablePath>;

    // Turn it into a SearchState
    gbwt::SearchState start_state = gbwt_graph.get_state(from_handle);
    
    if (start_state.empty()) {
        // No haplotypes even visit the first node. Have a 0-mapping dead end.
        limit_callback(ImmutablePath());
        return;
    }

    // The search state represents searching through the end of the node, so we have to consume that much search limit.

    // Tack on how much search limit distance we consume by going to the end of
    // the node. Our start position is a cut *between* bases, and we take everything after it.
    // If the cut is at the offset of the whole length of the node, we take 0 bases.
    // If it is at 0, we take all the bases in the node.
    size_t distance_to_node_end = gbwt_graph.get_length(from_handle) - from_offset;    
    
    // And make a Path that represents the part of the node we're on that goes out to the end.
    // This may be empty if the hit already stopped at the end of the node
    ImmutablePath path_to_end;
    if (distance_to_node_end != 0) {
        // We didn't hit the end of the node already.

        // Make a mapping that starts on the right side of the cut we started our search at.
        Mapping m;
        *m.mutable_position() = make_position(gbwt_graph.get_id(from_handle),
            gbwt_graph.get_is_reverse(from_handle), from_offset);
        m.mutable_position()->set_offset(m.position().offset());

        // Make it the requested length of perfect match.
        Edit* e = m.add_edit();
        e->set_from_length(distance_to_node_end);
        e->set_to_length(distance_to_node_end);
        
        // Put it in the list
        path_to_end = path_to_end.push_front(m);
    }
   
#ifdef debug   
    cerr << "Starting traversal with";
    for (auto& mapping : path_to_end) {
        cerr << " " << pb2json(mapping);
    }
    cerr << " from " << gbwt_graph.get_id(from_handle) << " " << gbwt_graph.get_is_reverse(from_handle)
        << " + " << from_offset << endl;
#endif
    
    // Glom these together into a traversal state and queue it up.

    // Holds a queue of search states to extend.
    list<traversal_state_t> queue{{start_state, path_to_end}};
    // Track queue size separately since getting it form the queue is O(n)
    size_t queue_size = 1;
    // Track queue high-water mark for debugging.
    size_t queue_max = 1;
    
    
    // Track how many times we hit the end/limit
    size_t limit_hits = 0;
    // And dead ends sweparately
    size_t dead_ends = 0;

    while (!queue.empty()) {
        // While there are things in the queue
        
        // Record max size
        queue_max = max(queue_max, queue_size);

        // Grab one
        traversal_state_t here(std::move(queue.front()));
        queue.pop_front();
        queue_size--;
        gbwt::SearchState& here_state = here.first;
        const ImmutablePath& here_path = here.second;
        
        // follow_paths on it
        bool got_anywhere = false;
        gbwt_graph.follow_paths(here_state, [&](const gbwt::SearchState& there_state) -> bool {
            if (there_state.empty()) {
                // Ignore places that no haplotypes go, and get the next place instead.
                return true;
            }
            
            // For each place it can go
            handle_t there_handle = gbwt_graph.node_to_handle(there_state.node);
            
#ifdef debug
            handle_t here_handle = gbwt_graph.node_to_handle(here_state.node);
            cerr << "Saw GBWT edge " << gbwt_graph.get_id(here_handle) << " " << gbwt_graph.get_is_reverse(here_handle)
                << " -> " << gbwt_graph.get_id(there_handle) << " " << gbwt_graph.get_is_reverse(there_handle) << endl;
            assert(gbwt_graph.has_edge(here_handle, there_handle));
#endif
            
            // Record that we got there
            got_anywhere = true;

            // Say we can go from here to there. Should we?
            bool continue_extending = visit_callback(here_path, there_handle);

            if (continue_extending) {
                // Generate the path we take if we take up all of that node.
                Mapping m;
                m.mutable_position()->set_node_id(gbwt_graph.get_id(there_handle));
                m.mutable_position()->set_is_reverse(gbwt_graph.get_is_reverse(there_handle));
                Edit* e = m.add_edit();
                e->set_from_length(gbwt_graph.get_length(there_handle));
                e->set_to_length(gbwt_graph.get_length(there_handle));
                
                ImmutablePath extended = here_path.push_front(m);

                // See if we can get to the end of the node without going outside the search length.
                if (immutable_path_from_length(here_path) + gbwt_graph.get_length(there_handle) <= walk_distance) {
                    // If so, continue the search
                    queue.emplace_back(there_state, extended);
                    queue_size++;
                } else {
                    // Report that, with this extension, we hit the limit.
                    limit_callback(extended);
                    limit_hits++;
                }
            }

            // Look at other possible haplotypes from where we came from
            return true;
        });
        
        if (!got_anywhere) {
            // We hit a dead end here. Report that.
            limit_callback(here_path);
            dead_ends++;
        }
    }
    
#ifdef debug
    cerr << "Queue max: " << queue_max << " Limit hits: " << limit_hits << " Dead ends: " << dead_ends << endl;
#endif
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
    
    // Holds the gbwt::SearchState we are at, and the distance we have consumed
    using traversal_state_t = pair<gbwt::SearchState, size_t>;

    // Turn from_handle into a SearchState.
    // TODO: Let a search state come in.
    gbwt::SearchState start_state = gbwt_graph.get_state(from_handle);
    
    if (start_state.empty()) {
        // No haplotypes even visit the first node. Stop.
        return;
    }

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
                
                if (there_state.empty()) {
                    // If it is empty, don't do it
                    return true;
                }
                
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


