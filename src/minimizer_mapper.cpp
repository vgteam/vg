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

MinimizerMapper::MinimizerMapper(const GBWTGraph& graph, const MinimizerIndex& minimizer_index,
     MinimumDistanceIndex& distance_index, const PathPositionHandleGraph* path_graph) :
    path_graph(path_graph), minimizer_index(minimizer_index),
    distance_index(distance_index), gbwt_graph(graph),
    extender(gbwt_graph, *(get_regular_aligner())), clusterer(distance_index) {
    
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
   
    if (track_provenance) {
        // Start the minimizer finding stage
        funnel.stage("minimizer");
    }
    
    // We will find all the seed hits
    vector<pos_t> seeds;
    
    // This will hold all the minimizers in the query
    vector<MinimizerIndex::minimizer_type> minimizers;
    // And either way this will map from seed to minimizer that generated it
    vector<size_t> seed_to_source;
    
    // Find minimizers in the query
    minimizers = minimizer_index.minimizers(aln.sequence());
    
    if (track_provenance) {
        // Record how many we found, as new lines.
        funnel.introduce(minimizers.size());
        
        // Start the minimizer locating stage
        funnel.stage("seed");
    }

    // Compute minimizer scores for all minimizers as 1 + ln(hard_hit_cap) - ln(hits).
    std::vector<double> minimizer_score(minimizers.size(), 0.0);
    double target_score = 0.0;
    for (size_t i = 0; i < minimizers.size(); i++) {
        size_t hits = minimizer_index.count(minimizers[i]);
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

        if (track_provenance) {
            // Say we're working on it
            funnel.processing_input(minimizer_num);
        }

        // Select the minimizer if it is informative enough or if the total score
        // of the selected minimizers is not high enough.
        size_t hits = minimizer_index.count(minimizers[minimizer_num]);

        //Take all the seeds

        // Locate the hits.
        for (auto& hit : minimizer_index.find(minimizers[minimizer_num])) {
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
        
        if (track_provenance) {
            // Record in the funnel that this minimizer gave rise to these seeds.
            funnel.expand(minimizer_num, hits);
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
    cerr << "Read " << aln.name() << ": " << aln.sequence() << endl;
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

    for (size_t i = 0; i < clusters.size(); i++) {
        // For each cluster
        auto& cluster = clusters[i];
        
        if (track_provenance) {
            // Say we're making it
            funnel.producing_output(i);
        }

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
        
        if (track_provenance) {
            // Record the cluster in the funnel as a group of the size of the number of items.
            funnel.merge_group(cluster.begin(), cluster.end());
            funnel.score(funnel.latest(), cluster_score[i]);
            
            // Say we made it.
            funnel.produced_output();
        }

        //TODO:
        //Get the cluster coverage
        // We set bits in here to true when query anchors cover them
        sdsl::bit_vector covered(aln.sequence().size(), 0);
        std::uint64_t k_bit_mask = sdsl::bits::lo_set[minimizer_index.k()];

        for (auto hit_index : cluster) {
            // For each hit in the cluster, work out what anchor sequence it is from.
            size_t source_index = seed_to_source[hit_index];

            // The offset of a reverse minimizer is the endpoint of the kmer
            size_t start_offset = minimizers[source_index].offset;
            if (minimizers[source_index].is_reverse) {
                start_offset = start_offset + 1 - minimizer_index.k();
            }

            // Set the k bits starting at start_offset.
            covered.set_int(start_offset, k_bit_mask, minimizer_index.k());
        }

        // Count up the covered positions
        size_t covered_count = sdsl::util::cnt_one_bits(covered);

        // Turn that into a fraction
        read_coverage_by_cluster.push_back(covered_count / (double) covered.size());


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
    std::sort(cluster_indexes_in_order.begin(), cluster_indexes_in_order.end(), 
        [&](const size_t& a, const size_t& b) -> bool {
            // Return true if a must come before b, and false otherwise
            return (read_coverage_by_cluster[a] > read_coverage_by_cluster[b]);
    });

    //Retain clusters only if their read coverage is better than this
    double cluster_coverage_cutoff = cluster_indexes_in_order.size() == 0 ? 0 : 
                                 read_coverage_by_cluster[cluster_indexes_in_order[0]]
                                    - cluster_coverage_threshold;
    //Retain clusters only if their score is better than this
    double cluster_score_cutoff = cluster_score.size() == 0 ? 0 :
                    *std::max_element(cluster_score.begin(), cluster_score.end()) - cluster_score_threshold;
    
    if (track_provenance) {
        // Now we go from clusters to gapless extensions
        funnel.stage("extend");
    }
    
    // These are the GaplessExtensions for all the clusters, in cluster_indexes_in_order order.
    vector<vector<GaplessExtension>> cluster_extensions;
    cluster_extensions.reserve(cluster_indexes_in_order.size());
    
    size_t num_extensions = 0;
    for (size_t i = 0; i < clusters.size(); i++) {
        // For each cluster, in sorted order
        size_t& cluster_num = cluster_indexes_in_order[i];
        num_extensions ++;
        
        if (track_provenance) {
            funnel.processing_input(cluster_num);
        }

        vector<size_t>& cluster = clusters[cluster_num];

#ifdef debug
        cerr << "Cluster " << cluster_num << " rank " << i << ": " << endl;
#endif
         
        // Pack the seeds for GaplessExtender.
        GaplessExtender::cluster_type seed_matchings;
        for (auto& seed_index : cluster) {
            // Insert the (graph position, read offset) pair.
            seed_matchings.insert(GaplessExtender::to_seed(seeds[seed_index], minimizers[seed_to_source[seed_index]].offset));
#ifdef debug
            cerr << "Seed read:" << minimizers[seed_to_source[seed_index]].offset << " = " << seeds[seed_index]
                << " from minimizer " << seed_to_source[seed_index] << "(" << minimizer_index.count(minimizers[seed_to_source[seed_index]]) << ")" << endl;
#endif
        }
        
        // Extend seed hits in the cluster into one or more gapless extensions
        vector<GaplessExtension> extensions = extender.extend(seed_matchings, aln.sequence());
        // Find the best scoring extension
        vector<GaplessExtension> filtered_extensions;
        int32_t best_extension_score = 0;
        for (GaplessExtension& extension : extensions) {
            best_extension_score = max(best_extension_score, extension.score);
        }
        //Keep only the extensions whose score is within extension_score_threshold
        //of the best scoring extension
        for (GaplessExtension& extension : extensions) {
            filtered_extensions.push_back(std::move(extension));
        }
        cluster_extensions.emplace_back(std::move(filtered_extensions));
        
        if (track_provenance) {
            // Record with the funnel that the previous group became a group of this size.
            // Don't bother recording the seed to extension matching...
            funnel.project_group(cluster_num, cluster_extensions.back().size());
            
            // Say we finished with this cluster, for now.
            funnel.processed_input();
        }
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
        // Count the matches suggested by the group and use that as a score.
        cluster_extension_scores.push_back(estimate_extension_group_score(aln, extensions));
        
        if (track_provenance) {
            // Record the score with the funnel
            funnel.score(i, cluster_extension_scores.back());
            funnel.produced_output();
        }
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

    //Retain cluster_extensions only if their score (coverage of the read) is at
    //least as good as this
    double extension_set_cutoff = cluster_extension_scores.size() == 0 ? 0 :
                              cluster_extension_scores[extension_indexes_in_order[0]]
                                    - extension_set_score_threshold;
    
    if (track_provenance) {
        funnel.stage("align");
    }
    
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
    int second_best_score = cluster_extension_scores.size() < 2 ? 0 :
                            cluster_extension_scores[extension_indexes_in_order[1]];
    for (size_t i = 0; i < extension_indexes_in_order.size() && i < max_alignments ; i++) {
        // Find the extension group we are talking about
        size_t& extension_num = extension_indexes_in_order[i];
        
        if (track_provenance) {
            funnel.processing_input(extension_num);
        }

        auto& extensions = cluster_extensions[extension_num];
        
        if (true) {
            // Always take the first and second.
            // For later ones, check if this score is significant relative to the running best and second best scores.
            
            // If so, get an Alignment out of it somehow, and throw it in.
            alignments.emplace_back(aln);
            Alignment& out = alignments.back();
            
            if (extensions.size() == 1 && extensions[0].full()) {
                // We got a full-length extension, so directly convert to an Alignment.
                
                if (track_provenance) {
                    funnel.substage("direct");
                }

                *out.mutable_path() = extensions.front().to_path(gbwt_graph, out.sequence());
                
                // The score estimate is exact.
                int alignment_score = cluster_extension_scores[extension_num];
                
                // Compute identity from mismatch count.
                size_t mismatch_count = extensions[0].mismatches();
                double identity = out.sequence().size() == 0 ? 0.0 : (out.sequence().size() - mismatch_count) / (double) out.sequence().size();
                
                // Fill in the score and identity
                out.set_score(alignment_score);
                out.set_identity(identity);
                
                if (track_provenance) {
                    // Stop the current substage
                    funnel.substage_stop();
                }
            } else if (do_dp) {
                // We need to do chaining.
                
                if (track_provenance) {
                    funnel.substage("chain");
                }
                
                // Do the DP and compute alignment into out 
                find_optimal_tail_alignments(aln, extensions, out);
                
                if (track_provenance) {
                    // We're done chaining. Next alignment may not go through this substage.
                    funnel.substage_stop();
                }
            } else {
                // We would do chaining but it is disabled.
                // Leave out unaligned
            }
            
            
            if (track_provenance) {
                // Record the Alignment and its score with the funnel
                funnel.project(extension_num);
                funnel.score(i, out.score());
                
                // We're done with this input item
                funnel.processed_input();
            }
        } else {
            // If this score is insignificant, nothing past here is significant.
            // Don't do any more.
            
            if (track_provenance) {
                funnel.kill_all(extension_indexes_in_order.begin() + i, extension_indexes_in_order.end());
                funnel.processed_input();
            }
            
            break;
        }
    }
    
    if (alignments.size() == 0) {
        // Produce an unaligned Alignment
        alignments.emplace_back(aln);
        
        if (track_provenance) {
            // Say it came from nowhere
            funnel.introduce();
        }
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
    
    if (track_provenance) {
        // Now say we are finding the winner(s)
        funnel.stage("winner");
    }
    
    vector<Alignment> mappings;
    mappings.reserve(min(alignments_in_order.size(), max_multimaps));
    for (size_t i = 0; i < alignments_in_order.size() && i < max_multimaps; i++) {
        // For each output slot, fill it with the alignment at that rank if available.
        size_t& alignment_num = alignments_in_order[i];
        mappings.emplace_back(std::move(alignments[alignment_num]));
        
        if (track_provenance) {
            // Tell the funnel
            funnel.project(alignment_num);
            funnel.score(alignment_num, mappings.back().score());
        }
    }

    if (track_provenance) {
        if (max_multimaps < alignments_in_order.size()) {
            // Some things stop here
            funnel.kill_all(alignments_in_order.begin() + max_multimaps, alignments_in_order.end());
        }
        
        funnel.substage("mapq");
    }
    
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
    double mapq = get_regular_aligner()->compute_mapping_quality(scores, false);

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
    
        // And with the number of results in play at each stage
        funnel.for_each_stage([&](const string& stage, const vector<size_t>& result_sizes) {
            // Save the number of items
            set_annotation(mappings[0], "stage_" + stage + "_results", (double)result_sizes.size());
            // Save the size of each item
            vector<double> converted;
            converted.reserve(result_sizes.size());
            std::copy(result_sizes.begin(), result_sizes.end(), std::back_inserter(converted));
            set_annotation(mappings[0], "stage_" + stage + "_sizes", converted);
        });

        if (track_correctness) {
            // And with the last stage at which we had any descendants of the correct seed hit locations
            set_annotation(mappings[0], "last_correct_stage", funnel.last_correct_stage());
        }
    }
    
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
        // This is a full length match. We already have the score.
        return extended_seeds.front().score;
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
                next_seed_start = extended_seeds[unentered].read_interval.first;
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
            
            while (unentered < extended_seeds.size() && extended_seeds[unentered].read_interval.first == sweep_line) {
                // Bring in anything that starts here
                end_heap.emplace_back(extended_seeds[unentered].read_interval.second, unentered);
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

void MinimizerMapper::find_optimal_tail_alignments(const Alignment& aln, const vector<GaplessExtension>& extended_seeds, Alignment& out) const {

#ifdef debug
    cerr << "Trying again find tail alignments for " << extended_seeds.size() << " extended seeds" << endl;
#endif

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
            aln.sequence().begin() + extended_seed.read_interval.first));
        // The position in the read it occurs at will be handled by the multipath topology.
        if (extended_seed.read_interval.first == 0) {
            // But if it occurs at the very start of the read we need to mark that now.
            mp.add_start(mp.subpath_size() - 1);
        }
    }
    
    // Handle left tails as a forest of trees
   
    // Get the forests of all left tails by extension they belong to
    auto tails_by_extension = get_tail_forests(extended_seeds, aln.sequence().size(), true);
   
    for (auto& kv : tails_by_extension) {
        // For each extension that has a nonempty tail
        auto& source = kv.first;
        auto& forest = kv.second;
   
        // Grab the part of the read sequence that comes before it
        string before_sequence = aln.sequence().substr(0, extended_seeds[source].read_interval.first);
        
        // Do right-pinned alignment
        pair<Path, int64_t> result = get_best_alignment_against_any_tree(forest, before_sequence,
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

    // We must have somewhere to start.
    assert(mp.start_size() > 0);

    // Handle right tails as a forest of trees

    // Get the forests of all right tails by extension they belong to
    tails_by_extension = get_tail_forests(extended_seeds, aln.sequence().size(), false);
    
    for (auto& kv : tails_by_extension) {
        // For each extension that has a nonempty tail
        auto& from = kv.first;
        auto& forest = kv.second;
        
#ifdef debug
        cerr << "Consider right tails for extension " << from << " with interval "
            << extended_seeds[from].read_interval.first << " - " << extended_seeds[from].read_interval.second << endl;
#endif
        
        // Find the sequence
        string trailing_sequence = aln.sequence().substr(extended_seeds[from].read_interval.second);
        
        // There should be actual trailing sequence to align on this escape path
        assert(!trailing_sequence.empty());
        
        // Do left-pinned alignment
        pair<Path, int64_t> result = get_best_alignment_against_any_tree(forest, trailing_sequence,
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
            cerr << "\tPath: " << pb2json(path) << endl;
#else
            cerr << endl;
#endif
#endif
            
            // Align, accounting for full length bonus
            
            if (pinned) {
                get_regular_aligner()->get_xdrop()->align_pinned(current_alignment, subgraph, subgraph.get_topological_order(), pin_left);
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

unordered_map<size_t, vector<TreeSubgraph>> MinimizerMapper::get_tail_forests(const vector<GaplessExtension>& extended_seeds,
    size_t read_length, bool left_tails) const {

    // We will fill this in with all the trees we return, by parent extension.
    unordered_map<size_t, vector<TreeSubgraph>> to_return;

    for (size_t extension_number = 0; extension_number < extended_seeds.size(); extension_number++) {
        // Now for each extension that can have tails, walk the GBWT in the appropriate direction
        
#ifdef debug
        cerr << "Look for " << (left_tails ? "left" : "right") << " tails from extension " << extension_number << endl;
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
            from = extended_seeds[extension_number].starting_position(gbwt_graph);
            // And then flip to look the other way at the prev base
            from = reverse(from, gbwt_graph.get_length(gbwt_graph.get_handle(from.node_id(), false)));
           
            // Use the search state going backward
            base_state = &extended_seeds[extension_number].state.backward;
           
            tail_length = extended_seeds[extension_number].read_interval.first;
        } else {
            // Look right from end
            from = extended_seeds[extension_number].tail_position(gbwt_graph);
            
            // Use the search state going forward
            base_state = &extended_seeds[extension_number].state.forward;
            
            tail_length = read_length - extended_seeds[extension_number].read_interval.second;
        }
    
        if (tail_length == 0) {
            // Don't go looking for places to put no tail.
            continue;
        }
    
        // Make sure we at least have an empty forest if we have any tail.
        to_return.emplace(extension_number, vector<TreeSubgraph>());
        
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
        
#ifdef debug
        if (to_return.count(extension_number)) {
            cerr << "Found " << to_return[extension_number].size() << " trees" << endl;
        } else {
            cerr << "Found no trees in no forest" << endl;
        }
#endif
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


