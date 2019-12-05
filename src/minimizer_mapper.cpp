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

MinimizerMapper::MinimizerMapper(const gbwtgraph::GBWTGraph& graph, const gbwtgraph::DefaultMinimizerIndex& minimizer_index,
    MinimumDistanceIndex& distance_index, const PathPositionHandleGraph* path_graph) :
    path_graph(path_graph), minimizer_index(minimizer_index),
    distance_index(distance_index), gbwt_graph(graph),
    extender(gbwt_graph, *(get_regular_aligner())), clusterer(distance_index) {
    
    // Nothing to do!
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
    vector<gbwtgraph::DefaultMinimizerIndex::minimizer_type> minimizers;
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
    double base_target_score = 0.0;
    for (size_t i = 0; i < minimizers.size(); i++) {
        size_t hits = minimizer_index.count(minimizers[i]);
        if (hits > 0) {
            if (hits <= hard_hit_cap) {
                minimizer_score[i] = 1.0 + std::log(hard_hit_cap) - std::log(hits);
            } else {
                minimizer_score[i] = 1.0;
            }
        }
        base_target_score += minimizer_score[i];
    }
    double target_score = (base_target_score * minimizer_score_fraction) + 0.000001;

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
        
#ifdef debug
        cerr << "Minimizer " << minimizer_num << " = " << minimizers[minimizer_num].key.decode(minimizer_index.k())
            << " has " << hits << " hits" << endl;
#endif
        
        if (hits == 0) {
            // A minimizer with no hits can't go on.
            if (track_provenance) {
                funnel.fail("any-hits", minimizer_num);
            }
        } else if (hits <= hit_cap || (hits <= hard_hit_cap && selected_score + minimizer_score[minimizer_num] <= target_score)) {
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
                funnel.pass("any-hits", minimizer_num);
                funnel.pass("hard-hit-cap", minimizer_num);
                funnel.pass("hit-cap||score-fraction", minimizer_num, selected_score  / base_target_score);
                funnel.expand(minimizer_num, hits);
            }
        } else if (hits <= hard_hit_cap) {
            // Passed hard hit cap but failed score fraction/normal hit cap
            rejected_count++;
            if (track_provenance) {
                funnel.pass("any-hits", minimizer_num);
                funnel.pass("hard-hit-cap", minimizer_num);
                funnel.fail("hit-cap||score-fraction", minimizer_num, (selected_score + minimizer_score[minimizer_num]) / base_target_score);
            }
        } else {
            // Failed hard hit cap
            rejected_count++;
            if (track_provenance) {
                funnel.pass("any-hits", minimizer_num);
                funnel.fail("hard-hit-cap", minimizer_num);
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
                                    
    // Retain clusters only if their score is better than this, in addition to the coverage cutoff
    double cluster_score_cutoff = cluster_score.size() == 0 ? 0 :
                    *std::max_element(cluster_score.begin(), cluster_score.end())
                                    - cluster_score_threshold;
    
    if (track_provenance) {
        // Now we go from clusters to gapless extensions
        funnel.stage("extend");
    }
    
    // These are the GaplessExtensions for all the clusters, in cluster_indexes_in_order order.
    vector<vector<GaplessExtension>> cluster_extensions;
    cluster_extensions.reserve(clusters.size());
    //For each cluster, what fraction of "equivalent" clusters did we keep?
    vector<double> probability_cluster_lost;
    //What is the score and coverage we are considering and how many reads
    size_t curr_coverage = 0;
    size_t curr_score = 0;
    size_t curr_kept = 0;
    size_t curr_count = 0;
    
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
            } else if (!read_coverage_by_cluster[cluster_num] == curr_coverage ||
                    !cluster_score[cluster_num] == curr_score) {
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
                return false;
            }
            

            //Only keep this cluster if we have few enough equivalent clusters

            vector<size_t>& cluster = clusters[cluster_num];

#ifdef debug
            cerr << "Cluster " << cluster_num << endl;
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
            cluster_extensions.emplace_back(std::move(extender.extend(seed_matchings, aln.sequence())));
            
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
        
        vector<GaplessExtension>& extensions = cluster_extensions[i];
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
    //probability_cluster_lost but ordered by alignment
    vector<double> probability_alignment_lost;
    probability_alignment_lost.reserve(cluster_extensions.size());

    
    // Clear any old refpos annotation and path
    aln.clear_refpos();
    aln.clear_path();
    aln.set_score(0);
    aln.set_identity(0);
    aln.set_mapping_quality(0);
    
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
                //If there is a second extension and its score is at least half of the best score
                alignments.push_back(std::move(second_best_alignment));
                probability_alignment_lost.push_back(probability_cluster_lost[extension_num]);

                if (track_provenance) {
    
                    funnel.project(extension_num);
                    funnel.score(alignments.size() - 1, alignments.back().score());
                    // We're done with this input item
                    funnel.processed_input();
                }
            }

            alignments.push_back(std::move(best_alignment));
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
    
    // Fill this in with the alignments we will output
    vector<Alignment> mappings;
    mappings.reserve(min(alignments.size(), max_multimaps));
    
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
    
    if (probability_mapping_lost.front() > 0) {
        mapq = min(mapq,round(prob_to_phred(probability_mapping_lost.front())));
    }
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
        set_annotation(mappings[0], "param_cluster-coverage", (double) cluster_coverage_threshold);
        set_annotation(mappings[0], "param_extension-set", (double) extension_set_score_threshold);
        set_annotation(mappings[0], "param_max-multimaps", (double) max_multimaps);
    }
    
    // Ship out all the aligned alignments
    alignment_emitter.emit_mapped_single(std::move(mappings));

#ifdef debug
    // Dump the funnel info graph.
    funnel.to_dot(cerr);
#endif
}

void MinimizerMapper::map_paired(Alignment& aln1, Alignment& aln2, AlignmentEmitter& alignment_emitter) {
    // For each input alignment
    
#ifdef debug
    cerr << "Read " << aln.name() << ": " << aln.sequence() << endl;
#endif

    // Make two new funnel instrumenters to watch us map this read pair.
    vector<Funnel> funnels;
    funnels.resize(2);
    // Start this alignment 
    funnels[0].start(aln1.name());
    funnels[1].start(aln2.name());
    
    // Annotate the original read with metadata
    if (!sample_name.empty()) {
        aln1.set_sample_name(sample_name);
        aln2.set_sample_name(sample_name);
    }
    if (!read_group.empty()) {
        aln1.set_read_group(read_group);
        aln2.set_read_group(read_group);
    }
   
    if (track_provenance) {
        // Start the minimizer finding stage
        funnels[0].stage("minimizer");
        funnels[1].stage("minimizer");
    }
    
    // We will find all the seed hits
    vector<vector<pos_t>> seeds;
    
    // This will hold all the minimizers in the query
    vector<vector<gbwtgraph::DefaultMinimizerIndex::minimizer_type>> minimizers;

    // And either way this will map from seed to minimizer that generated it
    vector<vector<size_t>> seed_to_source;
    vector<vector<double>> minimizer_score;
    
    for (size_t read_num = 0 ; read_num < 2 ; read_num++) {
        // Find minimizers in the query
        minimizers.emplace_back(read_num == 0 ? minimizer_index.minimizers(aln1.sequence()) : 
                                                minimizer_index.minimizers(aln2.sequence()) );

        seeds.emplace_back();
        seed_to_source.emplace_back();
        
        if (track_provenance) {
            // Record how many we found, as new lines.
            funnels[read_num].introduce(minimizers.back().size());
            
            // Start the minimizer locating stage
            funnels[read_num].stage("seed");
        }
    
        // Compute minimizer scores for all minimizers as 1 + ln(hard_hit_cap) - ln(hits).
        minimizer_score.emplace_back(minimizers.back().size(), 0.0);
        double base_target_score = 0.0;
        for (size_t i = 0; i < minimizers.back().size(); i++) {
            size_t hits = minimizer_index.count(minimizers.back()[i]);
            if (hits > 0) {
                if (hits <= hard_hit_cap) {
                    minimizer_score.back()[i] = 1.0 + std::log(hard_hit_cap) - std::log(hits);
                } else {
                    minimizer_score.back()[i] = 1.0;
                }
            }
            base_target_score += minimizer_score.back()[i];
        }
        double target_score = (base_target_score * minimizer_score_fraction) + 0.000001;
    
        // Sort the minimizers by score.
        std::vector<size_t> minimizers_in_order(minimizers.back().size());
        for (size_t i = 0; i < minimizers_in_order.size(); i++) {
            minimizers_in_order[i] = i;
        }
        std::sort(minimizers_in_order.begin(), minimizers_in_order.end(), [&minimizer_score](const size_t a, const size_t b) {
            return (minimizer_score.back()[a] > minimizer_score.back()[b]);
        });
    
        // Select the minimizers we use for seeds.
        size_t rejected_count = 0;
        double selected_score = 0.0;
        for (size_t i = 0; i < minimizers.back().size(); i++) {
            size_t minimizer_num = minimizers_in_order[i];
    
            if (track_provenance) {
                // Say we're working on it
                funnels[read_num].processing_input(minimizer_num);
            }
    
            // Select the minimizer if it is informative enough or if the total score
            // of the selected minimizers is not high enough.
            size_t hits = minimizer_index.count(minimizers.back()[minimizer_num]);
            
#ifdef debug
            cerr << "Minimizer " << minimizer_num << " in read " << read_num << " = " << minimizers.back()[minimizer_num].key.decode(minimizer_index.k())
                << " has " << hits << " hits" << endl;
#endif
            
            if (hits == 0) {
                // A minimizer with no hits can't go on.
                if (track_provenance) {
                    funnels[read_num].fail("any-hits", minimizer_num);
                }
            } else if (hits <= hit_cap || (hits <= hard_hit_cap && selected_score + minimizer_score.back()[minimizer_num] <= target_score)) {
                // Locate the hits.
                for (auto& hit : minimizer_index.find(minimizers.back()[minimizer_num])) {
                    // Reverse the hits for a reverse minimizer
                    if (minimizers.back()[minimizer_num].is_reverse) {
                        size_t node_length = gbwt_graph.get_length(gbwt_graph.get_handle(id(hit)));
                        hit = reverse_base_pos(hit, node_length);
                    }
                    // For each position, remember it and what minimizer it came from
                    seeds.back().push_back(hit);
                    seed_to_source.back().push_back(minimizer_num);
                }
                selected_score += minimizer_score.back()[minimizer_num];
                
                if (track_provenance) {
                    // Record in the funnel that this minimizer gave rise to these seeds.
                    funnels[read_num].pass("any-hits", minimizer_num);
                    funnels[read_num].pass("hard-hit-cap", minimizer_num);
                    funnels[read_num].pass("hit-cap||score-fraction", minimizer_num, selected_score  / base_target_score);
                    funnels[read_num].expand(minimizer_num, hits);
                }
            } else if (hits <= hard_hit_cap) {
                // Passed hard hit cap but failed score fraction/normal hit cap
                rejected_count++;
                if (track_provenance) {
                    funnels[read_num].pass("any-hits", minimizer_num);
                    funnels[read_num].pass("hard-hit-cap", minimizer_num);
                    funnels[read_num].fail("hit-cap||score-fraction", minimizer_num, (selected_score + minimizer_score.back()[minimizer_num]) / base_target_score);
                }
            } else {
                // Failed hard hit cap
                rejected_count++;
                if (track_provenance) {
                    funnels[read_num].pass("any-hits", minimizer_num);
                    funnels[read_num].fail("hard-hit-cap", minimizer_num);
                }
            }
            if (track_provenance) {
                // Say we're done with this input item
                funnels[read_num].processed_input();
            }
        }
#ifdef debug
        cerr << "For read " << read_num << ", found " << seeds.back().size() << " seeds from " << (minimizers.back().size() - rejected_count) << " minimizers, rejected " << rejected_count << endl;
#endif
        if (track_provenance && track_correctness) {
            // Tag seeds with correctness based on proximity along paths to the input read's refpos
            funnels[read_num].substage("correct");
          
            if (path_graph == nullptr) {
                cerr << "error[vg::MinimizerMapper] Cannot use track_correctness with no XG index" << endl;
                exit(1);
            }
            
            if (read_num == 0 ? (aln1.refpos_size() != 0) :  (aln2.refpos_size() != 0)) {
                // Take the first refpos as the true position.
                auto& true_pos = aln1.refpos(0);
                
                for (size_t i = 0; i < seeds[read_num].size(); i++) {
                    // Find every seed's reference positions. This maps from path name to pairs of offset and orientation.
                    auto offsets = algorithms::nearest_offsets_in_paths(path_graph, seeds[read_num][i], 100);
                    for (auto& hit_pos : offsets[path_graph->get_path_handle(true_pos.name())]) {
                        // Look at all the ones on the path the read's true position is on.
                        if (abs((int64_t)hit_pos.first - (int64_t) true_pos.offset()) < 200) {
                            // Call this seed hit close enough to be correct
                            funnels[read_num].tag_correct(i);
                        }
                    }
                }
            }
        }
    }

    if (track_provenance) {
        // Begin the clustering stage
        funnels[0].stage("cluster");
        funnels[1].stage("cluster");
    }
        
    // Cluster the seeds. Get sets of input seed indexes that go together.
    vector<vector<pair<vector<size_t>, size_t>>> read_clusters = clusterer.cluster_seeds(seeds, distance_limit, expected_fragment_length);
    
    if (track_provenance) {
        funnels[0].substage("score");
        funnels[1].substage("score");
    }

    //For each fragment cluster (cluster of clusters), for each read, a vector of all alignments
    vector<vector<vector<Alignment>>> alignments;

    for (size_t read_num = 0 ; read_num < 2 ; read_num++) {
        Alignment& aln = read_num == 0 ? aln1 : aln2;
        vector<pair<vector<size_t>, size_t>>& clusters = read_clusters[read_num];

        // Cluster score is the sum of minimizer scores.
        vector<double> cluster_score(clusters.size(), 0.0);
        vector<double> read_coverage_by_cluster;
        read_coverage_by_cluster.reserve(clusters.size());

        for (size_t i = 0; i < clusters.size(); i++) {
            // For each cluster
            auto& cluster = clusters[i].first;
            
            if (track_provenance) {
                // Say we're making it
                funnels[read_num].producing_output(i);
            }

            // Which minimizers are present in the cluster.
            vector<bool> present(minimizers[read_num].size(), false);
            for (auto hit_index : cluster) {
                present[seed_to_source.back()[hit_index]] = true;
            }

            // Compute the score.
            for (size_t j = 0; j < minimizers[read_num].size(); j++) {
                if (present[j]) {
                    cluster_score[i] += minimizer_score[read_num][j];
                }
            }
            
            if (track_provenance) {
                // Record the cluster in the funnel as a group of the size of the number of items.
                funnels[read_num].merge_group(cluster.begin(), cluster.end());
                funnels[read_num].score(funnels[read_num].latest(), cluster_score[i]);
                
                // Say we made it.
                funnels[read_num].produced_output();
            }

            //Get the cluster coverage
            // We set bits in here to true when query anchors cover them
            sdsl::bit_vector covered(read_num == 0 ? aln1.sequence().size() : aln2.sequence().size(), 0);
            std::uint64_t k_bit_mask = sdsl::bits::lo_set[minimizer_index.k()];

            for (auto hit_index : cluster) {
                // For each hit in the cluster, work out what anchor sequence it is from.
                size_t source_index = seed_to_source[read_num][hit_index];

                // The offset of a reverse minimizer is the endpoint of the kmer
                size_t start_offset = minimizers[read_num][source_index].offset;
                if (minimizers[read_num][source_index].is_reverse) {
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
        cerr << "Found " << clusters.size() << " clusters for read " << read_num << endl;
#endif
                                        
        // Retain clusters only if their score is better than this, in addition to the coverage cutoff
        double cluster_score_cutoff = cluster_score.size() == 0 ? 0 :
                        *std::max_element(cluster_score.begin(), cluster_score.end()) - cluster_score_threshold;
        
        if (track_provenance) {
            // Now we go from clusters to gapless extensions
            funnels[read_num].stage("extend");
        }
        
        // These are the GaplessExtensions for all the clusters (and fragment cluster assignments), in cluster_indexes_in_order order.
        vector<pair<vector<GaplessExtension>, size_t>> cluster_extensions;
        cluster_extensions.reserve(clusters.size());
        size_t max_fragment_num = 0;
        //For each cluster, what fraction of "equivalent" clusters did we keep?
        //TODO: Maybe put this back vector<vector<double>> probability_cluster_lost;
        //What is the score and coverage we are considering and how many reads
        //size_t curr_coverage = 0;
        //size_t curr_score = 0;
        //size_t curr_kept = 0;
        //size_t curr_count = 0;
        
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
                    funnels[read_num].pass("cluster-coverage", cluster_num, read_coverage_by_cluster[cluster_num]);
                    funnels[read_num].pass("max-extensions", cluster_num);
                }
                
                // First check against the additional score filter
                if (cluster_score_threshold != 0 && cluster_score[cluster_num] < cluster_score_cutoff) {
                    //If the score isn't good enough, ignore this cluster
                    if (track_provenance) {
                        funnels[read_num].fail("cluster-score", cluster_num, cluster_score[cluster_num]);
                    }
                    return false;
                }
                
                if (track_provenance) {
                    funnels[read_num].pass("cluster-score", cluster_num, cluster_score[cluster_num]);
                    funnels[read_num].processing_input(cluster_num);
                }

//TODO
//                if (read_coverage_by_cluster[cluster_num] == curr_coverage &&
//                    cluster_score[cluster_num] == curr_score &&
//                    curr_kept < max_extensions * 0.75) {
//                    curr_kept++;
//                    curr_count++;
//                } else if (!read_coverage_by_cluster[cluster_num] == curr_coverage ||
//                        !cluster_score[cluster_num] == curr_score) {
//                        //If this is a cluster that has scores different than the previous one
//                        for (size_t i = 0 ; i < curr_kept ; i++ ) {
//                            probability_cluster_lost[read_num].push_back(1.0 - (double(curr_kept) / double(curr_count)));
//                        }
//                        curr_coverage = read_coverage_by_cluster[cluster_num];
//                        curr_score = cluster_score[cluster_num];
//                        curr_kept = 1;
//                        curr_count = 1;
//                } else {
//                    //If this cluster is equivalent to the previous one and we already took enough
//                    //equivalent clusters
//                    curr_count ++;
//                    return false;
//                }
//                
//
//                //Only keep this cluster if we have few enough equivalent clusters
//
                vector<size_t>& cluster = clusters[cluster_num].first;

#ifdef debug
                cerr << "Cluster " << cluster_num << endl;
#endif
                 
                // Pack the seeds for GaplessExtender.
                GaplessExtender::cluster_type seed_matchings;
                for (auto& seed_index : cluster) {
                    // Insert the (graph position, read offset) pair.
                    seed_matchings.insert(GaplessExtender::to_seed(seeds[read_num][seed_index], minimizers[read_num][seed_to_source.back()[seed_index]].offset));
#ifdef debug
                    cerr << "Seed read:" << minimizers[read_num][seed_to_source.back()[seed_index]].offset << " = " << seeds[seed_index]
                        << " from minimizer " << seed_to_source.back()[seed_index] << "(" << minimizer_index.count(minimizers[read_num][seed_to_source.back()[seed_index]]) << ")" << endl;
#endif
                }
                
                // Extend seed hits in the cluster into one or more gapless extensions
                cluster_extensions.emplace_back(std::move(extender.extend(seed_matchings, aln.sequence())), 
                                                clusters[cluster_num].second);
                max_fragment_num = max(max_fragment_num, clusters[cluster_num].second);
                
                if (track_provenance) {
                    // Record with the funnel that the previous group became a group of this size.
                    // Don't bother recording the seed to extension matching...
                    funnels[read_num].project_group(cluster_num, cluster_extensions.back().first.size());
                    
                    // Say we finished with this cluster, for now.
                    funnels[read_num].processed_input();
                }
                return true;
                
            }, [&](size_t cluster_num) {
                // There are too many sufficiently good clusters
                if (track_provenance) {
                    funnels[read_num].pass("cluster-coverage", cluster_num, read_coverage_by_cluster[cluster_num]);
                    funnels[read_num].fail("max-extensions", cluster_num);
                }
//                if (read_coverage_by_cluster[cluster_num] == curr_coverage &&
//                    cluster_score[cluster_num] == curr_score) {
//                    curr_count ++;
//                } else {
//        
//                    //TODO:
//                    //for (size_t i = 0 ; i < curr_kept ; i++ ) {
//                    //    probability_cluster_lost[read_num].push_back(1.0 - (double(curr_kept) / double(curr_count)));
//                    //}
//                    curr_score = 0;
//                    curr_coverage = 0;
//                    curr_kept = 0;
//                    curr_count = 0;
//                }
            }, [&](size_t cluster_num) {
                // This cluster is not sufficiently good.
                if (track_provenance) {
                    funnels[read_num].fail("cluster-coverage", cluster_num, read_coverage_by_cluster[cluster_num]);
                }
//                //TODO:
//                for (size_t i = 0 ; i < curr_kept ; i++ ) {
//                    probability_cluster_lost[read_num].push_back(1.0 - (double(curr_kept) / double(curr_count)));
//                }
//                curr_kept = 0;
//                curr_count = 0;
//                curr_score = 0;
//                curr_coverage = 0;
            });
            
            //TODO
            //for (size_t i = 0 ; i < curr_kept ; i++ ) {
            //    probability_cluster_lost[read_num].push_back(1.0 - (double(curr_kept) / double(curr_count)));
            //}
        
        if (track_provenance) {
            funnels[read_num].substage("score");
        }

        // We now estimate the best possible alignment score for each cluster.
        vector<int> cluster_extension_scores;
        cluster_extension_scores.reserve(cluster_extensions.size());
        for (size_t i = 0; i < cluster_extensions.size(); i++) {
            // For each group of GaplessExtensions
            
            if (track_provenance) {
                funnels[read_num].producing_output(i);
            }
            
            vector<GaplessExtension>& extensions = cluster_extensions[i].first;
            // Work out the score of the best chain of extensions, accounting for gaps/overlaps in the read only
            cluster_extension_scores.push_back(score_extension_group(aln, extensions,
                get_regular_aligner()->gap_open, get_regular_aligner()->gap_extension));
            
            if (track_provenance) {
                // Record the score with the funnel
                funnels[read_num].score(i, cluster_extension_scores.back());
                funnels[read_num].produced_output();
            }
        }
        
        if (track_provenance) {
            funnels[read_num].stage("align");
        }
        
        // Now start the alignment step. Everything has to become an alignment.

        // We will fill this with all computed alignments in estimated score order.
        alignments.resize(max_fragment_num);

        
        // Clear any old refpos annotation and path
        aln.clear_refpos();
        aln.clear_path();
        aln.set_score(0);
        aln.set_identity(0);
        aln.set_mapping_quality(0);
        
        // Go through the gapless extension groups in score order.
        process_until_threshold(cluster_extensions, cluster_extension_scores,
            extension_set_score_threshold, 2, max_alignments,
            [&](size_t extension_num) {
                // This extension set is good enough.
                // Called in descending score order.
                
                if (track_provenance) {
                    funnels[read_num].pass("extension-set", extension_num, cluster_extension_scores[extension_num]);
                    funnels[read_num].pass("max-alignments", extension_num);
                    funnels[read_num].processing_input(extension_num);
                }
                
                auto& extensions = cluster_extensions[extension_num].first;
                
                // Get an Alignments best_ and second_best_alignment of it somehow, and throw it in.
                Alignment best_alignment = aln;
                Alignment second_best_alignment = aln;
                
                if (extensions[0].full()) {
                    // We got full-length extensions, so directly convert to an Alignment.
                    
                    if (track_provenance) {
                        funnels[read_num].substage("direct");
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
                        funnels[read_num].substage_stop();
                    }
                } else if (do_dp) {
                    // We need to do chaining.
                    
                    if (track_provenance) {
                        funnels[read_num].substage("chain");
                    }
                    
                    // Do the DP and compute alignment into best_alignment and
                    // second_best_alignment, if there is a second best
                    find_optimal_tail_alignments(aln, extensions, best_alignment, second_best_alignment);

                    
                    if (track_provenance) {
                        // We're done chaining. Next alignment may not go through this substage.
                        funnels[read_num].substage_stop();
                    }
                } else {
                    // We would do chaining but it is disabled.
                    // Leave best_alignment unaligned
                }
                
                
                
                if (second_best_alignment.score() != 0 && 
                    second_best_alignment.score() > best_alignment.score() * 0.8) {
                    //If there is a second extension and its score is at least half of the best score
                    alignments[cluster_extensions[extension_num].second ][read_num].push_back(std::move(second_best_alignment));

                    if (track_provenance) {
        
                        funnels[read_num].project(extension_num);
                        funnels[read_num].score(alignments[cluster_extensions[extension_num].second ][read_num].size() - 1, 
                                     alignments[cluster_extensions[extension_num].second ][read_num].back().score());
                        // We're done with this input item
                        funnels[read_num].processed_input();
                    }
                }

                alignments[cluster_extensions[extension_num].second ][read_num].push_back(std::move(best_alignment));

                if (track_provenance) {

                    funnels[read_num].project(extension_num);
                    funnels[read_num].score(alignments[cluster_extensions[extension_num].second ][read_num].size() - 1, 
                    alignments[cluster_extensions[extension_num].second ][read_num].back().score());
                    
                    // We're done with this input item
                    funnels[read_num].processed_input();
                }
                
                return true;
            }, [&](size_t extension_num) {
                // There are too many sufficiently good extensions
                if (track_provenance) {
                    funnels[read_num].pass("extension-set", extension_num, cluster_extension_scores[extension_num]);
                    funnels[read_num].fail("max-alignments", extension_num);
                }
            }, [&](size_t extension_num) {
                // This extension is not good enough.
                if (track_provenance) {
                    funnels[read_num].fail("extension-set", extension_num, cluster_extension_scores[extension_num]);
                }
            });
        
    }
    
    if (track_provenance) {
        // Now say we are finding the winner(s)
        funnels[0].stage("winner");
        funnels[1].stage("winner");
    }
    
    // Fill this in with the alignments we will output
    // Tuple of fragment index, index of first alignment, index of second alignment in alignments
    vector<tuple<size_t, size_t, size_t>> paired_alignments;
    // Grab all the scores in order for MAPQ computation.
    vector<double> paired_scores;
    paired_alignments.reserve(alignments.size());
    paired_scores.reserve(alignments.size());


    for (size_t fragment_num = 0 ; fragment_num < alignments.size() ; fragment_num ++ ) {
        //Get pairs of plausible alignments
        //TODO: Maybe don't keep all pairs per fragment cluster
        vector<vector<Alignment>>& fragment_alignments = alignments[fragment_num];
        if (!fragment_alignments[0].empty() && ! fragment_alignments[1].empty()) {
            for (size_t i1 = 0 ; i1 < fragment_alignments[0].size() ; i1++)  {
                Alignment& alignment1 = fragment_alignments[0][i1];
                for (size_t i2 = 0 ; i2 < fragment_alignments[1].size() ; i2++) {
                    Alignment& alignment2 = fragment_alignments[1][i2];
                    pos_t pos1 = initial_position(alignment1.path());
                    pos_t pos2 = final_position(alignment2.path());
                    int64_t fragment_distance = distance_index.minDistance(pos1, pos2); 
                    //TODO: Scoring of pairs
                    double score = alignment1.score() + alignment2.score() - (double)abs(expected_fragment_length-fragment_distance); 
                    paired_alignments.emplace_back(fragment_num, i1, i2);
                    paired_scores.emplace_back(score);
                }
            }
        }
    }
    
    // Fill this in with the alignments we will output
    // Same as paired_alignments: indexes into alignments
    vector<tuple<size_t, size_t, size_t>> mappings;
    // Grab all the scores in order for MAPQ computation.
    vector<double> scores;
    mappings.reserve(paired_alignments.size());
    scores.reserve(paired_alignments.size());
    process_until_threshold(paired_alignments, (std::function<double(size_t)>) [&](size_t i) -> double {
        return paired_scores[i];
    }, 0, 1, max_multimaps, [&](size_t alignment_num) {
        // This alignment makes it
        // Called in score order
        
        // Remember the score at its rank
        scores.emplace_back(paired_scores[alignment_num]);
        
        // Remember the output alignment
        mappings.emplace_back(paired_alignments[alignment_num]);
        
//TODO: After this point, the indices of each item in the funnel doesn't make much sense since we've re-ordered everything based on the fragment cluster it belongs to        
//        if (track_provenance) {
//            // Tell the funnel
//            funnels[0].pass("max-multimaps", std::get<1>(paired_alignments[alignment_num]));
//            funnels[0].project(alignment_num);
//            funnels[0].score(alignment_num, scores.back());
//        }
        
        return true;
    }, [&](size_t alignment_num) {
        // We already have enough alignments, although this one has a good score
        
        // Remember the score at its rank anyway
        scores.emplace_back(paired_scores[alignment_num]);
        
//        if (track_provenance) {
//            funnel.fail("max-multimaps", alignment_num);
//        }
    }, [&](size_t alignment_num) {
        // This alignment does not have a sufficiently good score
        // Score threshold is 0; this should never happen
        assert(false);
    });
    
//    if (track_provenance) {
//        funnel.substage("mapq");
//    }
    
#ifdef debug
    cerr << "For scores ";
    for (auto& score : scores) cerr << score << " ";
#endif

    size_t winning_index;
    // Compute MAPQ if not unmapped. Otherwise use 0 instead of the 50% this would give us.
    double mapq = (mappings.empty() || scores[0] == 0) ? 0 : 
        get_regular_aligner()->maximum_mapping_quality_exact(scores, &winning_index) / 2;
    
#ifdef debug
    cerr << "MAPQ is " << mapq << endl;
#endif
 
    pair<vector<Alignment>, vector<Alignment>> paired_mappings;
    if (mappings.empty()) {
        //If we didn't get an alignment, return empty alignments
        paired_mappings.first.emplace_back(aln1);
        paired_mappings.second.emplace_back(aln2);
        paired_mappings.first.back().clear_refpos();
        paired_mappings.first.back().clear_path();
        paired_mappings.first.back().set_score(0);
        paired_mappings.first.back().set_identity(0);
        paired_mappings.first.back().set_mapping_quality(0);
        paired_mappings.second.back().clear_refpos();
        paired_mappings.second.back().clear_path();
        paired_mappings.second.back().set_score(0);
        paired_mappings.second.back().set_identity(0);
        paired_mappings.second.back().set_mapping_quality(0);

    } else {
        for (size_t i = 0 ; i < mappings.size() ; i++) {
            paired_mappings.first.emplace_back( alignments[std::get<0>(mappings[i])][0][std::get<1>(mappings[i])]);
            paired_mappings.second.emplace_back(alignments[std::get<0>(mappings[i])][1][std::get<2>(mappings[i])]);
            // Make sure to clamp 0-60.
            // TODO: Maybe don't just give them the same mapq
            if (i == 0 ) {
                paired_mappings.first.back().set_mapping_quality(max(min(mapq, 60.0), 0.0));
                paired_mappings.second.back().set_mapping_quality(max(min(mapq, 60.0), 0.0));
            } else {
                paired_mappings.first.back().set_is_secondary(true);
                paired_mappings.second.back().set_is_secondary(true);
            }
        }
    }
   
    
    if (track_provenance) {
        funnels[0].substage_stop();
        funnels[1].substage_stop();
    }
    
    // Stop this alignment
    funnels[0].stop();
    funnels[1].stop();
    
    if (track_provenance) {
    
        // Annotate with the number of results in play at each stage
        for (size_t read_num = 0 ; read_num < 2 ; read_num++) {
            funnels[read_num].for_each_stage([&](const string& stage, const vector<size_t>& result_sizes) {
                // Save the number of items
                set_annotation(paired_mappings.first [0], "stage_" + stage + "_results", (double)result_sizes.size());
                set_annotation(paired_mappings.second[0], "stage_" + stage + "_results", (double)result_sizes.size());
            });
            // Annotate with the performances of all the filters
            // We need to track filter number
            size_t filter_num = 0;
            funnels[read_num].for_each_filter([&](const string& stage, const string& filter,
                const Funnel::FilterPerformance& by_count, const Funnel::FilterPerformance& by_size,
                const vector<double>& filter_statistics_correct, const vector<double>& filter_statistics_non_correct) {
                
                string filter_id = to_string(filter_num) + "_" + filter + "_" + stage;
                
                if (read_num == 0) {
                    // Save the stats
                    set_annotation(paired_mappings.first[0], "filter_" + filter_id + "_passed_count_total", (double) by_count.passing);
                    set_annotation(paired_mappings.first[0], "filter_" + filter_id + "_failed_count_total", (double) by_count.failing);
                    set_annotation(paired_mappings.first[0], "filter_" + filter_id + "_passed_size_total", (double) by_size.passing);
                    set_annotation(paired_mappings.first[0], "filter_" + filter_id + "_failed_size_total", (double) by_size.failing);
                    
                    if (track_correctness) {
                        set_annotation(paired_mappings.first[0], "filter_" + filter_id + "_passed_count_correct", (double) by_count.passing_correct);
                        set_annotation(paired_mappings.first[0], "filter_" + filter_id + "_failed_count_correct", (double) by_count.failing_correct);
                        set_annotation(paired_mappings.first[0], "filter_" + filter_id + "_passed_size_correct", (double) by_size.passing_correct);
                        set_annotation(paired_mappings.first[0], "filter_" + filter_id + "_failed_size_correct", (double) by_size.failing_correct);
                    }
                    
                    // Save the correct and non-correct filter statistics, even if
                    // everything is non-correct because correctness isn't computed
                    set_annotation(paired_mappings.first[0], "filterstats_" + filter_id + "_correct", filter_statistics_correct);
                    set_annotation(paired_mappings.first[0], "filterstats_" + filter_id + "_noncorrect", filter_statistics_non_correct);
                
                } else {
                    // Save the stats
                    set_annotation(paired_mappings.second[0], "filter_" + filter_id + "_passed_count_total", (double) by_count.passing);
                    set_annotation(paired_mappings.second[0], "filter_" + filter_id + "_failed_count_total", (double) by_count.failing);
                    set_annotation(paired_mappings.second[0], "filter_" + filter_id + "_passed_size_total", (double) by_size.passing);
                    set_annotation(paired_mappings.second[0], "filter_" + filter_id + "_failed_size_total", (double) by_size.failing);
                    
                    if (track_correctness) {
                        set_annotation(paired_mappings.second[0], "filter_" + filter_id + "_passed_count_correct", (double) by_count.passing_correct);
                        set_annotation(paired_mappings.second[0], "filter_" + filter_id + "_failed_count_correct", (double) by_count.failing_correct);
                        set_annotation(paired_mappings.second[0], "filter_" + filter_id + "_passed_size_correct", (double) by_size.passing_correct);
                        set_annotation(paired_mappings.second[0], "filter_" + filter_id + "_failed_size_correct", (double) by_size.failing_correct);
                    }
                    
                    // Save the correct and non-correct filter statistics, even if
                    // everything is non-correct because correctness isn't computed
                    set_annotation(paired_mappings.second[0], "filterstats_" + filter_id + "_correct", filter_statistics_correct);
                    set_annotation(paired_mappings.second[0], "filterstats_" + filter_id + "_noncorrect", filter_statistics_non_correct);
                }
                filter_num++;
            });
        }
//        
//        if (track_correctness) {
//            // And with the last stage at which we had any descendants of the correct seed hit locations
//            set_annotation(mappings[0], "last_correct_stage", funnel.last_correct_stage());
//        }
//        
//        
        // Annotate with parameters used for the filters.
        set_annotation(paired_mappings.first [0], "param_hit-cap", (double) hit_cap);
        set_annotation(paired_mappings.first [0], "param_hard-hit-cap", (double) hard_hit_cap);
        set_annotation(paired_mappings.first [0], "param_score-fraction", (double) minimizer_score_fraction);
        set_annotation(paired_mappings.first [0], "param_max-extensions", (double) max_extensions);
        set_annotation(paired_mappings.first [0], "param_max-alignments", (double) max_alignments);
        set_annotation(paired_mappings.first [0], "param_cluster-score", (double) cluster_score_threshold);
        set_annotation(paired_mappings.first [0], "param_cluster-coverage", (double) cluster_coverage_threshold);
        set_annotation(paired_mappings.first [0], "param_extension-set", (double) extension_set_score_threshold);
        set_annotation(paired_mappings.first [0], "param_max-multimaps", (double) max_multimaps);
        set_annotation(paired_mappings.second[0], "param_hit-cap", (double) hit_cap);
        set_annotation(paired_mappings.second[0], "param_hard-hit-cap", (double) hard_hit_cap);
        set_annotation(paired_mappings.second[0], "param_score-fraction", (double) minimizer_score_fraction);
        set_annotation(paired_mappings.second[0], "param_max-extensions", (double) max_extensions);
        set_annotation(paired_mappings.second[0], "param_max-alignments", (double) max_alignments);
        set_annotation(paired_mappings.second[0], "param_cluster-score", (double) cluster_score_threshold);
        set_annotation(paired_mappings.second[0], "param_cluster-coverage", (double) cluster_coverage_threshold);
        set_annotation(paired_mappings.second[0], "param_extension-set", (double) extension_set_score_threshold);
        set_annotation(paired_mappings.second[0], "param_max-multimaps", (double) max_multimaps);
    }
    
    // Ship out all the aligned alignments
    alignment_emitter.emit_mapped_pair(std::move(paired_mappings.first), std::move(paired_mappings.second));

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


