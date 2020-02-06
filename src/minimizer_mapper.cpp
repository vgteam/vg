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

//#define debug

namespace vg {

using namespace std;

MinimizerMapper::MinimizerMapper(const gbwtgraph::GBWTGraph& graph,
    const std::vector<std::unique_ptr<gbwtgraph::DefaultMinimizerIndex>>& minimizer_indexes,
    MinimumDistanceIndex& distance_index, const PathPositionHandleGraph* path_graph) :
    path_graph(path_graph), minimizer_indexes(minimizer_indexes),
    distance_index(distance_index), gbwt_graph(graph),
    extender(gbwt_graph, *(get_regular_aligner())), clusterer(distance_index),
    fragment_length_distr(1000,1000,0.95){
        //TODO: Picked fragment_length_distr params from mpmap
    
    // Nothing to do!
}

void MinimizerMapper::map(Alignment& aln, AlignmentEmitter& alignment_emitter) {
    // Ship out all the aligned alignments
    alignment_emitter.emit_mapped_single(map(aln));
}

vector<Alignment> MinimizerMapper::map(Alignment& aln) {
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
    
    // Minimizers from each index and scores as 1 + ln(hard_hit_cap) - ln(hits).
    struct Minimizer {
        typename gbwtgraph::DefaultMinimizerIndex::minimizer_type value;
        size_t hits;
        const typename gbwtgraph::DefaultMinimizerIndex::code_type* occs;
        size_t origin; // From minimizer_indexes[origin].
        double score;

        // Sort the minimizers in descending order by score.
        bool operator< (const Minimizer& another) const {
            return (this->score > another.score);
        }
    };
    std::vector<Minimizer> minimizers;

    // Find the minimizers and score them.
    double base_target_score = 0.0;
    double base_score = 1.0 + std::log(hard_hit_cap);
    for (size_t i = 0; i < minimizer_indexes.size(); i++) {
        auto current_minimizers = minimizer_indexes[i]->minimizers(aln.sequence());
        for (auto& minimizer : current_minimizers) {
            double score = 0.0;
            auto hits = minimizer_indexes[i]->count_and_find(minimizer);
            if (hits.first > 0) {
                if (hits.first <= hard_hit_cap) {
                    score = base_score - std::log(hits.first);
                } else {
                    score = 1.0;
                }
            }
            minimizers.push_back({ minimizer, hits.first, hits.second, i, score });
            base_target_score += score;
        }
    }
    double target_score = (base_target_score * minimizer_score_fraction) + 0.000001;
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
        
        if (minimizer.hits == 0) {
            // A minimizer with no hits can't go on.
            if (track_provenance) {
                funnel.fail("any-hits", i);
            }
        } else if (minimizer.hits <= hit_cap || (minimizer.hits <= hard_hit_cap && selected_score + minimizer.score <= target_score)) {
            // Locate the hits.
            for (size_t j = 0; j < minimizer.hits; j++) {
                pos_t hit = gbwtgraph::DefaultMinimizerIndex::decode(minimizer.occs[j]);
                // Reverse the hits for a reverse minimizer
                if (minimizer.value.is_reverse) {
                    size_t node_length = gbwt_graph.get_length(gbwt_graph.get_handle(id(hit)));
                    hit = reverse_base_pos(hit, node_length);
                }
                // For each position, remember it and what minimizer it came from
                seeds.push_back(hit);
                seed_to_source.push_back(i);
            }
            selected_score += minimizer.score;
            
            if (track_provenance) {
                // Record in the funnel that this minimizer gave rise to these seeds.
                funnel.pass("any-hits", i);
                funnel.pass("hard-hit-cap", i);
                funnel.pass("hit-cap||score-fraction", i, selected_score  / base_target_score);
                funnel.expand(i, minimizer.hits);
            }
        } else if (minimizer.hits <= hard_hit_cap) {
            // Passed hard hit cap but failed score fraction/normal hit cap
            rejected_count++;
            if (track_provenance) {
                funnel.pass("any-hits", i);
                funnel.pass("hard-hit-cap", i);
                funnel.fail("hit-cap||score-fraction", i, (selected_score + minimizer.score) / base_target_score);
            }
        } else {
            // Failed hard hit cap
            rejected_count++;
            if (track_provenance) {
                funnel.pass("any-hits", i);
                funnel.fail("hard-hit-cap", i);
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
                seed_matchings.insert(GaplessExtender::to_seed(seeds[seed_index], minimizers[seed_to_source[seed_index]].value.offset));
#ifdef debug
                const Minimizer& minimizer = minimizers[seed_to_source[seed_index]];
                cerr << "Seed read:" << minimizer.value.offset << " = " << seeds[seed_index]
                    << " from minimizer " << seed_to_source[seed_index] << "(" << minimizer.hits << ")" << endl;
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
    
#ifdef debug
    // Dump the funnel info graph.
    funnel.to_dot(cerr);
#endif

    return mappings;
}
pair<vector<Alignment>, vector<Alignment>> MinimizerMapper::map_paired(Alignment& aln1, Alignment& aln2,
                                                      vector<pair<Alignment, Alignment>>& ambiguous_pair_buffer){

    if (fragment_length_distr.is_finalized()) {
        //If we know the fragment length distribution then we just map paired ended 
        return map_paired(aln1, aln2);

    } else {
        //If we don't know the fragment length distribution, map the reads single ended

        vector<Alignment> alns1(map(aln1));
        vector<Alignment> alns2(map(aln2));

        // TODO: How do we decide what an unambiguous mapping is? 
        // TODO: Double check what the actual best possible score is - which aligner are we actually using? 
        // Check if the separately-mapped ends are both sufficiently perfect and sufficiently unique
        int32_t max_score_aln_1 = get_regular_aligner()->score_exact_match(aln1, 0, aln1.sequence().size());
        int32_t max_score_aln_2 = get_regular_aligner()->score_exact_match(aln2, 0, aln2.sequence().size());
        if (!alns1.empty() && ! alns2.empty()  && 
            alns1.front().mapping_quality() == 60 && alns2.front().mapping_quality() == 60 &&
            alns1.front().score() >= max_score_aln_1 * 0.85 && alns2.front().score() >= max_score_aln_2 * 0.85) {

            //Flip the second alignment to get the proper fragment distance 
            reverse_complement_alignment_in_place(&alns2.front(), [&](vg::id_t node_id) {
                    return gbwt_graph.get_length(gbwt_graph.get_handle(node_id));
                    });           
            int64_t dist = distance_between(alns1.front(), alns2.front());
            // And that they have an actual pair distance and set of relative orientations

            if (dist == std::numeric_limits<int64_t>::max()) {
                //If the distance between them is ambiguous, leave them unmapped

                ambiguous_pair_buffer.emplace_back(aln1, aln2);
                pair<vector<Alignment>, vector<Alignment>> empty;
                return empty;
            }

            //If we're keeping this alignment, flip the second alignment back
            reverse_complement_alignment_in_place(&alns2.front(), [&](vg::id_t node_id) {
                    return gbwt_graph.get_length(gbwt_graph.get_handle(node_id));
                    });           
            // If that all checks out, say they're mapped, emit them, and register their distance and orientations
            fragment_length_distr.register_fragment_length(dist);

            pair<vector<Alignment>, vector<Alignment>> mapped_pair;
            mapped_pair.first.emplace_back(alns1.front());
            mapped_pair.second.emplace_back(alns2.front());
            return mapped_pair;

        } else {
            // Otherwise, discard the mappings and put them in the ambiguous buffer

            ambiguous_pair_buffer.emplace_back(aln1, aln2);
            pair<vector<Alignment>, vector<Alignment>> empty;
            return empty;
        }
    }
}
pair<vector<Alignment>, vector< Alignment>> MinimizerMapper::map_paired(Alignment& aln1, Alignment& aln2) {
    // For each input alignment
    
#ifdef debug
    cerr << "Read pair " << aln1.name() << ": " << aln1.sequence() << " and " << aln2.name() << ": " << aln2.sequence() << endl;
#endif

    // Assume reads are in inward orientations on input, and
    // convert to rightward orientations before mapping
    // and flip the second read back before output 

    aln2.clear_path();
    reverse_complement_alignment_in_place(&aln2, [&](vg::id_t node_id) {
        return gbwt_graph.get_length(gbwt_graph.get_handle(node_id));
    });


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
    vector<vector<pos_t>> seeds_by_read;
    
    // This will hold all the minimizers in the query, one vector of minimizers per read
    // Minimizers from each index and scores as 1 + ln(hard_hit_cap) - ln(hits).
    struct Minimizer {
        typename gbwtgraph::DefaultMinimizerIndex::minimizer_type value;
        size_t origin; // From minimizer_indexes[origin].
        double score;

        // Sort the minimizers in descending order by score.
        bool operator< (const Minimizer& another) const {
            return (this->score > another.score);
        }
    };
    pair<std::vector<Minimizer>, std::vector<Minimizer>> minimizers_by_read;


    // And either way this will map from seed to minimizer that generated it
    pair<vector<size_t>, vector<size_t>> seed_to_source_by_read;
    
    for (size_t read_num = 0 ; read_num < 2 ; read_num++) {
        std::vector<Minimizer>& minimizers = read_num == 0 ? minimizers_by_read.first : minimizers_by_read.second;
        seeds_by_read.emplace_back();
        vector<pos_t>& seeds = seeds_by_read.back();
        Funnel& funnel = funnels[read_num];
        Alignment& aln = read_num == 0 ? aln1 : aln2;

        vector<size_t>& seed_to_source = read_num == 0 ? seed_to_source_by_read.first : seed_to_source_by_read.second;
        
        if (track_provenance) {
            // Record how many we found, as new lines.
            funnel.introduce(minimizers.size());
            
            // Start the minimizer locating stage
            funnel.stage("seed");
        }

        // Find the minimizers and score them.
        double base_target_score = 0.0;
        double base_score = 1.0 + std::log(hard_hit_cap);
        for (size_t i = 0; i < minimizer_indexes.size(); i++) {
            auto current_minimizers = minimizer_indexes[i]->minimizers(aln.sequence());
            for (auto& minimizer : current_minimizers) {
                double score = 0.0;
                size_t hits = minimizer_indexes[i]->count(minimizer);
                if (hits > 0) {
                    if (hits <= hard_hit_cap) {
                        score = base_score - std::log(hits);
                    } else {
                        score = 1.0;
                    }
                }
                minimizers.push_back({ minimizer, i, score });
                base_target_score += score;
            }
        }
        double target_score = (base_target_score * minimizer_score_fraction) + 0.000001;
        std::sort(minimizers.begin(), minimizers.end());
    
        if (track_provenance) {
            // Record how many we found, as new lines.
            funnel.introduce(minimizers.size());
            
            // Start the minimizer locating stage
            funnel.stage("seed");
        }
    
    
        // Select the minimizers we use for seeds.
        size_t rejected_count = 0;
        double selected_score = 0.0;
        for (size_t i = 0; i < minimizers.size(); i++) {
            if (track_provenance) {
                // Say we're working on it
                funnel.processing_input(i);
            }
    
            // Select the minimizer if it is informative enough or if the total score
            // of the selected minimizers is not high enough.
            const Minimizer& minimizer = minimizers[i];
            size_t hits = minimizer_indexes[minimizer.origin]->count(minimizer.value);
            
#ifdef debug
            cerr << "Minimizer " << i << " = " << minimizer.value.key.decode(minimizer_indexes[minimizer.origin]->k())
                 << " has " << hits << " hits" << endl;
#endif
            
            if (hits == 0) {
                // A minimizer with no hits can't go on.
                if (track_provenance) {
                    funnel.fail("any-hits", i);
                }
            } else if (hits <= hit_cap || (hits <= hard_hit_cap && selected_score + minimizer.score <= target_score)) {
                // Locate the hits.
                for (auto& hit : minimizer_indexes[minimizer.origin]->find(minimizer.value)) {
                    // Reverse the hits for a reverse minimizer
                    if (minimizer.value.is_reverse) {
                        size_t node_length = gbwt_graph.get_length(gbwt_graph.get_handle(id(hit)));
                        hit = reverse_base_pos(hit, node_length);
                    }
                    // For each position, remember it and what minimizer it came from
                    seeds.push_back(hit);
                    seed_to_source.push_back(i);
                }
                selected_score += minimizer.score;
                
                if (track_provenance) {
                    // Record in the funnel that this minimizer gave rise to these seeds.
                    funnel.pass("any-hits", i);
                    funnel.pass("hard-hit-cap", i);
                    funnel.pass("hit-cap||score-fraction", i, selected_score  / base_target_score);
                    funnel.expand(i, hits);
                }
            } else if (hits <= hard_hit_cap) {
                // Passed hard hit cap but failed score fraction/normal hit cap
                rejected_count++;
                if (track_provenance) {
                    funnel.pass("any-hits", i);
                    funnel.pass("hard-hit-cap", i);
                    funnel.fail("hit-cap||score-fraction", i, (selected_score + minimizer.score) / base_target_score);
                }
            } else {
                // Failed hard hit cap
                rejected_count++;
                if (track_provenance) {
                    funnel.pass("any-hits", i);
                    funnel.fail("hard-hit-cap", i);
                }
            }
            if (track_provenance) {
                // Say we're done with this input item
                funnel.processed_input();
            }
        }
#ifdef debug
        cerr << "For read " << read_num <<  " found " << seeds.size() << " seeds from " << (minimizers.size() - rejected_count) << " minimizers, rejected " << rejected_count << endl;
#endif
    
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
    }

    if (track_provenance) {
        // Begin the clustering stage
        funnels[0].stage("cluster");
        funnels[1].stage("cluster");
    }

    // Cluster the seeds. Get sets of input seed indexes that go together.
    // If the fragment length distribution hasn't been fixed yet (if the expected fragment length = 0),
    // then everything will be in the same cluster and the best pair will be the two best independent mappings
    vector<vector<pair<vector<size_t>, size_t>>> all_clusters = clusterer.cluster_seeds(seeds_by_read, distance_limit, 
            fragment_length_distr.mean() + 2*fragment_length_distr.stdev());
            //TODO: Choose a good distance for the fragment distance limit ^

    //For each fragment cluster, determine if it has clusters from both reads
    size_t max_fragment_num = 0;
    for (pair<vector< size_t>, size_t>& cluster : all_clusters[0]) {
        max_fragment_num = std::max(max_fragment_num, cluster.second);
    }
    for (pair<vector< size_t>, size_t>& cluster : all_clusters[1]) {
        max_fragment_num = std::max(max_fragment_num, cluster.second);
    }
#ifdef debug
    cerr << "Found " << max_fragment_num << " fragment clusters" << endl;
#endif
    vector<bool> has_first_read (max_fragment_num+1, false);//For each fragment cluster, does it have a cluster for the first read
    vector<bool> fragment_cluster_has_pair (max_fragment_num+1, false);//Does a fragment cluster have both reads
    bool found_paired_cluster = false;
    for (pair<vector<size_t>, size_t>& cluster : all_clusters[0]) {
        size_t fragment_num = cluster.second;
        has_first_read[fragment_num] = true;
    }
    for (pair<vector<size_t>, size_t>& cluster : all_clusters[1]) {
        size_t fragment_num = cluster.second;
        fragment_cluster_has_pair[fragment_num] = has_first_read[fragment_num];
        if (has_first_read[fragment_num]) {
            found_paired_cluster = true;
#ifdef debug
            cerr << "Fragment cluster " << fragment_num << " has read clusters from both reads" << endl;
#endif
        }
    }

    if (track_provenance) {
        funnels[0].substage("score");
        funnels[1].substage("score");
    }

    //For each fragment cluster (cluster of clusters), for each read, a vector of all alignments + the order they were fed into the funnel 
    //so the funnel can track them
    vector<pair<vector<Alignment>, vector<Alignment>>> alignments;
    vector<pair<vector<size_t>, vector<size_t>>> alignment_indices;
    pair<int, int> best_alignment_scores (0, 0); // The best alignment score for each end


    //Scores and coverage of each of the clusters
    pair<vector<double>, vector<double>> cluster_scores;
    pair<vector<double>, vector<double>>  cluster_coverages;


    //Keep track of the best cluster score and coverage per end for each fragment cluster
    pair<vector<double>, vector<double>> cluster_score_by_fragment;
    cluster_score_by_fragment.first.resize(max_fragment_num + 1, 0.0);
    cluster_score_by_fragment.second.resize(max_fragment_num + 1, 0.0);
    pair<vector<double>, vector<double>> cluster_coverage_by_fragment;
    cluster_coverage_by_fragment.first.resize(max_fragment_num + 1, 0.0);
    cluster_coverage_by_fragment.second.resize(max_fragment_num + 1, 0.0);

    for (size_t read_num = 0 ; read_num < 2 ; read_num++) {

        Alignment& aln = read_num == 0 ? aln1 : aln2;
        vector<size_t>& seed_to_source = read_num == 0 ? seed_to_source_by_read.first : seed_to_source_by_read.second;
        vector<pair<vector<size_t>, size_t>>& clusters = all_clusters[read_num];
        std::vector<Minimizer>& minimizers = read_num == 0 ? minimizers_by_read.first : minimizers_by_read.second;
        vector<pos_t>& seeds = seeds_by_read[read_num];
        vector<double>& best_cluster_score = read_num == 0 ? cluster_score_by_fragment.first : cluster_score_by_fragment.second;
        vector<double>& best_cluster_coverage = read_num == 0 ? cluster_coverage_by_fragment.first : cluster_coverage_by_fragment.second;

        // Cluster score is the sum of minimizer scores.
        vector<double>& cluster_score = read_num == 0 ? cluster_scores.first : cluster_scores.second;
        cluster_score.resize(clusters.size(), 0.0);
        vector<double>& read_coverage_by_cluster = read_num == 0 ? cluster_coverages.first : cluster_coverages.second;
        read_coverage_by_cluster.reserve(clusters.size());

        for (size_t i = 0; i < clusters.size(); i++) {
            // For each cluster
            auto& cluster = clusters[i].first;

            if (track_provenance) {
                // Say we're making it
                funnels[read_num].producing_output(i);
            }

            // Which minimizers are present in the cluster.
            vector<bool> present(minimizers.size(), false);
            for (auto hit_index : cluster) {
                present[seed_to_source[hit_index]] = true;
            }

            // Compute the score.
            for (size_t j = 0; j < minimizers.size(); j++) {
                if (present[j]) {
                    cluster_score[i] += minimizers[j].score;
                }
            }
            best_cluster_score[clusters[i].second] = max(best_cluster_score[clusters[i].second], cluster_score[i]);

            if (track_provenance) {
                // Record the cluster in the funnel as a group of the size of the number of items.
                funnels[read_num].merge_group(cluster.begin(), cluster.end());
                funnels[read_num].score(funnels[read_num].latest(), cluster_score[i]);

                // Say we made it.
                funnels[read_num].produced_output();
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
            best_cluster_coverage[clusters[i].second] = max(best_cluster_coverage[clusters[i].second], read_coverage_by_cluster.back());

        }
    }

    //TODO: Might not want to use this
    //For each fragment cluster, we want to know how many equivalent or better clusters we found
    vector<size_t> fragment_cluster_indices_by_score (max_fragment_num + 1);
    for (size_t i = 0 ; i < fragment_cluster_indices_by_score.size() ; i++) {
        fragment_cluster_indices_by_score[i] = i;
    }
    std::sort(fragment_cluster_indices_by_score.begin(), fragment_cluster_indices_by_score.end(), [&](size_t a, size_t b) {
        return cluster_coverage_by_fragment.first[a] + cluster_coverage_by_fragment.second[a] + cluster_score_by_fragment.first[a] + cluster_score_by_fragment.second[a]  
            > cluster_coverage_by_fragment.first[b] + cluster_coverage_by_fragment.second[b] + cluster_score_by_fragment.first[b] + cluster_score_by_fragment.second[b];  
    });

    vector<size_t> better_cluster_count (max_fragment_num+1); // How many fragment clusters are at least as good as the one at each index
    for (int j = fragment_cluster_indices_by_score.size() - 1 ; j >= 0 ; j--) {
        size_t i = fragment_cluster_indices_by_score[j];
        if (j == fragment_cluster_indices_by_score.size()-1) {
            better_cluster_count[i] = j;
        } else {
            size_t i2 = fragment_cluster_indices_by_score[j+1];
            if(cluster_coverage_by_fragment.first[i] + cluster_coverage_by_fragment.second[i] + cluster_score_by_fragment.first[i] + cluster_score_by_fragment.second[i] 
                == cluster_coverage_by_fragment.first[i2] + cluster_coverage_by_fragment.second[i2] + cluster_score_by_fragment.first[i2] + cluster_score_by_fragment.second[i2]) {
                better_cluster_count[i] = better_cluster_count[i2];
            } else {
                better_cluster_count[i] = j;
            }
        }
    }

    //Now that we've scored each of the clusters, extend and align them
    for (size_t read_num = 0 ; read_num < 2 ; read_num++) {

        Alignment& aln = read_num == 0 ? aln1 : aln2;
        vector<size_t>& seed_to_source = read_num == 0 ? seed_to_source_by_read.first : seed_to_source_by_read.second;
        vector<pair<vector<size_t>, size_t>>& clusters = all_clusters[read_num];
        std::vector<Minimizer>& minimizers = read_num == 0 ? minimizers_by_read.first : minimizers_by_read.second;
        vector<pos_t>& seeds = seeds_by_read[read_num];
        vector<double>& cluster_score = read_num == 0 ? cluster_scores.first : cluster_scores.second;
        vector<double>& read_coverage_by_cluster = read_num == 0 ? cluster_coverages.first : cluster_coverages.second;

#ifdef debug
        cerr << "Found " << clusters.size() << " clusters for read " << read_num << endl;
#endif

        // Retain clusters only if their score is better than this, in addition to the coverage cutoff
        double cluster_score_cutoff = cluster_score.size() == 0 ? 0 :
            *std::max_element(cluster_score.begin(), cluster_score.end()) - cluster_score_threshold;
        double cluster_coverage_cutoff = read_coverage_by_cluster.size() == 0 ? 0 :
                    *std::max_element(read_coverage_by_cluster.begin(), read_coverage_by_cluster.end())
                                    - cluster_coverage_threshold;

        if (track_provenance) {
            // Now we go from clusters to gapless extensions
            funnels[read_num].stage("extend");
        }

        // These are the GaplessExtensions for all the clusters (and fragment cluster assignments), in cluster_indexes_in_order order.
        vector<pair<vector<GaplessExtension>, size_t>> cluster_extensions;
        cluster_extensions.reserve(clusters.size());
        //TODO: Maybe put this back 
        //For each cluster, what fraction of "equivalent" clusters did we keep?
        //vector<vector<double>> probability_cluster_lost;
        //What is the score and coverage we are considering and how many reads
        //size_t curr_coverage = 0;
        //size_t curr_score = 0;
        //size_t curr_kept = 0;
        //size_t curr_count = 0;
        
        //Could just turn off these filters
        //Process clusters sorted by both score and read coverage
        process_until_threshold(clusters, read_coverage_by_cluster,
            [&](size_t a, size_t b) -> bool {
                //Sort clusters first by whether it was paired, then by the best coverage and score of any pair in the fragment cluster, 
                //then by its coverage and score
                size_t fragment_a = clusters[a].second;
                size_t fragment_b = clusters[b].second;

                double coverage_a = cluster_coverage_by_fragment.first[fragment_a]+cluster_coverage_by_fragment.second[fragment_a];
                double coverage_b = cluster_coverage_by_fragment.first[fragment_b]+cluster_coverage_by_fragment.second[fragment_b];
                double score_a = cluster_score_by_fragment.first[fragment_a]+cluster_score_by_fragment.second[fragment_a];
                double score_b = cluster_score_by_fragment.first[fragment_b]+cluster_score_by_fragment.second[fragment_b];

                if (fragment_cluster_has_pair[fragment_a] != fragment_cluster_has_pair[fragment_b]) {
                    return fragment_cluster_has_pair[fragment_a];
                } else if (coverage_a != coverage_b){
                    return coverage_a > coverage_b;
                } else if (score_a != score_b) {
                    return score_a > score_b;
                } else if (read_coverage_by_cluster[a] != read_coverage_by_cluster[b]){
                    return read_coverage_by_cluster[a] > read_coverage_by_cluster[b];
                } else {
                    return cluster_score[a] > cluster_score[b];
                }
            },
            0, 1, max_extensions,
            [&](size_t cluster_num) {
                // Handle sufficiently good clusters 
                
                if (!found_paired_cluster || fragment_cluster_has_pair[clusters[cluster_num].second] || 
                    (read_coverage_by_cluster[cluster_num] == cluster_coverage_cutoff + cluster_coverage_threshold &&
                           cluster_score[cluster_num] == cluster_score_cutoff + cluster_score_threshold)) { 
                    //If this cluster has a pair or if we aren't looking at pairs
                    //Or if it is the best cluster
                    
                    // First check against the additional score filter
                    if (cluster_coverage_threshold != 0 && read_coverage_by_cluster[cluster_num] < cluster_coverage_cutoff) {
                        //If the score isn't good enough, ignore this cluster
                        if (track_provenance) {
                            funnels[read_num].fail("cluster-coverage", cluster_num, read_coverage_by_cluster[cluster_num]);
                        }
                        return false;
                    }
                    if (cluster_score_threshold != 0 && cluster_score[cluster_num] < cluster_score_cutoff) {
                        //If the score isn't good enough, ignore this cluster
                        if (track_provenance) {
                            funnels[read_num].pass("cluster-coverage", cluster_num, read_coverage_by_cluster[cluster_num]);
                            funnels[read_num].pass("max-extensions", cluster_num);
                            funnels[read_num].fail("cluster-score", cluster_num, cluster_score[cluster_num]);
                        }
                        return false;
                    }
                    if (track_provenance) {
                        funnels[read_num].pass("cluster-coverage", cluster_num, read_coverage_by_cluster[cluster_num]);
                        funnels[read_num].pass("max-extensions", cluster_num);
                        funnels[read_num].pass("cluster-score", cluster_num, cluster_score[cluster_num]);
                        funnels[read_num].pass("paired-clusters", cluster_num);

                        funnels[read_num].processing_input(cluster_num);
                    }
                    vector<size_t>& cluster = clusters[cluster_num].first;

#ifdef debug
                    cerr << "Cluster " << cluster_num << endl;
#endif
                     
                    // Pack the seeds for GaplessExtender.
                    GaplessExtender::cluster_type seed_matchings;
                    for (auto& seed_index : cluster) {
                        // Insert the (graph position, read offset) pair.
                        seed_matchings.insert(GaplessExtender::to_seed(seeds[seed_index], minimizers[seed_to_source[seed_index]].value.offset));
#ifdef debug
                        cerr << "Seed read:" << minimizers[seed_to_source[seed_index]].offset << " = " << seeds[seed_index]
                            << " from minimizer " << seed_to_source[seed_index] << "(" << minimizer_index.count(minimizers[seed_to_source[seed_index]]) << ")" << endl;
#endif
                    }
                    
                    // Extend seed hits in the cluster into one or more gapless extensions
                    cluster_extensions.emplace_back(std::move(extender.extend(seed_matchings, aln.sequence())), 
                                                    clusters[cluster_num].second);
                    
                    if (track_provenance) {
                        // Record with the funnel that the previous group became a group of this size.
                        // Don't bother recording the seed to extension matching...
                        funnels[read_num].project_group(cluster_num, cluster_extensions.back().first.size());
                        
                        // Say we finished with this cluster, for now.
                        funnels[read_num].processed_input();
                    }
                    return true;
                } else {
                    //We were looking for clusters in a paired fragment cluster but this one doesn't have any on the other end
                    if (track_provenance) {
                        funnels[read_num].pass("cluster-coverage", cluster_num, read_coverage_by_cluster[cluster_num]);
                        funnels[read_num].pass("max-extensions", cluster_num);
                        funnels[read_num].pass("cluster-score", cluster_num, cluster_score[cluster_num]);
                        funnels[read_num].fail("paired-clusters", cluster_num);
                    }
                    return false;
                }
                
            }, [&](size_t cluster_num) {
                // There are too many sufficiently good clusters
                if (track_provenance) {
                    funnels[read_num].pass("cluster-coverage", cluster_num, read_coverage_by_cluster[cluster_num]);
                    funnels[read_num].fail("max-extensions", cluster_num);
                }
            }, [&](size_t cluster_num) {
                // This cluster is not sufficiently good.
                // TODO: I'm not sure if this is what should be failing
                if (track_provenance) {
                    funnels[read_num].fail("cluster-coverage", cluster_num, read_coverage_by_cluster[cluster_num]);
                }
            });
            
        
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
        alignments.resize(max_fragment_num + 2);
        alignment_indices.resize(max_fragment_num + 2);


        
        // Clear any old refpos annotation and path
        aln.clear_refpos();
        aln.clear_path();
        aln.set_score(0);
        aln.set_identity(0);
        aln.set_mapping_quality(0);
        
        //Since we will lose the order in which we pass alignments to the funnel, use this to keep track
        size_t curr_funnel_index = 0;

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
                    double identity = best_alignment.sequence().size() == 0 ? 0.0 
                            : (best_alignment.sequence().size() - mismatch_count) / (double) best_alignment.sequence().size();
                    
                    // Fill in the score and identity
                    best_alignment.set_score(extensions.front().score);
                    best_alignment.set_identity(identity);

                    if (extensions.size() > 1) {
                        //Do the same thing for the second extension, if one exists
                        *second_best_alignment.mutable_path() = extensions.back().to_path(gbwt_graph, second_best_alignment.sequence());
                        size_t mismatch_count = extensions.back().mismatches();
                        double identity = second_best_alignment.sequence().size() == 0 ? 0.0 
                                : (second_best_alignment.sequence().size() - mismatch_count) / (double) second_best_alignment.sequence().size();
                        
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
                
                
                size_t fragment_num = cluster_extensions[extension_num].second;
                if (second_best_alignment.score() != 0 && 
                    second_best_alignment.score() > best_alignment.score() * 0.8) {
                    //If there is a second extension and its score is at least half of the best score
                    read_num == 0 ? best_alignment_scores.first = max(best_alignment_scores.first, second_best_alignment.score())
                                  : best_alignment_scores.second = max(best_alignment_scores.second, second_best_alignment.score());
                    read_num == 0 ? alignments[fragment_num ].first.emplace_back(std::move(second_best_alignment) ) :
                                    alignments[fragment_num ].second.emplace_back(std::move(second_best_alignment));
                    read_num == 0 ? alignment_indices[fragment_num].first.emplace_back(curr_funnel_index)
                                  : alignment_indices[fragment_num].second.emplace_back(curr_funnel_index);
                        curr_funnel_index++;

                    if (track_provenance) {
        
                        funnels[read_num].project(extension_num);
                        read_num == 0 ? funnels[read_num].score(extension_num, alignments[fragment_num ].first.back().score()) :
                                        funnels[read_num].score(extension_num, alignments[fragment_num].second.back().score());
                    }
                }

                read_num == 0 ? best_alignment_scores.first = max(best_alignment_scores.first, best_alignment.score())
                              : best_alignment_scores.second = max(best_alignment_scores.second, best_alignment.score());
                read_num == 0 ? alignments[fragment_num].first.emplace_back(std::move(best_alignment))
                              : alignments[fragment_num].second.emplace_back(std::move(best_alignment));
                read_num == 0 ? alignment_indices[fragment_num].first.emplace_back(curr_funnel_index)
                              : alignment_indices[fragment_num].second.emplace_back(curr_funnel_index);

                curr_funnel_index++; 

                if (track_provenance) {

                    funnels[read_num].project(extension_num);
                    read_num == 0 ? funnels[read_num].score(extension_num, alignments[fragment_num].first.back().score())
                                  : funnels[read_num].score(extension_num, alignments[fragment_num].second.back().score());
                    
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
        funnels[0].stage("pairing");
        funnels[1].stage("pairing");
    }
    // Fill this in with the pairs of alignments we will output
    // each alignment is stored as <fragment index, alignment index> into alignments
    vector<pair<pair<size_t, size_t>, pair<size_t, size_t>>> paired_alignments;
    paired_alignments.reserve(alignments.size());
    //For each alignment in alignments, which paired_alignment includes it. Follows structure of alignments
    vector<pair<vector<vector<size_t>>, vector<vector<size_t>>>> alignment_groups(alignments.size());

    // Grab all the scores in order for MAPQ computation.
    vector<double> paired_scores;
    paired_scores.reserve(alignments.size());
    vector<int64_t> fragment_distances;
    fragment_distances.reserve(alignments.size());

    //For each fragment cluster, get the fraction of equivalent or better clusters that got thrown away

    vector<size_t> better_cluster_count_alignment_pairs; 
    better_cluster_count_alignment_pairs.reserve(alignments.size());

    //Keep track of alignments with no pairs in the same fragment cluster
    bool found_pair = false;
    //Alignments that don't have a mate
    // <fragment index, alignment_index, true if its the first end> 
    vector<tuple<size_t, size_t, bool>> unpaired_alignments;

    for (size_t fragment_num = 0 ; fragment_num < alignments.size() ; fragment_num ++ ) {
        //Get pairs of plausible alignments
        alignment_groups[fragment_num].first.resize(alignments[fragment_num].first.size());
        alignment_groups[fragment_num].second.resize(alignments[fragment_num].second.size());
        
        pair<vector<Alignment>, vector<Alignment>>& fragment_alignments = alignments[fragment_num];
        if (!fragment_alignments.first.empty() && ! fragment_alignments.second.empty()) {
            //Only keep pairs of alignments that were in the same fragment cluster
            found_pair = true;
            for (size_t i1 = 0 ; i1 < fragment_alignments.first.size() ; i1++)  {
                Alignment& alignment1 = fragment_alignments.first[i1];
                size_t j1 = alignment_indices[fragment_num].first[i1];
                for (size_t i2 = 0 ; i2 < fragment_alignments.second.size() ; i2++) {
                    Alignment& alignment2 = fragment_alignments.second[i2];
                    size_t j2 = alignment_indices[fragment_num].second[i2];

                    //Get the likelihood of the fragment distance
                    int64_t fragment_distance = distance_between(alignment1, alignment2); 
                    double dev = fragment_distance - fragment_length_distr.mean();
                    double fragment_length_log_likelihood = -dev * dev / (2.0 * fragment_length_distr.stdev() * fragment_length_distr.stdev());
                    if (fragment_distance != std::numeric_limits<int64_t>::max() ) {
                        double score = alignment1.score() + alignment2.score() + (fragment_length_log_likelihood / get_aligner()->log_base);
                        alignment_groups[fragment_num].first[i1].emplace_back(paired_alignments.size());
                        alignment_groups[fragment_num].second[i2].emplace_back(paired_alignments.size());
                        paired_alignments.emplace_back(make_pair(fragment_num, i1), make_pair(fragment_num, i2));
                        paired_scores.emplace_back(score);
                        fragment_distances.emplace_back(fragment_distance);
                        better_cluster_count_alignment_pairs.emplace_back(better_cluster_count[fragment_num]);

#ifdef debug
        cerr << "Found pair of alignments from fragment " << fragment_num << " with scores " 
             << alignment1.score() << " " << alignment2.score() << " at distance " << fragment_distance 
             << " gets pair score " << score << endl;
        cerr << "Alignment 1: " << pb2json(alignment1) << endl << "Alignment 2: " << pb2json(alignment2) << endl;
#endif
                    }

                    if (track_provenance) {
                        funnels[0].processing_input(j1);
                        funnels[1].processing_input(j2);
                        funnels[0].substage("pair-clusters");
                        funnels[1].substage("pair-clusters");
                        funnels[0].pass("max-rescue-attempts", j1);
                        funnels[0].project(j1);
                        funnels[1].pass("max-rescue-attempts", j2);
                        funnels[1].project(j2);
                        funnels[0].substage_stop();
                        funnels[1].substage_stop();
                        funnels[0].processed_input();
                        funnels[1].processed_input();
                    }
                }
            }
        } else if (!fragment_alignments.first.empty()) {
#ifdef debug
            cerr << "Found unpaired alignments from fragment " << fragment_num << " for first read" << endl;
#endif
            for (size_t i = 0 ; i < fragment_alignments.first.size() ; i++) {
                unpaired_alignments.emplace_back(fragment_num, i, true);
#ifdef debug
                cerr << "\t" << pb2json(fragment_alignments.first[i]) << endl;
#endif
            }
        } else if (!fragment_alignments.second.empty()) {
#ifdef debug
            cerr << "Found unpaired alignments from fragment " << fragment_num << " for second read" << endl;
#endif
            for (size_t i = 0 ; i < fragment_alignments.second.size() ; i++) {
                unpaired_alignments.emplace_back(fragment_num, i, false);
#ifdef debug
                cerr << "\t" << pb2json(fragment_alignments.second[i]) << endl;
#endif
            }
        }
    }

    if (!unpaired_alignments.empty()) {
        //If we didn't find any pairs within the fragment clusters, return the best individual alignments
        if (!found_pair && max_rescue_attempts == 0 ) {
            //If we aren't attempting rescue, just return the best for each

#ifdef debug
            cerr << "Found no pairs and we aren't doing rescue: return best alignment for each read" << endl;
#endif
            Alignment& best_aln1 = aln1;
            Alignment& best_aln2 = aln2;

            best_aln1.clear_refpos();
            best_aln1.clear_path();
            best_aln1.set_score(0);
            best_aln1.set_identity(0);
            best_aln1.set_mapping_quality(0);

            best_aln2.clear_refpos();
            best_aln2.clear_path();
            best_aln2.set_score(0);
            best_aln2.set_identity(0);
            best_aln2.set_mapping_quality(0);

            for (tuple<size_t, size_t, bool> index : unpaired_alignments ) {
                Alignment& alignment = std::get<2>(index) ? alignments[std::get<0>(index)].first[std::get<1>(index)]
                                                          : alignments[std::get<0>(index)].second[std::get<1>(index)];


                if (std::get<2>(index)) {
                    if (alignment.score() > best_aln1.score()) {
                        best_aln1 = alignment;
                    }
                } else {
                    if (alignment.score() > best_aln2.score()) {
                        best_aln2 = alignment;
                    }
                }
               
            }
            set_annotation(best_aln1, "unpaired", true);
            set_annotation(best_aln2, "unpaired", true);

            pair<vector<Alignment>, vector<Alignment>> paired_mappings;
            paired_mappings.first.emplace_back(std::move(best_aln1));
            paired_mappings.second.emplace_back(std::move(best_aln2));
            // Flip aln2 back to input orientation
            reverse_complement_alignment_in_place(&paired_mappings.second.back(), [&](vg::id_t node_id) {
                return gbwt_graph.get_length(gbwt_graph.get_handle(node_id));
            });
            // TODO: Maybe don't just give them the same mapq

            paired_mappings.first.back().set_mapping_quality(1);
            paired_mappings.second.back().set_mapping_quality(1);

            // Stop this alignment
            funnels[0].stop();
            funnels[1].stop();
            
            if (track_provenance) {
            
                // Annotate with the number of results in play at each stage
                for (size_t read_num = 0 ; read_num < 2 ; read_num++) {
                    funnels[read_num].for_each_stage([&](const string& stage, const vector<size_t>& result_sizes) {
                        // Save the number of items
                        set_annotation(read_num == 0 ? paired_mappings.first[0] : paired_mappings.second[0], "stage_" + stage + "_results", (double)result_sizes.size());
                    });
                    if (track_correctness) {
                        // And with the last stage at which we had any descendants of the correct seed hit locations
                        set_annotation(read_num == 0 ? paired_mappings.first[0] : paired_mappings.second[0], "last_correct_stage", funnels[read_num].last_correct_stage());
                    }
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
            }
            return paired_mappings;
        } else {
            
            //Attempt rescue

            process_until_threshold(unpaired_alignments, (std::function<double(size_t)>) [&](size_t i) -> double{
                tuple<size_t, size_t, bool> index = unpaired_alignments.at(i);
                return (double) std::get<2>(index) ? alignments[std::get<0>(index)].first[std::get<1>(index)].score()
                                                   : alignments[std::get<0>(index)].second[std::get<1>(index)].score();
            }, 0, 1, max_rescue_attempts, [&](size_t i) {
                //TODO: How many rescue attempts or how do we decide which ones?
                tuple<size_t, size_t, bool> index = unpaired_alignments.at(i);
                bool found_first = std::get<2>(index); 
                size_t j = found_first ? alignment_indices[std::get<0>(index)].first[std::get<1>(index)]
                                        : alignment_indices[std::get<0>(index)].second[std::get<1>(index)];
                if (track_provenance) {
                    funnels[found_first ? 0 : 1].processing_input(j);
                    funnels[found_first ? 0 : 1].substage("rescue");
                }
                Alignment& mapped_aln = found_first ? alignments[std::get<0>(index)].first[std::get<1>(index)]
                                                    : alignments[std::get<0>(index)].second[std::get<1>(index)];
                Alignment rescued_aln = found_first ? aln2 : aln1;
                rescued_aln.clear_path();

                if (found_pair && (double) mapped_aln.score() < (double) (found_first ? best_alignment_scores.first : best_alignment_scores.second) * 0.9) {
                    //TODO: 0.9?
                    //If this is not the best alignment we found for this end, do nothing
                    return true;
                }

                found_first ? attempt_rescue(mapped_aln, rescued_aln, true ) : 
                                attempt_rescue(mapped_aln, rescued_aln, false); 

                int64_t fragment_dist = found_first ? distance_between(mapped_aln, rescued_aln) 
                                                      : distance_between(rescued_aln, mapped_aln);
                if (fragment_dist != std::numeric_limits<int64_t>::max()) {
                    bool duplicated = false;

                    double dev = fragment_dist - fragment_length_distr.mean();
                    double fragment_length_log_likelihood = -dev * dev / (2.0 * fragment_length_distr.stdev() * fragment_length_distr.stdev());
                    double score = mapped_aln.score() + rescued_aln.score() + (fragment_length_log_likelihood / get_aligner()->log_base);

                    set_annotation(mapped_aln, "rescuer", true);
                    set_annotation(rescued_aln, "rescued", true);
                    set_annotation(mapped_aln,  "fragment_length", (double)fragment_dist);
                    set_annotation(rescued_aln, "fragment_length", (double)fragment_dist);
                    pair<size_t, size_t> mapped_index (std::get<0>(index), std::get<1>(index)); 
                    pair<size_t, size_t> rescued_index (alignments.size() - 1, 
                                found_first ? alignments.back().second.size() : alignments.back().first.size());
                    found_first ? alignments.back().second.emplace_back(std::move(rescued_aln)) 
                                : alignments.back().first.emplace_back(std::move(rescued_aln));

                    found_first ? alignment_groups.back().second.emplace_back() : alignment_groups.back().first.emplace_back();
                    pair<pair<size_t, size_t>, pair<size_t, size_t>> index_pair = found_first ? 
                                make_pair(mapped_index, rescued_index) : make_pair(rescued_index, mapped_index);
                    paired_alignments.push_back(index_pair);
                    paired_scores.emplace_back(score);
                    fragment_distances.emplace_back(fragment_dist);
                    better_cluster_count_alignment_pairs.emplace_back(0);
                    if (track_provenance) {
                        funnels[found_first ? 0 : 1].pass("max-rescue-attempts", j);
                        funnels[found_first ? 0 : 1].project(j);
                        funnels[found_first ? 1 : 0].introduce();
                    }
                } 
                if (track_provenance) {
                    funnels[found_first ? 0 : 1].processed_input();
                    funnels[found_first ? 0 : 1].substage_stop();
                }
                return true;
            }, [&](size_t i) {
                if (track_provenance) {
                    tuple<size_t, size_t, bool> index = unpaired_alignments.at(i);
                    bool found_first = std::get<2>(index); 
                    size_t j = found_first ? alignment_indices[std::get<0>(index)].first[std::get<1>(index)]
                                            : alignment_indices[std::get<0>(index)].second[std::get<1>(index)];
                    funnels[found_first ? 0 : 1].fail("max-rescue-attempts", j);
                }
                return;
            }, [&] (size_t i) {
                if (track_provenance) {
                    tuple<size_t, size_t, bool> index = unpaired_alignments.at(i);
                    bool found_first = std::get<2>(index); 
                    size_t j = found_first ? alignment_indices[std::get<0>(index)].first[std::get<1>(index)]
                                            : alignment_indices[std::get<0>(index)].second[std::get<1>(index)];
                }
                return;
            });
        }
    }

    
    
    if (track_provenance) {
        // Now say we are finding the winner(s)
        funnels[0].stage("winner");
        funnels[1].stage("winner");
    }
    // Fill this in with the alignments we will output
    pair<vector<Alignment>, vector<Alignment>> mappings;
    // Grab all the scores in order for MAPQ computation.
    vector<double> scores;
    vector<double> scores_group_1;
    vector<double> scores_group_2;
    vector<int64_t> distances;
    mappings.first.reserve(paired_alignments.size());
    mappings.second.reserve(paired_alignments.size());
    scores.reserve(paired_scores.size());
    distances.reserve(fragment_distances.size());
    vector<size_t> better_cluster_count_mappings;
    better_cluster_count_mappings.reserve(better_cluster_count_alignment_pairs.size());

    process_until_threshold(paired_alignments, (std::function<double(size_t)>) [&](size_t i) -> double {
        return paired_scores[i];
    }, 0, 1, max_multimaps, [&](size_t alignment_num) {
        // This alignment makes it
        // Called in score order

        pair<pair<size_t, size_t>, pair<size_t, size_t>> index_pair = paired_alignments[alignment_num];
        
        // Remember the score at its rank
        scores.emplace_back(paired_scores[alignment_num]);
        distances.emplace_back(fragment_distances[alignment_num]);
        // Remember the output alignment
        mappings.first.emplace_back( alignments[index_pair.first.first].first[index_pair.first.second]);
        mappings.second.emplace_back(alignments[index_pair.second.first].second[index_pair.second.second]);

        better_cluster_count_mappings.emplace_back(better_cluster_count_alignment_pairs[alignment_num]);
        if (mappings.first.size() == 1 && found_pair) {
            //If this is the best pair of alignments that we're going to return and we didn't attempt rescue, 
            //get the group scores for mapq

            //Get the scores of 
            scores_group_1.push_back(paired_scores[alignment_num]);
            scores_group_2.push_back(paired_scores[alignment_num]);

            //The indices (into paired_alignments) of pairs with the same first read as this
            vector<size_t>& alignment_group_1 = alignment_groups[index_pair.first.first].first[index_pair.first.second];
            vector<size_t>& alignment_group_2 = alignment_groups[index_pair.second.first].second[index_pair.second.second];

            for (size_t other_alignment_num : alignment_group_1) {
                if (other_alignment_num != alignment_num) {
                    scores_group_1.push_back(paired_scores[other_alignment_num]);
                }
            }
            for (size_t other_alignment_num : alignment_group_2) {
                if (other_alignment_num != alignment_num) {
                    scores_group_2.push_back(paired_scores[other_alignment_num]);
                }
            }
        }



        // Flip aln2 back to input orientation
        reverse_complement_alignment_in_place(&mappings.second.back(), [&](vg::id_t node_id) {
            return gbwt_graph.get_length(gbwt_graph.get_handle(node_id));
        });
        if (mappings.first.size() > 1) {
            mappings.first.back().set_is_secondary(true);
            mappings.second.back().set_is_secondary(true);
        }
        
        if (track_provenance) {
            // Tell the funnel
            funnels[0].pass("max-multimaps", alignment_num);
            funnels[0].project(alignment_num);
            funnels[0].score(alignment_num, scores.back());
            funnels[1].pass("max-multimaps", alignment_num);
            funnels[1].project(alignment_num);
            funnels[1].score(alignment_num, scores.back());
        }
        
        return true;
    }, [&](size_t alignment_num) {
        // We already have enough alignments, although this one has a good score
        
        // Remember the score at its rank anyway
        scores.emplace_back(paired_scores[alignment_num]);
        distances.emplace_back(fragment_distances[alignment_num]);
        
        if (track_provenance) {
            funnels[0].fail("max-multimaps", alignment_num);
            funnels[1].fail("max-multimaps", alignment_num);
        }
    }, [&](size_t alignment_num) {
        // This alignment does not have a sufficiently good score
        // Score threshold is 0; this should never happen
        assert(false);
    });
    
    if (track_provenance) {
        funnels[0].substage("mapq");
        funnels[1].substage("mapq");
    }
 
    if (mappings.first.empty()) {
        //If we didn't get an alignment, return empty alignments
        mappings.first.emplace_back(aln1);
        mappings.second.emplace_back(aln2);

        // Flip aln2 back to input orientation
        reverse_complement_alignment_in_place(&mappings.second.back(), [&](vg::id_t node_id) {
                return gbwt_graph.get_length(gbwt_graph.get_handle(node_id));
                });

        mappings.first.back().clear_refpos();
        mappings.first.back().clear_path();
        mappings.first.back().set_score(0);
        mappings.first.back().set_identity(0);
        mappings.first.back().set_mapping_quality(0);
        mappings.first.back().clear_refpos();
        mappings.first.back().clear_path();
        mappings.first.back().set_score(0);
        mappings.first.back().set_identity(0);
        mappings.first.back().set_mapping_quality(0);

    } else {
    
    #ifdef debug
        cerr << "For scores ";
        for (auto& score : scores) cerr << score << " ";
    #endif

    
        size_t winning_index;
        // Compute MAPQ if not unmapped. Otherwise use 0 instead of the 50% this would give us.
        // If either of the mappings was duplicated in other pairs, use the group scores to determine mapq
        double mapq = (mappings.first.empty() || scores[0] == 0) ? 0 : 
            get_regular_aligner()->maximum_mapping_quality_exact(scores, &winning_index) / 2;

        //Cap mapq at probability that there was a better fragment cluster that had the correct mapping
        if (better_cluster_count_mappings.size() != 0 && better_cluster_count_mappings.front() > 0) {
            mapq = min(mapq,round(prob_to_phred((1.0 / (double) better_cluster_count_mappings.front()))));
        }

        double mapq_group1 = scores_group_1.size() <= 1 ? mapq : 
            min(mapq, get_regular_aligner()->maximum_mapping_quality_exact(scores_group_1, &winning_index) / 2);
        double mapq_group2 = scores_group_2.size() <= 1 ? mapq : 
            min(mapq, get_regular_aligner()->maximum_mapping_quality_exact(scores_group_2, &winning_index) / 2);
        
    #ifdef debug
        cerr << "MAPQ is " << mapq << ", group MAPQ scores are " << mapq_group1 << " and " << mapq_group2 << endl;
    #endif

        mappings.first.front().set_mapping_quality(max(min(mapq_group1, 60.0), 0.0)) ;
        mappings.second.front().set_mapping_quality(max(min(mapq_group2, 60.0), 0.0)) ;
    
        //Annotate top pair with its fragment distance, fragment length distrubution, and secondary scores
        set_annotation(mappings.first.front(), "fragment_length", (double) distances.front());
        set_annotation(mappings.second.front(), "fragment_length", (double) distances.front());
        string distribution = "-I " + to_string(fragment_length_distr.mean()) + " -D " + to_string(fragment_length_distr.stdev());
        set_annotation(mappings.first.front(),"fragment_length_distribution", distribution);
        set_annotation(mappings.second.front(),"fragment_length_distribution", distribution);
        set_annotation(mappings.first.front(),"secondary_scores", scores);
        set_annotation(mappings.second.front(),"secondary_scores", scores);
    
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
                set_annotation(read_num == 0 ? mappings.first[0] : mappings.second[0], "stage_" + stage + "_results", (double)result_sizes.size());
            });
            if (track_correctness) {
                // And with the last stage at which we had any descendants of the correct seed hit locations
                set_annotation(read_num == 0 ? mappings.first[0] : mappings.second[0], "last_correct_stage", funnels[read_num].last_correct_stage());
            }
            // Annotate with the performances of all the filters
            // We need to track filter number
            size_t filter_num = 0;
            funnels[read_num].for_each_filter([&](const string& stage, const string& filter,
                const Funnel::FilterPerformance& by_count, const Funnel::FilterPerformance& by_size,
                const vector<double>& filter_statistics_correct, const vector<double>& filter_statistics_non_correct) {
                
                string filter_id = to_string(filter_num) + "_" + filter + "_" + stage;
                
                if (read_num == 0) {
                    // Save the stats
                    set_annotation(mappings.first[0], "filter_" + filter_id + "_passed_count_total", (double) by_count.passing);
                    set_annotation(mappings.first[0], "filter_" + filter_id + "_failed_count_total", (double) by_count.failing);
                    set_annotation(mappings.first[0], "filter_" + filter_id + "_passed_size_total", (double) by_size.passing);
                    set_annotation(mappings.first[0], "filter_" + filter_id + "_failed_size_total", (double) by_size.failing);
                    
                    if (track_correctness) {
                        set_annotation(mappings.first[0], "filter_" + filter_id + "_passed_count_correct", (double) by_count.passing_correct);
                        set_annotation(mappings.first[0], "filter_" + filter_id + "_failed_count_correct", (double) by_count.failing_correct);
                        set_annotation(mappings.first[0], "filter_" + filter_id + "_passed_size_correct", (double) by_size.passing_correct);
                        set_annotation(mappings.first[0], "filter_" + filter_id + "_failed_size_correct", (double) by_size.failing_correct);
                    }
                    
                    // Save the correct and non-correct filter statistics, even if
                    // everything is non-correct because correctness isn't computed
                    set_annotation(mappings.first[0], "filterstats_" + filter_id + "_correct", filter_statistics_correct);
                    set_annotation(mappings.first[0], "filterstats_" + filter_id + "_noncorrect", filter_statistics_non_correct);
                
                } else {
                    // Save the stats
                    set_annotation(mappings.second[0], "filter_" + filter_id + "_passed_count_total", (double) by_count.passing);
                    set_annotation(mappings.second[0], "filter_" + filter_id + "_failed_count_total", (double) by_count.failing);
                    set_annotation(mappings.second[0], "filter_" + filter_id + "_passed_size_total", (double) by_size.passing);
                    set_annotation(mappings.second[0], "filter_" + filter_id + "_failed_size_total", (double) by_size.failing);
                    
                    if (track_correctness) {
                        set_annotation(mappings.second[0], "filter_" + filter_id + "_passed_count_correct", (double) by_count.passing_correct);
                        set_annotation(mappings.second[0], "filter_" + filter_id + "_failed_count_correct", (double) by_count.failing_correct);
                        set_annotation(mappings.second[0], "filter_" + filter_id + "_passed_size_correct", (double) by_size.passing_correct);
                        set_annotation(mappings.second[0], "filter_" + filter_id + "_failed_size_correct", (double) by_size.failing_correct);
                    }
                    
                    // Save the correct and non-correct filter statistics, even if
                    // everything is non-correct because correctness isn't computed
                    set_annotation(mappings.second[0], "filterstats_" + filter_id + "_correct", filter_statistics_correct);
                    set_annotation(mappings.second[0], "filterstats_" + filter_id + "_noncorrect", filter_statistics_non_correct);
                }
                filter_num++;
            });
        }
        
        
        
        // Annotate with parameters used for the filters.
        set_annotation(mappings.first[0] , "param_hit-cap", (double) hit_cap);
        set_annotation(mappings.first[0] , "param_hard-hit-cap", (double) hard_hit_cap);
        set_annotation(mappings.first[0] , "param_score-fraction", (double) minimizer_score_fraction);
        set_annotation(mappings.first[0] , "param_max-extensions", (double) max_extensions);
        set_annotation(mappings.first[0] , "param_max-alignments", (double) max_alignments);
        set_annotation(mappings.first[0] , "param_cluster-score", (double) cluster_score_threshold);
        set_annotation(mappings.first[0] , "param_cluster-coverage", (double) cluster_coverage_threshold);
        set_annotation(mappings.first[0] , "param_extension-set", (double) extension_set_score_threshold);
        set_annotation(mappings.first[0] , "param_max-multimaps", (double) max_multimaps);
        set_annotation(mappings.first[0] , "param_max-rescue-attempts", (double) max_rescue_attempts);
        set_annotation(mappings.second[0], "param_hit-cap", (double) hit_cap);
        set_annotation(mappings.second[0], "param_hard-hit-cap", (double) hard_hit_cap);
        set_annotation(mappings.second[0], "param_score-fraction", (double) minimizer_score_fraction);
        set_annotation(mappings.second[0], "param_max-extensions", (double) max_extensions);
        set_annotation(mappings.second[0], "param_max-alignments", (double) max_alignments);
        set_annotation(mappings.second[0], "param_cluster-score", (double) cluster_score_threshold);
        set_annotation(mappings.second[0], "param_cluster-coverage", (double) cluster_coverage_threshold);
        set_annotation(mappings.second[0], "param_extension-set", (double) extension_set_score_threshold);
        set_annotation(mappings.second[0], "param_max-multimaps", (double) max_multimaps);
        set_annotation(mappings.second[0] , "param_max-rescue-attempts", (double) max_rescue_attempts);

    }
    
    // Ship out all the aligned alignments
    return mappings;

#ifdef debug
    // Dump the funnel info graph.
    funnels[0].to_dot(cerr);
    funnels[1].to_dot(cerr);
#endif
}

void MinimizerMapper::attempt_rescue( const Alignment& aligned_read, Alignment& rescued_alignment,  bool rescue_forward) {
    

    //Get the subgraph of all nodes within a reasonable range from aligned_read
    SubHandleGraph sub_graph(path_graph);
    //TODO: maybe should be gbwt_graph
    //TODO: How big should the rescue subgraph be?
    int64_t min_distance = max(0.0, fragment_length_distr.mean() - rescued_alignment.sequence().size() - 4*fragment_length_distr.stdev());
    int64_t max_distance = fragment_length_distr.mean() + 4*fragment_length_distr.stdev();
    distance_index.subgraphInRange(aligned_read.path(), path_graph, min_distance, max_distance, sub_graph, rescue_forward); 


    //Convert subgraph to directed, acyclic graph
    //Borrowed heavily from mpmap
    bdsg::HashGraph align_graph;
    unordered_map<id_t, pair<id_t, bool> > node_trans = algorithms::split_strands(&sub_graph, &align_graph);
    if (!algorithms::is_directed_acyclic(&sub_graph)) {

        bdsg::HashGraph dagified;
        unordered_map<id_t, id_t> dagify_trans = algorithms::dagify(&align_graph, &dagified, rescued_alignment.sequence().size());
        align_graph = move(dagified);
        node_trans = overlay_node_translations(dagify_trans, node_trans);
    }

    //Align to the subgraph
    rescued_alignment.clear_path();
    get_regular_aligner()->align(rescued_alignment, align_graph, true, false);

    translate_oriented_node_ids(*rescued_alignment.mutable_path(), node_trans);

    //TODO: mpmap also checks the score here
    return;

}
int64_t MinimizerMapper::distance_between(const Alignment& aln1, const Alignment& aln2) {
    assert(aln1.path().mapping_size() != 0); 
    assert(aln2.path().mapping_size() != 0); 
     
    pos_t pos1 = initial_position(aln1.path()); 
    pos_t pos2 = final_position(aln2.path());

    int64_t min_dist = distance_index.minDistance(pos1, pos2);
    return min_dist == -1 ? numeric_limits<int64_t>::max() : min_dist;
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


