/**
 * \file minimizer_mapper.cpp
 * Defines the code for the minimizer-and-GBWT-based mapper.
 */

#include "minimizer_mapper.hpp"

#include "annotation.hpp"
#include "path_subgraph.hpp"
#include "multipath_alignment.hpp"
#include "split_strand_graph.hpp"

#include "algorithms/dagify.hpp"
#include "algorithms/dijkstra.hpp"

#include <bdsg/overlays/strand_split_overlay.hpp>
#include <gbwtgraph/algorithms.h>
#include <gbwtgraph/cached_gbwtgraph.h>

#include <iostream>
#include <algorithm>
#include <cmath>

//#define debug
//#define print_minimizers

namespace vg {

using namespace std;

MinimizerMapper::MinimizerMapper(const gbwtgraph::GBWTGraph& graph,
    const std::vector<gbwtgraph::DefaultMinimizerIndex*>& minimizer_indexes,
    MinimumDistanceIndex& distance_index, const PathPositionHandleGraph* path_graph) :
    path_graph(path_graph), minimizer_indexes(minimizer_indexes),
    distance_index(distance_index), gbwt_graph(graph),
    extender(gbwt_graph, *(get_regular_aligner())), clusterer(distance_index),
    fragment_length_distr(1000,1000,0.95) {

   
}

//-----------------------------------------------------------------------------

void MinimizerMapper::map(Alignment& aln, AlignmentEmitter& alignment_emitter) {
    // Ship out all the aligned alignments
    alignment_emitter.emit_mapped_single(map(aln));
}

vector<Alignment> MinimizerMapper::map(Alignment& aln) {
    
#ifdef debug
    cerr << "Read " << aln.name() << ": " << aln.sequence() << endl;
#endif
#ifdef print_minimizers
    cerr << aln.sequence() << "\t";
    for (char c : aln.quality()) {
        cerr << (char)(c+33);
    }
#endif

    // Make a new funnel instrumenter to watch us map this read.
    Funnel funnel;
    funnel.start(aln.name());
    
    // Minimizers sorted by score in descending order.
    std::vector<Minimizer> minimizers = this->find_minimizers(aln.sequence(), funnel);

    // Find the seeds and mark the minimizers that were located.
    std::vector<Seed> seeds = this->find_seeds(minimizers, aln, funnel);

    // Cluster the seeds. Get sets of input seed indexes that go together.
    if (track_provenance) {
        funnel.stage("cluster");
    }
    std::vector<Cluster> clusters = clusterer.cluster_seeds(seeds, get_distance_limit(aln.sequence().size()));

    // Determine the scores and read coverages for each cluster.
    // Also find the best and second-best cluster scores.
    if (this->track_provenance) {
        funnel.substage("score");
    }
    double best_cluster_score = 0.0, second_best_cluster_score = 0.0;
    for (size_t i = 0; i < clusters.size(); i++) {
        Cluster& cluster = clusters[i];
        this->score_cluster(cluster, i, minimizers, seeds, aln.sequence().length(), funnel);
        if (cluster.score > best_cluster_score) {
            second_best_cluster_score = best_cluster_score;
            best_cluster_score = cluster.score;
        } else if (cluster.score > second_best_cluster_score) {
            second_best_cluster_score = cluster.score;
        }
    }

#ifdef debug
    cerr << "Found " << clusters.size() << " clusters" << endl;
#endif
#ifdef print_minimizers
    cerr << "\t" << clusters.size();
#endif
    
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
        // Now we go from clusters to gapless extensions
        funnel.stage("extend");
    }
    
    // These are the GaplessExtensions for all the clusters.
    vector<vector<GaplessExtension>> cluster_extensions;
    cluster_extensions.reserve(clusters.size());
    // To compute the windows for explored minimizers, we need to get
    // all the minimizers that are explored.
    SmallBitset minimizer_explored(minimizers.size());
    //How many hits of each minimizer ended up in each extended cluster?
    vector<vector<size_t>> minimizer_extended_cluster_count; 
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
    
    //Process clusters sorted by both score and read coverage
    process_until_threshold_c<Cluster, double>(clusters, [&](size_t i) -> double {
            return clusters[i].coverage;
        }, [&](size_t a, size_t b) -> bool {
            return ((clusters[a].coverage > clusters[b].coverage) ||
                    (clusters[a].coverage == clusters[b].coverage && clusters[a].score > clusters[b].score));
        },
        cluster_coverage_threshold, 1, max_extensions,
        [&](size_t cluster_num) {
            // Handle sufficiently good clusters in descending coverage order
            
            Cluster& cluster = clusters[cluster_num];
            if (track_provenance) {
                funnel.pass("cluster-coverage", cluster_num, cluster.coverage);
                funnel.pass("max-extensions", cluster_num);
            }
            
            // First check against the additional score filter
            if (cluster_score_threshold != 0 && cluster.score < cluster_score_cutoff) {
                //If the score isn't good enough, ignore this cluster
                if (track_provenance) {
                    funnel.fail("cluster-score", cluster_num, cluster.score);
                }
                
                // Record MAPQ implications of not extending this cluster.
                unextended_clusters.push_back(cluster_num);
#ifdef debug
            cerr << "Cluster " << cluster_num << " fails cluster score cutoff" <<  endl;
            cerr << "Covers " << clusters[cluster_num].coverage << "/best-" << cluster_coverage_threshold << " of read" << endl;
            cerr << "Scores " << clusters[cluster_num].score << "/" << cluster_score_cutoff << endl;
#endif
                return false;
            }
            
            if (track_provenance) {
                funnel.pass("cluster-score", cluster_num, cluster.score);
                funnel.processing_input(cluster_num);
            }
            if (cluster.coverage == curr_coverage &&
                cluster.score == curr_score &&
                curr_kept < max_extensions * 0.75) {
                curr_kept++;
                curr_count++;
            } else if (cluster.coverage != curr_coverage ||
                    cluster.score != curr_score) {
                    //If this is a cluster that has scores different than the previous one
                    for (size_t i = 0 ; i < curr_kept ; i++ ) {
                        probability_cluster_lost.push_back(1.0 - (double(curr_kept) / double(curr_count)));
                    }
                    curr_coverage = cluster.coverage;
                    curr_score = cluster.score;
                    curr_kept = 1;
                    curr_count = 1;
            } else {
                //If this cluster is equivalent to the previous one and we already took enough
                //equivalent clusters
                curr_count ++;
                
                // TODO: shouldn't we fail something for the funnel here?
                
                // Record MAPQ implications of not extending this cluster.
                unextended_clusters.push_back(cluster_num);
#ifdef debug
            cerr << "Cluster " << cluster_num << " fails because we took too many identical clusters" <<  endl;
            cerr << "Covers " << clusters[cluster_num].coverage << "/best-" << cluster_coverage_threshold << " of read" << endl;
            cerr << "Scores " << clusters[cluster_num].score << "/" << cluster_score_cutoff << endl;
#endif
                return false;
            }
            

            //Only keep this cluster if we have few enough equivalent clusters

#ifdef debug
            cerr << "Cluster " << cluster_num << endl;
            cerr << "Covers " << cluster.coverage << "/best-" << cluster_coverage_threshold << " of read" << endl;
            cerr << "Scores " << cluster.score << "/" << cluster_score_cutoff << endl;
#endif
             
            minimizer_extended_cluster_count.emplace_back(minimizers.size(), 0);
            // Pack the seeds for GaplessExtender.
            GaplessExtender::cluster_type seed_matchings;
            for (auto seed_index : cluster.seeds) {
                // Insert the (graph position, read offset) pair.
                const Seed& seed = seeds[seed_index];
                seed_matchings.insert(GaplessExtender::to_seed(seed.pos, minimizers[seed.source].value.offset));
                minimizer_extended_cluster_count.back()[seed.source]++;
#ifdef debug
                const Minimizer& minimizer = minimizers[seed.source];
                cerr << "Seed read:" << minimizer.value.offset << " = " << seed.pos
                    << " from minimizer " << seed.source << "(" << minimizer.hits << ")" << endl;
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
            Cluster& cluster = clusters[cluster_num];
            if (track_provenance) {
                funnel.pass("cluster-coverage", cluster_num, cluster.coverage);
                funnel.fail("max-extensions", cluster_num);
            }
            if (cluster.coverage == curr_coverage &&
                cluster.score == curr_score) {
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
            unextended_clusters.push_back(cluster_num);
#ifdef debug
            cerr << "Cluster " << cluster_num << " passes cluster cutoffs but we have too many" <<  endl;
            cerr << "Covers " << cluster.coverage << "/best-" << cluster_coverage_threshold << " of read" << endl;
            cerr << "Scores " << cluster.score << "/" << cluster_score_cutoff << endl;
#endif
            
        }, [&](size_t cluster_num) {
            // This cluster is not sufficiently good.
            if (track_provenance) {
                funnel.fail("cluster-coverage", cluster_num, clusters[cluster_num].coverage);
            }
            for (size_t i = 0 ; i < curr_kept ; i++ ) {
                probability_cluster_lost.push_back(1.0 - (double(curr_kept) / double(curr_count)));
            }
            curr_kept = 0;
            curr_count = 0;
            curr_score = 0;
            curr_coverage = 0;
            
            // Record MAPQ implications of not extending this cluster.
            unextended_clusters.push_back(cluster_num);
#ifdef debug
            cerr << "Cluster " << cluster_num << " fails cluster coverage cutoffs" <<  endl;
            cerr << "Covers " << clusters[cluster_num].coverage << "/best-" << cluster_coverage_threshold << " of read" << endl;
            cerr << "Scores " << clusters[cluster_num].score << "/" << cluster_score_cutoff << endl;
#endif
        });
        
    for (size_t i = 0 ; i < curr_kept ; i++ ) {
        probability_cluster_lost.push_back(1.0 - (double(curr_kept) / double(curr_count)));
    }
    std::vector<int> cluster_extension_scores = this->score_extensions(cluster_extensions, aln, funnel);

    if (track_provenance) {
        funnel.stage("align");
    }

    //How many of each minimizer ends up in an extension set that actually gets turned into an alignment?
    vector<size_t> minimizer_extensions_count(minimizers.size(), 0);
    
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
    process_until_threshold_b(cluster_extensions, cluster_extension_scores,
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

            // Determine if we got full-length alignments or local alignments.
            bool full_length_extensions = (extensions.size() <= 2);
            for (auto& extension : extensions) {
                full_length_extensions &= extension.full();
            }

            if (full_length_extensions) {
                // We got full-length extensions, so directly convert to an Alignment.
                
                if (track_provenance) {
                    funnel.substage("direct");
                }

                //Fill in the best alignments from the extension
                this->extension_to_alignment(extensions.front(), best_alignment);

                if (extensions.size() > 1) {
                    //Do the same thing for the second extension, if one exists
                    this->extension_to_alignment(extensions.back(), second_best_alignment);
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
#ifdef debug
                cerr << "Found second best alignment from gapless extension " << extension_num << ": " << pb2json(second_best_alignment) << endl;
#endif
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
#ifdef debug
                cerr << "Found best alignment from gapless extension " << extension_num << ": " << pb2json(best_alignment) << endl;
#endif

            alignments.push_back(std::move(best_alignment));
            alignments_to_source.push_back(extension_num);
            probability_alignment_lost.push_back(probability_cluster_lost[extension_num]);

            if (track_provenance) {

                funnel.project(extension_num);
                funnel.score(alignments.size() - 1, alignments.back().score());
                
                // We're done with this input item
                funnel.processed_input();
            }

            for (size_t i = 0 ; i < minimizer_extended_cluster_count[extension_num].size() ; i++) {
                minimizer_extensions_count[i] += minimizer_extended_cluster_count[extension_num][i];
                if (minimizer_extended_cluster_count[extension_num][i] > 0) {
                    // This minimizer is in an extended cluster that gave rise
                    // to at least one alignment, so it is explored.
                    minimizer_explored.insert(i);
                }
            }
            
            return true;
        }, [&](size_t extension_num) {
            // There are too many sufficiently good extensions
            if (track_provenance) {
                funnel.pass("extension-set", extension_num, cluster_extension_scores[extension_num]);
                funnel.fail("max-alignments", extension_num);
            }
#ifdef debug
                cerr << "gapless extension " << extension_num << " failed because there were too many good extensions" << endl;
#endif
        }, [&](size_t extension_num) {
            // This extension is not good enough.
            if (track_provenance) {
                funnel.fail("extension-set", extension_num, cluster_extension_scores[extension_num]);
            }
#ifdef debug
                cerr << "gapless extension " << extension_num << " failed because its score was not good enough (score=" << cluster_extension_scores[extension_num] << ")" << endl;
#endif
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
    process_until_threshold_a(alignments, (std::function<double(size_t)>) [&](size_t i) -> double {
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
    cerr << "Picked best alignment " << pb2json(mappings[0]) << endl;
    cerr << "For scores ";
    for (auto& score : scores) cerr << score << " ";
#endif

    size_t winning_index;
    assert(!mappings.empty());
    // Compute MAPQ if not unmapped. Otherwise use 0 instead of the 50% this would give us.
    double mapq = (mappings.front().path().mapping_size() == 0) ? 0 : 
        get_regular_aligner()->maximum_mapping_quality_exact(scores, &winning_index) / 2;

#ifdef print_minimizers
double uncapped_mapq = mapq;
#endif
#ifdef debug
    cerr << "uncapped MAPQ is " << mapq << endl;
#endif
    
    if (probability_mapping_lost.front() > 0) {
        mapq = min(mapq,round(prob_to_phred(probability_mapping_lost.front())));
    }
    
    // TODO: give SmallBitset iterators so we can use it instead of an index vector.
    vector<size_t> explored_minimizers;
    for (size_t i = 0; i < minimizers.size(); i++) {
        if (minimizer_explored.contains(i)) {
            explored_minimizers.push_back(i);
        }
    }
    // Compute caps on MAPQ. TODO: avoid needing to pass as much stuff along.
    double mapq_explored_cap = 2 * faster_cap(minimizers, explored_minimizers, aln.sequence(), aln.quality());

    // Remember the uncapped MAPQ and the caps
    set_annotation(mappings.front(), "mapq_uncapped", mapq);
    set_annotation(mappings.front(), "mapq_explored_cap", mapq_explored_cap);

    // Apply the caps
    mapq = min(mapq, mapq_explored_cap);

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
    
#ifdef print_minimizers
    for (size_t i = 0 ; i < minimizers.size() ; i++) {
        auto& minimizer = minimizers[i];
        cerr << "\t"
             << minimizer.value.key.decode(minimizer.length) << "\t"
             << minimizer.forward_offset() << "\t"
             << minimizer.agglomeration_start << "\t"
             << minimizer.agglomeration_length << "\t"
             << minimizer.hits << "\t"
             << minimizer_extensions_count[i];
         if (minimizer_extensions_count[i]>0) {
             assert(minimizer.hits<=hard_hit_cap) ;
         }
    }
    cerr << "\t" << uncapped_mapq << "\t" << mapq_explored_cap << "\t" << probability_mapping_lost.front();
    if (track_correctness) {
        cerr << "\t" << funnel.last_correct_stage() << endl;
    } else {
        cerr << endl;
    }
#endif
#ifdef debug
    // Dump the funnel info graph.
    funnel.to_dot(cerr);
#endif

    return mappings;
}

//-----------------------------------------------------------------------------

pair<vector<Alignment>, vector<Alignment>> MinimizerMapper::map_paired(Alignment& aln1, Alignment& aln2,
                                                      vector<pair<Alignment, Alignment>>& ambiguous_pair_buffer){
    if (fragment_length_distr.is_finalized()) {
        //If we know the fragment length distribution then we just map paired ended 
        return map_paired(aln1, aln2);

    } else {
        //If we don't know the fragment length distribution, map the reads single ended

        vector<Alignment> alns1(map(aln1));
        vector<Alignment> alns2(map(aln2));

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
    
    // Minimizers for both reads, sorted by score in descending order.
    std::vector<std::vector<Minimizer>> minimizers_by_read(2);
    minimizers_by_read[0] = this->find_minimizers(aln1.sequence(), funnels[0]);
    minimizers_by_read[1] = this->find_minimizers(aln2.sequence(), funnels[1]);

    // Seeds for both reads, stored in separate vectors.
    std::vector<std::vector<Seed>> seeds_by_read(2);
    seeds_by_read[0] = this->find_seeds(minimizers_by_read[0], aln1, funnels[0]);
    seeds_by_read[1] = this->find_seeds(minimizers_by_read[1], aln2, funnels[1]);

    // Cluster the seeds. Get sets of input seed indexes that go together.
    // If the fragment length distribution hasn't been fixed yet (if the expected fragment length = 0),
    // then everything will be in the same cluster and the best pair will be the two best independent mappings
    if (track_provenance) {
        funnels[0].stage("cluster");
        funnels[1].stage("cluster");
    }
    std::vector<std::vector<Cluster>> all_clusters = clusterer.cluster_seeds(seeds_by_read, get_distance_limit(aln1.sequence().size()), 
            fragment_length_distr.mean() + paired_distance_stdevs * fragment_length_distr.std_dev());

    //For each fragment cluster, determine if it has clusters from both reads
    size_t max_fragment_num = 0;
    for (auto& cluster : all_clusters[0]) {
        max_fragment_num = std::max(max_fragment_num, cluster.fragment);
    }
    for (auto& cluster : all_clusters[1]) {
        max_fragment_num = std::max(max_fragment_num, cluster.fragment);
    }
#ifdef debug
    cerr << "Found " << max_fragment_num << " fragment clusters" << endl;
#endif
    vector<bool> has_first_read (max_fragment_num+1, false);//For each fragment cluster, does it have a cluster for the first read
    vector<bool> fragment_cluster_has_pair (max_fragment_num+1, false);//Does a fragment cluster have both reads
    bool found_paired_cluster = false;
    for (auto& cluster : all_clusters[0]) {
        size_t fragment_num = cluster.fragment;
        has_first_read[fragment_num] = true;
    }
    for (auto& cluster : all_clusters[1]) {
        size_t fragment_num = cluster.fragment;
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

    //Keep track of the best cluster score and coverage per end for each fragment cluster
    pair<vector<double>, vector<double>> cluster_score_by_fragment;
    cluster_score_by_fragment.first.resize(max_fragment_num + 1, 0.0);
    cluster_score_by_fragment.second.resize(max_fragment_num + 1, 0.0);
    pair<vector<double>, vector<double>> cluster_coverage_by_fragment;
    cluster_coverage_by_fragment.first.resize(max_fragment_num + 1, 0.0);
    cluster_coverage_by_fragment.second.resize(max_fragment_num + 1, 0.0);

    for (size_t read_num = 0 ; read_num < 2 ; read_num++) {
        Alignment& aln = read_num == 0 ? aln1 : aln2;
        std::vector<Cluster>& clusters = all_clusters[read_num];
        std::vector<Minimizer>& minimizers = minimizers_by_read[read_num];
        std::vector<Seed>& seeds = seeds_by_read[read_num];
        vector<double>& best_cluster_score = read_num == 0 ? cluster_score_by_fragment.first : cluster_score_by_fragment.second;
        vector<double>& best_cluster_coverage = read_num == 0 ? cluster_coverage_by_fragment.first : cluster_coverage_by_fragment.second;

        for (size_t i = 0; i < clusters.size(); i++) {
            // Deterimine cluster score and read coverage.
            Cluster& cluster = clusters[i];
            this->score_cluster(cluster, i, minimizers, seeds, aln.sequence().length(), funnels[read_num]);
            size_t fragment = cluster.fragment;
            best_cluster_score[fragment] = std::max(best_cluster_score[fragment], cluster.score);
            best_cluster_coverage[fragment] = std::max(best_cluster_coverage[fragment], cluster.coverage);
        }
    }

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

    // We track unextended clusters.
    vector<vector<size_t>> unextended_clusters_by_read(2);
    // To compute the windows that are explored, we need to get
    // all the minimizers that are explored.
    vector<SmallBitset> minimizer_explored_by_read(2);
    //How many hits of each minimizer ended up in each extended cluster?
    vector<vector<vector<size_t>>> minimizer_extended_cluster_count_by_read(2); 
    

    //Now that we've scored each of the clusters, extend and align them
    for (size_t read_num = 0 ; read_num < 2 ; read_num++) {
        Alignment& aln = read_num == 0 ? aln1 : aln2;
        std::vector<Cluster>& clusters = all_clusters[read_num];
        std::vector<Minimizer>& minimizers = minimizers_by_read[read_num];
        std::vector<Seed>& seeds = seeds_by_read[read_num];

#ifdef debug
        cerr << "Found " << clusters.size() << " clusters for read " << read_num << endl;
#endif

        // Retain clusters only if their score is better than this, in addition to the coverage cutoff
        double cluster_score_cutoff = 0.0, cluster_coverage_cutoff = 0.0;
        for (auto& cluster : clusters) {
            cluster_score_cutoff = std::max(cluster_score_cutoff, cluster.score);
            cluster_coverage_cutoff = std::max(cluster_coverage_cutoff, cluster.coverage);
        }
        cluster_score_cutoff -= cluster_score_threshold;
        cluster_coverage_cutoff -= cluster_coverage_threshold;

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

        unextended_clusters_by_read[read_num].reserve(clusters.size());
        minimizer_explored_by_read[read_num] = SmallBitset(minimizers.size());
        
        //Process clusters sorted by both score and read coverage
        process_until_threshold_c<Cluster, double>(clusters, [&](size_t i) -> double {
                return clusters[i].coverage;
            }, [&](size_t a, size_t b) -> bool {
                //Sort clusters first by whether it was paired, then by the best coverage and score of any pair in the fragment cluster, 
                //then by its coverage and score
                size_t fragment_a = clusters[a].fragment;
                size_t fragment_b = clusters[b].fragment;

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
                } else if (clusters[a].coverage != clusters[b].coverage){
                    return clusters[a].coverage > clusters[b].coverage;
                } else {
                    return clusters[a].score > clusters[b].score;
                }
            },
            0, 1, max_extensions,
            [&](size_t cluster_num) {
                // Handle sufficiently good clusters 
                Cluster& cluster = clusters[cluster_num];
                if (!found_paired_cluster || fragment_cluster_has_pair[cluster.fragment] || 
                    (cluster.coverage == cluster_coverage_cutoff + cluster_coverage_threshold &&
                           cluster.score == cluster_score_cutoff + cluster_score_threshold)) { 
                    //If this cluster has a pair or if we aren't looking at pairs
                    //Or if it is the best cluster
                    
                    // First check against the additional score filter
                    if (cluster_coverage_threshold != 0 && cluster.coverage < cluster_coverage_cutoff) {
                        //If the score isn't good enough, ignore this cluster
                        if (track_provenance) {
                            funnels[read_num].fail("cluster-coverage", cluster_num, cluster.coverage);
                        }
                        unextended_clusters_by_read[read_num].push_back(cluster_num);
                        return false;
                    }
                    if (cluster_score_threshold != 0 && cluster.score < cluster_score_cutoff) {
                        //If the score isn't good enough, ignore this cluster
                        if (track_provenance) {
                            funnels[read_num].pass("cluster-coverage", cluster_num, cluster.coverage);
                            funnels[read_num].pass("max-extensions", cluster_num);
                            funnels[read_num].fail("cluster-score", cluster_num, cluster.score);
                        }
                        unextended_clusters_by_read[read_num].push_back(cluster_num);
                        return false;
                    }
                    if (track_provenance) {
                        funnels[read_num].pass("cluster-coverage", cluster_num, cluster.coverage);
                        funnels[read_num].pass("max-extensions", cluster_num);
                        funnels[read_num].pass("cluster-score", cluster_num, cluster.score);
                        funnels[read_num].pass("paired-clusters", cluster_num);

                        funnels[read_num].processing_input(cluster_num);
                    }

#ifdef debug
                    cerr << "Cluster " << cluster_num << endl;
#endif
                    
                    //Count how many of each minimizer is in each cluster extension
                    minimizer_extended_cluster_count_by_read[read_num].emplace_back(minimizers.size(), 0);
                    // Pack the seeds for GaplessExtender.
                    GaplessExtender::cluster_type seed_matchings;
                    for (auto seed_index : cluster.seeds) {
                        // Insert the (graph position, read offset) pair.
                        const Seed& seed = seeds[seed_index];
                        seed_matchings.insert(GaplessExtender::to_seed(seed.pos, minimizers[seed.source].value.offset));
                        minimizer_extended_cluster_count_by_read[read_num].back()[seed.source]++;
#ifdef debug
                        cerr << "Seed read:" << minimizers[seed.source].value.offset << " = " << seed.pos
                            << " from minimizer " << seed.source << endl;
#endif
                    }
                    
                    // Extend seed hits in the cluster into one or more gapless extensions
                    cluster_extensions.emplace_back(std::move(extender.extend(seed_matchings, aln.sequence())), 
                                                    cluster.fragment);
                    
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
                        funnels[read_num].pass("cluster-coverage", cluster_num, cluster.coverage);
                        funnels[read_num].pass("max-extensions", cluster_num);
                        funnels[read_num].pass("cluster-score", cluster_num, cluster.score);
                        funnels[read_num].fail("paired-clusters", cluster_num);
                    }
                    unextended_clusters_by_read[read_num].push_back(cluster_num);
                    return false;
                }
                
            }, [&](size_t cluster_num) {
                // There are too many sufficiently good clusters
                if (track_provenance) {
                    funnels[read_num].pass("cluster-coverage", cluster_num, clusters[cluster_num].coverage);
                    funnels[read_num].fail("max-extensions", cluster_num);
                }
                unextended_clusters_by_read[read_num].push_back(cluster_num);
            }, [&](size_t cluster_num) {
                // This cluster is not sufficiently good.
                // TODO: I don't think it should ever get here unless we limit the scores of the fragment clusters we look at
                unextended_clusters_by_read[read_num].push_back(cluster_num);
            });
            
        
        // We now estimate the best possible alignment score for each cluster.
        std::vector<int> cluster_extension_scores = this->score_extensions(cluster_extensions, aln, funnels[read_num]);
        
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
        process_until_threshold_b(cluster_extensions, cluster_extension_scores,
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
                

                // Determine if we got full-length alignments or local alignments.
                bool full_length_extensions = (extensions.size() <= 2);
                for (auto& extension : extensions) {
                    full_length_extensions &= extension.full();
                }

                if (full_length_extensions) {
                    // We got full-length extensions, so directly convert to an Alignment.
                    
                    if (track_provenance) {
                        funnels[read_num].substage("direct");
                    }

                    //Fill in the best alignments from the extension
                    this->extension_to_alignment(extensions.front(), best_alignment);

                    if (extensions.size() > 1) {
                        //Do the same thing for the second extension, if one exists
                        this->extension_to_alignment(extensions.back(), second_best_alignment);
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

                
                for (size_t i = 0 ; i < minimizer_extended_cluster_count_by_read[read_num][extension_num].size() ; i++) {
                    if (minimizer_extended_cluster_count_by_read[read_num][extension_num][i] > 0) {
                        // This minimizer is in an extended cluster that gave rise
                        // to at least one alignment, so it is explored.
                        minimizer_explored_by_read[read_num].insert(i);
                    }
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
#ifdef print_minimizers
    vector<pair<bool, bool>> alignment_was_rescued;
#endif
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
    size_t unpaired_count_1 = 0;
    size_t unpaired_count_2 = 0;

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
                    double fragment_length_log_likelihood = -dev * dev / (2.0 * fragment_length_distr.std_dev() * fragment_length_distr.std_dev());
                    if (fragment_distance != std::numeric_limits<int64_t>::max() ) {
                        double score = alignment1.score() + alignment2.score() + (fragment_length_log_likelihood / get_aligner()->log_base);
                        alignment_groups[fragment_num].first[i1].emplace_back(paired_alignments.size());
                        alignment_groups[fragment_num].second[i2].emplace_back(paired_alignments.size());
                        paired_alignments.emplace_back(make_pair(fragment_num, i1), make_pair(fragment_num, i2));
                        paired_scores.emplace_back(score);
                        fragment_distances.emplace_back(fragment_distance);
                        better_cluster_count_alignment_pairs.emplace_back(better_cluster_count[fragment_num]);
#ifdef print_minimizers
                        alignment_was_rescued.emplace_back(false, false);
#endif

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
                unpaired_count_1++;
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
                unpaired_count_2++;
#ifdef debug
                cerr << "\t" << pb2json(fragment_alignments.second[i]) << endl;
#endif
            }
        }
    }
    size_t rescued_count_1 = 0;
    size_t rescued_count_2 = 0;
    vector<bool> rescued_from;

    if (!unpaired_alignments.empty()) {
        //If we found some clusters that don't belong to a fragment cluster
        if (!found_pair && max_rescue_attempts == 0 ) {
            //If we didn't find any pairs and we aren't attempting rescue, just return the best for each end

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
#ifdef print_minimizers
        cerr << aln1.sequence() << "\t";
        for (char c : aln1.quality()) {
            cerr << (char)(c+33);
        }
        cerr << "\t" << all_clusters.size() << "\t" << 0 << "\t" << 0 << "\t?";
        for (size_t i = 0 ; i < minimizers_by_read[0].size() ; i++) {
            auto& minimizer = minimizers_by_read[0][i];
            cerr << "\t"
                 << minimizer.value.key.decode(minimizer.length) << "\t"
                 << minimizer.forward_offset() << "\t"
                 << minimizer.agglomeration_start << "\t"
                 << minimizer.agglomeration_length << "\t"
                 << minimizer.hits << "\t"
                 << minimizer_explored_by_read[0].contains(i);
             if (minimizer_explored_by_read[0].contains(i)) {
                 assert(minimizer.hits<=hard_hit_cap) ;
             }
        }
        cerr << "\t" << 1 << "\t" << 1 << "\t" << 1 << "\t" << 1;  
        if (track_correctness) {
            cerr << "\t" << funnels[0].last_correct_stage() << endl;
        } else {
            cerr << "\t?" << endl;
        }

        cerr << aln2.sequence() << "\t";
        for (char c : aln2.quality()) {
            cerr << (char)(c+33);
        }
        cerr << "\t" << all_clusters.size() << "\t" << 0 << "\t" << 0 << "\t?";
        for (size_t i = 0 ; i < minimizers_by_read[1].size() ; i++) {
            auto& minimizer = minimizers_by_read[1][i];
            cerr << "\t"
                 << minimizer.value.key.decode(minimizer.length) << "\t"
                 << minimizer.forward_offset() << "\t"
                 << minimizer.agglomeration_start << "\t"
                 << minimizer.agglomeration_length << "\t"
                 << minimizer.hits << "\t"
                 << minimizer_explored_by_read[0].contains(i);
             if (minimizer_explored_by_read[0].contains(i)) {
                 assert(minimizer.hits<=hard_hit_cap) ;
             }
        }
        cerr << "\t" << 1 << "\t" << 1 << "\t" << 1 << "\t" << 1;  
        if (track_correctness) {
            cerr << "\t" << funnels[1].last_correct_stage() << endl;
        } else {
            cerr << "\t?" << endl;
        }
#endif   
            return paired_mappings;
        } else {
            
            //Attempt rescue on unpaired alignments if either we didn't find any pairs or if the unpaired alignments are very good

            process_until_threshold_a(unpaired_alignments, (std::function<double(size_t)>) [&](size_t i) -> double{
                tuple<size_t, size_t, bool> index = unpaired_alignments.at(i);
                return (double) std::get<2>(index) ? alignments[std::get<0>(index)].first[std::get<1>(index)].score()
                                                   : alignments[std::get<0>(index)].second[std::get<1>(index)].score();
            }, 0, 1, max_rescue_attempts, [&](size_t i) {
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

                if (found_pair && (double) mapped_aln.score() < (double) (found_first ? best_alignment_scores.first : best_alignment_scores.second) * paired_rescue_score_limit) {
                    //If we have already found paired clusters and this unpaired alignment is not good enough, do nothing
                    return true;
                }

                attempt_rescue(mapped_aln, rescued_aln, minimizers_by_read[(found_first ? 1 : 0)], found_first);

                int64_t fragment_dist = found_first ? distance_between(mapped_aln, rescued_aln) 
                                                      : distance_between(rescued_aln, mapped_aln);
                if (fragment_dist != std::numeric_limits<int64_t>::max()) {
                    bool duplicated = false;

                    double dev = fragment_dist - fragment_length_distr.mean();
                    double fragment_length_log_likelihood = -dev * dev / (2.0 * fragment_length_distr.std_dev() * fragment_length_distr.std_dev());
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
                    found_first ? rescued_count_1++ : rescued_count_2++;

                    found_first ? alignment_groups.back().second.emplace_back() : alignment_groups.back().first.emplace_back();
                    pair<pair<size_t, size_t>, pair<size_t, size_t>> index_pair = found_first ? 
                                make_pair(mapped_index, rescued_index) : make_pair(rescued_index, mapped_index);
                    paired_alignments.push_back(index_pair);
                    paired_scores.emplace_back(score);
                    fragment_distances.emplace_back(fragment_dist);
                    better_cluster_count_alignment_pairs.emplace_back(0);
                    rescued_from.push_back(found_first); 
#ifdef print_minimizers
                    alignment_was_rescued.emplace_back(!found_first, found_first);
#endif
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

    double estimated_multiplicity_from_1 = unpaired_count_1 > 0 ? (double) unpaired_count_1 / min(rescued_count_1, max_rescue_attempts) : 1.0;
    double estimated_multiplicity_from_2 = unpaired_count_2 > 0 ? (double) unpaired_count_2 / min(rescued_count_2, max_rescue_attempts) : 1.0;
    vector<double> paired_multiplicities;
    for (bool rescued_from_first : rescued_from) {
        paired_multiplicities.push_back(rescued_from_first ? estimated_multiplicity_from_1 : estimated_multiplicity_from_2);
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

#ifdef print_minimizers
vector<pair<bool, bool>> mapping_was_rescued;
#endif

    process_until_threshold_a(paired_alignments, (std::function<double(size_t)>) [&](size_t i) -> double {
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

#ifdef print_minimizers
        mapping_was_rescued.emplace_back(alignment_was_rescued[alignment_num]);
#endif
        
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
#ifdef print_minimizers
    double uncapped_mapq = 0.0;
    double fragment_cluster_cap = 0.0;
    vector<double> mapq_extend_caps;
    double mapq_score_group_1 = 0.0;
    double mapq_score_group_2 = 0.0;
#endif
 
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

        mappings.second.back().clear_refpos();
        mappings.second.back().clear_path();
        mappings.second.back().set_score(0);
        mappings.second.back().set_identity(0);
        mappings.second.back().set_mapping_quality(0);
#ifdef print_minimizers
        mapping_was_rescued.emplace_back(false, false);
#endif

    } else {
    
#ifdef debug
        cerr << "For scores ";
        for (auto& score : scores) cerr << score << " ";
#endif

        size_t winning_index;
        // Compute MAPQ if not unmapped. Otherwise use 0 instead of the 50% this would give us.
        // If either of the mappings was duplicated in other pairs, use the group scores to determine mapq
        const vector<double>* multiplicities = paired_multiplicities.size() == scores.size() ? &paired_multiplicities : nullptr; 
        // Compute MAPQ if not unmapped. Otherwise use 0 instead of the 50% this would give us.
        // If either of the mappings was duplicated in other pairs, use the group scores to determine mapq
        double mapq = scores[0] == 0 ? 0 : 
            get_regular_aligner()->maximum_mapping_quality_exact(scores, &winning_index, multiplicities) / 2;
#ifdef print_minimizers
        uncapped_mapq = mapq;
#endif

        //Cap mapq at 1 / # equivalent or better fragment clusters
        if (better_cluster_count_mappings.size() != 0 && better_cluster_count_mappings.front() > 0) {
            mapq = min(mapq,round(prob_to_phred((1.0 / (double) better_cluster_count_mappings.front()))));
        }
#ifdef print_minimizers
        if (better_cluster_count_mappings.size() != 0 && better_cluster_count_mappings.front() > 0) {
            fragment_cluster_cap = round(prob_to_phred((1.0 / (double) better_cluster_count_mappings.front())));
        }
#endif

        //If one alignment was duplicated in other pairs, cap the mapq for that alignment at the mapq
        //of the group of duplicated alignments
        double mapq_group1 = scores_group_1.size() <= 1 ? mapq : 
            min(mapq, get_regular_aligner()->maximum_mapping_quality_exact(scores_group_1, &winning_index) / 2);
        double mapq_group2 = scores_group_2.size() <= 1 ? mapq : 
            min(mapq, get_regular_aligner()->maximum_mapping_quality_exact(scores_group_2, &winning_index) / 2);

    
#ifdef print_minimizers
    mapq_score_group_1 = scores_group_1.size() <= 1 ? 0.0 : get_regular_aligner()->maximum_mapping_quality_exact(scores_group_1, &winning_index);
    mapq_score_group_2 = scores_group_2.size() <= 1 ? 0.0 : get_regular_aligner()->maximum_mapping_quality_exact(scores_group_2, &winning_index);
#endif
        // Compute one MAPQ cap across all the fragments
        double mapq_cap = -std::numeric_limits<float>::infinity();
        for (auto read_num : {0, 1}) {
            // For each fragment

            // Find the source read
            auto& aln = read_num == 0 ? aln1 : aln2;
            
            vector<size_t> explored_minimizers;
            for (size_t i = 0; i < minimizers_by_read[read_num].size(); i++) {
                if (minimizer_explored_by_read[read_num].contains(i)) {
                    explored_minimizers.push_back(i);
                }
            }
            // Compute caps on MAPQ. TODO: avoid needing to pass as much stuff along.
            double mapq_explored_cap = 2 * faster_cap(minimizers_by_read[read_num], explored_minimizers, aln.sequence(), aln.quality());

#ifdef print_minimizers
            mapq_extend_caps.emplace_back(mapq_explored_cap);
#endif

            // Remember the caps
            auto& to_annotate = (read_num == 0 ? mappings.first : mappings.second).front();
            set_annotation(to_annotate, "mapq_explored_cap", mapq_explored_cap);

            // Compute the cap. It should be the higher of the caps for the two reads 
            // (unless one has no minimizers, i.e. if it was rescued) 
            // The individual cap values are either actual numbers or +inf, so the cap can't stay as -inf.
            mapq_cap = mapq_explored_cap == numeric_limits<double>::infinity() ? mapq_cap : max(mapq_cap, mapq_explored_cap);
        }
        mapq_cap = mapq_cap == -std::numeric_limits<float>::infinity() ? numeric_limits<double>::infinity() : mapq_cap;

        for (auto read_num : {0, 1}) {
            // For each fragment
            // Find the source read
            auto& aln = read_num == 0 ? aln1 : aln2;

            // Find the MAPQ to cap
            auto& read_mapq = read_num == 0 ? mapq_group1 : mapq_group2;
            
            // Remember the uncapped MAPQ
            auto& to_annotate = (read_num == 0 ? mappings.first : mappings.second).front();
            set_annotation(to_annotate, "mapq_uncapped", read_mapq);
            // And the cap we actually applied (possibly from the pair partner)
            set_annotation(to_annotate, "mapq_applied_cap", mapq_cap);

            // Apply the cap
            read_mapq = min(read_mapq, mapq_cap);
        }
        
#ifdef debug
        cerr << "MAPQ is " << mapq << ", capped group MAPQ scores are " << mapq_group1 << " and " << mapq_group2 << endl;
#endif

        mappings.first.front().set_mapping_quality(max(min(mapq_group1, 60.0), 0.0)) ;
        mappings.second.front().set_mapping_quality(max(min(mapq_group2, 60.0), 0.0)) ;
    
        //Annotate top pair with its fragment distance, fragment length distrubution, and secondary scores
        set_annotation(mappings.first.front(), "fragment_length", (double) distances.front());
        set_annotation(mappings.second.front(), "fragment_length", (double) distances.front());
        string distribution = "-I " + to_string(fragment_length_distr.mean()) + " -D " + to_string(fragment_length_distr.std_dev());
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
 
#ifdef print_minimizers
    if (mapq_extend_caps.size() == 0) {
        mapq_extend_caps.resize(2,0.0);
    }
    if (distances.size() == 0) {
        distances.emplace_back(0);
    }
    cerr << aln1.sequence() << "\t";
    for (char c : aln1.quality()) {
        cerr << (char)(c+33);
    }
    cerr << "\t" << all_clusters.size() << "\t" << mapping_was_rescued[0].first << "\t" << mapping_was_rescued[0].second << "\t" << distances.front();
    for (size_t i = 0 ; i < minimizers_by_read[0].size() ; i++) {
        auto& minimizer = minimizers_by_read[0][i];
        cerr << "\t"
             << minimizer.value.key.decode(minimizer.length) << "\t"
             << minimizer.forward_offset() << "\t"
             << minimizer.agglomeration_start << "\t"
             << minimizer.agglomeration_length << "\t"
             << minimizer.hits << "\t"
             << minimizer_explored_by_read[0].contains(i);
         if (minimizer_explored_by_read[0].contains(i)) {
             assert(minimizer.hits<=hard_hit_cap) ;
         }
    } cerr << "\t" << uncapped_mapq << "\t" << fragment_cluster_cap << "\t" << mapq_score_group_1 << "\t" << mapq_extend_caps[0];  
    if (track_correctness) {
        cerr << "\t" << funnels[0].last_correct_stage() << endl;
    } else {
        cerr << "\t?" << endl;
    }

    cerr << aln2.sequence() << "\t";
    for (char c : aln2.quality()) {
        cerr << (char)(c+33);
    }
    cerr << "\t" << all_clusters.size() << "\t" << mapping_was_rescued[0].second << "\t" << mapping_was_rescued[0].first << "\t" << distances.front();
    for (size_t i = 0 ; i < minimizers_by_read[1].size() ; i++) {
        auto& minimizer = minimizers_by_read[1][i];
        cerr << "\t"
             << minimizer.value.key.decode(minimizer.length) << "\t"
             << minimizer.forward_offset() << "\t"
             << minimizer.agglomeration_start << "\t"
             << minimizer.agglomeration_length << "\t"
             << minimizer.hits << "\t"
             << minimizer_explored_by_read[1].contains(i);
         if (minimizer_explored_by_read[1].contains(i)) {
             assert(minimizer.hits<=hard_hit_cap) ;
         }
    }
    cerr << "\t" << uncapped_mapq << "\t" << fragment_cluster_cap << "\t" << mapq_score_group_1 << "\t" << mapq_extend_caps[1];
    if (track_correctness) {
        cerr << "\t" << funnels[1].last_correct_stage() << endl;
    } else {
        cerr << "\t?" << endl;
    }
#endif

    // Ship out all the aligned alignments
    return mappings;

#ifdef debug
    // Dump the funnel info graph.
    funnels[0].to_dot(cerr);
    funnels[1].to_dot(cerr);
#endif
}

//-----------------------------------------------------------------------------

double MinimizerMapper::compute_mapq_caps(const Alignment& aln, 
    const std::vector<Minimizer>& minimizers,
    const SmallBitset& explored) {

    // We need to cap MAPQ based on the likelihood of generating all the windows in the explored minimizers by chance, too.
#ifdef debug
    cerr << "Cap based on explored minimizers all being faked by errors..." << endl;
#endif

    // Convert our flag vector to a list of the minimizers actually in extended clusters
    vector<size_t> explored_minimizers;
    explored_minimizers.reserve(minimizers.size());
    for (size_t i = 0; i < minimizers.size(); i++) {
        if (explored.contains(i)) {
            explored_minimizers.push_back(i);
        }
    }
    double mapq_explored_cap = window_breaking_quality(minimizers, explored_minimizers, aln.sequence(), aln.quality());
    
    return mapq_explored_cap;
}

double MinimizerMapper::window_breaking_quality(const vector<Minimizer>& minimizers, vector<size_t>& broken,
    const string& sequence, const string& quality_bytes) {
    
#ifdef debug
    cerr << "Computing MAPQ cap based on " << broken.size() << "/" << minimizers.size() << " minimizers' windows" << endl;
#endif

    if (broken.empty() || quality_bytes.empty()) {
        // If we have no agglomerations or no qualities, bail
        return numeric_limits<double>::infinity();
    }
    
    assert(sequence.size() == quality_bytes.size());

    // Sort the agglomerations by start position
    std::sort(broken.begin(), broken.end(), [&](const size_t& a, const size_t& b) {
        return minimizers[a].agglomeration_start < minimizers[b].agglomeration_start;
    });

    // Have a priority queue for tracking the agglomerations that a base is currently overlapping.
    // Prioritize earliest-ending best, and then latest-starting best.
    // This is less, and greater priority is at the front of the queue (better).
    // TODO: do we really need to care about the start here?
    auto agglomeration_priority = [&](const size_t& a, const size_t& b) {
        // Returns true if a is worse (ends later, or ends at the same place and starts earlier).
        auto& ma = minimizers[a];
        auto& mb = minimizers[b];
        auto a_end = ma.agglomeration_start + ma.agglomeration_length;
        auto b_end = mb.agglomeration_start + mb.agglomeration_length;
        return (a_end > b_end) || (a_end == b_end && ma.agglomeration_start < mb.agglomeration_start);
    };
    // We maintain our own heap so we can iterate over it.
    vector<size_t> active_agglomerations;

    // A window in flight is a pair of start position, inclusive end position
    struct window_t {
        size_t first;
        size_t last;
    };

    // Have a priority queue of window starts and ends, prioritized earliest-ending best, and then latest-starting best.
    // This is less, and greater priority is at the front of the queue (better).
    auto window_priority = [&](const window_t& a, const window_t& b) {
        // Returns true if a is worse (ends later, or ends at the same place and starts earlier).
        return (a.last > b.last) || (a.last == b.last && a.first < b.first);
    };
    priority_queue<window_t, vector<window_t>, decltype(window_priority)> active_windows(window_priority);
    
    // Have a cursor for which agglomeration should come in next.
    auto next_agglomeration = broken.begin();
    
    // Have a DP table with the cost of the cheapest solution to the problem up to here, including a hit at this base.
    // Or numeric_limits<double>::infinity() if base cannot be hit.
    // We pre-fill it because I am scared to use end() if we change its size.
    vector<double> costs(sequence.size(), numeric_limits<double>::infinity());
    
    // Keep track of the latest-starting window ending before here. If none, this will be two numeric_limits<size_t>::max() values.
    window_t latest_starting_ending_before = { numeric_limits<size_t>::max(), numeric_limits<size_t>::max() };
    // And keep track of the min cost DP table entry, or end if not computed yet.
    auto latest_starting_ending_before_winner = costs.end();
    
    for (size_t i = 0; i < sequence.size(); i++) {
        // For each base in the read
        
#ifdef debug
        cerr << "At base " << i << endl;
#endif

        // Bring in new agglomerations and windows that start at this base.
        while (next_agglomeration != broken.end() && minimizers[*next_agglomeration].agglomeration_start == i) {
            // While the next agglomeration starts here

            // Record it
            active_agglomerations.push_back(*next_agglomeration);
            std::push_heap(active_agglomerations.begin(), active_agglomerations.end(), agglomeration_priority);

            // Look it up
            auto& minimizer = minimizers[*next_agglomeration];
            
            // Determine its window size from its index.
            size_t window_size = minimizer.window_size();
            
#ifdef debug
            cerr << "\tBegin agglomeration of " << (minimizer.agglomeration_length - window_size + 1)
                << " windows of " << window_size << " bp each for minimizer "
                << *next_agglomeration << " (" << minimizer.forward_offset()
                << "-" << (minimizer.forward_offset() + minimizer.length) << ")" << endl;
#endif

            // Make sure the minimizer instance isn't too far into the agglomeration to actually be conatined in the same k+w-1 window as the first base.
            assert(minimizer.agglomeration_start + minimizer.candidates_per_window >= minimizer.forward_offset());

            for (size_t start = minimizer.agglomeration_start;
                start + window_size - 1 < minimizer.agglomeration_start + minimizer.agglomeration_length;
                start++) {
                // Add all the agglomeration's windows to the queue, looping over their start bases in the read.
                window_t add = {start, start + window_size - 1};
#ifdef debug
                cerr << "\t\t" << add.first << " - " << add.last << endl;
#endif
                active_windows.push(add);
            }
            
            // And advance the cursor
            ++next_agglomeration;
        }
        
        // We have the start and end of the latest-starting window ending before here (may be none)
        
        if (isATGC(sequence[i]) &&
            (!active_windows.empty() || 
            (next_agglomeration != broken.end() && minimizers[*next_agglomeration].agglomeration_start == i))) {
            
            // This base is not N, and it is either covered by an agglomeration
            // that hasn't ended yet, or a new agglomeration starts here.
            
#ifdef debug
            cerr << "\tBase is acceptable (" << sequence[i] << ", " << active_agglomerations.size() << " active agglomerations, "
                << active_windows.size() << " active windows)" << endl;
#endif
            
            // Score mutating the base itself, thus causing all the windows it touches to be wrong.
            // TODO: account for windows with multiple hits not having to be explained at full cost here.
            // We start with the cost to mutate the base.
            double base_cost = quality_bytes[i];

#ifdef debug
            cerr << "\t\tBase base quality: " << base_cost << endl;
#endif

            for (const size_t& active_index : active_agglomerations) {
                // Then, we look at each agglomeration the base overlaps

                // Find the minimizer whose agglomeration we are looking at
                auto& active = minimizers[active_index];

                if (i >= active.forward_offset() &&
                    i < active.forward_offset() + active.length) {
                    // If the base falls in the minimizer, we don't do anything. Just mutating the base is enough to create this wrong minimizer.
                    
#ifdef debug
                    cerr << "\t\tBase in minimizer " << active_index << " (" << active.forward_offset()
                        << "-" << (active.forward_offset() + active.length) << ") for agglomeration "
                        << active.agglomeration_start << "-" << (active.agglomeration_start + active.agglomeration_length) << endl;
#endif
                    
                    continue;
                }

                // If the base falls outside the minimizer, it participates in
                // some number of other possible minimizers in the
                // agglomeration. Compute that, accounting for edge effects.
                size_t possible_minimizers = min((size_t) active.length, min(i - active.agglomeration_start + 1, (active.agglomeration_start + active.agglomeration_length) - i));

                // Then for each of those possible minimizers, we need P(would have beaten the current minimizer).
                // We approximate this as constant across the possible minimizers.
                // And since we want to OR together beating, we actually AND together not-beating and then not it.
                // So we track the probability of not beating.
                double any_beat_phred = phred_for_at_least_one(active.value.hash, possible_minimizers);

#ifdef debug
                cerr << "\t\tBase flanks minimizer " << active_index << " (" << active.forward_offset()
                    << "-" << (active.forward_offset() + active.length) << ") for agglomeration "
                    << active.agglomeration_start << "-" << (active.agglomeration_start + active.agglomeration_length) 
                    << " and has " << possible_minimizers
                    << " chances to have beaten it; cost to have beaten with any is Phred " << any_beat_phred << endl;
#endif

                // Then we AND (sum) the Phred of that in, as an additional cost to mutating this base and kitting all the windows it covers.
                // This makes low-quality bases outside of minimizers produce fewer low MAPQ caps, and accounts for the robustness of minimizers to some errors.
                base_cost += any_beat_phred;

            }
                
            // Now we know the cost of mutating this base, so we need to
            // compute the total cost of a solution through here, mutating this
            // base.

            // Look at the start of that latest-starting window ending before here.
            
            if (latest_starting_ending_before.first == numeric_limits<size_t>::max()) {
                // If there is no such window, this is the first base hit, so
                // record the cost of hitting it.
                
                costs[i] = base_cost;
                
#ifdef debug
                cerr << "\tFirst base hit, costs " << costs[i] << endl;
#endif
            } else {
                // Else, scan from that window's start to its end in the DP
                // table, and find the min cost.

                if (latest_starting_ending_before_winner == costs.end()) {
                    // We haven't found the min in the window we come from yet, so do that.
                    latest_starting_ending_before_winner = std::min_element(costs.begin() + latest_starting_ending_before.first,
                        costs.begin() + latest_starting_ending_before.last + 1);
                }
                double min_prev_cost = *latest_starting_ending_before_winner;
                    
                // Add the cost of hitting this base
                costs[i] = min_prev_cost + base_cost;
                
#ifdef debug
                cerr << "\tComes from prev base at " << (latest_starting_ending_before_winner - costs.begin()) << ", costs " << costs[i] << endl;
#endif
            }
            
        } else {
            // This base is N, or not covered by an agglomeration.
            // Leave infinity there to say we can't hit it.
            // Nothing to do!
        }
        
        // Now we compute the start of the latest-starting window ending here or before, and deactivate agglomerations/windows.
        
        while (!active_agglomerations.empty() &&
            minimizers[active_agglomerations.front()].agglomeration_start + minimizers[active_agglomerations.front()].agglomeration_length - 1 == i) {
            // Look at the queue to see if an agglomeration ends here.

#ifdef debug
            cerr << "\tEnd agglomeration " << active_agglomerations.front() << endl;
#endif
            
            std::pop_heap(active_agglomerations.begin(), active_agglomerations.end(), agglomeration_priority);
            active_agglomerations.pop_back();
        }

        while (!active_windows.empty() && active_windows.top().last == i) {
            // The look at the queue to see if a window ends here. This is
            // after windows are added so that we can handle 1-base windows.
            
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
                // And clear our cache of the lowest cost base to hit it.
                latest_starting_ending_before_winner = costs.end();
                
#ifdef debug
                cerr << "\t\t\tNow have: " << latest_starting_ending_before.first << " - " << latest_starting_ending_before.last << endl;
#endif
                
            }
            
            // And pop them all off.
            active_windows.pop();
        }
        // If not, use the latest-starting window ending at the previous base or before (i.e. do nothing).
        
        // Loop around; we will have the latest-starting window ending before the next here.
    }
    
#ifdef debug
    cerr << "Final window to scan: " << latest_starting_ending_before.first << " - " << latest_starting_ending_before.last << endl;
#endif
    
    // When we get here, all the agglomerations should have been handled
    assert(next_agglomeration == broken.end());
    // And all the windows should be off the queue.
    assert(active_windows.empty());
    
    // When we get to the end, we have the latest-starting window overall. It must exist.
    assert(latest_starting_ending_before.first != numeric_limits<size_t>::max());
    
    // Scan it for the best final base to hit and return the cost there.
    // Don't use the cache here because nothing can come after the last-ending window.
    auto min_cost_at = std::min_element(costs.begin() + latest_starting_ending_before.first,
        costs.begin() + latest_starting_ending_before.last + 1);
#ifdef debug
    cerr << "Overall min cost: " << *min_cost_at << " at base " << (min_cost_at - costs.begin()) << endl;
#endif
    return *min_cost_at;
}

double MinimizerMapper::faster_cap(const vector<Minimizer>& minimizers, vector<size_t>& minimizers_explored,
    const string& sequence, const string& quality_bytes) {

    // Sort minimizer subset so we go through minimizers in increasing order of start position
    std::sort(minimizers_explored.begin(), minimizers_explored.end(), [&](size_t a, size_t b) {
        // Return true if a must come before b, and false otherwise
        return minimizers[a].forward_offset() < minimizers[b].forward_offset();
    });
#ifdef debug
    cerr << "Sorted " << minimizers_explored.size() << " minimizers" << endl;
#endif

#ifdef debug
    // Dump read and minimizers
    int digits_needed = (int) ceil(log10(sequence.size()));
    for (int digit = digits_needed - 1; digit >= 0; digit--) {
        for (size_t i = 0; i < sequence.size(); i++) {
            // Output the correct digit for this place in this number
            cerr << (char) ('0' + (uint8_t) round(i % (int) round(pow(10, digit + 1)) / pow(10, digit)));
        }
        cerr << endl;
    }
    cerr << sequence << endl;
    
    for (auto& explored_index : minimizers_explored) {
        // For each explored minimizer
        auto& m = minimizers[explored_index];
        for (size_t i = 0; i < m.agglomeration_start; i++) {
            // Space until its agglomeration starts
            cerr << ' ';
        }
        
        for (size_t i = m.agglomeration_start; i < m.forward_offset(); i++) {
            // Do the beginnign of the agglomeration
            cerr << '-';
        }
        // Do the minimizer itself
        cerr << m.value.key.decode(m.length);
        for (size_t i = m.forward_offset() + m.length ; i < m.agglomeration_start + m.agglomeration_length; i++) {
            // Do the tail end of the agglomeration
            cerr << '-';
        }
        cerr << endl;
    }
#endif

    // Make a DP table holding the log10 probability of having an error disrupt each minimizer.
    // Entry i+1 is log prob of mutating minimizers 0, 1, 2, ..., i.
    // Make sure to have an extra field at the end to support this.
    // Initialize with -inf for unfilled.
    vector<double> c(minimizers_explored.size() + 1, -numeric_limits<double>::infinity());
    c[0] = 0.0;
    
    for_each_aglomeration_interval(minimizers, sequence, quality_bytes, minimizers_explored, [&](size_t left, size_t right, size_t bottom, size_t top) {
        // For each overlap range in the agglomerations
        
#ifdef debug
        cerr << "Consider overlap range " << left << " to " << right << " in minimizer ranks " << bottom << " to " << top << endl;
        cerr << "log10prob for bottom: " << c[bottom] << endl;
#endif
        
        // Calculate the probability of a disruption here
        double p_here = get_log10_prob_of_disruption_in_interval(minimizers, sequence, quality_bytes,
            minimizers_explored.begin() + bottom, minimizers_explored.begin() + top, left, right);

#ifdef debug
        cerr << "log10prob for here: " << p_here << endl;
#endif
        
        // Calculate prob of all intervals up to top being disrupted
        double p = c[bottom] + p_here;
        
#ifdef debug
        cerr << "log10prob overall: " << p << endl;
#endif

        for (size_t i = bottom + 1; i < top + 1; i++) {
            // Replace min-prob for minimizers in the interval
            if (c[i] < p) {
#ifdef debug
                cerr << "\tBeats " << c[i] << " at rank " << i-1 << endl;
#endif
                c[i] = p;
            } else {
#ifdef debug
                cerr << "\tBeaten by " << c[i] << " at rank " << i-1 << endl;
#endif
            }
        }
    });
    
    assert(!isinf(c.back()));
    // Conver to Phred.
    return -c.back() * 10;
}

void MinimizerMapper::for_each_aglomeration_interval(const vector<Minimizer>& minimizers,
    const string& sequence, const string& quality_bytes,
    const vector<size_t>& minimizer_indices,
    const function<void(size_t, size_t, size_t, size_t)>& iteratee) {
    
    if (minimizer_indices.empty()) {
        // Handle no item case
        return;
    }

    // Items currently being iterated over
    list<const Minimizer*> stack = {&minimizers[minimizer_indices.front()]};
    // The left end of an item interval
    size_t left = stack.front()->agglomeration_start;
    // The index of the first item in the interval in the sequence of selected items
    size_t bottom = 0;

    // Emit all intervals that precede a given point "right"
    auto emit_preceding_intervals = [&](size_t right) {
        while (left < right) {
            // Work out the end position of the top thing on the stack
            size_t stack_top_end = stack.front()->agglomeration_start + stack.front()->agglomeration_length;
            if (stack_top_end <= right) {
                // Case where the left-most item ends before the start of the new item
                iteratee(left, stack_top_end, bottom, bottom + stack.size());

                // If the stack contains only one item there is a gap between the item
                // and the new item, otherwise just shift to the end of the leftmost item
                left = stack.size() == 1 ? right : stack_top_end;

                bottom += 1;
                stack.pop_front();
            } else {
                // Case where the left-most item ends at or after the beginning of the new new item
                iteratee(left, right, bottom, bottom + stack.size());
                left = right;
            }
        }
    };

    for (auto it = minimizer_indices.begin() + 1; it != minimizer_indices.end(); ++it) {
        // For each item in turn
        auto& item = minimizers[*it];
        
        assert(stack.size() > 0);

        // For each new item we return all intervals that
        // precede its start
        emit_preceding_intervals(item.agglomeration_start);

        // Add the new item for the next loop
        stack.push_back(&item);
    }

    // Intervals of the remaining intervals on the stack
    emit_preceding_intervals(sequence.size());
}

double MinimizerMapper::get_log10_prob_of_disruption_in_interval(const vector<Minimizer>& minimizers,
    const string& sequence, const string& quality_bytes,
    const vector<size_t>::iterator& disrupt_begin, const vector<size_t>::iterator& disrupt_end,
    size_t left, size_t right) {
    
    // We don't allow zero length intervals
    assert(left < right);
    
    // Start with the first column
    double p = get_log10_prob_of_disruption_in_column(minimizers, sequence, quality_bytes, disrupt_begin, disrupt_end, left);
    for(size_t i = left + 1 ; i < right; i++) {
        // OR up probability of all the other columns
        p = add_log10(get_log10_prob_of_disruption_in_column(minimizers, sequence, quality_bytes, disrupt_begin, disrupt_end, i), p);
    }
    
    // Return the result
    return p;
 
}

double MinimizerMapper::get_log10_prob_of_disruption_in_column(const vector<Minimizer>& minimizers,
    const string& sequence, const string& quality_bytes,
    const vector<size_t>::iterator& disrupt_begin, const vector<size_t>::iterator& disrupt_end,
    size_t index) {
    
    // Base cost is quality.
    double p = -(uint8_t)quality_bytes[index] / 10;
    for (auto it = disrupt_begin; it != disrupt_end; ++it) {
        // For each minimizer to disrupt
        auto m = minimizers[*it];
        
        if (!(m.forward_offset() <= index && index < m.forward_offset() + m.length)) {
            // Index is out of range of the minimizer
        }
        
        // How many new possible minimizers would an error here create in this agglomeration,
        // to compete with its minimizer?
        // No more than one per position in a minimizer sequence.
        // No more than 1 per base from the start of the agglomeration to here, inclusive.
        // No more than 1 per base from here to the last base of the agglomeration, inclusive.
        size_t possible_minimizers = min((size_t) m.length,
                                         min(index - m.agglomeration_start + 1,
                                         (m.agglomeration_start + m.agglomeration_length) - index));

        // Account for at least one of them beating the minimizer.
        p += phred_for_at_least_one(m.value.hash, possible_minimizers);
        
        // TODO: handle N somehow??? It can occur outside the minimizer itself, here in the flank.
    }
    
    return p;
}

//-----------------------------------------------------------------------------

void MinimizerMapper::attempt_rescue(const Alignment& aligned_read, Alignment& rescued_alignment, const std::vector<Minimizer>& minimizers, bool rescue_forward ) {

    if (this->rescue_algorithm == rescue_none) { return; }

    // We are traversing the same small subgraph repeatedly, so it's better to use a cache.
    gbwtgraph::CachedGBWTGraph cached_graph(this->gbwt_graph);

    // Find all nodes within a reasonable range from aligned_read.
    std::unordered_set<id_t> rescue_nodes;
    int64_t min_distance = max(0.0, fragment_length_distr.mean() - rescued_alignment.sequence().size() - rescue_subgraph_stdevs * fragment_length_distr.std_dev());
    int64_t max_distance = fragment_length_distr.mean() + rescue_subgraph_stdevs * fragment_length_distr.std_dev();
    distance_index.subgraph_in_range(aligned_read.path(), &cached_graph, min_distance, max_distance, rescue_nodes, rescue_forward);

    // Remove node ids that do not exist in the GBWTGraph from the subgraph.
    // We may be using the distance index of the original graph, and nodes
    // not visited by any thread are missing from the GBWTGraph.
    for (auto iter = rescue_nodes.begin(); iter != rescue_nodes.end(); ) {
        if (!cached_graph.has_node(*iter)) {
            iter = rescue_nodes.erase(iter);
        } else {
            ++iter;
        }
    }

    // Get rid of the old path.
    rescued_alignment.clear_path();

    // Find all seeds in the subgraph and try to get a full-length extension.
    GaplessExtender::cluster_type seeds = this->seeds_in_subgraph(minimizers, rescue_nodes);
    std::vector<GaplessExtension> extensions = this->extender.extend(seeds, rescued_alignment.sequence(), &cached_graph);
    size_t best = extensions.size();
    for (size_t i = 0; i < extensions.size(); i++) {
        if (best >= extensions.size() || extensions[i].score > extensions[best].score) {
            best = i;
        }
    }

    // If we have a full-length extension, use it as the rescued alignment.
    if (best < extensions.size() && extensions[best].full()) {
        this->extension_to_alignment(extensions[best], rescued_alignment);
        return;
    }

    // The haplotype-based algorithm is a special case.
    if (this->rescue_algorithm == rescue_haplotypes) {
        // Find and unfold the local haplotypes in the subgraph.
        std::vector<std::vector<handle_t>> haplotype_paths;
        bdsg::HashGraph align_graph;
        this->extender.unfold_haplotypes(rescue_nodes, haplotype_paths, align_graph);

        // Align to the subgraph.
        this->get_regular_aligner()->align_xdrop(rescued_alignment, align_graph,
                                                 std::vector<MaximalExactMatch>(), false);
        this->fix_dozeu_score(rescued_alignment, align_graph, std::vector<handle_t>());

        // Get the corresponding alignment to the original graph.
        this->extender.transform_alignment(rescued_alignment, haplotype_paths);
        return;
    }

    // Use the best extension as a seed for dozeu.
    // Also ensure that the entire extension is in the subgraph.
    std::vector<MaximalExactMatch> dozeu_seed;
    if (best < extensions.size()) {
        const GaplessExtension& extension = extensions[best];
        for (handle_t handle : extension.path) {
            rescue_nodes.insert(cached_graph.get_id(handle));
        }
        dozeu_seed.emplace_back();
        dozeu_seed.back().begin = rescued_alignment.sequence().begin() + extension.read_interval.first;
        dozeu_seed.back().end = rescued_alignment.sequence().begin() + extension.read_interval.second;
        nid_t id = cached_graph.get_id(extension.path.front());
        bool is_reverse = cached_graph.get_is_reverse(extension.path.front());
        gcsa::node_type node = gcsa::Node::encode(id, extension.offset, is_reverse);
        dozeu_seed.back().nodes.push_back(node);
    }

    // GSSW and dozeu assume that the graph is a DAG.
    std::vector<handle_t> topological_order = gbwtgraph::topological_order(cached_graph, rescue_nodes);
    if (!topological_order.empty()) {
        if (rescue_algorithm == rescue_dozeu) {
            get_regular_aligner()->align_xdrop(rescued_alignment, cached_graph, topological_order,
                                               dozeu_seed, false);
            this->fix_dozeu_score(rescued_alignment, cached_graph, topological_order);
        } else {
            get_regular_aligner()->align(rescued_alignment, cached_graph, topological_order);
        }
        return;
    }

    // Build a subgraph overlay.
    SubHandleGraph sub_graph(&cached_graph);
    for (id_t id : rescue_nodes) {
        sub_graph.add_handle(cached_graph.get_handle(id));
    }

    // Create an overlay where each strand is a separate node.
    StrandSplitGraph split_graph(&sub_graph);

    // Dagify the subgraph.
    bdsg::HashGraph dagified;
    std::unordered_map<id_t, id_t> dagify_trans =
        algorithms::dagify(&split_graph, &dagified, rescued_alignment.sequence().size());

    // Align to the subgraph.
    // TODO: Map the seed to the dagified subgraph.
    if (this->rescue_algorithm == rescue_dozeu) {
        get_regular_aligner()->align_xdrop(rescued_alignment, dagified, std::vector<MaximalExactMatch>(), false);
        this->fix_dozeu_score(rescued_alignment, dagified, std::vector<handle_t>());
    } else if (this->rescue_algorithm == rescue_gssw) {
        get_regular_aligner()->align(rescued_alignment, dagified, true);
    }

    // Map the alignment back to the original graph.
    Path& path = *(rescued_alignment.mutable_path());
    for (size_t i = 0; i < path.mapping_size(); i++) {
        Position& pos = *(path.mutable_mapping(i)->mutable_position());
        id_t id = dagify_trans[pos.node_id()];
        handle_t handle = split_graph.get_underlying_handle(split_graph.get_handle(id));
        pos.set_node_id(sub_graph.get_id(handle));
        pos.set_is_reverse(sub_graph.get_is_reverse(handle));
    }
}

GaplessExtender::cluster_type MinimizerMapper::seeds_in_subgraph(const std::vector<Minimizer>& minimizers,
                                                                 const std::unordered_set<id_t>& subgraph) const {
    std::vector<id_t> sorted_ids(subgraph.begin(), subgraph.end());
    std::sort(sorted_ids.begin(), sorted_ids.end());
    GaplessExtender::cluster_type result;
    for (const Minimizer& minimizer : minimizers) {
        gbwtgraph::hits_in_subgraph(minimizer.hits, minimizer.occs, sorted_ids, [&](pos_t pos, gbwtgraph::payload_type) {
            if (minimizer.value.is_reverse) {
                size_t node_length = this->gbwt_graph.get_length(this->gbwt_graph.get_handle(id(pos)));
                pos = reverse_base_pos(pos, node_length);
            }
            result.insert(GaplessExtender::to_seed(pos, minimizer.value.offset));
        });
    }
    return result;
}

void MinimizerMapper::fix_dozeu_score(Alignment& rescued_alignment, const HandleGraph& rescue_graph,
                                      const std::vector<handle_t>& topological_order) const {

    const Aligner* aligner = this->get_regular_aligner();
    int32_t score = aligner->score_ungapped_alignment(rescued_alignment);
    if (score > 0) {
        rescued_alignment.set_score(score);
    } else {
        rescued_alignment.clear_path();
        if (topological_order.empty()) {
            aligner->align(rescued_alignment, rescue_graph, true);
        } else {
            aligner->align(rescued_alignment, rescue_graph, topological_order);
        }
    }
}

//-----------------------------------------------------------------------------

int64_t MinimizerMapper::distance_between(const Alignment& aln1, const Alignment& aln2) {
    assert(aln1.path().mapping_size() != 0); 
    assert(aln2.path().mapping_size() != 0); 
     
    pos_t pos1 = initial_position(aln1.path()); 
    pos_t pos2 = final_position(aln2.path());

    int64_t min_dist = distance_index.min_distance(pos1, pos2);
    return min_dist == -1 ? numeric_limits<int64_t>::max() : min_dist;
}

void MinimizerMapper::extension_to_alignment(const GaplessExtension& extension, Alignment& alignment) const {
    *(alignment.mutable_path()) = extension.to_path(this->gbwt_graph, alignment.sequence());
    alignment.set_score(extension.score);
    double identity = 0.0;
    if (!alignment.sequence().empty()) {
        size_t len = alignment.sequence().length();
        identity = (len - extension.mismatches()) / static_cast<double>(len);
    }
    alignment.set_identity(identity);
}

//-----------------------------------------------------------------------------

std::vector<MinimizerMapper::Minimizer> MinimizerMapper::find_minimizers(const std::string& sequence, Funnel& funnel) const {

    if (this->track_provenance) {
        // Start the minimizer finding stage
        funnel.stage("minimizer");
    }

    std::vector<Minimizer> result;
    double base_score = 1.0 + std::log(this->hard_hit_cap);
    for (size_t i = 0; i < this->minimizer_indexes.size(); i++) {
        // Get minimizers and their window agglomeration starts and lengths
        vector<tuple<gbwtgraph::DefaultMinimizerIndex::minimizer_type, size_t, size_t>> current_minimizers = 
            minimizer_indexes[i]->minimizer_regions(sequence);
        for (auto& m : current_minimizers) {
            double score = 0.0;
            auto hits = this->minimizer_indexes[i]->count_and_find(get<0>(m));
            if (hits.first > 0) {
                if (hits.first <= this->hard_hit_cap) {
                    score = base_score - std::log(hits.first);
                } else {
                    score = 1.0;
                }
            }
            result.push_back({ std::get<0>(m), std::get<1>(m), std::get<2>(m), hits.first, hits.second,
                               (int32_t) minimizer_indexes[i]->k(), (int32_t) minimizer_indexes[i]->w(), score });
        }
    }
    std::sort(result.begin(), result.end());

    if (this->track_provenance) {
        // Record how many we found, as new lines.
        funnel.introduce(result.size());
    }

    return result;
}

std::vector<MinimizerMapper::Seed> MinimizerMapper::find_seeds(const std::vector<Minimizer>& minimizers, const Alignment& aln, Funnel& funnel) const {

    if (this->track_provenance) {
        // Start the minimizer locating stage
        funnel.stage("seed");
    }

    // One of the filters accepts minimizers until selected_score reaches target_score.
    double base_target_score = 0.0;
    for (const Minimizer& minimizer : minimizers) {
        base_target_score += minimizer.score;
    }
    double target_score = (base_target_score * this->minimizer_score_fraction) + 0.000001;
    double selected_score = 0.0;

    // In order to consistently take either all or none of the minimizers in
    // the read with a particular sequence, we track whether we took the
    // previous one.
    bool took_last = false;

    // Select the minimizers we use for seeds.
    size_t rejected_count = 0;
    std::vector<Seed> seeds;
    // Flag whether each minimizer in the read was located or not, for MAPQ capping.
    // We ignore minimizers with no hits (count them as not located), because
    // they would have to be created in the read no matter where we say it came
    // from, and because adding more of them should lower the MAPQ cap, whereas
    // locating more of the minimizers that are present and letting them pass
    // to the enxt stage should raise the cap.
    for (size_t i = 0; i < minimizers.size(); i++) {
        if (this->track_provenance) {
            // Say we're working on it
            funnel.processing_input(i);
        }

        // Select the minimizer if it is informative enough or if the total score
        // of the selected minimizers is not high enough.
        const Minimizer& minimizer = minimizers[i];

#ifdef debug
        std::cerr << "Minimizer " << i << " = " << minimizer.value.key.decode(minimizer.length)
             << " has " << minimizer.hits << " hits" << std::endl;
#endif

        if (minimizer.hits == 0) {
            // A minimizer with no hits can't go on.
            took_last = false;
            // We do not treat it as located for MAPQ capping purposes.
            if (this->track_provenance) {
                funnel.fail("any-hits", i);
            }
        } else if (minimizer.hits <= this->hit_cap ||
            (minimizer.hits <= this->hard_hit_cap && selected_score + minimizer.score <= target_score) ||
            (took_last && i > 0 && minimizer.value.key == minimizers[i - 1].value.key)) {
            
            // We should keep this minimizer instance because it is
            // sufficiently rare, or we want it to make target_score, or it is
            // the same sequence as the previous minimizer which we also took.

            // Locate the hits.
            for (size_t j = 0; j < minimizer.hits; j++) {
                pos_t hit = gbwtgraph::Position::decode(minimizer.occs[j].pos);
                // Reverse the hits for a reverse minimizer
                if (minimizer.value.is_reverse) {
                    size_t node_length = this->gbwt_graph.get_length(this->gbwt_graph.get_handle(id(hit)));
                    hit = reverse_base_pos(hit, node_length);
                }
                // Extract component id and offset in the root chain, if we have them for this seed.
                // TODO: Get all the seed values here
                tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info
                    (false, MIPayload::NO_VALUE, MIPayload::NO_VALUE, false, MIPayload::NO_VALUE, MIPayload::NO_VALUE, MIPayload::NO_VALUE, MIPayload::NO_VALUE, false );
                if (minimizer.occs[j].payload != MIPayload::NO_CODE) {
                    chain_info = MIPayload::decode(minimizer.occs[j].payload);
                }
                seeds.push_back({ hit, i, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info), 
                    std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info) });
            }
            
            if (!(took_last && i > 0 && minimizer.value.key == minimizers[i - 1].value.key)) {
                // We did not also take a previous identical-sequence minimizer, so count this one towards the score.
                selected_score += minimizer.score;
            }

            // Remember that we took this minimizer
            took_last = true;

            if (this->track_provenance) {
                // Record in the funnel that this minimizer gave rise to these seeds.
                funnel.pass("any-hits", i);
                funnel.pass("hard-hit-cap", i);
                funnel.pass("hit-cap||score-fraction", i, selected_score  / base_target_score);
                funnel.expand(i, minimizer.hits);
            }
        } else if (minimizer.hits <= this->hard_hit_cap) {
            // Passed hard hit cap but failed score fraction/normal hit cap
            took_last = false;
            rejected_count++;
            if (this->track_provenance) {
                funnel.pass("any-hits", i);
                funnel.pass("hard-hit-cap", i);
                funnel.fail("hit-cap||score-fraction", i, (selected_score + minimizer.score) / base_target_score);
            }
            //Stop looking for more minimizers once we fail the score fraction
            target_score = selected_score; 
        } else {
            // Failed hard hit cap
            took_last = false;
            rejected_count++;
            if (this->track_provenance) {
                funnel.pass("any-hits", i);
                funnel.fail("hard-hit-cap", i);
            }
        }
        if (this->track_provenance) {
            // Say we're done with this input item
            funnel.processed_input();
        }
    }

    if (this->track_provenance && this->track_correctness) {
        // Tag seeds with correctness based on proximity along paths to the input read's refpos
        funnel.substage("correct");

        if (this->path_graph == nullptr) {
            cerr << "error[vg::MinimizerMapper] Cannot use track_correctness with no XG index" << endl;
            exit(1);
        }

        if (aln.refpos_size() != 0) {
            // Take the first refpos as the true position.
            auto& true_pos = aln.refpos(0);

            for (size_t i = 0; i < seeds.size(); i++) {
                // Find every seed's reference positions. This maps from path name to pairs of offset and orientation.
                auto offsets = algorithms::nearest_offsets_in_paths(this->path_graph, seeds[i].pos, 100);
                for (auto& hit_pos : offsets[this->path_graph->get_path_handle(true_pos.name())]) {
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
    std::cerr << "Found " << seeds.size() << " seeds from " << (minimizers.size() - rejected_count) << " minimizers, rejected " << rejected_count << std::endl;
#endif

    return seeds;
}

//-----------------------------------------------------------------------------

void MinimizerMapper::score_cluster(Cluster& cluster, size_t i, const std::vector<Minimizer>& minimizers, const std::vector<Seed>& seeds, size_t seq_length, Funnel& funnel) const {

    if (this->track_provenance) {
        // Say we're making it
        funnel.producing_output(i);
    }

    // Initialize the values.
    cluster.score = 0.0;
    cluster.coverage = 0.0;
    cluster.present = SmallBitset(minimizers.size());

    // Determine the minimizers that are present in the cluster.
    for (auto hit_index : cluster.seeds) {
        cluster.present.insert(seeds[hit_index].source);
#ifdef debug
        cerr << "Minimizer " << seeds[hit_index].source << " is present in cluster " << i << endl;
#endif
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
        // Record the cluster in the funnel as a group of the size of the number of items.
        funnel.merge_group(cluster.seeds.begin(), cluster.seeds.end());
        funnel.score(funnel.latest(), cluster.score);

        // Say we made it.
        funnel.produced_output();
    }
}

//-----------------------------------------------------------------------------

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

std::vector<int> MinimizerMapper::score_extensions(const std::vector<std::vector<GaplessExtension>>& extensions, const Alignment& aln, Funnel& funnel) const {

    // Extension scoring substage.
    if (this->track_provenance) {
        funnel.substage("score");
    }

    // We now estimate the best possible alignment score for each cluster.
    std::vector<int> result(extensions.size(), 0);
    for (size_t i = 0; i < extensions.size(); i++) {
        
        if (this->track_provenance) {
            funnel.producing_output(i);
        }
        
        result[i] = score_extension_group(aln, extensions[i], get_regular_aligner()->gap_open, get_regular_aligner()->gap_extension);
        
        // Record the score with the funnel.
        if (this->track_provenance) {
            funnel.score(i, result[i]);
            funnel.produced_output();
        }
    }

    return result;
}

std::vector<int> MinimizerMapper::score_extensions(const std::vector<std::pair<std::vector<GaplessExtension>, size_t>>& extensions, const Alignment& aln, Funnel& funnel) const {

    // Extension scoring substage.
    if (this->track_provenance) {
        funnel.substage("score");
    }

    // We now estimate the best possible alignment score for each cluster.
    std::vector<int> result(extensions.size(), 0);
    for (size_t i = 0; i < extensions.size(); i++) {
        
        if (this->track_provenance) {
            funnel.producing_output(i);
        }
        
        result[i] = score_extension_group(aln, extensions[i].first, get_regular_aligner()->gap_open, get_regular_aligner()->gap_extension);
        
        // Record the score with the funnel.
        if (this->track_provenance) {
            funnel.score(i, result[i]);
            funnel.produced_output();
        }
    }

    return result;
}

//-----------------------------------------------------------------------------

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
    process_until_threshold_b(extended_seeds, extension_path_scores,
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
            
            // X-drop align, accounting for full length bonus.
            // We *always* do left-pinned alignment internally, since that's the shape of trees we get.
            get_regular_aligner()->align_pinned(current_alignment, subgraph, true, true);
            
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


