/**
 * \file minimizer_mapper.cpp
 * Defines the code for the minimizer-and-GBWT-based mapper.
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

namespace vg {

using namespace std;

MinimizerMapper::MinimizerMapper(const gbwtgraph::GBWTGraph& graph,
    const gbwtgraph::DefaultMinimizerIndex& minimizer_index,
    MinimumDistanceIndex* old_distance_index, SnarlDistanceIndex* distance_index, 
    const PathPositionHandleGraph* path_graph) :
    path_graph(path_graph), minimizer_index(minimizer_index),
    distance_index(distance_index), old_distance_index(old_distance_index), 
    clusterer(distance_index, &graph), old_clusterer(old_distance_index),
    gbwt_graph(graph),
    extender(gbwt_graph, *(get_regular_aligner())),
    fragment_length_distr(1000,1000,0.95) {
    
    // The GBWTGraph needs a GBWT
    assert(graph.index != nullptr);
}

//-----------------------------------------------------------------------------

string MinimizerMapper::log_name() {
    return "T" + to_string(omp_get_thread_num()) + ":\t";
}

void MinimizerMapper::dump_debug_sequence(ostream& out, const string& sequence) {
    int digits_needed = (int) ceil(log10(sequence.size()));
    for (int digit = digits_needed - 1; digit >= 0; digit--) {
        out << log_name();
        for (size_t i = 0; i < sequence.size(); i++) {
            // Output the correct digit for this place in this number
            out << (char) ('0' + (uint8_t) floor(i % (int) round(pow(10, digit + 1)) / pow(10, digit)));
        }
        out << endl;
    }
    out << log_name() << sequence << endl;
}

void MinimizerMapper::dump_debug_extension_set(const HandleGraph& graph, const Alignment& aln, const vector<GaplessExtension>& extended_seeds) {
    dump_debug_sequence(cerr, aln.sequence());
    
    for (auto& ext : extended_seeds) {
        // For each extension
        
        cerr << log_name();
        
        for (size_t i = 0; i < ext.read_interval.first; i++) {
            // Space until it starts
            cerr << ' ';
        }
        
        for (size_t i = ext.read_interval.first; i < ext.read_interval.second; i++) {
            if (std::find(ext.mismatch_positions.begin(), ext.mismatch_positions.end(), i) != ext.mismatch_positions.end()) {
                // Has an error here
                cerr << "*";
            } else {
                // A match
                cerr << aln.sequence()[i];
            }
        }
        cerr << " @";
        for (const handle_t& h : ext.path) {
            cerr << " " << graph.get_id(h);
        }
        cerr << endl;
    }
}

void MinimizerMapper::dump_debug_minimizers(const vector<MinimizerMapper::Minimizer>& minimizers, const string& sequence, const vector<size_t>* to_include) {

    dump_debug_sequence(cerr, sequence);
    
    vector<size_t> all;
    if (to_include == nullptr) {
        // Synthesize a list of all minimizers
        to_include = &all;
        for (size_t i = 0; i < minimizers.size(); i++) {
            all.push_back(i);
        }
        
        // Sort minimizer subset so we go through minimizers in increasing order of start position
        std::sort(all.begin(), all.end(), [&](size_t a, size_t b) {
            // Return true if a must come before b, and false otherwise
            return minimizers[a].forward_offset() < minimizers[b].forward_offset();
        });
    }
    
    // Dump minimizers
    for (auto& index : *to_include) {
        // For each minimizer
        
        cerr << log_name();
        
        auto& m = minimizers[index];
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
        
        // Tag with metadata
        cerr << " (#" << index << ", " << m.hits << " hits)" << endl;
    }
}

//-----------------------------------------------------------------------------

void MinimizerMapper::map(Alignment& aln, AlignmentEmitter& alignment_emitter) {
    // Ship out all the aligned alignments
    alignment_emitter.emit_mapped_single(map(aln));
}

vector<Alignment> MinimizerMapper::map(Alignment& aln) {
    
    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "Read " << aln.name() << ": " << aln.sequence() << endl;

        }
    }
    
    // Make a new funnel instrumenter to watch us map this read.
    Funnel funnel;
    funnel.start(aln.name());
    
    // Prepare the RNG for shuffling ties, if needed
    LazyRNG rng([&]() {
        return aln.sequence();
    });


    // Minimizers sorted by score in descending order.
    std::vector<Minimizer> minimizers = this->find_minimizers(aln.sequence(), funnel);

    //Since there can be two different versions of a distance index, find seeds and clusters differently

    //One of these two will be filled
    std::vector<Seed> seeds; std::vector<OldSeed> old_seeds;
    std::vector<Cluster> clusters;
    if (distance_index != nullptr) {
        //If there is the new version of the distance index

        // Find the seeds and mark the minimizers that were located.
        seeds = this->find_seeds(minimizers, aln, funnel);

        // Cluster the seeds. Get sets of input seed indexes that go together.
        if (track_provenance) {
            funnel.stage("cluster");
        }

        clusters = clusterer.cluster_seeds(seeds, get_distance_limit(aln.sequence().size()));
    } else {
        //Old version of the distance index
        assert(old_distance_index != nullptr);

        // Find the seeds and mark the minimizers that were located.
        old_seeds = this->find_seeds_old(minimizers, aln, funnel);

        // Cluster the seeds. Get sets of input seed indexes that go together.
        if (track_provenance) {
            funnel.stage("cluster");
        }

        clusters = old_clusterer.cluster_seeds(old_seeds, get_distance_limit(aln.sequence().size()));
    }

    // Determine the scores and read coverages for each cluster.
    // Also find the best and second-best cluster scores.
    if (this->track_provenance) {
        funnel.substage("score");
    }
    double best_cluster_score = 0.0, second_best_cluster_score = 0.0;
    for (size_t i = 0; i < clusters.size(); i++) {
        Cluster& cluster = clusters[i];
        if (distance_index != nullptr) {
            this->score_cluster(cluster, i, minimizers, seeds, aln.sequence().length(), funnel);
        } else {
            this->score_cluster_old(cluster, i, minimizers, old_seeds, aln.sequence().length(), funnel);
        }
        if (cluster.score > best_cluster_score) {
            second_best_cluster_score = best_cluster_score;
            best_cluster_score = cluster.score;
        } else if (cluster.score > second_best_cluster_score) {
            second_best_cluster_score = cluster.score;
        }
    }

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

    size_t kept_cluster_count = 0;
    
    //Process clusters sorted by both score and read coverage
    process_until_threshold_c<Cluster, double>(clusters, [&](size_t i) -> double {
            return clusters[i].coverage;
        }, [&](size_t a, size_t b) -> bool {
            return ((clusters[a].coverage > clusters[b].coverage) ||
                    (clusters[a].coverage == clusters[b].coverage && clusters[a].score > clusters[b].score));
        }, cluster_coverage_threshold, min_extensions, max_extensions, rng, [&](size_t cluster_num) {
            // Handle sufficiently good clusters in descending coverage order
            
            Cluster& cluster = clusters[cluster_num];
            if (track_provenance) {
                funnel.pass("cluster-coverage", cluster_num, cluster.coverage);
                funnel.pass("max-extensions", cluster_num);
            }
            
            // First check against the additional score filter
            if (cluster_score_threshold != 0 && cluster.score < cluster_score_cutoff 
                && kept_cluster_count >= min_extensions) {
                //If the score isn't good enough and we already kept at least min_extensions clusters,
                //ignore this cluster
                if (track_provenance) {
                    funnel.fail("cluster-score", cluster_num, cluster.score);
                }
                if (show_work) {
                    #pragma omp critical (cerr)
                    {
                        cerr << log_name() << "Cluster " << cluster_num << " fails cluster score cutoff" <<  endl;
                        cerr << log_name() << "Covers " << clusters[cluster_num].coverage << "/best-" << cluster_coverage_threshold << " of read" << endl;
                        cerr << log_name() << "Scores " << clusters[cluster_num].score << "/" << cluster_score_cutoff << endl;
                    }
                }
                return false;
            }
            
            if (track_provenance) {
                funnel.pass("cluster-score", cluster_num, cluster.score);
                funnel.processing_input(cluster_num);
            }
            

            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Cluster " << cluster_num << endl;
                    cerr << log_name() << "Covers " << cluster.coverage << "/best-" << cluster_coverage_threshold << " of read" << endl;
                    cerr << log_name() << "Scores " << cluster.score << "/" << cluster_score_cutoff << endl;
                }
            }
             
            minimizer_extended_cluster_count.emplace_back(minimizers.size(), 0);
            // Pack the seeds for GaplessExtender.
            GaplessExtender::cluster_type seed_matchings;
            if (distance_index != nullptr) {
                for (auto seed_index : cluster.seeds) {
                    // Insert the (graph position, read offset) pair.
                    const Seed& seed = seeds[seed_index];
                    seed_matchings.insert(GaplessExtender::to_seed(seed.pos, minimizers[seed.source].value.offset));
                    minimizer_extended_cluster_count.back()[seed.source]++;
                    
                    if (show_work) {
                        #pragma omp critical (cerr)
                        {
                            const Minimizer& minimizer = minimizers[seed.source];
                            cerr << log_name() << "Seed read:" << minimizer.value.offset << " = " << seed.pos
                                << " from minimizer " << seed.source << "(" << minimizer.hits << "), #" << seed_index << endl;
                        }
                    }
                }
            } else {
                for (auto seed_index : cluster.seeds) {
                    // Insert the (graph position, read offset) pair.
                    const OldSeed& seed = old_seeds[seed_index];
                    seed_matchings.insert(GaplessExtender::to_seed(seed.pos, minimizers[seed.source].value.offset));
                    minimizer_extended_cluster_count.back()[seed.source]++;
                    
                    if (show_work) {
                        #pragma omp critical (cerr)
                        {
                            const Minimizer& minimizer = minimizers[seed.source];
                            cerr << log_name() << "Seed read:" << minimizer.value.offset << " = " << seed.pos
                                << " from minimizer " << seed.source << "(" << minimizer.hits << "), #" << seed_index << endl;
                        }
                    }
                }
            }
            
            // Extend seed hits in the cluster into one or more gapless extensions
            cluster_extensions.emplace_back(std::move(extender.extend(seed_matchings, aln.sequence())));

            kept_cluster_count ++;
            
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Extensions:" << endl;
                    for (auto& e : cluster_extensions.back()) {
                        cerr << log_name() << "\tRead " << e.read_interval.first << "-" << e.read_interval.second << " with " << e.mismatch_positions.size() << " mismatches:";
                        for (auto& pos : e.mismatch_positions) {
                            cerr << " " << pos;
                        }
                        cerr << endl;
                    }
                }
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
            Cluster& cluster = clusters[cluster_num];
            if (track_provenance) {
                funnel.pass("cluster-coverage", cluster_num, cluster.coverage);
                funnel.fail("max-extensions", cluster_num);
            }
            
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    
                    cerr << log_name() << "Cluster " << cluster_num << " passes cluster cutoffs but we have too many" <<  endl;
                    cerr << log_name() << "Covers " << cluster.coverage << "/best-" << cluster_coverage_threshold << " of read" << endl;
                    cerr << log_name() << "Scores " << cluster.score << "/" << cluster_score_cutoff << endl;
                }
            }
            
        }, [&](size_t cluster_num) {
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
    std::vector<int> cluster_extension_scores;
    // We optionally also get tracebacks for aligining each cluster into a base-level alignment
    std::vector<std::vector<size_t>> cluster_extension_tracebacks;
    if (align_from_chains) {
        // Chain all the extension groups and get scores and chains
        auto result = this->chain_extensions(cluster_extensions, aln, funnel);
        cluster_extension_scores = std::move(result.first);
        cluster_extension_tracebacks = std::move(result.second);
    } else {
        // Just score the extension groups, with a slithgly simpler algorithm; don't chain them
        cluster_extension_scores = this->score_extensions(cluster_extensions, aln, funnel);
    }
    
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
        extension_set_score_threshold, 2, max_alignments, rng, [&](size_t extension_num) {
            // This extension set is good enough.
            // Called in descending score order.
            
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "gapless extension group " << extension_num << " is good enough (score=" << cluster_extension_scores[extension_num] << ")" << endl;
                    if (track_correctness && funnel.was_correct(extension_num)) {
                        cerr << log_name() << "\tCORRECT!" << endl;
                        dump_debug_extension_set(gbwt_graph, aln, cluster_extensions[extension_num]);
                    }
                }
            }
            if (track_provenance) {
                funnel.pass("extension-set", extension_num, cluster_extension_scores[extension_num]);
                funnel.pass("max-alignments", extension_num);
                funnel.processing_input(extension_num);
            }

            
            auto& extensions = cluster_extensions[extension_num];
            
            // Collect the top alignments. Make sure we have at least one always, starting with unaligned.
            vector<Alignment> best_alignments(1, aln);

            if (GaplessExtender::full_length_extensions(extensions)) {
                // We got full-length extensions, so directly convert to an Alignment.
                
                if (track_provenance) {
                    funnel.substage("direct");
                }
                
                //Fill in the best alignments from the extension. We know the top one is always full length and exists.
                this->extension_to_alignment(extensions.front(), best_alignments.front());
                
                if (show_work) {
                    #pragma omp critical (cerr)
                    {
                        cerr << log_name() << "Produced alignment directly from full length gapless extension " << extension_num << endl;
                    }
                }
                
                for (auto next_ext_it = extensions.begin() + 1; next_ext_it != extensions.end() && next_ext_it->full(); ++next_ext_it) {
                    // For all subsequent full length extensions, make them into alignments too.
                    // We want them all to go on to the pairing stage so we don't miss a possible pairing in a tandem repeat.
                    best_alignments.emplace_back(aln);
                    this->extension_to_alignment(*next_ext_it, best_alignments.back());
                    
                    if (show_work) {
                        #pragma omp critical (cerr)
                        {
                            cerr << log_name() << "Produced additional alignment directly from full length gapless extension " << (next_ext_it - extensions.begin()) << endl;
                        }
                    }
                    
                }
                
                if (track_provenance) {
                    // Stop the current substage
                    funnel.substage_stop();
                }
            } else if (do_dp) {
                // We need to do base-level alignment.
                
                if (track_provenance) {
                    funnel.substage("align");
                }
                
                if (align_from_chains) {
                    // Do the DP between the extensions in the cluster as specified by the chain we got for it.
                    auto& chain_order = cluster_extension_tracebacks[extension_num];
                    
                    // TODO: implement!
                    throw std::runtime_error("Not implemented!");
                    
                } else {
                    // Do the DP and compute up to 2 alignments from the individual gapless extensions
                    best_alignments.emplace_back(aln);
                    find_optimal_tail_alignments(aln, extensions, rng, best_alignments[0], best_alignments[1]);
                }

                if (show_work) {
                    #pragma omp critical (cerr)
                    {
                        cerr << log_name() << "Did dynamic programming for gapless extension group " << extension_num << endl;
                    }
                }
                
                if (track_provenance) {
                    // We're done base-level alignment. Next alignment may not go through this substage.
                    funnel.substage_stop();
                }
            } else {
                // We would do base-level alignment but it is disabled.
                // Leave best_alignment unaligned
            }
           
            // Have a function to process the best alignments we obtained
            auto observe_alignment = [&](Alignment& aln) {
                alignments.emplace_back(std::move(aln));
                alignments_to_source.push_back(extension_num);

                if (track_provenance) {
    
                    funnel.project(extension_num);
                    funnel.score(alignments.size() - 1, alignments.back().score());
                }
                if (show_work) {
                    #pragma omp critical (cerr)
                    {
                        cerr << log_name() << "Produced alignment from gapless extension group " << extension_num
                            << " with score " << alignments.back().score() << ": " << pb2json(alignments.back()) << endl;
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
            
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "gapless extension group " << extension_num << " failed because there were too many good extensions (score=" << cluster_extension_scores[extension_num] << ")" << endl;
                    if (track_correctness && funnel.was_correct(extension_num)) {
                        cerr << log_name() << "\tCORRECT!" << endl;
                        dump_debug_extension_set(gbwt_graph, aln, cluster_extensions[extension_num]);
                    }
                }
            }
        }, [&](size_t extension_num) {
            // This extension is not good enough.
            if (track_provenance) {
                funnel.fail("extension-set", extension_num, cluster_extension_scores[extension_num]);
            }
            
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "gapless extension group " << extension_num << " failed because its score was not good enough (score=" << cluster_extension_scores[extension_num] << ")" << endl;
                    if (track_correctness && funnel.was_correct(extension_num)) {
                        cerr << log_name() << "\tCORRECT!" << endl;
                        dump_debug_extension_set(gbwt_graph, aln, cluster_extensions[extension_num]);
                    }
                }
            }
        });
    
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
    
    process_until_threshold_a(alignments, (std::function<double(size_t)>) [&](size_t i) -> double {
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
            cerr << log_name() << "Picked best alignment " << pb2json(mappings[0]) << endl;
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
            if (distance_index != nullptr) {
                annotate_with_minimizer_statistics(mappings[0], minimizers, seeds, funnel);
            } else {
                old_annotate_with_minimizer_statistics(mappings[0], minimizers, old_seeds, funnel);
            }
        }
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
             << minimizer_extensions_count[i];
         if (minimizer_extensions_count[i]>0) {
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

    if (track_provenance && show_work) {
        // Dump the funnel info graph.
        #pragma omp critical (cerr)
        {
            funnel.to_dot(cerr);
        }
    }

    return mappings;
}

//-----------------------------------------------------------------------------

void MinimizerMapper::pair_all(pair<vector<Alignment>, vector<Alignment>>& mappings) const {
    if (!mappings.first.empty()) {
        for (auto& next : mappings.second) {
            // Each read 2 needs to point to read 1 as its predecessor
            next.mutable_fragment_prev()->set_name(mappings.first.front().name());
        }
    }
    if (!mappings.second.empty()) {
        for (auto& prev : mappings.first) {
            // Each read 1 needs to point to read 2 as its successor
            prev.mutable_fragment_next()->set_name(mappings.second.front().name());
        }
    }
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

            if (dist == std::numeric_limits<int64_t>::max() ||
                dist >= max_fragment_length) {
                //If the distance between them is ambiguous or it it large enough that we don't think it's valid, leave them unmapped

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
            mapped_pair.first.emplace_back(std::move(alns1.front()));
            mapped_pair.second.emplace_back(std::move(alns2.front()));
            pair_all(mapped_pair);

            //TODO: This is "properly paired" if there is a path between the reads. Since we
            //haven't finalized the distribution we can't compare it
            bool properly_paired = dist != std::numeric_limits<int64_t>::max();
            set_annotation(mapped_pair.first.back(), "proper_pair", properly_paired);
            set_annotation(mapped_pair.second.back(), "proper_pair", properly_paired);


#ifdef debug_fragment_distr
            //Print stats about finalizing the fragment length distribution, copied from mpmap
            if (fragment_length_distr.is_finalized()) {
                cerr << "finalized read distribution with " << fragment_length_distr.max_sample_size() << " measurements" << endl;
                cerr << "mean: " << fragment_length_distr.mean() << endl;
                cerr << "std dev: " << fragment_length_distr.std_dev() << endl;
                cerr << "ambiguous buffer contains pairs:" << endl;
                for (pair<Alignment,Alignment>& aln_pair : ambiguous_pair_buffer) {
                    cerr << "\t" << aln_pair.first.name() << ", " << aln_pair.second.name() << endl;
                }
                cerr << "distance measurements:" << endl;
                auto iter = fragment_length_distr.measurements_begin();
                if (iter != fragment_length_distr.measurements_end()) {
                    cerr << *iter;
                    iter++;
                }
                for (; iter != fragment_length_distr.measurements_end(); iter++) {
                    cerr << ", " << *iter;
                }
                cerr << endl;
            }
#endif


            return mapped_pair;

        } else {
            // Otherwise, discard the mappings and put them in the ambiguous buffer

            ambiguous_pair_buffer.emplace_back(aln1, aln2);
            pair<vector<Alignment>, vector<Alignment>> empty;
            return empty;
        }
    }
}

pair<vector<Alignment>, vector<Alignment>> MinimizerMapper::map_paired(Alignment& aln1, Alignment& aln2) {
    
    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "Read pair " << aln1.name() << ": " << aln1.sequence() << " and " << aln2.name() << ": " << aln2.sequence() << endl;
        }
    }

    // Make sure we actually have a working fragment length distribution that the clusterer will accept.
    int64_t fragment_distance_limit = fragment_length_distr.mean() + paired_distance_stdevs * fragment_length_distr.std_dev();
    if (fragment_distance_limit < get_distance_limit(aln1.sequence().size())) {
        // We can't use this distribution
        
        if (!warned_about_bad_distribution.test_and_set()) {
            // We get to print the warning
            cerr << "warning[vg::giraffe]: Cannot cluster reads with a fragment distance smaller than read distance" << endl;
            cerr << "                      Fragment length distribution: mean=" << fragment_length_distr.mean() 
                 << ", stdev=" << fragment_length_distr.std_dev() << endl;
            cerr << "                      Fragment distance limit: " << fragment_distance_limit 
                 << ", read distance limit: " << get_distance_limit(aln1.sequence().size()) << endl;
            cerr << "warning[vg::giraffe]: Falling back on single-end mapping" << endl;
        }
        
        // Map single-ended and bail
        auto mapped_pair = make_pair(map(aln1), map(aln2));
        pair_all(mapped_pair);
        return mapped_pair;
    }


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
    
    // Prepare the RNG for shuffling ties, if needed
    LazyRNG rng([&]() {
        return aln1.sequence() + aln2.sequence();
    });
    
    // Minimizers for both reads, sorted by score in descending order.
    std::vector<std::vector<Minimizer>> minimizers_by_read(2);
    minimizers_by_read[0] = this->find_minimizers(aln1.sequence(), funnels[0]);
    minimizers_by_read[1] = this->find_minimizers(aln2.sequence(), funnels[1]);

    // Seeds for both reads, stored in separate vectors.
    std::vector<std::vector<Seed>> seeds_by_read(2);
    std::vector<std::vector<OldSeed>> old_seeds_by_read(2);
    if (distance_index != nullptr) {
        seeds_by_read[0] = this->find_seeds(minimizers_by_read[0], aln1, funnels[0]);
        seeds_by_read[1] = this->find_seeds(minimizers_by_read[1], aln2, funnels[1]);
    } else {
        old_seeds_by_read[0] = this->find_seeds_old(minimizers_by_read[0], aln1, funnels[0]);
        old_seeds_by_read[1] = this->find_seeds_old(minimizers_by_read[1], aln2, funnels[1]);
    }

    // Cluster the seeds. Get sets of input seed indexes that go together.
    if (track_provenance) {
        funnels[0].stage("cluster");
        funnels[1].stage("cluster");
    }
    std::vector<std::vector<Cluster>> all_clusters;
    if (distance_index != nullptr) {
        all_clusters = clusterer.cluster_seeds(seeds_by_read, get_distance_limit(aln1.sequence().size()), fragment_distance_limit);
    } else {
        all_clusters = old_clusterer.cluster_seeds(old_seeds_by_read, get_distance_limit(aln1.sequence().size()), fragment_distance_limit);
    }


    //Keep track of which fragment clusters (clusters of clusters) have read clusters from each end

    size_t max_fragment_num = 0;
    for (auto& cluster : all_clusters[0]) {
        max_fragment_num = std::max(max_fragment_num, cluster.fragment);
    }
    for (auto& cluster : all_clusters[1]) {
        max_fragment_num = std::max(max_fragment_num, cluster.fragment);
    }

    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "Found " << max_fragment_num + 1 << " fragment clusters" << endl;
        }
    }

    vector<bool> has_first_read (max_fragment_num+1, false);//For each fragment cluster, does it have a cluster for the first read
    vector<bool> fragment_cluster_has_pair (max_fragment_num+1, false);//Does a fragment cluster have both reads
    bool found_paired_cluster = false;
    for (auto& cluster : all_clusters[0]) {
        has_first_read[cluster.fragment] = true;
    }
    for (auto& cluster : all_clusters[1]) {
        size_t fragment_num = cluster.fragment;
        fragment_cluster_has_pair[fragment_num] = has_first_read[fragment_num];
        if (has_first_read[fragment_num]) {
            found_paired_cluster = true;
            
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Fragment cluster " << fragment_num << " has read clusters from both reads" << endl;
                }
            }
        }
    }

    if (track_provenance) {
        funnels[0].substage("score");
        funnels[1].substage("score");
    }


    //Keep track of the best cluster score and coverage per end for each fragment cluster
    pair<vector<double>, vector<double>> cluster_score_by_fragment;
    cluster_score_by_fragment.first.resize(max_fragment_num + 1, 0.0);
    cluster_score_by_fragment.second.resize(max_fragment_num + 1, 0.0);
    pair<vector<double>, vector<double>> cluster_coverage_by_fragment;
    cluster_coverage_by_fragment.first.resize(max_fragment_num + 1, 0.0);
    cluster_coverage_by_fragment.second.resize(max_fragment_num + 1, 0.0);

    //Get the scores of each cluster
    for (size_t read_num = 0 ; read_num < 2 ; read_num++) {
        Alignment& aln = read_num == 0 ? aln1 : aln2;
        std::vector<Cluster>& clusters = all_clusters[read_num];
        std::vector<Minimizer>& minimizers = minimizers_by_read[read_num];
        if (distance_index != nullptr) {
            std::vector<Seed>& seeds = seeds_by_read[read_num];
            vector<double>& best_cluster_score = read_num == 0 ? cluster_score_by_fragment.first 
                                                               : cluster_score_by_fragment.second;
            vector<double>& best_cluster_coverage = read_num == 0 ? cluster_coverage_by_fragment.first 
                                                                  : cluster_coverage_by_fragment.second;

            for (size_t i = 0; i < clusters.size(); i++) {
                // Determine cluster score and read coverage.
                Cluster& cluster = clusters[i];
                this->score_cluster(cluster, i, minimizers, seeds, aln.sequence().length(), funnels[read_num]);
                size_t fragment = cluster.fragment;
                best_cluster_score[fragment] = std::max(best_cluster_score[fragment], cluster.score);
                best_cluster_coverage[fragment] = std::max(best_cluster_coverage[fragment], cluster.coverage);
            }
        } else {
            std::vector<OldSeed>& seeds = old_seeds_by_read[read_num];
            vector<double>& best_cluster_score = read_num == 0 ? cluster_score_by_fragment.first 
                                                               : cluster_score_by_fragment.second;
            vector<double>& best_cluster_coverage = read_num == 0 ? cluster_coverage_by_fragment.first 
                                                                  : cluster_coverage_by_fragment.second;

            for (size_t i = 0; i < clusters.size(); i++) {
                // Determine cluster score and read coverage.
                Cluster& cluster = clusters[i];
                this->score_cluster_old(cluster, i, minimizers, seeds, aln.sequence().length(), funnels[read_num]);
                size_t fragment = cluster.fragment;
                best_cluster_score[fragment] = std::max(best_cluster_score[fragment], cluster.score);
                best_cluster_coverage[fragment] = std::max(best_cluster_coverage[fragment], cluster.coverage);
            }
        }
    }

    //For each fragment cluster, we want to know how many equivalent or better clusters we found
    //We consider two fragment clusters to be "equivalent" if the sum of the best score and
    //coverage for both ends is the same

    //Get a vector of the indices of fragment clusters so we can sort
    vector<size_t> fragment_cluster_indices_by_score (max_fragment_num + 1);
    for (size_t i = 0 ; i < fragment_cluster_indices_by_score.size() ; i++) {
        fragment_cluster_indices_by_score[i] = i;
    }

    //Sort by the sum of the score and coverage of the best cluster for each end
    sort_shuffling_ties(fragment_cluster_indices_by_score.begin(), fragment_cluster_indices_by_score.end(), [&](size_t a, size_t b) {
        return cluster_coverage_by_fragment.first[a] + cluster_coverage_by_fragment.second[a] 
               + cluster_score_by_fragment.first[a] + cluster_score_by_fragment.second[a]  
            > cluster_coverage_by_fragment.first[b] + cluster_coverage_by_fragment.second[b] 
                + cluster_score_by_fragment.first[b] + cluster_score_by_fragment.second[b];  
    }, rng);

    // How many fragment clusters are at least as good as the one at each index
    vector<size_t> better_cluster_count (max_fragment_num+1); 

    double prev_score_sum = 0.0;
    for (int rank = fragment_cluster_indices_by_score.size() - 1 ; rank >= 0 ; rank--) {
        //Go through fragment clusters in descending score order and count how many equivalent or
        //better clusters we found
        size_t fragment_num = fragment_cluster_indices_by_score[rank];
        if (rank == fragment_cluster_indices_by_score.size()-1) {
            better_cluster_count[fragment_num] = rank+1;
        } else {
            size_t prev_fragment_num = fragment_cluster_indices_by_score[rank+1];
            double curr_score_sum =  cluster_coverage_by_fragment.first[fragment_num] +
                cluster_coverage_by_fragment.second[fragment_num] +
                cluster_score_by_fragment.first[fragment_num] +  
                cluster_score_by_fragment.second[fragment_num];
            if (curr_score_sum == prev_score_sum) {
                //If this is the same as the last cluster, it has the same count
                better_cluster_count[fragment_num] = better_cluster_count[prev_fragment_num];
            } else {
                //Otherwise, its count is the index
                better_cluster_count[fragment_num] = rank+1;
                prev_score_sum = curr_score_sum;
            }
        }
    }
#ifdef debug
    for (size_t count : better_cluster_count) {
        assert(count >= 1);
    }
#endif

    // To compute the windows that are explored, we need to get
    // all the minimizers that are explored.
    vector<SmallBitset> minimizer_explored_by_read(2);
    vector<vector<size_t>> minimizer_aligned_count_by_read(2);
    //How many hits of each minimizer ended up in each extended cluster?
    vector<vector<vector<size_t>>> minimizer_extended_cluster_count_by_read(2); 
    
    // To compute the windows present in any extended cluster, we need to get
    // all the minimizers in any extended cluster.
 
    //For each fragment cluster (cluster of clusters), for each read, a vector of all alignments + the order they were fed into the funnel 
    //so the funnel can track them
    vector<pair<vector<Alignment>, vector<Alignment>>> alignments;
    vector<pair<vector<size_t>, vector<size_t>>> alignment_indices;
    pair<int, int> best_alignment_scores (0, 0); // The best alignment score for each end   

    // We will fill this with all computed alignments in estimated score order.
    // alignments has one entry for each fragment cluster and an extra for unpaired alignment 
    alignments.resize(max_fragment_num + 2);
    alignment_indices.resize(max_fragment_num + 2);

    //Now that we've scored each of the clusters, extend and align them
    for (size_t read_num = 0 ; read_num < 2 ; read_num++) {
        Alignment& aln = read_num == 0 ? aln1 : aln2;
        std::vector<Cluster>& clusters = all_clusters[read_num];
        std::vector<Minimizer>& minimizers = minimizers_by_read[read_num];
        std::vector<Seed>& seeds = seeds_by_read[read_num];
        std::vector<OldSeed>& old_seeds = old_seeds_by_read[read_num];

        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "Found " << clusters.size() << " clusters for read " << read_num << endl;
            }
        }

        // Retain clusters only if their score is better than this, in addition to the coverage cutoff
        double cluster_score_cutoff = 0.0, cluster_coverage_cutoff = 0.0, second_best_cluster_score = 0.0;

        //The score and coverage of the best cluster, "best" is determined first by coverage then score
        pair<double, double> best_cluster_coverage_score (0.0, 0.0);
        for (auto& cluster : clusters) {

            if (cluster.coverage > best_cluster_coverage_score.first) {
                //If this is the new best coverage, update both best coverage and best score
                best_cluster_coverage_score.first = cluster.coverage;
                best_cluster_coverage_score.second = cluster.score;
            } else if (cluster.coverage ==  best_cluster_coverage_score.first) {
                //If this is the same as the best coverage, get best score
                best_cluster_coverage_score.second = std::max(best_cluster_coverage_score.second, cluster.score);
            }

            cluster_coverage_cutoff = std::max(cluster_coverage_cutoff, cluster.coverage);

            if (cluster.score > cluster_score_cutoff) {
                second_best_cluster_score = cluster_score_cutoff;
                cluster_score_cutoff = cluster.score;
            } else if (cluster.score > second_best_cluster_score) {
                second_best_cluster_score = cluster.score;
            }
        }
        cluster_score_cutoff -= cluster_score_threshold;
        cluster_coverage_cutoff -= cluster_coverage_threshold;

        if (cluster_score_cutoff - pad_cluster_score_threshold < second_best_cluster_score) {
            cluster_score_cutoff = std::min(cluster_score_cutoff, second_best_cluster_score);
        }

        if (track_provenance) {
            // Now we go from clusters to gapless extensions
            funnels[read_num].stage("extend");
        }

        // These are the GaplessExtensions for all the clusters (and fragment cluster assignments), in cluster_indexes_in_order order.
        vector<pair<vector<GaplessExtension>, size_t>> cluster_extensions;
        cluster_extensions.reserve(clusters.size());

        minimizer_explored_by_read[read_num] = SmallBitset(minimizers.size());
        minimizer_aligned_count_by_read[read_num].resize(minimizers.size(), 0);
        size_t kept_cluster_count = 0;
        
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
            0, min_extensions, max_extensions, rng, [&](size_t cluster_num) {
                // Handle sufficiently good clusters 
                Cluster& cluster = clusters[cluster_num];
                if (!found_paired_cluster || fragment_cluster_has_pair[cluster.fragment] || 
                    (cluster.coverage == best_cluster_coverage_score.first &&
                     cluster.score    == best_cluster_coverage_score.second)) { 
                    //If this cluster has a pair or if we aren't looking at pairs
                    //Or if it is the best cluster
                    
                    // First check against the additional score filter
                    if (cluster_coverage_threshold != 0 && cluster.coverage < cluster_coverage_cutoff 
                            && kept_cluster_count >= min_extensions) {
                        //If the coverage isn't good enough, ignore this cluster
                        if (track_provenance) {
                            funnels[read_num].fail("cluster-coverage", cluster_num, cluster.coverage);
                        }
                        return false;
                    }
                    if (cluster_score_threshold != 0 && cluster.score < cluster_score_cutoff 
                            && kept_cluster_count >= min_extensions) {
                        //If the score isn't good enough, ignore this cluster
                        if (track_provenance) {
                            funnels[read_num].pass("cluster-coverage", cluster_num, cluster.coverage);
                            funnels[read_num].pass("max-extensions", cluster_num);
                            funnels[read_num].fail("cluster-score", cluster_num, cluster.score);
                        }
                        return false;
                    }
                    if (track_provenance) {
                        funnels[read_num].pass("cluster-coverage", cluster_num, cluster.coverage);
                        funnels[read_num].pass("max-extensions", cluster_num);
                        funnels[read_num].pass("cluster-score", cluster_num, cluster.score);
                        funnels[read_num].pass("paired-clusters", cluster_num);

                        funnels[read_num].processing_input(cluster_num);
                    }

                    if (show_work) {
                        #pragma omp critical (cerr)
                        {
                            cerr << log_name() << "Cluster " << cluster_num << endl;
                        }
                    }
                    
                    //Count how many of each minimizer is in each cluster extension
                    minimizer_extended_cluster_count_by_read[read_num].emplace_back(minimizers.size(), 0);
                    // Pack the seeds for GaplessExtender.
                    GaplessExtender::cluster_type seed_matchings;
                    if (distance_index != nullptr) {
                        for (auto seed_index : cluster.seeds) {
                            // Insert the (graph position, read offset) pair.
                            const Seed& seed = seeds[seed_index];
                            seed_matchings.insert(GaplessExtender::to_seed(seed.pos, minimizers[seed.source].value.offset));
                            minimizer_extended_cluster_count_by_read[read_num].back()[seed.source]++;
                            
                            if (show_work) {
                                #pragma omp critical (cerr)
                                {
                                    cerr << log_name() << "Seed read:" << minimizers[seed.source].value.offset << " = " << seed.pos
                                        << " from minimizer " << seed.source << ", #" << seed_index << endl;
                                }
                            }
                        }
                    } else {
                        for (auto seed_index : cluster.seeds) {
                            // Insert the (graph position, read offset) pair.
                            const OldSeed& seed = old_seeds[seed_index];
                            seed_matchings.insert(GaplessExtender::to_seed(seed.pos, minimizers[seed.source].value.offset));
                            minimizer_extended_cluster_count_by_read[read_num].back()[seed.source]++;
                            
                            if (show_work) {
                                #pragma omp critical (cerr)
                                {
                                    cerr << log_name() << "Seed read:" << minimizers[seed.source].value.offset << " = " << seed.pos
                                        << " from minimizer " << seed.source << ", #" << seed_index << endl;
                                }
                            }
                        }
                    }
                    
                    // Extend seed hits in the cluster into one or more gapless extensions
                    cluster_extensions.emplace_back(std::move(extender.extend(seed_matchings, aln.sequence())), 
                                                    cluster.fragment);
                    
                    kept_cluster_count++;
                    
                    if (show_work) {
                        #pragma omp critical (cerr)
                        {
                            cerr << log_name() << "Extensions:" << endl;
                            for (auto& e : cluster_extensions.back().first) {
                                cerr << log_name() << "\tRead " << e.read_interval.first
                                    << "-" << e.read_interval.second << " with "
                                    << e.mismatch_positions.size() << " mismatches:";
                                for (auto& pos : e.mismatch_positions) {
                                    cerr << " " << pos;
                                }
                                cerr << endl;
                            }
                        }
                    }
                    
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
                    return false;
                }
                
            }, [&](size_t cluster_num) {
                // There are too many sufficiently good clusters
                if (track_provenance) {
                    funnels[read_num].pass("cluster-coverage", cluster_num, clusters[cluster_num].coverage);
                    funnels[read_num].fail("max-extensions", cluster_num);
                }
            }, [&](size_t cluster_num) {
                // This cluster is not sufficiently good.
                // TODO: I don't think it should ever get here unless we limit the scores of the fragment clusters we look at
            });

        // We now estimate the best possible alignment score for each cluster.
        std::vector<int> cluster_extension_scores = this->score_extensions(cluster_extensions, aln, funnels[read_num]);
        
        if (track_provenance) {
            funnels[read_num].stage("align");
        }
        
        // Now start the alignment step. Everything has to become an alignment.



        
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
            extension_set_score_threshold, 2, max_alignments, rng, [&](size_t extension_num) {
                // This extension set is good enough.
                // Called in descending score order.
                
                if (track_provenance) {
                    funnels[read_num].pass("extension-set", extension_num, cluster_extension_scores[extension_num]);
                    funnels[read_num].pass("max-alignments", extension_num);
                    funnels[read_num].processing_input(extension_num);
                }
                
                auto& extensions = cluster_extensions[extension_num].first;
                
                // Collect the top alignments. Make sure we have at least one always, starting with unaligned.
                vector<Alignment> best_alignments(1, aln);
                
                if (GaplessExtender::full_length_extensions(extensions)) {
                    // We got full-length extensions, so directly convert to an Alignment.
                    
                    if (track_provenance) {
                        funnels[read_num].substage("direct");
                    }

                    //Fill in the best alignments from the extension. We know the top one is always full length and exists.
                    this->extension_to_alignment(extensions.front(), best_alignments.front());
                    
                    
                    if (show_work) {
                        #pragma omp critical (cerr)
                        {
                            cerr << log_name() << "Produced alignment directly from full length gapless extension " << extension_num << endl;
                        }
                    }
                    
                    for (auto next_ext_it = extensions.begin() + 1; next_ext_it != extensions.end() && next_ext_it->full(); ++next_ext_it) {
                        // For all subsequent full length extensions, make them into alignments too.
                        // We want them all to go on to the pairing stage so we don't miss a possible pairing in a tandem repeat.
                        best_alignments.emplace_back(aln);
                        this->extension_to_alignment(*next_ext_it, best_alignments.back());
                        
                        if (show_work) {
                            #pragma omp critical (cerr)
                            {
                                cerr << log_name() << "Produced additional alignment directly from full length gapless extension " << (next_ext_it - extensions.begin()) << endl;
                            }
                        }
                        
                    }

                    if (track_provenance) {
                        // Stop the current substage
                        funnels[read_num].substage_stop();
                    }
                } else if (do_dp) {
                    // We need to do base-level alignment.
                    
                    if (track_provenance) {
                        funnels[read_num].substage("align");
                    }
                    
                    // Do the DP and compute up to 2 alignments
                    best_alignments.emplace_back(aln);
                    find_optimal_tail_alignments(aln, extensions, rng, best_alignments[0], best_alignments[1]);

                    
                    if (track_provenance) {
                        // We're done base-level alignment. Next alignment may not go through this substage.
                        funnels[read_num].substage_stop();
                    }
                } else {
                    // We would do base-level alignment but it is disabled.
                    // Leave best_alignments unaligned
                }
                
                
                size_t fragment_num = cluster_extensions[extension_num].second;
                    
                // Have a function to process the best alignments we obtained
                auto observe_alignment = [&](Alignment& aln) {
                    auto& best_score = read_num == 0 ? best_alignment_scores.first : best_alignment_scores.second;
                    best_score = max(best_score, aln.score());
                    
                    auto& alignment_list = read_num == 0 ? alignments[fragment_num].first 
                                                         : alignments[fragment_num].second;
                    alignment_list.emplace_back(std::move(aln));
                    
                    auto& indices_list = read_num == 0 ? alignment_indices[fragment_num].first 
                                                       : alignment_indices[fragment_num].second;
                    indices_list.emplace_back(curr_funnel_index);
                    curr_funnel_index++;

                    if (track_provenance) {
                        funnels[read_num].project(extension_num);
                        funnels[read_num].score(funnels[read_num].latest(), alignment_list.back().score());
                    }
                    
                    if (show_work) {
                        #pragma omp critical (cerr)
                        {
                            cerr << log_name() << "Produced fragment option " << fragment_num << " end " << read_num
                                << " alignment with score " << alignment_list.back().score() << ": " << pb2json(alignment_list.back()) << endl;
                        }
                    }
                };
                
                for(auto aln_it = best_alignments.begin() ; aln_it != best_alignments.end() && aln_it->score() != 0 && aln_it->score() >= best_alignments[0].score() * 0.8; ++aln_it) {
                    //For each additional extension with score at least 0.8 of the best score
                    observe_alignment(*aln_it);
                }

                if (track_provenance) {
                    // We're done with this input item
                    funnels[read_num].processed_input();
                }
                
                for (size_t i = 0 ; i < minimizer_extended_cluster_count_by_read[read_num][extension_num].size() ; i++) {
                    if (minimizer_extended_cluster_count_by_read[read_num][extension_num][i] > 0) {
                        // This minimizer is in an extended cluster that gave rise
                        // to at least one alignment, so it is explored.
                        minimizer_explored_by_read[read_num].insert(i);
                        minimizer_aligned_count_by_read[read_num][i] += minimizer_extended_cluster_count_by_read[read_num][extension_num][i];
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


    //Now that we have alignments, figure out how to pair them up
    
    if (track_provenance) {
        // Now say we are finding the pairs
        funnels[0].stage("pairing");
        funnels[1].stage("pairing");
    }
    // Fill this in with the indexes of pairs of alignments we will output
    // each alignment is stored as <fragment index, alignment index> into alignments
    // fragment_index should be the same for both ends, unless one was rescued
    vector<pair<pair<size_t, size_t>, pair<size_t, size_t>>> paired_alignments;
    paired_alignments.reserve(alignments.size());


#ifdef print_minimizer_table
    vector<pair<bool, bool>> alignment_was_rescued;
#endif

    //For each alignment in alignments, which paired_alignment includes it. Follows structure of alignments
    vector<pair<vector<vector<size_t>>, vector<vector<size_t>>>> alignment_groups(alignments.size());

    // Grab all the scores in order for MAPQ computation.
    vector<double> paired_scores;
    paired_scores.reserve(alignments.size());
    vector<int64_t> fragment_distances;
    fragment_distances.reserve(alignments.size());

    //for each alignment pair, what type of pair is it
    enum PairType {paired, unpaired, rescued_from_first, rescued_from_second};
    vector<PairType> pair_types; 
    
    //For each pair of alignments in paired_alignments, how many equivalent or better fragment clusters
    //did we find 
    vector<size_t> better_cluster_count_by_pairs;


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
            for (size_t aln_index1 = 0 ; aln_index1 < fragment_alignments.first.size() ; aln_index1++)  {
                Alignment& alignment1 = fragment_alignments.first[aln_index1];
                size_t funnel_index1 = alignment_indices[fragment_num].first[aln_index1];
                for (size_t aln_index2 = 0 ; aln_index2 < fragment_alignments.second.size() ; aln_index2++) {
                    Alignment& alignment2 = fragment_alignments.second[aln_index2];
                    size_t funnel_index2 = alignment_indices[fragment_num].second[aln_index2];

                    //Get the likelihood of the fragment distance
                    int64_t fragment_distance = distance_between(alignment1, alignment2); 

                    double score = score_alignment_pair(alignment1, alignment2, fragment_distance);
                    alignment_groups[fragment_num].first[aln_index1].emplace_back(paired_alignments.size());
                    alignment_groups[fragment_num].second[aln_index2].emplace_back(paired_alignments.size());
                    paired_alignments.emplace_back(make_pair(fragment_num, aln_index1), make_pair(fragment_num, aln_index2));
                    paired_scores.emplace_back(score);
                    fragment_distances.emplace_back(fragment_distance);
                    better_cluster_count_by_pairs.emplace_back(better_cluster_count[fragment_num]);
                    pair_types.emplace_back(paired);
#ifdef print_minimizer_table
                    alignment_was_rescued.emplace_back(false, false);
#endif

                    if (show_work) {
                        #pragma omp critical (cerr)
                        {
                            cerr << log_name() << "Found pair of alignments from fragment " << fragment_num << " with scores " 
                                 << alignment1.score() << " " << alignment2.score() << " at distance " << fragment_distance 
                                 << " gets pair score " << score << endl;
                            cerr << log_name() << "Alignment 1: " << pb2json(alignment1) << endl << "Alignment 2: " << pb2json(alignment2) << endl;
                        }
                    }

                    if (track_provenance) {
                        funnels[0].processing_input(funnel_index1);
                        funnels[1].processing_input(funnel_index2);
                        funnels[0].substage("pair-clusters");
                        funnels[1].substage("pair-clusters");
                        funnels[0].pass("max-rescue-attempts", funnel_index1);
                        funnels[0].project(funnel_index1);
                        funnels[0].score(funnels[0].latest(), score);
                        funnels[1].pass("max-rescue-attempts", funnel_index2);
                        funnels[1].project(funnel_index2);
                        funnels[1].score(funnels[1].latest(), score);
                        funnels[0].substage_stop();
                        funnels[1].substage_stop();
                        funnels[0].processed_input();
                        funnels[1].processed_input();
                    }
                }
            }
        } else if (!fragment_alignments.first.empty()) {
            //If this fragment cluster has only alignments from the first read
            
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Found unpaired alignments from fragment " << fragment_num << " for first read" << endl;
                }
            }
            for (size_t i = 0 ; i < fragment_alignments.first.size() ; i++) {
                unpaired_alignments.emplace_back(fragment_num, i, true);
                unpaired_count_1++;
                if (show_work) {
                    #pragma omp critical (cerr)
                    {
                        cerr << log_name() << "\t" << pb2json(fragment_alignments.first[i]) << endl;
                    }
                }
            }
        } else if (!fragment_alignments.second.empty()) {
            //If this fragment cluster has only alignments from the second read
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Found unpaired alignments from fragment " << fragment_num << " for second read" << endl;
                }
            }
            for (size_t i = 0 ; i < fragment_alignments.second.size() ; i++) {
                unpaired_alignments.emplace_back(fragment_num, i, false);
                unpaired_count_2++;
                
                if (show_work) {
                    #pragma omp critical (cerr)
                    {
                        cerr << log_name() << "\t" << pb2json(fragment_alignments.second[i]) << endl;
                    }
                }
            }
        }
    }
    size_t rescued_count_1 = 0;
    size_t rescued_count_2 = 0;
    vector<vector<double>> unpaired_scores(2) ;
    unpaired_scores[0].reserve(unpaired_alignments.size());
    unpaired_scores[1].reserve(unpaired_alignments.size());


    if (!unpaired_alignments.empty()) {
        //If we found some clusters that had no pair in a fragment cluster
        if (!found_pair) {
            //If we didn't find any pairs find the best alignment for each end

            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Found no pairs and we aren't doing rescue: return best alignment for each read" << endl;
                }
            }
            tuple<size_t, size_t, size_t> best_index_1(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
            tuple<size_t, size_t, size_t> best_index_2(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
            int32_t best_score_1 = 0;
            int32_t best_score_2 = 0;

            for (tuple<size_t, size_t, bool> index : unpaired_alignments ) {
                Alignment& alignment = std::get<2>(index) ? alignments[std::get<0>(index)].first[std::get<1>(index) ]
                                                          : alignments[std::get<0>(index)].second[std::get<1>(index)];
                if (std::get<2>(index)) {
                    unpaired_scores[0].emplace_back(alignment.score());
                    if (deterministic_beats(alignment.score(), best_score_1, rng)) {
                        best_index_1 = index;
                        best_score_1 = alignment.score();
                    }
                } else {
                    unpaired_scores[1].emplace_back(alignment.score());
                    if (deterministic_beats(alignment.score(), best_score_2, rng)) {
                        best_index_2 = index;
                        best_score_2 = alignment.score();
                    }
                }
            }
            if (max_rescue_attempts == 0 ) { 
                //If we aren't attempting rescue, just return the best alignment from each end
                Alignment& best_aln1 = aln1;
                Alignment& best_aln2 = aln2;
    
                if (std::get<0>(best_index_1) != std::numeric_limits<size_t>::max()) {
                    //If there was a best alignment for 1, use it
                    best_aln1 = alignments[std::get<0>(best_index_1)].first[std::get<1>(best_index_1)]; 
                } else {
                    //Otherwise return an empty alignment
                    best_aln1.clear_refpos();
                    best_aln1.clear_path();
                    best_aln1.set_score(0);
                    best_aln1.set_identity(0);
                    best_aln1.set_mapping_quality(0);
                }

                if (std::get<0>(best_index_2) != std::numeric_limits<size_t>::max()) {
                    //If there was a best alignment for 2, use it
                    best_aln2 = alignments[std::get<0>(best_index_2)].second[std::get<2>(best_index_2)]; 
                } else {
                    //Otherwise return an empty alignment
                    best_aln2.clear_refpos();
                    best_aln2.clear_path();
                    best_aln2.set_score(0);
                    best_aln2.set_identity(0);
                    best_aln2.set_mapping_quality(0);
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
                
                // Annotate with whatever's in the funnel
                funnels[0].annotate_mapped_alignment(paired_mappings.first[0], track_correctness);
                funnels[0].annotate_mapped_alignment(paired_mappings.second[0], track_correctness);
                
                return paired_mappings;
            } else if (best_score_1 != 0 and best_score_2 != 0) {
                //We are attempting rescue, but we still want to keep the best alignments as a potential (unpaired) pair
                
                size_t funnel_index1, funnel_index2;
                if (track_provenance) {
                    // Work out what the paired-up alignments are numbered in the funnel.
                    // TODO: can we flatten these lookup paths or change tuples
                    // to structs to be more understandable?
                    funnel_index1 = alignment_indices[std::get<0>(best_index_1)].first[std::get<1>(best_index_1)];
                    funnel_index2 = alignment_indices[std::get<0>(best_index_2)].second[std::get<1>(best_index_2)];
                    funnels[0].processing_input(funnel_index1);
                    funnels[1].processing_input(funnel_index2);
                }
                
                pair<pair<size_t, size_t>, pair<size_t, size_t>> index_pair =  make_pair(make_pair(std::get<0>(best_index_1), std::get<1>(best_index_1)), 
                                                                                         make_pair(std::get<0>(best_index_2), std::get<1>(best_index_2)));

                Alignment& aln1 = alignments[std::get<0>(best_index_1)].first[std::get<1>(best_index_1)];
                Alignment& aln2 = alignments[std::get<0>(best_index_2)].second[std::get<1>(best_index_2)];
                paired_alignments.push_back(index_pair);
                //Assume the distance between them is infinite
                double pair_score = score_alignment_pair(aln1, aln2, std::numeric_limits<int64_t>::max());
                paired_scores.emplace_back(pair_score);
                fragment_distances.emplace_back(std::numeric_limits<int64_t>::max());
                better_cluster_count_by_pairs.emplace_back(0);
                pair_types.emplace_back(unpaired);
                
                if (track_provenance) {
                    funnels[0].substage("pair-clusters");
                    funnels[1].substage("pair-clusters");
                    funnels[0].project(funnel_index1);
                    funnels[0].score(funnels[0].latest(), pair_score);
                    funnels[1].project(funnel_index2);
                    funnels[1].score(funnels[1].latest(), pair_score);
                    funnels[0].substage_stop();
                    funnels[1].substage_stop();
                    funnels[0].processed_input();
                    funnels[1].processed_input();
                }
            }
        }

        if (max_rescue_attempts != 0) {
            //Attempt rescue on unpaired alignments if either we didn't find any pairs or if the unpaired alignments are very good

            process_until_threshold_a(unpaired_alignments, (std::function<double(size_t)>) [&](size_t i) -> double{
                tuple<size_t, size_t, bool>& index = unpaired_alignments.at(i);
                return (double) std::get<2>(index) ? alignments[std::get<0>(index)].first[std::get<1>(index)].score()
                                                   : alignments[std::get<0>(index)].second[std::get<1>(index)].score();
            }, 0, 1, max_rescue_attempts, rng, [&](size_t i) {
                tuple<size_t, size_t, bool>& index = unpaired_alignments.at(i);
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

                //Rescue the alignment
                attempt_rescue(mapped_aln, rescued_aln, minimizers_by_read[(found_first ? 1 : 0)], found_first);

                if (rescued_aln.path().mapping_size() != 0) {
                    //If we actually found an alignment

                    int64_t fragment_dist = found_first ? distance_between(mapped_aln, rescued_aln) 
                                                      : distance_between(rescued_aln, mapped_aln);

                    double score = score_alignment_pair(mapped_aln, rescued_aln, fragment_dist);

                    set_annotation(mapped_aln, "rescuer", true);
                    set_annotation(rescued_aln, "rescued", true);
                    set_annotation(mapped_aln,  "fragment_length", distance_to_annotation(fragment_dist));
                    set_annotation(rescued_aln, "fragment_length", distance_to_annotation(fragment_dist));
                    bool properly_paired = fragment_dist == std::numeric_limits<int64_t>::max() ? false :
                        (std::abs(fragment_dist-fragment_length_distr.mean()) <= 6.0*fragment_length_distr.std_dev()) ;
                    set_annotation(mapped_aln, "proper_pair", properly_paired);
                    set_annotation(rescued_aln, "proper_pair", properly_paired);

                    //Since we're still accumulating a list of indexes of pairs of alignments,
                    //add the new alignment to the list of alignments 
                    //(in a separate "fragment cluster" vector for rescued alignments) and keep track of its index
                    //
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
                    pair_types.push_back(found_first ? rescued_from_first : rescued_from_second); 
                    better_cluster_count_by_pairs.emplace_back(better_cluster_count[mapped_index.first]);

#ifdef print_minimizer_table
                    alignment_was_rescued.emplace_back(!found_first, found_first);
#endif
                    if (track_provenance) {
                        funnels[found_first ? 0 : 1].pass("max-rescue-attempts", j);
                        funnels[found_first ? 0 : 1].project(j);
                        funnels[found_first ? 1 : 0].introduce();
                        funnels[0].score(funnels[0].latest(), score);
                        funnels[1].score(funnels[1].latest(), score);
                    }
                }
                if (track_provenance) {
                    funnels[found_first ? 0 : 1].processed_input();
                    funnels[found_first ? 0 : 1].substage_stop();
                }
                return true;
            }, [&](size_t i) {
                //This alignment is good enough but we already rescued enough
                if (track_provenance) {
                    tuple<size_t, size_t, bool> index = unpaired_alignments.at(i);
                    bool found_first = std::get<2>(index); 
                    size_t j = found_first ? alignment_indices[std::get<0>(index)].first[std::get<1>(index)]
                                            : alignment_indices[std::get<0>(index)].second[std::get<1>(index)];
                    funnels[found_first ? 0 : 1].fail("max-rescue-attempts", j);
                }
                return;
            }, [&] (size_t i) {
                //This alignment is insufficiently good
                if (track_provenance) {
                    //TODO: Fail something here
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
    vector<PairType> types; 
    mappings.first.reserve(paired_alignments.size());
    mappings.second.reserve(paired_alignments.size());
    scores.reserve(paired_scores.size());
    distances.reserve(fragment_distances.size());
    types.reserve(pair_types.size());

    //For each pair of alignments in mappings, how many equivalent or better fragment clusters were there
    vector<size_t> better_cluster_count_by_mappings;

#ifdef print_minimizer_table
    vector<pair<bool, bool>> mapping_was_rescued;
    vector<pair<pair<size_t, size_t>, pair<size_t, size_t>>> pair_indices;
#endif

    process_until_threshold_a(paired_alignments, (std::function<double(size_t)>) [&](size_t i) -> double {
        return paired_scores[i];
    }, 0, 1, max_multimaps, rng, [&](size_t alignment_num) {
        // This alignment makes it
        // Called in score order

        pair<pair<size_t, size_t>, pair<size_t, size_t>> index_pair = paired_alignments[alignment_num];
        
        
        // Remember the score at its rank
        scores.emplace_back(paired_scores[alignment_num]);
        distances.emplace_back(fragment_distances[alignment_num]);
        types.emplace_back(pair_types[alignment_num]);
        better_cluster_count_by_mappings.emplace_back(better_cluster_count_by_pairs[alignment_num]);
        // Remember the output alignment
        mappings.first.emplace_back(alignments[index_pair.first.first].first[index_pair.first.second]);
        mappings.second.emplace_back(alignments[index_pair.second.first].second[index_pair.second.second]);

        if (mappings.first.size() == 1 && found_pair) {
            //If this is the best pair of alignments that we're going to return and we didn't attempt rescue, 
            //get the group scores for mapq

            //Get the scores of this pair 
            scores_group_1.push_back(paired_scores[alignment_num]);
            scores_group_2.push_back(paired_scores[alignment_num]);

            //The indices (into paired_alignments) of pairs with the same first read as this
            vector<size_t>& alignment_group_1 = alignment_groups[index_pair.first.first].first[index_pair.first.second];
            //And second read
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

#ifdef print_minimizer_table
        mapping_was_rescued.emplace_back(alignment_was_rescued[alignment_num]);
        pair_indices.push_back(index_pair);
#endif
        
        if (track_provenance) {
            // Tell the funnel
            funnels[0].pass("max-multimaps", alignment_num);
            funnels[0].project(alignment_num);
            funnels[0].score(funnels[0].latest(), scores.back());
            funnels[1].pass("max-multimaps", alignment_num);
            funnels[1].project(alignment_num);
            funnels[1].score(funnels[1].latest(), scores.back());
        }
        
        return true;
    }, [&](size_t alignment_num) {
        // We already have enough alignments, although this one has a good score
        
        // Remember the score at its rank anyway
        scores.emplace_back(paired_scores[alignment_num]);
        distances.emplace_back(fragment_distances[alignment_num]);
        types.emplace_back(pair_types[alignment_num]);
        better_cluster_count_by_mappings.emplace_back(better_cluster_count_by_pairs[alignment_num]);

 
#ifdef print_minimizer_table
        pair<pair<size_t, size_t>, pair<size_t, size_t>> index_pair = paired_alignments[alignment_num];
        pair_indices.push_back(index_pair);
#endif       
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

    
    // Compute raw explored caps (with 2.0 scaling, like for single-end) and raw group caps.
    // Non-capping caps stay at infinity.
    vector<double> mapq_explored_caps(2, std::numeric_limits<float>::infinity());
    vector<double> mapq_score_groups(2, std::numeric_limits<float>::infinity());
    // We also have one fragment_cluster_cap across both ends.
    double fragment_cluster_cap = std::numeric_limits<float>::infinity();
    // And one base uncapped MAPQ
    double uncapped_mapq = 0;
    double new_cluster_cap = numeric_limits<double>::infinity();

    // Store multiplicities, if we fill them in
    vector<double> paired_multiplicities;

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
#ifdef print_minimizer_table
        mapping_was_rescued.emplace_back(false, false);
        pair_indices.emplace_back(make_pair(std::numeric_limits<size_t>::infinity(), std::numeric_limits<size_t>::infinity()), 
                        make_pair(std::numeric_limits<size_t>::infinity(), std::numeric_limits<size_t>::infinity()));
#endif

    } else {
    
        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "For scores";
                for (auto& score : scores) cerr << " " << score << ":" << endl;
            }
        }

        //Get the multiplicities for mapq calculation
        //We're only using multiplicities if the alignments were rescued
        double estimated_multiplicity_from_1 = unpaired_count_1 > 0 ? (double) unpaired_count_1 / min(rescued_count_1, max_rescue_attempts) : 1.0;
        double estimated_multiplicity_from_2 = unpaired_count_2 > 0 ? (double) unpaired_count_2 / min(rescued_count_2, max_rescue_attempts) : 1.0;
        bool all_rescued = true;
        for (PairType type : types) {
            switch (type){
                case paired:
                    paired_multiplicities.push_back(1.0);
                    all_rescued=false;
                case unpaired:
                    paired_multiplicities.push_back(1.0);
                case rescued_from_first:
                    paired_multiplicities.push_back(estimated_multiplicity_from_1);
                case rescued_from_second:
                    paired_multiplicities.push_back(estimated_multiplicity_from_2);
            }
        }
        const vector<double>* multiplicities = all_rescued ? &paired_multiplicities : nullptr; 
        // Compute base MAPQ if not unmapped. Otherwise use 0 instead of the 50% this would give us.
        // If all of the alignment pairs were found with rescue, use the multiplicities to determine mapq
        // Use exact mapping quality
        uncapped_mapq = scores[0] == 0 ? 0 : 
            get_regular_aligner()->compute_max_mapping_quality(scores, false, multiplicities);

        //Cap mapq at 1 - 1 / # equivalent or better fragment clusters, including self
        if (better_cluster_count_by_mappings.front() > 1) {
            // TODO: why is this a sensible cap?
            fragment_cluster_cap = prob_to_phred(1.0 - (1.0 / (double) better_cluster_count_by_mappings.front()));
            // Leave zeros in here and don't round.
        }

        //If one alignment was duplicated in other pairs, cap the mapq for that alignment at the mapq
        //of the group of duplicated alignments. Always compute this even if not quite sensible.
        mapq_score_groups[0] = get_regular_aligner()->compute_max_mapping_quality(scores_group_1, false);
        mapq_score_groups[1] = get_regular_aligner()->compute_max_mapping_quality(scores_group_2, false);
        
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
            // Compute exploration cap on MAPQ. TODO: avoid needing to pass as much stuff along.
            double mapq_explored_cap =  faster_cap(minimizers_by_read[read_num], explored_minimizers, aln.sequence(), aln.quality());

            mapq_explored_caps[read_num] = mapq_explored_cap;

            // Remember the caps
            auto& to_annotate = (read_num == 0 ? mappings.first : mappings.second).front();
            set_annotation(to_annotate, "mapq_explored_cap", mapq_explored_cap);
            set_annotation(to_annotate, "mapq_score_group", mapq_score_groups[read_num]);
        }
        
        // Have a function to transform interesting cap values to uncapped.
        auto preprocess_cap = [&](double cap) {
            return (cap != -numeric_limits<double>::infinity()) ? cap : numeric_limits<double>::infinity();
        };
        
        for (auto read_num : {0, 1}) {
            // For each fragment

            // Compute the overall cap for just this read, now that the individual cap components are filled in for both reads.
            double escape_bonus = uncapped_mapq < std::numeric_limits<int32_t>::max() ? 1.0 : 2.0;
            double mapq_cap = std::min(fragment_cluster_cap, ((mapq_explored_caps[0] + mapq_explored_caps[1])*escape_bonus) );

            //TODO: How to cap mapq when the reads were unpaired
            if (types.front() == unpaired) {
                //If this pair came from two different fragment cluster, then cap mapq at the mapq
                //from only unpaired alignments of this read
                mapq_cap = std::min(mapq_cap, (double)get_regular_aligner()->compute_max_mapping_quality(unpaired_scores[read_num], false));
            }
            
            // Find the MAPQ to cap
            double read_mapq = uncapped_mapq;
            
            // Remember the uncapped MAPQ
            auto& to_annotate = (read_num == 0 ? mappings.first : mappings.second).front();
            set_annotation(to_annotate, "mapq_uncapped", read_mapq);
            // And the cap we actually applied (possibly from the pair partner)
            set_annotation(to_annotate, "mapq_applied_cap", mapq_cap);

            // Apply the cap, and limit to 0-60
            double capped_mapq = min(mapq_cap, read_mapq); 
            if (distances.front() == std::numeric_limits<int64_t>::max()) {
                //If the two reads are not reachable, lower cap
                capped_mapq = capped_mapq / 2.0;
            }
            read_mapq = max(min(capped_mapq, 120.0) / 2.0, 0.0);
            
            // Save the MAPQ
            to_annotate.set_mapping_quality(read_mapq);
            
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "MAPQ for read " << read_num << " is " << read_mapq << ", was " << uncapped_mapq
                        << " capped by fragment cluster cap " << fragment_cluster_cap
                        << ", score group cap " << (mapq_score_groups[read_num] / 2.0)
                        << ", combined explored cap " << ((mapq_explored_caps[0] + mapq_explored_caps[1]) / 2.0)  << endl;
                }
            }
        }
        
        //Annotate top pair with its fragment distance, fragment length distrubution, and secondary scores
        set_annotation(mappings.first.front(), "fragment_length", distance_to_annotation(distances.front()));
        set_annotation(mappings.second.front(), "fragment_length", distance_to_annotation(distances.front()));
        bool properly_paired = distances.front() == std::numeric_limits<int64_t>::max() ? false :
            (std::abs(distances.front()-fragment_length_distr.mean()) <= 6.0*fragment_length_distr.std_dev()) ;
        set_annotation(mappings.first.front(), "proper_pair", properly_paired);
        set_annotation(mappings.second.front(), "proper_pair", properly_paired);

        string distribution = "-I " + to_string(fragment_length_distr.mean()) + " -D " + to_string(fragment_length_distr.std_dev());
        set_annotation(mappings.first.front(),"fragment_length_distribution", distribution);
        set_annotation(mappings.second.front(),"fragment_length_distribution", distribution);
        set_annotation(mappings.first.front(),"secondary_scores", scores);
        set_annotation(mappings.second.front(),"secondary_scores", scores);
    
    }
    
    // Make sure pair partners reference each other
    pair_all(mappings);
        
    if (track_provenance) {
        funnels[0].substage_stop();
        funnels[1].substage_stop();
    }
    
    // Stop this alignment
    funnels[0].stop();
    funnels[1].stop();
    
    // Annotate with whatever's in the funnel.
    funnels[0].annotate_mapped_alignment(mappings.first[0], track_correctness);
    funnels[1].annotate_mapped_alignment(mappings.second[0], track_correctness);
    
    if (track_provenance) {
        if (track_correctness) {
            if (distance_index != nullptr) {
                annotate_with_minimizer_statistics(mappings.first[0], minimizers_by_read[0], seeds_by_read[0], funnels[0]);
                annotate_with_minimizer_statistics(mappings.second[0], minimizers_by_read[1], seeds_by_read[1], funnels[1]);
            } else {
                old_annotate_with_minimizer_statistics(mappings.first[0], minimizers_by_read[0], old_seeds_by_read[0], funnels[0]);
                old_annotate_with_minimizer_statistics(mappings.second[0], minimizers_by_read[1], old_seeds_by_read[1], funnels[1]);
            }
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
 
#ifdef print_minimizer_table

    if (distances.size() == 0) {
        distances.emplace_back(0);
    }
    cerr << aln1.sequence() << "\t";
    for (char c : aln1.quality()) {
        cerr << (char)(c+33);
    }
    cerr << "\t" << max_fragment_num << "\t" << mapping_was_rescued[0].first << "\t" << mapping_was_rescued[0].second 
         << "\t" << distances.front();
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
    cerr << "\t" << uncapped_mapq << "\t" << fragment_cluster_cap << "\t" << mapq_score_groups[0] << "\t" 
         << mapq_explored_caps[0] << "\t" << new_cluster_cap << "\t" << mappings.first.front().mapping_quality() << "\t";  
    for (size_t i = 0 ; i < scores.size() ; i++) {
        pair<pair<size_t, size_t>, pair<size_t, size_t>> indices = pair_indices[i];
        Alignment& aln_1 = alignments[indices.first.first].first[indices.first.second];
        Alignment& aln_2 = alignments[indices.second.first].second[indices.second.second];

        int64_t dist = distances[i];
        assert(dist == distance_between(aln_1, aln_2)); 

        assert(scores[i] == score_alignment_pair(aln_1, aln_2, dist));

        double multiplicity = paired_multiplicities.size() == scores.size() ? paired_multiplicities[i] : 1.0;

        cerr << aln_1.score() << "," 
             << aln_2.score() << "," 
             << multiplicity << "," 
             << scores[i] << ";";
    }

    
    if (track_correctness) {
        cerr << "\t" << funnels[0].last_correct_stage() << endl;
    } else {
        cerr << "\t?" << endl;
    }

    cerr << aln2.sequence() << "\t";
    for (char c : aln2.quality()) {
        cerr << (char)(c+33);
    }
    cerr << "\t" << max_fragment_num << "\t" << mapping_was_rescued[0].second << "\t" << mapping_was_rescued[0].first 
         << "\t" << distances.front();
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
    cerr << "\t" << uncapped_mapq << "\t" << fragment_cluster_cap << "\t" << mapq_score_groups[1] << "\t" 
         << mapq_explored_caps[1] << "\t" << new_cluster_cap << "\t" << mappings.second.front().mapping_quality() << "\t";

    for (size_t i = 0 ; i < scores.size() ; i++) {
        pair<pair<size_t, size_t>, pair<size_t, size_t>> indices = pair_indices[i];
        Alignment& aln_1 = alignments[indices.first.first].first[indices.first.second];
        Alignment& aln_2 = alignments[indices.second.first].second[indices.second.second];

        int64_t dist = distances[i];
        assert(dist == distance_between(aln_1, aln_2)); 

        assert(scores[i] == score_alignment_pair(aln_1, aln_2, dist));

        double multiplicity = paired_multiplicities.size() == scores.size() ? paired_multiplicities[i] : 1.0;

        cerr << aln_1.score() << "," 
             << aln_2.score() << "," 
             << multiplicity << "," 
             << scores[i] << ";";
    }

    if (track_correctness) {
        cerr << "\t" << funnels[1].last_correct_stage() << endl;
    } else {
        cerr << "\t?" << endl;
    }
#endif

    if (track_provenance && show_work) {
        // Dump the funnel info graph.
        #pragma omp critical (cerr)
        {
            funnels[0].to_dot(cerr);
            funnels[1].to_dot(cerr);
        }
    }

    // Ship out all the aligned alignments
    return mappings;
}

//-----------------------------------------------------------------------------

double MinimizerMapper::faster_cap(const vector<Minimizer>& minimizers, vector<size_t>& minimizers_explored,
    const string& sequence, const string& quality_bytes) {

    // TODO: Maybe we should to something smarter if we do not have base qualities.
    if (quality_bytes.empty()) {
        return numeric_limits<double>::infinity();
    }

    // Sort minimizer subset so we go through minimizers in increasing order of start position
    std::sort(minimizers_explored.begin(), minimizers_explored.end(), [&](size_t a, size_t b) {
        // Return true if a must come before b, and false otherwise
        return minimizers[a].forward_offset() < minimizers[b].forward_offset();
    });
#ifdef debug
    cerr << "Sorted " << minimizers_explored.size() << " minimizers" << endl;
#endif

#ifdef debug
    cerr << "Explored minimizers:" << endl;
    dump_debug_minimizers(minimizers, sequence, &minimizers_explored);
#endif

    // Make a DP table holding the log10 probability of having an error disrupt each minimizer.
    // Entry i+1 is log prob of mutating minimizers 0, 1, 2, ..., i.
    // Make sure to have an extra field at the end to support this.
    // Initialize with -inf for unfilled.
    vector<double> c(minimizers_explored.size() + 1, -numeric_limits<double>::infinity());
    c[0] = 0.0;
    
    for_each_agglomeration_interval(minimizers, sequence, quality_bytes, minimizers_explored, [&](size_t left, size_t right, size_t bottom, size_t top) {
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
    
#ifdef debug
    cerr << "log10prob after all minimizers is " << c.back() << endl;
#endif
    
    assert(!isinf(c.back()));
    // Conver to Phred.
    double result = -c.back() * 10;
    return result;
}

void MinimizerMapper::for_each_agglomeration_interval(const vector<Minimizer>& minimizers,
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
    
#ifdef debug
    cerr << "Compute log10 probability in interval " << left << "-" << right << endl;
#endif
    
    if (left == right) {
        // 0-length intervals need no disruption.
        return 0;
    }
   
    // We want an OR over all the columns, but some of the probabilities are tiny.
    // So instead of NOT(AND(NOT())), which also would assume independence the
    // way we calculate AND by multiplication, we just assume independence and
    // compute OR as (p1 + p2 - (p1 * p2)).
    // Start with the first column.
    double p = get_prob_of_disruption_in_column(minimizers, sequence, quality_bytes, disrupt_begin, disrupt_end, left);
#ifdef debug
    cerr << "\tProbability disrupted at column " << left << ": " << p << endl;
#endif
    for(size_t i = left + 1 ; i < right; i++) {
        // OR up probability of all the other columns
        double col_p = get_prob_of_disruption_in_column(minimizers, sequence, quality_bytes, disrupt_begin, disrupt_end, i);
#ifdef debug
        cerr << "\tProbability disrupted at column " << i << ": " << col_p << endl;
#endif
        p = (p + col_p - (p * col_p));
#ifdef debug
        cerr << "\tRunning OR of disrupted anywhere: " << p << endl;
#endif
    }
    
    // Convert to log10prob.
    return log10(p);
 
}

double MinimizerMapper::get_prob_of_disruption_in_column(const vector<Minimizer>& minimizers,
    const string& sequence, const string& quality_bytes,
    const vector<size_t>::iterator& disrupt_begin, const vector<size_t>::iterator& disrupt_end,
    size_t index) {
    
#ifdef debug
    cerr << "\tCompute probability at column " << index << endl;
#endif
    
    // Base cost is quality. Make sure to compute a non-integral answer.
    double p = phred_to_prob((uint8_t)quality_bytes[index]);
#ifdef debug
    cerr << "\t\tBase probability from quality: " << p << endl;
#endif
    for (auto it = disrupt_begin; it != disrupt_end; ++it) {
        // For each minimizer to disrupt
        auto& m = minimizers[*it];
        
#ifdef debug
        cerr << "\t\tRelative rank " << (it - disrupt_begin) << " is minimizer " << m.value.key.decode(m.length) << endl;
#endif
        
        if (!(m.forward_offset() <= index && index < m.forward_offset() + m.length)) {
            // Index is out of range of the minimizer itself. We're in the flank.
#ifdef debug
            cerr << "\t\t\tColumn " << index << " is in flank." << endl;
#endif
            // How many new possible minimizers would an error here create in this agglomeration,
            // to compete with its minimizer?
            // No more than one per position in a minimizer sequence.
            // No more than 1 per base from the start of the agglomeration to here, inclusive.
            // No more than 1 per base from here to the last base of the agglomeration, inclusive.
            size_t possible_minimizers = min((size_t) m.length,
                                             min(index - m.agglomeration_start + 1,
                                             (m.agglomeration_start + m.agglomeration_length) - index));

            // Account for at least one of them beating the minimizer.
            double any_beat_prob = prob_for_at_least_one(m.value.hash, possible_minimizers);
            
#ifdef debug
            cerr << "\t\t\tBeat hash " << m.value.hash << " at least 1 time in " << possible_minimizers << " gives probability: " << any_beat_prob << endl;
#endif
            
            p *= any_beat_prob;
            
            // TODO: handle N somehow??? It can occur outside the minimizer itself, here in the flank.
        }
#ifdef debug
        cerr << "\t\t\tRunning AND prob: " << p << endl;
#endif
    }
    
    return p;
}

//-----------------------------------------------------------------------------

void MinimizerMapper::attempt_rescue(const Alignment& aligned_read, Alignment& rescued_alignment, const std::vector<Minimizer>& minimizers, bool rescue_forward ) {

    // Get rid of the old path.
    rescued_alignment.clear_path();

    if (this->rescue_algorithm == rescue_none) { return; }

    // We are traversing the same small subgraph repeatedly, so it's better to use a cache.
    gbwtgraph::CachedGBWTGraph cached_graph(this->gbwt_graph);

    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "Attempt rescue from: " << pb2json(aligned_read) << endl;
        }
    }

    // Find all nodes within a reasonable range from aligned_read.
    std::unordered_set<id_t> rescue_nodes;
    int64_t min_distance = max(0.0, fragment_length_distr.mean() - rescued_alignment.sequence().size() - rescue_subgraph_stdevs * fragment_length_distr.std_dev());
    int64_t max_distance = fragment_length_distr.mean() + rescue_subgraph_stdevs * fragment_length_distr.std_dev();

    if (distance_index != nullptr) {
        subgraph_in_distance_range(*distance_index, aligned_read.path(), &cached_graph, min_distance, max_distance, rescue_nodes, rescue_forward);
    } else {
        old_distance_index->subgraph_in_range(aligned_read.path(), &cached_graph, min_distance, max_distance, rescue_nodes, rescue_forward);
    }

    if (rescue_nodes.size() == 0) {
        //If the rescue subgraph is empty
        return;
    }
    
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

    // Find all seeds in the subgraph and try to get a full-length extension.
    GaplessExtender::cluster_type seeds = this->seeds_in_subgraph(minimizers, rescue_nodes);
    if (seeds.size() > this->rescue_seed_limit) {
        return;
    }
    std::vector<GaplessExtension> extensions = this->extender.extend(seeds, rescued_alignment.sequence(), &cached_graph);

    // If we have a full-length extension, use it as the rescued alignment.
    if (GaplessExtender::full_length_extensions(extensions)) {
        this->extension_to_alignment(extensions.front(), rescued_alignment);
        return;
    }

    // Determine the best extension.
    size_t best = extensions.size();
    for (size_t i = 0; i < extensions.size(); i++) {
        if (best >= extensions.size() || extensions[i].score > extensions[best].score) {
            best = i;
        }
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
        
        size_t rescue_subgraph_bases = 0;
        for (auto& h : topological_order) {
            rescue_subgraph_bases += cached_graph.get_length(h);
        }
        if (rescue_subgraph_bases * rescued_alignment.sequence().size() > max_dozeu_cells) {
            if (!warned_about_rescue_size.test_and_set()) {
                cerr << "warning[vg::giraffe]: Refusing to perform too-large rescue alignment of "
                    << rescued_alignment.sequence().size() << " bp against "
                    << rescue_subgraph_bases << " bp ordered subgraph for read " << rescued_alignment.name()
                    << " which would use more than " << max_dozeu_cells
                    << " cells and might exhaust Dozeu's allocator; suppressing further warnings." << endl;
            }
            return; 
        }
    
        if (rescue_algorithm == rescue_dozeu) {
            size_t gap_limit = this->get_regular_aligner()->longest_detectable_gap(rescued_alignment);
            get_regular_aligner()->align_xdrop(rescued_alignment, cached_graph, topological_order,
                                               dozeu_seed, false, gap_limit);
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
        handlealgs::dagify(&split_graph, &dagified, rescued_alignment.sequence().size());

    size_t rescue_subgraph_bases = dagified.get_total_length();
    if (rescue_subgraph_bases * rescued_alignment.sequence().size() > max_dozeu_cells) {
        if (!warned_about_rescue_size.test_and_set()) {
            cerr << "warning[vg::giraffe]: Refusing to perform too-large rescue alignment of "
                << rescued_alignment.sequence().size() << " bp against "
                << rescue_subgraph_bases << " bp dagified subgraph for read " << rescued_alignment.name()
                << " which would use more than " << max_dozeu_cells
                << " cells and might exhaust Dozeu's allocator; suppressing further warnings." << endl;
        }
        return; 
    }
    
    // Align to the subgraph.
    // TODO: Map the seed to the dagified subgraph.
    if (this->rescue_algorithm == rescue_dozeu) {
        size_t gap_limit = this->get_regular_aligner()->longest_detectable_gap(rescued_alignment);
        get_regular_aligner()->align_xdrop(rescued_alignment, dagified, std::vector<MaximalExactMatch>(), false, gap_limit);
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
    
    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "Rescue result: " << pb2json(rescued_alignment) << endl;
        }
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
    int32_t score = aligner->score_contiguous_alignment(rescued_alignment);
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

    int64_t min_dist;
    if (distance_index != nullptr) {
        min_dist = distance_index->minimum_distance(get_id(pos1), get_is_rev(pos1), get_offset(pos1), get_id(pos2), get_is_rev(pos1), get_offset(pos2), false, &gbwt_graph);
    } else {
        min_dist = old_distance_index->min_distance(pos1, pos2);
    }
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
    // Get minimizers and their window agglomeration starts and lengths
    // Starts and lengths are all 0 if we are using syncmers.
    vector<tuple<gbwtgraph::DefaultMinimizerIndex::minimizer_type, size_t, size_t>> minimizers =
        this->minimizer_index.minimizer_regions(sequence);
    for (auto& m : minimizers) {
        double score = 0.0;
        auto hits = this->minimizer_index.count_and_find(get<0>(m));
        if (hits.first > 0) {
            if (hits.first <= this->hard_hit_cap) {
                score = base_score - std::log(hits.first);
            } else {
                score = 1.0;
            }
        }
        
        // Length of the match from this minimizer or syncmer
        int32_t match_length = (int32_t) minimizer_index.k();
        // Number of candidate kmers that this minimizer is minimal of
        int32_t candidate_count = this->minimizer_index.uses_syncmers() ? 1 : (int32_t) minimizer_index.w();
        
        auto& value = std::get<0>(m);
        size_t agglomeration_start = std::get<1>(m);
        size_t agglomeration_length = std::get<2>(m);
        if (this->minimizer_index.uses_syncmers()) {
            // The index says the start and length are 0. Really they should be where the k-mer is.
            // So start where the k-mer is on the forward strand
            agglomeration_start = value.is_reverse ? (value.offset - (match_length - 1)) : value.offset;
            // And run for the k-mer length
            agglomeration_length = match_length;
        }
        
        result.push_back({ value, agglomeration_start, agglomeration_length, hits.first, hits.second,
                            match_length, candidate_count, score });
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

    // We group all all occurrences of the same minimizer in the read together
    // and either take all of them (if the total number of hits is low enough)
    // or skip all of them. Such minimizers are expensive to process, because
    // they tend to have many hits and each hit in the graph is matched with
    // each occurrence in the read.
    bool took_last = false;
    size_t start = 0, limit = 0;
    size_t run_hits = 0;

    if (show_work) {
        #pragma omp critical (cerr)
        {
            std::cerr << log_name() << "All minimizers:" << std::endl;
            dump_debug_minimizers(minimizers, aln.sequence());
        }
    }

    // bit vector length of read to check for overlaps
    size_t num_minimizers = 0;
    size_t read_len = aln.sequence().size();
    std::vector<bool> read_bit_vector (read_len, false);

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

        // Find the next run of identical minimizers.
        if (i >= limit) {
            start = i; limit = i + 1;
            run_hits = minimizers[i].hits;
            for (size_t j = i + 1; j < minimizers.size() && minimizers[j].value.key == minimizers[i].value.key; j++) {
                limit++;
                run_hits += minimizers[j].hits;
            }
        }

        // Select the minimizer if it is informative enough or if the total score
        // of the selected minimizers is not high enough.
        const Minimizer& minimizer = minimizers[i];

        // minimizer information
        size_t min_start_index = minimizer.forward_offset();
        size_t min_len = minimizer.length;
        bool overlapping = false;
        if (this->exclude_overlapping_min) {
          if (read_bit_vector[min_start_index] == true || read_bit_vector[min_start_index + min_len] == true) {
            overlapping = true;
          }
        }

        if (minimizer.hits == 0) {
            // A minimizer with no hits can't go on.
            took_last = false;
            // We do not treat it as located for MAPQ capping purposes.
            if (this->track_provenance) {
                funnel.fail("any-hits", i);
            }
        } else if (  // passes reads
              ((minimizer.hits <= this->hit_cap) ||
              (run_hits <= this->hard_hit_cap && selected_score + minimizer.score <= target_score) ||
              (took_last && i > start)) &&
              (num_minimizers < this->max_unique_min) &&
              (overlapping == false)
            ) {

            // set minimizer overlap as a reads
            num_minimizers += 1;    // tracking number of minimizers selected
            if (this->exclude_overlapping_min) {
              for (size_t i = min_start_index; i < min_start_index + min_len; i++) {
                read_bit_vector[i] = true;
              }
            }

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
                // TODO: Don't use the seed payload anymore
                tuple<size_t, size_t, size_t, size_t, bool> chain_info (MIPayload::NO_VALUE, MIPayload::NO_VALUE, MIPayload::NO_VALUE, MIPayload::NO_VALUE, false);
                if (minimizer.occs[j].payload != MIPayload::NO_CODE) {
                    chain_info = MIPayload::decode(minimizer.occs[j].payload);
                }
                seeds.push_back({ hit, i, chain_info});
            }

            // Remember that we took this minimizer
            selected_score += minimizer.score;
            took_last = true;

            if (this->track_provenance) {
                // Record in the funnel that this minimizer gave rise to these seeds.
                funnel.pass("any-hits", i);
                funnel.pass("hard-hit-cap", i);
                funnel.pass("hit-cap||score-fraction", i, selected_score  / base_target_score);
                funnel.expand(i, minimizer.hits);
            }
        } else if (run_hits <= this->hard_hit_cap) {
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
            for (size_t i = 0; i < seeds.size(); i++) {
                // Find every seed's reference positions. This maps from path name to pairs of offset and orientation.
                auto offsets = algorithms::nearest_offsets_in_paths(this->path_graph, seeds[i].pos, 100);
                
                for (auto& true_pos : aln.refpos()) {
                    // For every annotated true position
                    for (auto& hit_pos : offsets[this->path_graph->get_path_handle(true_pos.name())]) {
                        // Look at all the hit positions on the path the read's true position is on.
                        if (abs((int64_t)hit_pos.first - (int64_t) true_pos.offset()) < 200) {
                            // Call this seed hit close enough to be correct
                            funnel.tag_correct(i);
                        }
                    }
                }
            }
        }
    }

    if (show_work) {
        #pragma omp critical (cerr)
        {
            std::cerr << log_name() << "Found " << seeds.size() << " seeds from "
                << (minimizers.size() - rejected_count) << " minimizers, rejected "
                << rejected_count << std::endl;
        }
    }

    return seeds;
}
std::vector<MinimizerMapper::OldSeed> MinimizerMapper::find_seeds_old(const std::vector<Minimizer>& minimizers, const Alignment& aln, Funnel& funnel) const {

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

    // We group all all occurrences of the same minimizer in the read together
    // and either take all of them (if the total number of hits is low enough)
    // or skip all of them. Such minimizers are expensive to process, because
    // they tend to have many hits and each hit in the graph is matched with
    // each occurrence in the read.
    bool took_last = false;
    size_t start = 0, limit = 0;
    size_t run_hits = 0;

    if (show_work) {
        #pragma omp critical (cerr)
        {
            std::cerr << log_name() << "All minimizers:" << std::endl;
            dump_debug_minimizers(minimizers, aln.sequence());
        }
    }

    // Select the minimizers we use for seeds.
    size_t rejected_count = 0;
    std::vector<OldSeed> seeds;
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

        // Find the next run of identical minimizers.
        if (i >= limit) {
            start = i; limit = i + 1;
            run_hits = minimizers[i].hits;
            for (size_t j = i + 1; j < minimizers.size() && minimizers[j].value.key == minimizers[i].value.key; j++) {
                limit++;
                run_hits += minimizers[j].hits;
            }
        }

        // Select the minimizer if it is informative enough or if the total score
        // of the selected minimizers is not high enough.
        const Minimizer& minimizer = minimizers[i];

        if (minimizer.hits == 0) {
            // A minimizer with no hits can't go on.
            took_last = false;
            // We do not treat it as located for MAPQ capping purposes.
            if (this->track_provenance) {
                funnel.fail("any-hits", i);
            }
        } else if (minimizer.hits <= this->hit_cap ||
            (run_hits <= this->hard_hit_cap && selected_score + minimizer.score <= target_score) ||
            (took_last && i > start)) {
            
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
                // TODO: Don't use the seed payload anymore
                tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info 
                    (false, MinimumDistanceIndex::MIPayload::NO_VALUE, MinimumDistanceIndex::MIPayload::NO_VALUE, false, MinimumDistanceIndex::MIPayload::NO_VALUE, MinimumDistanceIndex::MIPayload::NO_VALUE, MinimumDistanceIndex::MIPayload::NO_VALUE, MinimumDistanceIndex::MIPayload::NO_VALUE, false );
                if (minimizer.occs[j].payload != MinimumDistanceIndex::MIPayload::NO_CODE) {
                    chain_info = MinimumDistanceIndex::MIPayload::decode(minimizer.occs[j].payload);
                }
                seeds.push_back({ hit, i, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info), 
                    std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info) });

            }

            // Remember that we took this minimizer
            selected_score += minimizer.score;
            took_last = true;

            if (this->track_provenance) {
                // Record in the funnel that this minimizer gave rise to these seeds.
                funnel.pass("any-hits", i);
                funnel.pass("hard-hit-cap", i);
                funnel.pass("hit-cap||score-fraction", i, selected_score  / base_target_score);
                funnel.expand(i, minimizer.hits);
            }
        } else if (run_hits <= this->hard_hit_cap) {
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
            cerr << "error[vg::MinimizerMapper] Cannot use track_correctness with no graph with paths" << endl;
            exit(1);
        }

        if (aln.refpos_size() != 0) {
            // Compose path handle to position and strand data structure for the annotated reference positions
            algorithms::path_offset_collection_t ref_positions;
            for (auto& true_pos : aln.refpos()) {
               ref_positions[this->path_graph->get_path_handle(true_pos.name())].emplace_back(true_pos.offset(), true_pos.is_reverse());
            }
            // Make sure it is sorted
            algorithms::sort_path_offsets(ref_positions);
            
            // We need to know whether each seed is correct or not, and the
            // seeds probably have about one position each. So we look up the
            // seeds' positions in the sorted ref positions, and take: 
            // O(#seed positions * log(#ref positions)) 
            
            for (size_t i = 0; i < seeds.size(); i++) {
                // Find every seed's reference positions. This maps from path name to pairs of offset and orientation.
                algorithms::path_offset_collection_t seed_positions = algorithms::nearest_offsets_in_paths(this->path_graph, seeds[i].pos, 100);
                // Check against the annotated truth position set
                if (algorithms::intersect_path_offsets(ref_positions, seed_positions, 200)) {
                    // It's within range of a truth position.
                    // Call this seed hit close enough to be correct.
                    funnel.tag_correct(i);
                }
            }
        }
    }

    if (show_work) {
        #pragma omp critical (cerr)
        {
            std::cerr << log_name() << "Found " << seeds.size() << " seeds from "
                << (minimizers.size() - rejected_count) << " minimizers, rejected "
                << rejected_count << std::endl;
        }
    }

    return seeds;
}

void MinimizerMapper::old_annotate_with_minimizer_statistics(Alignment& target, const std::vector<Minimizer>& minimizers, const std::vector<OldSeed>& seeds, const Funnel& funnel) const {
    // Annotate with fraction covered by correct (and necessarily located) seed hits.
    
    // First make the set of minimizers that got correct seeds
    std::unordered_set<size_t> seeded;
    for (size_t i = 0; i < seeds.size(); i++) {
        // We need to get correctness out of the funnel, since we don't tag the seed or minimizer.
        // Correctness is assessed per seed, not per minimizer.
        // We know seed finding was always stage 1.
        if (funnel.was_correct(1, "seed", i)) {
            seeded.insert(seeds[i].source);
        }
    }
    
    // Then we make a table of all the ranges covered by correct minimizers
    std::vector<std::pair<size_t, size_t>> bounds;
    bounds.reserve(seeded.size());
    for(auto& minimizer_number : seeded) {
        // For each minimizer with correct seeds
        auto& minimizer = minimizers[minimizer_number];
        // Cover the minimizer k-mer itself.
        size_t start = minimizer.forward_offset();
        bounds.emplace_back(start, start + minimizer.length);
    }
    // Then we count the positions covered
    size_t covered_count = algorithms::count_covered(bounds);
    // And turn it into a fraction
    double covered_fraction = (double) covered_count / target.sequence().size();
    // And add the annotation
    set_annotation(target, "correct-minimizer-coverage", covered_fraction);
}
void MinimizerMapper::annotate_with_minimizer_statistics(Alignment& target, const std::vector<Minimizer>& minimizers, const std::vector<Seed>& seeds, const Funnel& funnel) const {
    // Annotate with fraction covered by correct (and necessarily located) seed hits.
    
    // First make the set of minimizers that got correct seeds
    std::unordered_set<size_t> seeded;
    for (size_t i = 0; i < seeds.size(); i++) {
        // We need to get correctness out of the funnel, since we don't tag the seed or minimizer.
        // Correctness is assessed per seed, not per minimizer.
        // We know seed finding was always stage 1.
        if (funnel.was_correct(1, "seed", i)) {
            seeded.insert(seeds[i].source);
        }
    }
    
    // Then we make a table of all the ranges covered by correct minimizers
    std::vector<std::pair<size_t, size_t>> bounds;
    bounds.reserve(seeded.size());
    for(auto& minimizer_number : seeded) {
        // For each minimizer with correct seeds
        auto& minimizer = minimizers[minimizer_number];
        // Cover the minimizer k-mer itself.
        size_t start = minimizer.forward_offset();
        bounds.emplace_back(start, start + minimizer.length);
    }
    // Then we count the positions covered
    size_t covered_count = algorithms::count_covered(bounds);
    // And turn it into a fraction
    double covered_fraction = (double) covered_count / target.sequence().size();
    // And add the annotation
    set_annotation(target, "correct-minimizer-coverage", covered_fraction);
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
        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "Minimizer " << seeds[hit_index].source << " is present in cluster " << i << endl;
            }
        }
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
void MinimizerMapper::score_cluster_old(Cluster& cluster, size_t i, const std::vector<Minimizer>& minimizers, const std::vector<OldSeed>& seeds, size_t seq_length, Funnel& funnel) const {

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
        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "Minimizer " << seeds[hit_index].source << " is present in cluster " << i << endl;
            }
        }
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

//#define debug_chaining

// For doing scores with backtracing, we use this type, which is a
// score and a number for the place we came from to get it.
// It can be maxed with std::max() as long as the destinations are all filled in right.
using traced_score_t = pair<int, size_t>;

/// Accessors for attributes of a type when used as a (possibly provenance-traced) DP table score.
template<typename S>
struct score_traits {
    // We don't actually define the interface here, because templates don't do any sort of inheritance.
};

/// This is how to use a traced_score_t as a DP score, and all the operations
/// you can do with it.
template<>
struct score_traits<traced_score_t> {
    
    /// What is the sentinel for an empty provenance?
    static size_t NOWHERE;
    
    /// What's the default value for an empty table cell?
    static traced_score_t UNSET;
    
    /// Pack a score together with its provenance
    static traced_score_t annotate(int points, size_t from) {
        return {points, from};
    }

    /// Accessor to get the score
    static int& score(traced_score_t& item) {
        return item.first;
    }
    
    /// Accessor to get the score, when const
    static const int& score(const traced_score_t& item) {
        return item.first;
    }
    
    /// Accessor to get the source
    static size_t& source(traced_score_t& item) {
        return item.second;
    }
    
    /// Accessor to get the source, when const
    static const size_t& source(const traced_score_t& item) {
        return item.second;
    }
    
    /// Max in a score from a DP table. If it wins, record provenance.
    static void max_in(traced_score_t& dest, const vector<traced_score_t>& options, size_t option_number) {
        auto& option = options[option_number];
        if (score(option) > score(dest) || source(dest) == NOWHERE) {
            // This is the new winner.
            score(dest) = score(option);
            source(dest) = option_number;
        }
    }
    
    /// Get a score from a table and record provenance in it.
    static traced_score_t score_from(const vector<traced_score_t>& options, size_t option_number) {
        traced_score_t got = options[option_number];
        source(got) = option_number;
        return got;
    }
    
    /// Add (or remove) points along a route to somewhere. Return a modified copy.
    static traced_score_t add_points(const traced_score_t& item, int adjustment) {
        return annotate(score(item) + adjustment, source(item));
    }
};
// Give constants a compilation unit
size_t score_traits<traced_score_t>::NOWHERE = numeric_limits<size_t>::max();
traced_score_t score_traits<traced_score_t>::UNSET = {0, score_traits<traced_score_t>::NOWHERE};


/// This is how to use an int as a DP score, and all the operations you can do
/// with it.
template<>
struct score_traits<int> {
    
    /// What is the sentinel for an empty provenance?
    static size_t NOWHERE;
    
    /// What's the default value for an empty table cell?
    static int UNSET;

    /// Pack a score together with its provenance
    static int annotate(int points, size_t from) {
        return points;
    }

    /// Accessor to get the score.
    static int& score(int& item) {
        return item;
    }
    
    /// Accessor to get the score, when const.
    static const int& score(const int& item) {
        return item;
    }
    
    /// Accessor to get the source, when const (always NOWHERE)
    static size_t source(const int& item) {
        return NOWHERE;
    }
    // Note that source can't be set.
    
    /// Max in a score from a DP table. If it wins, record provenance.
    static void max_in(int& dest, const vector<int>& options, size_t option_number) {
        dest = std::max(dest, options[option_number]);
    }
    
    /// Get a score from a table.
    static int score_from(const vector<int>& options, size_t option_number) {
        return options[option_number];
    }
    
    /// Add (or remove) points along a route to somewhere. Return a modified copy.
    static int add_points(const int& item, int adjustment) {
        return item + adjustment;
    }
};
// Give constants a compilation unit
size_t score_traits<int>::NOWHERE = numeric_limits<size_t>::max();
int score_traits<int>::UNSET = 0;

/// Print operator
static ostream& operator<<(ostream& out, const traced_score_t& value) {
    if (score_traits<traced_score_t>::source(value) == score_traits<traced_score_t>::NOWHERE) {
        return out << score_traits<traced_score_t>::score(value) << " from nowhere";
    }
    return out << score_traits<traced_score_t>::score(value) << " from " << score_traits<traced_score_t>::source(value);
}

size_t MinimizerMapper::get_graph_overlap(const GaplessExtension& left,
                                          const GaplessExtension& right,
                                          const HandleGraph* graph) {
                                          
    
    // We can't just use GaplessExtension::overlap() because we need to worry
    // about containment, or the left thing actually starting after the right
    // thing.
    
    // We also don't really care if they take different overall paths in the graph. We find overlap if a start node is shared.
    
    // This is the handle the right extension starts on
    handle_t right_starts_on = right.path.front();
    
    // Scan back along the left path until we find the handle the right path starts on 
    auto last_right_start_in_left = left.path.rbegin();
    while (last_right_start_in_left != left.path.rend() && *last_right_start_in_left != right_starts_on) {
        last_right_start_in_left++;
    }
    
    if (last_right_start_in_left != left.path.rend()) {
        // The right extension starts at some handle where the left extension also is.
        
        // Work out how many bases of the left extension are before that handle.
        // We start with all the bases on the path.
        size_t left_bases_before_shared_handle = 0;
        auto it = last_right_start_in_left;
        ++it;
        if (it != left.path.rend()) {
            // There are some left bases before the shared handle.
            while (it != left.path.rend()) {
                left_bases_before_shared_handle += graph->get_length(*it);
                ++it;
            }
            // Since we have some bases, we back out the offset.
            left_bases_before_shared_handle -= left.offset;
            
            // Work out how many bases of the left extension are on or after that handle
            size_t left_bases_remaining = left.length() - left_bases_before_shared_handle;
            
            // Work out how many bases into the handle the right extension starts
            size_t right_offset = right.offset;
            
            // Compute how many bases of the left extension the right extension
            // starts at or before and return that.
            if (left_bases_remaining <= right_offset) {
                // They don't actually overlap
                return 0;
            } else {
                // They overlap by however many bases the left extension goes into
                // this handle past where the right extension starts.
                return left_bases_remaining - right_offset;
            }
                
        } else {
            // The left extension actually starts on the shared handle where the right extension also starts. It might be completely before or completely after the right one actually.
            size_t left_past_end = left.offset + left.length();
            if (left_past_end <= right.offset) {
                // Actually there's no overlap, left will end before right starts.
                return 0;
            } else {
                // There's some overlap
                return left_past_end - right.offset;
            }
        }
    } else {
        // The handle the right extension starts on doesn't actually exist in the left path.
        // So the right extension either doesn't actually overlap, or else contains the left one.
        // Or else the right extension starts somewhere else and reads into the left extension's path.
        
        // This is the handle the left extension starts on
        handle_t left_starts_on = left.path.front();
        
        // So, find the earliest place the left extension's path could start in the right one
        auto first_left_start_in_right = right.path.begin();
        // And how many bases come before it along the left extension
        size_t right_bases_before_shared_handle = 0;
        while (first_left_start_in_right != right.path.end() && *first_left_start_in_right != left_starts_on) {
            right_bases_before_shared_handle += graph->get_length(*first_left_start_in_right);
            first_left_start_in_right++;
        }
        // Since we know the right path's start handle does not appear on the
        // left path, we know we went through at least one handle and we can
        // safely back out the offset.
        right_bases_before_shared_handle -= right.offset;
        
        if (first_left_start_in_right != right.path.end()) {
            // The left path actually does have its starting handle somewhere in the right path.
            
            // We know how many right path bases are before that handle.
            // Then we add to that overlap the distance along the handle to the
            // start of the left extension, and also the length of the left
            // extension.
            return right_bases_before_shared_handle + left.offset + left.length();
        } else {
            // Neither path has the start handle of the other.
            // We could maybe come up with some notion of overlap between points along the paths, but we can't do it for the extensions as a whole.
            // TODO: is this going to make the overall algorithm work well when these e.g. start/end on SNP nodes?
            return 0;
        }
    }
    
}

size_t MinimizerMapper::get_graph_distance(const GaplessExtension& left,
                                           const GaplessExtension& right,
                                           const SnarlDistanceIndex* distance_index,
                                           const HandleGraph* graph) {
    // Find where the left extension past-ends
    Position left_past_end = left.tail_position(*graph);
    // Find where the right one starts
    Position right_start = right.starting_position(*graph);
    // Get the oriented minimum distance from the index
    size_t distance = distance_index->minimum_distance(
        left_past_end.node_id(), left_past_end.is_reverse(), left_past_end.offset(),
        right_start.node_id(), right_start.is_reverse(), right_start.offset(),
        false, graph);
    // And return it
    return distance;
}

template<typename S>
S MinimizerMapper::chain_extension_group_dp(vector<S>& best_chain_score,
                                            const Alignment& aln,
                                            const vector<GaplessExtension>& extended_seeds,
                                            int gap_open_penalty,
                                            int gap_extend_penalty,
                                            const SnarlDistanceIndex* distance_index,
                                            const HandleGraph* graph) {
                    
    // Grab the traits into a short name so we can use the accessors concisely.
    using ST = score_traits<S>;

#ifdef debug_chaining
    cerr << "Chaining extension group of " << extended_seeds.size() << " extensions" << endl;
#endif
    
    
    // This is a collection of one or more non-full-length extended seeds.
    
    // We use a sweep line algorithm to find relevant points along the read: extension starts or ends.
    // This records the last base to be covered by the current sweep line.
    int64_t sweep_line = 0;
    // This records the first base not covered by the last sweep line.
    int64_t next_unswept = 0;
    
    // And we track the next unentered gapless extension
    size_t unentered = 0;
    
    // Extensions we are in are in this min-heap of past-end position and gapless extension number.
    using ending_at_t = pair<size_t, size_t>;
    priority_queue<ending_at_t, vector<ending_at_t>, std::greater<ending_at_t>> end_heap;
    
    // We track the best score for a chain reaching the position before this one and ending in a gap.
    // We never let it go below 0.
    // Will be 0 when there's no gap that can be open
    S best_gap_score = ST::UNSET;
    
    // We track the score for the best chain ending with each gapless extension
    best_chain_score.clear();
    best_chain_score.resize(extended_seeds.size(), ST::UNSET);
    
    // And we're after the best score overall that we can reach when an extension ends
    S best_past_ending_score_ever = ST::UNSET;
    
    // Overlaps are more complicated.
    // We need a heap of all the extensions for which we have seen the
    // start and that we can thus overlap.
    // We filter things at the top of the heap if their past-end positions
    // have occurred.
    // So we store pairs of score we get backtracking to the current
    // position, and past-end position for the thing we are backtracking
    // from.
    // We can just use the standard max-heap comparator.
    priority_queue<pair<S, size_t>> overlap_heap;
    
    // We encode the score relative to a counter that we increase by the
    // gap extend every base we go through, so we don't need to update and
    // re-sort the heap.
    int overlap_score_offset = 0;
    
    while(next_unswept <= aln.sequence().size()) {
        // We are processed through the position before next_unswept.
        
        // Find a place for sweep_line to go
        
        // Find the next seed start
        int64_t next_seed_start = numeric_limits<int64_t>::max();
        if (unentered < extended_seeds.size()) {
            next_seed_start = extended_seeds[unentered].read_interval.first;
        }
        
        // Find the next seed end
        int64_t next_seed_end = numeric_limits<int64_t>::max();
        if (!end_heap.empty()) {
            next_seed_end = end_heap.top().first;
        }
        
        // Whichever is closer between those points and the end, do that.
        sweep_line = min(min(next_seed_end, next_seed_start), (int64_t) aln.sequence().size());
        
        // So now we're only interested in things that happen at sweep_line.
#ifdef debug_chaining
        cerr << "Sweep to " << sweep_line << endl;
#endif
        
        // Compute the distance from the previous sweep line position
        // Make sure to account for next_unswept's semantics as the next unswept base.
        int sweep_distance = sweep_line - next_unswept + 1;
        
        // We need to track the score of the best thing that past-ended here
        S best_past_ending_score_here = ST::UNSET;
        
        while(!end_heap.empty() && end_heap.top().first == sweep_line) {
            // Find anything that past-ends here
            size_t past_ending = end_heap.top().second;
            
            // Mix it into the score
            ST::max_in(best_past_ending_score_here, best_chain_score, past_ending);
            
            // Remove it from the end-tracking heap
            end_heap.pop();
        }
#ifdef debug_chaining
        cerr << "Best score of an extension past-ending here: " << best_past_ending_score_here << endl;
#endif
        

        // Pick between that and the best score overall
        best_past_ending_score_ever = std::max(best_past_ending_score_ever, best_past_ending_score_here);
        
        if (sweep_line == aln.sequence().size()) {
            // We don't need to think about gaps or backtracking anymore since everything has ended
            break;
        }
        
        // Update the overlap score offset by removing some gap extends from it.
        overlap_score_offset += sweep_distance * gap_extend_penalty;
        
        // The best way to backtrack to here is whatever is on top of the heap, if anything, that doesn't past-end here.
        S best_overlap_score = ST::UNSET;
        while (!overlap_heap.empty()) {
            // While there is stuff on the heap
            if (overlap_heap.top().second <= sweep_line) {
                // We are already past this thing, so drop it
                overlap_heap.pop();
            } else {
                // This is at the top of the heap and we aren't past it
                // Decode and use its score offset if we only backtrack to here.
                best_overlap_score = ST::add_points(overlap_heap.top().first, overlap_score_offset);
                // Stop looking in the heap
                break;
            }
        }
#ifdef debug_chaining
        cerr << "Best score of overlapping back to here: " << best_overlap_score << endl;
#endif
        
        // The best way to end 1 before here in a gap is either:
        
        if (best_gap_score != ST::UNSET) {
            // Best way to end 1 before our last sweep line position with a gap, plus distance times gap extend penalty
            ST::score(best_gap_score) -= sweep_distance * gap_extend_penalty;
        }
        
        // Best way to end 1 before here with an actual extension, plus the gap open part of the gap open penalty.
        // (Will never be taken over an actual adjacency)
        best_gap_score = std::max(ST::UNSET, std::max(best_gap_score, ST::add_points(best_past_ending_score_here, -(gap_open_penalty - gap_extend_penalty))));
#ifdef debug_chaining
        cerr << "Best score here but in a gap: " << best_gap_score << endl;
#endif
        
        while (unentered < extended_seeds.size() && extended_seeds[unentered].read_interval.first == sweep_line) {
            // For each thing that starts here
            
            // We want to compute how much we should charge to come from
            // each place we could come from, and then come from the best
            // one.
            S overlap_option = best_overlap_score;
            S gap_option = best_gap_score;
            S here_option = best_past_ending_score_here;
            
            if (distance_index && graph) {
                // We have a distance index, so we might charge different
                // amounts for different transitions depending on the graph
                // distance between the extensions.
                
                // So we adjust the options we have bgased on where we are
                // coming from and where we are going to.
                
                // We could take:
                // An overlap from source(overlap_option):
                if (ST::source(overlap_option) != ST::NOWHERE) {
                    // See whether and how much they overlap in the graph
                    size_t graph_overlap = get_graph_overlap(extended_seeds[ST::source(overlap_option)], extended_seeds[unentered], graph);
                    // And also get the overlap in the read.
                    size_t read_overlap = extended_seeds[ST::source(overlap_option)].read_interval.second - extended_seeds[unentered].read_interval.first;
                    if (graph_overlap > read_overlap) {
                        // They overlap in the graph more than they do in the read.
                        // We don't want to charge for the overlap in the read, we just want to charge for the additional overlap in the graph.
                        // We have (read overlap - 1) extend penalties, and we want (graph overlap - read overlap - 1) of them.
                        // (graph overlap - read overlap - 1) - (read overlap - 1) = graph overlap - 2 * read overlap
                        // Which is how many penalties to charge
                        overlap_option = ST::add_points(overlap_option, -gap_extend_penalty * graph_overlap + gap_extend_penalty * 2 * read_overlap);
                        // TODO: We need to not double-count the matches right?
                        // When we go to synthesize an actual alignment we can cut out the double-matched area and match it only one place.
                    } else if (graph_overlap == read_overlap) {
                        // They overlap in the graph exactly as much as they do in the read.
                        // Back out the whole read gap
                        overlap_option = ST::add_points(overlap_option, gap_open_penalty + gap_extend_penalty * (read_overlap - 1));
                        // TODO: We need to not double-count the matches right?
                    } else if (graph_overlap > 0) {
                        // They overlap in the graph but less than they do in the read.
                        // Back out some of the read gap
                        overlap_option = ST::add_points(overlap_option, gap_extend_penalty * (read_overlap - graph_overlap));
                    } else {
                        // They don't overlap in the graph at all; they might abut or be separated or be unreachable.
                        
                        // Compute the graph distance if they don't overlap in the graph
                        size_t graph_distance = get_graph_distance(extended_seeds[ST::source(overlap_option)],
                                                                   extended_seeds[unentered],
                                                                   distance_index,
                                                                   graph);
                        
                        if (graph_distance == numeric_limits<size_t>::max()) {
                            // We can't actually get between the relevant graph areas to do the overlap.
                            // TODO: check if we can get between the points on either side of the overlapped area instead?
                            // Right now just say we can't make this connection.
                            overlap_option = ST::UNSET;
                        } else {
                            // We need to charge for this extra graph distance being deleted, in addition to what we charge for the existing overlap gap.
                            overlap_option = ST::add_points(overlap_option, -gap_extend_penalty * graph_distance);
                        }
                    }
                }
                
                // A gap in the read from ST::source(gap_option):
                if (ST::source(gap_option) != ST::NOWHERE) {
                    // Compute the graph distance
                    size_t graph_distance = get_graph_distance(extended_seeds[ST::source(gap_option)],
                                                               extended_seeds[unentered],
                                                               distance_index,
                                                               graph);
                    if (graph_distance == numeric_limits<size_t>::max()) {
                        // This is actually unreachable in the graph.
                        // So say we can't actually do this.
                        gap_option = ST::UNSET;
                    } else {
                        // Compare to the read distance
                        size_t read_distance = extended_seeds[unentered].read_interval.first - extended_seeds[ST::source(gap_option)].read_interval.second;
                        if (graph_distance == read_distance) {
                            // This is actually even length.
                            // Charge nothing for the gap by backing out the penalty.
                        } else if (graph_distance > read_distance) {
                            // Graph distance is longer.
                            // We have (read distance - 1) extend penalties, and we want (graph distance - read distance - 1) of them.
                            // (graph distance - read distance - 1) - (read distance - 1) = graph distance - 2 * read distance
                            // Which is how many penalties to charge
                            gap_option = ST::add_points(gap_option, -gap_extend_penalty * graph_distance + gap_extend_penalty * 2 * read_distance);
                        } else {
                            // Read distance is longer.
                            // Back out extend penalties from the part of
                            // the gap that the graph also has, leaving
                            // only the remaining surplus read gap.
                            gap_option = ST::add_points(gap_option, gap_extend_penalty * graph_distance);
                        }
                    }
                }
                
                // An adjacency in the read from ST::source(here_option):
                if (ST::source(here_option) != ST::NOWHERE) {
                    // Compute the graph distance
                    size_t graph_distance = get_graph_distance(extended_seeds[ST::source(here_option)],
                                                               extended_seeds[unentered],
                                                               distance_index,
                                                               graph);
                    if (graph_distance == numeric_limits<size_t>::max()) {
                        // This is actually unreachable in the graph.
                        // So say we can't actually do this.
                        here_option = ST::UNSET;
                    } else if (graph_distance > 0) { 
                        // There's a gap in the graph but not in the read.
                        // Charge for the gap
                        here_option = ST::add_points(here_option, -gap_open_penalty - gap_extend_penalty * (graph_distance - 1));
                    } else {
                        // These abut in the read and the graph, so do nothing
                    }
                }
                
            }
            
            // Compute its chain score, allowing us to just start here also.
            best_chain_score[unentered] = ST::add_points(
                std::max(std::max(ST::UNSET,
                                  overlap_option),
                         std::max(gap_option,
                                  here_option)),
                extended_seeds[unentered].score);
#ifdef debug_chaining
            cerr << "Best score of chain ending in extension " << unentered << ": " << best_chain_score[unentered] << endl;
#endif
            
            // Compute its backtrack-to-here score and add it to the backtracking heap
            // We want how far we would have had to have backtracked to be
            // able to preceed the base we are at now, where this thing
            // starts. That's a gap open for the last base, and then an
            // extend for each base before it.
            size_t extension_length = extended_seeds[unentered].read_interval.second - extended_seeds[unentered].read_interval.first - 1;
            S raw_overlap_score = ST::add_points(
                ST::score_from(best_chain_score, unentered),
                -gap_open_penalty - gap_extend_penalty * extension_length);
            S encoded_overlap_score = ST::add_points(raw_overlap_score, -overlap_score_offset);
            
            // Stick it in the overlap heap
            overlap_heap.emplace(encoded_overlap_score, extended_seeds[unentered].read_interval.second);
#ifdef debug_chaining
            cerr << "Overlap encoded as " << encoded_overlap_score << " available until " << extended_seeds[unentered].read_interval.second << endl;
#endif
            
            // Add it to the end finding heap
            end_heap.emplace(extended_seeds[unentered].read_interval.second, unentered);
            
            // Advance and check the next thing to start
            unentered++;
        }
        
        // Move next_unswept to sweep_line.
        // We need to add 1 since next_unswept is the next *un*included base
        next_unswept = sweep_line + 1;
    }
    
    return best_past_ending_score_ever; 

    
}

int MinimizerMapper::score_extension_group(const Alignment& aln, const vector<GaplessExtension>& extended_seeds,
                                           int gap_open_penalty, int gap_extend_penalty) {
    
    if (extended_seeds.empty()) {
        // TODO: We should never see an empty group of extensions
        return 0;
    } else if (GaplessExtender::full_length_extensions(extended_seeds)) {
        // These are full-length matches. We already have the score.
        return extended_seeds.front().score;
    } else if (aln.sequence().size() == 0) {
        // No score here
        return extended_seeds.front().score;
    } else {
        // Do the DP but without the traceback
        vector<int> best_chain_score;
        return chain_extension_group_dp(best_chain_score,
                                        aln,
                                        extended_seeds,
                                        gap_open_penalty,
                                        gap_extend_penalty);
    }
}

pair<int, vector<size_t>> MinimizerMapper::chain_extension_group(const Alignment& aln,
                                                                 const vector<GaplessExtension>& extended_seeds,
                                                                 int gap_open_penalty,
                                                                 int gap_extend_penalty,
                                                                 const SnarlDistanceIndex* distance_index,
                                                                 const HandleGraph* graph) {
                                                                 
    if (extended_seeds.empty()) {
        // TODO: We should never see an empty group of extensions
        return std::make_pair(0, vector<size_t>());
    } else if (GaplessExtender::full_length_extensions(extended_seeds)) {
        // These are full-length matches. We already have the score.
        return std::make_pair(extended_seeds.front().score, vector<size_t>{0});
    } else if (aln.sequence().size() == 0) {
        // No score here
        return std::make_pair(extended_seeds.front().score, vector<size_t>());
    } else {
        // We actually need to do DP
        vector<traced_score_t> best_chain_score;
        traced_score_t best_past_ending_score_ever = chain_extension_group_dp(best_chain_score,
                                                                              aln,
                                                                              extended_seeds,
                                                                              gap_open_penalty,
                                                                              gap_extend_penalty,
                                                                              distance_index,
                                                                              graph);
        
        // Now we need to trace back.
        vector<size_t> traceback;
        size_t here = score_traits<traced_score_t>::source(best_past_ending_score_ever);
#ifdef debug_chaining
        cerr << "Chain ends at #" << here << " at " << extended_seeds[here].read_interval.first << "-" << extended_seeds[here].read_interval.second << " with score " << best_past_ending_score_ever << endl;
#endif
        while(here != score_traits<traced_score_t>::NOWHERE) {
            traceback.push_back(here);
#ifdef debug_chaining
            cerr << "Which gets score " << best_chain_score[here] << endl;
#endif
            here = score_traits<traced_score_t>::source(best_chain_score[here]);
#ifdef debug_chaining
            if (here != score_traits<traced_score_t>::NOWHERE) {
                cerr << "And comes after #" << here << " at " << extended_seeds[here].read_interval.first << "-" << extended_seeds[here].read_interval.second << endl;
            } else {
                cerr << "And is first" << endl;
            }
#endif
        }
        // Flip it around front-ways
        std::reverse(traceback.begin(), traceback.end());
        
#ifdef debug_chaining
        cerr << "Best score of chain overall: " << best_past_ending_score_ever << endl;
#endif
        return std::make_pair(score_traits<traced_score_t>::score(best_past_ending_score_ever), traceback);
    }
}

std::pair<std::vector<int>, std::vector<std::vector<size_t>>> MinimizerMapper::chain_extensions(
    const std::vector<std::vector<GaplessExtension>>& extensions, const Alignment& aln, Funnel& funnel) const {

    std::pair<std::vector<int>, std::vector<std::vector<size_t>>> result;
    result.first.reserve(extensions.size());
    result.second.reserve(extensions.size());

    // Extension chaining substage.
    if (this->track_provenance) {
        funnel.substage("chain");
    }

    // We now estimate the best possible alignment score for each cluster.
    for (size_t i = 0; i < extensions.size(); i++) {
        
        if (this->track_provenance) {
            funnel.producing_output(i);
        }
        
        // Do the chaining and move the answers
        auto chained = chain_extension_group(aln, extensions[i], get_regular_aligner()->gap_open, get_regular_aligner()->gap_extension, distance_index, &gbwt_graph);
        result.first.emplace_back(std::move(chained.first));
        result.second.emplace_back(std::move(chained.second));
        
        // Record the score with the funnel.
        if (this->track_provenance) {
            funnel.score(i, result.first.back());
            funnel.produced_output();
        }
    }

    return result;
}

// TODO: Combine the two score_extensions overloads into one template when we get constexpr if.

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

// (value, cost)
typedef std::pair<uint32_t, int32_t> pareto_point;

void find_pareto_frontier(std::vector<pareto_point>& v) {
    if(v.empty()) {
        return;
    }
    std::sort(v.begin(), v.end(), [](pareto_point a, pareto_point b) {
        return (a.second < b.second || (a.second == b.second && a.first > b.first));
    });
    size_t tail = 1;
    for (size_t i = 1; i < v.size(); i++) {
        if (v[i].first <= v[tail - 1].first) {
            continue;
        }
        v[tail] = v[i];
        tail++;
    }
    v.resize(tail);
    std::sort(v.begin(), v.end());
}

// Positive gap penalty if there is a gap.
int32_t gap_penalty(size_t length, const Aligner* aligner) {
    return (length == 0 ? 0 : aligner->gap_open + (length - 1) * aligner->gap_extension);
}

// Positive penalty for a number of mismatches.
int32_t mismatch_penalty(size_t n, const Aligner* aligner) {
    return n * (aligner->match + aligner->mismatch);
}

// Positive gap penalty, assuming that there is always a gap.
int32_t gap_penalty(size_t start, size_t limit, const Aligner* aligner) {
    return (start >= limit ? aligner->gap_open : aligner->gap_open + (limit - start - 1) * aligner->gap_extension);
}

// Positive flank penalty based on taking a gap to the end or to the Pareto frontier.
int32_t flank_penalty(size_t length, const std::vector<pareto_point>& frontier, const Aligner* aligner) {
    int32_t result = gap_penalty(length, aligner);
    for (size_t i = 0; i < frontier.size(); i++) {
        int32_t candidate = frontier[i].second + gap_penalty(frontier[i].first, length, aligner);
        result = std::min(result, candidate);
        if (frontier[i].first >= length) {
            break;
        }
    }
    return result;
}

void MinimizerMapper::find_optimal_tail_alignments(const Alignment& aln, const vector<GaplessExtension>& extended_seeds, LazyRNG& rng, Alignment& best, Alignment& second_best) const {

    // This assumes that full-length extensions have the highest scores.
    // We want to align at least two extensions and at least one
    // partial extension. However, we do not want to align more than one
    // partial extension, unless the score is very close to the best
    // extension or the extension looks very promising.
    size_t min_extensions = 1;
    for (const GaplessExtension& extension : extended_seeds) {
        if (extension.full()) {
            min_extensions++;
        }
    }
    if (min_extensions < 2) {
        min_extensions = 2;
    }

    /*
      (length, penalty) pairs sorted by length. Pareto frontiers for the
      number of bp we can align at each end and the corresponding alignment
      score penalty. We use three types of points:
      1. A gap from the start/end of the read to the start/end of the
         extension followed by the entire extension.
      2. A gap from the start/end of the read to the start/end of the
         extension followed by the extension until the first mismatch.
      3. An all-windows-length - 1 bp exact match at the start/end of the read.
    */
    const Aligner* aligner = this->get_regular_aligner();
    std::vector<pareto_point> left_frontier, right_frontier;
    {
        size_t seq_len = aln.sequence().length();
        for (const GaplessExtension& extension : extended_seeds) {
            if (extension.full()) {
                continue;
            }
            int32_t left_penalty = gap_penalty(extension.read_interval.first, aligner);
            int32_t mid_penalty = mismatch_penalty(extension.mismatches(), aligner);
            int32_t right_penalty = gap_penalty(seq_len - extension.read_interval.second, aligner);
            left_frontier.push_back(pareto_point(extension.read_interval.second, mid_penalty + left_penalty));
            right_frontier.push_back(pareto_point(seq_len - extension.read_interval.first, mid_penalty + right_penalty));
            if (extension.mismatches() > 0) {
                left_frontier.push_back(pareto_point(extension.mismatch_positions.front(), left_penalty));
                right_frontier.push_back(pareto_point(seq_len - extension.mismatch_positions.back() - 1, right_penalty));
            }
        }
        size_t window_length = this->minimizer_index.uses_syncmers() ? this->minimizer_index.k() : (this->minimizer_index.k() + this->minimizer_index.w() - 1);
        left_frontier.push_back(pareto_point(window_length - 1, 0));
        right_frontier.push_back(pareto_point(window_length - 1, 0));
    }
    find_pareto_frontier(left_frontier);
    find_pareto_frontier(right_frontier);

    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "Trying to find " << min_extensions << " tail alignments for "
                << extended_seeds.size() << " extended seeds" << endl;
        }
    }
    
    // We will keep the winning alignment here, in pieces
    Path winning_left;
    Path winning_middle;
    Path winning_right;
    int32_t winning_score = 0;

    Path second_left;
    Path second_middle;
    Path second_right;
    int32_t second_score = 0;
    
    // Handle each extension in the set
    bool partial_extension_aligned = false;
    int32_t threshold = -1;
    process_until_threshold_a<GaplessExtension, double>(extended_seeds,
        [&](size_t extended_seed_num) -> double {
            return static_cast<double>(extended_seeds[extended_seed_num].score);
        }, extension_score_threshold, min_extensions, max_local_extensions, rng, [&](size_t extended_seed_num) -> bool {
       
            // This extended seed looks good enough.
            const GaplessExtension& extension = extended_seeds[extended_seed_num];

            // Extensions with score at most this will not be aligned,
            // unless we do not have enough alignments.
            if (threshold < 0) {
                threshold = extension.score - extension_score_threshold;
            }

            // Identify the special case: We already have aligned a partial
            // extension and the current score is too far below the best
            // extension. We do not want to align further partial extensions,
            // unless they look very promising.
            // The estimate is based on taking a gap to read end or to another
            // extension on the Pareto frontier, for both ends.
            if (!extension.full()) {
                if (partial_extension_aligned && extension.score <= threshold) {
                    int32_t score_estimate = aln.sequence().length() * aligner->match + 2 * aligner->full_length_bonus -
                        mismatch_penalty(extension.mismatches(), aligner);
                    if (!extension.left_full) {
                        score_estimate -= flank_penalty(extension.read_interval.first, left_frontier, aligner);
                    }
                    if (!extension.right_full) {
                        score_estimate -= flank_penalty(aln.sequence().length() - extension.read_interval.second,
                            right_frontier, aligner);
                    }
                    if (score_estimate <= winning_score) {
                        return true;
                    }
                }
                partial_extension_aligned = true;
            }
            
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
            
            if (!extension.left_full) {
                // There is a left tail
                
                // Have scratch for the longest detectable gap
                size_t longest_detectable_gap;
    
                // Get the forest of all left tail placements
                auto forest = get_tail_forest(extension, aln.sequence().size(), true, &longest_detectable_gap);
           
                // Grab the part of the read sequence that comes before the extension
                string before_sequence = aln.sequence().substr(0, extension.read_interval.first);
                
                // Do right-pinned alignment
                left_tail_result = std::move(get_best_alignment_against_any_tree(forest, before_sequence,
                    extension.starting_position(gbwt_graph), false, longest_detectable_gap, rng));
            }
            
            if (!extension.right_full) {
                // There is a right tail
                
                // Have scratch for the longest detectable gap
                size_t longest_detectable_gap;
                
                // Get the forest of all right tail placements
                auto forest = get_tail_forest(extension, aln.sequence().size(), false, &longest_detectable_gap);
            
                // Find the sequence
                string trailing_sequence = aln.sequence().substr(extension.read_interval.second);
        
                // Do left-pinned alignment
                right_tail_result = std::move(get_best_alignment_against_any_tree(forest, trailing_sequence,
                    extension.tail_position(gbwt_graph), true, longest_detectable_gap, rng));
            }
            
            // Compute total score
            int32_t total_score = extension.score + left_tail_result.second + right_tail_result.second;
            
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Extended seed " << extended_seed_num << " has left tail of "
                        << extension.read_interval.first << "bp and right tail of "
                        << (aln.sequence().size() - extension.read_interval.second)
                        << "bp for total score " << total_score << endl;
                }
            }

            // Get the node ids of the beginning and end of each alignment
            id_t winning_start = winning_score == 0 ? 0 : (winning_left.mapping_size() == 0
                                          ? winning_middle.mapping(0).position().node_id()
                                          : winning_left.mapping(0).position().node_id());
            id_t current_start = left_tail_result.first.mapping_size() == 0
                                     ? gbwt_graph.get_id(extension.path.front())
                                     : left_tail_result.first.mapping(0).position().node_id();
            id_t winning_end = winning_score == 0 ? 0 : (winning_right.mapping_size() == 0
                                  ? winning_middle.mapping(winning_middle.mapping_size() - 1).position().node_id()
                                  : winning_right.mapping(winning_right.mapping_size()-1).position().node_id());
            id_t current_end = right_tail_result.first.mapping_size() == 0
                                ? gbwt_graph.get_id(extension.path.back())
                                : right_tail_result.first.mapping(right_tail_result.first.mapping_size()-1).position().node_id();

            // Is this left tail different from the currently winning left tail?
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
                winning_middle = extension.to_path(gbwt_graph, aln.sequence());
                winning_right = std::move(right_tail_result.first);

            } else if ((total_score > second_score || second_score == 0) && different_left && different_right) {
                // This is the new second best alignment seen so far and it is 
                // different from the best alignment.
                
                // Save the score
                second_score = total_score;
                // And the path parts
                second_left = std::move(left_tail_result.first);
                second_middle = extension.to_path(gbwt_graph, aln.sequence());
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

//-----------------------------------------------------------------------------

pair<Path, size_t> MinimizerMapper::get_best_alignment_against_any_tree(const vector<TreeSubgraph>& trees,
    const string& sequence, const Position& default_position, bool pin_left, size_t longest_detectable_gap, LazyRNG& rng) const {
   
    // We want the best alignment, to the base graph, done against any target path
    Path best_path;
    // And its score
    int32_t best_score = 0;
    
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
        
        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "First best alignment: " << pb2json(best_path) << " score " << best_score << endl;
            }
        }
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
            
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Align " << pb2json(current_alignment) << " pinned left" << endl;
                }
            }

#ifdef debug_dump_graph
            cerr << "Vs graph:" << endl;
            subgraph.for_each_handle([&](const handle_t& here) {
                cerr << subgraph.get_id(here) << " (" << subgraph.get_sequence(here) << "): " << endl;
                subgraph.follow_edges(here, true, [&](const handle_t& there) {
                    cerr << "\t" << subgraph.get_id(there) << " (" << subgraph.get_sequence(there) << ") ->" << endl;
                });
                subgraph.follow_edges(here, false, [&](const handle_t& there) {
                    cerr << "\t-> " << subgraph.get_id(there) << " (" << subgraph.get_sequence(there) << ")" << endl;
                });
            });
#endif

            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Limit gap length to " << longest_detectable_gap << " bp" << endl;
                }
            }
            
            size_t tail_subgraph_bases = subgraph.get_total_length();
            if (tail_subgraph_bases * sequence.size() > max_dozeu_cells) {
                if (!warned_about_tail_size.test_and_set()) {
                    cerr << "warning[vg::giraffe]: Refusing to perform too-large tail alignment of "
                        << sequence.size() << " bp against "
                        << tail_subgraph_bases << " bp tree which would use more than " << max_dozeu_cells
                        << " cells and might exhaust Dozeu's allocator; suppressing further warnings." << endl;
                }
            } else {
                // X-drop align, accounting for full length bonus.
                // We *always* do left-pinned alignment internally, since that's the shape of trees we get.
                // Make sure to pass through the gap length limit so we don't just get the default.
                get_regular_aligner()->align_pinned(current_alignment, subgraph, true, true, longest_detectable_gap);
            }
            
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "\tScore: " << current_alignment.score() << endl;
                }
            }
            
            if (current_alignment.path().mapping_size() > 0 && deterministic_beats(current_alignment.score(), best_score, rng)) {
                // This is a new best alignment, and it is nonempty.
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
                
                if (show_work) {
                    #pragma omp critical (cerr)
                    {
                        cerr << log_name() << "New best alignment is "
                            << pb2json(best_path) << " score " << best_score << endl;
                    }
                }
            }
        }
    }

    return make_pair(best_path, best_score);
}

vector<TreeSubgraph> MinimizerMapper::get_tail_forest(const GaplessExtension& extended_seed,
    size_t read_length, bool left_tails, size_t* longest_detectable_gap) const {

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
    
    // Make sure we have a place to store the longest detectable gap
    size_t gap_limit;
    if (!longest_detectable_gap) {
        longest_detectable_gap = &gap_limit;
    }
    
    // Work it out because we need it for the limit of our search distance
    *longest_detectable_gap = get_regular_aligner()->longest_detectable_gap(read_length, tail_length);

#ifdef debug
    cerr << "Tail length: " << tail_length << " Read length: " << read_length << " Longest detectable gap: " << *longest_detectable_gap << endl;
#endif
    
    // How long should we search? It should be the longest detectable gap plus the remaining sequence.
    size_t search_limit = *longest_detectable_gap + tail_length;
    
#ifdef debug
    cerr << "Search limit: now " << search_limit << endl;
#endif

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
    cerr << "DFS starting at offset " << from_offset << " on node "
        << gbwt_graph.get_id(from_handle) << " " << gbwt_graph.get_is_reverse(from_handle) << " of length "
        << gbwt_graph.get_length(from_handle) << " leaving " << distance_to_node_end << " bp" << endl;
#endif


    // Have a recursive function that does the DFS. We fire the enter and exit
    // callbacks, and the user can keep their own stack.
    function<void(const gbwt::SearchState&, size_t, bool, bool)> recursive_dfs = [&](const gbwt::SearchState& here_state,
        size_t used_distance, bool hide_root, bool root_measured_already) {
        
        handle_t here_handle = gbwt_graph.node_to_handle(here_state.node);
        
        if (!hide_root) {
            // Enter this handle if there are any bases on it to visit
            
#ifdef debug
            cerr << "Enter handle " << gbwt_graph.get_id(here_handle) << " " << gbwt_graph.get_is_reverse(here_handle) << endl;
#endif
            
            enter_handle(here_handle);
        }
        
        if (!root_measured_already) {
            // Up the used distance with our length
            used_distance += gbwt_graph.get_length(here_handle);

#ifdef debug
            cerr << "Node was " << gbwt_graph.get_length(here_handle) << " bp; Used " << used_distance << "/" << walk_distance << " bp distance" << endl;
#endif
        } else {
#ifdef debug
            cerr << "Node was already measured; Used " << used_distance << "/" << walk_distance << " bp distance" << endl;
#endif
        }
        
        if (used_distance < walk_distance) {
            // If we haven't used up all our distance yet
            
            gbwt_graph.follow_paths(here_state, [&](const gbwt::SearchState& there_state) -> bool {
                // For each next state
                
                // Otherwise, do it with the new distance value.
                // Don't hide the root on any child subtrees; only the top root can need hiding.
                recursive_dfs(there_state, used_distance, false, false);
                
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
    // Make sure we don't count the length of the root node inside the DFS,
    // since we are already feeding it in.
    recursive_dfs(start_state, distance_to_node_end, distance_to_node_end == 0, true);

}

//-----------------------------------------------------------------------------

double MinimizerMapper::score_alignment_pair(Alignment& aln1, Alignment& aln2, int64_t fragment_distance) {
    //Score a pair of alignments

    double dev = fragment_distance - fragment_length_distr.mean();
    double fragment_length_log_likelihood = (-dev * dev / (2.0 * fragment_length_distr.std_dev() * fragment_length_distr.std_dev()))/ get_aligner()->log_base;
    double score = aln1.score() + aln2.score() +fragment_length_log_likelihood ;

    //Don't let the fragment length log likelihood bring score down below the score of the best alignment
    double worse_score = std::min(aln1.score(), aln2.score());

    return std::max(score, worse_score);;
}

double MinimizerMapper::distance_to_annotation(int64_t distance) const {
    // We use numeric_limits<int64_t>::max() to represent no distance. But that
    // can't convert to double (which rounds up) and then safely back to int64.
    // We also aren't allowed inf or nan in a Protobuf double Value. So change
    // the sentinel value to 0 which is probably not a fragment length.
    if (distance == numeric_limits<int64_t>::max()) {
        distance = 0;
    }
    
    // Make sure we can't generate any >64 bit integers in the double cast by
    // clamping to the doubles that are also integers.
    static_assert(DBL_MANT_DIG <= 64, "We assume doubles have <64 bits of mantissa");
    double max_int_double = (double)((int64_t)1 << DBL_MANT_DIG);
    return max(min((double) distance, max_int_double), -max_int_double);
}

}


