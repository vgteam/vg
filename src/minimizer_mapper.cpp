/**
 * \file minimizer_mapper.cpp
 * Defines the code for the minimizer-and-GBWT-based mapper.
 */

#include "minimizer_mapper.hpp"
#include "annotation.hpp"
#include "path_subgraph.hpp"
#include "multipath_alignment.hpp"

#include <chrono>
#include <iostream>

namespace vg {

using namespace std;

MinimizerMapper::MinimizerMapper(const xg::XG* xg_index, const gbwt::GBWT* gbwt_index, const MinimizerIndex* minimizer_index,
    SnarlManager* snarl_manager, DistanceIndex* distance_index) :
    xg_index(xg_index), gbwt_index(gbwt_index), minimizer_index(minimizer_index),
    snarl_manager(snarl_manager), distance_index(distance_index), gbwt_graph(*gbwt_index, *xg_index),
    extender(gbwt_graph) {
    
    // Nothing to do!
}

void MinimizerMapper::map(Alignment& aln, AlignmentEmitter& alignment_emitter) {
    // For each input alignment
        
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
        
    // We will find all the seed hits
    vector<pos_t> seeds;
    
    // This will hold all the minimizers in the query
    vector<MinimizerIndex::minimizer_type> minimizers;
    // And either way this will map from seed to minimizer that generated it
    vector<size_t> seed_to_source;
    
    // Find minimizers in the query
    minimizers = minimizer_index->minimizers(aln.sequence());

    size_t rejected_count = 0;
    
    for (size_t i = 0; i < minimizers.size(); i++) {
        // For each minimizer
        if (hit_cap == 0 || minimizer_index->count(minimizers[i].first) <= hit_cap) {
            // The minimizer is infrequent enough to be informative, so feed it into clustering
            
            // Locate it in the graph
            for (auto& hit : minimizer_index->find(minimizers[i].first)) {
                // For each position, remember it and what minimizer it came from
                seeds.push_back(hit);
                seed_to_source.push_back(i);
            }
        } else {
            // The minimizer is too frequent
            rejected_count++;
        }
    }

#ifdef debug
    cerr << "Read " << aln.name() << ": " << aln.sequence() << endl;
    cerr << "Found " << seeds.size() << " seeds from " << (minimizers.size() - rejected_count) << " minimizers, rejected " << rejected_count << endl;
#endif
        
    // Cluster the seeds. Get sets of input seed indexes that go together.
    vector<hash_set<size_t>> clusters = clusterer.cluster_seeds(seeds, distance_limit, *snarl_manager, *distance_index);
    
    // Compute the covered portion of the read represented by each cluster.
    // TODO: Put this and sorting into the clusterer to deduplicate with vg cluster.
    vector<double> read_coverage_by_cluster;
    for (auto& cluster : clusters) {
        // We set bits in here to true when query anchors cover them
        vector<bool> covered(aln.sequence().size());
        // We use this to convert iterators to indexes
        auto start = aln.sequence().begin();
        
        for (auto& hit_index : cluster) {
            // For each hit in the cluster, work out what anchor sequence it is from.
            size_t source_index = seed_to_source.at(hit_index);
            
            for (size_t i = minimizers[source_index].second; i < minimizers[source_index].second + minimizer_index->k(); i++) {
                // Set all the bits in read space for that minimizer.
                // Each minimizr is a length-k exact match starting at a position
                covered[i] = true;
            }
        }
        
        // Count up the covered positions
        size_t covered_count = 0;
        for (auto bit : covered) {
            covered_count += bit;
        }
        
        // Turn that into a fraction
        read_coverage_by_cluster.push_back(covered_count / (double) covered.size());
    }

#ifdef debug
    cerr << "Found " << clusters.size() << " clusters" << endl;
#endif
    
    // Make a vector of cluster indexes to sort
    vector<size_t> cluster_indexes_in_order;
    for (size_t i = 0; i < clusters.size(); i++) {
        cluster_indexes_in_order.push_back(i);
    }

    // Put the most covering cluster's index first
    std::sort(cluster_indexes_in_order.begin(), cluster_indexes_in_order.end(), [&](const size_t& a, const size_t& b) -> bool {
        // Return true if a must come before b, and false otherwise
        return read_coverage_by_cluster.at(a) > read_coverage_by_cluster.at(b);
    });
    
    // We will fill this with the output alignments (primary and secondaries) in score order.
    vector<Alignment> aligned;
    aligned.reserve(cluster_indexes_in_order.size());
    
    // Annotate the original read with metadata before copying
    if (!sample_name.empty()) {
        aln.set_sample_name(sample_name);
    }
    if (!read_group.empty()) {
        aln.set_read_group(read_group);
    }
    
    for (size_t i = 0; i < max(min(max_alignments, cluster_indexes_in_order.size()), (size_t)1); i++) {
        // For each output alignment we will produce (always at least 1,
        // and possibly up to our alignment limit or the cluster count)
        
        // Produce an output Alignment
        aligned.emplace_back(aln);
        Alignment& out = aligned.back();
        // Clear any old refpos annotation
        out.clear_refpos();
        
        if (i < clusters.size()) {
            // We have a cluster; it actually mapped

#ifdef debug
            cerr << "Cluster " << cluster_indexes_in_order[i] << " rank " << i << ": " << endl;
#endif
        
            // For each cluster
            hash_set<size_t>& cluster = clusters[cluster_indexes_in_order[i]];
            
            // Pack the seeds into (read position, graph position) pairs.
            vector<pair<size_t, pos_t>> seed_matchings;
            seed_matchings.reserve(cluster.size());
            for (auto& seed_index : cluster) {
                // For each seed in the cluster, generate its matching pair
                seed_matchings.emplace_back(minimizers[seed_to_source[seed_index]].second, seeds[seed_index]);
#ifdef debug
                cerr << "Seed read:" << minimizers[seed_to_source[seed_index]].second << " = " << seeds[seed_index]
                    << " from minimizer " << seed_to_source[seed_index] << "(" << minimizer_index->count(minimizers[seed_to_source[seed_index]].first) << ")" << endl;
#endif
            }
            
            // Extend seed hits in the cluster into a real alignment path and mismatch count.
            std::pair<Path, size_t> extended = extender.extend_seeds(seed_matchings, aln.sequence());
            auto& path = extended.first;
            auto& mismatch_count = extended.second;

#ifdef debug
            cerr << "Produced path with " << path.mapping_size() << " mappings and " << mismatch_count << " mismatches" << endl;
#endif

            if (path.mapping_size() != 0) {
                // We have a mapping
                
                // Compute a score based on the sequence length and mismatch count.
                // Alignments will only contain matches and mismatches.
                int alignment_score = default_match * (aln.sequence().size() - mismatch_count) - default_mismatch * extended.second;
                
                if (path.mapping().begin()->edit_size() != 0 && edit_is_match(*path.mapping().begin()->edit().begin())) {
                    // Apply left full length bonus based on the first edit
                    alignment_score += default_full_length_bonus;
                }
                if (path.mapping().rbegin()->edit_size() != 0 && edit_is_match(*path.mapping().rbegin()->edit().rbegin())) {
                    // Apply right full length bonus based on the last edit
                    alignment_score += default_full_length_bonus;
                }
               
                // Compute identity from mismatch count.
                double identity = aln.sequence().size() == 0 ? 0.0 : (aln.sequence().size() - mismatch_count) / (double) aln.sequence().size();
                
                // Fill in the extension info
                *out.mutable_path() = path;
                out.set_score(alignment_score);
                out.set_identity(identity);
                
                // Read mapped successfully!
                continue;
            } else {
                // We need to generate some sub-full-length, maybe-extended seeds.
                // Call back into the extender and get the unambiguous perfect match extensions of the seeds in the cluster.
                vector<pair<Path, size_t>> extended_seeds = extender.maximal_extensions(seed_matchings, aln.sequence());

#ifdef debug
                cerr << "Trying again to chain " << extended_seeds.size() << " extended seeds" << endl;
#endif

                // TODO: split extended seeds when they overlap in the read, so
                // they either don't overlap or completely overlap (and are
                // thus mutually exclusive).

                // Then we need to find all the haplotypes between each pair of seeds that can connect.

                // Sort the extended seeds by read start position.
                // We won't be able to match them back to the minimizers anymore but we won't need to.
                std::sort(extended_seeds.begin(), extended_seeds.end(), [&](const pair<Path, size_t>& a, const pair<Path, size_t>& b) -> bool {
                    // Return true if a needs to come before b.
                    // This will happen if a is earlier in the read than b.
                    return a.second < b.second;
                });

                // Find the paths between pairs of extended seeds that agree with haplotypes.
                // We don't actually need the read sequence for this, just the read length for longest gap computation.
                // The paths in the seeds know the hit length.
                // We assume all overlapping hits are exclusive.
                unordered_map<size_t, unordered_map<size_t, vector<Path>>> paths_between_seeds = find_connecting_paths(extended_seeds,
                    aln.sequence().size());

                // Make a MultipathAlignment and feed in all the extended seeds as subpaths
                MultipathAlignment mp;
                mp.set_sequence(aln.sequence());
                mp.set_quality(aln.quality());
                for (auto& extended_seed : extended_seeds) {
                    Subpath* s = mp.add_subpath();
                    // Copy in the path.
                    *s->mutable_path() = extended_seed.first;
                    // Score it
                    s->set_score(get_regular_aligner()->score_partial_alignment(aln, gbwt_graph, extended_seed.first,
                        aln.sequence().begin() + extended_seed.second));
                    // The position in the read it occurs at will be handled by the multipath topology.
                    if (extended_seed.second == 0) {
                        // But if it occurs at the very start of the read we need to mark that now.
                        mp.add_start(mp.subpath_size() - 1);
                    }
                }

                for (auto& kv : paths_between_seeds[numeric_limits<size_t>::max()]) {
                    // For each source extended seed
                    const size_t& source = kv.first;
                    
                    // Grab the part of the read sequence that comes before it
                    string before_sequence = aln.sequence().substr(0, extended_seeds[source].second); 
                    
#ifdef debug
                    cerr << "There is a path into source extended seed " << source
                        << ": \"" << before_sequence << "\" against " << kv.second.size() << " haplotypes" << endl;
#endif

                    // We want the best alignment, to the base graph, done against any target path
                    Path best_path;
                    // And its score
                    int64_t best_score = numeric_limits<int64_t>::min();

                    // We can align it once per target path
                    for (auto& path : kv.second) {
                        // For each path we can take to get to the source
                        
                        if (path.mapping_size() == 0) {
                            // We might have extra read before where the graph starts. Handle leading insertions.
                            // We consider a pure softclip.
                            // We don't consider an empty sequence because if that were the case
                            // we would not have any paths_between_seeds entries for the dangling-left-sequence sentinel.
                            if (best_score < 0) {
                                best_score = 0;
                                best_path.clear_mapping();
                                Mapping* m = best_path.add_mapping();
                                Edit* e = m->add_edit();
                                e->set_from_length(0);
                                e->set_to_length(before_sequence.size());
                                e->set_sequence(before_sequence);
                                // Since the softclip consumes no graph, we place it on the node we are going to.
                                *m->mutable_position() = extended_seeds[source].first.mapping(0).position();
                                
#ifdef debug
                                cerr << "New best alignment: " << pb2json(best_path) << endl;
#endif
                            }
                        } else {

                            // Make a subgraph.
                            // TODO: don't copy the path
                            PathSubgraph subgraph(&gbwt_graph, path);
                            
                            // Do right-pinned alignment to the path subgraph with GSSWAligner.
                            Alignment before_alignment;
                            before_alignment.set_sequence(before_sequence);
                            // TODO: pre-make the topological order
                            
#ifdef debug
                            cerr << "Align " << pb2json(before_alignment) << " pinned right vs:" << endl;
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

                            get_regular_aligner()->align_pinned(before_alignment, subgraph, false);
                            // TODO: This should assign full length bonus! Does it?

                            if (before_alignment.score() > best_score) {
                                // This is a new best alignment. Translate from subgraph into base graph and keep it
                                best_path = subgraph.translate_down(before_alignment.path());
                                best_score = before_alignment.score();
                                
#ifdef debug
                                cerr << "New best alignment against: " << pb2json(path) << " is " << pb2json(best_path) << endl;
#endif
                            }
                        }
                    }
                    
                    // We really should have gotten something
                    assert(best_path.mapping_size() != 0);

                    // Put it in the MultipathAlignment
                    Subpath* s = mp.add_subpath();
                    *s->mutable_path() = std::move(best_path);
                    s->set_score(best_score);
                    
                    // And make the edge from it to the correct source
                    s->add_next(source);
                    
#ifdef debug
                    cerr << "Resulting source subpath: " << pb2json(*s) << endl;
#endif
                    
                    // And mark it as a start subpath
                    mp.add_start(mp.subpath_size() - 1);
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
                    size_t from_end = extended_seeds[from].second + path_from_length(extended_seeds[from].first);
                    
                    for (auto& to_and_paths : from_and_edges.second) {
                        const size_t& to = to_and_paths.first;
                        // For all the edges to other extended seeds
                        
#ifdef debug
                        cerr << "Consider " << to_and_paths.second.size() << " paths between extended seeds " << to << " and " << from << endl;
#endif

                        if (to == numeric_limits<size_t>::max()) {
                            // Do a bunch of left pinned alignments for the tails.
                            
                            // Find the sequence
                            string trailing_sequence = aln.sequence().substr(from_end);
                            
                            if (!trailing_sequence.empty()) {
                                // There is actual trailing sequence to align on this escape path

                                // Find the best path in backing graph space
                                Path best_path;
                                // And its score
                                int64_t best_score = numeric_limits<int64_t>::min();

                                // We can align it once per target path
                                for (auto& path : to_and_paths.second) {
                                    // For each path we can take to leave the "from" sink
                                    
                                    if (path.mapping_size() == 0) {
                                        // Consider the case of a nonempty trailing
                                        // softclip that bumped up against the end
                                        // of the underlying graph.
                                        
                                        if (best_score < 0) {
                                            best_score = 0;
                                            best_path.clear_mapping();
                                            Mapping* m = best_path.add_mapping();
                                            Edit* e = m->add_edit();
                                            e->set_from_length(0);
                                            e->set_to_length(trailing_sequence.size());
                                            e->set_sequence(trailing_sequence);
                                            // We need to set a position at the end of where we are coming from.
                                            const Mapping& prev_mapping = extended_seeds[from].first.mapping(
                                                extended_seeds[from].first.mapping_size() - 1);
                                            const Position& coming_from = prev_mapping.position();
                                            size_t last_node_length = gbwt_graph.get_length(gbwt_graph.get_handle(coming_from.node_id()));
                                            m->mutable_position()->set_node_id(coming_from.node_id());
                                            m->mutable_position()->set_is_reverse(coming_from.is_reverse());
                                            m->mutable_position()->set_offset(last_node_length);
                                            
                                            // We should only have this case if we are coming from the end of a node.
                                            assert(mapping_from_length(prev_mapping) + coming_from.offset() == last_node_length);
                                        }
                                    } else {

                                        // Make a subgraph.
                                        // TODO: don't copy the path
                                        PathSubgraph subgraph(&gbwt_graph, path);
                                        
                                        // Do left-pinned alignment to the path subgraph
                                        Alignment after_alignment;
                                        after_alignment.set_sequence(trailing_sequence);
                                        // TODO: pre-make the topological order

    #ifdef debug
                                        cerr << "Align " << pb2json(after_alignment) << " pinned left vs:" << endl;
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

                                        get_regular_aligner()->align_pinned(after_alignment, subgraph, true);

                                        if (after_alignment.score() > best_score) {
                                            // This is a new best alignment. Translate from subgraph into base graph and keep it
                                            best_path = subgraph.translate_down(after_alignment.path());
                                            best_score = after_alignment.score();
                                        }
                                    }
                                }
                                
                                // We need to come after from with this path

                                // We really should have gotten something
                                assert(best_path.mapping_size() != 0);

                                // Put it in the MultipathAlignment
                                Subpath* s = mp.add_subpath();
                                *s->mutable_path() = std::move(best_path);
                                s->set_score(best_score);
                                
                                // And make the edge to hook it up
                                mp.mutable_subpath(from)->add_next(mp.subpath_size() - 1);
                            }
                            
                            // If there's no sequence to align on the path going off to nowhere, don't do anything.
                            
                        } else {
                            // Do alignments between from and to

                            // Find the sequence
                            assert(extended_seeds[to].second >= from_end);
                            string intervening_sequence = aln.sequence().substr(from_end, extended_seeds[to].second - from_end); 

                            // Find the best path in backing graph space (which may be empty)
                            Path best_path;
                            // And its score
                            int64_t best_score = numeric_limits<int64_t>::min();

                            // We can align it once per target path
                            for (auto& path : to_and_paths.second) {
                                // For each path we can take to get to the source
                                
                                if (path.mapping_size() == 0) {
                                    // We're aligning against nothing
                                    if (intervening_sequence.empty()) {
                                        // Consider the nothing to nothing alignment, score 0
                                        if (best_score < 0) {
                                            best_score = 0;
                                            best_path.clear_mapping();
                                        }
                                    } else {
                                        // Consider the something to nothing alignment.
                                        // We can't use the normal code path because the BandedGlobalAligner 
                                        // wouldn't be able to generate a position form an empty graph.
                                        
                                        // We know the extended seeds we are between won't start/end with gaps, so we own the gap open.
                                        int64_t score = get_regular_aligner()->score_gap(intervening_sequence.size());
                                        if (score > best_score) {
                                            best_path.clear_mapping();
                                            Mapping* m = best_path.add_mapping();
                                            Edit* e = m->add_edit();
                                            e->set_from_length(0);
                                            e->set_to_length(intervening_sequence.size());
                                            e->set_sequence(intervening_sequence);
                                            // We can copy the position of where we are going to, since we consume no graph.
                                            *m->mutable_position() = extended_seeds[to].first.mapping(0).position();
                                        }
                                    }
                                } else {

                                    // Make a subgraph.
                                    // TODO: don't copy the path
                                    PathSubgraph subgraph(&gbwt_graph, path);
                                    
                                    // Do global alignment to the path subgraph
                                    Alignment between_alignment;
                                    between_alignment.set_sequence(intervening_sequence);
                                    
#ifdef debug
                                    cerr << "Align " << pb2json(between_alignment) << " global vs:" << endl;
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
#endif
                                    
                                    get_regular_aligner()->align_global_banded(between_alignment, subgraph, 5, true);
                                    
                                    if (between_alignment.score() > best_score) {
                                        // This is a new best alignment. Translate from subgraph into base graph and keep it
                                        best_path = subgraph.translate_down(between_alignment.path());
#ifdef debug
                                        cerr << "\tNew best: " << pb2json(best_path) << endl;
#endif
                                        
                                        best_score = between_alignment.score();
                                    }
                                }
                                
                            }
                            
                            // We may have an empty path. That's fine.

                            if (best_path.mapping_size() == 0 && intervening_sequence.empty()) {
                                // We just need an edge from from to to
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
                            }

                        }

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
                
                // Fill in the path and score
                optimal_alignment(mp, out, true);
                // Compute the identity from the path.
                out.set_identity(identity(out.path()));

                // Then continue so we don't emit the unaligned Alignment
                continue;
            }
        }
        
        // If we get here, either there was no cluster or the cluster produced no extension
        
        // Read was not able to be mapped.
        // Make our output alignment un-aligned.
        out.clear_path();
        out.set_score(0);
        out.set_identity(0);
    }
    
    // Sort again by actual score instead of cluster coverage
    std::sort(aligned.begin(), aligned.end(), [](const Alignment& a, const Alignment& b) -> bool {
        // Return true if a must come before b (i.e. it has a larger score)
        return a.score() > b.score();
    });
    
    if (aligned.size() > max_multimaps) {
        // Drop the lowest scoring alignments
        aligned.resize(max_multimaps);
    }
    
    for (size_t i = 0; i < aligned.size(); i++) {
        // For each output alignment in score order
        auto& out = aligned[i];
        
        // Assign primary and secondary status
        out.set_is_secondary(i > 0);
        out.set_mapping_quality(0);
    }
    
    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    
    if (!aligned.empty()) {
        // Annotate the primary alignment with mapping runtime
        set_annotation(aligned[0], "map_seconds", elapsed_seconds.count());
    }
    
    // Ship out all the aligned alignments
    alignment_emitter.emit_mapped_single(std::move(aligned));
}

unordered_map<size_t, unordered_map<size_t, vector<Path>>>
MinimizerMapper::find_connecting_paths(const vector<pair<Path, size_t>>& extended_seeds, size_t read_length) const {

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
        // Where does it start?
        auto& pos = extended_seeds[i].first.mapping(0).position();

        // Get the handle it is on
        handle_t handle = gbwt_graph.get_handle(pos.node_id(), pos.is_reverse());

        // Record that this extension starts at this offset along that handle
        extensions_by_handle[handle].emplace_back(pos.offset(), i);

        // Assume it is a source
        sources.insert(i);

#ifdef debug
        cerr << "Extended seed " << i << " starts on node " << pos.node_id() << " " << pos.is_reverse()
            << " at offset " << pos.offset() << endl;
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
        auto& last_mapping = extended_seeds[i].first.mapping(extended_seeds[i].first.mapping_size() - 1);
        Position cut_pos_graph = last_mapping.position();
        cut_pos_graph.set_offset(cut_pos_graph.offset() + mapping_from_length(last_mapping));
        // And the read?
        size_t cut_pos_read = extended_seeds[i].second + path_to_length(extended_seeds[i].first);

#ifdef debug
        cerr << "Extended seed " << i << ": cut after on node " << cut_pos_graph.node_id() << " " << cut_pos_graph.is_reverse()
            << " at point " << cut_pos_graph.offset() << endl;
#endif

        // Get a handle in the GBWTGraph
        handle_t start_handle = gbwt_graph.get_handle(cut_pos_graph.node_id(), cut_pos_graph.is_reverse());

        // Decide if we need to actually do GBWT search, or if we can find a destination on the same node we ended on
        bool do_gbwt_search = true;

        // Look on the same graph node.
        // See if we hit any other extensions on this node.
        auto same_node_found = extensions_by_handle.find(start_handle);
        if (same_node_found != extensions_by_handle.end()) {
            // If we have extended seeds
            
            for (auto& next_offset_and_index : same_node_found->second) {
                // Scan them in order.
                // TODO: Skip to after ourselves.
                
                if (extended_seeds[next_offset_and_index.second].second >= cut_pos_read && next_offset_and_index.first >= cut_pos_graph.offset()) { 
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

#ifdef debug
                    cerr << "Found path on node between seed at "
                    << pb2json(extended_seeds[i].first.mapping(extended_seeds[i].first.mapping_size() - 1).position())
                    << " and seed at "
                    << pb2json(extended_seeds[next_offset_and_index.second].first.mapping(0).position())
                    << ":" << endl << "\t" << pb2json(connecting) <<  endl;
#endif

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

        // Search everything in the GBWT graph right from the end of the start extended seed, up to the limit.
        explore_gbwt(cut_pos_graph, search_limit, [&](const Path& here_path, const handle_t& there_handle) -> bool {
            // When we encounter a new handle visited by haplotypes extending off of the last node in a Path

            // See if we hit any other extensions on this next node
            auto found = extensions_by_handle.find(there_handle);
            if (found != extensions_by_handle.end()) {
                // If we do
                
                for (auto& next_offset_and_index : found->second) {
                    // Look at them in order along the node
                    
                    if (extended_seeds[next_offset_and_index.second].second >= cut_pos_read) { 
                        // As soon as we find one that starts in the read after our start extended seed ended

                        // Extend the Path to connect to it.
                        // TODO: Make these shared tail lists for better algorithmics
                        Path extended = here_path;
                        
                        if (next_offset_and_index.first > 0) {
                            // There is actual material on this new node before the extended seed we have to hit.
                            Mapping* m = extended.add_mapping();
                            m->mutable_position()->set_node_id(gbwt_graph.get_id(there_handle));
                            m->mutable_position()->set_is_reverse(gbwt_graph.get_is_reverse(there_handle));
                            Edit* e = m->add_edit();
                            // Make sure it runs through the last base *before* the extended seed we are going for
                            e->set_from_length(next_offset_and_index.first);
                            e->set_to_length(next_offset_and_index.first);
                        }

#ifdef debug
                        cerr << "Found graph path between seed at "
                            << pb2json(extended_seeds[i].first.mapping(extended_seeds[i].first.mapping_size() - 1).position())
                            << " and seed at "
                            << pb2json(extended_seeds[next_offset_and_index.second].first.mapping(0).position())
                            << ":" << endl << "\t" << pb2json(extended) <<  endl;
#endif

                        // And emit that connection
                        to_return[i][next_offset_and_index.second].emplace_back(std::move(extended));

                        // Record that the destination is not a source
                        sources.erase(next_offset_and_index.second);

                        // Don't look at any more destinations, or any extensions past this node of this search state.
                        return false;
                    }
                }
            }

            // Otherwise we didn't hit anything we can stop at. Keep extending.
            return true;
        }, [&](const Path& limit_path) {
            // When we blow past the walk distance limit or hit a dead end

            // We have a way to escape.
            // Save that as a path.
            to_return[i][numeric_limits<size_t>::max()].emplace_back(limit_path);
            
            // If we end up with paths anywhere else we will destroy it, so we
            // will only keep it for sinks. Note that we might not have any
            // sequence to escape *with*! The caller has to look out for that
            // case and not waste time trying to align nothing, since we don't
            // have access to the read length here.
        });
    }

    for (auto& kv : to_return) {
        // For each seed, if it can reach anything *other* than numeric_limits<size_t>::max(), erase anything going to numeric_limits<size_t>::max()
        auto found = kv.second.find(numeric_limits<size_t>::max());
        if (kv.second.size() > 1 && found != kv.second.end()) {
            // We have a going-off-to-nothing path and also a path to somewhere else.
            // Don't go off to nothing.
            kv.second.erase(found);
        }
    }

    // Now we need the paths *from* numeric_limits<size_t>::max() to sources.
    // Luckily we know the sources.
    for (const size_t& i : sources) {
        // For each source
        
#ifdef debug
        cerr << "Extended seed " << i << " is a source" << endl;
#endif
        
        if (extended_seeds[i].second > 0) {
#ifdef debug
            cerr << "\tIt is not at the start of the read, so there is a left tail" << endl;
#endif

            // Find its start
            Position start = extended_seeds[i].first.mapping(0).position();
            
#ifdef debug
            cerr << "\tPosition read-forward to search left before: " << pb2json(start) << endl;
#endif

            // Flip it around to face left
            start = reverse(start, gbwt_graph.get_length(gbwt_graph.get_handle(start.node_id())));
            
#ifdef debug
            cerr << "\tPosition read-reverse to search right after: " << pb2json(start) << endl;
#endif

            // Now the search limit is all the read *before* the seed, plus the detectable gap
            size_t search_limit = get_regular_aligner()->longest_detectable_gap(read_length, extended_seeds[i].second) +
                extended_seeds[i].second;

            // Start another search, but going left.
            explore_gbwt(start, search_limit, [&](const Path& here_path, const handle_t& there_handle) -> bool {
                // If we weren't reachable from anyone, nobody should be reachable from us going the other way.
                // So always keep going.
                return true;
            }, [&](const Path& limit_path) {
                // We have this path going right from start and hitting the walk limit or the edge of the graph.

                for (auto& mapping : limit_path.mapping()) {
                    // Make sure nothing has a negative offset before flipping.
                    assert(mapping.position().offset() >= 0);
                }

                // Flip the path around
                Path flipped = reverse_complement_path(limit_path, [&](id_t id) -> size_t {
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

void MinimizerMapper::explore_gbwt(const Position& from, size_t walk_distance,
    const function<bool(const Path&, const handle_t&)>& visit_callback,
    const function<void(const Path&)>& limit_callback) const {
    
#ifdef debug
    cerr << "Exploring GBWT out from " << pb2json(from) << " to distance " << walk_distance << endl;
#endif
    
    // Holds the gbwt::SearchState we are at, and the Path from the end of the starting
    // seed up through the end of the node we just searched.
    // The from_length of the path tracks our consumption of distance limit.
    using traversal_state_t = pair<gbwt::SearchState, Path>;
    
    // Get a handle to the node the from position is on, in its forward orientation
    handle_t start_handle = gbwt_graph.get_handle(from.node_id(), from.is_reverse());

    // Turn it into a SearchState
    gbwt::SearchState start_state = gbwt_graph.get_state(start_handle);

    // The search state represents searching through the end of the node, so we have to consume that much search limit.

    // Tack on how much search limit distance we consume by going to the end of
    // the node. Our start position is a cut *between* bases, and we take everything after it.
    // If the cut is at the offset of the whole length of the node, we take 0 bases.
    // If it is at 0, wer take all the bases in the node.
    size_t distance_to_node_end = gbwt_graph.get_length(start_handle) - from.offset();    
    
    // And make a Path that represents the part of the node we're on that goes out to the end.
    // This may be empty if the hit already stopped at the end of the node
    Path path_to_end;
    if (distance_to_node_end != 0) {
        // We didn't hit the end of the node already.

        // Make a mapping that starts on the right side of the cut we started our search at.
        Mapping* m = path_to_end.add_mapping();
        *m->mutable_position() = from;
        m->mutable_position()->set_offset(m->position().offset());

        // Make it the requested length of perfect match.
        Edit* e = m->add_edit();
        e->set_from_length(distance_to_node_end);
        e->set_to_length(distance_to_node_end);
    }
   
#ifdef debug   
    cerr << "Starting traversal with " << pb2json(path_to_end) << " from " << pb2json(from) << endl;
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
    // And track end node IDs
    unordered_set<id_t> end_ids;

    while (!queue.empty()) {
        // While there are things in the queue
        
        // Record max size
        queue_max = max(queue_max, queue_size);

        // Grab one
        traversal_state_t here(std::move(queue.front()));
        queue.pop_front();
        queue_size--;
        gbwt::SearchState& here_state = here.first;
        Path& here_path = here.second;
        
        // follow_paths on it
        bool got_anywhere = false;
        gbwt_graph.follow_paths(here_state, [&](const gbwt::SearchState& there_state) -> bool {
            if (there_state.empty()) {
                // Ignore places that no haplotypes go, and get the next place instead.
                return true;
            }
            
            // For each place it can go
            handle_t there_handle = gbwt_graph.node_to_handle(there_state.node);
            
            // Record that we got there
            got_anywhere = true;

            // Say we can go from here to there. Should we?
            bool continue_extending = visit_callback(here_path, there_handle);

            if (continue_extending) {
                // Generate the path we take if we take up all of that node.
                Path extended = here_path;
                Mapping* m = extended.add_mapping();
                m->mutable_position()->set_node_id(gbwt_graph.get_id(there_handle));
                m->mutable_position()->set_is_reverse(gbwt_graph.get_is_reverse(there_handle));
                Edit* e = m->add_edit();
                e->set_from_length(gbwt_graph.get_length(there_handle));
                e->set_to_length(gbwt_graph.get_length(there_handle));


                // See if we can get to the end of the node without going outside the search length.
                if (path_from_length(here_path) + gbwt_graph.get_length(there_handle) <= walk_distance) {
                    // If so, continue the search
                    queue.emplace_back(there_state, extended);
                    queue_size++;
                } else {
                    // Report that, with this extension, we hit the limit.
                    limit_callback(extended);
                    limit_hits++;
                    end_ids.insert(gbwt_graph.get_id(there_handle));
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
    
    cerr << "Queue max: " << queue_max << " Limit hits: " << limit_hits << " Dead ends: " << dead_ends << " End IDs: " << end_ids.size() << endl;
}


}


