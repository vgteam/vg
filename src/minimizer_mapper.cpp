/**
 * \file minimizer_mapper.cpp
 * Defines the code for the minimizer-and-GBWT-based mapper.
 */

#include "minimizer_mapper.hpp"

#include "crash.hpp"
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
// Make sure by-index references are correct
//#define debug_validate_index_references

namespace vg {

using namespace std;

MinimizerMapper::MinimizerMapper(const gbwtgraph::GBWTGraph& graph,
    const gbwtgraph::DefaultMinimizerIndex& minimizer_index,
    SnarlDistanceIndex* distance_index, 
    const vector<ZipCode>* zipcodes, 
    const PathPositionHandleGraph* path_graph) :
    path_graph(path_graph), minimizer_index(minimizer_index),
    distance_index(distance_index),  
    zipcodes(zipcodes),
    clusterer(distance_index, &graph),
    zip_clusterer(distance_index, &graph),
    gbwt_graph(graph),
    extender(gbwt_graph, *(get_regular_aligner())),
    fragment_length_distr(1000,1000,0.95) {
    
    // The GBWTGraph needs a GBWT
    crash_unless(graph.index != nullptr);
}

//-----------------------------------------------------------------------------

string MinimizerMapper::log_name() {
    return "T" + to_string(omp_get_thread_num()) + ":\t";
}

string MinimizerMapper::log_alignment(const Alignment& aln) {
    if (aln.sequence().size() < LONG_LIMIT && aln.path().mapping_size() < LONG_LIMIT/32 ) {
        // Log as a short alignment
        return pb2json(aln);
    } else {
        // Log as a long alignment
        
        stringstream ss;
        ss << log_alignment(aln.path(), true);
        ss << " score " << aln.score();
        return ss.str();
    }
}

string MinimizerMapper::log_alignment(const Path& path, bool force_condensed) {
    if (path.mapping_size() < LONG_LIMIT/32 && !force_condensed) {
        // Log as a short alignment
        return pb2json(path);
    } else {
        // Log as a long alignment
        
        // Turn it into one big CIGAR string
        vector<pair<int, char>> cigar;
        for (auto& mapping : path.mapping()) {
            mapping_cigar(mapping, cigar);
        }
        
        // And then put that
        stringstream ss;
        if (!cigar.empty()) {
            ss << cigar_string(cigar);
        } else {
            ss << "<empty alignment>";
        }
        if (path.mapping_size() > 0) {
            ss << "@" << path.mapping(0).position().node_id()
                << (path.mapping(0).position().is_reverse() ? '-' : '+');
        }
        return ss.str();
        
        // TODO: Include visited nodes as collapsed ranges of consecutive IDs?
    }
}

string MinimizerMapper::log_bits(const std::vector<bool>& bits) {
    stringstream ss;
    
    // We're going to take groups of 8 bits into here, and then emit them as
    // 8-dot Braille characters.
    unsigned char scratch = 0;
    for (size_t i = 0; i < bits.size(); i++) {
        // We want earlier bits to be top and left aligned, and this is easier
        // if we make later bits the high bits and push earlier bits to the low
        // bits.
        //
        // So shift everything down and put in at the highest place.
        scratch = (scratch >> 1) | ((unsigned char)bits[i] << 7);
        if (i % 8 == 7 || i + 1 == bits.size()) {
            // Time to ship it
            
            // Make sure to finish aligning the first bit down if this is a
            // partial byte.
            scratch = scratch >> (7 - i % 8);
            
            // We want our layout to be (numbering LSB to MSB):
            //
            // 0 4
            // 1 5
            // 2 6
            // 3 7
            //
            // Unicode Braille is laid out:
            //
            // 0 3
            // 1 4
            // 2 5
            // 6 7
            //
            // So we have to:
            // - Keep 012 and 7 in place
            // - Move 3 up by 3 to become 6
            // - Move 456 down by 1 to become 345
            scratch = (scratch & 0b10000111) |
                     ((scratch & 0b00001000) << 3) |
                     ((scratch & 0b01110000) >> 1);
            
            // Braille block codepoints would be 0x2800 + scratch
            
            // Encode as UTF-8. These are all above U+0800 and below U+FFFF so
            // we need 3 bytes.
            // First 4 bits are always 0xe.
            // First 4 data bits (0x2) go in low nibble of first byte.
            ss << (char) 0xe2;
            // Next two bits are always 0b10.
            // Next 4 data bits (0x8) and high 2 bits of scratch go in next byte.
            ss << (char) ((0b10 << 6) | (0x80 >> 2) | (scratch >> 6));
            // Next two bits are also always 0b10.
            // Low 6 bits of scratch go in last byte.
            ss << (char) ((0b10 << 6) | (scratch & 0b00111111));
            
            // Clear out since the last cycle may be partial and we don't want
            // stray bits.
            scratch = 0;
        }
    }
    return ss.str();
}

void MinimizerMapper::dump_chaining_problem(const std::vector<algorithms::Anchor>& anchors, const std::vector<size_t>& cluster_seeds_sorted, const HandleGraph& graph) {
    ProblemDumpExplainer exp;
    
    // We need to keep track of all the points we want in our problem subgraph.
    std::vector<pos_t> seed_positions;
    seed_positions.reserve(cluster_seeds_sorted.size());
    
    exp.object_start();
    
    // Save all the items
    exp.key("items");
    exp.array_start();
    for (auto& index : cluster_seeds_sorted) {
        exp.object_start();
        exp.key("read_start");
        exp.value(anchors[index].read_start());
        exp.key("read_end");
        exp.value(anchors[index].read_end());
        
        pos_t graph_start = anchors[index].graph_start();
        seed_positions.push_back(graph_start);
        exp.key("graph_start");
        exp.value(graph_start);
        exp.key("graph_end");
        exp.value(anchors[index].graph_end());
        
        exp.object_end();
    }
    exp.array_end();
    
    // Get the subgraph for the cluster
    HashGraph subgraph;
    algorithms::extract_containing_graph(&graph, &subgraph, seed_positions, 10000);
    exp.key("subgraph");
    exp.value(subgraph);
    
    exp.object_end();
}

void MinimizerMapper::dump_debug_sequence(ostream& out, const string& sequence, size_t start_offset, size_t length_limit) {
    int digits_needed = (int) ceil(log10(std::min(sequence.size(), start_offset + length_limit)));
    for (int digit = digits_needed - 1; digit >= 0; digit--) {
        out << log_name();
        for (size_t i = start_offset; i < std::min(sequence.size(), start_offset + length_limit); i++) {
            // Output the correct digit for this place in this number
            out << (char) ('0' + (uint8_t) floor(i % (int) round(pow(10, digit + 1)) / pow(10, digit)));
        }
        out << endl;
    }
    out << log_name() << sequence.substr(start_offset, std::min(sequence.size() - start_offset, length_limit)) << endl;
}

void MinimizerMapper::dump_debug_extension_set(const HandleGraph& graph, const Alignment& aln, const vector<GaplessExtension>& extended_seeds) {
    
    if (aln.sequence().size() >= LONG_LIMIT) {
        // Describe the extensions, because the read is huge
        cerr << log_name() << "<" << extended_seeds.size() << " extensions>" << endl;
    } else {
        // Show a diagram
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
}

void MinimizerMapper::dump_debug_minimizers(const VectorView<MinimizerMapper::Minimizer>& minimizers, const string& sequence, const vector<size_t>* to_include, size_t start_offset, size_t length_limit) {

    // Work out what region we are actually looking at
    size_t region_start = start_offset;
    size_t region_length = std::min(sequence.size() - start_offset, length_limit);

    if (region_length >= LONG_LIMIT) {
        // Describe the minimizers, because the read is huge
        size_t minimizer_count = to_include ? to_include->size() : minimizers.size();
        if (minimizer_count < MANY_LIMIT) {
            auto print_minimizer = [&](size_t i) {
                cerr << log_name() << "Minimizer " << i << ": " << minimizers[i].forward_sequence() << "@" << minimizers[i].forward_offset() << " with " << minimizers[i].hits << " hits" << endl;
            };
            
            if (to_include) {
                for (auto& i : *to_include) {
                    print_minimizer(i);
                }
            } else {
                for (size_t i = 0; i < minimizers.size(); i++) {
                     print_minimizer(i);
                }
            }
        } else {
            if (region_start == 0 && length_limit == sequence.size()) {
                // Report as if we have a count
                cerr << log_name() << "<" << minimizer_count << " minimizers>" << endl;
            } else {
                // We don't know how many minimizers are actually in the region
                cerr << log_name() << "<Minimizers from " << region_start << " to " << (region_start + region_length) << ">" << endl;
            }
        }
    } else {
        // Draw a diagram
        dump_debug_sequence(cerr, sequence, region_start, region_length);
        
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
            auto& m = minimizers[index];
            
            if (m.forward_offset() < region_start || m.forward_offset() - region_start + m.length > region_length) {
                // Minimizer itself reaches out of bounds, so hide it
                continue;
            }
            
            // Compute its traits relative to the region we are interested in showing
            size_t relative_agglomeration_start = (m.agglomeration_start < region_start) ? (size_t)0 : m.agglomeration_start - region_start;
            size_t relative_forward_offset = m.forward_offset() - region_start;
            size_t capped_agglomeration_length = (relative_agglomeration_start + m.agglomeration_length > region_length) ? (region_length - relative_agglomeration_start) : m.agglomeration_length;
            
            cerr << log_name();
            
            
            for (size_t i = 0; i < relative_agglomeration_start; i++) {
                // Space until its agglomeration starts
                cerr << ' ';
            }
            
            for (size_t i = relative_agglomeration_start; i < relative_forward_offset; i++) {
                // Do the beginnign of the agglomeration
                cerr << '-';
            }
            // Do the minimizer itself
            cerr << m.value.key.decode(m.length);
            for (size_t i = relative_forward_offset + m.length ; i < relative_agglomeration_start + capped_agglomeration_length; i++) {
                // Do the tail end of the agglomeration
                cerr << '-';
            }
            
            // Tag with metadata
            cerr << " (#" << index << ", " << m.hits << " hits)" << endl;
        }
    }
}

void MinimizerMapper::dump_debug_clustering(const Cluster& cluster, size_t cluster_number, const VectorView<Minimizer>& minimizers, const std::vector<Seed>& seeds) {
    if (minimizers.size() < MANY_LIMIT) {
        // There are a few minimizers overall, so describe each in the cluster individually.
        for (auto hit_index : cluster.seeds) {
            cerr << log_name() << "Minimizer " << seeds[hit_index].source << " is present in cluster " << cluster_number << endl;
        }
    } else {
        // Describe the seeds in aggregate
        std::vector<bool> presence(minimizers.size());
        for (auto hit_index : cluster.seeds) {
            presence[seeds[hit_index].source] = true;
        }
        cerr << log_name() << "Minimizers in cluster " << cluster_number << "\t" << log_bits(presence) << endl;
    }
}
bool MinimizerMapper::validate_clusters(const std::vector<std::vector<Cluster>>& clusters, const std::vector<std::vector<Seed>>& seeds, size_t read_limit, size_t fragment_limit) const {
    vector<vector<pos_t>> fragment_clusters;
    for (size_t read_num = 0 ; read_num < clusters.size() ; read_num ++) {
        auto& one_read_clusters = clusters[read_num];
        if (one_read_clusters.size() > 0) {
            for (size_t cluster_num1 = 0; cluster_num1 < one_read_clusters.size(); cluster_num1++) {
                // For each cluster -cluster this cluster to ensure that there is only one
                vector<size_t> clust = one_read_clusters[cluster_num1].seeds;
                size_t fragment_cluster = one_read_clusters[cluster_num1].fragment;
                if (fragment_cluster >= fragment_clusters.size()) {
                    fragment_clusters.resize(fragment_cluster+1);
                }
                
                structures::UnionFind new_clusters (clust.size(), false);
                for (size_t i1 = 0 ; i1 < clust.size() ; i1++) {
                    pos_t pos1 = seeds[read_num][clust[i1]].pos;
                    fragment_clusters[fragment_cluster].emplace_back(pos1);
                    
                    for (size_t cluster_num2 = 0 ; cluster_num2 < one_read_clusters.size() ; cluster_num2++) {
                        if (cluster_num2 != cluster_num1) {
                            //For each other cluster, make sure that the seeds in the other cluster are far away
                            vector<size_t> clust2 = one_read_clusters[cluster_num2].seeds;
                            for (size_t i2 = 0 ; i2 < clust2.size() ; i2++) {
                                //And each position in each other cluster,
                                //make sure that this position is far away from i1
                                pos_t pos2 = seeds[read_num][clust2[i2]].pos;;
                                size_t distance_between = unoriented_distance_between(pos1, pos2);
                                
                                if ( distance_between != std::numeric_limits<int64_t>::max() && 
                                     distance_between <= read_limit) {
                                    cerr << "These should have been in the same read cluster: " ;
                                    cerr << pos1 << " and " << pos2 << " with distance between them " << distance_between << endl;
                                    return false;
                                }
                    
                            }
                        }
                    }
                    for (size_t i2 = 0 ; i2 < clust.size() ; i2++) {
                        //For each position in the same cluster, make sure it's close to something
                        pos_t pos2 = seeds[read_num][clust[i2]].pos;
                        size_t distance_between = unoriented_distance_between(pos1, pos2);
                        if ( distance_between != std::numeric_limits<int64_t>::max() && 
                            distance_between <= read_limit) {
                            new_clusters.union_groups(i1, i2);
                        }
                
                    }
                }
                auto actual_clusters = new_clusters.all_groups();
                if (actual_clusters.size() != 1) {
                    cerr << "These should be different read clusters: " << endl;
                    for (auto c : actual_clusters) {
                        cerr << "cluster: " ;
                        for (size_t i : c) {
                            cerr << seeds[read_num][clust[i]].pos << " ";
                        }
                        cerr << endl;
                    }
                    return false;
                }
            }
        }
    }
    //Now check the fragment clusters
    if (fragment_limit != 0) {
        for (size_t cluster_num1 = 0; cluster_num1 < fragment_clusters.size(); cluster_num1++) {
            // For each cluster -cluster this cluster to ensure that 
            // there is only one
            vector<pos_t> clust = fragment_clusters[cluster_num1];
        
            structures::UnionFind new_clusters (clust.size(), false);
        
            for (size_t i1 = 0 ; i1 < clust.size() ; i1++) {
                pos_t pos1 = clust[i1];
        
                for (size_t cluster_num2 = 0 ; cluster_num2 < fragment_clusters.size() ; cluster_num2++) {
                    if (cluster_num2 != cluster_num1) {
                        //For each other cluster
                        vector<pos_t> clust2 = fragment_clusters[cluster_num2];
                        for (size_t i2 = 0 ; i2 < clust2.size() ; i2++) {
                            //And each position in each other cluster,
                            //make sure that this position is far away from i1
                            pos_t pos2 = clust2[i2];
                            int64_t distance_between = unoriented_distance_between(pos1, pos2);
                            if ( distance_between != std::numeric_limits<int64_t>::max() && 
                                distance_between <= fragment_limit) {
                                cerr << "These should have been in the same fragment cluster: " ;
                                cerr << pos1 << " and " << pos2 << " with distance between " << distance_between << " and distance limit " << fragment_limit << endl;
                                return false;
                                
                            }
                
                        }
                    }
                }
                for (size_t i2 = 0 ; i2 < clust.size() ; i2++) {
                    //For each position in the same cluster
                    pos_t pos2 = clust[i2];
                    int64_t distance_between = unoriented_distance_between(pos1, pos2);
                    if ( distance_between != std::numeric_limits<int64_t>::max() && 
                        distance_between <= fragment_limit) {
                        new_clusters.union_groups(i1, i2);
                    }
            
                }
            }
            auto actual_clusters = new_clusters.all_groups();
            if (actual_clusters.size() != 1) {
                cerr << "These should be different fragment clusters with distance limit " << fragment_limit << endl;
                for (auto c : actual_clusters) {
                    cerr << "cluster: " ;
                    for (size_t i1 : c) {
                        cerr << clust[i1] << " ";
                    }
                    cerr << endl;
                }
                return false;
            }
        }
    }
    return true;
}



void MinimizerMapper::dump_debug_seeds(const VectorView<Minimizer>& minimizers, const std::vector<Seed>& seeds, const std::vector<size_t>& selected_seeds) {
    if (selected_seeds.size() < MANY_LIMIT) {
        // There are a few seeds so describe them individually.
        for (auto seed_index : selected_seeds) {
            const Seed& seed = seeds[seed_index];
            const Minimizer& minimizer = minimizers[seed.source];
            cerr << log_name() << "Seed read:" << minimizer.value.offset << (minimizer.value.is_reverse ? '-' : '+') << " = " << seed.pos
                << " from minimizer " << seed.source << "(" << minimizer.hits << "), #" << seed_index << endl;
        }
    } else {
        // Describe the seeds in aggregate
        size_t min_read_pos = std::numeric_limits<size_t>::max();
        size_t max_read_pos = 0;
        size_t min_hits = std::numeric_limits<size_t>::max();
        size_t max_hits = 0;
        for (auto seed_index : selected_seeds) {
            // Compute statistics over all the seeds
            const Seed& seed = seeds[seed_index];
            const Minimizer& minimizer = minimizers[seed.source];
            min_read_pos = std::min(min_read_pos, (size_t)minimizer.value.offset);
            max_read_pos = std::max(max_read_pos, (size_t)minimizer.value.offset);
            min_hits = std::min(min_hits, minimizer.hits);
            max_hits = std::max(max_hits, minimizer.hits);
        }
        cerr << log_name() << selected_seeds.size() << " seeds in read:" << min_read_pos << "-" << max_read_pos << " with " << min_hits << "-" << max_hits << " hits" << endl;
    }
}

void MinimizerMapper::dump_debug_query(const Alignment& aln) {
    cerr << log_name() << "Read " << aln.name() << ": " ;
    if (aln.sequence().size() < LONG_LIMIT) {
        cerr << aln.sequence();
    } else {
        cerr << "<" << aln.sequence().size() << " bp>";
    }
    cerr << endl;
}

void MinimizerMapper::dump_debug_query(const Alignment& aln1, const Alignment& aln2) {
    cerr << log_name() << "Read pair " << aln1.name() << ": ";
    if (aln1.sequence().size() < LONG_LIMIT) {
        cerr << aln1.sequence();
    } else {
        cerr << "<" << aln1.sequence().size() << " bp>";
    }
    cerr << " and " << aln2.name() << ": ";
    if (aln2.sequence().size() < LONG_LIMIT) {
        cerr << aln2.sequence();
    } else {
        cerr << "<" << aln2.sequence().size() << " bp>";
    }
    cerr << endl;
}


//-----------------------------------------------------------------------------

void MinimizerMapper::map(Alignment& aln, AlignmentEmitter& alignment_emitter) {
    // Ship out all the aligned alignments
    alignment_emitter.emit_mapped_single(map(aln));
}

vector<Alignment> MinimizerMapper::map(Alignment& aln) {
    if (align_from_chains) {
        return map_from_chains(aln);
    } else {
        return map_from_extensions(aln);
    }
}

vector<Alignment> MinimizerMapper::map_from_extensions(Alignment& aln) {
    
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

    // Find the seeds and mark the minimizers that were located.
    vector<Seed> seeds = this->find_seeds(minimizers_in_read, minimizers, aln, funnel);

    // Cluster the seeds. Get sets of input seed indexes that go together.
    if (track_provenance) {
        funnel.stage("cluster");
    }

    // Find the clusters
    std::vector<Cluster> clusters = clusterer.cluster_seeds(seeds, get_distance_limit(aln.sequence().size()));
    
#ifdef debug_validate_clusters
    vector<vector<Cluster>> all_clusters;
    all_clusters.emplace_back(clusters);
    vector<vector<Seed>> all_seeds;
    all_seeds.emplace_back(seeds);
    validate_clusters(all_clusters, all_seeds, get_distance_limit(aln.sequence().size()), 0);
#endif

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
        this->score_cluster(cluster, i, minimizers, seeds, aln.sequence().length());
        if (this->track_provenance) {
            // Record the cluster in the funnel as a group of the size of the number of items.
            funnel.merge_group(cluster.seeds.begin(), cluster.seeds.end());
            funnel.score(funnel.latest(), cluster.score);

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
    process_until_threshold_c<double>(clusters.size(), [&](size_t i) -> double {
            return clusters[i].coverage;
        }, [&](size_t a, size_t b) -> bool {
            return ((clusters[a].coverage > clusters[b].coverage) ||
                    (clusters[a].coverage == clusters[b].coverage && clusters[a].score > clusters[b].score));
        }, cluster_coverage_threshold, min_extensions, max_extensions, rng, [&](size_t cluster_num) -> bool {
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
            }
            

            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Cluster " << cluster_num << endl;
                    cerr << log_name() << "Covers " << cluster.coverage << "/best-" << cluster_coverage_threshold << " of read" << endl;
                    cerr << log_name() << "Scores " << cluster.score << "/" << cluster_score_cutoff << endl;
                }
            }
            
            // Extend seed hits in the cluster into one or more gapless extensions
            cluster_extensions.emplace_back(this->extend_cluster(
                cluster,
                cluster_num,
                minimizers,
                seeds,
                aln.sequence(),
                minimizer_extended_cluster_count,
                funnel));
            
            kept_cluster_count ++;
            
            return true;
            
        }, [&](size_t cluster_num) -> void {
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
    auto discard_processed_cluster_by_score = [&](size_t extension_num) -> void {
        // This extension is not good enough.
        if (track_provenance) {
            funnel.fail("extension-set", extension_num, cluster_extension_scores[extension_num]);
        }
        
        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "processed cluster " << extension_num << " failed because its score was not good enough (score=" << cluster_extension_scores[extension_num] << ")" << endl;
                if (track_correctness && funnel.was_correct(extension_num)) {
                    cerr << log_name() << "\tCORRECT!" << endl;
                    dump_debug_extension_set(gbwt_graph, aln, cluster_extensions[extension_num]);
                }
            }
        }
    };
    
    // Go through the gapless extension groups in score order.
    process_until_threshold_b<int>(cluster_extension_scores,
        extension_set_score_threshold, min_extension_sets, max_alignments, rng, [&](size_t extension_num) -> bool {
            // This extension set is good enough.
            // Called in descending score order.
            
            if (cluster_extension_scores[extension_num] < extension_set_min_score) {
                // Actually discard by score
                discard_processed_cluster_by_score(extension_num);
                return false;
            }
            
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
            
                // Do the DP and compute up to 2 alignments from the individual gapless extensions
                best_alignments.emplace_back(aln);
                find_optimal_tail_alignments(aln, extensions, rng, best_alignments[0], best_alignments[1]);
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

            for (size_t i = 0 ; i < minimizer_extended_cluster_count[extension_num].size() ; i++) {
                minimizer_extensions_count[i] += minimizer_extended_cluster_count[extension_num][i];
                if (minimizer_extended_cluster_count[extension_num][i] > 0) {
                    // This minimizer is in an extended cluster that gave rise
                    // to at least one alignment, so it is explored.
                    minimizer_explored.insert(i);
                }
            }
            
            return true;
        }, [&](size_t extension_num) -> void {
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
        crash_unless(false);
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

    crash_unless(!mappings.empty());
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
            annotate_with_minimizer_statistics(mappings[0], minimizers, seeds, seeds.size(), 0, funnel);
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
             crash_unless(minimizer.hits<=hard_hit_cap) ;
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

    if (track_provenance && show_work && aln.sequence().size() < LONG_LIMIT) {
        // Dump the funnel info graph.
        // TODO: Add a new flag for this.
        #pragma omp critical (cerr)
        {
            funnel.to_dot(cerr);
        }
    }

    return mappings;
}

// For MinimizerMapper::map_from_chains() see ./minimizer_mapper_from_chains.cpp

//-----------------------------------------------------------------------------

void MinimizerMapper::pair_all(std::array<vector<Alignment>, 2>& mappings) const {
    for (auto r : {0, 1}) {
        if (!mappings[r].empty()) {
            for (auto& other : mappings[1 - r]) {
                // Each read 1 needs to point to read 2 as its successor
                // Each read 2 needs to point to read 1 as its predecessor
                // If we started with the first read, other is the second read, so it needs to point to the first read as prev.
                Alignment* ref = r == 0 ? other.mutable_fragment_prev() : other.mutable_fragment_next();
                ref->set_name(mappings[r].front().name());
            }
        }
    }
}

pair<vector<Alignment>, vector<Alignment>> MinimizerMapper::map_paired(Alignment& aln1, Alignment& aln2,
                                                      vector<pair<Alignment, Alignment>>& ambiguous_pair_buffer){
    if (fragment_length_distr.is_finalized()) {

        //If we know the fragment length distribution then we just map paired ended 
        return map_paired(aln1, aln2);
    } else {
        std::array<Alignment*, 2> alns {&aln1, &aln2};
        
        std::array<vector<Alignment>, 2> single;
        std::array<int32_t, 2> max_score_aln;
        bool both_perfect_unique = true;
        for (auto r : {0, 1}) {
            //If we don't know the fragment length distribution, map the reads single ended
            single[r] = std::move(map(*alns[r]));
            // Check if the separately-mapped ends are both sufficiently perfect and sufficiently unique
            max_score_aln[r] = get_regular_aligner()->score_exact_match(*alns[r], 0, alns[r]->sequence().size());
            both_perfect_unique = both_perfect_unique && !single[r].empty() && single[r].front().mapping_quality() == 60 && single[r].front().score() >= max_score_aln[r] * 0.85;
        }
        
        if (both_perfect_unique) {
            //Flip the second alignment to get the proper fragment distance 
            reverse_complement_alignment_in_place(&single[1].front(), [&](vg::id_t node_id) {
                    return gbwt_graph.get_length(gbwt_graph.get_handle(node_id));
                    });           
            int64_t dist = distance_between(single[0].front(), single[1].front());
            // And that they have an actual pair distance and set of relative orientations

            if (dist == std::numeric_limits<int64_t>::max() ||
                dist >= max_fragment_length) {
                //If the distance between them is ambiguous or it it large enough that we don't think it's valid, leave them unmapped

                ambiguous_pair_buffer.emplace_back(aln1, aln2);
                pair<vector<Alignment>, vector<Alignment>> empty;
                return empty;
            }

            //If we're keeping this alignment, flip the second alignment back
            reverse_complement_alignment_in_place(&single[1].front(), [&](vg::id_t node_id) {
                    return gbwt_graph.get_length(gbwt_graph.get_handle(node_id));
                    });           
            // If that all checks out, say they're mapped, emit them, and register their distance and orientations
            fragment_length_distr.register_fragment_length(dist);

            std::array<vector<Alignment>, 2> mapped_pair;
            for (auto r : {0, 1}) {
                mapped_pair[r].emplace_back(std::move(single[r].front()));
            }
            pair_all(mapped_pair);

            //TODO: This is "properly paired" if there is a path between the reads. Since we
            //haven't finalized the distribution we can't compare it
            bool properly_paired = dist != std::numeric_limits<int64_t>::max();
            for (auto r : {0, 1}) {
                set_annotation(mapped_pair[r].back(), "proper_pair", properly_paired);
            }

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


            return {std::move(mapped_pair[0]), std::move(mapped_pair[1])};

        } else {
            // Otherwise, discard the mappings and put them in the ambiguous buffer

            ambiguous_pair_buffer.emplace_back(aln1, aln2);
            return {};
        }
    }
}

// For paired-end alignment we use a bunch of structures indexed by clustering
// fragment, then read, then alignment of the read. So we define some index
// types and lookup functions to deal with these.
// TODO: Make these local classes when C++ learns to let you use template members in local classes.

/// Type to point to an alignment of a known read
struct read_alignment_index_t {
    size_t fragment;
    size_t alignment;
    
    /// Allow looking up this index in a structure organized by fragment,
    /// read, and alignment, such as in alignments, or alignment_indices.
    template<typename NestedArray>
    auto lookup_for_read_in(bool read, NestedArray& a) const -> typename std::add_lvalue_reference<decltype(a[0][0][0])>::type {
        return a[fragment][read][alignment]; 
    }
    
    /// Make sure that this index actually points to an alignment of the given
    /// read in the given structure. Throws if not.
    template<typename NestedArray>
    void check_for_read_in(bool read, NestedArray& a) const {
        a.at(fragment).at(read).at(alignment); 
    }
    
    // Allow comparison
    inline bool operator==(const read_alignment_index_t& other) {
        return fragment == other.fragment && alignment == other.alignment;
    };
    inline bool operator!=(const read_alignment_index_t& other) {
        return !(*this == other);
    };
};
/// Represents an unset index
const read_alignment_index_t NO_READ_INDEX = {std::numeric_limits<size_t>::infinity(), std::numeric_limits<size_t>::infinity()};

/// Type to point to an alignment of either read
/// <fragment index, alignment_index, read number (0 or 1)>
struct alignment_index_t {
    size_t fragment;
    size_t alignment;
    bool read;
    
    /// Drop the read field
    inline read_alignment_index_t without_read() {
        return read_alignment_index_t {fragment, alignment};
    }
    
    /// Allow looking up this index in a structure organized by fragment,
    /// read, and alignment, such as in alignments, or alignment_indices.
    template<typename NestedArray>
    auto lookup_in(NestedArray& a) const -> typename std::add_lvalue_reference<decltype(a[0][0][0])>::type {
        return a[fragment][read][alignment]; 
    }
    
    // Allow comparison
    inline bool operator==(const alignment_index_t& other) {
        return fragment == other.fragment && alignment == other.alignment && read == other.read;
    };
    inline bool operator!=(const alignment_index_t& other) {
        return !(*this == other);
    };
};
/// Represents an unset index
const alignment_index_t NO_INDEX {std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), std::numeric_limits<bool>::max()};

pair<vector<Alignment>, vector<Alignment>> MinimizerMapper::map_paired(Alignment& aln1, Alignment& aln2) {
    
    if (show_work) {
        #pragma omp critical (cerr)
        dump_debug_query(aln1, aln2);
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
        std::array<vector<Alignment>, 2> mapped_pair = {map(aln1), map(aln2)};
        pair_all(mapped_pair);
        return {std::move(mapped_pair[0]), std::move(mapped_pair[1])};
    }


    // Assume reads are in inward orientations on input, and
    // convert to rightward orientations before mapping
    // and flip the second read back before output 

    aln2.clear_path();
    reverse_complement_alignment_in_place(&aln2, [&](vg::id_t node_id) {
        return gbwt_graph.get_length(gbwt_graph.get_handle(node_id));
    });
    
    // Lay out the alignments for looping
    std::array<Alignment*, 2> alns{&aln1, &aln2};

    // Make two new funnel instrumenters to watch us map this read pair.
    std::array<Funnel, 2> funnels;
    // Start this alignment 
    for (auto r : {0, 1}) {
        funnels[r].start(alns[r]->name());
    }
    
    // Annotate the original read with metadata
    for (auto r : {0, 1}) {
        if (!sample_name.empty()) {
            alns[r]->set_sample_name(sample_name);
        }
        if (!read_group.empty()) {
            alns[r]->set_read_group(read_group);
        }
    }
    
    // Prepare the RNG for shuffling ties, if needed
    LazyRNG rng([&]() {
        return aln1.sequence() + aln2.sequence();
    });
    
    // Minimizers for both reads, sorted by read position.
    std::array<std::vector<Minimizer>, 2> minimizers_in_read_by_read;
    // Indexes of minimizers for both reads, sorted into score order, best score first
    std::array<std::vector<size_t>, 2> minimizer_score_order_by_read;
    // Minimizers for both reads, sorted by best score first.
    std::array<VectorView<Minimizer>, 2> minimizers_by_read;
    for (auto r : {0, 1}) {
        minimizers_in_read_by_read[r] = this->find_minimizers(alns[r]->sequence(), funnels[r]);
        minimizer_score_order_by_read[r] = sort_minimizers_by_score(minimizers_in_read_by_read[r]);
        minimizers_by_read[r] = {minimizers_in_read_by_read[r], minimizer_score_order_by_read[r]};
    }

    // Seeds for both reads, stored in separate vectors.
    // These *MUST* be std::vector, because the clusterer's internal data
    // structures pass around pointers to std::vector<std::vector<seed type>>.
    // TODO: Let the clusterer use something else?
    std::vector<std::vector<Seed>> seeds_by_read(2);
    for (auto r : {0, 1}) {
        seeds_by_read[r] = this->find_seeds(minimizers_in_read_by_read[r], minimizers_by_read[r], *alns[r], funnels[r]);
    }

    // Cluster the seeds. Get sets of input seed indexes that go together.
    if (track_provenance) {
        for (auto r : {0, 1}) {
            funnels[r].stage("cluster");
        }
    }

    std::vector<std::vector<Cluster>> all_clusters = clusterer.cluster_seeds(seeds_by_read, get_distance_limit(aln1.sequence().size()), fragment_distance_limit);
#ifdef debug_validate_clusters
    validate_clusters(all_clusters, seeds_by_read, get_distance_limit(aln1.sequence().size()), fragment_distance_limit);

#endif


    //Keep track of which fragment clusters (clusters of clusters) have read clusters from each end

    size_t max_fragment_num = 0;
    for (auto r : {0, 1}) {
        for (auto& cluster : all_clusters[r]) {
            max_fragment_num = std::max(max_fragment_num, cluster.fragment);
        }
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
        for (auto r : {0, 1}) {
            funnels[r].substage("score");
        }
    }


    //Keep track of the best cluster score and coverage per end for each fragment cluster
    std::array<vector<double>, 2> cluster_score_by_fragment;
    std::array<vector<double>, 2> cluster_coverage_by_fragment;
    for (auto r : {0, 1}) {
        cluster_score_by_fragment[r].resize(max_fragment_num + 1, 0.0);
        cluster_coverage_by_fragment[r].resize(max_fragment_num + 1, 0.0);
    }
    
    // We have accessors for totals across both reads, which we use a lot later.
    // TODO: Make a macro?
    auto total_score = [&](size_t fragment_num) -> double {
        return cluster_score_by_fragment[0][fragment_num] + cluster_score_by_fragment[1][fragment_num];
    };
    auto total_coverage = [&](size_t fragment_num) -> double {
        return cluster_coverage_by_fragment[0][fragment_num] + cluster_coverage_by_fragment[1][fragment_num];
    };
    
    //Get the scores of each cluster
    for (auto r : {0, 1}) {
        auto& aln = *alns[r];
        std::vector<Cluster>& clusters = all_clusters[r];
        const VectorView<Minimizer>& minimizers = minimizers_by_read[r];
        vector<double>& best_cluster_score = cluster_score_by_fragment[r];
        vector<double>& best_cluster_coverage = cluster_coverage_by_fragment[r];

        for (size_t i = 0; i < clusters.size(); i++) {
            // Determine cluster score and read coverage.
            Cluster& cluster = clusters[i];
            
            if (this->track_provenance) {
                // Say we're making it
                funnels[r].producing_output(i);
            }
            this->score_cluster(cluster, i, minimizers, seeds_by_read[r], aln.sequence().length());
            if (this->track_provenance) {
                // Record the cluster in the funnel as a group of the size of the number of items.
                funnels[r].merge_group(cluster.seeds.begin(), cluster.seeds.end());
                funnels[r].score(funnels[r].latest(), cluster.score);

                // Say we made it.
                funnels[r].produced_output();
            }
            size_t fragment = cluster.fragment;
            best_cluster_score[fragment] = std::max(best_cluster_score[fragment], cluster.score);
            best_cluster_coverage[fragment] = std::max(best_cluster_coverage[fragment], cluster.coverage);
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
        return total_coverage(a) + total_score(a) > total_coverage(b) + total_score(b);  
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
            double curr_score_sum =  total_coverage(fragment_num) + total_score(fragment_num);
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
        crash_unless(count >= 1);
    }
#endif

    // To compute the windows that are explored, we need to get
    // all the minimizers that are explored.
    std::array<SmallBitset, 2> minimizer_explored_by_read;
    std::array<vector<size_t>, 2> minimizer_aligned_count_by_read;
    //How many hits of each minimizer ended up in each cluster that was kept?
    std::array<vector<vector<size_t>>, 2> minimizer_kept_cluster_count_by_read; 
    
    // To compute the windows present in any extended cluster, we need to get
    // all the minimizers in any extended cluster.
 
    //For each fragment cluster (cluster of clusters), for each read, a vector of all alignments + the order they were fed into the funnel 
    //so the funnel can track them
    vector<std::array<vector<Alignment>, 2>> alignments;
    vector<std::array<vector<size_t>, 2>> alignment_indices;
    std::array<int, 2> best_alignment_scores {0, 0}; // The best alignment score for each end   

    // We will fill this with all computed alignments in estimated score order.
    // alignments has one entry for each fragment cluster and an extra for unpaired alignment 
    alignments.resize(max_fragment_num + 2);
    alignment_indices.resize(max_fragment_num + 2);

    //Now that we've scored each of the clusters, extend and align them
    for (size_t read_num = 0 ; read_num < 2 ; read_num++) {
        Alignment& aln = *alns[read_num];
        std::vector<Cluster>& clusters = all_clusters[read_num];
        const VectorView<Minimizer>& minimizers = minimizers_by_read[read_num];
        std::vector<Seed>& seeds = seeds_by_read[read_num];

        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "Found " << clusters.size() << " clusters for read " << (read_num + 1) << endl;
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
        process_until_threshold_c<double>(clusters.size(), [&](size_t i) -> double {
                return clusters[i].coverage;
            }, [&](size_t a, size_t b) -> bool {
                //Sort clusters first by whether it was paired, then by the best coverage and score of any pair in the fragment cluster, 
                //then by its coverage and score
                size_t fragment_a = clusters[a].fragment;
                size_t fragment_b = clusters[b].fragment;

                double coverage_a = cluster_coverage_by_fragment[0][fragment_a]+cluster_coverage_by_fragment[1][fragment_a];
                double coverage_b = cluster_coverage_by_fragment[0][fragment_b]+cluster_coverage_by_fragment[1][fragment_b];
                double score_a = cluster_score_by_fragment[0][fragment_a]+cluster_score_by_fragment[1][fragment_a];
                double score_b = cluster_score_by_fragment[0][fragment_b]+cluster_score_by_fragment[1][fragment_b];

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
            0, min_extensions, max_extensions, rng, [&](size_t cluster_num) -> bool {
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
                    }

                    if (show_work) {
                        #pragma omp critical (cerr)
                        {
                            cerr << log_name() << "Cluster " << cluster_num << endl;
                        }
                    }
                    
                    // Extend seed hits in the cluster into one or more gapless extensions
                    cluster_extensions.emplace_back(std::move(this->extend_cluster(
                        cluster,
                        cluster_num,
                        minimizers,
                        seeds,
                        aln.sequence(),
                        minimizer_kept_cluster_count_by_read[read_num],
                        funnels[read_num])), cluster.fragment);
                    
                    kept_cluster_count ++;
                    
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
                
            }, [&](size_t cluster_num) -> void {
                // There are too many sufficiently good clusters
                if (track_provenance) {
                    funnels[read_num].pass("cluster-coverage", cluster_num, clusters[cluster_num].coverage);
                    funnels[read_num].fail("max-extensions", cluster_num);
                }
            }, [&](size_t cluster_num) -> void {
                // This cluster is not sufficiently good.
                // TODO: I don't think it should ever get here unless we limit the scores of the fragment clusters we look at
            });

        // We now estimate the best possible alignment score for each cluster.
        std::vector<int> cluster_alignment_score_estimates = this->score_extensions(cluster_extensions, aln, funnels[read_num]);
        
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

        // Go through the processed clusters in estimated-score order.
        process_until_threshold_b(cluster_alignment_score_estimates,
            extension_set_score_threshold, 2, max_alignments, rng, [&](size_t processed_num) {
                // This processed cluster is good enough.
                // Called in descending score order.
                
                if (track_provenance) {
                    funnels[read_num].pass("extension-set", processed_num, cluster_alignment_score_estimates[processed_num]);
                    funnels[read_num].pass("max-alignments", processed_num);
                    funnels[read_num].processing_input(processed_num);
                }
                
                auto& extensions = cluster_extensions[processed_num].first;
                
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
                            cerr << log_name() << "Produced alignment directly from full length gapless extension " << processed_num << endl;
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
                
                
                size_t fragment_num = cluster_extensions[processed_num].second;
                    
                // Have a function to process the best alignments we obtained
                auto observe_alignment = [&](Alignment& aln) {
                    auto& best_score = best_alignment_scores[read_num];
                    best_score = max(best_score, aln.score());
                    
                    auto& alignment_list = alignments[fragment_num][read_num]; 
                    alignment_list.emplace_back(std::move(aln));
                    
                    auto& indices_list = alignment_indices[fragment_num][read_num];
                    indices_list.emplace_back(curr_funnel_index);
                    curr_funnel_index++;

                    if (track_provenance) {
                        funnels[read_num].project(processed_num);
                        funnels[read_num].score(funnels[read_num].latest(), alignment_list.back().score());
                    }
                    
                    if (show_work) {
                        #pragma omp critical (cerr)
                        {
                            cerr << log_name() << "Produced fragment option " << fragment_num << " read " << (read_num + 1)
                                << " alignment with score " << alignment_list.back().score() << ": " << log_alignment(alignment_list.back()) << endl;
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
                
                for (size_t i = 0 ; i < minimizer_kept_cluster_count_by_read[read_num][processed_num].size() ; i++) {
                    if (minimizer_kept_cluster_count_by_read[read_num][processed_num][i] > 0) {
                        // This minimizer is in a cluster that gave rise
                        // to at least one alignment, so it is explored.
                        minimizer_explored_by_read[read_num].insert(i);
                        minimizer_aligned_count_by_read[read_num][i] += minimizer_kept_cluster_count_by_read[read_num][processed_num][i];
                    }
                }

                
                return true;
            }, [&](size_t processed_num) {
                // There are too many sufficiently good processed clusters
                if (track_provenance) {
                    funnels[read_num].pass("extension-set", processed_num, cluster_alignment_score_estimates[processed_num]);
                    funnels[read_num].fail("max-alignments", processed_num);
                }
            }, [&](size_t processed_num) {
                // This processed cluster is not good enough.
                if (track_provenance) {
                    funnels[read_num].fail("extension-set", processed_num, cluster_alignment_score_estimates[processed_num]);
                }
            });
        
    }


    //Now that we have alignments, figure out how to pair them up
    
    if (track_provenance) {
        // Now say we are finding the pairs
        for (auto r : {0, 1}) {
            funnels[r].stage("pairing");
        }
    }
    
    // Fill this in with the indexes of pairs of alignments we will output
    // each alignment is stored as <fragment index, alignment index> into alignments
    // fragment_index should be the same for both ends, unless one was rescued
    vector<std::array<read_alignment_index_t, 2>> paired_alignments;
    paired_alignments.reserve(alignments.size());


#ifdef print_minimizer_table
    vector<pair<bool, bool>> alignment_was_rescued;
#endif

    //For each alignment in alignments, which paired_alignment includes it. Follows structure of alignments
    vector<std::array<vector<vector<size_t>>, 2>> alignment_groups(alignments.size());

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
    vector<alignment_index_t> unpaired_alignments;
    std::array<size_t, 2> unpaired_count {0, 0};

    for (size_t fragment_num = 0 ; fragment_num < alignments.size() ; fragment_num ++ ) {
        //Get pairs of plausible alignments
        for (auto r : {0, 1}) {
            alignment_groups[fragment_num][r].resize(alignments[fragment_num][r].size());
        }
        
        std::array<vector<Alignment>, 2>& fragment_alignments = alignments[fragment_num];
        if (!fragment_alignments[0].empty() && ! fragment_alignments[1].empty()) {
            //Only keep pairs of alignments that were in the same fragment cluster
            found_pair = true;
            std::array<size_t, 2> aln_index;
            std::array<Alignment*, 2> alignment;
            std::array<size_t, 2> funnel_index;
            
            // TODO: replace the nested all-combinations loop with something like a cross product?
            for (aln_index[0] = 0 ; aln_index[0] < fragment_alignments[0].size() ; aln_index[0]++)  {
                alignment[0] = &fragment_alignments[0][aln_index[0]];
                funnel_index[0] = alignment_indices[fragment_num][0][aln_index[0]];
                for (aln_index[1] = 0 ; aln_index[1] < fragment_alignments[1].size() ; aln_index[1]++) {
                    alignment[1] = &fragment_alignments[1][aln_index[1]];
                    funnel_index[1] = alignment_indices[fragment_num][1][aln_index[1]];

                    //Get the likelihood of the fragment distance
                    int64_t fragment_distance = distance_between(*alignment[0], *alignment[1]); 
                    double score = score_alignment_pair(*alignment[0], *alignment[1], fragment_distance);
                    
                    for (auto r : {0, 1}) {
                        alignment_groups[fragment_num][r][aln_index[r]].emplace_back(paired_alignments.size());
                    }
                    paired_alignments.emplace_back();
                    for (auto r : {0, 1}) {
                        paired_alignments.back()[r] = read_alignment_index_t {fragment_num, aln_index[r]};
                    }
#ifdef debug_validate_index_references
                    for (auto r : {0, 1}) {
                        // Make sure we refer to things that exist.
                        paired_alignments.back().at(r).check_for_read_in(r, alignments);
                    }
#endif
                    
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
                                 << alignment[0]->score() << " " << alignment[1]->score() << " at distance " << fragment_distance 
                                 << " gets pair score " << score << endl;
                            cerr << log_name() << "Alignment 1: " << log_alignment(*alignment[0]) << endl << "Alignment 2: " << log_alignment(*alignment[1]) << endl;
                        }
                    }

                    if (track_provenance) {
                        for (auto r : {0, 1}) {
                            funnels[r].processing_input(funnel_index[r]);
                            funnels[r].substage("pair-clusters");
                            funnels[r].pass("max-rescue-attempts", funnel_index[r]);
                            funnels[r].project(funnel_index[r]);
                            funnels[r].score(funnels[r].latest(), score);
                            funnels[r].substage_stop();
                            funnels[r].processed_input();
                        }
                    }
                }
            }
        } else {
            // At most one of these is set
            for (auto r : {0, 1}) {
                if (!fragment_alignments[r].empty()) {
                    //If this fragment cluster has only alignments from this
                
                    if (show_work) {
                        #pragma omp critical (cerr)
                        {
                            cerr << log_name() << "Found unpaired alignments from fragment " << fragment_num << " for read " << (r + 1) << endl;
                        }
                    }
                    for (size_t i = 0 ; i < fragment_alignments[r].size() ; i++) {
                        unpaired_alignments.push_back({fragment_num, i, (bool)r});
                        unpaired_count[r]++;
                        if (show_work) {
                            #pragma omp critical (cerr)
                            {
                                cerr << log_name() << "\t" << log_alignment(fragment_alignments[r][i]) << endl;
                            }
                        }
                    }
                }
            }
        }
    }
    std::array<size_t, 2> rescued_count {0, 0};
    std::array<vector<double>, 2> unpaired_scores;
    for (auto r : {0, 1}) {
        unpaired_scores[r].reserve(unpaired_alignments.size());
    }


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
            std::array<alignment_index_t, 2> best_index {NO_INDEX, NO_INDEX};
            std::array<int32_t, 2> best_score {0, 0};

            for (auto& index : unpaired_alignments ) {
                const Alignment& alignment = index.lookup_in(alignments);
                unpaired_scores[index.read].emplace_back(alignment.score());
                if (deterministic_beats(alignment.score(), best_score[index.read], rng)) {
                    best_index[index.read] = index;
                    best_score[index.read] = alignment.score();
                }
            }
            if (max_rescue_attempts == 0 ) { 
                // If we aren't attempting rescue, just return the best alignment from each end.
                // By default, use argument alignments for scratch.
                std::array<Alignment*, 2> best_aln {alns[0], alns[1]};
    
                for (auto r : {0, 1}) {
                    if (best_index[r] != NO_INDEX) {
                        //If there was a best alignment, use it
                        best_aln[r] = &best_index[r].lookup_in(alignments); 
                    } else {
                        //Otherwise return an empty alignment
                        best_aln[r]->clear_refpos();
                        best_aln[r]->clear_path();
                        best_aln[r]->set_score(0);
                        best_aln[r]->set_identity(0);
                        best_aln[r]->set_mapping_quality(0); // TODO: this gets clobbered
                    }
                    set_annotation(*best_aln[r], "unpaired", true);
                }
                
                std::array<vector<Alignment>, 2> paired_mappings;
                for (auto r : {0, 1}) {
                    paired_mappings[r].emplace_back(std::move(*best_aln[r]));
                }
                // Now we may have moved from our arguments so we can't use them anymore.
                
                // Flip second read back to input orientation
                reverse_complement_alignment_in_place(&paired_mappings[1].back(), [&](vg::id_t node_id) {
                    return gbwt_graph.get_length(gbwt_graph.get_handle(node_id));
                });

                for (auto r : {0, 1}) {
                    paired_mappings[r].back().set_mapping_quality(1);

                    // Stop this alignment
                    funnels[r].stop();
                
                    // Annotate with whatever's in the funnel
                    funnels[r].annotate_mapped_alignment(paired_mappings[r].back(), track_correctness);
                }
                
                return {std::move(paired_mappings[0]), std::move(paired_mappings[1])};
            } else if (best_score[0] != 0 && best_score[1] != 0) {
                //We are attempting rescue, but we still want to keep the best alignments as a potential (unpaired) pair
                
                std::array<size_t, 2> funnel_index;
                if (track_provenance) {
                    // Work out what the paired-up alignments are numbered in the funnel.
                    // TODO: can we flatten these lookup paths or change tuples
                    // to structs to be more understandable?
                    for (auto r : {0, 1}) {
                        funnel_index[r] = best_index[r].lookup_in(alignment_indices);
                        funnels[r].processing_input(funnel_index[r]);
                    }
                }
                
                
                std::array<Alignment*, 2> winners;
                paired_alignments.emplace_back();
                for (auto r : {0, 1}) {
                    paired_alignments.back()[r] = best_index[r].without_read();
                    winners[r] = &best_index[r].lookup_in(alignments);
                }
#ifdef debug_validate_index_references
                for (auto r : {0, 1}) {
                    // Make sure we refer to things that exist.
                    paired_alignments.back().at(r).check_for_read_in(r, alignments);
                }
#endif
                
                //Assume the distance between them is infinite
                double pair_score = score_alignment_pair(*winners[0], *winners[1], std::numeric_limits<int64_t>::max());
                paired_scores.emplace_back(pair_score);
                fragment_distances.emplace_back(std::numeric_limits<int64_t>::max());
                better_cluster_count_by_pairs.emplace_back(0);
                pair_types.emplace_back(unpaired);
                
                if (track_provenance) {
                    for (auto r : {0, 1}) {
                        funnels[r].substage("pair-clusters");
                        funnels[r].project(funnel_index[0]);
                        funnels[r].score(funnels[r].latest(), pair_score);
                        funnels[r].substage_stop();
                        funnels[r].processed_input();
                    }
                }
            }
        }

        if (max_rescue_attempts != 0) {
            //Attempt rescue on unpaired alignments if either we didn't find any pairs or if the unpaired alignments are very good

            process_until_threshold_a(unpaired_alignments.size(), (std::function<double(size_t)>) [&](size_t i) -> double{
                return (double) unpaired_alignments.at(i).lookup_in(alignments).score();
            }, 0, 1, max_rescue_attempts, rng, [&](size_t i) {
                auto& index = unpaired_alignments.at(i);
                size_t j = index.lookup_in(alignment_indices);
                if (track_provenance) {
                    funnels[index.read].processing_input(j);
                    funnels[index.read].substage("rescue");
                }
                Alignment& mapped_aln = index.lookup_in(alignments);
                Alignment rescued_aln = *alns[1 - index.read];
                rescued_aln.clear_path();
                
                if (found_pair && (double) mapped_aln.score() < (double) best_alignment_scores[index.read] * paired_rescue_score_limit) {
                    //If we have already found paired clusters and this unpaired alignment is not good enough, do nothing
                    return true;
                }

                //Rescue the alignment
                attempt_rescue(mapped_aln, rescued_aln, minimizers_by_read[1 - index.read], index.read == 0);

                if (rescued_aln.path().mapping_size() != 0) {
                    //If we actually found an alignment

                    int64_t fragment_dist = index.read == 0 ? distance_between(mapped_aln, rescued_aln) 
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
                    read_alignment_index_t mapped_index = index.without_read();
                    read_alignment_index_t rescued_index {alignments.size() - 1, alignments.back()[1 - index.read].size()};
                    alignments.back()[1 - index.read].emplace_back(std::move(rescued_aln));
                    rescued_count[index.read]++;
                    
                    alignment_groups.back()[1 - index.read].emplace_back();
                    std::array<read_alignment_index_t, 2> index_pair;
                    index_pair[index.read] = mapped_index;
                    index_pair[1 - index.read] = rescued_index;
                    
                    paired_alignments.emplace_back(std::move(index_pair));
#ifdef debug_validate_index_references
                    for (auto r : {0, 1}) {
                        // Make sure we refer to things that exist.
                        paired_alignments.back().at(r).check_for_read_in(r, alignments);
                    }
#endif
                    
                    paired_scores.emplace_back(score);
                    fragment_distances.emplace_back(fragment_dist);
                    pair_types.push_back(index.read == 0 ? rescued_from_first : rescued_from_second); 
                    better_cluster_count_by_pairs.emplace_back(better_cluster_count[mapped_index.fragment]);

#ifdef print_minimizer_table
                    alignment_was_rescued.emplace_back(index.read == 1, index.read == 0);
#endif
                    if (track_provenance) {
                        funnels[index.read].pass("max-rescue-attempts", j);
                        funnels[index.read].project(j);
                        funnels[1 - index.read].introduce();
                        for (auto r : {0, 1}) {
                            funnels[r].score(funnels[r].latest(), score);
                        }
                    }
                }
                if (track_provenance) {
                    funnels[index.read].processed_input();
                    funnels[index.read].substage_stop();
                }
                return true;
            }, [&](size_t i) {
                //This alignment is good enough but we already rescued enough
                if (track_provenance) {
                    auto& index = unpaired_alignments.at(i);
                    size_t j = index.lookup_in(alignment_indices);
                    funnels[index.read].fail("max-rescue-attempts", j);
                }
                return;
            }, [&] (size_t i) {
                //This alignment is insufficiently good
                if (track_provenance) {
                    //TODO: Fail something here
                    auto& index = unpaired_alignments.at(i);
                    size_t j = index.lookup_in(alignment_indices);
                }
                return;
            });
        }
    }

    
    
    if (track_provenance) {
        // Now say we are finding the winner(s)
        for (auto r : {0, 1}) {
            funnels[r].stage("winner");
        }
    }

    // Fill this in with the alignments we will output.
    std::array<vector<Alignment>, 2> mappings;
    // Grab all the scores in order for MAPQ computation.
    vector<double> scores;
    std::array<vector<double>, 2> scores_group;
    vector<int64_t> distances;
    vector<PairType> types;
    for (auto r : {0, 1}) {
        mappings[r].reserve(paired_alignments.size());
    }
    scores.reserve(paired_scores.size());
    distances.reserve(fragment_distances.size());
    types.reserve(pair_types.size());

    //For each pair of alignments in mappings, how many equivalent or better fragment clusters were there
    vector<size_t> better_cluster_count_by_mappings;

#ifdef print_minimizer_table
    vector<pair<bool, bool>> mapping_was_rescued;
    vector<std::array<read_alignment_index_t, 2>> pair_indices;
#endif

    process_until_threshold_a(paired_alignments.size(), (std::function<double(size_t)>) [&](size_t i) -> double {
        return paired_scores[i];
    }, 0, 1, max_multimaps, rng, [&](size_t alignment_num) {
        // This alignment makes it
        // Called in score order

        const std::array<read_alignment_index_t, 2>& index_pair = paired_alignments[alignment_num];
        
        
        // Remember the score at its rank
        scores.emplace_back(paired_scores[alignment_num]);
        distances.emplace_back(fragment_distances[alignment_num]);
        types.emplace_back(pair_types[alignment_num]);
        better_cluster_count_by_mappings.emplace_back(better_cluster_count_by_pairs[alignment_num]);
        // Remember the output alignment
        for (auto r : {0, 1}) {
            mappings[r].emplace_back(index_pair[r].lookup_for_read_in(r, alignments));
        }

        if (mappings[0].size() == 1 && found_pair) {
            //If this is the best pair of alignments that we're going to return and we didn't attempt rescue, 
            //get the group scores for mapq

            //Get the scores of this pair
            for (auto r : {0, 1}) {
                scores_group[r].push_back(paired_scores[alignment_num]);
            }

            //The indices (into paired_alignments) of pairs with the same first/second read as this
            std::array<vector<size_t>*, 2> alignment_group;
            for (auto r : {0, 1}) {
                alignment_group[r] = &index_pair[r].lookup_for_read_in(r, alignment_groups);
            }

            for (auto r : {0, 1}) {
                for (size_t other_alignment_num : *alignment_group[r]) {
                    if (other_alignment_num != alignment_num) {
                        scores_group[r].push_back(paired_scores[other_alignment_num]);
                    }
                }
            }
        }

        // Flip second alignment back to input orientation
        reverse_complement_alignment_in_place(&mappings[1].back(), [&](vg::id_t node_id) {
            return gbwt_graph.get_length(gbwt_graph.get_handle(node_id));
        });
        
        if (mappings[0].size() > 1) {
            // Mark pair as secondary alignments
            for (auto r : {0, 1}) {
                mappings[r].back().set_is_secondary(true);
            }
        }

#ifdef print_minimizer_table
        mapping_was_rescued.emplace_back(alignment_was_rescued[alignment_num]);
        pair_indices.push_back(index_pair);
#endif
        
        if (track_provenance) {
            // Tell the funnel
            for (auto r : {0, 1}) {
                funnels[r].pass("max-multimaps", alignment_num);
                funnels[r].project(alignment_num);
                funnels[r].score(funnels[r].latest(), scores.back());
            }
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
        std::array<read_alignment_index_t, 2> index_pair = paired_alignments[alignment_num];
        pair_indices.push_back(index_pair);
#endif       
        if (track_provenance) {
            for (auto r : {0, 1}) {
                funnels[r].fail("max-multimaps", alignment_num);
            }
        }
    }, [&](size_t alignment_num) {
        // This alignment does not have a sufficiently good score
        // Score threshold is 0; this should never happen
        crash_unless(false);
    });

    if (track_provenance) {
        for (auto r : {0, 1}) {
            funnels[r].substage("mapq");
        }
    }

    
    // Compute raw explored caps (with 2.0 scaling, like for single-end) and raw group caps.
    // Non-capping caps stay at infinity.
    std::array<double, 2> mapq_explored_caps {
        std::numeric_limits<float>::infinity(),
        std::numeric_limits<float>::infinity()
    };
    std::array<double, 2> mapq_score_groups {
        std::numeric_limits<float>::infinity(),
        std::numeric_limits<float>::infinity()
    };
    // We also have one fragment_cluster_cap across both ends.
    double fragment_cluster_cap = std::numeric_limits<float>::infinity();
    // And one base uncapped MAPQ
    double uncapped_mapq = 0;
    double new_cluster_cap = numeric_limits<double>::infinity();

    // Store multiplicities, if we fill them in
    vector<double> paired_multiplicities;

    if (mappings[0].empty()) {
        //If we didn't get an alignment, return empty alignments
        for (auto r : {0, 1}) {
            mappings[r].emplace_back(*alns[r]);
        }

        // Flip second read back to input orientation
        reverse_complement_alignment_in_place(&mappings[1].back(), [&](vg::id_t node_id) {
                return gbwt_graph.get_length(gbwt_graph.get_handle(node_id));
                });

        for (auto r : {0, 1}) {
            mappings[r].back().clear_refpos();
            mappings[r].back().clear_path();
            mappings[r].back().set_score(0);
            mappings[r].back().set_identity(0);
            mappings[r].back().set_mapping_quality(0);
        }

#ifdef print_minimizer_table
        mapping_was_rescued.emplace_back(false, false);
        pair_indices.emplace_back(NO_READ_INDEX, NO_READ_INDEX);
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
        std::array<double, 2> estimated_multiplicity_from;
        for (auto r : {0, 1}) {
            estimated_multiplicity_from[r] = unpaired_count[r] > 0 ? (double) unpaired_count[r] / min(rescued_count[r], max_rescue_attempts) : 1.0;
        }
        bool all_rescued = true;
        // TODO: This switch used to fall-through on all cases.
        for (PairType type : types) {
            switch (type){
                case paired:
                    paired_multiplicities.push_back(1.0);
                    all_rescued=false;
                    break;
                case unpaired:
                    paired_multiplicities.push_back(1.0);
                    break;
                case rescued_from_first:
                    paired_multiplicities.push_back(estimated_multiplicity_from[0]);
                    break;
                case rescued_from_second:
                    paired_multiplicities.push_back(estimated_multiplicity_from[1]);
                    break;
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

        for (auto r : {0, 1}) {
            // For each fragment
        
            //If one alignment was duplicated in other pairs, cap the mapq for that alignment at the mapq
            //of the group of duplicated alignments. Always compute this even if not quite sensible.
            mapq_score_groups[r] = get_regular_aligner()->compute_max_mapping_quality(scores_group[r], false);
        
            vector<size_t> explored_minimizers;
            for (size_t i = 0; i < minimizers_by_read[r].size(); i++) {
                if (minimizer_explored_by_read[r].contains(i)) {
                    explored_minimizers.push_back(i);
                }
            }
            // Compute exploration cap on MAPQ. TODO: avoid needing to pass as much stuff along.
            double mapq_explored_cap =  faster_cap(minimizers_by_read[r], explored_minimizers, alns[r]->sequence(), alns[r]->quality());

            mapq_explored_caps[r] = mapq_explored_cap;

            // Remember the caps
            set_annotation(mappings[r].front(), "mapq_explored_cap", mapq_explored_cap);
            set_annotation(mappings[r].front(), "mapq_score_group", mapq_score_groups[r]);
        }
        
        // Have a function to transform interesting cap values to uncapped.
        auto preprocess_cap = [](double cap) {
            return (cap != -numeric_limits<double>::infinity()) ? cap : numeric_limits<double>::infinity();
        };
        
        for (auto r : {0, 1}) {
            // For each fragment

            // Compute the overall cap for just this read, now that the individual cap components are filled in for both reads.
            double escape_bonus = uncapped_mapq < std::numeric_limits<int32_t>::max() ? 1.0 : 2.0;
            double mapq_cap = std::min(fragment_cluster_cap, ((mapq_explored_caps[0] + mapq_explored_caps[1])*escape_bonus) );

            //TODO: How to cap mapq when the reads were unpaired
            if (types.front() == unpaired) {
                //If this pair came from two different fragment cluster, then cap mapq at the mapq
                //from only unpaired alignments of this read
                mapq_cap = std::min(mapq_cap, (double)get_regular_aligner()->compute_max_mapping_quality(unpaired_scores[r], false));
            }
            
            // Find the MAPQ to cap
            double read_mapq = uncapped_mapq;
            
            // Remember the uncapped MAPQ
            set_annotation(mappings[r].front(), "mapq_uncapped", read_mapq);
            // And the cap we actually applied (possibly from the pair partner)
            set_annotation(mappings[r].front(), "mapq_applied_cap", mapq_cap);

            // Apply the cap, and limit to 0-60
            double capped_mapq = min(mapq_cap, read_mapq); 
            if (distances.front() == std::numeric_limits<int64_t>::max()) {
                //If the two reads are not reachable, lower cap
                capped_mapq = capped_mapq / 2.0;
            }
            read_mapq = max(min(capped_mapq, 120.0) / 2.0, 0.0);
            
            // Save the MAPQ
            mappings[r].front().set_mapping_quality(read_mapq);
            
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "MAPQ for read " << (r + 1) << " is " << read_mapq << ", was " << uncapped_mapq
                        << " capped by fragment cluster cap " << fragment_cluster_cap
                        << ", score group cap " << (mapq_score_groups[r] / 2.0)
                        << ", combined explored cap " << ((mapq_explored_caps[0] + mapq_explored_caps[1]) / 2.0)  << endl;
                }
            }
        }
        
        //Annotate top pair with its fragment distance, properly-paired-ness, fragment length distrubution, and secondary scores
        bool properly_paired = distances.front() == std::numeric_limits<int64_t>::max() ? false :
            (std::abs(distances.front()-fragment_length_distr.mean()) <= 6.0*fragment_length_distr.std_dev()) ;
        string distribution = "-I " + to_string(fragment_length_distr.mean()) + " -D " + to_string(fragment_length_distr.std_dev());
        for (auto r : {0, 1}) {
            set_annotation(mappings[r].front(), "fragment_length", distance_to_annotation(distances.front()));
            set_annotation(mappings[r].front(), "proper_pair", properly_paired);
            set_annotation(mappings[r].front(),"fragment_length_distribution", distribution);
            set_annotation(mappings[r].front(),"secondary_scores", scores);
        }
    }
    
    // Make sure pair partners reference each other
    pair_all(mappings);
        
    
    
    for (auto r : {0, 1}) {
        if (track_provenance) {
            funnels[r].substage_stop();
        }
        // Stop this alignment
        funnels[r].stop();
    }
    
    for (auto r : {0, 1}) {
        // Annotate with whatever's in the funnel.
        funnels[r].annotate_mapped_alignment(mappings[r].front(), track_correctness);
    
        if (track_provenance) {
            if (track_correctness) {
                annotate_with_minimizer_statistics(mappings[r].front(), minimizers_by_read[r], seeds_by_read[r], seeds_by_read[r].size(), 0, funnels[r]);
            }
            // Annotate with parameters used for the filters.
            set_annotation(mappings[r].front(), "param_hit-cap", (double) hit_cap);
            set_annotation(mappings[r].front(), "param_hard-hit-cap", (double) hard_hit_cap);
            set_annotation(mappings[r].front(), "param_score-fraction", (double) minimizer_score_fraction);
            set_annotation(mappings[r].front(), "param_max-extensions", (double) max_extensions);
            set_annotation(mappings[r].front(), "param_max-alignments", (double) max_alignments);
            set_annotation(mappings[r].front(), "param_cluster-score", (double) cluster_score_threshold);
            set_annotation(mappings[r].front(), "param_cluster-coverage", (double) cluster_coverage_threshold);
            set_annotation(mappings[r].front(), "param_extension-set", (double) extension_set_score_threshold);
            set_annotation(mappings[r].front(), "param_max-multimaps", (double) max_multimaps);
            set_annotation(mappings[r].front(), "param_max-rescue-attempts", (double) max_rescue_attempts);
        }
    }
 
#ifdef print_minimizer_table

    if (distances.size() == 0) {
        distances.emplace_back(0);
    }
    for (auto r : {0, 1}) {
        cerr << alns[r]->sequence() << "\t";
        for (char c : alns[r]->quality()) {
            cerr << (char)(c+33);
        }
        cerr << "\t" << max_fragment_num << "\t" << mapping_was_rescued[0].first << "\t" << mapping_was_rescued[0].second 
             << "\t" << distances.front();
        for (size_t i = 0 ; i < minimizers_by_read[r].size() ; i++) {
            auto& minimizer = minimizers_by_read[r][i];
            cerr << "\t"
                 << minimizer.value.key.decode(minimizer.length) << "\t"
                 << minimizer.forward_offset() << "\t"
                 << minimizer.agglomeration_start << "\t"
                 << minimizer.agglomeration_length << "\t"
                 << minimizer.hits << "\t"
                 << minimizer_explored_by_read[r].contains(i);
             if (minimizer_explored_by_read[r].contains(i)) {
                 crash_unless(minimizer.hits<=hard_hit_cap) ;
             }
        }
        cerr << "\t" << uncapped_mapq << "\t" << fragment_cluster_cap << "\t" << mapq_score_groups[0] << "\t" 
             << mapq_explored_caps[0] << "\t" << new_cluster_cap << "\t" << mappings[r].front().mapping_quality() << "\t";  
        for (size_t i = 0 ; i < scores.size() ; i++) {
            std::array<read_alignment_index_t, 2>& indices = pair_indices[i];
            std::array<Alignment*, 2> paired_alignments;
            for (auto r2 : {0, 1}) {
                paired_alignments[r2] = &indices[r2].lookup_for_read_in(r2, alignments);
            }

            int64_t dist = distances[i];
            crash_unless(dist == distance_between(paired_alignments[0], paired_alignments[1])); 

            crash_unless(scores[i] == score_alignment_pair(paired_alignments[0], paired_alignments[1], dist));

            double multiplicity = paired_multiplicities.size() == scores.size() ? paired_multiplicities[i] : 1.0;

            cerr << paired_alignments[0].score() << "," 
                 << paired_alignments[1].score() << "," 
                 << multiplicity << "," 
                 << scores[i] << ";";
        }

        
        if (track_correctness) {
            cerr << "\t" << funnels[r].last_correct_stage() << endl;
        } else {
            cerr << "\t?" << endl;
        }
    }
#endif

    if (track_provenance && show_work) {
        // Dump the funnel info graph.
        #pragma omp critical (cerr)
        {
            for (auto r : {0, 1}) {
                funnels[r].to_dot(cerr);
            }
        }
    }

    // Ship out all the aligned alignments
    return {std::move(mappings[0]), std::move(mappings[1])};
}

//-----------------------------------------------------------------------------

double MinimizerMapper::faster_cap(const VectorView<Minimizer>& minimizers, vector<size_t>& minimizers_explored,
    const string& sequence, const string& quality_bytes) {

    // TODO: Maybe we should to something smarter if we do not have base qualities.
    if (quality_bytes.empty()) {
        return numeric_limits<double>::infinity();
    }

    // Sort minimizer subset so we go through minimizers in increasing order of
    // agglomeration end position, and then increasing order of agglomeration
    // start position
    std::sort(minimizers_explored.begin(), minimizers_explored.end(), [&](size_t a, size_t b) {
        // Return true if a must come before b, and false otherwise
        size_t a_end = minimizers[a].agglomeration_start + minimizers[a].agglomeration_length;
        size_t b_end = minimizers[b].agglomeration_start + minimizers[b].agglomeration_length;
        return a_end < b_end || (a_end == b_end && minimizers[a].agglomeration_start < minimizers[b].agglomeration_start);
    });
#ifdef debug
    cerr << "Sorted " << minimizers_explored.size() << " minimizers" << endl;
#endif

#ifdef debug
    cerr << "Explored minimizers:" << endl;
    dump_debug_minimizers(minimizers, sequence, &minimizers_explored);
#endif

    for (auto it = minimizers_explored.begin(); it != minimizers_explored.end(); ++it) {
        if (minimizers[*it].length == 0) {
            #pragma omp critical (cerr)
            {
                std::cerr << "error[MinimizerMapper::faster_cap]: Minimizer with no sequence found in read with sequence " << sequence << std::endl;
                dump_debug_minimizers(minimizers, sequence, &minimizers_explored);
                for (size_t i = 0 ; i < minimizers_explored.size() ; i++) {
                    auto& m = minimizers[minimizers_explored[i]];
                    std::cerr << "Mininizer " << minimizers_explored[i] << " agg start " << m.agglomeration_start << " length " << m.agglomeration_length
                              << " core start " << m.value.offset << " length " << m.length << std::endl;
                }
                std::cerr << "Read sequence: " << sequence << std::endl;
                std::cerr << "Read quality: ";
                for (char q : quality_bytes) {
                    std::cerr << (char) (33 + (int)q);
                }
                std::cerr << std::endl;
                exit(1);
            }
        }
    }

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

        if (isinf(p_here)) {
            #pragma omp critical (cerr)
            {
                std::cerr << "error[MinimizerMapper::faster_cap]: Minimizers seem impossible to disrupt in region " << left << " " << right << " " << bottom << " " << top << std::endl;
                dump_debug_minimizers(minimizers, sequence, &minimizers_explored);
                for (size_t i = 0 ; i < minimizers_explored.size() ; i++) {
                    auto& m = minimizers[minimizers_explored[i]];
                    std::cerr << "Mininizer " << minimizers_explored[i] << " agg start " << m.agglomeration_start << " length " << m.agglomeration_length
                              << " core start " << m.value.offset << " length " << m.length << std::endl;
                }
                std::cerr << "Read sequence: " << sequence << std::endl;
                std::cerr << "Read quality: ";
                for (char q : quality_bytes) {
                    std::cerr << (char) (33 + (int)q);
                }
                std::cerr << std::endl;
            }
            exit(1);
        }
        
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
    
    if (isinf(c.back())) {
        #pragma omp critical (cerr)
        {
            std::cerr << "error[MinimizerMapper::faster_cap]: Minimizers seem impossible to disrupt!" << std::endl;
            dump_debug_minimizers(minimizers, sequence, &minimizers_explored);
            for (size_t i = 0 ; i < minimizers_explored.size() ; i++) {
                auto& m = minimizers[minimizers_explored[i]];
                std::cerr << "Mininizer " << minimizers_explored[i] << " agg start " << m.agglomeration_start << " length " << m.agglomeration_length
                          << " core start " << m.value.offset << " length " << m.length << std::endl;
            }
            std::cerr << "Read sequence: " << sequence << std::endl;
            std::cerr << "Read quality: ";
            for (char q : quality_bytes) {
                std::cerr << (char) (33 + (int)q);
            }
            std::cerr << std::endl;
        }
        exit(1);
    }
    
    // Conver to Phred.
    double result = -c.back() * 10;
    return result;
}

void MinimizerMapper::for_each_agglomeration_interval(const VectorView<Minimizer>& minimizers,
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
                
                if (stack_top_end < left) {
                    // Something is wrong with the order we are visiting these in.
                    #pragma omp critical (cerr)
                    {
                        std::cerr << "error[MinimizerMapper::faster_cap]: Minimizers not sorted properly for read with sequence " << sequence << "! Agglomeration on stack ends at " << stack_top_end << " but we are at " << left << " from a previous agglomeration" << std::endl;
                        exit(1);
                    }
                }
                
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
        
        if (stack.size() == 0) {
            // Something is wrong with our stacking algorithm
            #pragma omp critical (cerr)
            {
                std::cerr << "error[MinimizerMapper::faster_cap]: Minimizers not stacked up properly for read with sequence " << sequence << "!" << std::endl;
                exit(1);
            }
        }

        // For each new item we return all intervals that
        // precede its start
        emit_preceding_intervals(item.agglomeration_start);

        // Add the new item for the next loop
        stack.push_back(&item);
    }

    // Intervals of the remaining intervals on the stack
    emit_preceding_intervals(sequence.size());
}

double MinimizerMapper::get_log10_prob_of_disruption_in_interval(const VectorView<Minimizer>& minimizers,
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

double MinimizerMapper::get_prob_of_disruption_in_column(const VectorView<Minimizer>& minimizers,
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

#ifdef debug
            cerr << "\t\t\tBeat hash " << m.value.hash << " at least 1 time in " << possible_minimizers << endl;
#endif

            // Account for at least one of them beating the minimizer.
            double any_beat_prob = prob_for_at_least_one(m.value.hash, possible_minimizers);
            
#ifdef debug
            cerr << "\t\t\t\tGives probability: " << any_beat_prob << endl;
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

void MinimizerMapper::attempt_rescue(const Alignment& aligned_read, Alignment& rescued_alignment, const VectorView<Minimizer>& minimizers, bool rescue_forward ) {

    // Get rid of the old path.
    rescued_alignment.clear_path();

    if (this->rescue_algorithm == rescue_none) { return; }

    // We are traversing the same small subgraph repeatedly, so it's better to use a cache.
    gbwtgraph::CachedGBWTGraph cached_graph(this->gbwt_graph);

    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "Attempt rescue from: " << log_alignment(aligned_read) << endl;
        }
    }

    // Find all nodes within a reasonable range from aligned_read.
    std::unordered_set<id_t> rescue_nodes;
    int64_t min_distance = max(0.0, fragment_length_distr.mean() - rescued_alignment.sequence().size() - rescue_subgraph_stdevs * fragment_length_distr.std_dev());
    int64_t max_distance = fragment_length_distr.mean() + rescue_subgraph_stdevs * fragment_length_distr.std_dev();

    subgraph_in_distance_range(*distance_index, aligned_read.path(), &cached_graph, min_distance, max_distance, rescue_nodes, rescue_forward);

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
            cerr << log_name() << "Rescue result: " << log_alignment(rescued_alignment) << endl;
        }
    }
}

GaplessExtender::cluster_type MinimizerMapper::seeds_in_subgraph(const VectorView<Minimizer>& minimizers,
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


int64_t MinimizerMapper::distance_between(const pos_t& pos1, const pos_t& pos2) {
    int64_t min_dist= distance_index->minimum_distance(id(pos1), is_rev(pos1), offset(pos1),
                                                       id(pos2), is_rev(pos2), offset(pos2),
                                                       false, &gbwt_graph);
    return min_dist;
}

int64_t MinimizerMapper::unoriented_distance_between(const pos_t& pos1, const pos_t& pos2) const {

    int64_t min_dist = distance_index->minimum_distance(id(pos1), is_rev(pos1), offset(pos1),
                                                        id(pos2), is_rev(pos2), offset(pos2),
                                                        true, &gbwt_graph);
    return min_dist;
}

int64_t MinimizerMapper::distance_between(const Alignment& aln1, const Alignment& aln2) {
    crash_unless(aln1.path().mapping_size() != 0); 
    crash_unless(aln2.path().mapping_size() != 0); 
     
    pos_t pos1 = initial_position(aln1.path()); 
    pos_t pos2 = final_position(aln2.path());

    return distance_between(pos1, pos2);
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
    
    if (this->track_provenance) {
        // Record how many we found, as new lines.
        // THey are going to be numbered in score order, not read order. Probably...
        funnel.introduce(result.size());
    }

    return result;
}

std::vector<size_t> MinimizerMapper::sort_minimizers_by_score(const std::vector<Minimizer>& minimizers) const {
    // We defined operator< so the minimizers always sort descening by score by default.
    return sort_permutation(minimizers.begin(), minimizers.end());
}

std::vector<MinimizerMapper::Seed> MinimizerMapper::find_seeds(const std::vector<Minimizer>& minimizers_in_read_order, const VectorView<Minimizer>& minimizers, const Alignment& aln,  Funnel& funnel) const {

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
    size_t start = 0, limit = 0;
    size_t run_hits = 0;
    bool taking_run = false;

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
    size_t num_min_by_read_len = read_len / this->num_bp_per_min;
    std::vector<bool> read_bit_vector (read_len, false);

    // Select the minimizers we use for seeds.
    size_t rejected_count = 0;
    std::vector<Seed> seeds;

    // Prefilter and downsample the minimizers with a sliding window.
    // We do all this in *read* order!
    auto minimizer_worse_than = [&](const Minimizer& min_a, const Minimizer& min_b) {
        // Return true if minimizer a is worse than minimizer b.
        // The worse minimizer is the one that doesn't match the reference, or
        // if both match the reference it is the one that has less score.
        return (min_a.hits == 0 && min_b.hits > 0) || (min_a.hits != 0 && min_b.hits != 0 && min_a.score < min_b.score);
    };
    // We keep a set of the minimizers that pass downsampling.
    // We later need to filter given a minimizer reference and that makes it hard to use a bit vector here.
    // TODO: change how the filters work!
    std::unordered_set<const Minimizer*> downsampled;
    if (this->minimizer_downsampling_window_size != 0) {
        // This will hold all the minimizers in the sliding window of bases
        std::deque<size_t> queue;
        for (size_t i = 0; i < minimizers_in_read_order.size(); i++) {
            // We use the minimizer sampling algorithm again, as described in the Winnowmap paper (Jain et al., 2020).
            
            auto& new_minimizer = minimizers_in_read_order[i];
            size_t new_window_end = new_minimizer.forward_offset() + new_minimizer.length;

            while (!queue.empty() && minimizer_worse_than(minimizers_in_read_order.at(queue.back()), new_minimizer)) { 
                // Drop minimizers off the end of the queue until it is empty
                // or we find one that is at least as good as the new
                // minimizer.
                queue.pop_back();
            }

            // Add the new minimizer.
            queue.push_back(i);

            while(!queue.empty() && minimizers_in_read_order[queue.front()].forward_offset() + this->minimizer_downsampling_window_size < new_window_end) {
                // Drop minimizers off the front of the queue until it is empty
                // or we find one that is in-window.
                queue.pop_front();
            }
            
            if (queue.empty()) {
                // We removed the minimizer we just added. The window is probably too small.
                #pragma omp critical (cerr)
                std::cerr << "error:[vg::MinimizerMapper] no minimizer found in downsampling window. Make sure that the downsampling window is at least " << new_minimizer.length << " bp" << std::endl;
                exit(1);
            }

            // Since we never add a better minimizer after a worse one, the first thing in the queue is the best minimizer in the window.
            downsampled.insert(&minimizers_in_read_order[queue.front()]);
        }
        if (show_work) {
            #pragma omp critical (cerr)
            std::cerr << log_name() << "Downsampled to " << downsampled.size() << " minimizers" << std::endl;
        }
    }
    
    // Define the filters for minimizers.
    //
    // Each has a name, a function to check if a minimizer passes, a function
    // to check the minimizer's statistic, and a pass and a fail function to
    // update state.
    //
    // Filters that keep state and don't want to treat multiple duplicate
    // minimizers differently should check taking_run to see if something in
    // the run already passed all filters (and therefore this filter).
    //
    // Functions are only allowed to be called after run_hits is set and the
    // run logic has determined run membership for the minimizer.
    using filter_t = std::tuple<const char*, std::function<bool(const Minimizer&)>, std::function<double(const Minimizer&)>, std::function<void(const Minimizer&)>, std::function<void(const Minimizer&)>>;
    std::vector<filter_t> minimizer_filters;
    minimizer_filters.reserve(5);
    if (this->minimizer_downsampling_window_size != 0) {
        // Drop minimizers if we cleared their downsampling flag. Sneakily go back from minimizer itself to index in the array.
        minimizer_filters.emplace_back(
            "window-downsampling", 
            [&](const Minimizer& m) { return downsampled.empty() || downsampled.count(&m); },
            [](const Minimizer& m) { return nan(""); },
            [](const Minimizer& m) {},
            [](const Minimizer& m) {}
        );
    }
    minimizer_filters.emplace_back(
        "any-hits", 
        [&](const Minimizer& m) { return m.hits > 0; },
        [](const Minimizer& m) { return nan(""); },
        [](const Minimizer& m) {},
        [](const Minimizer& m) {}
    );
    minimizer_filters.emplace_back(
        "hard-hit-cap",
        [&](const Minimizer& m) { return run_hits <= this->hard_hit_cap; },
        [&](const Minimizer& m) { return (double)run_hits; },
        [](const Minimizer& m) {},
        [](const Minimizer& m) {}
    );
    if (this->exclude_overlapping_min) {
        minimizer_filters.emplace_back(
            "exclude-overlapping-min",
            [&](const Minimizer& m) {
                // TODO: Is this matching abutting minimizers?
                // If at both ehnds the flag isn't set yet, we pass. 
                return !read_bit_vector[m.forward_offset()] && 
                !read_bit_vector[m.forward_offset() + m.length];
            },
            [](const Minimizer& m) { return nan(""); },
            [&](const Minimizer& m) {
              for (size_t i = m.forward_offset(); i < m.forward_offset() + m.length; i++) {
                read_bit_vector[i] = true;
              }
            },
            [](const Minimizer& m) {}
        );
    }
    if (this->max_unique_min != 0) {
        minimizer_filters.emplace_back(
            "max-unique-min||num-bp-per-min",
            [&](const Minimizer& m) {
                return num_minimizers < std::max(this->max_unique_min, num_min_by_read_len);
            },
            [](const Minimizer& m) { return nan(""); },
            [](const Minimizer& m) {},
            [](const Minimizer& m) {}
        );
    }
    minimizer_filters.emplace_back(
        "hit-cap||score-fraction",
        [&](const Minimizer& m) {
            return (m.hits <= this->hit_cap) || // We pass if we are under the soft hit cap
            (run_hits <= this->hard_hit_cap && selected_score + m.score <= target_score) || // Or the run as a whole is under the hard hit cap and we need the score
            (taking_run); // Or we already took one duplicate and we want to finish out the run 
        },
        [&](const Minimizer& m) {
            return (selected_score + m.score) / base_target_score;
        },
        [&](const Minimizer& m) {
            // Remember that we took this minimizer for evaluating later ones
            selected_score += m.score;
        },
        [&](const Minimizer& m) {
            //Stop looking for more minimizers once we fail the score fraction
            target_score = selected_score; 
        }
    );
     
    
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
            // We are starting a new run
            start = i; limit = i + 1;
            run_hits = minimizers[i].hits;
            for (size_t j = i + 1; j < minimizers.size() && minimizers[j].value.key == minimizers[i].value.key; j++) {
                limit++;
                run_hits += minimizers[j].hits;
            }
            // We haven't taken the first thing in the run yet.
            taking_run = false;
        }

        // Select the minimizer if it is informative enough or if the total score
        // of the selected minimizers is not high enough.
        const Minimizer& minimizer = minimizers[i];

        // Evaluate the filters
        bool passing = true;
        for (auto& filter : minimizer_filters) {
            auto& filter_name = std::get<0>(filter);
            auto& filter_function = std::get<1>(filter);
            auto& filter_stat_function = std::get<2>(filter);
            auto& filter_pass_function = std::get<3>(filter);
            auto& filter_fail_function = std::get<4>(filter);
            
            passing = filter_function(minimizer);
            if (passing) {
                // Pass this filter
                if (this->track_provenance) {
                    funnel.pass(filter_name, i, filter_stat_function(minimizer));
                }
                filter_pass_function(minimizer);
            } else {
                // Fail this filter.
                if (this->track_provenance) {
                    funnel.fail(filter_name, i, filter_stat_function(minimizer));
                }
                filter_fail_function(minimizer);
                // Don't do later filters
                break;
            }
        }
        
        if (passing) {
            // We passed all filters.
            // So we are taking this item and ought to take the others in the same run in most cases.
            taking_run = true;
            // Track number of minimizers selected.
            num_minimizers++;
            
            // We should keep this minimizer instance because it is
            // sufficiently rare, or we want it to make target_score, or it is
            // the same sequence as a previous minimizer in this run of identical
            // minimizers which we also took.

            // Locate the hits.
            for (size_t j = 0; j < minimizer.hits; j++) {
                pos_t hit = gbwtgraph::Position::decode(minimizer.occs[j].pos);
                // Reverse the hits for a reverse minimizer
                if (minimizer.value.is_reverse) {
                    size_t node_length = this->gbwt_graph.get_length(this->gbwt_graph.get_handle(id(hit)));
                    hit = reverse_base_pos(hit, node_length);
                }
                // Extract component id and offset in the root chain, if we have them for this seed.
                seeds.emplace_back();
                seeds.back().pos = hit;
                seeds.back().source = i;

                //Get the zipcode
                if (minimizer.occs[j].payload == MIPayload::NO_CODE) {
                    //If the zipcocde wasn't saved, then calculate it
                    seeds.back().zipcode.fill_in_zipcode(*(this->distance_index), hit);
                } else if (minimizer.occs[j].payload.first == 0) {
                    //If the minimizer stored the index into a list of zipcodes
                    if (this->zipcodes != nullptr) {
                        //If we have the oversized zipcodes
                        seeds.back().zipcode = zipcodes->at(minimizer.occs[j].payload.second);
                    } else {
                        //If we don't have the oversized payloads, then fill in the zipcode using the pos
                        seeds.back().zipcode.fill_in_zipcode(*(this->distance_index), hit);
                    }
                } else {
                    //If the zipcode was saved in the payload
                    seeds.back().zipcode.fill_in_zipcode_from_payload(minimizer.occs[j].payload);
                }
                ZipCodeDecoder* decoder = new ZipCodeDecoder(&seeds.back().zipcode);
                seeds.back().zipcode_decoder.reset(decoder);

            }
            
            if (this->track_provenance) {
                // Record in the funnel that this minimizer gave rise to these seeds.
                funnel.expand(i, minimizer.hits);
            }
                
        } else {
            // We failed a filter.
            // Track number of minimizers rejected.
            rejected_count++;
        }
        
        if (this->track_provenance) {
            // Say we're done with this input item
            funnel.processed_input();
        }
    }

    if (this->track_provenance) {
        if (this->track_correctness) {
            // Tag seeds with correctness 
            funnel.substage("correct");
        } else {
            // We're just tagging them with read positions
            funnel.substage("placed");
        }
        this->tag_seeds(aln, seeds.cbegin(), seeds.cend(), minimizers, 0, funnel);
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

void MinimizerMapper::tag_seeds(const Alignment& aln, const std::vector<Seed>::const_iterator& begin, const std::vector<Seed>::const_iterator& end, const VectorView<Minimizer>& minimizers, size_t funnel_offset, Funnel& funnel) const { 
    if (this->track_correctness && this->path_graph == nullptr) {
        cerr << "error[vg::MinimizerMapper] Cannot use track_correctness with no XG index" << endl;
        exit(1);
    }
    
    // Track the index of each seed in the funnel
    size_t funnel_index = funnel_offset;
    for (std::vector<Seed>::const_iterator it = begin; it != end; ++it) {
        
        // We know the seed is placed somewhere.
        Funnel::State tag = Funnel::State::PLACED;
        if (this->track_correctness && aln.refpos_size() != 0) {
            // It might also be correct
            // Find every seed's reference positions. This maps from path name to pairs of offset and orientation.
            auto offsets = algorithms::nearest_offsets_in_paths(this->path_graph, it->pos, 100);
            
            for (auto& true_pos : aln.refpos()) {
                // For every annotated true position
                for (auto& hit_pos : offsets[this->path_graph->get_path_handle(true_pos.name())]) {
                    // Look at all the hit positions on the path the read's true position is on.
                    if (abs((int64_t)hit_pos.first - (int64_t) true_pos.offset()) < 200) {
                        // We're close enough to be correct
                        tag = Funnel::State::CORRECT;
                        break;
                    }
                }
                if (tag == Funnel::State::CORRECT) {
                    break;
                }
            }
        }
                
        // Tag this seed as making some of the read space placed or even correct.
        funnel.tag(funnel_index, tag, minimizers[it->source].forward_offset(), minimizers[it->source].length);
        
        // Look at the next seed
        funnel_index++;
    }
}

void MinimizerMapper::annotate_with_minimizer_statistics(Alignment& target, const VectorView<Minimizer>& minimizers, const std::vector<Seed>& seeds, size_t old_seed_count, size_t new_seed_offset, const Funnel& funnel) const {
    // Annotate with fraction covered by correct (and necessarily located) seed hits.
    
    // First make the set of minimizers that got correct seeds
    std::unordered_set<size_t> seeded;
    for (size_t i = 0; i < seeds.size(); i++) {
        // We need to get correctness out of the funnel, since we don't tag the seed or minimizer.
        // Correctness is assessed per seed, not per minimizer.
        
        if (i < old_seed_count) {
            // This seed was created at the seed stage.
            // We know seed finding was always stage 1.
            if (funnel.was_correct(1, "seed", i)) {
                seeded.insert(seeds[i].source);
            }
        } else {
            // This seed must have been created at the reseed stage in the
            // align_from_chains codepath. We happen to know the magic stager
            // number for it.
            // Make sure to translate seed number to funnel number in the reseed stage.
            // TODO: This is a tightly coupled hack, do we even need this annotation anymore? And does it really need to include the reseeded seeds?
            if (funnel.was_correct(3, "reseed", i - old_seed_count + new_seed_offset)) {
                seeded.insert(seeds[i].source);
            }
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

void MinimizerMapper::score_cluster(Cluster& cluster, size_t i, const VectorView<Minimizer>& minimizers, const std::vector<Seed>& seeds, size_t seq_length) const {

    // Initialize the values.
    cluster.score = 0.0;
    cluster.coverage = 0.0;
    cluster.present = SmallBitset(minimizers.size());

    // Determine the minimizers that are present in the cluster.
    for (auto hit_index : cluster.seeds) {
        cluster.present.insert(seeds[hit_index].source);
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
}

//-----------------------------------------------------------------------------

vector<GaplessExtension> MinimizerMapper::extend_cluster(const Cluster& cluster,
    size_t cluster_num,
    const VectorView<Minimizer>& minimizers,
    const std::vector<Seed>& seeds,
    const string& sequence,
    vector<vector<size_t>>& minimizer_kept_cluster_count,
    Funnel& funnel) const {

    if (track_provenance) {
        // Say we're working on this cluster
        funnel.processing_input(cluster_num);
    }

    // Count how many of each minimizer is in each cluster that we kept
    minimizer_kept_cluster_count.emplace_back(minimizers.size(), 0);
    // Pack the seeds for GaplessExtender.
    GaplessExtender::cluster_type seed_matchings;
    for (auto seed_index : cluster.seeds) {
        // Insert the (graph position, read offset) pair.
        auto& seed = seeds[seed_index];
        seed_matchings.insert(GaplessExtender::to_seed(seed.pos, minimizers[seed.source].value.offset));
        minimizer_kept_cluster_count.back()[seed.source]++;
        
        if (show_work) {
            #pragma omp critical (cerr)
            {
                dump_debug_seeds(minimizers, seeds, cluster.seeds);
            }
        }
    }
    
    vector<GaplessExtension> cluster_extension = extender.extend(seed_matchings, sequence);

    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "Extensions:" << endl;
            for (auto& e : cluster_extension) {
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
        funnel.project_group(cluster_num, cluster_extension.size());
        // Say we finished with this cluster, for now.
        funnel.processed_input();
    }
    
    return cluster_extension;
}

//-----------------------------------------------------------------------------

int MinimizerMapper::score_extension_group(const Alignment& aln, const vector<GaplessExtension>& extended_seeds,
    int gap_open_penalty, int gap_extend_penalty) {
        
    if (extended_seeds.empty()) {
        // TODO: We should never see an empty group of extensions
        return 0;
    } else if (GaplessExtender::full_length_extensions(extended_seeds)) {
        // These are full-length matches. We already have the score.
        return extended_seeds.front().score;
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

static void find_pareto_frontier(std::vector<pareto_point>& v) {
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
static int32_t gap_penalty(size_t length, const Aligner* aligner) {
    return (length == 0 ? 0 : aligner->gap_open + (length - 1) * aligner->gap_extension);
}

// Positive penalty for a number of mismatches.
static int32_t mismatch_penalty(size_t n, const Aligner* aligner) {
    return n * (aligner->match + aligner->mismatch);
}

// Positive gap penalty, assuming that there is always a gap.
static int32_t gap_penalty(size_t start, size_t limit, const Aligner* aligner) {
    return (start >= limit ? aligner->gap_open : aligner->gap_open + (limit - start - 1) * aligner->gap_extension);
}

// Positive flank penalty based on taking a gap to the end or to the Pareto frontier.
static int32_t flank_penalty(size_t length, const std::vector<pareto_point>& frontier, const Aligner* aligner) {
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
    size_t min_tails = 1;
    for (const GaplessExtension& extension : extended_seeds) {
        if (extension.full()) {
            min_tails++;
        }
    }
    if (min_tails < 2) {
        min_tails = 2;
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
            cerr << log_name() << "Trying to find " << min_tails << " tail alignments for "
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
    process_until_threshold_a<double>(extended_seeds.size(),
        [&](size_t extended_seed_num) -> double {
            return static_cast<double>(extended_seeds[extended_seed_num].score);
        }, extension_score_threshold, min_tails, max_local_extensions, rng, [&](size_t extended_seed_num) -> bool {
       
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
        }, [&](size_t extended_seed_num) -> void {
            // This extended seed is good enough by its own score, but we have too many.
            // Do nothing
        }, [&](size_t extended_seed_num) -> void {
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
                crash_unless(mapping.position().node_id() == best.path().mapping(best.path().mapping_size() - 1).position().node_id());

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
                crash_unless(mapping.position().node_id() == second_best.path().mapping(second_best.path().mapping_size() - 1).position().node_id());

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
                cerr << log_name() << "First best alignment: " << log_alignment(best_path) << " score " << best_score << endl;
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
                    cerr << log_name() << "Align " << log_alignment(current_alignment) << " pinned left" << endl;
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
                            << log_alignment(best_path) << " score " << best_score << endl;
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
    size_t remaining_root = gbwt_graph.get_length(from_handle) - from_offset;
    
#ifdef debug
    cerr << "DFS starting at offset " << from_offset << " on node "
        << gbwt_graph.get_id(from_handle) << " " << gbwt_graph.get_is_reverse(from_handle) << " of length "
        << gbwt_graph.get_length(from_handle) << " leaving " << remaining_root << " bp" << endl;
#endif
    
    // We would do a recursive formulation, but if attempted on long reads that can run out of stack.
    // See <https://github.com/vgteam/vg/issues/3733>.
    // So we go the iterative route. Which means we need a stack frame for our simulated stack.
    struct Frame {
        gbwt::SearchState here_state;
        size_t used_distance;
        bool visit;
    };
    std::stack<Frame> stack;
    
    // Start the DFS with our stating node, consuming the distance from our
    // offset to its end. Don't show the root state to the user if we don't
    // actually visit any bases on that node.
    // Make sure we don't count the length of the root node inside the DFS,
    // since we are already feeding it in.
    stack.push({start_state, 0, false});
    
    while (!stack.empty()) {
        auto& frame = stack.top();
        
        handle_t here_handle = gbwt_graph.node_to_handle(frame.here_state.node);
        
        // If we're at the root node we don't measure it the same way
        bool is_root = (stack.size() == 1);
        // If we're at the root node but there's no bases on it left, we don't show it.
        bool is_hidden_root = (is_root && remaining_root == 0);
        
        if (frame.visit == false) {
            // First visit
            frame.visit = true;
            
            if (!is_hidden_root) {
                // Enter this handle if there are any bases on it to visit
#ifdef debug
                cerr << "Enter handle " << gbwt_graph.get_id(here_handle) << " " << gbwt_graph.get_is_reverse(here_handle) << endl;
#endif
                enter_handle(here_handle);
            }
            
            // Use the length of the node, or the remaining root if we're the root.
            size_t node_length = is_root ? remaining_root : gbwt_graph.get_length(here_handle);
            frame.used_distance += node_length;
#ifdef debug
            cerr << "Node was " << node_length << " bp; Used " << frame.used_distance << "/" << walk_distance << " bp distance" << endl;
#endif
            if (frame.used_distance < walk_distance) {
                // If we haven't used up all our distance yet
                
                // Stack up all the paths onto the stack to be processed.
                gbwt_graph.follow_paths(frame.here_state, [&](const gbwt::SearchState& there_state) -> bool {
                    // For each next state
                    
                    // Do it with the new distance value.
                    stack.push({there_state, frame.used_distance, 0});
                    
                    return true;
                });
                
                // Jump to handling the new top of the stack. When we come back
                // to us, we'll take the second visit branch.
                continue;
            }
        }
        
        // Second visit, or first visit didn't expand the stack.
        
        if (!is_hidden_root) {
            // Exit this handle if we entered it
            
#ifdef debug
            cerr << "Exit handle " << gbwt_graph.get_id(here_handle) << " " << gbwt_graph.get_is_reverse(here_handle) << endl;
#endif
            
            exit_handle();
        }
        
        // Remove ourselves from the stack.
        stack.pop();
    }
}

//-----------------------------------------------------------------------------

double MinimizerMapper::score_alignment_pair(Alignment& aln1, Alignment& aln2, int64_t fragment_distance) {
    //Score a pair of alignments

    double dev = fragment_distance - fragment_length_distr.mean();
    double fragment_length_log_likelihood = (-dev * dev / (2.0 * fragment_length_distr.std_dev() * fragment_length_distr.std_dev()))/ get_regular_aligner()->log_base;
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


