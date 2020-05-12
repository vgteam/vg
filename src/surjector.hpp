#ifndef VG_SURJECTOR_HPP_INCLUDED
#define VG_SURJECTOR_HPP_INCLUDED

/** \file
 *
 *  A class to hold surjection algorithms that do lossy realignment restricted to paths in the graph
 */

#include <set>

#include "alignment.hpp"
#include "aligner.hpp"
#include "vg.hpp"
#include "translator.hpp"
#include "utility.hpp"
#include <vg/vg.pb.h>
#include "multipath_alignment_graph.hpp"
#include "memoizing_graph.hpp"
#include "split_strand_graph.hpp"

#include "algorithms/topological_sort.hpp"

#include "bdsg/hash_graph.hpp"

namespace vg {

using namespace std;

    class Surjector : public AlignerClient {
    public:
        
        Surjector(const PathPositionHandleGraph* graph);
        
        /// Extract the portions of an alignment that are on a chosen set of paths and try to
        /// align realign the portions that are off of the chosen paths to the intervening
        /// path segments to obtain an alignment that is fully restricted to the paths.
        ///
        /// Also returns the path name, position, and strand of the new alignment.
        ///
        /// Optionally either allow softclips so that the alignment has a nonnegative score on
        /// the path or require the full-length alignment, possibly creating a negative score.
        ///
        /// Also optionally leaves deletions against the reference path in the final alignment
        /// (useful for splicing).
        Alignment surject(const Alignment& source,
                          const set<string>& path_names,
                          string& path_name_out,
                          int64_t& path_pos_out,
                          bool& path_rev_out,
                          bool allow_negative_scores = false,
                          bool preserve_deletions = false) const;
                          
        /// Extract the portions of an alignment that are on a chosen set of
        /// paths and try to align realign the portions that are off of the
        /// chosen paths to the intervening path segments to obtain an
        /// alignment that is fully restricted to the paths.
        ///
        /// Replaces the alignment's refpos with the path name, position, and
        /// strand the alignment has been surjected to.
        ///
        /// Optionally either allow softclips so that the alignment has a
        /// nonnegative score on the path or require the full-length alignment,
        /// possibly creating a negative score.
        ///
        /// Also optionally leaves deletions against the reference path in the final
        /// alignment (useful for splicing).
        Alignment surject(const Alignment& source,
                          const set<string>& path_names,
                          bool allow_negative_scores = false,
                          bool preserve_deletions = false) const;
        
        /// a local type that represents a read interval matched to a portion of the alignment path
        using path_chunk_t = pair<pair<string::const_iterator, string::const_iterator>, Path>;
        
    private:
        
        Alignment
        realigning_surject(const PathPositionHandleGraph* graph, const Alignment& source,
                           const path_handle_t& path_handle, const vector<path_chunk_t>& path_chunks,
                           bool allow_negative_scores, bool preserve_N_alignments = false,
                           bool preserve_tail_indel_anchors = false) const;
        
        Alignment
        spliced_surject(const PathPositionHandleGraph* path_position_graph, const Alignment& source,
                        const path_handle_t& path_handle, const vector<path_chunk_t>& path_chunks,
                        const vector<pair<step_handle_t, step_handle_t>>& ref_chunks,
                        bool allow_negative_scores) const;
        
        ///////////////////////
        // Support methods for the realigning surject algorithm
        ///////////////////////
        
        /// get the chunks of the alignment path that follow the given reference paths
        unordered_map<path_handle_t, pair<vector<path_chunk_t>, vector<pair<step_handle_t, step_handle_t>>>>
        extract_overlapping_paths(const PathPositionHandleGraph* graph, const Alignment& source,
                                  const unordered_set<path_handle_t>& surjection_paths) const;
        
        /// compute the widest interval of path positions that the realigned sequence could align to
        pair<size_t, size_t>
        compute_path_interval(const PathPositionHandleGraph* graph, const Alignment& source, path_handle_t path_handle,
                              const vector<path_chunk_t>& path_chunks) const;
        
        /// make a linear graph that corresponds to a path interval, possibly duplicating nodes in case of cycles
        unordered_map<id_t, pair<id_t, bool>>
        extract_linearized_path_graph(const PathPositionHandleGraph* graph, MutableHandleGraph* into,
                                      path_handle_t path_handle, size_t first, size_t last) const;
        
        
        /// associate a path position and strand to a surjected alignment against this path
        void set_path_position(const PathPositionHandleGraph* graph, const Alignment& surjected,
                               path_handle_t best_path_handle,
                               string& path_name_out, int64_t& path_pos_out, bool& path_rev_out) const;
        
        /// another  algorithm that doesn't blow up if the alignment preserved deletions
        void set_path_position_inexact(const PathPositionHandleGraph* graph, const Alignment& surjected,
                                       path_handle_t best_path_handle,
                                       string& path_name_out, int64_t& path_pos_out, bool& path_rev_out) const;
        
        ///////////////////////
        // Support methods for the spliced surject algorithm
        ///////////////////////
        
        /// reverses an adjacency list
        vector<vector<size_t>> reverse_adjacencies(const vector<vector<size_t>>& adj) const;
        
        /// returns a vector assignming each node to a connectd component, requires both the forward and reverse adjacency
        /// lists. optionally also returns the total number of components
        vector<size_t> connected_components(const vector<vector<size_t>>& adj,
                                            const vector<vector<size_t>>& rev_adj,
                                            size_t* num_comps_out) const;
        
        /// returns the transitive reduction of a topologically sorted DAG's adjacency list
        vector<vector<size_t>> transitive_reduction(const vector<vector<size_t>>& adj) const;
        
        /// returns the nodes in an adjacency list that have an outgoing edge through which all of their component's
        /// source-to-sink paths flow
        vector<size_t> find_constriction_edges(const vector<vector<size_t>>& adj) const;
        
        /// make a sentinel meant to indicate an unmapped read
        static Alignment make_null_alignment(const Alignment& source);
        
        /// the graph we're surjecting onto
        const PathPositionHandleGraph* graph = nullptr;
        
        /// the minimum length deletion that the spliced algorithm will interpret as a splice event
        size_t min_splice_length = 20;
    };
}

#endif
