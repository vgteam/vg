#ifndef VG_SURJECTOR_HPP_INCLUDED
#define VG_SURJECTOR_HPP_INCLUDED

/** \file
 *
 *  A class to hold surjection algorithms that do lossy realignment restricted to paths in the graph
 */

#include <set>
#include <atomic>
#include <sstream>
#include <algorithm>
#include <functional>
#include <limits>

#include "aligner.hpp"
#include "handle.hpp"
#include "path.hpp"
#include <vg/vg.pb.h>
#include "multipath_alignment.hpp"
#include "algorithms/pad_band.hpp"


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
                          const unordered_set<path_handle_t>& paths,
                          string& path_name_out,
                          int64_t& path_pos_out,
                          bool& path_rev_out,
                          bool allow_negative_scores = false,
                          bool preserve_deletions = false) const;
        
        /// Same as above, but include alignments to all paths instead of only the optimal one
        vector<Alignment> multi_surject(const Alignment& source,
                                        const unordered_set<path_handle_t>& paths,
                                        vector<tuple<string, int64_t, bool>>& positions_out,
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
                          const unordered_set<path_handle_t>& paths,
                          bool allow_negative_scores = false,
                          bool preserve_deletions = false) const;
        
        /// Same as above, but include alignments to all paths instead of only the optimal one
        vector<Alignment> multi_surject(const Alignment& source,
                                        const unordered_set<path_handle_t>& paths,
                                        bool allow_negative_scores = false,
                                        bool preserve_deletions = false) const;
        
        /// Same semantics as with alignments except that connections are always
        /// preserved as splices. The output consists of a multipath alignment with
        /// a single path, separated by splices (either from large deletions or from
        /// connections)
        multipath_alignment_t surject(const multipath_alignment_t& source,
                                      const unordered_set<path_handle_t>& paths,
                                      string& path_name_out, int64_t& path_pos_out,
                                      bool& path_rev_out,
                                      bool allow_negative_scores = false,
                                      bool preserve_deletions = false) const;
        
        /// Same as above, but include alignments to all paths instead of only the optimal one
        vector<multipath_alignment_t> multi_surject(const multipath_alignment_t& source,
                                                    const unordered_set<path_handle_t>& paths,
                                                    vector<tuple<string, int64_t, bool>>& positions_out,
                                                    bool allow_negative_scores = false,
                                                    bool preserve_deletions = false) const;
        
        /// a local type that represents a read interval matched to a portion of the alignment path
        using path_chunk_t = pair<pair<string::const_iterator, string::const_iterator>, path_t>;
        
        /// the minimum length deletion that the spliced algorithm will interpret as a splice event
        int64_t min_splice_length = 20;
        
        int64_t dominated_path_chunk_diff = 10;

        /// the minimum length apparent intron that we will try to repair
        int64_t min_splice_repair_length = 250;

        /// the maximum length of a tail that we will try to align
        size_t max_tail_length = 10000;
        
        /// We have a different default max_subgraph_bases_per_read_base to use for spliced alignment.
        static constexpr double SPLICED_DEFAULT_SUBGRAPH_LIMIT = 16 * 1024 * 1024 / 125.0;
        /// And an accessible default max_subgraph_bases_per_read_base for normal alignment.
        static constexpr double DEFAULT_SUBGRAPH_LIMIT = 100 * 1024 / 125.0;
        /// How big of a graph (in graph bases per read base) should we ever try to align against for realigning surjection?
        double max_subgraph_bases_per_read_base = DEFAULT_SUBGRAPH_LIMIT;
       
        
        /// in spliced surject, downsample if the base-wise average coverage by chunks is this high
        int64_t min_fold_coverage_for_downsample = 8;
        /// while downsampling, try to get down to this coverage on each base
        int64_t downsample_coverage = 16;
        
        int64_t min_shift_for_prune = 32 * 1024;
        int64_t shift_prune_diff = 16 * 1024;
        
        /// And have we complained about hitting it?
        mutable atomic_flag warned_about_subgraph_size = ATOMIC_FLAG_INIT;
        
        bool prune_suspicious_anchors = false;
        int64_t max_tail_anchor_prune = 4;
        static constexpr int64_t DEFAULT_MAX_SLIDE = 6;
        /// Declare an anchor suspicious if it appears again at any offset up
        /// to this limit or the anchor length.
        int64_t max_slide = DEFAULT_MAX_SLIDE;
        double low_complexity_p_value = .0075;
        int64_t max_low_complexity_anchor_prune = 40;
        int64_t max_low_complexity_anchor_trim = 65;
        /// When examining anchors for low complexity, try and make them at
        /// least this long. To ensure orientation symmetry, we will make
        /// anchors with the oppsite parity (even if this is odd, or odd if
        /// this is even) 1bp longer.
        int64_t pad_suspicious_anchors_to_length = 12;
        
        // A function for computing band padding
        std::function<size_t(const Alignment&, const HandleGraph&)> choose_band_padding;
        
        /// How many anchors (per path) will we use when surjecting using
        /// anchors?
        /// Excessive anchors will be pruned away.
        size_t max_anchors = std::numeric_limits<size_t>::max();
        
        bool annotate_with_all_path_scores = false;
        
    protected:
        
        void surject_internal(const Alignment* source_aln, const multipath_alignment_t* source_mp_aln,
                              vector<Alignment>* alns_out, vector<multipath_alignment_t>* mp_alns_out,
                              const unordered_set<path_handle_t>& paths,
                              vector<tuple<string, int64_t, bool>>& positions_out, bool all_paths,
                              bool allow_negative_scores, bool preserve_deletions) const;
        
        Alignment
        realigning_surject(const PathPositionHandleGraph* graph, const Alignment& source,
                           const path_handle_t& path_handle, bool rev_strand,
                           const vector<path_chunk_t>& path_chunks,
                           const vector<pair<step_handle_t, step_handle_t>>& ref_chunks,
                           pair<step_handle_t, step_handle_t>& path_range_out,
                           bool allow_negative_scores,
                           bool preserve_N_alignments = false,
                           bool sinks_are_anchors = false,
                           bool sources_are_anchors = false,
                           vector<pair<step_handle_t, step_handle_t>>* all_path_ranges_out = nullptr) const;
        
        multipath_alignment_t
        spliced_surject(const PathPositionHandleGraph* path_position_graph,
                        const string& src_sequence, const string& src_quality,
                        const int32_t src_mapping_quality,
                        const path_handle_t& path_handle, bool rev_strand,
                        vector<path_chunk_t>& path_chunks,
                        vector<pair<step_handle_t, step_handle_t>>& ref_chunks,
                        vector<tuple<size_t, size_t, int32_t>>& connections,
                        pair<step_handle_t, step_handle_t>& path_range_out,
                        bool allow_negative_scores, bool deletions_as_splices) const;
        
        ///////////////////////
        // Support methods for the realigning surject algorithm
        ///////////////////////
        
        /// get the chunks of the alignment path that follow the given reference paths
        unordered_map<pair<path_handle_t, bool>, pair<vector<path_chunk_t>, vector<pair<step_handle_t, step_handle_t>>>>
        extract_overlapping_paths(const PathPositionHandleGraph* graph, const Alignment& source,
                                  const unordered_set<path_handle_t>& surjection_paths) const;
        
        /// same semantics except for a multipath alignment
        unordered_map<pair<path_handle_t, bool>, pair<vector<path_chunk_t>, vector<pair<step_handle_t, step_handle_t>>>>
        extract_overlapping_paths(const PathPositionHandleGraph* graph,
                                  const multipath_alignment_t& source,
                                  const unordered_set<path_handle_t>& surjection_paths,
                                  unordered_map<pair<path_handle_t, bool>, vector<tuple<size_t, size_t, int32_t>>>& connections_out) const;
        
        /// remove any path chunks and corresponding ref chunks that are identical to a longer
        /// path chunk over the region where they overlap
        void filter_redundant_path_chunks(bool path_rev, vector<path_chunk_t>& path_chunks,
                                          vector<pair<step_handle_t, step_handle_t>>& ref_chunks,
                                          vector<tuple<size_t, size_t, int32_t>>& connections) const;
        
        void prune_and_trim_anchors(const string& sequence, vector<path_chunk_t>& path_chunks,
                                    vector<pair<step_handle_t, step_handle_t>>& step_ranges) const;
        
        /// Compute the widest end-inclusive interval of path positions that
        /// the realigned sequence could align to, or an interval where start >
        /// end if there are no path chunks.
        pair<size_t, size_t>
        compute_path_interval(const PathPositionHandleGraph* graph, const Alignment& source, path_handle_t path_handle,
                              bool rev_strand, const vector<path_chunk_t>& path_chunks,
                              const vector<pair<step_handle_t, step_handle_t>>& ref_chunks,
                              bool no_left_expansion, bool no_right_expansion) const;
        
        /// use the graph position bounds and the path range bounds to assign a path position to a surjected read
        void set_path_position(const PathPositionHandleGraph* graph, const pos_t& init_surj_pos,
                               const pos_t& final_surj_pos,
                               const step_handle_t& range_begin, const step_handle_t& range_end,
                               bool rev_strand, string& path_name_out, int64_t& path_pos_out, bool& path_rev_out) const;
        
        template<class AlnType>
        string path_score_annotations(const unordered_map<pair<path_handle_t, bool>, pair<AlnType, pair<step_handle_t, step_handle_t>>>& surjections) const;
        
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
        
        /// eliminate any path chunks that have the exact same colinearities as another but are much shorter
        vector<vector<size_t>> remove_dominated_chunks(const string& src_sequence,
                                                       const vector<vector<size_t>>& adj,
                                                       vector<path_chunk_t>& path_chunks,
                                                       vector<pair<step_handle_t, step_handle_t>>& ref_chunks,
                                                       vector<tuple<size_t, size_t, int32_t>>& connections) const;
        
        /// if any anchors overlap each other, cut the second at the implied overlap position
        void cut_anchors(bool rev_strand, vector<path_chunk_t>& path_chunks,
                         vector<pair<step_handle_t, step_handle_t>>& ref_chunks,
                         vector<tuple<size_t, size_t, int32_t>>& connections) const;
        
        /// if there are too many chunks, downsample to a given level
        void downsample_chunks(const string& src_sequence,
                               vector<path_chunk_t>& path_chunks,
                               vector<pair<step_handle_t, step_handle_t>>& ref_chunks,
                               vector<tuple<size_t, size_t, int32_t>>& connections) const;
        
        /// returns all sets of chunks such that 1) all of chunks on the left set abut all of the chunks on the right
        /// set on the read, 2) all source-to-sink paths in the connected component go through an edge between
        /// the left and right sides, 3) all of the chunks that do not have a connection between them are fully
        /// connected (i.e. form a biclique)
        vector<pair<vector<size_t>, vector<size_t>>> find_constriction_bicliques(const vector<vector<size_t>>& adj,
                                                                                 const string& src_sequence,
                                                                                 const string& src_quality,
                                                                                 vector<path_chunk_t>& path_chunks,
                                                                                 vector<pair<step_handle_t, step_handle_t>>& ref_chunks,
                                                                                 const vector<tuple<size_t, size_t, int32_t>>& connections) const;
        
        void prune_unconnectable(vector<vector<size_t>>& adj,
                                 vector<vector<tuple<size_t, int32_t, bool>>>& splice_adj,
                                 vector<size_t>& component,
                                 vector<vector<size_t>>& comp_groups,
                                 vector<path_chunk_t>& path_chunks,
                                 vector<pair<step_handle_t, step_handle_t>>& ref_chunks) const;
        
        /// make a sentinel meant to indicate an unmapped read
        static Alignment make_null_alignment(const Alignment& source);
        
        static multipath_alignment_t make_null_mp_alignment(const string& src_sequence,
                                                            const string& src_quality);
        
        template<class AlnType>
        static int32_t get_score(const AlnType& aln);
        
        /// the graph we're surjecting onto
        const PathPositionHandleGraph* graph = nullptr;
    };


    template<class AlnType>
    string Surjector::path_score_annotations(const unordered_map<pair<path_handle_t, bool>, pair<AlnType, pair<step_handle_t, step_handle_t>>>& surjections) const {
        
        vector<tuple<int32_t, string, bool>> paths;
        for (const auto& surjection : surjections) {
            paths.emplace_back(get_score(surjection.second.first), graph->get_path_name(surjection.first.first), surjection.first.second);
        }
        sort(paths.begin(), paths.end(), greater<tuple<int32_t, string, bool>>());
        
        stringstream sstrm;
        
        for (size_t i = 0; i < paths.size(); ++i) {
            if (i != 0) {
                sstrm << ',';
            }
            sstrm << get<1>(paths[i]);
            sstrm << (get<2>(paths[i]) ? '-' : '+');
            sstrm << get<0>(paths[i]);
        }
        
        return sstrm.str();
    }
}

#endif
