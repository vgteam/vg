#ifndef VG_SURJECTOR_HPP_INCLUDED
#define VG_SURJECTOR_HPP_INCLUDED

/** \file
 *
 *  A class to hold surjection algorithms that do lossy realignment restricted to paths in the graph
 */

#include <set>
#include <map>
#include <atomic>
#include <sstream>
#include <algorithm>
#include <functional>
#include <limits>
#include <queue>

#include "aligner.hpp"
#include "handle.hpp"
#include "path.hpp"
#include <vg/vg.pb.h>
#include "multipath_alignment.hpp"
#include "algorithms/pad_band.hpp"


namespace vg {

using namespace std;
    
    /**
     * Widget to surject alignments down to linear paths in the graph.
     *
     * Assumes the alignments actually go with the graph; the caller is
     * repsonsible for ensuring that e.g. all nodes referenced by the
     * alignments actually exist.
     */
    class Surjector : public AlignerClient {
    public:
        
        Surjector(const PathPositionHandleGraph* graph);
        
        /// Override alignment score setting to let DP use slightly adjusted scores.
        /// The provided scores will be adjusted to increase gap open cost
        /// slightly when doing DP, but output scores will be scored with the
        /// provided parameters.
        virtual void set_alignment_scores(const int8_t* score_matrix, int8_t gap_open, int8_t gap_extend, int8_t full_length_bonus);
        
        /// Extract the portions of an alignment that are on a chosen set of paths and try to
        /// align realign the portions thaet are off of the chosen paths to the intervening
        /// path segments to obtain an alignment that is fully restricted to the paths.
        ///
        /// Also returns the path name, position, and strand of the new alignment.
        ///
        /// Optionally either allow softclips so that the alignment has a nonnegative score on
        /// the path or require the full-length alignment, possibly creating a negative score.
        ///
        /// Also optionally leaves deletions against the reference path in the final alignment
        /// (useful for splicing).
        vector<Alignment> surject(const Alignment& source,
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
        vector<Alignment> surject(const Alignment& source,
                                  const unordered_set<path_handle_t>& paths,
                                  bool allow_negative_scores = false,
                                  bool preserve_deletions = false) const;
        
        /// Same semantics as with alignments except that connections are always
        /// preserved as splices. The output consists of a multipath alignment with
        /// a single path, separated by splices (either from large deletions or from
        /// connections)
        vector<multipath_alignment_t> surject(const multipath_alignment_t& source,
                                              const unordered_set<path_handle_t>& paths,
                                              vector<tuple<string, int64_t, bool>>& positions_out,
                                              bool allow_negative_scores = false,
                                              bool preserve_deletions = false) const;
        
        /// a local type that represents a read interval matched to a portion of the alignment path
        using path_chunk_t = pair<pair<string::const_iterator, string::const_iterator>, path_t>;

        /// When doing DP alignments, we use slightly adjusted alignment scores, which need to live in their own Aligner
        std::unique_ptr<Aligner> dp_aligner;
        
        /// When doing DP alignments, we multiply scores by this factor before slightly increasing gap open. (Changes won't take effect until set_alignment_scores() is called.)
        int8_t dp_score_scale = 10;
        /// When doing DP alignment, we increase gap open score by this much. (Changes won't take effect until set_alignment_scores() is called.)
        int8_t dp_gap_open_extra_cost = 0;

        

        /// the minimum length deletion that the spliced algorithm will interpret as a splice event
        int64_t min_splice_length = 20;
        
        int64_t dominated_path_chunk_diff = 10;

        /// the minimum length apparent intron that we will try to repair
        int64_t min_splice_repair_length = 250;

        /// the maximum length of a tail that we will try to align
        size_t max_tail_length = 10000;

        // the maximum number of estimated band cells that we are willing to try to fill when connecting anchors
        uint64_t max_band_cells = 8000000000;
        
        /// We have a different default max_subgraph_bases_per_read_base to use for spliced alignment.
        static constexpr double SPLICED_DEFAULT_SUBGRAPH_LIMIT = 16 * 1024 * 1024 / 125.0;
        /// And an accessible default max_subgraph_bases_per_read_base for normal alignment.
        static constexpr double DEFAULT_SUBGRAPH_LIMIT = 100 * 1024 / 125.0;
        /// How big of a graph (in graph bases per read base) should we ever try to align against for realigning surjection?
        double max_subgraph_bases_per_read_base = DEFAULT_SUBGRAPH_LIMIT;
        /// Don't refuse to align  (graph size) * (read size) is at least this size (overrides max_subgraph_bases_per_read_base)
        int64_t min_absolute_align_size_to_refuse = 1024;
        
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

        /// Return alignments to overlapping paths as secondaries instead of just the best one
        bool multimap_to_all_paths = false;
        
        // A function for computing band padding
        std::function<size_t(const Alignment&, const HandleGraph&)> choose_band_padding;
        
        /// How many anchors (per path) will we use when surjecting using
        /// anchors?
        /// Excessive anchors will be pruned away.
        size_t max_anchors = std::numeric_limits<size_t>::max();
        
        /// Should we report supplementary alignments?
        bool report_supplementary = true;
        /// What fraction of a read length separation should we expect before separating into supplementaries?
        double read_length_prop_disjoint_gap = 1.0;
        /// The maximum gap size before we always separate into supplementary alignments
        size_t max_disjoint_gap = 10000;
        /// The minimum gap size we should never separate into supplementary alignments
        size_t min_disjoint_gap = 50;
        /// The max size overlap on the read between supplementaries before we designate them as non-colinear
        int64_t disjoint_interval_allowable_overlap = 4;
        
        bool annotate_with_all_path_scores = false;
        bool annotate_with_graph_alignment = false;
        
    protected:

        /// Do the extra score setup for the DP-only Aligner.
        void set_dp_alignment_scores(const int8_t* score_matrix, int8_t gap_open, int8_t gap_extend, int8_t full_length_bonus);
        
        void surject_internal(const Alignment* source_aln, const multipath_alignment_t* source_mp_aln,
                              vector<Alignment>* alns_out, vector<multipath_alignment_t>* mp_alns_out,
                              const unordered_set<path_handle_t>& paths,
                              vector<tuple<string, int64_t, bool>>& positions_out,
                              bool allow_negative_scores, bool preserve_deletions) const;
        
        vector<pair<Alignment, pair<step_handle_t, step_handle_t>>>
        realigning_surject(const PathPositionHandleGraph* graph, const Alignment& source,
                           const path_handle_t& path_handle, bool rev_strand,
                           const vector<path_chunk_t>& path_chunks,
                           const vector<pair<step_handle_t, step_handle_t>>& ref_chunks,
                           bool allow_negative_scores,
                           bool preserve_N_alignments = false,
                           bool sinks_are_anchors = false,
                           bool sources_are_anchors = false,
                           vector<vector<pair<step_handle_t, step_handle_t>>>* all_path_ranges_out = nullptr,
                           size_t override_read_length = 0) const;
        
        vector<pair<multipath_alignment_t, pair<step_handle_t, step_handle_t>>>
        spliced_surject(const PathPositionHandleGraph* path_position_graph,
                        const string& src_sequence, const string& src_quality,
                        const int32_t src_mapping_quality,
                        const path_handle_t& path_handle, bool rev_strand,
                        vector<path_chunk_t>& path_chunks,
                        vector<pair<step_handle_t, step_handle_t>>& ref_chunks,
                        vector<tuple<size_t, size_t, int32_t>>& connections,
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
        
        /// Compute a set of end-inclusive intervals of path positions that the
        /// realign sequence could align to, and associate with each one a vector
        /// indexes into path_chunks to indicate which path chunks it contains.
        vector<tuple<size_t, size_t, vector<size_t>>>
        compute_disjoint_path_intervals(const PathPositionHandleGraph* graph, const Alignment& source, path_handle_t path_handle,
                                        bool rev_strand, const vector<path_chunk_t>& path_chunks,
                                        const vector<pair<step_handle_t, step_handle_t>>& ref_chunks,
                                        bool no_left_expansion, bool no_right_expansion, size_t max_gap) const;
        
        /// use the graph position bounds and the path range bounds to assign a path position to a surjected read
        void set_path_position(const PathPositionHandleGraph* graph, const pos_t& init_surj_pos,
                               const pos_t& final_surj_pos,
                               const step_handle_t& range_begin, const step_handle_t& range_end,
                               bool rev_strand, string& path_name_out, int64_t& path_pos_out, bool& path_rev_out) const;
        
        template<class AlnType>
        string path_score_annotations(const unordered_map<pair<path_handle_t, bool>, vector<pair<AlnType, pair<step_handle_t, step_handle_t>>>>& surjections) const;
        
        // helpers to choose one among supplementary alignments to be the primary
        template<class AlnType>
        void choose_primary_internal(vector<pair<AlnType, pair<step_handle_t, step_handle_t>>>& surjections,
                                     const function<void(AlnType&)>& annotate_supplementary) const;
        void choose_primary(vector<pair<Alignment, pair<step_handle_t, step_handle_t>>>& surjections) const {
            function<void(Alignment&)> annotate_supplementary = [](Alignment& aln) { set_annotation<bool>(aln, "supplementary", true); };
            choose_primary_internal(surjections, annotate_supplementary);
        }
        void choose_primary(vector<pair<multipath_alignment_t, pair<step_handle_t, step_handle_t>>>& surjections) const {
            function<void(multipath_alignment_t&)> annotate_supplementary = [](multipath_alignment_t& mp_aln) { mp_aln.set_annotation("supplementary", true); };
            choose_primary_internal(surjections, annotate_supplementary);
        }

        // helpers to add the SA tag
        template<class AlnType>
        void add_SA_tag_internal(vector<AlnType>& surjected, const vector<tuple<string, int64_t, bool>>& positions,
                                 const PathPositionHandleGraph& graph, bool spliced, const function<void(AlnType&,const string&)>& update_sa) const;
        static string update_tag_string_for_SA(const string& tags, const string& sa_val);
        void add_SA_tag(vector<Alignment>& surjected, const vector<tuple<string, int64_t, bool>>& positions,
                        const PathPositionHandleGraph& graph, bool spliced) const {
            function<void(Alignment&,const string&)> update_sa = [](Alignment& aln,const string& val) {
                string annotation;
                if (has_annotation(aln, "tags")) {
                    annotation = std::move(get_annotation<string>(aln, "tags"));
                }
                annotation = update_tag_string_for_SA(annotation, val);
                set_annotation(aln, "tags", annotation);
            };
            add_SA_tag_internal(surjected, positions, graph, spliced, update_sa);
        }
        void add_SA_tag(vector<multipath_alignment_t>& surjected, const vector<tuple<string, int64_t, bool>>& positions,
                        const PathPositionHandleGraph& graph, bool spliced) const {
            function<void(multipath_alignment_t&,const string&)> update_sa = [](multipath_alignment_t& mp_aln,const string& val) {
                string annotation;
                if (mp_aln.has_annotation("tags")) {
                    annotation = *((const string*) mp_aln.get_annotation("tags").second);
                }
                annotation = update_tag_string_for_SA(annotation, val);
                mp_aln.set_annotation("tags", annotation);
            };
            add_SA_tag_internal(surjected, positions, graph, spliced, update_sa);
        }
        
        // calculate the total length of overlaps between alignments, measured in terms of read sequence
        template<class AlnType>
        int64_t total_overlap(const vector<pair<AlnType, pair<step_handle_t, step_handle_t>>>& surjections) const;
        
        // choose a primary strand, and if it is possible to cover the read with disjoint alignments, follow it with supplementary strands
        template<class AlnType>
        vector<pair<path_handle_t, bool>> 
        supplementary_cover(const unordered_map<pair<path_handle_t, bool>, vector<pair<AlnType, pair<step_handle_t, step_handle_t>>>>& surjections) const;
        
        // identify the strand that has the best alignment to choose it as the primary, optionally restrict to a subset of strands
        template<class AlnType>
        pair<path_handle_t, bool> choose_primary_strand(const unordered_map<pair<path_handle_t, bool>, vector<pair<AlnType, pair<step_handle_t, step_handle_t>>>>& surjections,
                                                        const unordered_set<pair<path_handle_t, bool>>* among_strands = nullptr) const;

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
        
        void annotate_graph_cigar(vector<Alignment>& surjections, const Alignment& source, 
                                  const vector<tuple<string, int64_t, bool>>& positions) const;
        
        void annotate_graph_cigar(vector<multipath_alignment_t>& surjections, const multipath_alignment_t& source, 
                                  const vector<tuple<string, int64_t, bool>>& positions) const;
        
        template<class AlnType>
        static int32_t get_score(const AlnType& aln);

        template<class AlnType>
        static size_t count_mismatches(const AlnType& surjected);

        template<class AlnType>
        vector<pair<int, char>> get_cigar(const AlnType& surjected, const tuple<string, int64_t, bool>& position, const PathPositionHandleGraph& graph, bool spliced) const;
        
        /// the graph we're surjecting onto
        const PathPositionHandleGraph* graph = nullptr;
    };

    template<class AlnType>
    string Surjector::path_score_annotations(const unordered_map<pair<path_handle_t, bool>, vector<pair<AlnType, pair<step_handle_t, step_handle_t>>>>& surjections) const {
        
        vector<tuple<int32_t, string, bool>> paths;
        for (const auto& path_surjections : surjections) {
            for (const auto& surjection : path_surjections.second) {
                if (!is_supplementary(surjection.first)) {
                    paths.emplace_back(get_score(surjection.first),
                                       graph->get_path_name(path_surjections.first.first), path_surjections.first.second);
                }
            }
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

    template<class AlnType>
    void Surjector::choose_primary_internal(vector<pair<AlnType, pair<step_handle_t, step_handle_t>>>& surjections,
                                            const function<void(AlnType&)>& annotate_supplementary) const {
        if (surjections.size() > 1) {
            size_t opt_idx = 0;
            int32_t opt_score = get_score(surjections.front().first);
            for (size_t i = 1; i < surjections.size(); ++i) {
                int32_t score = get_score(surjections[i].first);
                if (score > opt_score) {
                    opt_score = score;
                    opt_idx = i;
                }
            }
            
            for (size_t i = 0; i < surjections.size(); ++i) {
                if (i != opt_idx) {
                    annotate_supplementary(surjections[i].first);
                }
            }
        }
    }

    template<class AlnType>
    int64_t Surjector::total_overlap(const vector<pair<AlnType, pair<step_handle_t, step_handle_t>>>& surjections) const {
        
        vector<pair<int64_t, int64_t>> intervals(surjections.size());
        for (size_t i = 0; i < surjections.size(); ++i) {
            intervals[i] = aligned_interval(surjections[i].first);
        }
        sort(intervals.begin(), intervals.end());
        
        priority_queue<int64_t, vector<int64_t>, greater<int64_t>> ends;
        int64_t overlap = 0;
        int64_t previous = -1;
        for (size_t i = 0; i < intervals.size(); ++i) {
            
            const auto& interval = intervals[i];
            // add overlaps from constant-depth sub-intervals that end before this interval starts
            while (!ends.empty() && ends.top() <= interval.first) {
                overlap += (ends.size() - 1) * (ends.top() - previous);
                previous = ends.top();
                ends.pop();
            }
            
            if (!ends.empty()) {
                // add the constant-depth sub-interval that ends at this interval's beginning
                overlap += (ends.size() - 1) * (interval.first - previous);
            }
            previous = interval.first;
            // queue up the point at which this interval is subtracted from the depth
            ends.push(interval.second);
        }
        // clear out the queue
        while (!ends.empty()) {
            overlap += (ends.size() - 1) * (ends.top() - previous);
            previous = ends.top();
            ends.pop();
        }
        
        return overlap;
    }

    template<class AlnType>
    void Surjector::add_SA_tag_internal(vector<AlnType>& surjected, const vector<tuple<string, int64_t, bool>>& positions, 
                                        const PathPositionHandleGraph& graph, bool spliced, const function<void(AlnType&,const string&)>& update_sa) const {
        
        if (surjected.size() <= 1) {
            // there are no other reads to get supplementary information from
            return;
        }
        size_t primary_idx = -1;
        bool has_supplementary = false;
        for (size_t i = 0; i < surjected.size() && (!has_supplementary || primary_idx == -1); ++i) {
            if (!is_supplementary(surjected[i])) {
                if (primary_idx == -1) {
                    primary_idx = i;
                }
            }
            else {
                has_supplementary = true;
            }
        }
        if (primary_idx == -1 || !has_supplementary) {
            // there is no primary or no supplementaries
            return;
        }

        // gather the information from each alignment that will go into the other alignment's SA tag
        vector<string> SA_entry(surjected.size());
        for (size_t i = 0; i < surjected.size(); ++i) {
            if (i == primary_idx || is_supplementary(surjected[i])) {
                
                const auto& surj = surjected[i];

                stringstream strm;

                const auto& pos = positions[i];
                // note: convert to 1-based
                strm << (get<0>(pos).empty() ? string("*") : get<0>(pos)) << ',' << (get<1>(pos) >= 0 ? get<1>(pos) + 1 : (int64_t) 0) << ',' << (get<2>(pos) ? '-' : '+') << ',';
                auto cigar = get_cigar(surj, pos, graph, spliced);
                if (cigar.empty()) {
                    strm << '*';
                }
                else {
                    for (const auto& op : cigar)  {
                        strm << op.first << op.second;
                    }
                }
                strm << ',' << surj.mapping_quality() << ',' << count_mismatches(surj) << ';';

                SA_entry[i] = std::move(strm.str());
            }
        }

        // add the SA tag onto the primary and the supplementaries
        for (size_t i = 0; i < surjected.size(); ++i) {
            if (i != primary_idx && !is_supplementary(surjected[i])) {
                continue;
            }
            // this is a primary or a supplementary that needs the SA tag
            stringstream strm;
            for (size_t j = 0; j < surjected.size(); ++j) {
                if (j == i || (j != primary_idx && !is_supplementary(surjected[j]))) {
                    continue;
                }
                strm << SA_entry[j];
            }

            update_sa(surjected[i], strm.str());
        }
    }
    
    template<class AlnType>
    pair<path_handle_t, bool> 
    Surjector::choose_primary_strand(const unordered_map<pair<path_handle_t, bool>, vector<pair<AlnType, pair<step_handle_t, step_handle_t>>>>& surjections,
                                     const unordered_set<pair<path_handle_t, bool>>* among_strands) const {

        pair<path_handle_t, bool> best_path_strand;
        int32_t score = numeric_limits<int32_t>::min();
        for (const auto& strand_surjections : surjections) {
            if (among_strands && !among_strands->count(strand_surjections.first)) {
                continue;
            }
            int32_t total_score = 0;
            for (const auto& surjection : strand_surjections.second) {
                total_score += get_score(surjection.first);
            }
            // approximate the score from overlapping alignments
            // TODO: it would be possible to make this exact
            total_score -= total_overlap(strand_surjections.second) * get_aligner()->match;
            
            if (total_score >= score) {
#ifdef debug_anchored_surject
                cerr << "surjection against path " << graph->get_path_name(surjections.first.first) << " strand " << surjections.first.second << " achieves highest score (so far) of " << total_score << endl;
#endif
                score = total_score;
                best_path_strand = strand_surjections.first;
            }
        }
        return best_path_strand;
    }

    template<class AlnType>
    vector<pair<path_handle_t, bool>> 
    Surjector::supplementary_cover(const unordered_map<pair<path_handle_t, bool>, vector<pair<AlnType, pair<step_handle_t, step_handle_t>>>>& strand_surjections) const {

        int64_t read_len = strand_surjections.begin()->second.front().first.sequence().size();

        // identify each strand with the interval of the read it aligns
        vector<tuple<int64_t, int64_t, path_handle_t, bool, size_t>> intervals;
        for (const auto& strand_surjection : strand_surjections) {            
            for (size_t i = 0; i < strand_surjection.second.size(); ++i) {
                const auto& surjection = strand_surjection.second[i];
                auto interval = aligned_interval(surjection.first);
                if (interval.first < interval.second) {
                    intervals.emplace_back(interval.first, interval.second, strand_surjection.first.first, strand_surjection.first.second, i);
                }
            }
        }

        sort(intervals.begin(), intervals.end());
        
        // try to find a set of path strands that completely partition the read

        // sparse DP structure, records of (end index -> (bases covered, interval idx, backpointer))
        map<int64_t, tuple<int64_t, size_t, int64_t>> dp;
        // boundary condition
        dp[0] = tuple<int64_t, size_t, int64_t>(0, -1, -1);
        for (size_t i = 0; i < intervals.size(); ++i) {

            const auto& interval = intervals[i];
            
            auto here = dp.find(get<1>(interval));

            // TODO: with a sparse RMQ, I could guarantee O(log n) time here
            int64_t lower_lim = get<0>(interval) - disjoint_interval_allowable_overlap;
            int64_t upper_lim = min<int64_t>(get<0>(interval) + disjoint_interval_allowable_overlap, get<1>(interval));
            for (auto it = dp.lower_bound(lower_lim); it != dp.end() && it->first < upper_lim; ++it) {
                int64_t cov = get<0>(it->second) + (get<1>(interval) - max<int64_t>(it->first, get<0>(interval)));

                if (here == dp.end()) {
                    here = dp.insert(make_pair(get<1>(interval), tuple<int64_t, size_t, int64_t>(cov, i, it->first))).first;
                }
                else if (cov > get<0>(here->second)) {
                    here->second = tuple<int64_t, size_t, int64_t>(cov, i, it->first);
                }
            }
        }

        // choose a feasible optimum
        auto opt = dp.end();
        for (auto it = dp.lower_bound(max<int64_t>(read_len - disjoint_interval_allowable_overlap, 0)); it != dp.end(); ++it) {
            if (opt == dp.end() || get<0>(it->second) > get<0>(opt->second)) {
                opt = it;
            }
        }

        // traceback
        // FIXME: unless we implement a system to select a subset of supplementaries from each strand, it's possible that
        // we select non-disjoint alignment sets here
        unordered_set<pair<path_handle_t, bool>> strands;
        if (opt != dp.end()) {
            auto here = opt;
            while (get<2>(here->second) >= 0) {
                const auto& interval = intervals[get<1>(here->second)];
                strands.emplace(get<2>(interval), get<3>(interval));
                here = dp.find(get<2>(here->second));
            }
        }

        // convert into the output
        vector<pair<path_handle_t, bool>> return_val;
        if (!strands.empty()) {
            return_val.insert(return_val.end(), strands.begin(), strands.end());
            sort(return_val.begin(), return_val.end());
            auto primary = choose_primary_strand(strand_surjections, &strands);
            auto it = find(return_val.begin(), return_val.end(), primary);
            std::swap(*it, return_val.front());
        }
        else {
            return_val.emplace_back(choose_primary_strand(strand_surjections));
        }
        return return_val;
    }

}

#endif
