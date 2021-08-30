/**
 * \file splicing.hpp
 *
 * Defines splicing related objects and functions
 *
 */
#ifndef VG_SPLICE_REGION_HPP_INCLUDED
#define VG_SPLICE_REGION_HPP_INCLUDED

#include "dinucleotide_machine.hpp"
#include "incremental_subgraph.hpp"
#include "aligner.hpp"
#include "multipath_alignment.hpp"
#include "statistics.hpp"

namespace vg {

using namespace std;

/*
 * Object that represents:
 * 1. A table of the acceptable splice motifs and their scores
 * 2. A distribution of intron lengths and their scores
 */
class SpliceStats {
public:
    // Defaults to the three canonical human splicing motifs with the frequencies
    // reported in Burset, Seledstov, and Solovyev (2000)
    // - GT-AG: 0.9924
    // - GC-AG: 0.0069
    // - AT-AC: 0.0005
    // and a default intron length distribution trained on GENCODE v.29.
    SpliceStats(const GSSWAligner& scorer);
    
    // Construct with triples of (5' dinucleotide, 3' dinucleotide, frequency).
    // Frequencies may sum to < 1.0, but may not sum to > 1.0.
    // Lognormal parameters should be given in pairs of (mu, sigma).
    SpliceStats(const vector<tuple<string, string, double>>& motifs,
                const vector<double>& lognormal_mixture_weights,
                const vector<pair<double, double>>& lognormal_component_params,
                const GSSWAligner& scorer);
    
    // the number of splicing motifs
    size_t motif_size() const;
    // the dinucleotide motif on one side in the order that the nucleotides
    // are encountered when traversing into the intron
    const string& oriented_motif(size_t motif_num, bool left_side) const;
    // true if the corresponding motif would be encountered on the reverse
    // strand of a transcript
    bool motif_is_reverse(size_t motif_num) const;
    // the dinucleotide motif on one side ins its standard representation
    string unoriented_motif(size_t motif_num, bool left_side) const;
    // the score associated with a splicing motif
    int32_t motif_score(size_t motif_num) const;
    // the score associated with an intron length
    int32_t intron_length_score(int64_t length) const;
    // must be called if scoring parameters are changed
    void update_scoring(const GSSWAligner& scorer);
private:
    
    // internal function for the constructor
    void init(const vector<tuple<string, string, double>>& motifs,
              const vector<double>& lognormal_mixture_weights,
              const vector<pair<double, double>>& lognormal_component_params,
              const GSSWAligner& scorer);
    
    double intron_length_log_likelihood(int64_t length) const;
    
    // the splice table stored in its most useful form
    vector<tuple<string, string, int32_t>> motif_data;
    // the original input data, in case we need to update the scoring
    vector<tuple<string, string, double>> unaltered_motif_data;
    // intron length distribution component weights
    vector<double> mixture_weights;
    // intron length distribution component parameters
    vector<pair<double, double>> component_params;
    // the maximum possible log likelihood for an intron length
    double mode_log_likelihood;
    // a copy of the aligner's log base (it's convenient to have it here)
    double log_base;
};


/*
 * Object that identifies possible splice sites in a small region of
 * the graph and answers queries about them.
 */
class SpliceRegion {
public:
    
    SpliceRegion(const pos_t& seed_pos, bool search_left, int64_t search_dist,
                 const HandleGraph& graph,
                 const DinucleotideMachine& dinuc_machine,
                 const SpliceStats& splice_stats);
    SpliceRegion() = default;
    ~SpliceRegion() = default;
    
    // returns the locations in the of a splice motif and the distance of those location
    // from the search position. crashes if given a motif that was not provided
    // to the constructor.
    const vector<tuple<handle_t, size_t, int64_t>>& candidate_splice_sites(size_t motif_num) const;
    
    // the underlying subgraph we've extracted
    const IncrementalSubgraph& get_subgraph() const;
    
    // the position that extraction began from
    const pair<handle_t, size_t>& get_seed_pos() const;
    
private:

    IncrementalSubgraph subgraph;
    
    pair<handle_t, size_t> seed;
        
    vector<vector<tuple<handle_t, size_t, int64_t>>> motif_matches;
    
};

/*
 * HandleGraph implementation that fuses two incremental subgraphs at
 * an identified splice junction
 */
class JoinedSpliceGraph : public HandleGraph {
public:
    
    JoinedSpliceGraph(const HandleGraph& parent_graph,
                      const IncrementalSubgraph& left_subgraph,
                      handle_t left_splice_node, size_t left_splice_offset,
                      const IncrementalSubgraph& right_subgraph,
                      handle_t right_splice_node, size_t right_splice_offset);
    
    JoinedSpliceGraph() = default;
    ~JoinedSpliceGraph() = default;
    
    //////////////////////////
    /// Specialized interface
    //////////////////////////
    
    /// Translates an aligned path to the underlying graph of the two incremental subgraphs
    /// Also returns the indexes of the left and right splice points. If either is not present,
    /// the corresponding index is set to numeric_limits<size_t>::max()
    pair<size_t, size_t> translate_node_ids(Path& path) const;
    
    /// Get the handle corresponding to the seed node of the left subgraph
    handle_t left_seed_node() const;
    
    /// Get the handle corresponding to the seed node of the right subgraph
    handle_t right_seed_node() const;
    
    /// Get the handle corresponding to the left side of the splice join
    handle_t left_splice_node() const;
    
    /// Get the handle corresponding to the right side of the splice join
    handle_t right_splice_node() const;
    
    /// Get the minimum length of a walk from the left seed to the right seed
    int64_t min_link_length() const;
    
    /// Get the maximum length of a walk from the left seed to the right seed
    int64_t max_link_length() const;
    
    //////////////////////////
    /// HandleGraph interface
    //////////////////////////
    
    /// Method to check if a node exists by ID
    bool has_node(id_t node_id) const;
    
    /// Look up the handle for the node with the given ID in the given orientation
    handle_t get_handle(const id_t& node_id, bool is_reverse = false) const;
    
    /// Get the ID from a handle
    id_t get_id(const handle_t& handle) const;
    
    /// Get the orientation of a handle
    bool get_is_reverse(const handle_t& handle) const;
    
    /// Invert the orientation of a handle (potentially without getting its ID)
    handle_t flip(const handle_t& handle) const;
    
    /// Get the length of a node
    size_t get_length(const handle_t& handle) const;
    
    /// Get the sequence of a node, presented in the handle's local forward
    /// orientation.
    string get_sequence(const handle_t& handle) const;
    
    /// Loop over all the handles to next/previous (right/left) nodes. Passes
    /// them to a callback which returns false to stop iterating and true to
    /// continue. Returns true if we finished and false if we stopped early.
    bool follow_edges_impl(const handle_t& handle, bool go_left, const function<bool(const handle_t&)>& iteratee) const;
    
    /// Loop over all the nodes in the graph in their local forward
    /// orientations, in their internal stored order. Stop if the iteratee
    /// returns false. Can be told to run in parallel, in which case stopping
    /// after a false return value is on a best-effort basis and iteration
    /// order is not defined.
    bool for_each_handle_impl(const function<bool(const handle_t&)>& iteratee, bool parallel = false) const;
    
    /// Return the number of nodes in the graph
    size_t get_node_count() const;
    
    /// Return the smallest ID in the graph, or some smaller number if the
    /// smallest ID is unavailable. Return value is unspecified if the graph is empty.
    id_t min_node_id() const;
    
    /// Return the largest ID in the graph, or some larger number if the
    /// largest ID is unavailable. Return value is unspecified if the graph is empty.
    id_t max_node_id() const;
    
    ///////////////////////////////////
    /// Optional HandleGraph interface
    ///////////////////////////////////
    
    /// Returns one base of a handle's sequence, in the orientation of the
    /// handle.
    char get_base(const handle_t& handle, size_t index) const;
    
    /// Returns a substring of a handle's sequence, in the orientation of the
    /// handle. If the indicated substring would extend beyond the end of the
    /// handle's sequence, the return value is truncated to the sequence's end.
    string get_subsequence(const handle_t& handle, size_t index, size_t size) const;
    
private:
    
    // the index subsequence of the underlying sequence that this handle corresponds
    // to (on the same strand as the handle)
    pair<size_t, size_t> underlying_interval(const handle_t& handle) const;
    
    const HandleGraph* parent_graph;
    
    const IncrementalSubgraph* left_subgraph;
    
    const IncrementalSubgraph* right_subgraph;
    
    vector<int64_t> left_handle_trans;
    
    vector<int64_t> right_handle_trans;
    
    size_t left_splice_offset;
    
    size_t right_splice_offset;
    
    // handles from the underlying subgraphs
    vector<size_t> handle_idxs;
    
    // the size of the prefix of the handles vector that come from the
    // left subgraph
    size_t num_left_handles;
};

multipath_alignment_t from_hit(const Alignment& alignment, const HandleGraph& graph,
                               const pos_t& hit_pos, const MaximalExactMatch& mem,
                               const GSSWAligner& scorer);

// return the position of the base when trimming a given length from either the start
// or the end of the alignment. softclips do not contributed to the total, and extra
// will be trimmed to avoid splitting an insertion. also returns the total amount of
// sequence trimmed (including softclip) and the score of the sub-alignment that would
// be removed.
tuple<pos_t, int64_t, int32_t> trimmed_end(const Alignment& aln, int64_t len, bool from_end,
                                           const HandleGraph& graph, const GSSWAligner& aligner);

// remove the portion of the path either to the left or right of a given position along it,
// described in terms of the indexes in the path object. returns true if any edits were
// actually removed.
bool trim_path(path_t* path, bool from_left, int64_t mapping_idx, int64_t edit_idx, int64_t base_idx);

// consumes two multipath alignments and returns a spliced multipath alignment, which is connected
// by the splice segment. the splice segment may be discontiguous at the splice point, which is
// is indicated by giving the index of the first mapping that is past the splice junctinon.
// the splice score refers to the score value of th splice itself
multipath_alignment_t&& fuse_spliced_alignments(const Alignment& alignment,
                                                multipath_alignment_t&& left_mp_aln, multipath_alignment_t&& right_mp_aln,
                                                int64_t left_bridge_point, const Alignment& splice_segment,
                                                int64_t splice_junction_idx, int32_t splice_score, const GSSWAligner& scorer,
                                                const HandleGraph& graph);  

}

#endif // VG_SPLICE_REGION_HPP_INCLUDED
