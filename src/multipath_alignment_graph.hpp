//
//  multipath_alignment_graph.hpp
//
// Contains class that computes multipath alignments to a small graph given
// a set of anchoring paths
//

#ifndef multipath_alignment_graph_hpp
#define multipath_alignment_graph_hpp

#include <stdio.h>
#include <vector>
#include <unordered_map>
#include <algorithm>

#include "structures/rank_pairing_heap.hpp"

#include <vg/vg.pb.h>
#include "vg.hpp"
#include "snarls.hpp"
#include "multipath_mapper.hpp"

namespace vg {
    
    
    // TODO: put in MultipathAlignmentGraph namespace
    class PathNode {
    public:
        string::const_iterator begin;
        string::const_iterator end;
        path_t path;
        
        // pairs of (target index, path length)
        vector<pair<size_t, size_t>> edges;
    };
    
    class MultipathAlignmentGraph {
    public:
        
        /// Create the constant injection translation data, which maps a node
        /// in the original graph to every one of its occurrences in the
        /// dagified graph, by reversing the projection translation. This data
        /// is needed to construct the MultipathAlignmentGraph, and to perform
        /// some other operations on it, but is big enough that it is worth not
        /// making it a member.
        static unordered_multimap<id_t, pair<id_t, bool>> create_injection_trans(const unordered_map<id_t, pair<id_t, bool>>& projection_trans);
        
        /// Create the constant injection translation data from a function instead
        /// of a map
        static unordered_multimap<id_t, pair<id_t, bool>> create_injection_trans(const HandleGraph& graph,
                                                                                 const function<pair<id_t, bool>(id_t)>& project);
        
        /// Create an identity projection translation from a DAG that did not
        /// need to be modified during dagification.
        static unordered_map<id_t, pair<id_t, bool>> create_identity_projection_trans(const HandleGraph& graph);
        
        /// Create a lambda function that projects using a map that projects
        static function<pair<id_t, bool>(id_t)> create_projector(const unordered_map<id_t, pair<id_t, bool>>& projection_trans);
        
        /// Construct a graph of the reachability between MEMs in a DAG-ified
        /// graph. If a GCSA is specified, use it to collapse MEMs whose
        /// lengths bump up against the GCSA's order limit on MEM length.
        /// Produces a graph with reachability edges. Assumes that the cluster
        /// is sorted by primarily length and secondarily lexicographically by
        /// read interval.
        /// If a hit fails to be walked ouut in the graph, it is removed from hits.
        MultipathAlignmentGraph(const HandleGraph& graph, MultipathMapper::memcluster_t& hits,
                                const function<pair<id_t, bool>(id_t)>& project,
                                const unordered_multimap<id_t, pair<id_t, bool>>& injection_trans,
                                vector<size_t>& path_node_provenance,
                                size_t max_branch_trim_length = 0, gcsa::GCSA* gcsa = nullptr,
                                const MultipathMapper::match_fanouts_t* fanout_breaks = nullptr);
                                
        /// Same as the previous constructor, but construct injection_trans implicitly and temporarily.
        MultipathAlignmentGraph(const HandleGraph& graph, MultipathMapper::memcluster_t& hits,
                                const unordered_map<id_t, pair<id_t, bool>>& projection_trans,
                                vector<size_t>& path_node_provenance,
                                size_t max_branch_trim_length = 0, gcsa::GCSA* gcsa = nullptr,
                                const MultipathMapper::match_fanouts_t* fanout_breaks = nullptr);
        
        /// Same as the previous constructor, but construct injection_trans implicitly and temporarily
        /// using a lambda for a projector
        MultipathAlignmentGraph(const HandleGraph& graph, MultipathMapper::memcluster_t& hits,
                                const function<pair<id_t, bool>(id_t)>& project,
                                vector<size_t>& path_node_provenance,
                                size_t max_branch_trim_length = 0, gcsa::GCSA* gcsa = nullptr,
                                const MultipathMapper::match_fanouts_t* fanout_breaks = nullptr);
        
        /// Construct a graph of the reachability between aligned chunks in a linearized
        /// path graph. Produces a graph with reachability edges.
        MultipathAlignmentGraph(const HandleGraph& graph, const vector<pair<pair<string::const_iterator, string::const_iterator>, Path>>& path_chunks,
                                const Alignment& alignment, const function<pair<id_t, bool>(id_t)>& project,
                                const unordered_multimap<id_t, pair<id_t, bool>>& injection_trans, bool realign_Ns = true,
                                bool preserve_tail_anchors = false);
       
        /// Same as the previous constructor, but construct injection_trans implicitly and temporarily
        MultipathAlignmentGraph(const HandleGraph& graph, const vector<pair<pair<string::const_iterator, string::const_iterator>, Path>>& path_chunks,
                                const Alignment& alignment, const unordered_map<id_t, pair<id_t, bool>>& projection_trans, bool realign_Ns = true,
                                bool preserve_tail_anchors = false);
        
        /// Same as the previous constructor, but construct injection_trans implicitly and temporarily
        /// and using a lambda for a projector
        MultipathAlignmentGraph(const HandleGraph& graph, const vector<pair<pair<string::const_iterator, string::const_iterator>, Path>>& path_chunks,
                                const Alignment& alignment, const function<pair<id_t, bool>(id_t)>& project, bool realign_Ns = true,
                                bool preserve_tail_anchors = false);
        
        /// Make a multipath alignment graph using the path of a single-path alignment. Only
        /// one of snarl_manager and dist_index need be supplied.
        MultipathAlignmentGraph(const HandleGraph& graph, const Alignment& alignment, SnarlManager* snarl_manager,
                                MinimumDistanceIndex* dist_index, size_t max_snarl_cut_size,
                                const function<pair<id_t, bool>(id_t)>& project,
                                const unordered_multimap<id_t, pair<id_t, bool>>& injection_trans);
        
        /// Same as the previous constructor, but construct injection_trans implicitly and temporarily
        MultipathAlignmentGraph(const HandleGraph& graph, const Alignment& alignment, SnarlManager* snarl_manager,
                                MinimumDistanceIndex* dist_index, size_t max_snarl_cut_size,
                                const unordered_map<id_t, pair<id_t, bool>>& projection_trans);
        
        /// Same as the previous constructor, but construct injection_trans implicitly and temporarily
        /// using a function instead of a map
        MultipathAlignmentGraph(const HandleGraph& graph, const Alignment& alignment, SnarlManager* snarl_manager,
                                MinimumDistanceIndex* dist_index, size_t max_snarl_cut_size,
                                const function<pair<id_t, bool>(id_t)>& project);
        
        ~MultipathAlignmentGraph();
        
        /// Fills input vector with node indices of a topological sort. 
        /// Reachability edges must be in the graph.
        void topological_sort(vector<size_t>& order_out);
        
        /// Removes non-softclip indels from path nodes. Does not update edges--should be called
        /// prior to adding computing edges.  If preserve tail anchors is true, then a null anchor (no
        /// bases and no path) will be preserved if the read segment orresponds to the beginning or
        /// end of the alignment sequence.
        void trim_hanging_indels(const Alignment& alignment, bool trim_Ns = true, bool preserve_tail_anchors = false);
        
        /// Removes all transitive edges from graph (reduces to minimum equivalent graph),
        /// except for edges between path nodes that abut either on the graph or read. These
        /// edges often correspond to overlap breaks in low complexity sequence, so retaining
        /// them improves alignment in low-complexity regions like STR expansions.
        /// Note: reorders internal representation of adjacency lists.
        /// Reachability edges must be in the graph.
        void remove_transitive_edges(const vector<size_t>& topological_order);
        
        /// Removes nodes and edges that are not part of any path that has an estimated score
        /// within some amount of the highest scoring path. Reachability edges must be present.
        void prune_to_high_scoring_paths(const Alignment& alignment, const GSSWAligner* aligner,
                                         double max_suboptimal_score_ratio, const vector<size_t>& topological_order,
                                         function<pair<id_t, bool>(id_t)>& translator,
                                         vector<size_t>& path_node_provenance);
        
        /// Clear reachability edges, so that add_reachability_edges can be run
        /// (possibly after modifying the graph).
        void clear_reachability_edges();
        
        /// Remove the ends of paths, up to a maximum length, if they cause the path
        /// to extend past a branch point in the graph.
        void trim_to_branch_points(const HandleGraph* graph, size_t max_trim_length = 1);
        
        /// Cut the interior of snarls out of anchoring paths (and split
        /// alignment nodes accordingly) unless they are longer than the max
        /// cut size. Snarls can be stored either in a SnarlManager or a
        /// MinimumDistanceIndex (only one need be supplied).
        void resect_snarls_from_paths(SnarlManager* cutting_snarls, MinimumDistanceIndex* dist_index,
                                      const function<pair<id_t, bool>(id_t)>& project, int64_t max_snarl_cut_size = 5);
        
        /// Do some exploratory alignments of the tails of the graph, outside
        /// the outermost existing anchors, and define new anchoring paths from
        /// them. After this, you can call resect_snarls_from_paths, in order
        /// to get better coverage of possible combinations of snarl traversals
        /// in parts of the alignment that didn't originally have anchors.
        /// Produces *only* perfect match anchors, so it is still safe to use
        /// score_anchors_as_matches. The Alignment passed *must* be the same
        /// Alignment that owns the sequence into which iterators were passed
        /// when the MultipathAlignmentGraph was constructed! TODO: Shouldn't
        /// the class hold a reference to the Alignment then?
        void synthesize_tail_anchors(const Alignment& alignment, const HandleGraph& align_graph, const GSSWAligner* aligner,
                                     size_t min_anchor_size, size_t max_alt_alns, bool dynamic_alt_alns, size_t max_gap,
                                     double pessimistic_tail_gap_multiplier);
        
        /// Add edges between reachable nodes and split nodes at overlaps
        void add_reachability_edges(const HandleGraph& vg,
                                    const function<pair<id_t, bool>(id_t)>& project,
                                    const unordered_multimap<id_t, pair<id_t, bool>>& injection_trans,
                                    vector<size_t>* path_node_provenance = nullptr);
                                    
        /// Do intervening and tail alignments between the anchoring paths and
        /// store the result in a multipath_alignment_t. Reachability edges must
        /// be in the graph. The Alignment passed *must* be the same Alignment
        /// that owns the sequence into which iterators were passed when the
        /// MultipathAlignmentGraph was constructed! TODO: Shouldn't the class
        /// hold a reference to the Alignment then?
        ///
        /// Note that the output alignment may NOT be in topologically-sorted
        /// order, even if this MultipathAlignmentGraph is. You MUST sort it
        /// with topologically_order_subpaths() before trying to run DP on it.
        void align(const Alignment& alignment, const HandleGraph& align_graph, const GSSWAligner* aligner, bool score_anchors_as_matches,
                   size_t max_alt_alns, bool dynamic_alt_alns, size_t max_gap, double pessimistic_tail_gap_multiplier, size_t band_padding,
                   multipath_alignment_t& multipath_aln_out, bool allow_negative_scores = false);
        
        /// Do intervening and tail alignments between the anchoring paths and
        /// store the result in a multipath_alignment_t. Reachability edges must
        /// be in the graph. Also, choose the band padding dynamically as a
        /// function of the inter-MEM sequence and graph. The Alignment passed
        /// *must* be the same Alignment that owns the sequence into which
        /// iterators were passed when the MultipathAlignmentGraph was
        /// constructed! TODO: Shouldn't the class hold a reference to the
        /// Alignment then?
        ///
        /// Note that the output alignment may NOT be in topologically-sorted
        /// order, even if this MultipathAlignmentGraph is. You MUST sort it
        /// with topologically_order_subpaths() before trying to run DP on it.
        void align(const Alignment& alignment, const HandleGraph& align_graph, const GSSWAligner* aligner, bool score_anchors_as_matches,
                   size_t max_alt_alns, bool dynamic_alt_alns, size_t max_gap, double pessimistic_tail_gap_multiplier,
                   function<size_t(const Alignment&,const HandleGraph&)> band_padding_function,
                   multipath_alignment_t& multipath_aln_out, bool allow_negative_scores = false);
        
        /// Converts a MultipathAlignmentGraph to a GraphViz Dot representation, output to the given ostream.
        /// If given the Alignment query we are working on, can produce information about subpath iterators.
        void to_dot(ostream& out, const Alignment* alignment = nullptr) const;
        
        /// Get lists of the vg node IDs that participate in each connected component in the MultipathAlignmentGraph
        vector<vector<id_t>> get_connected_components() const;
        
        /// Does the multipath alignment graph have any nodes?
        bool empty() const;
        
        size_t size() const;
        
    private:
        
        /// Nodes representing walked MEMs in the graph
        vector<PathNode> path_nodes;
        
        /// We keep a flag for whether the reachability edges are set. This is
        /// for error checking, and is kind of a forgery (you should just check
        /// the actual edge records), but the state system we use is confusing
        /// so we want to make sure we have lots of asserts to enforce it.
        /// If this is set and you want it unset, use clear_reachability_edges().
        /// If this is unset and you want it set, use add_reachability_edges().
        bool has_reachability_edges = false;
        
        /// Trim down the given PathNode of everything except softclips.
        /// Return true if it all gets trimmed away and should be removed.
        /// Fills in removed_start_from_length and/or removed_end_from_length
        /// with the bases in the graph removed from the path on each end
        /// during trimming, if set. If preserve tail anchors is true, then a null
        /// anchor (no bases and no path) will be preserved if the read segment
        /// corresponds to the beginning or end of the alignment sequence.
        static bool trim_and_check_for_empty(const Alignment& alignment, bool trim_Ns, PathNode& path_node,
                                             bool preserve_tail_anchors, int64_t* removed_start_from_length = nullptr,
                                             int64_t* removed_end_from_length = nullptr);
        
        /// Add the path chunks as nodes to the connectivity graph
        void create_path_chunk_nodes(const HandleGraph& graph, const vector<pair<pair<string::const_iterator, string::const_iterator>, Path>>& path_chunks,
                                     const Alignment& alignment, const function<pair<id_t, bool>(id_t)>& project,
                                     const unordered_multimap<id_t, pair<id_t, bool>>& injection_trans);
        
        /// Walk out MEMs into match nodes and filter out redundant sub-MEMs
        void create_match_nodes(const HandleGraph& graph, MultipathMapper::memcluster_t& hits,
                                const function<pair<id_t, bool>(id_t)>& project,
                                const unordered_multimap<id_t, pair<id_t, bool>>& injection_trans,
                                vector<size_t>& path_node_provenance,
                                int64_t max_branch_trim_length,
                                const MultipathMapper::match_fanouts_t* fanout_breaks);
        
        /// If path nodes partially overlap, merge the sections that overlap into a single path node
        void merge_partially_redundant_match_nodes(const unordered_map<int64_t, vector<int64_t>>& node_matches,
                                                   vector<size_t>& path_node_provenance);
        
        void jitter_homopolymer_ends(const HandleGraph& graph,
                                     vector<size_t>& path_node_provenance,
                                     const MultipathMapper::memcluster_t& hits,
                                     int64_t max_branch_trim_length);
        
        /// Identifies runs of exact matches that are sub-maximal because they hit the order of the GCSA
        /// index and merges them into a single node, assumes that match nodes are sorted by length and
        /// then lexicographically by read interval, does not update edges
        void collapse_order_length_runs(const HandleGraph& graph, gcsa::GCSA* gcsa,
                                        vector<size_t>& path_node_provenance);
        
        /// Reorders adjacency list representation of edges so that they follow the indicated
        /// ordering of their target nodes
        void reorder_adjacency_lists(const vector<size_t>& order);
        
        /// Return the pessimistic gap length corresponding to a certain tail length and multiplier (proportional to
        /// the square root of the tail length)
        int64_t pessimistic_tail_gap(int64_t tail_length, double multiplier);
        
        /// Returns true if we're pointing into a snarl that we want to cut out of paths
        bool into_cutting_snarl(id_t node_id, bool is_rev,
                                SnarlManager* snarl_manager, MinimumDistanceIndex* dist_index);
        
        /// Generate alignments of the tails of the query sequence, beyond the
        /// sources and sinks. The Alignment passed *must* be the one that owns
        /// the sequence we are working on. Returns a map from tail
        /// (left=false, right=true), to a map from subpath number to all the
        /// Alignments of the tail off of that subpath. Also computes the
        /// source subpaths and adds their numbers to the given set if not
        /// null.
        /// If dynamic alignment count is also selected, can indicate a minimum number
        /// of paths that must be in the extending graph in order to do an alignment
        unordered_map<bool, unordered_map<size_t, vector<Alignment>>>
        align_tails(const Alignment& alignment, const HandleGraph& align_graph, const GSSWAligner* aligner,
                    size_t max_alt_alns, bool dynamic_alt_alns, size_t max_gap, double pessimistic_tail_gap_multiplier,
                    size_t min_paths, unordered_set<size_t>* sources = nullptr);
        
        /// If a list of aligned subsequences are identifical in a prefix/suffix, remove that
        /// prefix/suffix from all of the alignments and return it as a separate alignment.
        /// If there is no shared prefix/suffix, returns an empty path with 0 score.
        static pair<path_t, int32_t> zip_alignments(vector<pair<path_t, int32_t>>& alt_alns, bool from_left,
                                                    const Alignment& alignment, const HandleGraph& align_graph,
                                                    string::const_iterator begin, const GSSWAligner* aligner);
        
        /// Memo for the transcendental pessimistic tail gap function (thread local to maintain thread-safety)
        static thread_local unordered_map<double, vector<int64_t>> pessimistic_tail_gap_memo;
        
        /// The largest size we will memoize up to
        static const size_t tail_gap_memo_max_size;
    };
}


#endif /* multipath_alignment_graph_hpp */




