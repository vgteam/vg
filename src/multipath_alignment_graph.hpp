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

#include "vg.pb.h"
#include "vg.hpp"
#include "snarls.hpp"
#include "multipath_mapper.hpp"

namespace vg {
    
    
    // TODO: put in MultipathAlignmentGraph namespace
    class PathNode {
    public:
        string::const_iterator begin;
        string::const_iterator end;
        Path path;
        
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
        
        /// Create an identity projection translation from a DAG that did not
        /// need to be modified during dagification.
        static unordered_map<id_t, pair<id_t, bool>> create_identity_projection_trans(const VG& vg);
        
        /// Construct a graph of the reachability between MEMs in a DAG-ified
        /// graph. If a GCSA is specified, use it to collapse MEMs whose
        /// lengths bump up against the GCSA's order limit on MEM length.
        /// Produces a graph with reachability edges. Assumes that the cluster
        /// is sorted by primarily length and secondarily lexicographically by
        /// read interval.
        MultipathAlignmentGraph(VG& vg, const MultipathMapper::memcluster_t& hits,
                                const unordered_map<id_t, pair<id_t, bool>>& projection_trans,
                                const unordered_multimap<id_t, pair<id_t, bool>>& injection_trans,
                                gcsa::GCSA* gcsa = nullptr);
                                
        /// Same as the previous constructor, but construct injection_trans implicitly and temporarily.
        MultipathAlignmentGraph(VG& vg, const MultipathMapper::memcluster_t& hits,
                                const unordered_map<id_t, pair<id_t, bool>>& projection_trans,
                                gcsa::GCSA* gcsa = nullptr);
        
        /// Construct a graph of the reachability between MEMs in a linearized
        /// path graph. Produces a graph with reachability edges.
        MultipathAlignmentGraph(VG& vg, const vector<pair<pair<string::const_iterator, string::const_iterator>, Path>>& path_chunks,
                                const Alignment& alignment, const unordered_map<id_t, pair<id_t, bool>>& projection_trans,
                                const unordered_multimap<id_t, pair<id_t, bool>>& injection_trans);
       
        /// Same as the previous constructor, but construct injection_trans implicitly and temporarily
        MultipathAlignmentGraph(VG& vg, const vector<pair<pair<string::const_iterator, string::const_iterator>, Path>>& path_chunks,
                                const Alignment& alignment, const unordered_map<id_t, pair<id_t, bool>>& projection_trans);
        
        /// Make a multipath alignment graph using the path of a single-path alignment
        MultipathAlignmentGraph(VG& vg, const Alignment& alignment, SnarlManager& snarl_manager, size_t max_snarl_cut_size,
                                const unordered_map<id_t, pair<id_t, bool>>& projection_trans,
                                const unordered_multimap<id_t, pair<id_t, bool>>& injection_trans);
        
        /// Same as the previous constructor, but construct injection_trans implicitly and temporarily
        MultipathAlignmentGraph(VG& vg, const Alignment& alignment, SnarlManager& snarl_manager, size_t max_snarl_cut_size,
                                const unordered_map<id_t, pair<id_t, bool>>& projection_trans);
        
        ~MultipathAlignmentGraph();
        
        /// Fills input vector with node indices of a topological sort. 
        /// Reachability edges must be in the graph.
        void topological_sort(vector<size_t>& order_out);
        
        /// Removes non-softclip indels from path nodes. Does not update edges--should be called
        /// prior to computing edges.
        void trim_hanging_indels(const Alignment& alignment);
        
        /// Removes all transitive edges from graph (reduces to minimum equivalent graph).
        /// Note: reorders internal representation of adjacency lists.
        /// Reachability edges must be in the graph.
        void remove_transitive_edges(const vector<size_t>& topological_order);
        
        /// Removes nodes and edges that are not part of any path that has an estimated score
        /// within some amount of the highest scoring path. Reachability edges must be present.
        void prune_to_high_scoring_paths(const Alignment& alignment, const BaseAligner* aligner,
                                         double max_suboptimal_score_ratio, const vector<size_t>& topological_order);
        
        /// Clear reachability edges, so that add_reachability_edges can be run
        /// (possibly after modifying the graph).
        void clear_reachability_edges();
        
        /// Cut all PathNodes where they cross over decision points in the
        /// graph. Reachability edges must be cleared. They will not be updated
        /// to account for cutting of PathNodes.
        void cut_at_forks(VG& align_graph);
        
        /// Find all possible anchors out from the existing anchors without
        /// going through more than max_mismatches mismatches. The Alignment
        /// passed *must* be the same Alignment that owns the sequence into
        /// which iterators were passed when the MultipathAlignmentGraph was
        /// constructed!
        ///
        /// If clear_originals is set, all the original subpaths are removed,
        /// and the search is set to re-find them, assuming they are perfect
        /// matches.
        ///
        /// Does *not* create new reachability edges; they will have to be
        /// computed afterwards.
        void synthesize_anchors_by_search(const Alignment& alignment, VG& align_graph,
            size_t max_mismatches, bool clear_originals = false);
        
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
        void synthesize_tail_anchors(const Alignment& alignment, VG& align_graph, BaseAligner* aligner,
                                     size_t max_alt_alns, bool dynamic_alt_alns);
                                     
        /// Add edges between reachable nodes and split nodes at overlaps
        void add_reachability_edges(VG& vg,
                                    const unordered_map<id_t, pair<id_t, bool>>& projection_trans,
                                    const unordered_multimap<id_t, pair<id_t, bool>>& injection_trans);
                                    
        /// Cut the interior of snarls out of anchoring paths (and split
        /// alignment nodes accordingly) unless they are longer than the max
        /// cut size. Reachability edges should be present and will be updated. 
        void resect_snarls_from_paths(SnarlManager* cutting_snarls, const unordered_map<id_t, pair<id_t, bool>>& projection_trans,
                                      int64_t max_snarl_cut_size = 5);
                                    
        /// Do intervening and tail alignments between the anchoring paths and
        /// store the result in a MultipathAlignment. Reachability edges must
        /// be in the graph. The Alignment passed *must* be the same Alignment
        /// that owns the sequence into which iterators were passed when the
        /// MultipathAlignmentGraph was constructed! TODO: Shouldn't the class
        /// hold a reference to the Alignment then?
        ///
        /// Note that the output alignment may NOT be in topologically-sorted
        /// order, even if this MultipathAlignmentGraph is. You MUST sort it
        /// with topologically_order_subpaths() before trying to run DP on it.
        void align(const Alignment& alignment, VG& align_graph, BaseAligner* aligner, bool score_anchors_as_matches,
                   size_t max_alt_alns, bool dynamic_alt_alns, size_t band_padding, MultipathAlignment& multipath_aln_out, const bool allow_negative_scores = false);
        
        /// Do intervening and tail alignments between the anchoring paths and
        /// store the result in a MultipathAlignment. Reachability edges must
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
        void align(const Alignment& alignment, VG& align_graph, BaseAligner* aligner, bool score_anchors_as_matches,
                   size_t max_alt_alns, bool dynamic_alt_alns,
                   function<size_t(const Alignment&,const HandleGraph&)> band_padding_function,
                   MultipathAlignment& multipath_aln_out, const bool allow_negative_scores = false);
        
        /// Converts a MultipathAlignmentGraph to a GraphViz Dot representation, output to the given ostream.
        /// If given the Alignment query we are working on, can produce information about subpath iterators.
        void to_dot(ostream& out, const Alignment* alignment = nullptr) const;
        
        /// Get lists of the vg node IDs that participate in each connected component in the MultipathAlignmentGraph
        vector<vector<id_t>> get_connected_components() const;
        
        /// Does the multipath alignment xgraph have any nodes?
        bool empty();
        
    protected:
        
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
        /// during trimming, if set.
        static bool trim_and_check_for_empty(const Alignment& alignment, PathNode& path_node,
            int64_t* removed_start_from_length = nullptr, int64_t* removed_end_from_length = nullptr);
        
        /// Add the path chunks as nodes to the connectivity graph
        void create_path_chunk_nodes(VG& vg, const vector<pair<pair<string::const_iterator, string::const_iterator>, Path>>& path_chunks,
                                     const Alignment& alignment, const unordered_map<id_t, pair<id_t, bool>>& projection_trans,
                                     const unordered_multimap<id_t, pair<id_t, bool>>& injection_trans);
        
        /// Walk out MEMs into match nodes and filter out redundant sub-MEMs
        void create_match_nodes(VG& vg, const MultipathMapper::memcluster_t& hits,
                                const unordered_map<id_t, pair<id_t, bool>>& projection_trans,
                                const unordered_multimap<id_t, pair<id_t, bool>>& injection_trans);
        
        /// Identifies runs of exact matches that are sub-maximal because they hit the order of the GCSA
        /// index and merges them into a single node, assumes that match nodes are sorted by length and
        /// then lexicographically by read interval, does not update edges
        void collapse_order_length_runs(VG& vg, gcsa::GCSA* gcsa);
        
        /// Reorders adjacency list representation of edges so that they follow the indicated
        /// ordering of their target nodes
        void reorder_adjacency_lists(const vector<size_t>& order);
        
        /// Reorders the nodes of a Protobuf graph in topological order, flips doubly reversing edges,
        /// and removes empty sequence nodes (invariants required for gssw alignment)
        /// TODO: this is duplicative with VG::lazy_sort, but I don't want to construct a VG here
        void groom_graph_for_gssw(Graph& graph);
        
        /// Generate alignments of the tails of the query sequence, beyond the
        /// sources and sinks. The Alignment passed *must* be the one that owns
        /// the sequence we are working on. Returns a map from tail
        /// (left=false, right=true), to a map from subpath number to all the
        /// Alignments of the tail off of that subpath. Also computes the
        /// source subpaths and adds their numbers to the given set if not
        /// null.
        unordered_map<bool, unordered_map<size_t, vector<Alignment>>>
        align_tails(const Alignment& alignment, VG& align_graph, BaseAligner* aligner,
                    size_t max_alt_alns, bool dynamic_alt_alns, unordered_set<size_t>* sources = nullptr);
    };
}


#endif /* multipath_alignment_graph_hpp */




