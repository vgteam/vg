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
        
        ~MultipathAlignmentGraph();
        
        /// Fills input vector with node indices of a topological sort. 
        /// Reachability edges must be in the graph.
        void topological_sort(vector<size_t>& order_out);
        
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
        
        /// Cut the interior of snarls out of anchoring paths (and split
        /// alignment nodes accordingly) unless they are longer than the max
        /// cut size. Reachability edges must be cleared.
        void resect_snarls_from_paths(SnarlManager* cutting_snarls, const unordered_map<id_t, pair<id_t, bool>>& projection_trans,
                                      int64_t max_snarl_cut_size = 5);
        
        /// Add edges between reachable nodes and split nodes at overlaps
        void add_reachability_edges(VG& vg,
                                    const unordered_map<id_t, pair<id_t, bool>>& projection_trans,
                                    const unordered_multimap<id_t, pair<id_t, bool>>& injection_trans);
                                    
        /// Do intervening and tail alignments between the anchoring paths and store the result
        /// in a MultipathAlignment. Reachability edges must be in the graph.
        /// Reachability edges must be in the graph.
        void align(const Alignment& alignment, VG& align_graph, BaseAligner* aligner, bool score_anchors_as_matches,
                   size_t num_alt_alns, size_t band_padding, MultipathAlignment& multipath_aln_out);
                   
        /// Converts a MultipathAlignmentGraph to a GraphViz Dot representation, output to the given ostream.
        void to_dot(ostream& out) const;
        
        /// Get lists of the vg node IDs that participate in each connected component in the MultipathAlignmentGraph
        vector<vector<id_t>> get_connected_components() const;
        
        /// Does the multipath alignment xgraph have any nodes?
        bool empty();
        
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
        
        // Reorders the nodes of a Protobuf graph in topological order, flips doubly reversing edges,
        // and removes empty sequence nodes (invariants required for gssw alignment)
        // TODO: this is duplicative with VG::lazy_sort, but I don't want to construct a VG here
        void groom_graph_for_gssw(Graph& graph);
    };
}


#endif /* multipath_alignment_graph_hpp */




