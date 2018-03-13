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
        
        /// Construct a graph of the reachability between MEMs in a DAG-ified graph. Removes redundant
        /// sub-MEMs. Assumes that the cluster is sorted by primarily length and secondarily lexicographically
        /// by read interval. Optionally cuts snarl interiors from the paths and splits nodes accordingly
        MultipathAlignmentGraph(VG& vg, const MultipathMapper::memcluster_t& hits,
                                const unordered_map<id_t, pair<id_t, bool>>& projection_trans, gcsa::GCSA* gcsa = nullptr,
                                SnarlManager* cutting_snarls = nullptr, int64_t max_snarl_cut_size = 5);
        
        /// Construct a graph of the reachability between MEMs in a linearized path graph.
        MultipathAlignmentGraph(VG& vg, const vector<pair<pair<string::const_iterator, string::const_iterator>, Path>>& path_chunks,
                                const unordered_map<id_t, pair<id_t, bool>>& projection_trans);
        
        ~MultipathAlignmentGraph();
        
        /// Fills input vector with node indices of a topological sort
        void topological_sort(vector<size_t>& order_out);
        
        /// Removes all transitive edges from graph (reduces to minimum equivalent graph)
        /// Note: reorders internal representation of adjacency lists
        void remove_transitive_edges(const vector<size_t>& topological_order);
        
        /// Removes nodes and edges that are not part of any path that has an estimated score
        /// within some amount of the highest scoring path
        void prune_to_high_scoring_paths(const Alignment& alignment, const BaseAligner* aligner,
                                         double max_suboptimal_score_ratio, const vector<size_t>& topological_order);
        
        /// Do intervening and tail alignments between the anchoring paths and store the result
        /// in a MultipathAlignment
        void align(const Alignment& alignment, VG& align_graph, BaseAligner* aligner, bool score_anchors_as_matches,
                   size_t num_alt_alns, size_t band_padding, MultipathAlignment& multipath_aln_out);
        
    private:
        
        /// Nodes representing walked MEMs in the graph
        vector<PathNode> path_nodes;
        
        /// Walk out MEMs into match nodes and filter out redundant sub-MEMs
        void create_match_nodes(VG& vg, const MultipathMapper::memcluster_t& hits,
                                const unordered_map<id_t, pair<id_t, bool>>& projection_trans,
                                const unordered_multimap<id_t, pair<id_t, bool>>& injection_trans);
        
        /// Identifies runs of exact matches that are sub-maximal because they hit the order of the GCSA
        /// index and merges them into a single node, assumes that match nodes are sorted by length and
        /// then lexicographically by read interval, does not update edges
        void collapse_order_length_runs(VG& vg, gcsa::GCSA* gcsa);
        
        /// Cut the interior of snarls out of anchoring paths unless they are longer than the
        /// max cut size
        void resect_snarls_from_paths(SnarlManager* cutting_snarls, const unordered_map<id_t, pair<id_t, bool>>& projection_trans,
                                      int64_t max_snarl_cut_size);
        
        /// Add edges between reachable nodes and split nodes at overlaps
        void add_reachability_edges(VG& vg, const MultipathMapper::memcluster_t& hits,
                                    const unordered_map<id_t, pair<id_t, bool>>& projection_trans,
                                    const unordered_multimap<id_t, pair<id_t, bool>>& injection_trans);
        
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




