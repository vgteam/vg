//
//  multipath_mapper.hpp
//
//
//

#ifndef multipath_mapper_hpp
#define multipath_mapper_hpp

#include <stdio.h>

#include "hash_map.hpp"
#include "suffix_tree.hpp"
#include "mapper.hpp"
#include "gssw_aligner.hpp"
#include "types.hpp"
#include "utility.hpp"
#include "xg.hpp"
#include "vg.pb.h"
#include "position.hpp"
#include "path.hpp"
#include "vg_algorithms.hpp"

using namespace std;

namespace vg {
    
    class MultipathAligner  {
    public:
        
        MultipathAlignment multipath_align(const Alignment& alignment,
                                           const vector<MaximalExactMatch>& mems);
        
        
    private:
        
        xg::XG xgindex;
        LRUCache<id_t, Node> node_cache;
        
        int8_t full_length_bonus = 0;
        size_t max_strand_dist_probes = 2;
        size_t max_expected_dist_approx_error = 8;
    };
    
    
    class MultipathClusterer {
    public:
        MultipathClusterer(const Alignment& alignment,
                           const vector<MaximalExactMatch>& mems,
                           const BaseAligner& aligner,
                           xg::XG& xgindex,
                           int8_t full_length_bonus = 0,
                           size_t max_strand_dist_probes = 2,
                           size_t max_expected_dist_approx_error = 8);
        
        /// Returns a vector of clusters. Each cluster is represented a vector of MEM hits. Each hit
        /// contains a pointer to the original MEM and the position of that particular hit in the graph.
        vector<vector<pair<MaximalExactMatch* const, pos_t>>> clusters(int32_t max_qual_score = 60);
        
    private:
        class MPCNode;
        class MPCEdge;
        
        /// Fills input vector with node indices of a topological sort
        void topological_order(vector<size_t>& order_out);
        
        /// Fills input vectors with indices of source and sink nodes
        void identify_sources_and_sinks(vector<size_t>& sources_out, vector<size_t>& sinks_out);
        
        /// Identify weakly connected components in the graph
        void connected_components(vector<vector<size_t>>& components_out);
        
        /// Perform dynamic programming
        void perform_dp();
        
        /// Computes the number of bases a traceback through the MEM graph covers
        int64_t trace_read_coverage(const vector<size_t>& trace);
        
        vector<MPCNode> nodes;
        
        const BaseAligner& aligner;
    };
    
    class MultipathClusterer::MPCNode {
    public:
        MPCNode(const MaximalExactMatch& mem, pos_t start_pos, int32_t score) :
                mem(&mem), start_pos(start_pos), score(score) {}
        MPCNode() = default;
        ~MPCNode() = default;
        
        MaximalExactMatch* const mem;
        
        /// Position of GCSA hit in the graph
        pos_t start_pos;
        /// Score of the exact match this node represents
        int32_t score;
        
        /// Score during dynamic programming
        int32_t dp_score;
        
        /// Edges from this node that are colinear with the read
        vector<MPCEdge> edges_from;
        
        /// Edges to this node that are colinear with the read
        vector<MPCEdge> edges_to;
    };
    
    class MultipathClusterer::MPCEdge {
    public:
        MPCEdge(size_t to_idx, int32_t weight) :
                to_idx(to_idx), weight(weight) {}
        MPCEdge() = default;
        ~MPCEdge() = default;
        
        /// Index of the node that the edge points to
        size_t to_idx;
        
        /// Weight for dynamic programming
        int32_t weight;
    };
    
    struct MultipathClusterer::DPScoreComparator {
    private:
        const vector<MPCNode>& nodes;
    public:
        DPScoreComparator() = delete;
        DPScoreComparator(const vector<MPCNode>& nodes) : nodes(nodes) {}
        ~DPScoreComparator() {}
        inline bool operator()(const size_t i, const size_t j) {
            return nodes[i].dp_score < nodes[j].dp_score;
        }
    };
}



#endif /* multipath_mapper_hpp */
