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
    class MultipathClusterer {
        
    public:
        MultipathClusterer(const Alignment& alignment,
                           const vector<MaximalExactMatch>& mems,
                           const QualAdjAligner& aligner,
                           xg::XG& xgindex,
                           int8_t full_length_bonus = 0,
                           size_t num_pruning_tracebacks = 5);
        
        /// Returns a vector of clusters. Each cluster is represented by a boolean indicating whether
        /// it is from the forward (false) or reverse complement (true) strand of the read, and a
        /// vector of MEM hits. Each hit contains a pointer to the original MEM and the position of
        /// that particular hit in the graph.
        vector<pair<bool, vector<pair<const MaximalExactMatch*, pos_t>>>> clusters();
        
    private:
        class MPCNode;
        class MPCEdge;
        class TracebackManager;
        
        struct DPScoreComparator;
        
        void init_mem_graph();
        
        /// The longest gap detectable from one side of a MEM without soft-clipping
        inline size_t longest_detectable_gap(const string::const_iterator& read_pos);
        
        vector<MPCNode> nodes;
        
        const Alignment& alignment;
        const vector<MaximalExactMatch>& mems;
        const QualAdjAligner& aligner;
        int8_t full_length_bonus;
        
    };
    
    class MultipathClusterer::MPCNode {
    public:
        MPCNode(string::const_iterator begin, string::const_iterator end, pos_t start_pos,
                int32_t score) :
                begin(begin), end(end), start_pos(start_pos), score(score) {}
        MPCNode() = default;
        ~MPCNode() = default;
        
        /// Beginning position on the read
        string::const_iterator begin;
        /// End position on the read
        string::const_iterator end;
        /// Position of GCSA hit in the graph
        pos_t start_pos;
        /// Score of the exact match this node represents
        int32_t score;
        
        /// Score during dynamic programming
        int32_t forward_dp_score;
        int32_t reverse_dp_score;
        
        /// Edges from this node that are colinear with the read
        vector<MPCEdge> forward_edges_from;
        
        /// Edges to this node that are colinear with the read
        vector<MPCEdge> forward_edges_to;
        
        /// Edges from this node that are anti-colinear with the read
        vector<MPCEdge> reverse_edges_from;
        
        /// Edges to this node that are anti-colinear with the read
        vector<MPCEdge> reverse_edges_to;
    };
    
    class MultipathClusterer::MPCEdge {
    public:
        MPCEdge(size_t to_idx, int32_t weight, size_t overlap = 0) :
                to_idx(to_idx), weight(weight), overlap(overlap) {}
        MPCEdge() = default;
        ~MPCEdge() = default;
        
        /// Index of the node that the edge points to
        size_t to_idx;
        
        /// Weight for dynamic programming
        int32_t weight;
        
        /// If an overlap is required to explain this edge, how long the overlap is
        size_t overlap;
    };
    
    struct MultipathClusterer::DPScoreComparator {
    private:
        const vector<MPCNode>& nodes;
        bool reverse;
    public:
        DPScoreComparator() = delete;
        DPScoreComparator(const vector<MPCNode>& nodes, bool reverse) : nodes(nodes), reverse(reverse) {}
        ~DPScoreComparator() {}
        inline bool operator()(const size_t i, const size_t j) {
            return reverse ? nodes[i].reverse_dp_score < nodes[j].reverse_dp_score
                           : nodes[i].forward_dp_score < nodes[j].forward_dp_score;
        }
    };
}



#endif /* multipath_mapper_hpp */
