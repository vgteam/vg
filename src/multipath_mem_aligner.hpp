//
//  multipath_mem_aligner.hpp
//  
//
//

#ifndef multipath_mem_aligner_hpp
#define multipath_mem_aligner_hpp

#include <stdio.h>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <priority_queue>
#include <list>

#include "hash_map.hpp"
#include "mapper.hpp"
#include "gssw_aligner.hpp"
#include "types.hpp"
#include "utility.hpp"

using namespace std;

namespace vg {
    
    class MultipathMEMAligner {
    public:
        MultipathMEMAligner(const Alignment& alignment,
                            const vector<MaximalExactMatch>& mems,
                            const QualAdjAligner& aligner,
                            const function<int64_t(pos_t,pos_t)>& distance,
                            int8_t full_length_bonus);
        
    private:
        class MultipathMEMNode;
        class MultipathMEMEdge;
        class TracebackManager;
        
        struct DPScoreComparator;
        
        /// The longest gap detectable from one side of a MEM without soft-clipping
        inline size_t longest_detectable_gap(const Alignment& alignment,
                                             const MaximalExactMatch& mem,
                                             const int32_t mem_score,
                                             const bool gap_to_right,
                                             const QualAdjAligner& aligner,
                                             const int8_t full_length_bonus);
        
        
        /// Fills input vector with node indices of a topological sort
        void topological_order(vector<size_t>& order_out);
        
        /// Fills input vectors with indices of source and sink nodes
        void identify_sources_and_sinks(vector<size_t>& sources_out, vector<size_t>& sinks_out);
        
        /// Identify connected components in the graph
        void connected_components(vector<vector<size_t>>& components_out);
        
        /// Perform dynamic programming
        void fill_dp();
        
        /// Returns the 
        unordered_set<size_t> find_nodes_on_tracebacks(size_t num_tracebacks);
        
        vector<MultipathMEMNode> nodes;
    };
    
    class MultipathMEMAligner::MultipathMEMNode {
    public:
        MultipathMEMNode(string::const_iterator begin,
                         string::const_iterator end,
                         pos_t start_pos,
                         size_t offset,
                         int32_t score);
        ~MultipathMEMNode() = default;
        
        /// Beginning position on the read
        string::const_iterator begin;
        /// End position on the read
        string::const_iterator end;
        /// Position of GCSA hit in the graph
        pos_t start_pos;
        /// Distance away from the GCSA hit along the MEM
        size_t offset;
        /// Score of the exact match this node represents
        int32_t score;
        
        /// Score during dynamic programming
        int32_t dp_score;
        
        /// Edges from this node
        vector<MultipathMEMEdge> edges_from;
        
        /// Edges to this node
        vector<MultipathMEMEdge> edges_to;
    };
    
    class MultipathMEMAligner::MultipathMEMEdge {
    public:
        MultipathMEMEdge(size_t to_idx, int32_t weight);
        ~MultipathMEMEdge() = default;
        
        const size_t to_idx;
        const int32_t weight;
    };
    
    class MultipathMEMAligner::TracebackManager {
        
        struct Deflection;
        
        TracebackManager(vector<MultipathMEMNode>& nodes, size_t max_num_tracebacks);
        
        inline void next();
        inline bool has_next();
        inline size_t get_traceback_start();
        inline void mark_traced(const size_t node_idx);
        
        inline bool at_next_deflection(const size_t node_idx);
        inline int32_t curr_traceback_score();
        inline void propose_deflection(const size_t from, const size_t to, const int32_t score);
        inline void insert_traceback(const vector<Deflection>& traceback_prefix,
                                     const size_t from, const size_t to, const int32_t score);
        
        inline size_t deflect();
        
        inline void find_local_traceback_end();
        
        const vector<MultipathMEMNode>& nodes;
        
        size_t max_num_tracebacks;
        list<pair<vector<Deflection>>, int32_t> alt_tracebacks;
        typename list<pair<vector<Deflection>>, int32_t> curr_traceback;
        typename vector<Deflection>::iterator curr_deflxn;
        
        priority_queue<size_t, vector<size_t>, DPScoreComparator> traceback_end_candidates;
        unordered_set<size_t> traced_node_idxs;
        
        
    };
    
    struct MultipathMEMAligner::TracebackManager::Deflection {
        Deflection(const size_t from, const size_t to) : from_idx(from), to_idx(to) {}
        const size_t from_idx;
        const size_t to_idx;
    }
    
    struct MultipathMEMAligner::DPScoreComparator {
    private:
        const vector<MultipathMEMNode>& nodes;
    public:
        DPScoreComparator(const vector<MultipathMEMNode>& nodes) : nodes(nodes) {}
        inline bool operator()(const size_t i, const size_t j) {
            return nodes[i].dp_score < nodes[j].dp_score;
        }
        
    };
}

#endif /* multipath_mem_aligner_hpp */
