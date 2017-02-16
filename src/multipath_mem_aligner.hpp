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
#include <queue>
#include <list>

#include "hash_map.hpp"
#include "suffix_tree.hpp"
#include "mapper.hpp"
#include "gssw_aligner.hpp"
#include "types.hpp"
#include "utility.hpp"
#include "xg.hpp"
#include "vg.pb.h"
#include "position.hpp"

using namespace std;

namespace vg {
    
    class MultipathMEMAligner {
    public:
        MultipathMEMAligner(const Alignment& alignment,
                            const vector<MaximalExactMatch>& mems,
                            const QualAdjAligner& aligner,
                            xg::XG& xgindex,
                            LRUCache<id_t, Node>& node_cache,
                            SnarlManager& snarl_manager,
                            int8_t full_length_bonus,
                            size_t num_pruning_tracebacks = 5);
        
    private:
        class MultipathMEMNode;
        class MultipathMEMEdge;
        class TracebackManager;
        
        struct DPScoreComparator;
        
        /// The longest gap detectable from one side of a MEM without soft-clipping
        inline size_t longest_detectable_gap(const Alignment& alignment,
                                             const MaximalExactMatch& mem,
                                             const int32_t mem_score,
                                             const int8_t match_score,
                                             const bool gap_to_right,
                                             const QualAdjAligner& aligner,
                                             const int8_t full_length_bonus);
        
        
        void init_first_pass_graph(const Alignment& alignment,
                                   const vector<MaximalExactMatch>& mems,
                                   const QualAdjAligner& aligner,
                                   const xg::XG& xgindex,
                                   int8_t full_length_bonus);
        
        /// Fills input vector with node indices of a topological sort
        void topological_order(vector<size_t>& order_out);
        
        /// Fills input vectors with indices of source and sink nodes
        void identify_sources_and_sinks(vector<size_t>& sources_out, vector<size_t>& sinks_out);
        
        /// Identify connected components in the graph
        void connected_components(vector<vector<size_t>>& components_out);
        
        /// Perform dynamic programming
        void perform_dp();
        
        /// Returns the indices of nodes and edges on a given number of the top scoring tracebacks
        void prune_to_nodes_on_tracebacks(size_t num_tracebacks);
        
        /// Subsets the graph to the nodes and edges with the indices given
        void prune_graph(const unordered_set<size_t>& node_idxs,
                         const unordered_set<pair<size_t, size_t>>& edge_idxs);
        
        /// Adds Paths to nodes and prunes away any nodes that turn out to be redundant sub-MEMs
        void query_node_matches(const Alignment& alignment, xg::XG& xgindex, LRUCache<id_t, Node>& node_cache);
        
        /// Cut out any part of a match node in the interior of an ultrabubble
        void remove_snarls(SnarlManager& snarl_manager);
        
        /// Master list of the nodes in the exact match graph
        vector<MultipathMEMNode> nodes;
        
        /// Index of the connect subgraphs for each each
        unordered_map<pair<size_t, size_t>, Graph> edge_subgraphs;
    };
    
    class MultipathMEMAligner::MultipathMEMNode {
    public:
        MultipathMEMNode(string::const_iterator begin, string::const_iterator end, pos_t start_pos,
                         size_t offset, int32_t score) :
                         begin(begin), end(end), start_pos(start_pos), offset(offset), score(score) {}
        MultipathMEMNode() = default;
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
        
        /// The walk of the exact match through the graph
        Path path;
    };
    
    class MultipathMEMAligner::MultipathMEMEdge {
    public:
        MultipathMEMEdge(size_t to_idx, int32_t weight) : to_idx(to_idx), weight(weight) {}
        MultipathMEMEdge() = default;
        ~MultipathMEMEdge() = default;
        
        /// Index of the node that the edge points to
        size_t to_idx;
        
        /// Weight for dynamic programming
        int32_t weight;
        
    };
    
    struct MultipathMEMAligner::DPScoreComparator {
    private:
        const vector<MultipathMEMNode>& nodes;
    public:
        DPScoreComparator() = delete;
        DPScoreComparator(const vector<MultipathMEMNode>& nodes) : nodes(nodes) {}
        inline bool operator()(const size_t i, const size_t j) {
            return nodes[i].dp_score < nodes[j].dp_score;
        }
        
    };
    
    class MultipathMEMAligner::TracebackManager {
    public:
        TracebackManager(const vector<MultipathMEMNode>& nodes, size_t max_num_tracebacks);
        
        /// Advance to the next traceback
        inline void next();
        /// Returns true if there are more tracebacks in the stack
        inline bool has_next();
        
        /// Get the starting node of the next traceback (call only once per traceback)
        inline size_t get_traceback_start();
        /// The score of the current traceback
        inline int32_t curr_traceback_score();
        /// Mark a node as being included in a traceback
        inline void mark_traced(const size_t node_idx);
        /// Returns true if this is the node where the current traceback takes its next deflection
        /// from the optimal traceback.
        inline bool at_next_deflection(const size_t node_idx);
        /// Get the node to take at the next deflection (call only once per deflection)
        inline size_t deflect();
        /// Check if a local alternate traceback should be inserted into the stack
        inline void propose_deflection(const size_t from, const size_t to, const int32_t score);
        
        // not crazy about exposing this, but it lets me save recomputing it
        unordered_set<size_t> traced_node_idxs;
        
    private:
        struct Deflection;
        
        /// Add the optimal traceback to the stack
        inline void init_alt_traceback_stack();
        /// Choose insert an alternate traceback into the stack, maintaining invariants
        inline void insert_traceback(const vector<Deflection>& traceback_prefix,
                                     const size_t from, const size_t to, const int32_t score);
        /// Check whether any tracebacks that do not branch from previous tracebacks should be added
        /// to the stack before the current traceback
        inline void find_local_traceback_end();
        /// Ensure that the stack does not grow beyond the maximum number of alternate tracebacks
        inline void trim_traceback_stack();
        
        /// Reference to the nodes in the graph
        const vector<MultipathMEMNode>& nodes;
        
        /// The highest number of alternate tracebacks we will compute
        size_t max_num_tracebacks;
        /// Stack in score-sorted order of alternate tracebacks
        list<pair<vector<Deflection>, int32_t>> alt_tracebacks;
        /// The element of the stack we are currently tracing
        typename list<pair<vector<Deflection>, int32_t>>::iterator curr_traceback;
        /// The next place we will differ from the optimal traceback path
        typename vector<Deflection>::iterator curr_deflxn;
        /// Priority queue of nodes ordered by DP score for finding non-branching tracebacks
        priority_queue<size_t, vector<size_t>, DPScoreComparator> traceback_end_candidates;
    };
    
    struct MultipathMEMAligner::TracebackManager::Deflection {
        Deflection(const size_t from, const size_t to) : from_idx(from), to_idx(to) {}
        /// Node index where traceback differs from optimal traceback
        const size_t from_idx;
        /// Node index that traceback goes to instead of optimal traceback
        const size_t to_idx;
    };}

#endif /* multipath_mem_aligner_hpp */
