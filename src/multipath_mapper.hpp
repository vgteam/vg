//
//  multipath_mapper.hpp
//
//
//

#ifndef multipath_mapper_hpp
#define multipath_mapper_hpp

#include <stdio.h>
#include <cmath>

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
#include "snarls.hpp"
#include "algorithms/vg_algorithms.hpp"

using namespace std;

namespace vg {

    
    
    class MultipathMapper : public BaseMapper  {
    public:
        MultipathMapper(xg::XG* xg_index, gcsa::GCSA* gcsa_index, gcsa::LCPArray* lcp_array,
                        SnarlManager* snarl_manager = nullptr);
        ~MultipathMapper();
        
        /// Map read in alignment to graph and make multipath alignments.
        void multipath_map(const Alignment& alignment,
                           list<MultipathAlignment>& multipath_alns_out,
                           size_t max_alt_alns);
        
        
    private:
        
        /// Extracts a subgraph around each cluster of MEMs that encompasses any graph position
        /// reachable with local alignment anchored at the MEMs. If any subgraphs overlap, they
        /// are merged into one subgraph. Also returns a map from each cluster subgraph to a
        /// vector containing the MEM hits in that graph.
        void query_cluster_graphs(const Alignment& alignment,
                                  const vector<MaximalExactMatch>& mems,
                                  const vector<vector<pair<const MaximalExactMatch*, pos_t>>>& clusters,
                                  vector<VG*>& cluster_graphs_out,
                                  unordered_map<VG*, vector<pair<const MaximalExactMatch*, pos_t>>>& cluster_graph_mems_out);
        
        /// Make a multipath alignment of the read against the indicated graph and add it to
        /// the list of multimappings.
        void multipath_align(const Alignment& alignment, VG* vg,
                             vector<pair<const MaximalExactMatch*, pos_t>>& graph_mems,
                             list<MultipathAlignment>& multipath_alns_out);
        
        /// Computes the number of read bases a cluster of MEM hits covers.
        int64_t read_coverage(const vector<pair<const MaximalExactMatch*, pos_t>>& mem_hits);
        
        /// Computes the Z-score of the number of matches against an equal length random DNA string.
        double read_coverage_z_score(int64_t coverage, const Alignment& alignment);
        
        SnarlManager* snarl_manager;
        
        //double z_score_cutoff = -1.0;
        size_t max_expected_dist_approx_error = 8;
        int32_t num_alt_alns = 4;
        double mem_coverage_min_ratio = 0.5;
        int32_t max_suboptimal_path_score_diff = 20;
        int64_t max_snarl_cut_size = 5;
        int32_t band_padding = 2;
    };
    
    // TODO: put in MultipathAlignmentGraph namespace
    class ExactMatchNode {
    public:
        string::const_iterator begin;
        string::const_iterator end;
        Path path;
        
        // pairs of (target index, path length)
        vector<pair<size_t, size_t>> edges;
    };
    
    // TODO: put in MultipathMapper namespace
    class MultipathAlignmentGraph {
    public:
        // removes duplicate sub-MEMs contained in parent MEMs
        MultipathAlignmentGraph(VG& vg, const vector<pair<const MaximalExactMatch*, pos_t>>& hits,
                                const unordered_multimap<id_t, pair<id_t, bool>>& injection_trans,
                                const unordered_map<id_t, pair<id_t, bool>>& projection_trans,
                                SnarlManager* cutting_snarls = nullptr, int64_t max_snarl_cut_size = 5);
        
        /// Fills input vector with node indices of a topological sort
        void topological_sort(vector<size_t>& order_out);
        
        /// Removes all transitive edges from graph (reduces to minimum equivalent graph)
        /// Note: reorders internal representation of adjacency lists
        void remove_transitive_edges(const vector<size_t>& topological_order);
        
        /// Removes nodes and edges that are not part of any path that has an estimated score
        /// within some amount of the highest scoring path
        void prune_to_high_scoring_paths(const BaseAligner& aligner, int32_t max_suboptimal_score_diff,
                                         const vector<size_t>& topological_order);
        
        /// Reorders adjacency list representation of edges so that they follow the indicated
        /// ordering of their target nodes
        void reorder_adjacency_lists(const vector<size_t>& order);
        
        vector<ExactMatchNode> match_nodes;
    };
}



#endif /* multipath_mapper_hpp */
