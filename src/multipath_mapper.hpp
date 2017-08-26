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
#include "multipath_alignment.hpp"
#include "utility.hpp"
#include "xg.hpp"
#include "vg.pb.h"
#include "position.hpp"
#include "path.hpp"
#include "edit.hpp"
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
                           vector<MultipathAlignment>& multipath_alns_out,
                           size_t max_alt_mappings);
                           
        /// Map a paired read to the graph and make paired multipath alignments. Assumes reads are on the
        /// same strand of the DNA/RNA molecule. If the fragment length distribution is still being estimated
        /// and the pair cannot be mapped unambiguously, adds the reads to a buffer for ambiguous pairs and
        /// does not output any multipath alignments.
        void multipath_map_paired(const Alignment& alignment1, const Alignment& alignment2,
                                  vector<pair<MultipathAlignment, MultipathAlignment>>& multipath_aln_pairs_out,
                                  vector<pair<Alignment, Alignment>>& ambiguous_pair_buffer,
                                  size_t max_alt_mappings);
        
        /// Prune path anchors from a multipath alignment if their approximate likelihood is this much
        /// lower than the optimal anchors
        void set_suboptimal_path_likelihood_ratio(double maximum_acceptable_ratio = 10000.0);
        
        /// Debugging function to check that multipath alignment meets the formalism's basic
        /// invariants. Returns true if multipath alignment is valid, else false. Does not
        /// validate alignment score.
        bool validate_multipath_alignment(const MultipathAlignment& multipath_aln) const;
        
        // parameters
        
        int64_t max_snarl_cut_size = 5;
        int32_t band_padding = 2;
        size_t max_expected_dist_approx_error = 8;
        int32_t num_alt_alns = 4;
        double mem_coverage_min_ratio = 0.5;
        int32_t max_suboptimal_path_score_diff = 20;
        size_t num_mapping_attempts = 1;
        double log_likelihood_approx_factor = 1.0;
        
        //static size_t PRUNE_COUNTER;
        //static size_t SUBGRAPH_TOTAL;
        
    private:
        
        /// Wrapped internal function that allows some code paths to circumvent the current
        /// mapping quality method option.
        void multipath_map_internal(const Alignment& alignment,
                                    MappingQualityMethod mapq_method,
                                    vector<MultipathAlignment>& multipath_alns_out,
                                    size_t max_alt_mappings);
        
        /// After clustering MEMs, extracting graphs, and assigning hits to cluster graphs, perform
        /// multipath alignment
        void align_to_cluster_graphs(const Alignment& alignment,
                                     MappingQualityMethod mapq_method,
                                     vector<tuple<VG*, vector<pair<const MaximalExactMatch*, pos_t>>, size_t>>& cluster_graphs,
                                     vector<MultipathAlignment>& multipath_alns_out,
                                     size_t max_alt_mappings);
        
        /// Extracts a subgraph around each cluster of MEMs that encompasses any
        /// graph position reachable with local alignment anchored at the MEMs.
        /// If any subgraphs overlap, they are merged into one subgraph. Also
        /// returns a map from each cluster subgraph to a vector containing the
        /// MEM hits in that graph. Also keeps track of which output graph each
        /// input cluster contributed to, in merged_into_out.
        void query_cluster_graphs(const Alignment& alignment,
                                  const vector<MaximalExactMatch>& mems,
                                  const vector<vector<pair<const MaximalExactMatch*, pos_t>>>& clusters,
                                  vector<tuple<VG*, vector<pair<const MaximalExactMatch*, pos_t>>, size_t>>& cluster_graphs_out);
        
        /// Make a multipath alignment of the read against the indicated graph and add it to
        /// the list of multimappings.
        void multipath_align(const Alignment& alignment, VG* vg,
                             vector<pair<const MaximalExactMatch*, pos_t>>& graph_mems,
                             MultipathAlignment& multipath_aln_out) const;
        
        /// Computes the number of read bases a cluster of MEM hits covers.
        int64_t read_coverage(const vector<pair<const MaximalExactMatch*, pos_t>>& mem_hits) const;
        
        /// Reorders the Subpaths in the MultipathAlignment to be in topological order (required by .proto specifications)
        void topologically_order_subpaths(MultipathAlignment& multipath_aln) const;
        
        /// Remove the full length bonus from all source or sink subpaths that received it
        void strip_full_length_bonuses(MultipathAlignment& mulipath_aln) const;
        
        /// Sorts mappings by score and store mapping quality of the optimal alignment in the MultipathAlignment object
        void sort_and_compute_mapping_quality(vector<MultipathAlignment>& multipath_alns, MappingQualityMethod mapq_method) const;
        
        /// Sorts mappings by score and store mapping quality of the optimal alignment in the MultipathAlignment object
        void sort_and_compute_mapping_quality(vector<pair<MultipathAlignment, MultipathAlignment>>& multipath_aln_pairs) const;
        
        /// Computes the Z-score of the number of matches against an equal length random DNA string.
        double read_coverage_z_score(int64_t coverage, const Alignment& alignment) const;
        
        SnarlManager* snarl_manager;
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
