//
//  multipath_mapper.hpp
//
//
//

#ifndef multipath_mapper_hpp
#define multipath_mapper_hpp

#include "hash_map.hpp"
#include "mapper.hpp"
#include "gssw_aligner.hpp"
#include "types.hpp"
#include "multipath_alignment.hpp"
#include "xg.hpp"
#include "vg.pb.h"
#include "position.hpp"
#include "nodeside.hpp"
#include "path.hpp"
#include "edit.hpp"
#include "snarls.hpp"
#include "haplotypes.hpp"

#include <structures/suffix_tree.hpp>
#include <gbwt/gbwt.h>

using namespace std;
using namespace haplo;

namespace vg {

    
    
    class MultipathMapper : public BaseMapper  {
    public:
    
        ////////////////////////////////////////////////////////////////////////
        // Interface
        ////////////////////////////////////////////////////////////////////////
    
        MultipathMapper(xg::XG* xg_index, gcsa::GCSA* gcsa_index, gcsa::LCPArray* lcp_array,
                        haplo::ScoreProvider* haplo_score_provider = nullptr, SnarlManager* snarl_manager = nullptr);
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
        
        /// Debugging function to check that multipath alignment meets the formalism's basic
        /// invariants. Returns true if multipath alignment is valid, else false. Does not
        /// validate alignment score.
        bool validate_multipath_alignment(const MultipathAlignment& multipath_aln) const;
        
        /// Sets the minimum clustering MEM length to the approximate length that a MEM would have to be to
        /// have at most the given probability of occurring in random sequence of the same size as the graph
        void set_automatic_min_clustering_length(double random_mem_probability = 0.5);
        
        /// Map random sequences against the graph to calibrate a parameterized distribution that detects
        /// when mappings are likely to have occurred by chance
        void calibrate_mismapping_detection(size_t num_simulations = 1000, size_t simulated_read_length = 150);
        
        // parameters
        
        int64_t max_snarl_cut_size = 5;
        int32_t band_padding = 2;
        size_t max_expected_dist_approx_error = 8;
        int32_t num_alt_alns = 4;
        double mem_coverage_min_ratio = 0.5;
        double max_suboptimal_path_score_ratio = 2.0;
        size_t num_mapping_attempts = 1;
        double log_likelihood_approx_factor = 1.0;
        size_t min_clustering_mem_length = 0;
        size_t max_p_value_memo_size = 500;
        double pseudo_length_multiplier = 1.65;
        double max_mapping_p_value = 0.00001;
        bool unstranded_clustering = true;
        size_t max_rescue_attempts = 32;
        size_t secondary_rescue_attempts = 4;
        double secondary_rescue_score_diff = 1.0;
        double mapq_scaling_factor = 1.0 / 4.0;
        // There must be a ScoreProvider provided if this is true
        bool use_population_mapqs = false;
        size_t population_max_paths = 1;
        // Note that, like the haplotype scoring code, we work with recombiantion penalties in exponent form.
        double recombination_penalty = 9 * 2.3;
        size_t rescue_only_min = 128;
        size_t rescue_only_anchor_max = 16;
        size_t order_length_repeat_hit_max = 0;
        
        //static size_t PRUNE_COUNTER;
        //static size_t SUBGRAPH_TOTAL;
        
        /// We often pass around clusters of MEMs and their graph positions.
        using memcluster_t = vector<pair<const MaximalExactMatch*, pos_t>>;
        
        /// This represents a graph for a cluster, and holds a pointer to the
        /// actual extracted graph, a list of assigned MEMs, and the number of
        /// bases of read coverage that that MEM cluster provides (which serves
        /// as a priority).
        using clustergraph_t = tuple<VG*, memcluster_t, size_t>;
        
    protected:
        
        /// Wrapped internal function that allows some code paths to circumvent the current
        /// mapping quality method option.
        void multipath_map_internal(const Alignment& alignment,
                                    MappingQualityMethod mapq_method,
                                    vector<MultipathAlignment>& multipath_alns_out,
                                    size_t max_alt_mappings);
        
        /// Before the fragment length distribution has been estimated, look for an unambiguous mapping of
        /// the reads using the single ended routine. If we find one record the fragment length and report
        /// the pair, if we don't find one, add the read pair to a buffer instead of the output vector.
        void attempt_unpaired_multipath_map_of_pair(const Alignment& alignment1, const Alignment& alignment2,
                                                    vector<pair<MultipathAlignment, MultipathAlignment>>& multipath_aln_pairs_out,
                                                    vector<pair<Alignment, Alignment>>& ambiguous_pair_buffer);
        
        /// Extracts a section of graph at a distance from the MultipathAlignment based on the fragment length
        /// distribution and attempts to align the other paired read to it. If rescuing forward, assumes the
        /// provided MultipathAlignment is the first read and vice versa if rescuing backward. Rescue constructs
        /// a conventional local alignment with gssw and converts the Alignment to a MultipathAlignment. The
        /// MultipathAlignment will be stored in the object passed by reference as an argument.
        bool attempt_rescue(const MultipathAlignment& multipath_aln, const Alignment& other_aln,
                            bool rescue_forward, MultipathAlignment& rescue_multipath_aln);
        
        
        /// After clustering MEMs, extracting graphs, and assigning hits to cluster graphs, perform
        /// multipath alignment
        void align_to_cluster_graphs(const Alignment& alignment,
                                     MappingQualityMethod mapq_method,
                                     vector<clustergraph_t>& cluster_graphs,
                                     vector<MultipathAlignment>& multipath_alns_out,
                                     size_t max_alt_mappings);
        
        /// After clustering MEMs, extracting graphs, assigning hits to cluster graphs, and determining
        /// which cluster graph pairs meet the fragment length distance constraints, perform multipath
        /// alignment
        void align_to_cluster_graph_pairs(const Alignment& alignment1, const Alignment& alignment2,
                                          vector<clustergraph_t>& cluster_graphs1,
                                          vector<clustergraph_t>& cluster_graphs2,
                                          vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs,
                                          vector<pair<MultipathAlignment, MultipathAlignment>>& multipath_aln_pairs_out,
                                          size_t max_alt_mappings,
                                          OrientedDistanceClusterer::paths_of_node_memo_t* paths_of_node_memo = nullptr,
                                          OrientedDistanceClusterer::oriented_occurences_memo_t* oriented_occurences_memo = nullptr,
                                          OrientedDistanceClusterer::handle_memo_t* handle_memo = nullptr);
        
        /// Align the read ends independently, but also try to form rescue alignments for each from
        /// the other. Return true if output obeys pair consistency and false otherwise.
        bool align_to_cluster_graphs_with_rescue(const Alignment& alignment1, const Alignment& alignment2,
                                                 vector<clustergraph_t>& cluster_graphs1,
                                                 vector<clustergraph_t>& cluster_graphs2,
                                                 bool block_rescue_from_1, bool block_rescue_from_2,
                                                 vector<pair<MultipathAlignment, MultipathAlignment>>& multipath_aln_pairs_out,
                                                 vector<pair<pair<size_t, size_t>, int64_t>>& pair_distances,
                                                 size_t max_alt_mappings);
        
        /// Use the rescue routine on strong suboptimal clusters to see if we can find a good secondary
        void attempt_rescue_for_secondaries(const Alignment& alignment1, const Alignment& alignment2,
                                            vector<clustergraph_t>& cluster_graphs1,
                                            vector<clustergraph_t>& cluster_graphs2,
                                            vector<pair<MultipathAlignment, MultipathAlignment>>& multipath_aln_pairs_out,
                                            vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs);
        
        /// Cluster and extract subgraphs for (possibly) only one end, meant to be a non-repeat, and use them to rescue
        /// an alignment for the other end, meant to be a repeat
        void attempt_rescue_of_repeat_from_non_repeat(const Alignment& alignment1, const Alignment& alignment2,
                                                      const vector<MaximalExactMatch>& mems1, const vector<MaximalExactMatch>& mems2,
                                                      bool do_repeat_rescue_from_1, bool do_repeat_rescue_from_2,
                                                      vector<memcluster_t>& clusters1, vector<memcluster_t>& clusters2,
                                                      vector<clustergraph_t>& cluster_graphs1, vector<clustergraph_t>& cluster_graphs2,
                                                      vector<pair<MultipathAlignment, MultipathAlignment>>& multipath_aln_pairs_out,
                                                      vector<pair<pair<size_t, size_t>, int64_t>>& pair_distances, size_t max_alt_mappings,
                                                      OrientedDistanceClusterer::paths_of_node_memo_t* paths_of_node_memo = nullptr,
                                                      OrientedDistanceClusterer::oriented_occurences_memo_t* oriented_occurences_memo = nullptr,
                                                      OrientedDistanceClusterer::handle_memo_t* handle_memo = nullptr);
        
        /// Merge the rescued mappings into the output vector and deduplicate pairs
        void merge_rescued_mappings(vector<pair<MultipathAlignment, MultipathAlignment>>& multipath_aln_pairs_out,
                                    vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs,
                                    vector<pair<MultipathAlignment, MultipathAlignment>>& rescued_multipath_aln_pairs,
                                    vector<pair<pair<size_t, size_t>, int64_t>>& rescued_cluster_pairs) const;
        
        /// Extracts a subgraph around each cluster of MEMs that encompasses any
        /// graph position reachable (according to the Mapper's aligner) with
        /// local alignment anchored at the MEMs. If any subgraphs overlap, they
        /// are merged into one subgraph. Returns a vector of all the merged
        /// cluster subgraphs, their MEMs assigned from the mems vector
        /// according to the MEMs' hits, and their read coverages in bp. The
        /// caller must delete the VG objects produced!
        vector<clustergraph_t> query_cluster_graphs(const Alignment& alignment,
                                                    const vector<MaximalExactMatch>& mems,
                                                    const vector<memcluster_t>& clusters);
        
        /// If there are any MultipathAlignments with multiple connected components, split them
        /// up and add them to the return vector
        void split_multicomponent_alignments(vector<MultipathAlignment>& multipath_alns_out) const;
        
        /// If there are any MultipathAlignments with multiple connected components, split them
        /// up and add them to the return vector, also measure the distance between them and add
        /// a record to the cluster pairs vector
        void split_multicomponent_alignments(vector<pair<MultipathAlignment, MultipathAlignment>>& multipath_aln_pairs_out,
                                             vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs) const;
        
        
        /// Make a multipath alignment of the read against the indicated graph and add it to
        /// the list of multimappings.
        void multipath_align(const Alignment& alignment, VG* vg,
                             memcluster_t& graph_mems,
                             MultipathAlignment& multipath_aln_out) const;
        
        /// Remove the full length bonus from all source or sink subpaths that received it
        void strip_full_length_bonuses(MultipathAlignment& mulipath_aln) const;
        
        /// Sorts mappings by score and store mapping quality of the optimal alignment in the MultipathAlignment object
        void sort_and_compute_mapping_quality(vector<MultipathAlignment>& multipath_alns, MappingQualityMethod mapq_method) const;
        
        /// Sorts mappings by score and store mapping quality of the optimal alignment in the MultipathAlignment object
        /// If there are ties between scores, breaks them by the expected distance between pairs as computed by the
        /// OrientedDistanceClusterer::cluster_pairs function (modified cluster_pairs vector)
        void sort_and_compute_mapping_quality(vector<pair<MultipathAlignment, MultipathAlignment>>& multipath_aln_pairs,
                                              vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs) const;
        
        /// Computes the log-likelihood of a given fragment length in the trained distribution
        double fragment_length_log_likelihood(int64_t length) const;
        
        /// Computes the number of read bases a cluster of MEM hits covers.
        static int64_t read_coverage(const memcluster_t& mem_hits);
        
        /// Would an alignment this good be expected against a graph this big by chance alone
        bool likely_mismapping(const MultipathAlignment& multipath_aln);
        
        /// A scaling of a score so that it approximately follows the distribution of the longest match in p-value test
        size_t pseudo_length(const MultipathAlignment& multipath_aln) const;
        
        /// The approximate p-value for a match length of the given size against the current graph
        double random_match_p_value(size_t match_length, size_t read_length);
        
        /// Compute the approximate distance between two multipath alignments
        int64_t distance_between(const MultipathAlignment& multipath_aln_1, const MultipathAlignment& multipath_aln_2,
                                 bool full_fragment = false, bool forward_strand = false) const;
        
        /// Are two multipath alignments consistently placed based on the learned fragment length distribution?
        bool are_consistent(const MultipathAlignment& multipath_aln_1, const MultipathAlignment& multipath_aln_2) const;
        
        /// Is this a consistent inter-pair distance based on the learned fragment length distribution?
        bool is_consistent(int64_t distance) const;
        
        /// Computes the Z-score of the number of matches against an equal length random DNA string.
        double read_coverage_z_score(int64_t coverage, const Alignment& alignment) const;
        
        /// Return true if any of the initial positions of the source Subpaths are shared between the two
        /// multipath alignments
        bool share_start_position(const MultipathAlignment& multipath_aln_1, const MultipathAlignment& multipath_aln_2) const;
        
        /// Get a thread_local RRMemo with these parameters
        haploMath::RRMemo& get_rr_memo(double recombination_penalty, size_t population_size) const;
        
        /// Detects if each pair can be assigned to a consistent strand of a path, and if not removes them. Also
        /// inverts the distances in the cluster pairs vector according to the strand
        void establish_strand_consistency(vector<pair<MultipathAlignment, MultipathAlignment>>& multipath_aln_pairs,
                                          vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs,
                                          OrientedDistanceClusterer::paths_of_node_memo_t* paths_of_node_memo = nullptr,
                                          OrientedDistanceClusterer::oriented_occurences_memo_t* oriented_occurences_memo = nullptr,
                                          OrientedDistanceClusterer::handle_memo_t* handle_memo = nullptr);
        
        SnarlManager* snarl_manager;
        
        /// Memos used by population model
        static thread_local unordered_map<pair<double, size_t>, haploMath::RRMemo> rr_memos;
        
        // a memo for the transcendental p-value function (thread local to maintain threadsafety)
        static thread_local unordered_map<pair<size_t, size_t>, double> p_value_memo;
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
        
        /// Construct a graph of the reachability between MEMs in a DAG-ified graph. Removes redundant
        /// sub-MEMs. Assumes that the cluster is sorted by primarily length and secondarily lexicographically
        /// by read interval. Optionally cuts snarl interiors from the paths and splits nodes accordingly
        MultipathAlignmentGraph(VG& vg, const MultipathMapper::memcluster_t& hits,
                                const unordered_map<id_t, pair<id_t, bool>>& projection_trans, gcsa::GCSA* gcsa = nullptr,
                                SnarlManager* cutting_snarls = nullptr, int64_t max_snarl_cut_size = 5);
        
        ~MultipathAlignmentGraph();
        
        /// Fills input vector with node indices of a topological sort
        void topological_sort(vector<size_t>& order_out);
        
        /// Removes all transitive edges from graph (reduces to minimum equivalent graph)
        /// Note: reorders internal representation of adjacency lists
        void remove_transitive_edges(const vector<size_t>& topological_order);
        
        /// Removes nodes and edges that are not part of any path that has an estimated score
        /// within some amount of the highest scoring path
        void prune_to_high_scoring_paths(const Alignment& alignment, const BaseAligner* aligner,
                                         double MultipathAlignmentGraph, const vector<size_t>& topological_order);
        
        /// Reorders adjacency list representation of edges so that they follow the indicated
        /// ordering of their target nodes
        void reorder_adjacency_lists(const vector<size_t>& order);
        
        /// Nodes representing walked MEMs in the graph
        vector<ExactMatchNode> match_nodes;
        
    private:
        
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
        
    };
}



#endif /* multipath_mapper_hpp */
