/**
 * \file multipath_mapper.hpp
 *
 * Defines the MultipathMapper class
 */

#ifndef multipath_mapper_hpp
#define multipath_mapper_hpp

#include <algorithm>
#include <vg/vg.pb.h>
#include <structures/union_find.hpp>
#include <gbwt/gbwt.h>
#include <vg/io/edit.hpp>
#include <bdsg/hash_graph.hpp>

#include "hash_map.hpp"
#include "mapper.hpp"
#include "aligner.hpp"
#include "types.hpp"
#include "multipath_alignment.hpp"
#include "position.hpp"
#include "nodeside.hpp"
#include "path.hpp"
#include "snarls.hpp"
#include "haplotypes.hpp"
#include "min_distance.hpp"
#include "utility.hpp"
#include "annotation.hpp"
#include "path_component_index.hpp"
#include "memoizing_graph.hpp"
#include "statistics.hpp"
#include "splicing.hpp"

#include "identity_overlay.hpp"
#include "reverse_graph.hpp"
#include "split_strand_graph.hpp"
#include "dagified_graph.hpp"

#include "algorithms/extract_containing_graph.hpp"
#include "algorithms/extract_connecting_graph.hpp"
#include "algorithms/extract_extending_graph.hpp"
#include "algorithms/locally_expand_graph.hpp"
#include "algorithms/jump_along_path.hpp"
#include "algorithms/ref_path_distance.hpp"


// note: only activated for single end mapping
//#define mpmap_instrument_mem_statistics

using namespace std;
using namespace haplo;
using namespace structures;

namespace vg {
    
    class MultipathMapper : public BaseMapper  {
    public:
    
        ////////////////////////////////////////////////////////////////////////
        // Interface
        ////////////////////////////////////////////////////////////////////////
    
        MultipathMapper(PathPositionHandleGraph* graph, gcsa::GCSA* gcsa_index, gcsa::LCPArray* lcp_array,
                        haplo::ScoreProvider* haplo_score_provider = nullptr, SnarlManager* snarl_manager = nullptr,
                        MinimumDistanceIndex* distance_index = nullptr);
        ~MultipathMapper();
        
        /// Map read in alignment to graph and make multipath alignments.
        void multipath_map(const Alignment& alignment,
                           vector<multipath_alignment_t>& multipath_alns_out);
                           
        /// Map a paired read to the graph and make paired multipath alignments. Assumes reads are on the
        /// same strand of the DNA/RNA molecule. If the fragment length distribution is still being estimated
        /// and the pair cannot be mapped unambiguously, adds the reads to a buffer for ambiguous pairs and
        /// does not output any multipath alignments.
        /// Returns true if the output is properly paired, or false if it is independent end mappings.
        bool multipath_map_paired(const Alignment& alignment1, const Alignment& alignment2,
                                  vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs_out,
                                  vector<pair<Alignment, Alignment>>& ambiguous_pair_buffer);
                                  
        /// Given a mapped multipath_alignment_t, reduce it to up to
        /// max_number + 1 nonoverlapping single path alignments, with
        /// mapping qualities accounting for positional uncertainty between
        /// them.
        /// Even if the read is unmapped, there will always be at least one (possibly score 0) output alignment.
        void reduce_to_single_path(const multipath_alignment_t& multipath_aln, vector<Alignment>& alns_out, size_t max_number) const;
        
        /// Sets the minimum clustering MEM length to the approximate length that a MEM would have to be to
        /// have at most the given probability of occurring in random sequence of the same size as the graph
        void set_automatic_min_clustering_length(double random_mem_probability = 0.5);
        
        /// Map random sequences against the graph to calibrate a parameterized distribution that detects
        /// when mappings are likely to have occurred by chance
        void calibrate_mismapping_detection(size_t num_simulations, const vector<size_t>& simulated_read_lengths);
        
        /// Should be called once after construction, or any time the band padding multiplier is changed
        void init_band_padding_memo();
        
        /// Set all the aligner scoring parameters and create the stored aligner instances.
        void set_alignment_scores(int8_t match, int8_t mismatch, int8_t gap_open, int8_t gap_extend, int8_t full_length_bonus);
        
        /// Set the algner scoring parameters and create the stored aligner instances. The
        /// stream should contain a 4 x 4 whitespace-separated substitution matrix (in the
        /// order ACGT)
        void set_alignment_scores(std::istream& matrix_stream, int8_t gap_open, int8_t gap_extend, int8_t full_length_bonus);
        
        /// Set the algner scoring parameters and create the stored aligner instances. The
        /// score matrix should by a 4 x 4 array in the order (ACGT)
        void set_alignment_scores(const int8_t* score_matrix, int8_t gap_open, int8_t gap_extend, int8_t full_length_bonus);
        
        /// How big of a softclip should lead us to attempt spliced alignment?
        void set_min_softclip_length_for_splice(size_t length);
        
        /// Decide how long of a tail alignment we want before we allow its subpath to be merged
        void set_max_merge_supression_length();
        
        // parameters
        
        size_t max_branch_trim_length = 1;
        bool agglomerate_multipath_alns = false;
        int64_t max_snarl_cut_size = 5;
        size_t max_tail_merge_supress_length = 4;
        bool suppress_tail_anchors = false;
        size_t min_tail_anchor_length = 3;
        double band_padding_multiplier = 1.0;
        bool use_pessimistic_tail_alignment = false;
        double pessimistic_gap_multiplier = 0.0;
        bool restrained_graph_extraction = false;
        size_t max_expected_dist_approx_error = 8;
        int32_t num_alt_alns = 4;
        double mem_coverage_min_ratio = 0.5;
        double unused_cluster_multiplicity_mq_limit = 7.0;
        double max_suboptimal_path_score_ratio = 2.0;
        size_t num_mapping_attempts = 48;
        double log_likelihood_approx_factor = 1.0;
        size_t min_clustering_mem_length = 0;
        bool use_stripped_match_alg = false;
        size_t stripped_match_alg_strip_length = 16;
        size_t stripped_match_alg_max_length = 0;
        size_t stripped_match_alg_target_count = 5;
        bool use_fanout_match_alg = false;
        int max_fanout_base_quality = 20;
        int max_fans_out = 5;
        size_t max_p_value_memo_size = 500;
        size_t band_padding_memo_size = 2000;
        double max_exponential_rate_intercept = 0.612045;
        double max_exponential_rate_slope = 0.000555181;
        double max_exponential_shape_intercept = 12.136;
        double max_exponential_shape_slope = 0.0113637;
        double max_mapping_p_value = 0.0001;
        double max_splice_p_value = 0.001;
        double max_rescue_p_value = 0.1;
        size_t max_alt_mappings = 1;
        size_t max_single_end_mappings_for_rescue = 64;
        size_t max_rescue_attempts = 32;
        size_t plausible_rescue_cluster_coverage_diff = 5;
        size_t secondary_rescue_attempts = 4;
        double secondary_rescue_score_diff = 1.0;
        bool get_rescue_graph_from_paths = true;
        double rescue_graph_std_devs = 6.0;
        double mapq_scaling_factor = 1.0;
        bool report_group_mapq = false;
        // There must be a ScoreProvider provided, and a positive population_max_paths, if this is true
        bool use_population_mapqs = false;
        // If this is nonzero, it takes precedence over any haplotype count
        // available from the score provider or the XG index. If neither of
        // those has a haplotype count, this must be set for haplotype scoring
        // to work.
        size_t force_haplotype_count = 0;
        // If this is set, use_population_mapqs must be set, and we will always
        // try to compute population scores, even if there is nothing to
        // disambiguate. This lets us get an accurate count of scorable reads.
        bool always_check_population = false;
        size_t population_max_paths = 10;
        size_t population_paths_hard_cap = 1000;
        bool top_tracebacks = false;
        // Note that, like the haplotype scoring code, we work with recombiantion penalties in exponent form.
        double recombination_penalty = 20.7; // 20.7 = 9 * 2.3
        size_t rescue_only_min = 128;
        size_t rescue_only_anchor_max = 16;
        size_t order_length_repeat_hit_max = 0;
        int32_t secondary_rescue_subopt_diff = 10;
        size_t min_median_mem_coverage_for_split = 0;
        bool suppress_cluster_merging = false;
        bool suppress_multicomponent_splitting = false;
        size_t alt_anchor_max_length_diff = 5;
        bool dynamic_max_alt_alns = false;
        bool simplify_topologies = false;
        bool use_tvs_clusterer = false;
        bool use_min_dist_clusterer = false;
        bool greedy_min_dist = false;
        bool component_min_dist = false;
        bool no_clustering = false;
        // length of reversing walks during graph extraction
        size_t reversing_walk_length = 0;
        bool suppress_p_value_memoization = false;
        size_t fragment_length_warning_factor = 0;
        size_t max_alignment_gap = 5000;
        bool suppress_mismapping_detection = false;
        bool do_spliced_alignment = false;
        int64_t max_softclip_overlap = 8;
        int64_t max_splice_overhang = 3;
        // about 250k
        int64_t max_intron_length = 1 << 18;
        int64_t min_splice_ref_search_length = 6;
        int64_t max_splice_ref_search_length = 32;
        
        //static size_t PRUNE_COUNTER;
        //static size_t SUBGRAPH_TOTAL;
        //static size_t SECONDARY_RESCUE_COUNT;
        //static size_t SECONDARY_RESCUE_ATTEMPT;
        //static size_t SECONDARY_RESCUE_TOTAL;
        
        /// We often pass around clusters of MEMs and their graph positions, paired with a multiplicity
        using memcluster_t = pair<vector<pair<const MaximalExactMatch*, pos_t>>, double>;
        
        /// This represents a graph for a cluster, and holds a pointer to the
        /// actual extracted graph, a list of assigned MEMs, and the number of
        /// bases of read coverage that that MEM cluster provides (which serves
        /// as a priority).
        using clustergraph_t = tuple<unique_ptr<bdsg::HashGraph>, memcluster_t, size_t>;
        
        /// Represents the mismatches that were allowed in "MEMs" from the fanout
        /// match algorithm
        using match_fanouts_t = unordered_map<const MaximalExactMatch*, deque<pair<string::const_iterator, char>>>;
        
        /// Unique identifier for an unaligned splicing candidate. Specified by:
        /// - Cluster candidate: (is read 1, cluster index, nullptr, pos_t())
        /// - Hit candidate: (is read 1, -1, MEM, position)
        using candidate_id_t = tuple<bool, int64_t, const MaximalExactMatch*, pos_t>;
        
    protected:
        
        /// Enum for the strand of a splice alignment's splice motifs
        enum SpliceStrand {Undetermined, Forward, Reverse};
        
        /// Wrapped internal function that allows some code paths to circumvent the current
        /// mapping quality method option.
        void multipath_map_internal(const Alignment& alignment,
                                    MappingQualityMethod mapq_method,
                                    vector<multipath_alignment_t>& multipath_alns_out);
        
        /// Before the fragment length distribution has been estimated, look for an unambiguous mapping of
        /// the reads using the single ended routine. If we find one record the fragment length and report
        /// the pair, if we don't find one, add the read pair to a buffer instead of the output vector.
        /// Returns true if we successfully find a measurable pair.
        bool attempt_unpaired_multipath_map_of_pair(const Alignment& alignment1, const Alignment& alignment2,
                                                    vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs_out,
                                                    vector<pair<Alignment, Alignment>>& ambiguous_pair_buffer);
        
        /// Extracts a section of graph at a distance from the multipath_alignment_t based on the fragment length
        /// distribution and attempts to align the other paired read to it. If rescuing forward, assumes the
        /// provided multipath_alignment_t is the first read and vice versa if rescuing backward. Rescue constructs
        /// a conventional local alignment with gssw and converts the Alignment to a multipath_alignment_t. The
        /// multipath_alignment_t will be stored in the object passed by reference as an argument.
        bool attempt_rescue(const multipath_alignment_t& multipath_aln, const Alignment& other_aln,
                            bool rescue_forward, multipath_alignment_t& rescue_multipath_aln);
        
        /// Use the algorithm implied by the mapper settings to extract a subgraph to perform a rescue alignment against
        void extract_rescue_graph(const multipath_alignment_t& multipath_aln, const Alignment& other_aln,
                                  bool rescue_forward, MutableHandleGraph* rescue_graph) const;
        
        /// After clustering MEMs, extracting graphs, and assigning hits to cluster graphs, perform
        /// multipath alignment.
        /// Produces topologically sorted multipath_alignment_ts.
        void align_to_cluster_graphs(const Alignment& alignment,
                                     MappingQualityMethod mapq_method,
                                     vector<clustergraph_t>& cluster_graphs,
                                     vector<multipath_alignment_t>& multipath_alns_out,
                                     vector<double>& multiplicities_out,
                                     size_t num_mapping_attempts,
                                     const match_fanouts_t* fanouts = nullptr,
                                     vector<size_t>* cluster_idxs = nullptr);
        
        /// After clustering MEMs, extracting graphs, assigning hits to cluster graphs, and determining
        /// which cluster graph pairs meet the fragment length distance constraints, perform multipath
        /// alignment
        /// Produces topologically sorted multipath_alignment_ts.
        void align_to_cluster_graph_pairs(const Alignment& alignment1, const Alignment& alignment2,
                                          vector<clustergraph_t>& cluster_graphs1,
                                          vector<clustergraph_t>& cluster_graphs2,
                                          vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs_out,
                                          vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs,
                                          vector<double>& pair_multiplicities,
                                          vector<pair<size_t, size_t>>& duplicate_pairs_out,
                                          const match_fanouts_t* fanouts1, const match_fanouts_t* fanouts2);
        
        /// Align the read ends independently, but also try to form rescue alignments for each from
        /// the other. Return true if output obeys pair consistency and false otherwise.
        /// Produces topologically sorted multipath_alignment_ts.
        bool align_to_cluster_graphs_with_rescue(const Alignment& alignment1, const Alignment& alignment2,
                                                 vector<clustergraph_t>& cluster_graphs1,
                                                 vector<clustergraph_t>& cluster_graphs2,
                                                 vector<MaximalExactMatch>& mems1,
                                                 vector<MaximalExactMatch>& mems2,
                                                 vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs_out,
                                                 vector<pair<pair<size_t, size_t>, int64_t>>& pair_distances_out,
                                                 vector<double>& pair_multiplicities_out,
                                                 const match_fanouts_t* fanouts1, const match_fanouts_t* fanouts2);
        
        /// Use the rescue routine on strong suboptimal clusters to see if we can find a good secondary.
        /// Produces topologically sorted multipath_alignment_ts.
        void attempt_rescue_for_secondaries(const Alignment& alignment1, const Alignment& alignment2,
                                            vector<clustergraph_t>& cluster_graphs1,
                                            vector<clustergraph_t>& cluster_graphs2,
                                            vector<pair<size_t, size_t>>& duplicate_pairs,
                                            vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs_out,
                                            vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs,
                                            vector<double>& pair_multiplicities,
                                            const match_fanouts_t* fanouts1, const match_fanouts_t* fanouts2);
        
        /// Merge the rescued mappings into the output vector and deduplicate pairs
        void merge_rescued_mappings(vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs_out,
                                    vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs,
                                    vector<double>& pair_multiplicities,
                                    vector<pair<multipath_alignment_t, multipath_alignment_t>>& rescued_multipath_aln_pairs,
                                    vector<pair<pair<size_t, size_t>, int64_t>>& rescued_cluster_pairs,
                                    vector<double>& rescued_multiplicities) const;
        
        /// Use the oriented distance clusterer or the TVS clusterer to cluster MEMs depending on parameters.
        /// If using oriented distance cluster, must alo provide an oriented distance measurer.
        vector<memcluster_t> get_clusters(const Alignment& alignment, const vector<MaximalExactMatch>& mems,
                                          OrientedDistanceMeasurer* distance_measurer = nullptr,
                                          const match_fanouts_t* fanouts = nullptr) const;
        
        /// Use the oriented distance clusterer or the TVS clusterer to cluster pairs of clusters. Assumes that
        /// the fragment length distribution has been estimated and fixed.
        vector<pair<pair<size_t, size_t>, int64_t>> get_cluster_pairs(const Alignment& alignment1,
                                                                      const Alignment& alignment2,
                                                                      vector<clustergraph_t>& cluster_graphs1,
                                                                      vector<clustergraph_t>& cluster_graphs2,
                                                                      OrientedDistanceMeasurer* distance_measurer = nullptr);
        
        /// Extracts a subgraph around each cluster of MEMs that encompasses any
        /// graph position reachable (according to the Mapper's aligner) with
        /// local alignment anchored at the MEMs. If any subgraphs overlap, they
        /// are merged into one subgraph. Returns a vector of all the merged
        /// cluster subgraphs, their MEMs assigned from the mems vector
        /// according to the MEMs' hits, and their read coverages in bp. The
        /// caller must delete the VG objects produced!
        vector<clustergraph_t> query_cluster_graphs(const Alignment& alignment,
                                                    const vector<MaximalExactMatch>& mems,
                                                    const vector<memcluster_t>& clusters) const;
        
        /// Return a graph (on the heap) that contains a cluster. The paired bool
        /// indicates whether the graph is known to be connected (but it is possible
        /// for the graph to be connected and have it return false)
        pair<unique_ptr<bdsg::HashGraph>, bool> extract_cluster_graph(const Alignment& alignment,
                                                                      const memcluster_t& mem_cluster) const;
        
        /// Extract a graph that is guaranteed to contain all local alignments that include
        /// the MEMs of the cluster.  The paired bool indicates whether the graph is
        /// known to be connected (but it is possible for the graph to be connected and have
        /// it return false)
        pair<unique_ptr<bdsg::HashGraph>, bool> extract_maximal_graph(const Alignment& alignment,
                                                                      const memcluster_t& mem_cluster) const;
        
        /// Extract a graph with an algorithm that tries to extract not much more than what
        /// is required to contain the cluster in a single connected component (can be slower
        /// than the maximal algorithm for alignments that require large indels),  The paired bool
        /// indicates whether the graph is known to be connected (but it is possible
        /// for the graph to be connected and have it return false)
        pair<unique_ptr<bdsg::HashGraph>, bool> extract_restrained_graph(const Alignment& alignment,
                                                                         const memcluster_t& mem_cluster) const;
        
        /// Returns the union of the intervals on the read that a cluster cover in sorted order
        vector<pair<int64_t, int64_t>> covered_intervals(const Alignment& alignment,
                                                         const clustergraph_t& cluster) const;
        
        /// If there are any multipath_alignment_ts with multiple connected components, split them
        /// up and add them to the return vector.
        /// Properly handles multipath_alignment_ts that are unmapped.
        /// Does not depend on or guarantee topological order in the multipath_alignment_ts.
        void split_multicomponent_alignments(vector<multipath_alignment_t>& multipath_alns_out,
                                             const Alignment* alignment = nullptr,
                                             vector<clustergraph_t>* cluster_graphs = nullptr,
                                             vector<size_t>* cluster_idxs = nullptr,
                                             vector<double>* multiplicities = nullptr) const;
        
        /// If there are any multipath_alignment_ts with multiple connected components, split them
        /// up and add them to the return vector, also measure the distance between them and add
        /// a record to the cluster pairs vector.
        /// Properly handles multipath_alignment_ts that are unmapped.
        /// Does not depend on or guarantee topological order in the multipath_alignment_ts.
        void split_multicomponent_alignments(const Alignment& alignment1, const Alignment& alignment2,
                                             vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs_out,
                                             vector<clustergraph_t>& cluster_graphs1,
                                             vector<clustergraph_t>& cluster_graphs2,
                                             vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs,
                                             vector<double>& multiplicities) const;
        
        /// Helper function to be called by split_multicomponent_alignments to reassign hits to the
        /// split clusters
        void reassign_split_clusters(const Alignment& alignment,
                                     vector<clustergraph_t>& cluster_graphs,
                                     const vector<const multipath_alignment_t*>& split_mp_alns,
                                     const vector<size_t*>& cluster_assignments,
                                     const vector<size_t*>& all_cluster_assignments) const;
        
        /// Combine all of the significant alignments into one. Requires alignments to be sorted by
        /// significance already
        void agglomerate_alignments(vector<multipath_alignment_t>& multipath_alns_out,
                                    vector<double>* multiplicities = nullptr) const;
        
        /// Combine all of the significant alignments into one pair. Requires alignments to be sorted by
        /// significance already
        void agglomerate_alignment_pairs(vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs_out,
                                         vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs,
                                         vector<double>& multiplicities) const;
        
        /// The internal agglomeration procedure
        void agglomerate(size_t idx, multipath_alignment_t& agglomerating, const multipath_alignment_t& multipath_aln,
                         vector<size_t>& agglomerated_group, unordered_set<pos_t>& agg_start_positions,
                         unordered_set<pos_t>& agg_end_positions) const;
        
        /// Look for spliced alignments among the results of various stages in the mapping algorithm
        /// Returns true if any spliced alignments were made
        bool find_spliced_alignments(const Alignment& alignment, vector<multipath_alignment_t>& multipath_alns_out,
                                     vector<double>& multiplicities, vector<size_t>& cluster_idxs,
                                     const vector<MaximalExactMatch>& mems, vector<clustergraph_t>& cluster_graphs,
                                     const match_fanouts_t* fanouts = nullptr);
        
        /// Look for spliced alignments among the results of various stages in the mapping algorithm for pairs
        /// Returns true if any spliced alignments were made
        bool find_spliced_alignments(const Alignment& alignment1, const Alignment& alignment2,
                                     vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs_out,
                                     vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs,
                                     vector<double>& pair_multiplicities,
                                     const vector<MaximalExactMatch>& mems1, const vector<MaximalExactMatch>& mems2,
                                     vector<clustergraph_t>& cluster_graphs1, vector<clustergraph_t>& cluster_graphs2,
                                     const match_fanouts_t* fanouts = nullptr);
        
        /// Find candidates for spliced alignment sections for a given multipath alignment among the
        /// aligned clusters
        void identify_aligned_splice_candidates(const Alignment& alignment, bool search_left,
                                                const pair<int64_t, int64_t>& primary_interval,
                                                const vector<multipath_alignment_t>& multipath_alns,
                                                const vector<size_t>& cluster_idxs,
                                                const vector<int64_t>& current_index, int64_t anchor,
                                                unordered_set<size_t>& clusters_used_out,
                                                vector<size_t>& mp_aln_candidates_out) const;

        /// Find candidates for spliced alignment sections for a given multipath alignment among the
        /// aligned cluster pairs
        void identify_aligned_splice_candidates(const Alignment& alignment, bool read_1, bool search_left,
                                                const pair<int64_t, int64_t>& primary_interval,
                                                const vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs,
                                                const vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs,
                                                const vector<int64_t>& current_index, int64_t anchor,
                                                unordered_set<size_t>& clusters_used_out,
                                                vector<size_t>& mp_aln_candidates_out) const;
        
        /// Find candidates for spliced alignment sections for a given multipath alignment among the
        /// unaligned clusters and MEMs
        void identify_unaligned_splice_candidates(const Alignment& alignment, bool search_left,
                                                  const pair<int64_t, int64_t>& primary_interval,
                                                  const vector<MaximalExactMatch>& mems,
                                                  const vector<clustergraph_t>& cluster_graphs,
                                                  const unordered_set<size_t>& clusters_already_used,
                                                  vector<size_t>& cluster_candidates_out,
                                                  vector<pair<const MaximalExactMatch*, pos_t>>& hit_candidates_out) const;
        
        /// Make alignments for the splice alignment cancidates from MEMs and unaligned clusters
        void align_to_splice_candidates(const Alignment& alignment,
                                        vector<clustergraph_t>& cluster_graphs,
                                        const vector<MaximalExactMatch>& mems,
                                        const vector<size_t>& cluster_candidates,
                                        const vector<pair<const MaximalExactMatch*, pos_t>>& hit_candidates,
                                        const pair<int64_t, int64_t>& primary_interval,
                                        bool searching_left,
                                        bool is_read_1,
                                        unordered_map<candidate_id_t, pair<multipath_alignment_t, double>>& unaligned_candidate_bank,
                                        vector<candidate_id_t>& candidates_out,
                                        const match_fanouts_t* mem_fanouts = nullptr) const;
        
        /// Check whether splice segment candidates can form a statistically significant spliced
        /// alignment. Returns true if a spliced alignment is made
        bool test_splice_candidates(const Alignment& alignment, bool searching_left,
                                    multipath_alignment_t& anchor_mp_aln, double& anchor_multiplicity,
                                    SpliceStrand& strand, int64_t num_candidates,
                                    const function<const multipath_alignment_t&(int64_t)>& get_candidate,
                                    const function<double(int64_t)>& get_multiplicity,
                                    const function<multipath_alignment_t&&(int64_t)>& consume_candidate);
        
        /// Check if any of the unpaired spliced alignments can make pairs now
        /// If any pairs are identified, can invalidate the input alignments
        bool retry_pairing_spliced_alignments(const Alignment& alignment1, const Alignment& alignment2,
                                              vector<multipath_alignment_t>& multipath_alns_1,
                                              vector<multipath_alignment_t>& multipath_alns_2,
                                              const vector<size_t>& cluster_idxs_1,
                                              const vector<size_t>& cluster_idxs_2,
                                              const vector<double>& multiplicities_1,
                                              const vector<double>& multiplicities_2,
                                              vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs_out,
                                              vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs_out,
                                              vector<double>& pair_multiplicities_out) const;
        

        /// Make a multipath alignment of the read against the indicated graph and add it to
        /// the list of multimappings.
        /// Does NOT necessarily produce a multipath_alignment_t in topological order.
        void multipath_align(const Alignment& alignment,
                             clustergraph_t& cluster_graph,
                             multipath_alignment_t& multipath_aln_out,
                             const match_fanouts_t* fanouts) const;
        
        /// If any softclips could have arisen because not enough graph was extracted, extract
        /// extra graph in those areas. Returns true if the graph was expanded.
        bool expand_for_softclips(clustergraph_t& cluster_graph,
                                  const multipath_alignment_t& multipath_aln) const;
        
        /// Removes the sections of an Alignment's path within snarls and re-aligns them with multiple traceback
        /// to create a multipath alignment with non-trivial topology.
        /// Guarantees that the resulting multipath_alignment_t is in topological order.
        void make_nontrivial_multipath_alignment(const Alignment& alignment, const HandleGraph& subgraph,
                                                 const function<pair<id_t, bool>(id_t)>& translator,
                                                 multipath_alignment_t& multipath_aln_out) const;
        
        /// Remove the full length bonus from all source or sink subpaths that received it
        void strip_full_length_bonuses(multipath_alignment_t& multipath_aln) const;
        
        /// Returns a vector of log-likelihoods for each mapping
        vector<double> mapping_likelihoods(vector<multipath_alignment_t>& multipath_alns) const;
        
        /// Returns a vector of log-likelihoods for each pair mapping
        vector<double> pair_mapping_likelihoods(vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs,
                                                const vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs) const;
        
        /// Compute a mapping quality from a list of scores, using the selected method.
        /// Optionally considers non-present duplicates of the scores encoded as multiplicities
        int32_t compute_raw_mapping_quality_from_scores(const vector<double>& scores, MappingQualityMethod mapq_method,
                                                        bool have_qualities, const vector<double>* multiplicities = nullptr) const;
        
        
        /// Sorts mappings by score and store mapping quality of the optimal alignment in the multipath_alignment_t object
        /// Optionally also sorts a vector of indexes to keep track of the cluster-of-origin
        /// Allows multipath alignments where the best single path alignment is leaving the read unmapped.
        /// multipath_alignment_ts MUST be topologically sorted.
        void sort_and_compute_mapping_quality(vector<multipath_alignment_t>& multipath_alns, MappingQualityMethod mapq_method,
                                              vector<size_t>* cluster_idxs = nullptr, vector<double>* multiplicities = nullptr) const;
        
        /// Sorts mappings by score and store mapping quality of the optimal alignment in the multipath_alignment_t object
        /// If there are ties between scores, breaks them by the expected distance between pairs as computed by the
        /// OrientedDistanceClusterer::cluster_pairs function (modified cluster_pairs vector)
        /// Allows multipath alignments where the best single path alignment is leaving the read unmapped.
        /// multipath_alignment_ts MUST be topologically sorted.
        /// Optionally considers non-present duplicates of the scores encoded as multiplicities
        void sort_and_compute_mapping_quality(vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs,
                                              vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs,
                                              vector<pair<size_t, size_t>>* duplicate_pairs_out = nullptr,
                                              vector<double>* pair_multiplicities = nullptr) const;
        
        /// Estimates the number of equivalent mappings (including this one), which we may not have seen due to
        /// unexplored rescues.
        double estimate_missed_rescue_multiplicity(size_t which_pair,
                                                   const vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs,
                                                   const vector<clustergraph_t>& cluster_graphs1,
                                                   const vector<clustergraph_t>& cluster_graphs2,
                                                   bool from_secondary_rescue) const;
        
        /// Estimates the number of equivalent mappings (including this one), which we may not have seen due to
        /// limits on the numbers of hits returns for a MEM
        double cluster_multiplicity(const memcluster_t& cluster) const;
        
        /// Estimates the number of equivalent pair mappings (including this one), which we may not have seen due to
        /// limits on the numbers of hits returns for a MEM
        double pair_cluster_multiplicity(const memcluster_t& cluster_1, const memcluster_t& cluster_2) const;
        
        /// Computes the log-likelihood of a given fragment length in the trained distribution
        double fragment_length_log_likelihood(int64_t length) const;
        
        /// Computes the number of read bases a cluster of MEM hits covers.
        static void set_read_coverage(clustergraph_t& cluster_graph);
        
        /// Would an alignment this good be expected against a graph this big by chance alone
        bool likely_mismapping(const multipath_alignment_t& multipath_aln);
        
        /// Would an alignment this good be expected against a graph this big by chance alone
        bool likely_misrescue(const multipath_alignment_t& multipath_aln);
        
        /// A scaling of a score so that it approximately follows the distribution of the longest match in p-value test
        int64_t pseudo_length(const multipath_alignment_t& multipath_aln) const;
        
        /// The approximate p-value for a match length of the given size against the current graph
        double random_match_p_value(int64_t match_length, size_t read_length);
        
        /// Reorganizes the fan-out breaks into the format that MultipathAlignmentGraph wants it in
        match_fanouts_t record_fanouts(const vector<MaximalExactMatch>& mems,
                                       vector<deque<pair<string::const_iterator, char>>>& fanouts) const;
        
        /// Get a distance measurer based on the configuartion of the mapper
        unique_ptr<OrientedDistanceMeasurer> get_distance_measurer(MemoizingGraph& memoizing_graph) const;
        
        /// Compute the approximate distance between two multipath alignments
        /// If either is unmapped, or the distance cannot be obtained, returns numeric_limits<int64_t>::max()
        int64_t distance_between(const multipath_alignment_t& multipath_aln_1, const multipath_alignment_t& multipath_aln_2,
                                 bool full_fragment = false, bool forward_strand = false) const;
        
        int64_t distance(const pos_t& pos_1, const pos_t& pos_2) const;
        
        /// Are two multipath alignments consistently placed based on the learned fragment length distribution?
        bool are_consistent(const multipath_alignment_t& multipath_aln_1, const multipath_alignment_t& multipath_aln_2) const;
        
        /// Is this a consistent inter-pair distance based on the learned fragment length distribution?
        bool is_consistent(int64_t distance) const;
        
        /// Return true if any of the initial positions of the source Subpaths are shared between the two
        /// multipath alignments
        bool share_terminal_positions(const multipath_alignment_t& multipath_aln_1, const multipath_alignment_t& multipath_aln_2) const;
        
        /// Get a thread_local RRMemo with these parameters
        haploMath::RRMemo& get_rr_memo(double recombination_penalty, size_t population_size) const;
        
        /// Detects if each pair can be assigned to a consistent strand of a path, and if not removes them. Also
        /// inverts the distances in the cluster pairs vector according to the strand
        void establish_strand_consistency(vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs,
                                          vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs);
        
        /// A restrained estimate of the amount of gap we would like to align for a read tail
        int64_t pessimistic_gap(int64_t length, double multiplier) const;

        /// Return exact matches according to the object's parameters
        /// If using the fan-out algorithm, we can optionally leave fan-out MEMs in tact and
        /// return a vector of their breaks.
        vector<MaximalExactMatch> find_mems(const Alignment& alignment,
                                            vector<deque<pair<string::const_iterator, char>>>* mem_fanout_breaks = nullptr);
        
        
        int64_t min_softclip_length_for_splice = 20;
        int64_t min_softclipped_score_for_splice = 25;
        
        DinucleotideMachine dinuc_machine;
        SpliceMotifs splice_motifs;
        SnarlManager* snarl_manager;
        MinimumDistanceIndex* distance_index;
        unique_ptr<PathComponentIndex> path_component_index;
        
        static const size_t RESCUED;
        
        /// Memos used by population model
        static thread_local unordered_map<pair<double, size_t>, haploMath::RRMemo> rr_memos;
        
        // a memo for the transcendental p-value function (thread local to maintain threadsafety)
        static thread_local unordered_map<pair<int64_t, size_t>, double> p_value_memo;
        
        // a memo for the transcendental restrained extraction function (thread local to maintain threadsafety)
        static thread_local unordered_map<double, vector<int64_t>> pessimistic_gap_memo;
        static const size_t gap_memo_max_size;
        
        // a memo for transcendental band padidng function (gets initialized at construction)
        vector<size_t> band_padding_memo;
        
#ifdef mpmap_instrument_mem_statistics
    public:
        ofstream _mem_stats;
        bool _wrote_mem_stats_header = false;
#endif
    };
        
}



#endif /* multipath_mapper_hpp */
