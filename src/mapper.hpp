#ifndef VG_MAPPER_HPP_INCLUDED
#define VG_MAPPER_HPP_INCLUDED

#include <iostream>
#include <map>
#include <chrono>
#include <ctime>
#include "omp.h"
#include "vg.hpp"
#include "bdsg/hash_graph.hpp"
#include <gcsa/gcsa.h>
#include <gcsa/lcp.h>
#include <gbwt/gbwt.h>
#include "position.hpp"
#include "vg/io/json2pb.h"
#include "aligner.hpp"
#include "mem.hpp"
#include "mem_accelerator.hpp"
#include "cluster.hpp"
// TODO: pull out ScoreProvider into its own file
#include "haplotypes.hpp"

// #define BENCH
// #include "bench.h"

namespace vg {

// uncomment to make vg map --debug very interesting
//#define debug_mapper

using namespace std;
    
enum MappingQualityMethod { Approx, Exact, Adaptive, None };

class Mapper;

// for banded long read alignment resolution

class AlignmentChainModelVertex {
public:
    Alignment* aln;
    vector<pair<AlignmentChainModelVertex*, double> > next_cost; // for forward
    vector<pair<AlignmentChainModelVertex*, double> > prev_cost; // for backward
    double weight;
    double score;
    unordered_map<path_handle_t, vector<pair<size_t, bool> > > positions;
    int band_begin;
    int band_idx;
    AlignmentChainModelVertex* prev;
    AlignmentChainModelVertex(void) = default;                                      // Copy constructor
    AlignmentChainModelVertex(const AlignmentChainModelVertex&) = default;               // Copy constructor
    AlignmentChainModelVertex(AlignmentChainModelVertex&&) = default;                    // Move constructor
    AlignmentChainModelVertex& operator=(const AlignmentChainModelVertex&) & = default;  // AlignmentChainModelVertexopy assignment operator
    AlignmentChainModelVertex& operator=(AlignmentChainModelVertex&&) & = default;       // Move assignment operator
    virtual ~AlignmentChainModelVertex() { }                     // Destructor
};

class AlignmentChainModel {
public:
    vector<AlignmentChainModelVertex> model;
    unordered_map<path_handle_t, map<int64_t, vector<vector<AlignmentChainModelVertex>::iterator> > > positions;
    set<vector<AlignmentChainModelVertex>::iterator> redundant_vertexes;
    vector<Alignment> unaligned_bands;
    AlignmentChainModel(
        vector<vector<Alignment> >& bands,
        Mapper* mapper,
        const function<double(const Alignment&, const Alignment&, const unordered_map<path_handle_t, vector<pair<size_t, bool> > >&, const unordered_map<path_handle_t, vector<pair<size_t, bool> > >&, int64_t)>& transition_weight,
        int vertex_band_width = 10,
        int position_depth = 1,
        int max_connections = 30);
    void score(const unordered_set<AlignmentChainModelVertex*>& exclude);
    AlignmentChainModelVertex* max_vertex(void);
    vector<Alignment> traceback(const Alignment& read, int alt_alns, bool paired, bool debug);
    void display(ostream& out);
    void display_dot(ostream& out, vector<AlignmentChainModelVertex*> vertex_trace);
    void clear_scores(void);
};

/*
 * A class that keeps a running estimation of a fragment length distribution
 * using a robust estimation formula in order to be insensitive to outliers.
 */
class FragmentLengthDistribution {
public:
    
    /// Initialize distribution
    ///
    /// Args:
    ///  maximum_sample_size         sample size at which reestimation stops
    ///  reestimation_frequency      update running estimate after this many samples
    ///  robust_estimation_fraction  robustly estimate using this fraction of samples
    FragmentLengthDistribution(size_t maximum_sample_size,
                               size_t reestimation_frequency,
                               double robust_estimation_fraction);
    FragmentLengthDistribution(void);
    ~FragmentLengthDistribution();
    
    
    /// Instead of estimating anything, just use these parameters.
    void force_parameters(double mean, double stddev);
    
    /// Record an observed fragment length
    void register_fragment_length(int64_t length);

    /// Robust mean of the distribution observed so far
    double mean() const;
    
    /// Robust standard deviation of the distribution observed so far
    double std_dev() const;
    
    /// Returns true if the maximum sample size has been reached, which finalizes the
    /// distribution estimate
    bool is_finalized() const;
    
    /// Returns the max sample size up to which the distribution will continue to reestimate
    /// parameters
    size_t max_sample_size() const;
    
    /// Returns the number of samples that have been collected so far
    size_t curr_sample_size() const;
    
    /// Begin iterator to the measurements that the distribution has used to estimate the
    /// parameters
    multiset<double>::const_iterator measurements_begin() const;
    
    /// End iterator to the measurements that the distribution has used to estimate the
    /// parameters
    multiset<double>::const_iterator measurements_end() const;
    
private:
    multiset<double> lengths;
    bool is_fixed = false;
    
    double robust_estimation_fraction;
    size_t maximum_sample_size;
    size_t reestimation_frequency;
    
    double mu = 0.0;
    double sigma = 1.0;
    
    void estimate_distribution();
};

/**
 * A class trait/mixin which defines a mapper's paired end distribution support.
 *
 * Doesn't actually define the paired-end mapping methods.
 */
class PairedEndMapper {
public:

    /// Set parameters for estimating fragment length distribution.
    /// TODO: setting alignment threads after this could mess up the internal memory for how many threads to reset to
    void set_fragment_length_distr_params(size_t maximum_sample_size = 1000, size_t reestimation_frequency = 1000,
                                          double robust_estimation_fraction = 0.95);

    /// Returns true if fragment length distribution has been fixed
    bool has_fixed_fragment_length_distr();
    
    /// Use the given fragment length distribution parameters instead of
    /// estimating them.
    void force_fragment_length_distr(double mean, double stddev);

protected:
    /// Holds the actual fragment length distribution and estimation information
    FragmentLengthDistribution fragment_length_distr;

};

/**
 * Base class for basic mapping functionality shared between the Mapper, MultipathMapper, etc.
 * Handles holding on to the random access and text indexes needed for mapping operations.
 */
class BaseMapper : public AlignerClient, public PairedEndMapper {
    
public:
    /**
     * Make a BaseMapper that pulls from an PathPositionHandleGraph succinct graph and a GCSA2 kmer
     * index + LCP array, and which can score reads against haplotypes using
     * the given ScoreProvider.
     *
     * If the GCSA and LCPArray are null, cannot do search, only alignment.
     */
    BaseMapper(PathPositionHandleGraph* xidex, gcsa::GCSA* g, gcsa::LCPArray* a, haplo::ScoreProvider* haplo_score_provider = nullptr);
    BaseMapper(void);
    
    /// We need to be able to estimate the GC content from the GCSA index in the constructor.
    /// The given index may be null.
    static double estimate_gc_content(const gcsa::GCSA* gcsa);
    
    int random_match_length(double chance_random);
                         
    
    /// Override alignment score setting to support haplotype consistency exponent
    void set_alignment_scores(int8_t match, int8_t mismatch, int8_t gap_open, int8_t gap_extend, int8_t full_length_bonus,
                              double haplotype_consistency_exponent = 1);
    
    /// Same, but loading a 4x4 substitution score matrix from a stream
    void set_alignment_scores(istream& matrix_stream, int8_t gap_open, int8_t gap_extend, int8_t full_length_bonus,
                              double haplotype_consistency_exponent = 1);
        
    // MEM-based mapping
    // find maximal exact matches
    // These are SMEMs by definition when shorter than the max_mem_length or GCSA2 order.
    // Designating reseed_length returns minimally-more-frequent sub-MEMs in addition to SMEMs when SMEM is >= reseed_length.
    // Minimally-more-frequent sub-MEMs are MEMs contained in an SMEM that have occurrences outside of the SMEM.
    // SMEMs and sub-MEMs will be automatically filled with the nodes they contain
    
    vector<MaximalExactMatch>
    find_mems_deep(string::const_iterator seq_begin,
                   string::const_iterator seq_end,
                   double& lcp_avg,
                   double& fraction_filtered,
                   int max_mem_length = 0,
                   int min_mem_length = 1,
                   int reseed_length = 0,
                   bool use_lcp_reseed_heuristic = false,
                   bool use_diff_based_fast_reseed = false,
                   bool include_parent_in_sub_mem_count = false,
                   bool record_max_lcp = false,
                   int reseed_below_count = 0);
    
    // Use the GCSA2 index to find super-maximal exact matches.
    vector<MaximalExactMatch>
    find_mems_simple(string::const_iterator seq_begin,
                     string::const_iterator seq_end,
                     int max_mem_length = 0,
                     int min_mem_length = 1,
                     int reseed_length = 0);
    
    vector<MaximalExactMatch>
    find_stripped_matches(string::const_iterator seq_begin,
                          string::const_iterator seq_end,
                          size_t strip_length,
                          size_t max_match_length,
                          size_t target_count);
    
    // finds MEMs where a pre-specified number of low-quality bases are
    // allowed to be any base. if the optional vector is provided, then it
    // will be filled to include all of the places that each returned MEM
    // mismatches the graph sequence. otherwise, the MEMs are walked out
    // and split into exact matches (can be expensive)
    vector<MaximalExactMatch>
    find_fanout_mems(string::const_iterator seq_begin,
                     string::const_iterator seq_end,
                     string::const_iterator qual_begin,
                     int max_fans_out,
                     char max_fanout_base_quality,
                     vector<deque<pair<string::const_iterator, char>>>* mem_fanout_breaks = nullptr);
    
    vector<pos_t> walk_fanout_path(string::const_iterator begin,
                                   string::const_iterator end,
                                   const deque<pair<string::const_iterator, char>>& fanout_breaks,
                                   gcsa::node_type pos);
    
    /// identifies tracts of order-length MEMs that were unfilled because their hit count was above the max
    /// and fills one MEM in the tract (the one with the smallest hit count), assumes MEMs are lexicographically
    /// ordered by read index
    void rescue_high_count_order_length_mems(vector<MaximalExactMatch>& mems,
                                             size_t max_rescue_hit_count);
    
    /// identifies hits for order-length MEMs that are actually part of longer MEMs above the GCSA's limit and
    /// merges them. for speed's sake, can have false negatives but no false positives
    void precollapse_order_length_runs(string::const_iterator seq_begin,
                                       vector<MaximalExactMatch>& mems);
    
    /// identifies hits for sub-MEMs that are redundant hits to the parent MEMs and removes them
    /// from the hit lists. for speed's sake, can have false negatives but no false positives
    void prefilter_redundant_sub_mems(vector<MaximalExactMatch>& mems,
                                      vector<pair<int, vector<size_t>>>& sub_mem_containment_graph);
    
    int sub_mem_thinning_burn_in = 16; // start counting at this many bases to verify sub-MEM count
    int sub_mem_count_thinning = 4; // count every this many bases to verify sub-MEM count
    int min_mem_length; // a mem must be >= this length
    int mem_reseed_length; // the length above which we reseed MEMs to get potentially missed hits
    bool fast_reseed = true; // use the fast reseed algorithm
    double fast_reseed_length_diff = 0.45; // how much smaller than its parent a sub-MEM can be in the fast reseed algorithm
    bool adaptive_reseed_diff = true; // use an adaptive length difference algorithm in reseed algorithm
    double adaptive_diff_exponent = 0.065; // exponent that describes limiting behavior of adaptive diff algorithm
    int hit_max = 0;       // only query at most this many hits for a MEM (0 for no limit)
    int hard_hit_max = 0; // don't query any hits for MEMs with this many occurrences or more (0 for no limit)
    bool use_approx_sub_mem_count = false;
    bool prefilter_redundant_hits = true;
    int max_sub_mem_recursion_depth = 2;
    bool use_greedy_mem_restarts = false;
    int greedy_restart_min_length = 40;
    int greedy_restart_max_count = 2;
    int greedy_restart_max_lcp = 0; // 0 for no max
    bool greedy_restart_assume_substitution = false;
    bool filter_short_mems = false;
    double short_mem_filter_factor = 0.45;
    int unpaired_penalty = 17;
    bool precollapse_order_length_hits = true;
    double avg_node_length = 0;
    size_t total_seq_length = 0;
    int fanout_length_threshold = 0;
    
    // The recombination rate (negative log per-base recombination probability) for haplotype-aware mapping
    double recombination_penalty = 20.7; // 9 * 2.3 = 20.7
    
    // Remove any bonuses used by the aligners from the final reported scores.
    // Does NOT (yet) remove the haplotype consistency bonus.
    bool strip_bonuses; 
    bool assume_acyclic; // the indexed graph is acyclic
    
    
    MappingQualityMethod mapping_quality_method; // how to compute mapping qualities
    int max_mapping_quality; // the cap for mapping quality
    
    bool exclude_unaligned = false;

    /// Set to enable debugging messages to cerr from the mapper, so a user can understand why a read maps the way it does.
    bool debug = false;
    
    /// Locate the sub-MEMs contained in the last MEM of the mems vector that have ending positions
    /// before the end the next SMEM, label each of the sub-MEMs with the indices of all of the SMEMs
    /// that contain it
    void find_sub_mems(const vector<MaximalExactMatch>& mems,
                       int parent_layer_begin,
                       int parent_layer_end,
                       int mem_idx,
                       string::const_iterator next_mem_end,
                       int min_mem_length,
                       vector<pair<MaximalExactMatch, vector<size_t>>>& sub_mems_out);
    
    /// Provides same semantics as find_sub_mems but with a different algorithm. This algorithm uses the
    /// min_mem_length as a pruning tool instead of the LCP index. It can be expected to be faster when both
    /// the min_mem_length reasonably large relative to the reseed_length (e.g. 1/2 of SMEM size or similar).
    void find_sub_mems_fast(const vector<MaximalExactMatch>& mems,
                            int parent_layer_begin,
                            int parent_layer_end,
                            int mem_idx,
                            string::const_iterator leftmost_guaranteed_disjoint_bound,
                            string::const_iterator leftmost_seeding_bound,
                            int min_sub_mem_length,
                            vector<pair<MaximalExactMatch, vector<size_t>>>& sub_mems_out);
    
    /// If possible, use the MEMAcclerator to get the initial range for a MEM and update the cursor
    /// accordingly. If this is not possible, return the full GCSA2 range and leave the cursor unaltered.
    gcsa::range_type accelerate_mem_query(string::const_iterator begin,
                                          string::const_iterator& cursor) const;
    
    // Use the GCSA index to look up the sequence
    set<pos_t> sequence_positions(const string& seq);
    
    // Algorithm for choosing an adaptive reseed length based on the length of the parent MEM
    size_t get_adaptive_min_reseed_length(size_t parent_mem_length);
    
    /// Score all of the alignments in the vector for haplotype consistency. If
    /// all of them can be scored (i.e. none of them visit nodes/edges with no
    /// haplotypes), adjust all of their scores to reflect haplotype
    /// consistency. If one or more cannot be scored for haplotype consistency,
    /// leave the alignment scores alone.
    void apply_haplotype_consistency_scores(const vector<Alignment*>& alns);
    
    // thread_local to allow alternating reads/writes
    thread_local static vector<size_t> adaptive_reseed_length_memo;
    
    // xg index
    PathPositionHandleGraph* xindex = nullptr;
    
    // GCSA index and its LCP array
    gcsa::GCSA* gcsa = nullptr;
    gcsa::LCPArray* lcp = nullptr;
    MEMAccelerator* accelerator = nullptr;
    
    // Haplotype score provider, if any, for determining haplotype concordance
    haplo::ScoreProvider* haplo_score_provider = nullptr;
    
    // The exponent for the haplotype consistency score.
    // 0 = no haplotype consistency scoring done.
    // 1 = multiply in haplotype likelihood once when computing alignment score
    double haplotype_consistency_exponent = 1;
};

/**
 * Keeps track of statistics about fragment length within the Mapper class.
 * Belongs to a single thread.
 */
class FragmentLengthStatistics {
public:

    void record_fragment_configuration(const Alignment& aln1, const Alignment& aln2, Mapper* mapper);

    string fragment_model_str(void);
    void save_frag_lens_to_alns(Alignment& aln1, Alignment& aln2,
                                const unordered_map<path_handle_t, int64_t>& approx_frag_lengths,
                                PathPositionHandleGraph* xindex,
                                bool is_consistent);
    
    // These functions are the authorities on the estimated parameters
    double fragment_length_stdev(void);
    double fragment_length_mean(void);
    double fragment_length_pdf(double length);
    double fragment_length_pval(double length);
    bool fragment_orientation(void);
    bool fragment_direction(void);
    
    // These cached versions of the parameters are updated periodically
    double cached_fragment_length_mean = 0;
    double cached_fragment_length_stdev = 0;
    bool cached_fragment_orientation_same = 0;
    bool cached_fragment_direction = 1;
    
    // These variables are used to manage the periodic updates
    int64_t since_last_fragment_length_estimate = 0;
    int64_t fragment_model_update_interval = 100;
    
    // These deques are used for the periodic running estimation of the fragment length distribution
    deque<double> fragment_lengths;
    deque<bool> fragment_orientations;
    deque<bool> fragment_directions;

    int64_t fragment_max = 10000; // the maximum length fragment which we will consider when estimating fragment lengths
    int64_t fragment_size = 0; // Used to bound clustering of MEMs during paired end mapping, also acts as sentinel to determine
                       // if consistent pairs should be reported; dynamically estimated at runtime
    double fragment_sigma = 10; // the number of times the standard deviation above the mean to set the fragment_size
    int64_t fragment_length_cache_size = 10000;
    float perfect_pair_identity_threshold = 0.9;
    bool fixed_fragment_model = true;
    
    
    
    
    
};

class Mapper : public BaseMapper {


private:
    
    vector<Alignment> align_multi_internal(bool compute_unpaired_qualities,
                                           const Alignment& aln,
                                           int kmer_size,
                                           int stride,
                                           int max_mem_length,
                                           int band_width,
                                           int band_overlap,
                                           double& cluster_mq,
                                           int keep_multimaps = 0,
                                           int additional_multimaps = 0,
                                           vector<MaximalExactMatch>* restricted_mems = nullptr,
                                           bool xdrop_alignment = false);
    void compute_mapping_qualities(vector<Alignment>& alns, double cluster_mq, double mq_estimate, double mq_cap);
    void compute_mapping_qualities(pair<vector<Alignment>, vector<Alignment>>& pair_alns, double cluster_mq, double mq_estmate1, double mq_estimate2, double mq_cap1, double mq_cap2);
    vector<Alignment> score_sort_and_deduplicate_alignments(vector<Alignment>& all_alns, const Alignment& original_alignment);
    void filter_and_process_multimaps(vector<Alignment>& all_alns, int total_multimaps);
    // Return the one best banded alignment.
    vector<Alignment> align_banded(const Alignment& read,
                                   int kmer_size = 0,
                                   int stride = 0,
                                   int max_mem_length = 0,
                                   int band_width = 1000,
                                   int band_overlap = 500,
                                   bool xdrop_alignment = false);
    // alignment based on the MEM approach
//    vector<Alignment> align_mem_multi(const Alignment& alignment, vector<MaximalExactMatch>& mems, double& cluster_mq, double lcp_avg, int max_mem_length, int additional_multimaps = 0);
    // uses approximate-positional clustering based on embedded paths in the xg index to find and align against alignment targets
    vector<Alignment> align_mem_multi(const Alignment& aln,
                                      vector<MaximalExactMatch>& mems,
                                      double& cluster_mq,
                                      double lcp_avg,
                                      double fraction_filtered,
                                      int max_mem_length,
                                      int keep_multimaps,
                                      int additional_multimaps,
                                      bool xdrop_alignment = false);
    
protected:
    Alignment align_to_graph(const Alignment& aln,
                             HandleGraph& graph,
                             bool do_flip,
                             bool traceback,
                             bool pinned_alignment = false,
                             bool pin_left = false,
                             bool banded_global = false,
                             bool keep_bonuses = true);
    Alignment align_to_graph(const Alignment& aln,
                             HandleGraph& graph,
                             const vector<MaximalExactMatch>& mems,
                             bool do_flip,
                             bool traceback,
                             bool pinned_alignment = false,
                             bool pin_left = false,
                             bool banded_global = false,
                             int xdrop_alignment = 0,           // 0 to disable X-drop, otherwise enabled (1 for forward, 2 for reverse-complemented; keep compatble to align_maybe_flip)
                             bool keep_bonuses = true);
    
    // make the bands used in banded alignment
    vector<Alignment> make_bands(const Alignment& read, int band_width, int band_overlap, vector<pair<int, int>>& to_strip);
public:
    // Make a Mapper that pulls from an PathPositionHandleGraph succinct graph, a GCSA2 kmer index +
    // LCP array, and an optional haplotype score provider.
    Mapper(PathPositionHandleGraph* xidex, gcsa::GCSA* g, gcsa::LCPArray* a, haplo::ScoreProvider* haplo_score_provider = nullptr);
    Mapper(void);
    ~Mapper(void);

    // a collection of read pairs which we'd like to realign once we have estimated the fragment_size
    vector<pair<Alignment, Alignment> > imperfect_pairs_to_retry;

    // Return true of the two alignments are consistent for paired reads, and false otherwise
    bool alignments_consistent(const map<string, double>& pos1,
                               const map<string, double>& pos2,
                               int fragment_size_bound);

    /// use the fragment length annotations to assess if the pair is consistent or not
    bool pair_consistent(Alignment& aln1, // may modify the alignments to store the reference positions
                         Alignment& aln2,
                         double pval);

    /// use the fragment configuration statistics to rescue more precisely
    pair<bool, bool> pair_rescue(Alignment& mate1, Alignment& mate2, bool& tried1, bool& tried2, int match_score, int full_length_bonus, bool traceback, bool xdrop_alignment);

    // uses heuristic clustering based on node id ranges to find alignment targets, and aligns
    vector<Alignment> mems_id_clusters_to_alignments(const Alignment& alignment, vector<MaximalExactMatch>& mems, int additional_multimaps);

    // use mapper parameters to determine which clusters we should drop
    set<const vector<MaximalExactMatch>* > clusters_to_drop(const vector<vector<MaximalExactMatch> >& clusters);

    /// Use the scoring provided by the internal aligner to re-score the
    /// alignment, scoring gaps between nodes using graph distance from the PathPositionHandleGraph
    /// index. Can use either approximate or exact (with approximate fallback)
    /// PathPositionHandleGraph-based distance estimation. Will strip out bonuses if the appropriate
    /// Mapper flag is set.
    /// Does not apply a haplotype consistency bonus, as this function is intended for alignments with large gaps.
    int32_t score_alignment(const Alignment& aln, bool use_approx_distance = false);
    
    /// Given an alignment scored with full length bonuses on, subtract out the full length bonus if it was applied.
    void remove_full_length_bonuses(Alignment& aln);
    
    // run through the alignment and attempt to align unaligned parts of the alignment to the graph in the region where they are anchored
    Alignment patch_alignment(const Alignment& aln, int max_patch_length, bool trim_internal_deletions = true, bool xdrop_alignment = false);
    // for aligning to a particular MEM cluster
    Alignment align_cluster(const Alignment& aln, const vector<MaximalExactMatch>& mems, bool traceback, bool xdrop_alignment = false);
    // compute the uniqueness metric based on the MEMs in the cluster
    double compute_uniqueness(const Alignment& aln, const vector<MaximalExactMatch>& mems);
    // wraps align_to_graph with flipping
    Alignment align_maybe_flip(const Alignment& base, HandleGraph& graph, bool flip, bool traceback, bool banded_global = false, bool xdrop_alignment = false);
    Alignment align_maybe_flip(const Alignment& base, HandleGraph& graph, const vector<MaximalExactMatch>& mems, bool flip, bool traceback, bool banded_global = false, bool xdrop_alignment = false);

    bool adjacent_positions(const Position& pos1, const Position& pos2);
    int64_t get_node_length(int64_t node_id);
    bool check_alignment(const Alignment& aln);
    
    // Align the given string and return an Alignment.
    Alignment align(const string& seq,
                    int kmer_size = 0,
                    int stride = 0,
                    int max_mem_length = 0,
                    int band_width = 1000,
                    int band_overlap = 500,
                    bool xdrop_alignment = false);

    // Align the given read and return an aligned copy. Does not modify the input Alignment.
    Alignment align(const Alignment& read,
                    int kmer_size = 0,
                    int stride = 0,
                    int max_mem_length = 0,
                    int band_width = 1000,
                    int band_overlap = 500,
                    bool xdrop_alignment = false);

    // Align the given read with multi-mapping. Returns the alignments in score
    // order, up to multimaps (or max_multimaps if multimaps is 0). Does not update the alignment passed in.
    // If the sequence is longer than the band_width, will only produce a single best banded alignment.
    // All alignments but the first are marked as secondary.
    vector<Alignment> align_multi(const Alignment& aln,
                                  int kmer_size = 0,
                                  int stride = 0,
                                  int max_mem_length = 0,
                                  int band_width = 1000,
                                  int band_overlap = 500,
                                  bool xdrop_alignment = false);
    
    // paired-end based
    
    // Both vectors of alignments will be sorted in order of increasing score.
    // All alignments but the first in each vector are marked as secondary.
    // Alignments at corresponding positions in the two vectors may or may not
    // be corresponding paired alignments. If a read does not map, its vector
    // will be empty.
    // If only_top_scoring_pair is set, then the vectors will be empty unless
    // the primary pair of alignments each have top scores individually as well. 
    // align the pair as a single component using MEM threading and patching on the pair simultaneously
    pair<vector<Alignment>, vector<Alignment>> 
        align_paired_multi(const Alignment& read1,
                           const Alignment& read2,
                           bool& queued_resolve_later,
                           int max_mem_length = 0,
                           bool only_top_scoring_pair = false,
                           bool retrying = false,
                           bool xdrop_alignment = false);
    
    // compute a mapping quality component based only on the MEMs we've obtained
    double compute_cluster_mapping_quality(const vector<vector<MaximalExactMatch> >& clusters, int read_length);
    // use an average length of an LCP to a parent in the suffix tree to estimate a mapping quality
    double estimate_max_possible_mapping_quality(int length, double min_diffs, double next_min_diffs);
    // absolute max possible mq
    double max_possible_mapping_quality(int length);
    // takes the min of graph_distance, approx_distance, and xindex->min_approx_path_distance()
    int64_t graph_mixed_distance_estimate(pos_t pos1, pos_t pos2, int64_t maximum);
    // use the offset in the sequence array to give an approximate distance
    int64_t approx_distance(pos_t pos1, pos_t pos2);
    // use the offset in the sequence array to get an approximate position
    int64_t approx_position(pos_t pos);
    // get the approximate position of the alignment or return -1 if it can't be had
    int64_t approx_alignment_position(const Alignment& aln);
    // get the end position of the alignment
    Position alignment_end_position(const Alignment& aln);
    // get the approximate distance between the starts of the alignments or return -1 if undefined
    int64_t approx_fragment_length(const Alignment& aln1, const Alignment& aln2);
    // get a set of positions that are likely based on the fragment model and the embedded paths
    vector<pos_t> likely_mate_positions(const Alignment& aln, bool is_first);

    // fargment length estimation
    unordered_map<path_handle_t, int64_t> min_pair_fragment_length(const Alignment& aln1, const Alignment& aln2);
    
    // mem mapper parameters
    //
    //int max_mem_length; // a mem must be <= this length
    int min_cluster_length; // a cluster needs this much sequence in it for us to consider it
    int max_attempts;  // maximum number of times to try to increase sensitivity or use a lower-hit subgraph
    int thread_extension; // add this many nodes in id space to the end of the thread when building thread into a subgraph
    int max_target_factor; // the maximum multiple of the read length we'll try to align to

    // multimapping
    int max_multimaps;
    // soft clip resolution
    int softclip_threshold; // if more than this many bp are clipped, try extension algorithm
    int max_softclip_iterations; // Extend no more than this many times (while softclips are getting shorter)
    float min_identity; // require that alignment identity is at least this much to accept alignment
    int min_banded_mq; // when aligning banded, treat bands with MQ < this as unaligned
    // paired-end consistency enforcement
    int extra_multimaps; // Extra mappings considered
    int min_multimaps; // Minimum number of multimappings
    int band_multimaps; // the number of multimaps for to attempt for each band in a banded alignment
    bool patch_alignments; // should we attempt alignment patching to resolve unaligned regions in banded alignment
    
    double maybe_mq_threshold; // quality below which we let the estimated mq kick in
    int max_cluster_mapping_quality; // the cap for cluster mapping quality
    bool use_cluster_mq; // should we use the cluster-based mapping quality component
    double identity_weight; // scale mapping quality by the alignment score identity to this power

    bool always_rescue; // Should rescue be attempted for all imperfect alignments?
    bool include_full_length_bonuses;
    int max_xdrop_gap_length;
    
    bool simultaneous_pair_alignment;
    int max_band_jump; // the maximum length edit we can detect via banded alignment
    float drop_chain; // drop chains shorter than this fraction of the longest overlapping chain
    float mq_overlap; // consider as alternative mappings any alignment with this overlap with our best
    int mate_rescues;

    double pair_rescue_hang_threshold;
    double pair_rescue_retry_threshold;
    
    // Keep track of fragment length distribution statistics
    FragmentLengthStatistics frag_stats;

    atomic<int> warned_about_chunking{0};
    
    // bench_t bench[4];
    // uint64_t counter[6];
};

// utility
const vector<string> balanced_kmers(const string& seq, int kmer_size, int stride);

set<pos_t> gcsa_nodes_to_positions(const vector<gcsa::node_type>& nodes);

int sub_overlaps_of_first_aln(const vector<Alignment>& alns, float overlap_fraction);

}

#endif
