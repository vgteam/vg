#ifndef VG_MAPPER_HPP_INCLUDED
#define VG_MAPPER_HPP_INCLUDED

#include <iostream>
#include <map>
#include <chrono>
#include <ctime>
#include "omp.h"
#include "vg.hpp"
#include "xg.hpp"
#include "index.hpp"
#include "gcsa/gcsa.h"
#include "gcsa/lcp.h"
#include "alignment.hpp"
#include "path.hpp"
#include "position.hpp"
#include "xg_position.hpp"
#include "lru_cache.h"
#include "json2pb.h"
#include "entropy.hpp"
#include "gssw_aligner.hpp"
#include "mem.hpp"
#include "cluster.hpp"
#include "graph.hpp"

namespace vg {

// uncomment to make vg map --debug very interesting
//#define debug_mapper

using namespace std;
    
enum MappingQualityMethod { Approx, Exact, None };

class Mapper;

// for banded long read alignment resolution

class AlignmentChainModelVertex {
public:
    Alignment* aln;
    vector<pair<AlignmentChainModelVertex*, double> > next_cost; // for forward
    vector<pair<AlignmentChainModelVertex*, double> > prev_cost; // for backward
    double weight;
    double score;
    map<string, double> positions;
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
    map<string, map<int64_t, vector<vector<AlignmentChainModelVertex>::iterator> > > positions;
    set<vector<AlignmentChainModelVertex>::iterator> redundant_vertexes;
    vector<Alignment> unaligned_bands;
    AlignmentChainModel(
        vector<vector<Alignment> >& bands,
        Mapper* mapper,
        const function<double(const Alignment&, const Alignment&, const map<string, double>&, const map<string, double>&)>& transition_weight,
        int vertex_width = 10,
        int position_depth = 1,
        int max_connections = 30);
    void score(const set<AlignmentChainModelVertex*>& exclude);
    AlignmentChainModelVertex* max_vertex(void);
    vector<Alignment> traceback(const Alignment& read, int alt_alns, bool paired, bool debug);
    void display(ostream& out);
    void clear_scores(void);
};

/*
 * A threadsafe class that keeps a running estimation of a fragment length distribution
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
    double stdev() const;
    
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
    
class BaseMapper : public Progressive {
    
public:
    // Make a Mapper that pulls from an XG succinct graph and a GCSA2 kmer index + LCP array
    BaseMapper(xg::XG* xidex, gcsa::GCSA* g, gcsa::LCPArray* a);
    BaseMapper(void);
    ~BaseMapper(void);
    
    double estimate_gc_content(void);
    
    int random_match_length(double chance_random);
    
    void set_alignment_scores(int8_t match, int8_t mismatch, int8_t gap_open, int8_t gap_extend, int8_t full_length_bonus);
    
    // TODO: setting alignment threads could mess up the internal memory for how many threads to reset to
    void set_fragment_length_distr_params(size_t maximum_sample_size = 1000, size_t reestimation_frequency = 1000,
                                          double robust_estimation_fraction = 0.95);
    
    /// Set the alignment thread count, updating internal data structures that
    /// are per thread. Note that this resets aligner scores to their default values!
    void set_alignment_threads(int new_thread_count);
    
    void set_cache_size(int new_cache_size);
    
    /// Returns true if fragment length distribution has been fixed
    bool has_fixed_fragment_length_distr();
    
    /// Use the given fragment length distribution parameters instead of
    /// estimating them.
    void force_fragment_length_distr(double mean, double stddev);
    
    // MEM-based mapping
    // find maximal exact matches
    // These are SMEMs by definition when shorter than the max_mem_length or GCSA2 order.
    // Designating reseed_length returns minimally-more-frequent sub-MEMs in addition to SMEMs when SMEM is >= reseed_length.
    // Minimally-more-frequent sub-MEMs are MEMs contained in an SMEM that have occurrences outside of the SMEM.
    // SMEMs and sub-MEMs will be automatically filled with the nodes they contain, which the occurrences of the sub-MEMs
    // that are inside SMEM hits filtered out. (filling sub-MEMs currently requires an XG index)
    
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
                   int reseed_below_count = 0);
    
    // Use the GCSA2 index to find super-maximal exact matches.
    vector<MaximalExactMatch>
    find_mems_simple(string::const_iterator seq_begin,
                     string::const_iterator seq_end,
                     int max_mem_length = 0,
                     int min_mem_length = 1,
                     int reseed_length = 0);
    
    
    int min_mem_length; // a mem must be >= this length
    int mem_reseed_length; // the length above which we reseed MEMs to get potentially missed hits
    bool fast_reseed; // use the fast reseed algorithm
    double fast_reseed_length_diff; // how much smaller than its parent a sub-MEM can be in the fast reseed algorithm
    int hit_max;       // ignore or MEMs with more than this many hits
    
    bool strip_bonuses; // remove any bonuses used by the aligners from the final reported scores
    bool assume_acyclic; // the indexed graph is acyclic
    bool adjust_alignments_for_base_quality; // use base quality adjusted alignments
    
    MappingQualityMethod mapping_quality_method; // how to compute mapping qualities
    int max_mapping_quality; // the cap for mapping quality
    
protected:
    /// Locate the sub-MEMs contained in the last MEM of the mems vector that have ending positions
    /// before the end the next SMEM, label each of the sub-MEMs with the indices of all of the SMEMs
    /// that contain it
    void find_sub_mems(const vector<MaximalExactMatch>& mems,
                       int mem_idx,
                       string::const_iterator next_mem_end,
                       int min_mem_length,
                       vector<pair<MaximalExactMatch, vector<size_t>>>& sub_mems_out);
    
    /// Provides same semantics as find_sub_mems but with a different algorithm. This algorithm uses the
    /// min_mem_length as a pruning tool instead of the LCP index. It can be expected to be faster when both
    /// the min_mem_length reasonably large relative to the reseed_length (e.g. 1/2 of SMEM size or similar).
    void find_sub_mems_fast(const vector<MaximalExactMatch>& mems,
                            int mem_idx,
                            string::const_iterator next_mem_end,
                            int min_sub_mem_length,
                            vector<pair<MaximalExactMatch, vector<size_t>>>& sub_mems_out);
    
    /// finds the nodes of sub MEMs that do not occur inside parent MEMs, each sub MEM should be associated
    /// with a vector of the indices of the SMEMs that contain it in the parent MEMs vector
    void fill_nonredundant_sub_mem_nodes(vector<MaximalExactMatch>& parent_mems,
                                         vector<pair<MaximalExactMatch, vector<size_t> > >::iterator sub_mem_records_begin,
                                         vector<pair<MaximalExactMatch, vector<size_t> > >::iterator sub_mem_records_end);
    
    /// fills a vector where each element contains the set of positions in the graph that the
    /// MEM touches at that index for the first MEM hit in the GCSA array
    void first_hit_positions_by_index(MaximalExactMatch& mem,
                                      vector<set<pos_t>>& positions_by_index_out);
    
    /// fills a vector where each element contains the set of positions in the graph that the
    /// MEM touches at that index starting at a given hit
    void mem_positions_by_index(MaximalExactMatch& mem, pos_t hit_pos,
                                vector<set<pos_t>>& positions_by_index_out);
    
    // use the xg index to get a character at a particular position (rc or foward)
    char pos_char(pos_t pos);
    
    // the next positions and their characters following the same strand of the graph
    map<pos_t, char> next_pos_chars(pos_t pos);
    
    // get the positions some specific distance from the given position (in the forward direction)
    set<pos_t> positions_bp_from(pos_t pos, int distance, bool rev);
    
    // Use the GCSA index to look up the sequence
    set<pos_t> sequence_positions(const string& seq);
    
    // debugging, checking of mems using find interface to gcsa
    void check_mems(const vector<MaximalExactMatch>& mems);
    
    int alignment_threads; // how many threads will *this* mapper use. Should not be set directly.
    /*
    int cache_size;
    
    // match walking support to prevent repeated calls to the xg index for the same node
    vector<LRUCache<id_t, Node>* > node_cache;
    LRUCache<id_t, Node>& get_node_cache(void);
    void init_node_cache(void);
    
    // node start cache for fast approximate position estimates
    vector<LRUCache<id_t, int64_t>* > node_start_cache;
    LRUCache<id_t, int64_t>& get_node_start_cache(void);
    void init_node_start_cache(void);
    
    // match node traversals to path positions
    vector<LRUCache<gcsa::node_type, map<string, vector<size_t> > >* > node_pos_cache;
    LRUCache<gcsa::node_type, map<string, vector<size_t> > >& get_node_pos_cache(void);
    void init_node_pos_cache(void);
    
    vector<LRUCache<id_t, vector<Edge> >* > edge_cache;
    LRUCache<id_t, vector<Edge> >& get_edge_cache(void);
    void init_edge_cache(void);
    */
    
    void init_aligner(int8_t match, int8_t mismatch, int8_t gap_open, int8_t gap_extend, int8_t full_length_bonus);
    void clear_aligners(void);
    
    // xg index
    xg::XG* xindex = nullptr;
    
    // GCSA index and its LCP array
    gcsa::GCSA* gcsa = nullptr;
    gcsa::LCPArray* lcp = nullptr;
    
    FragmentLengthDistribution fragment_length_distr;

    /// Get the appropriate aligner to use, based on
    /// adjust_alignments_for_base_quality. By setting have_qualities to false,
    /// you can force the non-quality-adjusted aligner, for reads that lack
    /// quality scores.
    BaseAligner* get_aligner(bool have_qualities = true) const;
    
    // Sometimes you really do need the two kinds of aligners, to pass to code
    // that expects one or the other.
    QualAdjAligner* get_qual_adj_aligner() const;
    Aligner* get_regular_aligner() const;

private:
    // GSSW aligners
    QualAdjAligner* qual_adj_aligner = nullptr;
    Aligner* regular_aligner = nullptr;    

};

/**
 * Keeps track of statistics about fragment length within the Mapper class.
 * Belongs to a single thread.
 */
class FragmentLengthStatistics {
public:

    void record_fragment_configuration(int length, const Alignment& aln1, const Alignment& aln2);
    
    string fragment_model_str(void);
    void save_frag_lens_to_alns(Alignment& aln1, Alignment& aln2, const map<string, int64_t>& approx_frag_lengths, bool is_consistent);
    
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
    bool cached_fragment_orientation = 0;
    bool cached_fragment_direction = 1;
    
    // These variables are used to manage the periodic updates
    int64_t since_last_fragment_length_estimate = 0;
    int64_t fragment_model_update_interval = 10;
    
    // These deques are used for the periodic running estimation of the fragment length distribution
    deque<double> fragment_lengths;
    deque<bool> fragment_orientations;
    deque<bool> fragment_directions;

    int64_t fragment_max = 10000; // the maximum length fragment which we will consider when estimating fragment lengths
    int64_t fragment_size = 0; // Used to bound clustering of MEMs during paired end mapping, also acts as sentinel to determine
                       // if consistent pairs should be reported; dynamically estimated at runtime
    double fragment_sigma = 10; // the number of times the standard deviation above the mean to set the fragment_size
    int64_t fragment_length_cache_size = 1000;
    float perfect_pair_identity_threshold = 0.95;
    bool fixed_fragment_model = true;
    
    
    
    
    
};

class Mapper : public BaseMapper {


private:
    
    Alignment align_to_graph(const Alignment& aln,
                             Graph& graph,
                             size_t max_query_graph_ratio,
                             bool traceback,
                             bool pinned_alignment = false,
                             bool pin_left = false,
                             bool global = false);
    vector<Alignment> align_multi_internal(bool compute_unpaired_qualities,
                                           const Alignment& aln,
                                           int kmer_size,
                                           int stride,
                                           int max_mem_length,
                                           int band_width,
                                           double& cluster_mq,
                                           int keep_multimaps = 0,
                                           int additional_multimaps = 0,
                                           vector<MaximalExactMatch>* restricted_mems = nullptr);
    void compute_mapping_qualities(vector<Alignment>& alns, double cluster_mq, double mq_estimate, double mq_cap);
    void compute_mapping_qualities(pair<vector<Alignment>, vector<Alignment>>& pair_alns, double cluster_mq, double mq_estmate1, double mq_estimate2, double mq_cap1, double mq_cap2);
    vector<Alignment> score_sort_and_deduplicate_alignments(vector<Alignment>& all_alns, const Alignment& original_alignment);
    void filter_and_process_multimaps(vector<Alignment>& all_alns, int total_multimaps);
    // make the bands used in banded alignment
    vector<Alignment> make_bands(const Alignment& read, int band_width, vector<pair<int, int>>& to_strip);
    // Return the one best banded alignment.
    vector<Alignment> align_banded(const Alignment& read,
                                   int kmer_size = 0,
                                   int stride = 0,
                                   int max_mem_length = 0,
                                   int band_width = 1000);
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
                                      int additional_multimaps);
    
public:
    // Make a Mapper that pulls from an XG succinct graph and a GCSA2 kmer index + LCP array
    Mapper(xg::XG* xidex, gcsa::GCSA* g, gcsa::LCPArray* a);
    Mapper(void);
    ~Mapper(void);

    map<string, vector<size_t> > node_positions_in_paths(gcsa::node_type node);
    
    // a collection of read pairs which we'd like to realign once we have estimated the fragment_size
    vector<pair<Alignment, Alignment> > imperfect_pairs_to_retry;

    double graph_entropy(void);

    // use the xg index to get the mean position of the nodes in the alignent for each reference that it corresponds to
    map<string, double> alignment_mean_path_positions(const Alignment& aln, bool first_hit_only = true);
    void annotate_with_mean_path_positions(vector<Alignment>& alns);
    
    // use the xg index to get the first position of an alignment on a reference path
    map<string, size_t> alignment_initial_path_positions(const Alignment& aln);
    void annotate_with_initial_path_positions(Alignment& aln);

    // Return true of the two alignments are consistent for paired reads, and false otherwise
    bool alignments_consistent(const map<string, double>& pos1,
                               const map<string, double>& pos2,
                               int fragment_size_bound);

    /// use the fragment length annotations to assess if the pair is consistent or not
    bool pair_consistent(const Alignment& aln1,
                         const Alignment& aln2,
                         double pval);

    /// use the fragment configuration statistics to rescue more precisely
    pair<bool, bool> pair_rescue(Alignment& mate1, Alignment& mate2, int match_score, int full_length_bonus, bool traceback);

    /// assuming the read has only been score-aligned, realign from the end position backwards
    Alignment realign_from_start_position(const Alignment& aln, int extra, int iteration);
    
    set<MaximalExactMatch*> resolve_paired_mems(vector<MaximalExactMatch>& mems1,
                                                vector<MaximalExactMatch>& mems2);

    // uses heuristic clustering based on node id ranges to find alignment targets, and aligns
    vector<Alignment> mems_id_clusters_to_alignments(const Alignment& alignment, vector<MaximalExactMatch>& mems, int additional_multimaps);

    // use mapper parameters to determine which clusters we should drop
    set<const vector<MaximalExactMatch>* > clusters_to_drop(const vector<vector<MaximalExactMatch> >& clusters);

    // takes the input alignment (with seq, etc) so we have reference to the base sequence
    // for reconstruction the alignments from the SMEMs
    Alignment mems_to_alignment(const Alignment& aln, vector<MaximalExactMatch>& mems);
    Alignment mem_to_alignment(MaximalExactMatch& mem);
    
    /// Use the scoring provided by the internal aligner to re-score the
    /// alignment, scoring gaps between nodes using graph distance from the XG
    /// index. Can use either approximate or exact (with approximate fallback)
    /// XG-based distance estimation. Will strip out bonuses if the appropriate
    /// Mapper flag is set.
    int32_t score_alignment(const Alignment& aln, bool use_approx_distance = false);
    
    /// Given an alignment scored with full length bonuses on, subtract out the full length bonus if it was applied.
    int32_t remove_full_length_bonus(const Alignment& aln);
    
    // run through the alignment and attempt to align unaligned parts of the alignment to the graph in the region where they are anchored
    Alignment patch_alignment(const Alignment& aln, int max_patch_length);
    // Get the graph context of a particular cluster, not expanding beyond the middles of MEMs.
    VG cluster_subgraph_strict(const Alignment& aln, const vector<MaximalExactMatch>& mems);
    // for aligning to a particular MEM cluster
    Alignment align_cluster(const Alignment& aln, const vector<MaximalExactMatch>& mems, bool traceback);
    // compute the uniqueness metric based on the MEMs in the cluster
    double compute_uniqueness(const Alignment& aln, const vector<MaximalExactMatch>& mems);
    // wraps align_to_graph with flipping
    Alignment align_maybe_flip(const Alignment& base, Graph& graph, bool flip, bool traceback, bool banded_global = false);

    bool adjacent_positions(const Position& pos1, const Position& pos2);
    int64_t get_node_length(int64_t node_id);
    bool check_alignment(const Alignment& aln);
    VG alignment_subgraph(const Alignment& aln, int context_size = 1);
    
    // Align the given string and return an Alignment.
    Alignment align(const string& seq,
                    int kmer_size = 0,
                    int stride = 0,
                    int max_mem_length = 0,
                    int band_width = 1000);

    // Align the given read and return an aligned copy. Does not modify the input Alignment.
    Alignment align(const Alignment& read,
                    int kmer_size = 0,
                    int stride = 0,
                    int max_mem_length = 0,
                    int band_width = 1000);

    // Align the given read with multi-mapping. Returns the alignments in score
    // order, up to multimaps (or max_multimaps if multimaps is 0). Does not update the alignment passed in.
    // If the sequence is longer than the band_width, will only produce a single best banded alignment.
    // All alignments but the first are marked as secondary.
    vector<Alignment> align_multi(const Alignment& aln,
                                  int kmer_size = 0,
                                  int stride = 0,
                                  int max_mem_length = 0,
                                  int band_width = 1000);
    
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
                           bool retrying = false);

    // lossily project an alignment into a particular path space of a graph
    // the resulting alignment is equivalent to a SAM record against the chosen path
    Alignment surject_alignment(const Alignment& source,
                                set<string>& path_names,
                                string& path_name,
                                int64_t& path_pos,
                                bool& path_reverse,
                                int window);

    
    // compute a mapping quality component based only on the MEMs we've obtained
    double compute_cluster_mapping_quality(const vector<vector<MaximalExactMatch> >& clusters, int read_length);
    // use an average length of an LCP to a parent in the suffix tree to estimate a mapping quality
    double estimate_max_possible_mapping_quality(int length, double min_diffs, double next_min_diffs);
    // walks the graph one base at a time from pos1 until we find pos2
    int64_t graph_distance(pos_t pos1, pos_t pos2, int64_t maximum = 1e3);
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
    // use the cached fragment model to estimate the likely place we'll find the mate
    pos_t likely_mate_position(const Alignment& aln, bool is_first);
    // get the node approximately at the given offset relative to our position (offset may be negative)
    id_t node_approximately_at(int64_t approx_pos);
    // convert a single MEM hit into an alignment (by definition, a perfect one)
    Alignment walk_match(const string& seq, pos_t pos);
    vector<Alignment> walk_match(const Alignment& base, const string& seq, pos_t pos);
    // convert the set of hits of a MEM into a set of alignments
    vector<Alignment> mem_to_alignments(MaximalExactMatch& mem);

    // fargment length estimation
    map<string, int64_t> approx_pair_fragment_length(const Alignment& aln1, const Alignment& aln2);
    int64_t first_approx_pair_fragment_length(const Alignment& aln1, const Alignment& aln2);
    // uses the cached information about the graph in the xg index to get an approximate node length
    double average_node_length(void);
    
    bool debug;

    // mem mapper parameters
    //
    //int max_mem_length; // a mem must be <= this length
    int min_cluster_length; // a cluster needs this much sequence in it for us to consider it
    int context_depth; // how deeply the mapper will extend out the subgraph prior to alignment
    int max_attempts;  // maximum number of times to try to increase sensitivity or use a lower-hit subgraph
    int thread_extension; // add this many nodes in id space to the end of the thread when building thread into a subgraph
    int max_target_factor; // the maximum multiple of the read length we'll try to align to

    size_t max_query_graph_ratio;

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
    
    int maybe_mq_threshold; // quality below which we let the estimated mq kick in
    int max_cluster_mapping_quality; // the cap for cluster mapping quality
    bool use_cluster_mq; // should we use the cluster-based mapping quality component
    double identity_weight; // scale mapping quality by the alignment score identity to this power

    bool always_rescue; // Should rescue be attempted for all imperfect alignments?
    
    bool simultaneous_pair_alignment;
    int max_band_jump; // the maximum length edit we can detect via banded alignment
    float drop_chain; // drop chains shorter than this fraction of the longest overlapping chain
    float mq_overlap; // consider as alternative mappings any alignment with this overlap with our best
    int mate_rescues;

    double pair_rescue_hang_threshold;
    double pair_rescue_retry_threshold;
    
    // Keep track of fragment length distribution statistics
    FragmentLengthStatistics frag_stats;

};

// utility
const vector<string> balanced_kmers(const string& seq, int kmer_size, int stride);

set<pos_t> gcsa_nodes_to_positions(const vector<gcsa::node_type>& nodes);

int sub_overlaps_of_first_aln(const vector<Alignment>& alns, float overlap_fraction);

}

#endif
