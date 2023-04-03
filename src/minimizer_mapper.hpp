#ifndef VG_MINIMIZER_MAPPER_HPP_INCLUDED
#define VG_MINIMIZER_MAPPER_HPP_INCLUDED

/** 
 * \file minimizer_mapper.hpp
 * Defines a mapper that uses the minimizer index and GBWT-based extension.
 */

#include "algorithms/chain_items.hpp"
#include "algorithms/nearest_offsets_in_paths.hpp"
#include "aligner.hpp"
#include "vg/io/alignment_emitter.hpp"
#include "gbwt_extender.hpp"
#include "snarl_seed_clusterer.hpp"
#include "mapper.hpp"
#include "snarls.hpp"
#include "tree_subgraph.hpp"
#include "funnel.hpp"

#include <gbwtgraph/minimizer.h>
#include <structures/immutable_list.hpp>

#include <atomic>

namespace vg {

//#define debug_chaining

using namespace std;
using namespace vg::io;

class MinimizerMapper : public AlignerClient {
public:

    /**
     * Construct a new MinimizerMapper using the given indexes. The PathPositionhandleGraph can be nullptr,
     * as we only use it for correctness tracking.
     */

    MinimizerMapper(const gbwtgraph::GBWTGraph& graph,
         const gbwtgraph::DefaultMinimizerIndex& minimizer_index,
         SnarlDistanceIndex* distance_index,
         const vector<ZipCode>* zipcodes,
         const PathPositionHandleGraph* path_graph = nullptr);

    /**
     * Map the given read, and send output to the given AlignmentEmitter. May be run from any thread.
     * TODO: Can't be const because the clusterer's cluster_seeds isn't const.
     */
    void map(Alignment& aln, AlignmentEmitter& alignment_emitter);
    
    /**
     * Map the given read. Return a vector of alignments that it maps to, winner first.
     */
    vector<Alignment> map(Alignment& aln);
    
    /**
     * Map the given read using chaining of seeds. Return a vector of alignments that it maps to, winner first.
     */
    vector<Alignment> map_from_chains(Alignment& aln);
    
    /**
     * Map the given read using gapless extensions. Return a vector of alignments that it maps to, winner first.
     */
    vector<Alignment> map_from_extensions(Alignment& aln);
    
    // The idea here is that the subcommand feeds all the reads to the version
    // of map_paired that takes a buffer, and then empties the buffer by
    // iterating over it in parallel with the version that doesn't.
    // TODO: how will we warn about not having a pair distribution yet then?
    
    /**
     * Map the given pair of reads, where aln1 is upstream of aln2 and they are
     * oriented towards each other in the graph.
     *
     * If the reads are ambiguous and there's no fragment length distribution
     * fixed yet, they will be dropped into ambiguous_pair_buffer.
     *
     * Otherwise, at least one result will be returned for them (although it
     * may be the unmapped alignment).
     */
    pair<vector<Alignment>, vector<Alignment>> map_paired(Alignment& aln1, Alignment& aln2,
        vector<pair<Alignment, Alignment>>& ambiguous_pair_buffer);
        
    /**
     * Map the given pair of reads, where aln1 is upstream of aln2 and they are
     * oriented towards each other in the graph.
     *
     * If the fragment length distribution is not yet fixed, reads will be
     * mapped independently. Otherwise, they will be mapped according to the
     * fragment length distribution.
     */
    pair<vector<Alignment>, vector<Alignment>> map_paired(Alignment& aln1, Alignment& aln2);




    // Mapping settings.
    // TODO: document each

    /// Use all minimizers with at most hit_cap hits
    static constexpr size_t default_hit_cap = 10;
    size_t hit_cap = default_hit_cap;
    
    /// Ignore all minimizers with more than hard_hit_cap hits
    static constexpr size_t default_hard_hit_cap = 500;
    size_t hard_hit_cap = default_hard_hit_cap;
    
    /// Take minimizers between hit_cap and hard_hit_cap hits until this fraction
    /// of total score
    static constexpr double default_minimizer_score_fraction = 0.9;
    double minimizer_score_fraction = default_minimizer_score_fraction;

    /// Maximum number of distinct minimizers to take
    static constexpr size_t default_max_unique_min = 500;
    size_t max_unique_min = default_max_unique_min;
    
    /// Number of minimzers to select based on read_len/num_min_per_bp
    static constexpr size_t default_num_bp_per_min = 1000;
    size_t num_bp_per_min = default_num_bp_per_min;

    /// If set, exclude overlapping minimizers
    static constexpr bool default_exclude_overlapping_min = false;
    bool exclude_overlapping_min = default_exclude_overlapping_min;
    
    //////////////
    // Alignment-from-gapless-extension/short read Giraffe specific parameters:
    //////////////

    ///Accept at least this many clusters for gapless extension
    static constexpr size_t default_min_extensions = 2;
    size_t min_extensions = default_min_extensions;

    /// How many clusters should we produce gapless extensions for, max?
    static constexpr size_t default_max_extensions = 800;
    size_t max_extensions = default_max_extensions;

    //If an extension set's score is smaller than the best 
    //extension's score by more than this much, don't align it
    static constexpr double default_extension_set_score_threshold = 20;
    double extension_set_score_threshold = default_extension_set_score_threshold;

    //If an extension's score is smaller than the best extension's score by
    //more than this much, don't align it
    static constexpr int default_extension_score_threshold = 1;
    int extension_score_threshold = default_extension_score_threshold;
    
    /// Disregard the extension set score thresholds when they would give us
    /// fewer than this many extension sets.
    static constexpr int default_min_extension_sets = 2;
    int min_extension_sets = default_min_extension_sets;
    
    /// Even if we would have fewer than min_extension_sets results, don't
    /// process anything with a score smaller than this.
    static constexpr int default_extension_set_min_score = 20;
    int extension_set_min_score = default_extension_set_min_score;
    
    /////////////////
    // More shared parameters:
    /////////////////

    /// How many extended clusters should we align, max?
    static constexpr size_t default_max_alignments = 8;
    size_t max_alignments = default_max_alignments;
    
    /// How many extensions should we try as seeds within a mapping location?
    static constexpr size_t default_max_local_extensions = numeric_limits<size_t>::max();
    size_t max_local_extensions = default_max_local_extensions;

    /// If a cluster's score is smaller than the best score of any cluster by more than
    /// this much, then don't extend it
    static constexpr double default_cluster_score_threshold = 50;
    double cluster_score_threshold = default_cluster_score_threshold;
    
    /// If the second best cluster's score is no more than this many points below
    /// the cutoff set by cluster_score_threshold, snap that cutoff down to the
    /// second best cluster's score, to avoid throwing away promising
    /// secondaries.
    static constexpr double default_pad_cluster_score_threshold = 20;
    double pad_cluster_score_threshold = default_pad_cluster_score_threshold;

    /// If the read coverage of a cluster is less than the best coverage of any cluster
    /// by more than this much, don't extend it
    static constexpr double default_cluster_coverage_threshold = 0.3;
    double cluster_coverage_threshold = default_cluster_coverage_threshold;
    
    //////////////////
    // Alignment-from-chains/long read Giraffe specific parameters:
    //////////////////
    
    /// If true, produce alignments from extension sets by chaining gapless
    /// extensions up and aligning the sequences between them. If false,
    /// produce alignments by aligning the tails off of individual gapless
    /// extensions.
    static constexpr bool default_align_from_chains = false;
    bool align_from_chains = default_align_from_chains;
    
    /// What multiple of the read length should we use for bucketing (coarse clustering/preclustering)?
    static constexpr double default_bucket_scale = 2.0;
    double bucket_scale = default_bucket_scale;
    
    /// How many fragments should we try and make in every bucket?
    static constexpr size_t default_max_fragments_per_bucket = std::numeric_limits<size_t>::max();
    size_t max_fragments_per_bucket = default_max_fragments_per_bucket;
    
    /// How many bases should we look back when making fragments?
    static constexpr size_t default_fragment_max_lookback_bases = 400;
    size_t fragment_max_lookback_bases = default_fragment_max_lookback_bases;
    /// In fragments, how many sources should we make sure to consider regardless of distance?
    static constexpr size_t default_fragment_min_lookback_items = 0;
    size_t fragment_min_lookback_items = default_fragment_min_lookback_items;
    /// In fragments, how many sources should we allow ourselves to consider ever?
    static constexpr size_t default_fragment_lookback_item_hard_cap = 3;
    size_t fragment_lookback_item_hard_cap = default_fragment_lookback_item_hard_cap;
    
    /// If the read coverage of a fragment connection is less than the best of any
    /// by more than this much, don't extend it
    static constexpr double default_fragment_connection_coverage_threshold = 0.3;
    double fragment_connection_coverage_threshold = default_fragment_connection_coverage_threshold;
    
    /// How many connections between fragments should we reseed over, minimum?
    static constexpr size_t default_min_fragment_connections = 10;
    size_t min_fragment_connections = default_min_fragment_connections;
    
    /// How many connections between fragments should we reseed over, maximum?
    static constexpr size_t default_max_fragment_connections = 50;
    size_t max_fragment_connections = default_max_fragment_connections;
    
    /// When connecting subclusters for reseeding, how far should we search?
    static constexpr size_t default_reseed_search_distance = 10000;
    size_t reseed_search_distance = default_reseed_search_distance;
    
    /// What read-length-independent distance threshold do we want to use for final clustering?
    static constexpr size_t default_chaining_cluster_distance = 100;
    size_t chaining_cluster_distance = default_chaining_cluster_distance;
    
    /// How many clusters should we produce chains for, max?
    static constexpr size_t default_max_buckets_to_fragment = 2;
    size_t max_buckets_to_fragment = default_max_buckets_to_fragment;

    /// When converting chains to alignments, what's the longest gap between
    /// items we will actually try to align? Passing strings longer than ~100bp
    /// can cause WFAAligner to run for a pathologically long amount of time.
    /// May not be 0.
    static constexpr size_t default_max_chain_connection = 100;
    size_t max_chain_connection = default_max_chain_connection;
    /// Similarly, what is the maximum tail length we will try to align?
    static constexpr size_t default_max_tail_length = 100;
    size_t max_tail_length = default_max_tail_length;
    
    /// How many bases should we look back when chaining?
    static constexpr size_t default_max_lookback_bases = 10000;
    size_t max_lookback_bases = default_max_lookback_bases;
    /// How many chaining sources should we make sure to consider regardless of distance?
    static constexpr size_t default_min_lookback_items = 1;
    size_t min_lookback_items = default_min_lookback_items;
    /// How many chaining sources should we allow ourselves to consider ever?
    static constexpr size_t default_lookback_item_hard_cap = 15;
    size_t lookback_item_hard_cap = default_lookback_item_hard_cap;
    /// How many bases should we try to look back initially when chaining?
    static constexpr size_t default_initial_lookback_threshold = 10;
    size_t initial_lookback_threshold = default_initial_lookback_threshold;
    /// How much chould we increase lookback when we can't find anything good?
    static constexpr double default_lookback_scale_factor = 2.0;
    double lookback_scale_factor = default_lookback_scale_factor;
    /// How bad can a transition be per base before lookback accepts it?
    static constexpr double default_min_good_transition_score_per_base = -0.1;
    double min_good_transition_score_per_base = default_min_good_transition_score_per_base;
    /// How much of a bonus should we give to each item in chaining?
    static constexpr int default_item_bonus = 0;
    int item_bonus = default_item_bonus;
    /// How many bases of indel should we allow in chaining?
    static constexpr size_t default_max_indel_bases = 6000;
    size_t max_indel_bases = default_max_indel_bases;
    
    /// If a chain's score is smaller than the best 
    /// chain's score by more than this much, don't align it
    static constexpr double default_chain_score_threshold = 100;
    double chain_score_threshold = default_chain_score_threshold;
    
    /// Disregard the chain score thresholds when they would give us
    /// fewer than this many chains.
    static constexpr int default_min_chains = 1;
    int min_chains = default_min_chains;
    
    /// Even if we would have fewer than min_chains results, don't
    /// process anything with a score smaller than this.
    static constexpr int default_chain_min_score = 100;
    int chain_min_score = default_chain_min_score;
    
    /// How long of a DP can we do before GSSW crashes due to 16-bit score
    /// overflow?
    static constexpr int MAX_DP_LENGTH = 30000;
    
    /// How many DP cells should we be willing to do in GSSW for an end-pinned
    /// alignment? If we want to do more than this, just leave tail unaligned.
    static constexpr size_t default_max_dp_cells = 16UL * 1024UL * 1024UL;
    size_t max_dp_cells = default_max_dp_cells;
    
    /////////////////
    // More shared parameters:
    /////////////////
    
    static constexpr size_t default_max_multimaps = 1;
    size_t max_multimaps = default_max_multimaps;
    static constexpr size_t default_distance_limit = 200;
    size_t distance_limit = default_distance_limit;
    
    /// If false, skip computing base-level alignments.
    static constexpr bool default_do_dp = true;
    bool do_dp = default_do_dp;
    
    /// Track which internal work items came from which others during each
    /// stage of the mapping algorithm.
    static constexpr bool default_track_provenance = false;
    bool track_provenance = default_track_provenance;

    /// Guess which seed hits are correct by location in the linear reference
    /// and track if/when their descendants make it through stages of the
    /// algorithm. Only works if track_provenance is true.
    static constexpr bool default_track_correctness = false;
    bool track_correctness = default_track_correctness;
    
    /// If set, log what the mapper is thinking in its mapping of each read.
    static constexpr bool default_show_work = false;
    bool show_work = default_show_work;

    ////How many stdevs from fragment length distr mean do we cluster together?
    static constexpr double default_paired_distance_stdevs = 2.0;
    double paired_distance_stdevs = default_paired_distance_stdevs; 

    ///How close does an alignment have to be to the best alignment for us to rescue on it
    static constexpr double default_paired_rescue_score_limit = 0.9;
    double paired_rescue_score_limit = default_paired_rescue_score_limit;

    ///How many stdevs from the mean do we extract a subgraph from?
    static constexpr double default_rescue_subgraph_stdevs = 4.0;
    double rescue_subgraph_stdevs = default_rescue_subgraph_stdevs;

    /// Do not attempt rescue if there are more seeds in the rescue subgraph.
    static constexpr size_t default_rescue_seed_limit = 100;
    size_t rescue_seed_limit = default_rescue_seed_limit;

    /// For paired end mapping, how many times should we attempt rescue (per read)?
    static constexpr size_t default_max_rescue_attempts = 15;
    size_t max_rescue_attempts = default_max_rescue_attempts;
    
    /// How big of an alignment in POA cells should we ever try to do with Dozeu?
    /// TODO: Lift this when Dozeu's allocator is able to work with >4 MB of memory.
    /// Each cell is 16 bits in Dozeu, and we leave some room for the query and
    /// padding to full SSE registers. Note that a very chopped graph might
    /// still break this!
    static constexpr size_t default_max_dozeu_cells = (size_t)(1.5 * 1024 * 1024);
    size_t max_dozeu_cells = default_max_dozeu_cells;
    
    ///What is the maximum fragment length that we accept as valid for paired-end reads?
    static constexpr size_t default_max_fragment_length = 2000;
    size_t max_fragment_length = default_max_fragment_length;
    
    /// Implemented rescue algorithms: no rescue, dozeu, GSSW.
    enum RescueAlgorithm { rescue_none, rescue_dozeu, rescue_gssw };

    /// The algorithm used for rescue.
    RescueAlgorithm rescue_algorithm = rescue_dozeu;
    
    /// Apply this sample name
    string sample_name;
    /// Apply this read group name
    string read_group;
    
    /// Have we complained about hitting the size limit for rescue?
    atomic_flag warned_about_rescue_size = ATOMIC_FLAG_INIT;
    
    /// Have we complained about hitting the size limit for tails?
    mutable atomic_flag warned_about_tail_size = ATOMIC_FLAG_INIT;

    bool fragment_distr_is_finalized () {return fragment_length_distr.is_finalized();}
    void finalize_fragment_length_distr() {
        if (!fragment_length_distr.is_finalized()) {
            fragment_length_distr.force_parameters(fragment_length_distr.mean(), fragment_length_distr.std_dev());
        } 
    }
    void force_fragment_length_distr(double mean, double stdev) {
        fragment_length_distr.force_parameters(mean, stdev);
    }
    double get_fragment_length_mean() const { return fragment_length_distr.mean(); }
    double get_fragment_length_stdev() const {return fragment_length_distr.std_dev(); }
    size_t get_fragment_length_sample_size() const { return fragment_length_distr.curr_sample_size(); }

    /**
     * Get the distance limit for the given read length
     */
    size_t get_distance_limit(size_t read_length) const {
        return max(distance_limit, read_length + 50);
    }

    /// The information we store for each seed.
    typedef SnarlDistanceIndexClusterer::Seed Seed;
    
    /**
     * We define our own type for minimizers, to use during mapping and to pass around between our internal functions.
     * Also used to represent syncmers, in which case the only window, the "minimizer", and the agglomeration are all the same region.
     */
    struct Minimizer {
        typename gbwtgraph::DefaultMinimizerIndex::minimizer_type value;
        size_t agglomeration_start; // What is the start base of the first window this minimizer instance is minimal in?
        size_t agglomeration_length; // What is the length in bp of the region of consecutive windows this minimizer instance is minimal in?
        size_t hits; // How many hits does the minimizer have?
        const gbwtgraph::hit_type* occs;
        int32_t length; // How long is the minimizer (index's k)
        int32_t candidates_per_window; // How many minimizers compete to be the best (index's w), or 1 for syncmers.  
        double score; // Scores as 1 + ln(hard_hit_cap) - ln(hits).

        // Sort the minimizers in descending order by score and group identical minimizers together.
        inline bool operator< (const Minimizer& another) const {
            return (this->score > another.score || (this->score == another.score && this->value.key < another.value.key));
        }
        
        /// Get the starting position of the given minimizer on the forward strand.
        /// Use this instead of value.offset which can really be the last base for reverse strand minimizers.
        inline size_t forward_offset() const {
            if (this->value.is_reverse) {
                // We have the position of the last base and we need the position of the first base.
                return this->value.offset - (this->length - 1);
            } else {
                // We already have the position of the first base.
                return this->value.offset;
            }
        }
        
        /// How many bases are in a window for which a minimizer is chosen?
        inline size_t window_size() const {
            return length + candidates_per_window - 1;
        }
        
        /// How many different windows are in this minimizer's agglomeration?
        inline size_t agglomeration_window_count() const {
            // Work out the length of a whole window, and then from that and the window count get the overall length.
            return agglomeration_length - window_size() + 1;
        }
        
        /// What is the minimizer sequence, in read orientation?
        inline string forward_sequence() const {
            string sequence = value.key.decode(length);
            return value.is_reverse ? reverse_complement(sequence) : sequence;
        }
    };
    
protected:
    
    /// Convert an integer distance, with limits standing for no distance, to a
    /// double annotation that can safely be parsed back from JSON into an
    /// integer if it is integral.
    double distance_to_annotation(int64_t distance) const;
    
    /// How should we initialize chain info when it's not stored in the minimizer index?
    inline static gbwtgraph::payload_type no_chain_info() {
        return MIPayload::NO_CODE;  
    } 
    
    /// How do we convert chain info to an actual seed of the type we are using?
    /// Also needs to know the hit position, and the minimizer number.
    inline static Seed chain_info_to_seed(const pos_t& hit, size_t minimizer, const ZipCode& zip, ZipCodeDecoder* decoder) {
        return { hit, minimizer, zip, std::unique_ptr<ZipCodeDecoder>(decoder)};
    }
    
    /// Convert a collection of seeds to a collection of chaining anchors.
    std::vector<algorithms::Anchor> to_anchors(const Alignment& aln, const VectorView<Minimizer>& minimizers, const std::vector<Seed>& seeds) const;
    
    /// Convert a single seed to a single chaining anchor.
    algorithms::Anchor to_anchor(const Alignment& aln, const VectorView<Minimizer>& minimizers, const Seed& seed) const;
    
    /// Convert an Anchor to a WFAAlignment
    WFAAlignment to_wfa_alignment(const algorithms::Anchor& anchor) const; 

    /// The information we store for each cluster.
    typedef SnarlDistanceIndexClusterer::Cluster Cluster;

    // These are our indexes
    const PathPositionHandleGraph* path_graph; // Can be nullptr; only needed for correctness tracking.
    const gbwtgraph::DefaultMinimizerIndex& minimizer_index;
    SnarlDistanceIndex* distance_index;
    const vector<ZipCode>* zipcodes;
    /// This is our primary graph.
    const gbwtgraph::GBWTGraph& gbwt_graph;
    
    /// We have a gapless extender to extend seed hits in haplotype space.
    GaplessExtender extender;
    
    /// We have a clusterer
    SnarlDistanceIndexClusterer clusterer;

    /// We have a distribution for read fragment lengths that takes care of
    /// knowing when we've observed enough good ones to learn a good
    /// distribution.
    FragmentLengthDistribution fragment_length_distr;
    /// We may need to complain exactly once that the distribution is bad.
    atomic_flag warned_about_bad_distribution = ATOMIC_FLAG_INIT;

//-----------------------------------------------------------------------------

    // Stages of mapping.

    /**
     * Find the minimizers in the sequence using the minimizer index, and
     * return them sorted in read order.
     */
    std::vector<Minimizer> find_minimizers(const std::string& sequence, Funnel& funnel) const;
    
    /**
     * Return the indices of all the minimizers, sorted in descending order by theit minimizers' scores.
     */
    std::vector<size_t> sort_minimizers_by_score(const std::vector<Minimizer>& minimizers) const;

    /**
     * Find seeds for all minimizers passing the filters.
     */
    std::vector<Seed> find_seeds(const VectorView<Minimizer>& minimizers, const Alignment& aln, Funnel& funnel) const;
    
    /**
     * If tracking correctness, mark seeds that are correctly mapped as correct
     * in the funnel, based on proximity along paths to the input read's
     * refpos. Otherwise, tag just as placed, with the seed's read interval.
     * Assumes we are tracking provenance.
     */
    void tag_seeds(const Alignment& aln, const std::vector<Seed>::const_iterator& begin, const std::vector<Seed>::const_iterator& end, const VectorView<Minimizer>& minimizers, size_t funnel_offset, Funnel& funnel) const;

    /**
     * Determine cluster score, read coverage, and a vector of flags for the
     * minimizers present in the cluster. Score is the sum of the scores of
     * distinct minimizers in the cluster, while read coverage is the fraction
     * of the read covered by seeds in the cluster.
     *
     * Puts the cluster in the funnel as coming from its seeds.
     */
    void score_cluster(Cluster& cluster, size_t i, const VectorView<Minimizer>& minimizers, const std::vector<Seed>& seeds, size_t seq_length) const;
    
    /**
     * Determine cluster score, read coverage, and a vector of flags for the
     * minimizers present in the cluster. Score is the sum of the scores of
     * distinct minimizers in the cluster, while read coverage is the fraction
     * of the read covered by seeds in the cluster.
     *
     * Thinks of the cluster as being made out of some fragments and
     * some new seeds from the tail end of seeds, which are already in the
     * funnel, clusters first. seed_to_fragment maps from seed to the old
     * cluster it is part of, or std::numeric_limits<size_t>::max() if it isn't
     * from an old cluster.
     *
     */
    void score_merged_cluster(Cluster& cluster, size_t i, const VectorView<Minimizer>& minimizers, const std::vector<Seed>& seeds, size_t first_new_seed, const std::vector<size_t>& seed_to_fragment, const std::vector<Cluster>& fragments, size_t seq_length, Funnel& funnel) const;
    
    /**
     * Reseed between the given graph and read positions. Produces new seeds by asking the given callback for minimizers' occurrence positions.
     *  Up to one end of the graph region can be a read end, with a pos_t matching is_empty().
     * The read region always needs to be fully defined.
     */
    std::vector<Seed> reseed_between(
        size_t read_region_start,
        size_t read_region_end,
        pos_t left_graph_pos,
        pos_t right_graph_pos,
        const HandleGraph& graph,
        const VectorView<Minimizer>& minimizers,
        const std::function<void(const Minimizer&, const std::vector<nid_t>&, const std::function<void(const pos_t&)>&)>& for_each_pos_for_source_in_subgraph) const;
    
    /// Represents configuration for chaining. May need to be derived from
    /// different class parameters depending on the chaining pass.
    struct chain_config_t {
        // Lookback config
        size_t max_lookback_bases;
        size_t min_lookback_items;
        size_t lookback_item_hard_cap;
        size_t initial_lookback_threshold;
        double lookback_scale_factor;
        double min_good_transition_score_per_base;
        
        // Item and gap scoring
        int item_bonus;
        size_t max_indel_bases;
        
        // Limits on clusters to keep
        double cluster_score_cutoff;
        bool cluster_score_cutoff_enabled;
        double cluster_coverage_threshold;
        size_t min_clusters_to_chain;
        size_t max_clusters_to_chain;
        
        // Limits on chains to compute
        size_t max_chains_per_cluster;
    };
    
    /// Represents a chaining result.
    struct chain_set_t {
        /// These are the numbers of the clusters in the order explored/the
        /// order the lists of chains appear in.
        vector<size_t> cluster_nums;
        /// These are all the chains for all the clusters, as score and sequence of visited seeds.
        /// Organized by cluster, and then best chain first.
        vector<vector<pair<int, vector<size_t>>>> cluster_chains;
        /// What cluster seeds define the space for clusters' chosen chains?
        vector<vector<size_t>> cluster_chain_seeds;
        /// Chainable anchors in the same order as seeds
        vector<algorithms::Anchor> seed_anchors;
        /// To compute the windows for explored minimizers, we need to get
        /// all the minimizers that are explored.
        SmallBitset minimizer_explored;
        /// How many hits of each minimizer ended up in each cluster we kept?
        vector<vector<size_t>> minimizer_kept_cluster_count;
        /// How many clusters were kept?
        size_t kept_cluster_count;
    };
    
    /**
     * Run chaining on some clusters. Returns the chains and the context needed to interpret them.
     */
    chain_set_t chain_clusters(const Alignment& aln, const VectorView<Minimizer>& minimizers, const std::vector<Seed>& seeds, const std::vector<Cluster>& clusters, const chain_config_t& cfg, size_t old_seed_count, size_t new_seed_start, Funnel& funnel, size_t seed_stage_offset, size_t reseed_stage_offset, LazyRNG& rng) const;
    
    /**
     * Extends the seeds in a cluster into a collection of GaplessExtension objects.
     */
    vector<GaplessExtension> extend_cluster(
        const Cluster& cluster,
        size_t cluster_num,
        const VectorView<Minimizer>& minimizers,
        const std::vector<Seed>& seeds,
        const string& sequence,
        vector<vector<size_t>>& minimizer_kept_cluster_count,
        Funnel& funnel) const;
    
    /**
     * Score the given group of gapless extensions. Determines the best score
     * that can be obtained by chaining extensions together, using the given
     * gap open and gap extend penalties to charge for either overlaps or gaps
     * in coverage of the read.
     *
     * Enforces that overlaps cannot result in containment.
     *
     * Input extended seeds must be sorted by start position.
     */
    static int score_extension_group(const Alignment& aln, const vector<GaplessExtension>& extended_seeds,
        int gap_open_penalty, int gap_extend_penalty);
    
    /**
     * Score the set of extensions for each cluster using score_extension_group().
     * Return the scores in the same order as the extension groups.
     */
    std::vector<int> score_extensions(const std::vector<std::vector<GaplessExtension>>& extensions, const Alignment& aln, Funnel& funnel) const;
    /**
     * Score the set of extensions for each cluster using score_extension_group().
     * Return the scores in the same order as the extensions.
     *
     * This version allows the collections of extensions to be scored to come
     * with annotating read numbers, which are ignored.
     */
    std::vector<int> score_extensions(const std::vector<std::pair<std::vector<GaplessExtension>, size_t>>& extensions, const Alignment& aln, Funnel& funnel) const;
    
    /**
     * Turn a chain into an Alignment.
     *
     * Operating on the given input alignment, align the tails and intervening
     * sequences along the given chain of perfect-match seeds, and return an
     * optimal Alignment.
     */
    Alignment find_chain_alignment(const Alignment& aln, const VectorView<algorithms::Anchor>& to_chain, const std::vector<size_t>& chain) const;
     
     /**
     * Operating on the given input alignment, align the tails dangling off the
     * given extended perfect-match seeds and produce an optimal alignment into
     * the given output Alignment object, best, and the second best alignment
     * into second_best.
     *
     * Uses the given RNG to break ties.
     */
    void find_optimal_tail_alignments(const Alignment& aln, const vector<GaplessExtension>& extended_seeds, LazyRNG& rng, Alignment& best, Alignment& second_best) const; 

//-----------------------------------------------------------------------------

    // Rescue.

    /**
     * Given an aligned read, extract a subgraph of the graph within a distance range
     * based on the fragment length distribution and attempt to align the unaligned
     * read to it.
     * Rescue_forward is true if the aligned read is the first and false otherwise.
     * Assumes that both reads are facing the same direction.
     * TODO: This should be const, but some of the function calls are not.
     */
    void attempt_rescue(const Alignment& aligned_read, Alignment& rescued_alignment, const VectorView<Minimizer>& minimizers, bool rescue_forward);

    /**
     * Return the all non-redundant seeds in the subgraph, including those from
     * minimizers not used for mapping.
     */
    GaplessExtender::cluster_type seeds_in_subgraph(const VectorView<Minimizer>& minimizers, const std::unordered_set<nid_t>& subgraph) const;

    /**
     * When we use dozeu for rescue, the reported alignment score is incorrect.
     * 1) Dozeu only gives the full-length bonus once.
     * 2) There is no penalty for a softclip at the edge of the subgraph.
     * This function calculates the score correctly. If the score is <= 0,
     * we realign the read using GSSW.
     * TODO: This should be unnecessary.
     */
    void fix_dozeu_score(Alignment& rescued_alignment, const HandleGraph& rescue_graph,
                         const std::vector<handle_t>& topological_order) const;

//-----------------------------------------------------------------------------

    // Helper functions.

    /**
     * Get the distance between a pair of positions, or std::numeric_limits<int64_t>::max() if unreachable.
     */
    int64_t distance_between(const pos_t& pos1, const pos_t& pos2);

    /**
     * Get the distance between a pair of read alignments, or std::numeric_limits<int64_t>::max() if unreachable.
     */
    int64_t distance_between(const Alignment& aln1, const Alignment& aln2);

    /**
     * Get the unoriented distance between a pair of positions
     */
    int64_t unoriented_distance_between(const pos_t& pos1, const pos_t& pos2) const;

    /**
     * Convert the GaplessExtension into an alignment. This assumes that the
     * extension is a full-length alignment and that the sequence field of the
     * alignment has been set.
     */
    void extension_to_alignment(const GaplessExtension& extension, Alignment& alignment) const;
    
    /**
     * Convert a WFAAlignment into a vg Alignment. This assumes that the
     * WFAAlignment is a full-length alignment and that the sequence field of
     * the vg Alignment has been set.
     */
    void wfa_alignment_to_alignment(const WFAAlignment& wfa_alignment, Alignment& alignment) const;
   
    /**
     * Clip out the part of the graph between the given positions, and dagify
     * it from the perspective of the anchors. If a left anchor is set, all
     * heads should correspond to the left anchor, and if a right anchor is
     * set, all tails should correspond to the right anchor. At least one
     * anchor must be set.
     *
     * Calls the callback with an extracted, strand-split, dagified graph, and
     * a function that translates from handle in the dagified graph to node ID
     * and orientation in the base graph.
     */
    static void with_dagified_local_graph(const pos_t& left_anchor, const pos_t& right_anchor, size_t max_path_length, const HandleGraph& graph, const std::function<void(DeletableHandleGraph&, const std::function<std::pair<nid_t, bool>(const handle_t&)>&)>& callback);
   
    /**
     * Clip out the part of the graph between the given positions and
     * global-align the sequence of the given Alignment to it. Populate the
     * Alignment's path and score.
     *
     * Finds an alignment against a graph path if it is <= max_path_length, and uses <= max_dp_cells GSSW cells.
     *
     * If one of the anchor positions is empty, does pinned alighnment against
     * the other position.
     */
    static void align_sequence_between(const pos_t& left_anchor, const pos_t& right_anchor, size_t max_path_length, const HandleGraph* graph, const GSSWAligner* aligner, Alignment& alignment, size_t max_dp_cells = std::numeric_limits<size_t>::max());
    
    /**
     * Set pair partner references for paired mapping results.
     */
    void pair_all(std::array<vector<Alignment>, 2>& mappings) const;
    
    /**
     * Add annotations to an Alignment with statistics about the minimizers.
     *
     * old_seed_count is the number of seeds in the seed vector actually
     * created at the "seed" stage of the alignment process. new_seed_offset is
     * where the first of thos eseeds appears in the funnel at the reseed stage.
     */
    void annotate_with_minimizer_statistics(Alignment& target, const VectorView<Minimizer>& minimizers, const std::vector<Seed>& seeds, size_t old_seed_count, size_t new_seed_offset, const Funnel& funnel) const;

//-----------------------------------------------------------------------------

    /**
     * Compute MAPQ caps based on all minimizers that are explored, for some definition of explored.
     *
     * Needs access to the input alignment for sequence and quality
     * information.
     *
     * Returns only an "extended" cap at the moment.
     */
    double compute_mapq_caps(const Alignment& aln, const VectorView<Minimizer>& minimizers,
                             const SmallBitset& explored);

    /**
     * Compute a bound on the Phred score probability of having created the
     * agglomerations of the specified minimizers by base errors from the given
     * sequence, which was sequenced with the given qualities.
     *
     * No limit is imposed if broken is empty.
     *
     * Takes the collection of all minimizers found, and a vector of the
     * indices of minimizers we are interested in the agglomerations of. May
     * modify the order of that index vector.
     *
     * Also takes the sequence of the read (to avoid Ns) and the quality string
     * (interpreted as a byte array).
     *
     * Currently computes a lower-score-bound, upper-probability-bound,
     * suitable for use as a mapping quality cap, by assuming the
     * easiest-to-disrupt possible layout of the windows, and the lowest
     * possible qualities for the disrupting bases.
     */
    static double window_breaking_quality(const VectorView<Minimizer>& minimizers, vector<size_t>& broken,
        const string& sequence, const string& quality_bytes);
    
    /**
     * Compute a bound on the Phred score probability of a mapping beign wrong
     * due to base errors and unlocated minimizer hits prevented us from
     * finding the true alignment.
     *  
     * Algorithm uses a "sweep line" dynamic programming approach.
     * For a read with minimizers aligned to it:
     *
     *              000000000011111111112222222222
     *              012345678901234567890123456789
     * Read:        ******************************
     * Minimizer 1:    *****
     * Minimizer 2:       *****
     * Minimizer 3:                   *****
     * Minimizer 4:                      *****
     *
     * For each distinct read interval of overlapping minimizers, e.g. in the
     * example the intervals 3,4,5; 6,7; 8,9,10; 18,19,20; 21,22; and 23,24,25
     * we consider base errors that would result in the minimizers in the
     * interval being incorrect
     *
     * We use dynamic programming sweeping left-to-right over the intervals to
     * compute the probability of the minimum number of base errors needed to
     * disrupt all the minimizers.
     *
     * Will sort minimizers_explored (which is indices into minimizers) by
     * minimizer start position.
     */
    static double faster_cap(const VectorView<Minimizer>& minimizers, vector<size_t>& minimizers_explored, const string& sequence, const string& quality_bytes);
    
    /**
     * Given a collection of minimizers, and a list of the minimizers we
     * actually care about (as indices into the collection), iterate over
     * common intervals of overlapping minimizer agglomerations.
     *   
     * Calls the given callback with (left, right, bottom, top), where left is
     * the first base of the agglomeration interval (inclusive), right is the
     * last base of the agglomeration interval (exclusive), bottom is the index
     * of the first minimizer with an agglomeration in the interval and top is
     * the index of the last minimizer with an agglomeration in the interval
     * (exclusive).
     *
     * Note that bottom and top are offsets into minimizer_indices, **NOT**
     * minimizers itself. Only contiguous ranges in minimizer_indices actually
     * make sense.
     */
    static void for_each_agglomeration_interval(const VectorView<Minimizer>& minimizers,
        const string& sequence, const string& quality_bytes,
        const vector<size_t>& minimizer_indices,
        const function<void(size_t, size_t, size_t, size_t)>& iteratee);      
    
    /**
     * Gives the log10 prob of a base error in the given interval of the read,
     * accounting for the disruption of specified minimizers.
     * 
     * minimizers is the collection of all minimizers
     *
     * disrupt_begin and disrupt_end are iterators defining a sequence of
     * **indices** of minimizers in minimizers that are disrupted.
     *
     * left and right are the inclusive and exclusive bounds of the interval
     * of the read where the disruption occurs.
     */
    static double get_log10_prob_of_disruption_in_interval(const VectorView<Minimizer>& minimizers,
        const string& sequence, const string& quality_bytes,
        const vector<size_t>::iterator& disrupt_begin, const vector<size_t>::iterator& disrupt_end,
        size_t left, size_t right);
    
    /**
     * Gives the raw probability of a base error in the given column of the
     * read, accounting for the disruption of specified minimizers.
     * 
     * minimizers is the collection of all minimizers
     *
     * disrupt_begin and disrupt_end are iterators defining a sequence of
     * **indices** of minimizers in minimizers that are disrupted.
     *
     * index is the position in the read where the disruption occurs.
     */
    static double get_prob_of_disruption_in_column(const VectorView<Minimizer>& minimizers,
        const string& sequence, const string& quality_bytes,
        const vector<size_t>::iterator& disrupt_begin, const vector<size_t>::iterator& disrupt_end,
        size_t index);
    
    /**
     * Get all the trees defining tails off the specified side of the specified
     * gapless extension. Should only be called if a tail on that side exists,
     * or this is a waste of time.
     *
     * If the gapless extension starts or ends at a node boundary, there may be
     * multiple trees produced, each with a distinct root.
     *
     * If the gapless extension abuts the edge of the read, an empty forest
     * will be produced.
     *
     * Each tree is represented as a TreeSubgraph over our gbwt_graph.
     *
     * If left_tails is true, the trees read out of the left sides of the
     * gapless extension. Otherwise they read out of the right side.
     *
     * As a side effect, saves the length of the longest detectable gap in an
     * alignment of a tail to the forest into the provided location, if set.
     */
    vector<TreeSubgraph> get_tail_forest(const GaplessExtension& extended_seed,
        size_t read_length, bool left_tails, size_t* longest_detectable_gap = nullptr) const;
        
    /**
     * Find the best alignment of the given sequence against any of the trees
     * provided in trees, where each tree is a TreeSubgraph over the GBWT
     * graph. Each tree subgraph is rooted at the left in its own local
     * coordinate space, even if we are pinning on the right.
     *
     * If no mapping is possible (for example, because there are no trees),
     * produce a pure insert at default_position.
     *
     * Alignment is always pinned.
     *
     * If pin_left is true, pin the alignment on the left to the root of each
     * tree. Otherwise pin it on the right to the root of each tree.
     *
     * Limits the length of the longest gap to longest_detectable_gap.
     *
     * Returns alignments in gbwt_graph space.
     */
    pair<Path, size_t> get_best_alignment_against_any_tree(const vector<TreeSubgraph>& trees, const string& sequence,
        const Position& default_position, bool pin_left, size_t longest_detectable_gap, LazyRNG& rng) const;
        
    /// We define a type for shared-tail lists of Mappings, to avoid constantly
    /// copying Path objects.
    using ImmutablePath = structures::ImmutableList<Mapping>;
    
    /**
     * Get the from length of an ImmutabelPath.
     *
     * Can't be called path_from_length or it will shadow the one for Paths
     * instead of overloading.
     */
    static size_t immutable_path_from_length(const ImmutablePath& path);
    
    /**
     * Convert an ImmutablePath to a Path.
     */
    static Path to_path(const ImmutablePath& path);

    /**
     * Run a DFS on valid haplotypes in the GBWT starting from the given
     * Position, and continuing up to the given number of bases.
     *
     * Calls enter_handle when the DFS enters a haplotype visit to a particular
     * handle, and exit_handle when it exits a visit. These let the caller
     * maintain a stack and track the traversals.
     *
     * The starting node is only entered if its offset isn't equal to its
     * length (i.e. bases remain to be visited).
     *
     * Stopping early is not permitted.
     */
    void dfs_gbwt(const Position& from, size_t walk_distance,
        const function<void(const handle_t&)>& enter_handle, const function<void(void)> exit_handle) const;
     
    /**
     * The same as dfs_gbwt on a Position, but takes a handle in the
     * backing gbwt_graph and an offset from the start of the handle instead.
     */ 
    void dfs_gbwt(handle_t from_handle, size_t from_offset, size_t walk_distance,
        const function<void(const handle_t&)>& enter_handle, const function<void(void)> exit_handle) const;
        
    /**
     * The same as dfs_gbwt on a handle and an offset, but takes a
     * gbwt::SearchState that defines only some haplotypes on a handle to start
     * with.
     */ 
    void dfs_gbwt(const gbwt::SearchState& start_state, size_t from_offset, size_t walk_distance,
        const function<void(const handle_t&)>& enter_handle, const function<void(void)> exit_handle) const;
 
    /**
     * Score a pair of alignments given the distance between them
     */
    double score_alignment_pair(Alignment& aln1, Alignment& aln2, int64_t fragment_distance);
    
    /**
     * Given a count of items, a function to get the score of each, a
     * score-difference-from-the-best cutoff, a min and max processed item
     * count, and a function to get a sort-shuffling seed for breaking ties,
     * process items in descending score order by calling process_item with the
     * item's number, until min_count items are processed and either max_count
     * items are processed or the score difference threshold is hit (or we run
     * out of items).
     *
     * If process_item returns false, the item is skipped and does not count
     * against min_count or max_count.
     *
     * Call discard_item_by_count with the item's number for all remaining
     * items that would pass the score threshold.
     *
     * Call discard_item_by_score with the item's number for all remaining
     * items that would fail the score threshold.
     */
    template<typename Score = double>
    void process_until_threshold_a(size_t items, const function<Score(size_t)>& get_score,
        double threshold, size_t min_count, size_t max_count,
        LazyRNG& rng,
        const function<bool(size_t)>& process_item,
        const function<void(size_t)>& discard_item_by_count,
        const function<void(size_t)>& discard_item_by_score) const;
     
    /**
     * Same as the other process_until_threshold functions, except using a vector to supply scores.
     */
    template<typename Score = double>
    void process_until_threshold_b(const vector<Score>& scores,
        double threshold, size_t min_count, size_t max_count,
        LazyRNG& rng,
        const function<bool(size_t)>& process_item,
        const function<void(size_t)>& discard_item_by_count,
        const function<void(size_t)>& discard_item_by_score) const;
     
    /**
     * Same as the other process_until_threshold functions, except user supplies 
     * comparator to sort the items (must still be sorted by score).
     */
    template<typename Score = double>
    void process_until_threshold_c(size_t items, const function<Score(size_t)>& get_score,
        const function<bool(size_t, size_t)>& comparator,
        double threshold, size_t min_count, size_t max_count,
        LazyRNG& get_seed,
        const function<bool(size_t)>& process_item,
        const function<void(size_t)>& discard_item_by_count,
        const function<void(size_t)>& discard_item_by_score) const;
        
    // Internal debugging functions
    
    /// Get the thread identifier prefix for logging
    static string log_name();
    
    /// Turn an Alignment into a conveniently-sized string for logging
    static string log_alignment(const Alignment& aln);
    
    /// Turn an Path from an alignment into a conveniently-sized string for logging
    static string log_alignment(const Path& path, bool force_condensed = false);
    
    /// Turn a list of bit flags into a compact representation.
    static string log_bits(const std::vector<bool>& bits);
    
    /// Dump a whole chaining problem
    static void dump_chaining_problem(const std::vector<algorithms::Anchor>& anchors, const std::vector<size_t>& cluster_seeds_sorted, const HandleGraph& graph);
    
    /// Dump all the given minimizers, with optional subset restriction
    static void dump_debug_minimizers(const VectorView<Minimizer>& minimizers, const string& sequence,
                                      const vector<size_t>* to_include = nullptr, size_t start_offset = 0, size_t length_limit = std::numeric_limits<size_t>::max());
    
    /// Dump all the extansions in an extension set
    static void dump_debug_extension_set(const HandleGraph& graph, const Alignment& aln, const vector<GaplessExtension>& extended_seeds);
    
    /// Print a sequence with base numbering
    static void dump_debug_sequence(ostream& out, const string& sequence, size_t start_offset = 0, size_t length_limit = std::numeric_limits<size_t>::max());
    
    /// Print the seed content of a cluster.
    static void dump_debug_clustering(const Cluster& cluster, size_t cluster_number, const VectorView<Minimizer>& minimizers, const std::vector<Seed>& seeds);

    /// Do a brute check of the clusters. Print errors to stderr
    bool validate_clusters(const std::vector<std::vector<Cluster>>& clusters, const std::vector<std::vector<Seed>>& seeds, size_t read_limit, size_t fragment_limit) const;
    
    /// Print information about a selected set of seeds.
    static void dump_debug_seeds(const VectorView<Minimizer>& minimizers, const std::vector<Seed>& seeds, const std::vector<size_t>& selected_seeds);
    
    /// Print information about a read to be aligned
    static void dump_debug_query(const Alignment& aln);
    
    /// Print information about a read pair to be aligned
    static void dump_debug_query(const Alignment& aln1, const Alignment& aln2);
    
    /// Length at which we cut over to long-alignment logging.
    const static size_t LONG_LIMIT = 256;
    
    /// Count at which we cut over to summary logging.
    const static size_t MANY_LIMIT = 20;


    friend class TestMinimizerMapper;
};

template<typename Score>
void MinimizerMapper::process_until_threshold_a(size_t items, const function<Score(size_t)>& get_score,
    double threshold, size_t min_count, size_t max_count,
    LazyRNG& rng,
    const function<bool(size_t)>& process_item,
    const function<void(size_t)>& discard_item_by_count,
    const function<void(size_t)>& discard_item_by_score) const {

    process_until_threshold_c<Score>(items, get_score, [&](size_t a, size_t b) -> bool {
        return (get_score(a) > get_score(b));
    },threshold, min_count, max_count, rng, process_item, discard_item_by_count, discard_item_by_score);
}

template<typename Score>
void MinimizerMapper::process_until_threshold_b(const vector<Score>& scores,
    double threshold, size_t min_count, size_t max_count,
    LazyRNG& rng,
    const function<bool(size_t)>& process_item,
    const function<void(size_t)>& discard_item_by_count,
    const function<void(size_t)>& discard_item_by_score) const {
    
    process_until_threshold_c<Score>(scores.size(), [&](size_t i) -> Score {
        return scores[i];
    }, [&](size_t a, size_t b) -> bool {
        return (scores[a] > scores[b]);
    },threshold, min_count, max_count, rng, process_item, discard_item_by_count, discard_item_by_score);
}

template<typename Score>
void MinimizerMapper::process_until_threshold_c(size_t items, const function<Score(size_t)>& get_score,
        const function<bool(size_t, size_t)>& comparator,
        double threshold, size_t min_count, size_t max_count,
        LazyRNG& rng,
        const function<bool(size_t)>& process_item,
        const function<void(size_t)>& discard_item_by_count,
        const function<void(size_t)>& discard_item_by_score) const {

    // Sort item indexes by item score
    vector<size_t> indexes_in_order;
    indexes_in_order.reserve(items);
    for (size_t i = 0; i < items; i++) {
        indexes_in_order.push_back(i);
    }
    
    // Put the highest scores first, but shuffle top ties so reads spray evenly
    // across equally good mappings
    sort_shuffling_ties(indexes_in_order.begin(), indexes_in_order.end(), comparator, rng);

    // Retain items only if their score is at least as good as this
    double cutoff = items == 0 ? 0 : get_score(indexes_in_order[0]) - threshold;
    
    // Count up non-skipped items for min_count and max_count
    size_t unskipped = 0;
    
    // Go through the items in descending score order.
    for (size_t i = 0; i < indexes_in_order.size(); i++) {
        // Find the item we are talking about
        size_t& item_num = indexes_in_order[i];
        
        if (threshold != 0 && get_score(item_num) <= cutoff) {
            // Item would fail the score threshold
            
            if (unskipped < min_count) {
                // But we need it to make up the minimum number.
                
                // Go do it.
                // If it is not skipped by the user, add it to the total number
                // of unskipped items, for min/max number accounting.
                unskipped += (size_t) process_item(item_num);
            } else {
                // We will reject it for score
                discard_item_by_score(item_num);
            }
        } else {
            // The item has a good enough score
            
            if (unskipped < max_count) {
                // We have room for it, so accept it.
                
                // Go do it.
                // If it is not skipped by the user, add it to the total number
                // of unskipped items, for min/max number accounting.
                unskipped += (size_t) process_item(item_num);
            } else {
                // We are out of room! Reject for count.
                discard_item_by_count(item_num);
            }
        }
    }
}

}



#endif
