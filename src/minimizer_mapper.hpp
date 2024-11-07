#ifndef VG_MINIMIZER_MAPPER_HPP_INCLUDED
#define VG_MINIMIZER_MAPPER_HPP_INCLUDED

/** 
 * \file minimizer_mapper.hpp
 * Defines a mapper that uses the minimizer index and GBWT-based extension.
 */

#include "algorithms/chain_items.hpp"
#include "algorithms/nearest_offsets_in_paths.hpp"
#include "algorithms/pad_band.hpp"
#include "aligner.hpp"
#include "vg/io/alignment_emitter.hpp"
#include "gbwt_extender.hpp"
#include "snarl_seed_clusterer.hpp"
#include "zip_code_tree.hpp"
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
         const ZipCodeCollection* zipcodes,
         const PathPositionHandleGraph* path_graph = nullptr);

    using AlignerClient::set_alignment_scores;
    virtual void set_alignment_scores(const int8_t* score_matrix, int8_t gap_open, int8_t gap_extend, int8_t full_length_bonus);

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

    /// Window count for minimizer downsampling
    static constexpr size_t default_minimizer_downsampling_window_count = 0;
    size_t minimizer_downsampling_window_count = default_minimizer_downsampling_window_count;

    static constexpr size_t default_minimizer_downsampling_max_window_length = std::numeric_limits<size_t>::max();
    size_t minimizer_downsampling_max_window_length = default_minimizer_downsampling_max_window_length;

    //We allow additional seeds past the maximum number of seeds allowed if they cover a region of the read that
    //was not covered by accepted seeds.
    //The coverage of a seed is its sequence plus the seed_coverage_flank on either end
    static constexpr size_t default_minimizer_coverage_flank = 250;
    size_t minimizer_coverage_flank = default_minimizer_coverage_flank;


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

    // If a cluster's score is smaller than the best score of any cluster by more than
    /// this much, then don't extend it
    static constexpr double default_cluster_score_threshold = 50;
    double cluster_score_threshold = default_cluster_score_threshold;
    
    /// If the second best cluster's score is no more than this many points below
    /// the cutoff set by cluster_score_threshold, snap that cutoff down to the
    /// second best cluster's score, to avoid throwing away promising
    /// secondaries.
    static constexpr double default_pad_cluster_score_threshold = 20;
    double pad_cluster_score_threshold = default_pad_cluster_score_threshold;

    /// If the read coverage of a cluster is less than the best coverage of any tree
    /// by more than this much, don't extend it
    static constexpr double default_cluster_coverage_threshold = 0.3;
    double cluster_coverage_threshold = default_cluster_coverage_threshold;

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

    /// How many extensions should we try as seeds within a mapping location?
    static constexpr size_t default_max_local_extensions = numeric_limits<size_t>::max();
    size_t max_local_extensions = default_max_local_extensions;

    
    /////////////////
    // More shared parameters:
    /////////////////
    
    /// How many alignments should we make, max?
    static constexpr size_t default_max_alignments = 8;
    size_t max_alignments = default_max_alignments;

    /// How many mismatches should we allow in gapless extension (except for
    /// start node where the limit doesn't count)?
    static constexpr size_t default_max_extension_mismatches = GaplessExtender::MAX_MISMATCHES;
    size_t max_extension_mismatches = default_max_extension_mismatches;
    
    //////////////////
    // Alignment-from-chains/long read Giraffe specific parameters:
    //////////////////
    
    /// If true, produce alignments from extension sets by chaining gapless
    /// extensions up and aligning the sequences between them. If false,
    /// produce alignments by aligning the tails off of individual gapless
    /// extensions.
    static constexpr bool default_align_from_chains = false;
    bool align_from_chains = default_align_from_chains;

    /// When making zipcode trees, at what multiple of the read length should the trees
    /// be split?
    static constexpr double default_zipcode_tree_scale = 2.0;
    double zipcode_tree_scale = default_zipcode_tree_scale;

    /// How far do we want to go down looking at zip code trees to make fragments?
    static constexpr double default_zipcode_tree_score_threshold = 50;
    double zipcode_tree_score_threshold = default_zipcode_tree_score_threshold;

    /// If the second best tree's score is no more than this many points below
    /// the cutoff set by zipcode_tree_score_threshold, snap that cutoff down
    /// to the second best tree's score, to avoid throwing away promising
    /// secondaries.
    static constexpr double default_pad_zipcode_tree_score_threshold = 20;
    double pad_zipcode_tree_score_threshold = default_pad_zipcode_tree_score_threshold;

    /// If the read coverage of a tree is less than the best coverage of any tree
    /// by more than this much, don't extend it
    static constexpr double default_zipcode_tree_coverage_threshold = 0.3;
    double zipcode_tree_coverage_threshold = default_zipcode_tree_coverage_threshold;

    /// How many things should we produce fragments for, min?
    static constexpr size_t default_min_to_fragment = 4;
    size_t min_to_fragment = default_min_to_fragment;

    /// How many things should we produce fragments for, max?
    static constexpr size_t default_max_to_fragment = 10;
    size_t max_to_fragment = default_max_to_fragment;
    
    /// Do gapless extension to the seeds in each tree before fragmenting the tree if the 
    /// read length is less than the limit.
    static constexpr size_t default_gapless_extension_limit = 0;
    size_t gapless_extension_limit = default_gapless_extension_limit;

    /// How many bases should we look back when making fragments?
    static constexpr size_t default_fragment_max_lookback_bases = 300;
    size_t fragment_max_lookback_bases = default_fragment_max_lookback_bases;
    /// How many bases should we look back when making fragments, per base of read length?
    static constexpr double default_fragment_max_lookback_bases_per_base = 0.03;
    double fragment_max_lookback_bases_per_base = default_fragment_max_lookback_bases_per_base;
    /// How many fragments should we try and make when fragmenting something?
    static constexpr size_t default_max_fragments = std::numeric_limits<size_t>::max();
    size_t max_fragments = default_max_fragments;
    
    /// How much of a multiple should we apply to each transition's gap penalty
    /// at fragmenting?
    static constexpr double default_fragment_gap_scale = 1.0;
    double fragment_gap_scale = default_fragment_gap_scale;
    // How many points should we treat a non-gap connection base as producing, at fragmenting?
    static constexpr double default_fragment_points_per_possible_match = 0;
    double fragment_points_per_possible_match = default_fragment_points_per_possible_match;
    /// How many bases of indel should we allow in fragments?
    static constexpr size_t default_fragment_max_indel_bases = 2000;
    size_t fragment_max_indel_bases = default_fragment_max_indel_bases;
    /// How many bases of indel should we allow in fragments per base of read length?
    static constexpr double default_fragment_max_indel_bases_per_base = 0.2;
    double fragment_max_indel_bases_per_base = default_fragment_max_indel_bases_per_base;
    
    /// When converting chains to alignments, what's the longest gap between
    /// items we will try to WFA align? Passing strings longer than ~100bp
    /// can cause WFAAligner to run for a pathologically long amount of time.
    /// May not be 0.
    static constexpr size_t default_max_chain_connection = 100;
    size_t max_chain_connection = default_max_chain_connection;
    /// Similarly, what is the maximum tail length we will try to WFA align?
    static constexpr size_t default_max_tail_length = 100;
    size_t max_tail_length = default_max_tail_length;
    
    /// How good should a fragment be in order to keep it? Fragments with
    /// scores less than this fraction of the best fragment's score
    /// will not be used.
    static constexpr double default_fragment_score_fraction = 0.1;
    double fragment_score_fraction = default_fragment_score_fraction;
    
    /// How high should we get the score threshold based on the best fragment's score get?
    static constexpr double default_fragment_max_min_score = std::numeric_limits<double>::max();
    double fragment_max_min_score = default_fragment_max_min_score;

    /// What minimum score in points should a fragment have in order to keep
    /// it? Needs to be set to some kind of significance threshold.
    static constexpr double default_fragment_min_score = 60;
    double fragment_min_score = default_fragment_min_score;

    /// If a fragment set's score is smaller than the best 
    /// fragment set's score by more than this much, don't align it
    static constexpr double default_fragment_set_score_threshold = 0;
    double fragment_set_score_threshold = default_fragment_set_score_threshold;

    /// Disregard the fragment set score thresholds when they would give us
    /// fewer than this many chainign problems done.
    static constexpr int default_min_chaining_problems = 1;
    int min_chaining_problems = default_min_chaining_problems;
    
    /// Do no more than this many chaining problems.
    static constexpr int default_max_chaining_problems = std::numeric_limits<int>::max();
    int max_chaining_problems = default_max_chaining_problems;

    /// Sometimes we don't do chaining but instead turn fragments directly into chains
    /// If this is 0, then do chaining. Otherwise take up to this many fragments and turn them into chains
    static constexpr size_t default_max_direct_to_chain = 0;
    size_t max_direct_to_chain = default_max_direct_to_chain;

    /// How many bases should we look back when chaining?
    static constexpr size_t default_max_lookback_bases = 3000;
    size_t max_lookback_bases = default_max_lookback_bases;
    /// How many bases should we look back when chaining, per base of read length?
    static constexpr double default_max_lookback_bases_per_base = 0.3;
    double max_lookback_bases_per_base = default_max_lookback_bases_per_base;

    /// How much of a bonus should we give to each item in
    /// fragmenting/chaining?
    static constexpr int default_item_bonus = 0;
    int item_bonus = default_item_bonus;
    /// How much of a multiple should we apply to each item's non-bonus score
    /// in fragmenting/chaining?
    static constexpr double default_item_scale = 1.0;
    double item_scale = default_item_scale;
    /// How much of a multiple should we apply to each transition's gap penalty
    /// at chaining?
    static constexpr double default_gap_scale = 1.0;
    double gap_scale = default_gap_scale;
    // How many points should we treat a non-gap connection base as producing, at chaining?
    static constexpr double default_points_per_possible_match = 0;
    double points_per_possible_match = default_points_per_possible_match;
    /// How many bases of indel should we allow in chaining?
    static constexpr size_t default_max_indel_bases = 2000;
    size_t max_indel_bases = default_max_indel_bases;
    /// How many bases of indel should we allow in chaining, per base of read length?
    static constexpr double default_max_indel_bases_per_base = 0.2;
    double max_indel_bases_per_base = default_max_indel_bases_per_base;
    
    /// If a chain's score is smaller than the best 
    /// chain's score by more than this much, don't align it
    static constexpr double default_chain_score_threshold = 100;
    double chain_score_threshold = default_chain_score_threshold;
    
    /// Disregard the chain score thresholds when they would give us
    /// fewer than this many chains aligned.
    static constexpr int default_min_chains = 4;
    int min_chains = default_min_chains;

    /// Allow up to this many chains per tree
    static constexpr size_t default_max_chains_per_tree = 1;
    size_t max_chains_per_tree = default_max_chains_per_tree;
    
    /// Even if we would have fewer than min_chains results, don't
    /// process anything with a score smaller than this, per read base.
    static constexpr double default_min_chain_score_per_base = 0.01;
    double min_chain_score_per_base = default_min_chain_score_per_base;

    /// Limit the min chain score to no more than this.
    static constexpr int default_max_min_chain_score = 200;
    int max_min_chain_score = default_max_min_chain_score;
    
    /// How long of a DP can we do before Dozeu gets lost at traceback due to
    /// 16-bit score overflow?
    static constexpr size_t default_max_tail_dp_length = 30000;
    size_t max_tail_dp_length = default_max_tail_dp_length;
    /// How long of a DP can we do before something might go wrong with BandedGlobalAligner or the GBWT-based WFA?
    static constexpr size_t default_max_middle_dp_length = std::numeric_limits<int32_t>::max();
    size_t max_middle_dp_length = default_max_middle_dp_length;
    
    /// How many DP cells should we be willing to do for an end-pinned
    /// alignment? If we want to do more than this, just leave tail unaligned.
    static constexpr size_t default_max_dp_cells = std::numeric_limits<size_t>::max();
    size_t max_dp_cells = default_max_dp_cells;

    /// How many gap bases should we allow in a Dozeu tail alignment, max?
    static constexpr size_t default_max_tail_gap = std::numeric_limits<size_t>::max();
    size_t max_tail_gap = default_max_tail_gap;

    /// How many gap bases should we allow in a between-seed alignment, max?
    static constexpr size_t default_max_middle_gap = std::numeric_limits<size_t>::max();
    size_t max_middle_gap = default_max_middle_gap;
    
    /// How many mismatch bases (or equivalent score of indels) should we allow in WFA connections and tails?
    static constexpr int default_wfa_max_mismatches = 2;
    int wfa_max_mismatches = default_wfa_max_mismatches;
    /// How many mismatch bases (or equivalent score of indels) should we allow in WFA connections and tails per base of read sequence?
    static constexpr double default_wfa_max_mismatches_per_base= 0.1;
    double wfa_max_mismatches_per_base = default_wfa_max_mismatches_per_base;
    /// How many mismatch bases (or equivalent score of indels) should we allow in WFA connections and tails maximum, at any read length?
    static constexpr int default_wfa_max_max_mismatches = 20;
    int wfa_max_max_mismatches = default_wfa_max_max_mismatches;

    /// How far behind the leader should the WFA be allowed to get?
    static constexpr int default_wfa_distance = WFAExtender::ErrorModel::default_distance().min;
    int wfa_distance = default_wfa_distance;
    /// How far behind the leader should the WFA be allowed to get, per base of read sequence?
    static constexpr double default_wfa_distance_per_base = WFAExtender::ErrorModel::default_distance().per_base;
    double wfa_distance_per_base = default_wfa_distance_per_base;
    /// How far behind the leader should the WFA be allowed to get, at any read length?
    static constexpr int default_wfa_max_distance = WFAExtender::ErrorModel::default_distance().max;
    int wfa_max_distance = default_wfa_max_distance;

    /// Should alignments be ranked by chain score instead of base-level score?
    static constexpr bool default_sort_by_chain_score = false;
    bool sort_by_chain_score = default_sort_by_chain_score;

    /// How much of an alignment needs to be from distinct nodes to be a distinct alignment?
    static constexpr double default_min_unique_node_fraction = 0.0;
    double min_unique_node_fraction = default_min_unique_node_fraction;

    /// If set, cap mapping quality based on minimizer layout in the read. Only
    /// really likely to help for short reads.
    static constexpr bool default_use_explored_cap = false;
    bool use_explored_cap = default_use_explored_cap;
    /// What number of bp should we re-scale scores to for MAPQ, for calibration? 0 for off.
    static constexpr size_t default_mapq_score_window = 0;
    size_t mapq_score_window = default_mapq_score_window;
    /// How should we scale scores before mapq, for calibration
    static constexpr double default_mapq_score_scale = 1.0;
    double mapq_score_scale = default_mapq_score_scale;

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

    /// Set refpos field of alignments to positions on nodes they visit.
    static constexpr bool default_set_refpos = false;
    bool set_refpos = default_set_refpos;
    
    /// Track which internal work items came from which others during each
    /// stage of the mapping algorithm.
    static constexpr bool default_track_provenance = false;
    bool track_provenance = default_track_provenance;

    /// Guess which seed hits are correct by location in the linear reference
    /// and track if/when their descendants make it through stages of the
    /// algorithm. Only works if track_provenance is true.
    static constexpr bool default_track_correctness = false;
    bool track_correctness = default_track_correctness;

    /// Track linear reference position for placements in log output.
    static constexpr bool default_track_position = false;
    bool track_position = default_track_position;
    
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
        const typename gbwtgraph::DefaultMinimizerIndex::value_type* occs;
        int32_t length; // How long is the minimizer (index's k)
        int32_t candidates_per_window; // How many minimizers compete to be the best (index's w), or 1 for syncmers.  
        double score; // Scores as 1 + ln(hard_hit_cap) - ln(hits).
        bool is_repetitive; //Is this minimizer in a repetitive region of the read based on its neighbors

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

        /// Get the position on the read's sequence that corresponds to the
        /// located graph positions. For reverse-strand minimizers this will be
        /// at the end of the minimizer's interval in the read.
        inline size_t pin_offset() const {
            return this->value.offset;
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
    inline static gbwtgraph::Payload no_chain_info() {
        return MIPayload::NO_CODE;  
    } 
    
    /// How do we convert chain info to an actual seed of the type we are using?
    /// Also needs to know the hit position, and the minimizer number.
    inline static Seed chain_info_to_seed(const pos_t& hit, size_t minimizer, const ZipCode& zip) {
        return { hit, minimizer, zip};
    }
    
    /// Convert a collection of seeds to a collection of chaining anchors.
    std::vector<algorithms::Anchor> to_anchors(const Alignment& aln, const VectorView<Minimizer>& minimizers, std::vector<Seed>& seeds) const;
    
    /// Convert a single seed to a single chaining anchor.
    static algorithms::Anchor to_anchor(const Alignment& aln, const VectorView<Minimizer>& minimizers, std::vector<Seed>& seeds, size_t seed_number, const HandleGraph& graph, const Aligner* aligner);

    /// Convert a read region, and the seeds that that region covers the
    /// stapled bases of (sorted by stapled base), into a single chaining
    /// anchor. Takes an iterator range of positions within the base range that
    /// are mismatches.
    static algorithms::Anchor to_anchor(const Alignment& aln, size_t read_start, size_t read_end, const std::vector<size_t>& sorted_seeds, const std::vector<algorithms::Anchor>& seed_anchors, const std::vector<size_t>::const_iterator& mismatch_begin, const std::vector<size_t>::const_iterator& mismatch_end, const HandleGraph& graph, const Aligner* aligner);
    
    /// Convert an Anchor to a WFAAlignment, given the input read it is from and the Aligner to use for scoring.
    /// Accounts for fuill length bonuses if the anchor abuts the end of the read.
    WFAAlignment to_wfa_alignment(const algorithms::Anchor& anchor, const Alignment& aln, const Aligner* aligner) const; 

    /// The information we store for each cluster.
    typedef SnarlDistanceIndexClusterer::Cluster Cluster;

    // These are our indexes
    const PathPositionHandleGraph* path_graph; // Can be nullptr; only needed for correctness or position tracking.
    const gbwtgraph::DefaultMinimizerIndex& minimizer_index;
    SnarlDistanceIndex* distance_index;
    const ZipCodeCollection* zipcodes;
    /// This is our primary graph.
    const gbwtgraph::GBWTGraph& gbwt_graph;
    
    /// We have a gapless extender to extend seed hits in haplotype space.
    /// Because this needs a reference to an Aligner, and because changing the
    /// scoring parameters deletes all the alignmers, we need to keep this
    /// somewhere we can clear out.
    std::unique_ptr<GaplessExtender> extender;
    
    /// We have a clusterer
    SnarlDistanceIndexClusterer clusterer;

    /// We have a zip code tree for finding distances between seeds 
    ZipCodeForest zip_forest;

    /// We have a function for determinign band paddding for banded alignment
    /// when aligning from chains.
    std::function<size_t(const Alignment&, const HandleGraph&)> choose_band_padding;

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
     * Flag minimizers as being in repetitive regions of the read
     */
    void flag_repetitive_minimizers(std::vector<Minimizer>& minimizers_in_read_order) const;
    
    /**
     * Return the indices of all the minimizers, sorted in descending order by their minimizers' scores.
     */
    std::vector<size_t> sort_minimizers_by_score(const std::vector<Minimizer>& minimizers_in_read_order, LazyRNG& rng) const;

    /**
     * Find seeds for all minimizers passing the filters. Takes in minimizers
     * sorted in read order, and a view of them sorted in score order.
     */
    std::vector<Seed> find_seeds(const std::vector<Minimizer>& minimizers_in_read_order, const VectorView<Minimizer>& minimizers, const Alignment& aln, Funnel& funnel) const;
    
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
    void score_cluster(Cluster& cluster, size_t i, const VectorView<Minimizer>& minimizers, const std::vector<Seed>& seeds, size_t seq_length, Funnel& funnel) const;

    /**
     * Determine score and read coverage for a zip code tree. Score is the sum
     * of the scores of distinct minimizers in the tree, while read coverage is
     * the fraction of the read covered by seeds in the tree.
     *
     * Puts the tree in the funnel as coming from its seeds.
     */
    std::pair<double, double> score_tree(const ZipCodeForest& zip_code_forest, size_t i, const VectorView<Minimizer>& minimizers, const std::vector<Seed>& seeds, size_t seq_length, Funnel& funnel) const;
    
    /**
     * Extends the seeds in a cluster or other grouping into a collection of
     * GaplessExtension objects.
     *
     * If funnel is set, the group is intended to come from the previous funnel
     * stage and will be introduced in this one.
     *
     * If seeds_used is not null, it should be an empty vector that gets filled
     * with, for each gapless extension, the numbers of the seeds in seeds that
     * are subsumed into the extension. They will be sorted by the stapled base
     * (first base for forward strand, last base for reverse strand) in the
     * read.
     *
     * Note that multiple gapless extensions might cover each seed position or
     * use each seed.
     */
    vector<GaplessExtension> extend_seed_group(
        const std::vector<size_t>& seed_group,
        size_t source_num,
        const VectorView<Minimizer>& minimizers,
        const std::vector<Seed>& seeds,
        const string& sequence,
        size_t max_mismatches,
        vector<vector<size_t>>* minimizer_kept_count = nullptr,
        Funnel* funnel = nullptr,
        std::vector<std::vector<size_t>>* seeds_used = nullptr) const;
    
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
     * Get the fraction of read bases covered by the given chains/fragments of
     * seeds. A base is covered if it is between the first and last endpoints
     * in the read of any of the given lists of seeds. The lists of seeds are
     * each assumed to be colinear in the read.
     */
    double get_read_coverage(const Alignment& aln, const VectorView<std::vector<size_t>>& seed_sets, const std::vector<Seed>& seeds, const VectorView<Minimizer>& minimizers) const;
    
    /// Struct to represent per-DP-method stats. 
    struct aligner_stats_t {

        /// Collection of values you can +=
        struct stat_collection_t {
            std::vector<double> values;
            inline stat_collection_t& operator+=(const double& value) {
                values.push_back(value);
                return *this;
            }
            inline stat_collection_t& operator+=(const stat_collection_t& other) {
                std::copy(other.values.begin(), other.values.end(), std::back_inserter(values));
                return *this;
            }

            inline double total() const {
                return std::accumulate(values.begin(), values.end(), 0.0);
            }
        };
        
        /// Struct to represent counts of bases or seconds or invocations used by different aligners.
        struct stat_set_t {
            stat_collection_t wfa_tail;
            stat_collection_t wfa_middle;
            stat_collection_t dozeu_tail;
            stat_collection_t bga_middle;

            inline stat_set_t& operator+=(const stat_set_t& other) {
                this->wfa_tail += other.wfa_tail;
                this->wfa_middle += other.wfa_middle;
                this->dozeu_tail += other.dozeu_tail;
                this->bga_middle += other.bga_middle;

                return *this;
            }

            inline void add_annotations(Alignment& aln, const std::string& scope, const std::string& type) {
                set_annotation(aln, "aligner_stats.per_" + scope + ".tail." + type + ".wfa", wfa_tail.total());
                set_annotation(aln, "aligner_stats.per_" + scope + ".tail." + type + ".wfa_values", wfa_tail.values);
                set_annotation(aln, "aligner_stats.per_" + scope + ".tail." + type + ".dozeu", dozeu_tail.total());
                set_annotation(aln, "aligner_stats.per_" + scope + ".tail." + type + ".dozeu_values", dozeu_tail.values);
                set_annotation(aln, "aligner_stats.per_" + scope + ".tail." + type + ".total", wfa_tail.total() + dozeu_tail.total());

                set_annotation(aln, "aligner_stats.per_" + scope + ".middle." + type + ".wfa", wfa_middle.total());
                set_annotation(aln, "aligner_stats.per_" + scope + ".middle." + type + ".wfa_values", wfa_middle.values);
                set_annotation(aln, "aligner_stats.per_" + scope + ".middle." + type + ".bga", bga_middle.total());
                set_annotation(aln, "aligner_stats.per_" + scope + ".middle." + type + ".bga_values", bga_middle.values);
                set_annotation(aln, "aligner_stats.per_" + scope + ".middle." + type + ".total", wfa_middle.total() + bga_middle.total());
            }
        };

        stat_set_t bases;
        stat_set_t time;
        stat_set_t invocations;
        stat_set_t fallbacks;

        inline aligner_stats_t& operator+=(const aligner_stats_t& other) {
            this->bases += other.bases;
            this->time += other.time;
            this->invocations += other.invocations;
            this->fallbacks += other.fallbacks;

            return *this;
        }

        inline void add_annotations(Alignment& aln, const std::string& scope) {
            bases.add_annotations(aln, scope, "bases");
            time.add_annotations(aln, scope, "time");
            invocations.add_annotations(aln, scope, "invocations");
            fallbacks.add_annotations(aln, scope, "fallbacks");
        }
    };

    /**
     * Given a collection of zipcode trees, score the trees and do fragmenting on the best trees.
     * 
     * This will fill in the given vectors of fragments, fragment scores, etc.
     *
     * If we do gapless extension, turn good full-length gapless extensions into alignments and return them in alignments
     * Gapless extensions are considered good enough if they have fewer than default_max_extension_mismatches mismatches
     */
    void do_fragmenting_on_trees(Alignment& aln, const ZipCodeForest& zip_code_forest, const std::vector<Seed>& seeds, const VectorView<MinimizerMapper::Minimizer>& minimizers,
                                  const vector<algorithms::Anchor>& seed_anchors,
                                  std::vector<std::vector<size_t>>& fragments, std::vector<double>& fragment_scores,
                                  std::vector<algorithms::Anchor>& fragment_anchors, std::vector<size_t>& fragment_source_tree,
                                  std::vector<std::vector<size_t>>& minimizer_kept_fragment_count, std::vector<double>& multiplicity_by_fragment,
                                  std::vector<Alignment>& alignments, SmallBitset& minimizer_explored, vector<double>& multiplicity_by_alignment,
                                  LazyRNG& rng, Funnel& funnel) const;
    
    /**
     * Given a collection of fragments, filter down to the good ones and do chaining on them
     */
    void do_chaining_on_fragments(Alignment& aln, const ZipCodeForest& zip_code_forest, const std::vector<Seed>& seeds, const VectorView<MinimizerMapper::Minimizer>& minimizers, 
                                  const std::vector<std::vector<size_t>>& fragments, const std::vector<double>& fragment_scores, 
                                  const std::vector<algorithms::Anchor>& fragment_anchors, const std::vector<size_t>& fragment_source_tree,
                                  const std::vector<std::vector<size_t>>& minimizer_kept_fragment_count, const std::vector<double>& multiplicity_by_fragment,
                                  std::vector<std::vector<size_t>>& chains, std::vector<size_t>& chain_source_tree, 
                                  std::vector<int>& chain_score_estimates, std::vector<std::vector<size_t>>& minimizer_kept_chain_count, 
                                  std::vector<double>& multiplicity_by_chain, vector<double>& multiplicity_by_tree,
                                  std::unordered_map<size_t, std::vector<size_t>>& good_fragments_in,
                                  LazyRNG& rng, Funnel& funnel) const;

    /**
     * Collect stats about the best chains for annotating the final alignment
     */
    void get_best_chain_stats( Alignment& aln, const ZipCodeForest& zip_code_forest, const std::vector<Seed>& seeds, 
                               const VectorView<MinimizerMapper::Minimizer>& minimizers,
                               const std::vector<std::vector<size_t>>& fragments,
                               const std::unordered_map<size_t, std::vector<size_t>>& good_fragments_in,
                               const std::vector<std::vector<size_t>>& chains,
                               const std::vector<size_t>& chain_source_tree,
                               const vector<algorithms::Anchor>& seed_anchors,
                               const std::vector<int>& chain_score_estimates,
                               bool& best_chain_correct, double& best_chain_coverage, size_t& best_chain_longest_jump, 
                               double& best_chain_average_jump, size_t& best_chain_anchors, size_t& best_chain_anchor_length, 
                               Funnel& funnel) const ;

    void do_alignment_on_chains(Alignment& aln, const std::vector<Seed>& seeds, 
                               const VectorView<MinimizerMapper::Minimizer>& minimizers, 
                               const vector<algorithms::Anchor>& seed_anchors,
                               const std::vector<std::vector<size_t>>& chains, 
                               const std::vector<size_t>& chain_source_tree,
                               const std::vector<double>& multiplicity_by_chain,
                               const std::vector<int>& chain_score_estimates,
                               const std::vector<std::vector<size_t>>& minimizer_kept_chain_count,
                               vector<Alignment>& alignments, vector<double>& multiplicity_by_alignment,
                               vector<size_t>& alignments_to_source,
                               SmallBitset& minimizer_explored, aligner_stats_t& stats, bool& funnel_depleted, LazyRNG& rng, Funnel& funnel) const;

    void pick_mappings_from_alignments(Alignment& aln, const std::vector<Alignment>& alignments, 
                                       const std::vector<double>& multiplicity_by_alignment, const std::vector<size_t>& alignments_to_source, 
                                       const std::vector<int>& chain_score_estimates,
                                       std::vector<Alignment>& mappings,
                                       std::vector<double>& scores, std::vector<double>& multiplicity_by_mapping,
                                       bool& funnel_depleted, LazyRNG& rng, Funnel& funnel) const;

    

    /**
     * Turn a chain into an Alignment.
     *
     * Operating on the given input alignment, align the tails and intervening
     * sequences along the given chain of perfect-match seeds, and return an
     * optimal Alignment.
     *
     * If given base processing stats for bases and for time, adds aligned bases and consumed time to them.
     */
    Alignment find_chain_alignment(const Alignment& aln, const VectorView<algorithms::Anchor>& to_chain, const std::vector<size_t>& chain, aligner_stats_t* stats = nullptr) const;
     
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
    
    /**
     * When dozeu doesn't have any seeds, it's scan heuristic can lead to
     * inaccurate anchoring with the end result that one end of the alignment
     * has a deletion that doesn't connect to an aligned base. This function
     * removes those deletions
     */
    void fix_dozeu_end_deletions(Alignment& rescued_alignment) const;

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
     * anchor must be set. Both anchors may be on the same node.
     *
     * Calls the callback with an extracted, strand-split, dagified graph, and
     * a function that translates from handle in the dagified graph to node ID
     * and orientation in the base graph.
     */
    static void with_dagified_local_graph(const pos_t& left_anchor, const pos_t& right_anchor, size_t max_path_length, const HandleGraph& graph, const std::function<void(DeletableHandleGraph&, const std::function<std::pair<nid_t, bool>(const handle_t&)>&)>& callback);
    
    /**
     * Determine the gap limit to use when aligning the given range of sequence
     * bases for the given Alignment.
     *
     * Accounts for the lognest gap that could be detected anywhere in the
     * range, not just at the very beginning or the very end, or at a single
     * point like GSSWAligner::longest_detectable_gap().
     */
    static size_t longest_detectable_gap_in_range(const Alignment& aln, const std::string::const_iterator& sequence_begin, const std::string::const_iterator& sequence_end, const GSSWAligner* aligner);

    /**
     * Clip out the part of the graph between the given positions and
     * global-align the sequence of the given Alignment to it. Populate the
     * Alignment's path and score.
     *
     * Finds an alignment against a graph path if it is <= max_path_length.
     *
     * If one of the anchor positions is empty, does pinned alignment against
     * the other position.
     *
     * For pinned alignment, restricts the alignment to have gaps no longer
     * than max_gap_length, and to use <= max_dp_cells cells. If too many DP
     * cells would be used, produces a softclip alignment.
     *
     * For connecting alignment, restricts the alignment to use <= max_dp_cells
     * cells. If too many DP cells would be used, produces an Alignment with
     * and empty path.
     *
     * Returns the number of nodes and bases in the graph aligned against.
     */
    static std::pair<size_t, size_t> align_sequence_between(const pos_t& left_anchor, const pos_t& right_anchor, size_t max_path_length, size_t max_gap_length, const HandleGraph* graph, const GSSWAligner* aligner, Alignment& alignment, const std::string* alignment_name = nullptr, size_t max_dp_cells = std::numeric_limits<size_t>::max(), const std::function<size_t(const Alignment&, const HandleGraph&)>& choose_band_padding = algorithms::pad_band_random_walk());

    /**
     * Version of align_sequence_between() that guarantees that you get the
     * same answer (modulo reverse-complementation) no matter whether the
     * sequence and anchors are reverse-complemented or not.
     */
    static std::pair<size_t, size_t> align_sequence_between_consistently(const pos_t& left_anchor, const pos_t& right_anchor, size_t max_path_length, size_t max_gap_length, const HandleGraph* graph, const GSSWAligner* aligner, Alignment& alignment, const std::string* alignment_name = nullptr, size_t max_dp_cells = std::numeric_limits<size_t>::max(), const std::function<size_t(const Alignment&, const HandleGraph&)>& choose_band_padding = algorithms::pad_band_random_walk());
    
    /**
     * Produce a WFAAlignment of the given sequence between the given points
     * that will be the same (modulo reverse-complementation) no matter whether
     * the sequence and anchors are reverse-complemented or not.
     */
    static WFAAlignment connect_consistently(const std::string& sequence, const pos_t& left_anchor, const pos_t& right_anchor, const WFAExtender& wfa_extender); 

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
     * minimizer_indices must be sorted by agglomeration end, and then by
     * agglomeration start, so they can be decomposed into nice rectangles.
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
     * item's number and the number of other items with the same or better score,
     * until min_count items are processed and either max_count
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
        const function<bool(size_t, size_t)>& process_item,
        const function<void(size_t)>& discard_item_by_count,
        const function<void(size_t)>& discard_item_by_score) const;
     
    /**
     * Same as the other process_until_threshold functions, except using a vector to supply scores.
     */
    template<typename Score = double>
    void process_until_threshold_b(const vector<Score>& scores,
        double threshold, size_t min_count, size_t max_count,
        LazyRNG& rng,
        const function<bool(size_t, size_t)>& process_item,
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
        const function<bool(size_t, size_t)>& process_item,
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

    /// Dump dotplot information for seeds.
    /// Displays one or more named collections of runs of seeds.
    static void dump_debug_dotplot(const std::string& name, const VectorView<Minimizer>& minimizers, const std::vector<Seed>& seeds, const std::vector<std::pair<std::string, std::vector<std::vector<size_t>>>>& seed_sets, const PathPositionHandleGraph* path_graph);

    /// Dump a graph
    static void dump_debug_graph(const HandleGraph& graph);
    
    /// Length at which we cut over to long-alignment logging.
    const static size_t LONG_LIMIT = 256;
    
    /// Count at which we cut over to summary logging.
    const static size_t MANY_LIMIT = 10;


    friend class TestMinimizerMapper;
};

template<typename Score>
void MinimizerMapper::process_until_threshold_a(size_t items, const function<Score(size_t)>& get_score,
    double threshold, size_t min_count, size_t max_count,
    LazyRNG& rng,
    const function<bool(size_t, size_t)>& process_item,
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
    const function<bool(size_t, size_t)>& process_item,
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
        const function<bool(size_t, size_t)>& process_item,
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

    // Find how many items have a better or equal score
    vector<size_t> better_or_equal_count(items, items);
    for (int i = items-2 ; i >= 0 ; --i) {
        //Starting from the second to last item, use the comparator to determine if it has the same
        // or lower score than the item after it
        if (comparator(indexes_in_order[i], indexes_in_order[i+1])){
            //If the score is less than the item after it
            better_or_equal_count[i] = i+1;
        } else {
            //Otherwise, they must be equal since they are ordered
            better_or_equal_count[i] = better_or_equal_count[i+1];
        }
    }

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
                unskipped += (size_t) process_item(item_num, better_or_equal_count[i]);
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
                unskipped += (size_t) process_item(item_num, better_or_equal_count[i]);
            } else {
                // We are out of room! Reject for count.
                discard_item_by_count(item_num);
            }
        }
    }
}

}



#endif
