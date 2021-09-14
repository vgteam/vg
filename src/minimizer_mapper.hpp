#ifndef VG_MINIMIZER_MAPPER_HPP_INCLUDED
#define VG_MINIMIZER_MAPPER_HPP_INCLUDED

/** 
 * \file minimizer_mapper.hpp
 * Defines a mapper that uses the minimizer index and GBWT-based extension.
 */

#include "algorithms/nearest_offsets_in_paths.hpp"
#include "aligner.hpp"
#include "vg/io/alignment_emitter.hpp"
#include "gapless_extender.hpp"
#include "mapper.hpp"
#include "min_distance.hpp"
#include "seed_clusterer.hpp"
#include "snarls.hpp"
#include "tree_subgraph.hpp"
#include "funnel.hpp"

#include <gbwtgraph/minimizer.h>
#include <structures/immutable_list.hpp>

#include <atomic>

namespace vg {

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
         MinimumDistanceIndex& distance_index, const PathPositionHandleGraph* path_graph = nullptr);

    /**
     * Map the given read, and send output to the given AlignmentEmitter. May be run from any thread.
     * TODO: Can't be const because the clusterer's cluster_seeds isn't const.
     */
    void map(Alignment& aln, AlignmentEmitter& alignment_emitter);
    
    /**
     * Map the given read. Return a vector of alignments that it maps to, winner first.
     */
    vector<Alignment> map(Alignment& aln);
    
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
    size_t hit_cap = 10;

    /// Ignore all minimizers with more than hard_hit_cap hits
    size_t hard_hit_cap = 500;

    /// Take minimizers between hit_cap and hard_hit_cap hits until this fraction
    /// of total score
    double minimizer_score_fraction = 0.9;

    ///Accept at least this many clusters
    size_t min_extensions = 2;

    /// How many clusters should we align?
    size_t max_extensions = 800;

    /// How many extended clusters should we align, max?
    size_t max_alignments = 8;
    
    /// How many extensions should we try as seeds within a mapping location?
    size_t max_local_extensions = numeric_limits<size_t>::max();

    //If a cluster's score is smaller than the best score of any cluster by more than
    //this much, then don't extend it
    double cluster_score_threshold = 50;
    
    //If the second best cluster's score is no more than this many points below
    //the cutoff set by cluster_score_threshold, snap that cutoff down to the
    //second best cluster's score, to avoid throwing away promising
    //secondaries.
    double pad_cluster_score_threshold = 20;

    //If the read coverage of a cluster is less than the best coverage of any cluster
    //by more than this much, don't extend it
    double cluster_coverage_threshold = 0.3;

    //If an extension set's score is smaller than the best 
    //extension's score by more than this much, don't align it
    double extension_set_score_threshold = 20;

    //If an extension's score is smaller than the best extension's score by
    //more than this much, don't align it
    int extension_score_threshold = 1;

    size_t max_multimaps = 1;
    size_t distance_limit = 200;
    bool do_dp = true;
    string sample_name;
    string read_group;
    
    /// Track which internal work items came from which others during each
    /// stage of the mapping algorithm.
    bool track_provenance = false;

    /// Guess which seed hits are correct by location in the linear reference
    /// and track if/when their descendants make it through stages of the
    /// algorithm. Only works if track_provenance is true.
    bool track_correctness = false;
    
    /// If set, log what the mapper is thinking in its mapping of each read.
    bool show_work = false;

    ////How many stdevs from fragment length distr mean do we cluster together?
    double paired_distance_stdevs = 2.0; 

    ///How close does an alignment have to be to the best alignment for us to rescue on it
    double paired_rescue_score_limit = 0.9;

    ///How many stdevs from the mean do we extract a subgraph from?
    double rescue_subgraph_stdevs = 4.0;

    /// Do not attempt rescue if there are more seeds in the rescue subgraph.
    size_t rescue_seed_limit = 100;

    /// For paired end mapping, how many times should we attempt rescue (per read)?
    size_t max_rescue_attempts = 15;
    
    /// How big of an alignment in POA cells should we ever try to do with Dozeu?
    /// TODO: Lift this when Dozeu's allocator is able to work with >4 MB of memory.
    /// Each cell is 16 bits in Dozeu, and we leave some room for the query and
    /// padding to full SSE registers. Note that a very chopped graph might
    /// still break this!
    size_t max_dozeu_cells = (size_t)(1.5 * 1024 * 1024);
    
    /// And have we complained about hitting it for rescue?
    atomic_flag warned_about_rescue_size = ATOMIC_FLAG_INIT;
    
    /// And have we complained about hitting it for tails?
    mutable atomic_flag warned_about_tail_size = ATOMIC_FLAG_INIT;

    ///What is the maximum fragment length that we accept as valid for paired-end reads?
    size_t max_fragment_length = 2000;

    /// Implemented rescue algorithms: no rescue, dozeu, GSSW, dozeu on local haplotypes.
    enum RescueAlgorithm { rescue_none, rescue_dozeu, rescue_gssw, rescue_haplotypes };

    /// The algorithm used for rescue.
    RescueAlgorithm rescue_algorithm = rescue_dozeu;

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
protected:

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

        // Sort the minimizers in descending order by score.
        inline bool operator< (const Minimizer& another) const {
            return (this->score > another.score);
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
    };
    
    /// The information we store for each seed.
    typedef SnarlSeedClusterer::Seed Seed;

    /// The information we store for each cluster.
    typedef SnarlSeedClusterer::Cluster Cluster;

    // These are our indexes
    const PathPositionHandleGraph* path_graph; // Can be nullptr; only needed for correctness tracking.
    const gbwtgraph::DefaultMinimizerIndex& minimizer_index;
    MinimumDistanceIndex& distance_index;
    /// This is our primary graph.
    const gbwtgraph::GBWTGraph& gbwt_graph;
    
    /// We have a gapless extender to extend seed hits in haplotype space.
    GaplessExtender extender;
    
    /// We have a clusterer
    SnarlSeedClusterer clusterer;

    FragmentLengthDistribution fragment_length_distr;
    atomic_flag warned_about_bad_distribution = ATOMIC_FLAG_INIT;

//-----------------------------------------------------------------------------

    // Stages of mapping.

    /**
     * Find the minimizers in the sequence using all minimizer indexes and
     * return them sorted in descending order by score.
     */
    std::vector<Minimizer> find_minimizers(const std::string& sequence, Funnel& funnel) const;

    /**
     * Find seeds for all minimizers passing the filters.
     */
    std::vector<Seed> find_seeds(const std::vector<Minimizer>& minimizers, const Alignment& aln, Funnel& funnel) const;

    /**
     * Determine cluster score, read coverage, and a vector of flags for the
     * minimizers present in the cluster. Score is the sum of the scores of
     * distinct minimizers in the cluster, while read coverage is the fraction
     * of the read covered by seeds in the cluster.
     */
    void score_cluster(Cluster& cluster, size_t i, const std::vector<Minimizer>& minimizers, const std::vector<Seed>& seeds, size_t seq_length, Funnel& funnel) const;

    /**
     * Score the set of extensions for each cluster using score_extension_group().
     * Return the scores in the same order as the extensions.
     */
    std::vector<int> score_extensions(const std::vector<std::vector<GaplessExtension>>& extensions, const Alignment& aln, Funnel& funnel) const;

    /**
     * Score the set of extensions for each cluster using score_extension_group().
     * Return the scores in the same order as the extensions.
     */
    std::vector<int> score_extensions(const std::vector<std::pair<std::vector<GaplessExtension>, size_t>>& extensions, const Alignment& aln, Funnel& funnel) const;

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
    void attempt_rescue(const Alignment& aligned_read, Alignment& rescued_alignment, const std::vector<Minimizer>& minimizers, bool rescue_forward);

    /**
     * Return the all non-redundant seeds in the subgraph, including those from
     * minimizers not used for mapping.
     */
    GaplessExtender::cluster_type seeds_in_subgraph(const std::vector<Minimizer>& minimizers, const std::unordered_set<id_t>& subgraph) const;

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
     * Get the distance between a pair of read alignments
     */
    int64_t distance_between(const Alignment& aln1, const Alignment& aln2);

    /**
     * Convert the GaplessExtension into an alignment. This assumes that the
     * extension is a full-length alignment and that the sequence field of the
     * alignment has been set.
     */
    void extension_to_alignment(const GaplessExtension& extension, Alignment& alignment) const;
    
    /**
     * Set pair partner references for paired mapping results.
     */
    void pair_all(pair<vector<Alignment>, vector<Alignment>>& mappings) const;
    


//-----------------------------------------------------------------------------

    /**
     * Compute MAPQ caps based on all minimizers that are explored, for some definition of explored.
     *
     * Needs access to the input alignment for sequence and quality
     * information.
     *
     * Returns only an "extended" cap at the moment.
     */
    double compute_mapq_caps(const Alignment& aln, const std::vector<Minimizer>& minimizers,
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
    static double window_breaking_quality(const vector<Minimizer>& minimizers, vector<size_t>& broken,
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
    static double faster_cap(const vector<Minimizer>& minimizers, vector<size_t>& minimizers_explored, const string& sequence, const string& quality_bytes);
    
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
    static void for_each_agglomeration_interval(const vector<Minimizer>& minimizers,
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
    static double get_log10_prob_of_disruption_in_interval(const vector<Minimizer>& minimizers,
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
    static double get_prob_of_disruption_in_column(const vector<Minimizer>& minimizers,
        const string& sequence, const string& quality_bytes,
        const vector<size_t>::iterator& disrupt_begin, const vector<size_t>::iterator& disrupt_end,
        size_t index);
    
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
     * Operating on the given input alignment, align the tails dangling off the
     * given extended perfect-match seeds and produce an optimal alignment into
     * the given output Alignment object, best, and the second best alignment
     * into second_best.
     *
     * Uses the given RNG to break ties.
     */
    void find_optimal_tail_alignments(const Alignment& aln, const vector<GaplessExtension>& extended_seeds, LazyRNG& rng, Alignment& best, Alignment& second_best) const; 
    
    /**
     * Find for each pair of extended seeds all the haplotype-consistent graph
     * paths against which the intervening read sequence needs to be aligned.
     *
     * Limits walks from each extended seed end to the longest detectable gap
     * plus the remaining to-be-alinged sequence, both computed using the read
     * length.
     *
     * extended_seeds must be sorted by read start position. Any extended seeds
     * that overlap in the read will be precluded from connecting.
     *
     * numeric_limits<size_t>::max() is used to store sufficiently long Paths
     * ending before sources (which cannot be reached from other extended
     * seeds) and starting after sinks (which cannot reach any other extended
     * seeds). Only sources and sinks have these "tail" paths.
     *
     * Tail paths are only calculated if the MinimizerMapper has linear_tails
     * set to true.
     */
    unordered_map<size_t, unordered_map<size_t, vector<Path>>> find_connecting_paths(const vector<GaplessExtension>& extended_seeds,
        size_t read_length) const;
        
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
     * Given a vector of items, a function to get the score of each, a
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
    template<typename Item, typename Score = double>
    void process_until_threshold_a(const vector<Item>& items, const function<Score(size_t)>& get_score,
        double threshold, size_t min_count, size_t max_count,
        LazyRNG& rng,
        const function<bool(size_t)>& process_item,
        const function<void(size_t)>& discard_item_by_count,
        const function<void(size_t)>& discard_item_by_score) const;
     
    /**
     * Same as the other process_until_threshold functions, except using a vector to supply scores.
     */
    template<typename Item, typename Score = double>
    void process_until_threshold_b(const vector<Item>& items, const vector<Score>& scores,
        double threshold, size_t min_count, size_t max_count,
        LazyRNG& rng,
        const function<bool(size_t)>& process_item,
        const function<void(size_t)>& discard_item_by_count,
        const function<void(size_t)>& discard_item_by_score) const;
     
    /**
     * Same as the other process_until_threshold functions, except user supplies 
     * comparator to sort the items (must still be sorted by score).
     */
    template<typename Item, typename Score = double>
    void process_until_threshold_c(const vector<Item>& items, const function<Score(size_t)>& get_score,
        const function<bool(size_t, size_t)>& comparator,
        double threshold, size_t min_count, size_t max_count,
        LazyRNG& get_seed,
        const function<bool(size_t)>& process_item,
        const function<void(size_t)>& discard_item_by_count,
        const function<void(size_t)>& discard_item_by_score) const;
        
    // Internal debugging functions
    
    /// Dump all the given minimizers, with optional subset restriction
    static void dump_debug_minimizers(const vector<Minimizer>& minimizers, const string& sequence, const vector<size_t>* to_include = nullptr);
    
    /// Dump all the extansions in an extension set
    static void dump_debug_extension_set(const HandleGraph& graph, const Alignment& aln, const vector<GaplessExtension>& extended_seeds);
    
    /// Print a sequence with base numbering
    static void dump_debug_sequence(ostream& out, const string& sequence);
    
    /// Get the thread identifier prefix for logging
    static string log_name();

    friend class TestMinimizerMapper;
};

template<typename Item, typename Score>
void MinimizerMapper::process_until_threshold_a(const vector<Item>& items, const function<Score(size_t)>& get_score,
    double threshold, size_t min_count, size_t max_count,
    LazyRNG& rng,
    const function<bool(size_t)>& process_item,
    const function<void(size_t)>& discard_item_by_count,
    const function<void(size_t)>& discard_item_by_score) const {

    process_until_threshold_c<Item, Score>(items, get_score, [&](size_t a, size_t b) -> bool {
        return (get_score(a) > get_score(b));
    },threshold, min_count, max_count, rng, process_item, discard_item_by_count, discard_item_by_score);
}

template<typename Item, typename Score>
void MinimizerMapper::process_until_threshold_b(const vector<Item>& items, const vector<Score>& scores,
    double threshold, size_t min_count, size_t max_count,
    LazyRNG& rng,
    const function<bool(size_t)>& process_item,
    const function<void(size_t)>& discard_item_by_count,
    const function<void(size_t)>& discard_item_by_score) const {
    
    assert(scores.size() == items.size());
    
    process_until_threshold_c<Item, Score>(items, [&](size_t i) -> Score {
        return scores[i];
    }, [&](size_t a, size_t b) -> bool {
        return (scores[a] > scores[b]);
    },threshold, min_count, max_count, rng, process_item, discard_item_by_count, discard_item_by_score);
}

template<typename Item, typename Score>
void MinimizerMapper::process_until_threshold_c(const vector<Item>& items, const function<Score(size_t)>& get_score,
        const function<bool(size_t, size_t)>& comparator,
        double threshold, size_t min_count, size_t max_count,
        LazyRNG& rng,
        const function<bool(size_t)>& process_item,
        const function<void(size_t)>& discard_item_by_count,
        const function<void(size_t)>& discard_item_by_score) const {

    // Sort item indexes by item score
    vector<size_t> indexes_in_order;
    indexes_in_order.reserve(items.size());
    for (size_t i = 0; i < items.size(); i++) {
        indexes_in_order.push_back(i);
    }
    
    // Put the highest scores first, but shuffle top ties so reads spray evenly
    // across equally good mappings
    sort_shuffling_ties(indexes_in_order.begin(), indexes_in_order.end(), comparator, rng);

    // Retain items only if their score is at least as good as this
    double cutoff = items.size() == 0 ? 0 : get_score(indexes_in_order[0]) - threshold;
    
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
