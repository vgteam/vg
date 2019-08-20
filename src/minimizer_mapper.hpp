#ifndef VG_MINIMIZER_MAPPER_HPP_INCLUDED
#define VG_MINIMIZER_MAPPER_HPP_INCLUDED

/** 
 * \file minimizer_mapper.hpp
 * Defines a mapper that uses the minimizer index and GBWT-based extension.
 */

#include "aligner.hpp"
#include "minimizer.hpp"
#include "alignment_emitter.hpp"
#include "gapless_extender.hpp"
#include "snarls.hpp"
#include "min_distance.hpp"
#include "seed_clusterer.hpp"
#include "tree_subgraph.hpp"
#include "algorithms/nearest_offsets_in_paths.hpp"

#include <structures/immutable_list.hpp>

namespace vg {

using namespace std;

class MinimizerMapper : public AlignerClient {
public:

    /**
     * Construct a new MinimizerMapper using the given indexes. The PathPositionhandleGraph can be nullptr,
     * as we only use it for correctness tracking.
     */

    MinimizerMapper(const GBWTGraph& graph, const MinimizerIndex& minimizer_index,
         MinimumDistanceIndex& distance_index, const PathPositionHandleGraph* path_graph = nullptr);

    /**
     * Map the given read, and send output to the given AlignmentEmitter. May be run from any thread.
     * TODO: Can't be const because the clusterer's cluster_seeds isn't const.
     */
    void map(Alignment& aln, AlignmentEmitter& alignment_emitter);

    // Mapping settings.
    // TODO: document each

    /// Use all minimizers with at most hit_cap hits
    size_t hit_cap = 10;

    /// Ignore all minimizers with more than hard_hit_cap hits
    size_t hard_hit_cap = 300;

    /// Take minimizers between hit_cap and hard_hit_cap hits until this fraction
    /// of total score
    double minimizer_score_fraction = 0.6;

    /// How many clusters should we align?
    size_t max_extensions = 48;

    /// How many extended clusters should we align, max?
    size_t max_alignments = 8;
    
    /// How many extensions should we try as seeds within a mapping location?
    size_t max_local_extensions = numeric_limits<size_t>::max();

    //If a cluster's score is smaller than the best score of any cluster by more than
    //this much, then don't extend it
    double cluster_score_threshold = 0;

    //If the read coverage of a cluster is less than the best coverage of any cluster
    //by more than this much, don't extend it
    double cluster_coverage_threshold = 0;

    //If an extension set's score is smaller than the best 
    //extension's score by more than this much, don't align it
    double extension_set_score_threshold = 0;

    //If an extension's score is smaller than the best extension's score by
    //more than this much, don't align it
    int extension_score_threshold = 1;

    size_t max_multimaps = 1;
    size_t distance_limit = 1000;
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
    
protected:
    // These are our indexes
    const PathPositionHandleGraph* path_graph; // Can be nullptr; only needed for correctness tracking.
    const MinimizerIndex& minimizer_index;
    MinimumDistanceIndex& distance_index;

    /// This is our primary graph.
    const GBWTGraph& gbwt_graph;
    
    /// We have a gapless extender to extend seed hits in haplotype space.
    GaplessExtender extender;
    
    /// We have a clusterer
    SnarlSeedClusterer clusterer;
    
    /**
     * Estimate the score it may be possible to achieve using the given group of GaplessExtensions.
     * Supports single full-length extensions and groups that need chaining.
     * May reorder the input extended_seeds vector if it is not sorted in read space.
     * Is not always an overestimate of the actual score.
     */
    int estimate_extension_group_score(const Alignment& aln, vector<GaplessExtension>& extended_seeds) const;
    
    /**
     * Operating on the given input alignment, align the tails dangling off the
     * given extended perfect-match seeds and produce an optimal alignment into
     * the given output Alignment object, best, and the second best alignment
     * into second_best.
     */
    void find_optimal_tail_alignments(const Alignment& aln, const vector<GaplessExtension>& extended_seeds, Alignment& best, Alignment& second_best) const; 
    
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
     */
    vector<TreeSubgraph> get_tail_forest(const GaplessExtension& extended_seed,
        size_t read_length, bool left_tails) const;
        
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
     * Returns alingments in gbwt_graph space.
     */
    pair<Path, size_t> get_best_alignment_against_any_tree(const vector<TreeSubgraph>& trees, const string& sequence,
        const Position& default_position, bool pin_left) const;
        
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
     * Given a vector of items, a function to get the score of each, a
     * score-difference-from-the-best cutoff, and a min and max processed item
     * count, process items in descending score order by calling process_item
     * with the item's number, until min_count items are processed and either
     * max_count items are processed or the score difference threshold is hit
     * (or we run out of items).
     *
     * If process_item returns false, the item is skipped and does not count
     * against min_count or max_count.
     *
     * Call discard_item with the item's number for all remaining items.
     */
    template<typename Item, typename Score = double>
    void process_until_threshold(const vector<Item>& items, const function<Score(size_t)>& get_score,
        double threshold, size_t min_count, size_t max_count,
        const function<bool(size_t)>& process_item, const function<void(size_t)>& discard_item) const;
     
    /**
     * Same as the other process_until_threshold overload, except using a vector to supply scores.
     */
    template<typename Item, typename Score = double>
    void process_until_threshold(const vector<Item>& items, const vector<Score>& scores,
        double threshold, size_t min_count, size_t max_count,
        const function<bool(size_t)>& process_item, const function<void(size_t)>& discard_item) const;
     
};

template<typename Item, typename Score>
void MinimizerMapper::process_until_threshold(const vector<Item>& items, const function<Score(size_t)>& get_score,
    double threshold, size_t min_count, size_t max_count,
    const function<bool(size_t)>& process_item, const function<void(size_t)>& discard_item) const {

    // Sort item indexes by item score
    vector<size_t> indexes_in_order;
    indexes_in_order.reserve(items.size());
    for (size_t i = 0; i < items.size(); i++) {
        indexes_in_order.push_back(i);
    }
    
    // Put the highest scores first
    std::sort(indexes_in_order.begin(), indexes_in_order.end(), [&](const size_t& a, const size_t& b) -> bool {
        // Return true if a must come before b, and false otherwise
        return get_score(a) > get_score(b);
    });

    // Retain items only if their score is at least as good as this
    double cutoff = items.size() == 0 ? 0 : get_score(indexes_in_order[0]) - threshold;
    
    // Count up non-skipped items for min_count and max_count
    size_t unskipped = 0;
    
    // Go through the items in descending score order.
    for (size_t i = 0; i < indexes_in_order.size() && unskipped < max_count ; i++) {
        // Find the item we are talking about
        size_t& item_num = indexes_in_order[i];
        
        if (unskipped < min_count || (threshold == 0 || get_score(item_num) > cutoff)) {
            // Do we definitely want this one, or is its score sufficiently good?
            
            // If so, go do it.
            // If it is not skipped, add it to the total number of unskipped items
            unskipped += (size_t) process_item(item_num);
        } else {
            // If this score is insignificant, nothing past here is significant.
            
            // Say we aren't doing this one.
            // We'll also say we aren't doing the rest.
            discard_item(item_num);
        }
    }
}

template<typename Item, typename Score>
void MinimizerMapper::process_until_threshold(const vector<Item>& items, const vector<Score>& scores,
    double threshold, size_t min_count, size_t max_count,
    const function<bool(size_t)>& process_item, const function<void(size_t)>& discard_item) const {
    
    assert(scores.size() == items.size());
    
    process_until_threshold<Item, Score>(items, [&](size_t i) -> Score {
        return scores.at(i);
    }, threshold, min_count, max_count, process_item, discard_item);
    
}

}



#endif
