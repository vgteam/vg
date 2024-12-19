#ifndef VG_ALGORITHMS_CHAIN_ITEMS_HPP_INCLUDED
#define VG_ALGORITHMS_CHAIN_ITEMS_HPP_INCLUDED

/**
 * \file
 * Algorithms for chaining subalignments into larger alignments.
 *
 * To use these algorithms, decide on the type (Anchor) you want to chain up.
 *
 * Then, make a dynamic programming table: vector<TracedScore>.
 *
 * Then, call chain_items_dp() to fill in the dynamic programming table and get
 * the score of the best chain.
 *
 * You can use chain_items_traceback() to get a traceback of the chain's items
 * in order.
 *
 * Helper entry points are find_best_chain() and score_best_chain() which set
 * up the DP for you and do the traceback if appropriate.
 */
 
#include "extract_containing_graph.hpp"

#include "../gbwt_extender.hpp"
#include "../snarl_seed_clusterer.hpp"
#include "../zip_code_tree.hpp"
#include "../handle.hpp"
#include "../explainer.hpp"
#include "../utility.hpp"

#include <bdsg/hash_graph.hpp>

namespace vg {
namespace algorithms {

using namespace std;

// Make sure all of vg's print operators are available.
using vg::operator<<;

//#define debug_chaining

/**
 * Represents a piece fo a graph node matching to a piece of a read. Can be
 * chained together.
 */
class Anchor {
public:
    /// Get the start position in the read of this anchor's match.
    inline size_t read_start() const {
        return start;
    }

    /// Get the start position in the graph of this anchor's match
    inline const pos_t& graph_start() const {
        return start_pos;
    }

    /// Get the start position in the read of the part of the read that you
    /// can't have another anchor in if you take this one.
    ///
    /// We trimmed the anchors down from the minimizers to avoid having to deal
    /// with the tail ends of the minimizers going multiple places in the
    /// graph. But we don't want to let you take anchors from minimizers that
    /// overlapped.
    inline size_t read_exclusion_start() const {
        return read_start() - margin_before;
    }

    /// Get the length of this anchor's match
    inline size_t length() const {
        return size;
    }
    /// Get the alignment score of the anchor
    inline int score() const {
        return points;
    }
    
    /// Get the end position in the read of this anchor's match
    inline size_t read_end() const {
        return read_start() + length();
    }

    /// Get the end position in the graph of this anchor's match
    inline pos_t graph_end() const {
        return end_pos;
    }
    
    /// Get the end position in the read of the part of the read that you
    /// can't have another anchor in if you take this one.
    inline size_t read_exclusion_end() const {
        return read_end() + margin_after;
    }
    
    /// Get the number of the seed at the start of the anchor, or
    /// std::numeric_limits<size_t>::max() if not set.
    inline size_t seed_start() const {
        return start_seed;
    }
    
    /// Get the number of the seed at the end of the chain, or
    /// std::numeric_limits<size_t>::max() if not set.
    inline size_t seed_end() const {
        return end_seed;
    }

    /// Get the distance-finding hint information (i.e. "zip code") for
    /// accelerating distance queries to the start of this anchor, or null if
    /// none is set.
    inline ZipCode* start_hint() const {
        return start_zip;
    }

    /// Get the graph distance from wherever the start hint is positioned back
    /// to the actual start of the anchor.
    inline size_t start_hint_offset() const {
        return start_offset;
    }
    
    /// Get the distance-finding hint information (i.e. "zip code") for
    /// accelerating distance queries from the end of this anchor, or null if
    /// none is set.
    inline ZipCode* end_hint() const {
        return end_zip;
    }

    /// Get the graph distance from wherever the end hint is positioned forward
    /// to the actual end of the anchor.
    inline size_t end_hint_offset() const {
        return end_offset;
    }

    /// Get the length of the exclusion zone for a primary anchor, or the
    /// average such length of the anchors this anchor is made from for a
    /// composite anchor. This is used in gap scoring during chaining, to make
    /// sure gap scores don't get enormous for long composite anchors.
    inline size_t base_seed_length() const {
        return seed_length; 
    }

    /// Can this anchor be skipped when constructing an alignment?
    /// This will be true if the anchor is in a repetitive region of the read
    inline bool is_skippable() const {
        return skippable;
    }
    
    // Construction
    
    /// Compose a read start position, graph start position, and match length into an Anchor.
    /// Can also bring along a distance hint and a seed number.
    inline Anchor(size_t read_start, const pos_t& graph_start, size_t length, size_t margin_before, size_t margin_after, int score, size_t seed_number = std::numeric_limits<size_t>::max(), ZipCode* hint = nullptr, size_t hint_start = 0, bool skippable = false) : start(read_start), size(length), margin_before(margin_before), margin_after(margin_after), start_pos(graph_start), end_pos(advance(graph_start, length)), points(score), start_seed(seed_number), end_seed(seed_number), start_zip(hint), end_zip(hint), start_offset(hint_start), end_offset(length - hint_start), seed_length(margin_before + length + margin_after), skippable(skippable) {
        // Nothing to do!
    }
    
    /// Compose two Anchors into an Anchor that represents coming in through
    /// the first one and going out through the second, like a tunnel. Useful
    /// for representing chains as chainable items.
    inline Anchor(const Anchor& first, const Anchor& last, size_t extra_margin_before, size_t extra_margin_after, int score) : start(first.read_start()), size(last.read_end() - first.read_start()), margin_before(first.margin_before + extra_margin_before), margin_after(last.margin_after + extra_margin_after), start_pos(first.graph_start()), end_pos(last.graph_end()), points(score), start_seed(first.seed_start()), end_seed(last.seed_end()), start_zip(first.start_hint()), end_zip(last.end_hint()), start_offset(first.start_offset), end_offset(last.end_offset), seed_length((first.base_seed_length() + last.base_seed_length()) / 2), skippable(first.is_skippable() || last.is_skippable()) {
        // Nothing to do!
    }

    // Act like data
    Anchor() = default;
    Anchor(const Anchor& other) = default;
    Anchor& operator=(const Anchor& other) = default;
    Anchor(Anchor&& other) = default;
    Anchor& operator=(Anchor&& other) = default;
    
protected:
    size_t start;
    size_t size;
    size_t margin_before;
    size_t margin_after;
    pos_t start_pos;
    pos_t end_pos;
    int points;
    size_t start_seed;
    size_t end_seed;
    ZipCode* start_zip;
    ZipCode* end_zip;
    size_t start_offset;
    size_t end_offset;
    size_t seed_length;
    bool skippable;
};

/// Explain an Anchor to the given stream
ostream& operator<<(ostream& out, const Anchor& anchor);

// For doing scores with backtracing, we use this type, which is a
// score and a number for the place we came from to get it.
class TracedScore {
public:
    /// What is the sentinel for an empty provenance?
    /// Use a function instead of a constant because that's easier when we're just a header.
    inline static size_t nowhere() {
        return numeric_limits<size_t>::max();
    }
    
    /// What's the default value for an empty table cell?
    /// Use a function instead of a constant because that's easier when we're just a header.
    inline static TracedScore unset() {
        return {0, nowhere()};
    }
    
    /// Max in a score from a DP table. If it wins, record provenance.
    void max_in(const vector<TracedScore>& options, size_t option_number);
    
    /// Get a score from a table of scores and record provenance in it.
    static TracedScore score_from(const vector<TracedScore>& options, size_t option_number);
    
    /// Add (or remove) points along a route to somewhere. Return a modified copy.
    TracedScore add_points(int adjustment) const;
    
    /// Compare for equality
    inline bool operator==(const TracedScore& other) const {
        return score == other.score && source == other.source;
    }
    
    /// Compare for inequality
    inline bool operator!=(const TracedScore& other) const {
        return !(*this == other);
    }
    
    /// Compare for less-than
    inline bool operator<(const TracedScore& other) const {
        return score < other.score || (score == other.score && source < other.source);
    }
    
    /// Compare for greater-than
    inline bool operator>(const TracedScore& other) const {
        return score > other.score || (score == other.score && source > other.source);
    }
    
    /// Subtraction to yield a difference in points
    inline int operator-(const TracedScore& other) const {
        return score - other.score;
    }
    
    // Number of points
    int score;
    // Index of source score among possibilities/traceback pointer
    size_t source;
};

}

}

namespace std {
    /// Allow maxing TracedScore
    inline vg::algorithms::TracedScore max(const vg::algorithms::TracedScore& a, const vg::algorithms::TracedScore& b) {
        return a > b ? a : b;
    }
}

namespace vg {

namespace algorithms {

using namespace std;

// Make sure all of vg's print operators are available.
using vg::operator<<;

/// Print operator
ostream& operator<<(ostream& out, const TracedScore& value);

/**
 * Sort indexes in the given list by by read start position (and end position)
 * of the anchors they refer to.
 */
void sort_anchor_indexes(const std::vector<Anchor>& items, std::vector<size_t>& indexes);

/**
 * Iteratee function type which can be called with each transition between
 * anchors.
 * 
 * Takes two anchor numbers (source and destination), and their read and graph
 * distances, in that order.
 */
using transition_iteratee = std::function<void(size_t from_anchor, size_t to_anchor, size_t read_distance, size_t graph_distance)>;

/**
 * Iterator function type which lets you iterate over transitions between
 * items, by calling a callback.
 *
 * Implementation will go throuch all the anchors and call the given callback
 * with pairs of anchor numbers, and their read and graph distances.
 * 
 * Transitions are always between anchors earlier and later in the read.
 * 
 * Transitions are from the first anchor, to the second.
 * 
 * Transitions are visited in order: all transititions to an anchor are visited
 * before any transitions from it.
 * 
 * to_chain must be sorted by read start.
 */
using transition_iterator = std::function<void(const VectorView<Anchor>& to_chain, const SnarlDistanceIndex& distance_index, const HandleGraph& graph, size_t max_indel_bases, const transition_iteratee& callback)>;

/**
 * Return a transition iterator that iterates along the read and uses the given lookback control parameters to filter transitions.
 * Closes over the arguments by value.
 */
transition_iterator lookback_transition_iterator(size_t max_lookback_bases,
                                                 size_t min_lookback_items,
                                                 size_t lookback_item_hard_cap);

/**
 * Return a transition iterator that uses zip code tree iteration to select traversals.
 */
transition_iterator zip_tree_transition_iterator(const std::vector<SnarlDistanceIndexClusterer::Seed>& seeds, const ZipCodeTree& zip_code_tree, size_t max_lookback_bases);

/**
 * Fill in the given DP table for the explored chain scores ending with each
 * item. Returns the best observed score overall from that table, with
 * provenance to its location in the table, if tracked in the type. Assumes
 * some items exist.
 *
 * We keep all the options to allow us to do multiple tracebacks and find
 * multiple good (ideally disjoint) chains.
 *
 * Input items must be sorted by start position in the read.
 *
 * Takes the given per-item bonus for each item collected, and scales item scores by the given scale.
 *
 * Uses a finite lookback in items and in read bases when checking where we can
 * come from to reach an item. Also, once a given number of good-looking
 * predecessor items have been found, stop looking back.
 *
 * Limits transitions to those involving indels of the given size or less, to
 * avoid very bad transitions.
 */
TracedScore chain_items_dp(vector<TracedScore>& chain_scores,
                           const VectorView<Anchor>& to_chain,
                           const SnarlDistanceIndex& distance_index,
                           const HandleGraph& graph,
                           int gap_open,
                           int gap_extension,
                           const transition_iterator& for_each_transition = lookback_transition_iterator(150, 0, 100),
                           int item_bonus = 0,
                           double item_scale = 1.0,
                           double gap_scale = 1.0,
                           double points_per_possible_match = 0,
                           size_t max_indel_bases = 100,
                           bool show_work = false);

/**
 * Trace back through in the given DP table from the best chain score.
 *
 * Returns tracebacks that visit disjoint sets of items, in score order, along
 * with their penalties from the optimal score. The best_past_ending_score_ever
 * is *not* always the source of the first traceback, if there is a tie.
 *
 *  Tracebacks are constrained to be nonoverlapping by stopping each traceback
 *  when the optimum place to come from has already been used. The second-best
 *  place to come from is *not* considered. It might be possible that two
 *  returned tracebacks could be pasted together to get a higher score, but it
 *  won't be possible to recombine two tracebacks to get a higher score; no
 *  edges followed between items will ever need to be cut.
 */
vector<pair<vector<size_t>, int>> chain_items_traceback(const vector<TracedScore>& chain_scores,
                                                        const VectorView<Anchor>& to_chain,
                                                        const TracedScore& best_past_ending_score_ever,
                                                        int item_bonus = 0,
                                                        double item_scale = 1.0,
                                                        size_t max_tracebacks = 1);


/**
 * Chain up the given group of items. Determines the best scores and
 * tracebacks that can be obtained by chaining items together.
 *
 * Input items must be sorted by start position in the read.
 *
 * Returns the scores and the list of indexes of items visited to achieve
 * that score, in order, with multiple tracebacks in descending score order.
 */
vector<pair<int, vector<size_t>>> find_best_chains(const VectorView<Anchor>& to_chain,
                                                   const SnarlDistanceIndex& distance_index,
                                                   const HandleGraph& graph,
                                                   int gap_open,
                                                   int gap_extension,
                                                   size_t max_chains = 1,
                                                   const transition_iterator& for_each_transition = lookback_transition_iterator(150, 0, 100), 
                                                   int item_bonus = 0,
                                                   double item_scale = 1.0,
                                                   double gap_scale = 1.0,
                                                   double points_per_possible_match = 0,
                                                   size_t max_indel_bases = 100,
                                                   bool show_work = false);

/**
 * Chain up the given group of items. Determines the best score and
 * traceback that can be obtained by chaining items together.
 *
 * Input items must be sorted by start position in the read.
 *
 * Returns the score and the list of indexes of items visited to achieve
 * that score, in order.
 */
pair<int, vector<size_t>> find_best_chain(const VectorView<Anchor>& to_chain,
                                          const SnarlDistanceIndex& distance_index,
                                          const HandleGraph& graph,
                                          int gap_open,
                                          int gap_extension,
                                          const transition_iterator& for_each_transition = lookback_transition_iterator(150, 0, 100),
                                          int item_bonus = 0,
                                          double item_scale = 1.0,
                                          double gap_scale = 1.0,
                                          double points_per_possible_match = 0,
                                          size_t max_indel_bases = 100);
                                          
/**
 * Score the given group of items. Determines the best score that can be
 * obtained by chaining items together.
 *
 * Input items must be sorted by start position in the read.
 */
int score_best_chain(const VectorView<Anchor>& to_chain, const SnarlDistanceIndex& distance_index, const HandleGraph& graph, int gap_open, int gap_extension);


/// Score a chaining gap using the Minimap2 method. See
/// <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6137996/> near equation 2.
/// This produces a penalty (positive number).
int score_chain_gap(size_t distance_difference, size_t average_anchor_length);

/// Get distance in the graph, or std::numeric_limits<size_t>::max() if unreachable or beyond the limit.
size_t get_graph_distance(const Anchor& from, const Anchor& to, const SnarlDistanceIndex& distance_index, const HandleGraph& graph, size_t distance_limit = std::numeric_limits<size_t>::max());

/// Get distance in the read, or std::numeric_limits<size_t>::max() if unreachable.
size_t get_read_distance(const Anchor& from, const Anchor& to);
                     
}
}

#endif
