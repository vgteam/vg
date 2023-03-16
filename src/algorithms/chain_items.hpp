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
    // Set up with accessors in case we want to stop copying stuff so much later.

    // Base API:
    
    /// Get the start position in the read of this anchor's match.
    inline size_t read_start() const {
        return start;
    }
    /// Get the start position in the graph of this anchor's match
    inline const pos_t& graph_start() const {
        return pos;
    }
    /// Get the length of this anchor's match
    inline size_t length() const {
        return size;
    }
    /// Get the alignment score of the anchor
    inline int score() const {
        return points;
    }
    
    // Other API implemented on top of this
    
    /// Get the end position in the read of this anchor's match
    inline size_t read_end() const {
        return read_start() + length();
    }
    
    /// Get the end position in the graph of this anchor's match
    inline pos_t graph_end() const {
        pos_t p = graph_start();
        get_offset(p) += length();
        return p;
    }
    
    /// Get the distance-finding hint information (i.e. "zip code") for
    /// accelerating distance queries, or null if none is set.
    inline ZipCodeDecoder* hint() const {
        return decoder;
    };
    
    // Construction
    
    /// Compose a read start position, graph start position, and match length into an Anchor.
    /// Can also bring along a distance hint
    inline Anchor(size_t read_start, const pos_t& graph_start, size_t length, int score, ZipCodeDecoder* hint = nullptr) : start(read_start), size(length), pos(graph_start), points(score), decoder(hint) {
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
    pos_t pos;
    int points;
    ZipCodeDecoder* decoder;
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
 * Get rid of items that are shadowed or contained by (or are identical to) others.
 *
 * Erases items that didn't survive from indexes, and sorts them by read start
 * position.
 */
void sort_and_shadow(const std::vector<Anchor>& items, std::vector<size_t>& indexes);

/**
 * Get rid of items that are shadowed or contained by (or are identical to) others.
 *
 * Erases items that didn't survive from items, and sorts them by read start
 * position.
 */
void sort_and_shadow(std::vector<Anchor>& items);

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
 * Takes the given per-item bonus for each item collected.
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
                           size_t max_lookback_bases = 150,
                           size_t min_lookback_items = 0,
                           size_t lookback_item_hard_cap = 100,
                           size_t initial_lookback_threshold = 10,
                           double lookback_scale_factor = 2.0,
                           double min_good_transition_score_per_base = -0.1,
                           int item_bonus = 0,
                           size_t max_indel_bases = 100);

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
                                                   size_t max_lookback_bases = 150,
                                                   size_t min_lookback_items = 0,
                                                   size_t lookback_item_hard_cap = 100,
                                                   size_t initial_lookback_threshold = 10,
                                                   double lookback_scale_factor = 2.0,
                                                   double min_good_transition_score_per_base = -0.1,
                                                   int item_bonus = 0,
                                                   size_t max_indel_bases = 100);

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
                                          size_t max_lookback_bases = 150,
                                          size_t min_lookback_items = 0,
                                          size_t lookback_item_hard_cap = 100,
                                          size_t initial_lookback_threshold = 10,
                                          double lookback_scale_factor = 2.0,
                                          double min_good_transition_score_per_base = -0.1,
                                          int item_bonus = 0,
                                          size_t max_indel_bases = 100);
                                          
/**
 * Score the given group of items. Determines the best score that can be
 * obtained by chaining items together.
 *
 * Input items must be sorted by start position in the read.
 */
int score_best_chain(const VectorView<Anchor>& to_chain, const SnarlDistanceIndex& distance_index, const HandleGraph& graph, int gap_open, int gap_extension);

/// Get distance in the graph, or std::numeric_limits<size_t>::max() if unreachable or beyond the limit.
size_t get_graph_distance(const Anchor& from, const Anchor& to, const SnarlDistanceIndex& distance_index, const HandleGraph& graph, size_t distance_limit = std::numeric_limits<size_t>::max());

/// Get distance in the read, or std::numeric_limits<size_t>::max() if unreachable.
size_t get_read_distance(const Anchor& from, const Anchor& to);
                     
}
}

#endif
