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

/// Represents a vount of alignment operations, some of which may be unspecified.
/// Supports addition and score finding under a scoring regime.
struct Operations {
    int matches;
    int mismatches;
    int opens;
    int extends;
    int unknowns;

    /// Allow default construction as zero
    inline Operations(): matches(0), mismatches(0), opens(0), extends(0) {
        // Nothing to do
    }

    /// Allow construction from a bunch of counts
    inline Operations(int matches, int mismatches, int opens, int extends, int unknowns): matches(matches), mismatches(mismatches), opens(opens), extends(extends), unknowns(unknowns) {
        // Nothing to do
    }

    /// Allow copy and move
    inline Operations(const Operations& other) = default;
    inline Operations(Operations&& other) = default;
    inline Operations& operator=(const Operations& other) = default;
    inline Operations& operator=(Operations&& other) = default;
    
    /// Add one collection of operations into another
    inline Operations& operator+=(const Operations& other) {
        matches += other.matches;
        mismatches += other.mismatches;
        opens += other.opens;
        extends += other.extends;
        unknowns += other.unknowns;
        return *this;
    }
    
    /// Add one collection of operations to another
    inline Operations operator+(const Operations& other) const {
        Operations added(*this);
        added += other;
        return added;
    }

    /// Allow negating a collection of operations
    inline Operations operator-() const {
        Operations copy(*this);
        copy.matches = -copy.matches;
        copy.mismatches = -copy.mismatches;
        copy.opens = -copy.opens;
        copy.extends = -copy.extends;
        copy.unknowns = -copy.unknowns;
        return copy;
    }

    /// Allow subtracting a collection of operations from this one
    inline Operations& operator-=(const Operations& other) {
        return (*this) += -other;
    }

    /// Allow subtracting two collections of operations to get a difference
    inline Operations operator-(const Operations& other) const {
        Operations copy = -other;
        copy += *this;
        return copy;
    }

    /// Make a match operation
    inline static Operations match(int count) {
        return {count, 0, 0, 0, 0};
    }

    /// Make a mismatch operation
    inline static Operations mismatch(int count) {
        return {0, count, 0, 0, 0};
    }

    /// Make a gap open operation
    inline static Operations open(int count) {
        return {0, 0, count, 0, 0};
    }

    /// Make a gap extend operation
    inline static Operations extend(int count) {
        return {0, 0, 0, count, 0};
    }

    /// Make an unknown/not yet determined operation
    inline static Operations unknown(int count) {
        return {0, 0, 0, 0, count};
    }

    /// Rescore according to the given operation scores, with penalties
    /// negative. Returns the computed score and leaves the object unmodified.
    inline int score_under(int match, int mismatch, int open, int extend) const {
        return match * matches + mismatch * mismatches + open * opens + extend * extends;
    }

    /// Rescore according to the given operation scores, with penalties
    /// negative, and assuming all unknown read bases are matches. Returns the
    /// computed score and leaves the object unmodified.
    inline int max_score_under(int match, int mismatch, int open, int extend) const {
        return match * (matches + unknowns) + mismatch * mismatches + open * opens + extend * extends;
    }
};

/// Represents a set of alignment operations together with a precomputed score.
struct ScoredOperations: public Operations {
    int score;

    /// Allow default construction as zero
    inline ScoredOperations(): Operations(), score(0) {
        // Nothing to do
    }

    /// Allow construction from a score and a bunch of counts
    inline ScoredOperations(int score, int matches, int mismatches, int opens, int extends, int unknowns): Operations(matches, mismatches, opens, extends, unknowns), score(score) {
        // Nothing to do
    }

    /// Allow construction from a score and Operations
    inline ScoredOperations(int score, const Operations& operations): Operations(operations), score(score) {
        // Nothing to do
    }

    /// Allow copy and move
    inline ScoredOperations(const ScoredOperations& other) = default;
    inline ScoredOperations(ScoredOperations&& other) = default;
    inline ScoredOperations& operator=(const ScoredOperations& other) = default;
    inline ScoredOperations& operator=(ScoredOperations&& other) = default;

    /// Add one collection of scored operations into another
    inline ScoredOperations& operator+=(const ScoredOperations& other) {
        Operations::operator+=(other);
        score += other.score;
        return *this;
    }

    /// Allow adding points
    inline ScoredOperations& operator+=(int points) {
        score += points;
        return *this;
    }
    
    /// Add one collection of scored operations to another
    inline ScoredOperations operator+(const ScoredOperations& other) const {
        ScoredOperations added(*this);
        added += other;
        return added;
    }

    /// Allow adding points to us
    inline ScoredOperations operator+(int points) const {
        ScoredOperations added(*this);
        added += points;
        return added;
    }

    /// Allow negating a collection of operations
    inline ScoredOperations operator-() const {
        ScoredOperations copy(-score, -*(const Operations*)this);
        return copy;
    }

    /// Allow subtracting a collection of operations from this one
    inline ScoredOperations& operator-=(const ScoredOperations& other) {
        Operations::operator-=(other);
        score -= other.score;
        return *this;
    }

    /// Allow subtracting two collections of operations to get a difference
    inline ScoredOperations operator-(const ScoredOperations& other) const {
        ScoredOperations copy = -other;
        copy += *this;
        return copy;
    }

    /// Allow multiplying a scale into the points
    inline ScoredOperations& operator*=(double scale) {
        score *= scale;
        return *this;
    }

    /// Allow multiplying the points by a scale
    inline ScoredOperations operator*(double scale) const {
        ScoredOperations multiplied(*this);
        multiplied *= scale;
        return multiplied;
    }

    /// Compare equality based only on score
    inline bool operator==(const ScoredOperations& other) const {
        return score == other.score;
    }

    /// Compare inequality based only on score
    inline bool operator!=(const ScoredOperations& other) const {
        return score != other.score;
    }

    /// Compare less than based only on score
    inline bool operator<(const ScoredOperations& other) const {
        return score < other.score;
    }

    /// Compare greater than based only on score
    inline bool operator>(const ScoredOperations& other) const {
        return score > other.score;
    }
    
    /// Make a match operation
    inline static ScoredOperations match(int score, int count) {
        return ScoredOperations(score, Operations::match(count));
    }

    /// Make a mismatch operation
    inline static ScoredOperations mismatch(int score, int count) {
        return ScoredOperations(score, Operations::mismatch(count));
    }

    /// Make a gap open operation
    inline static ScoredOperations open(int score, int count) {
        return ScoredOperations(score, Operations::open(count));
    }

    /// Make a gap extend operation
    inline static ScoredOperations extend(int score, int count) {
        return ScoredOperations(score, Operations::extend(count));
    }

    /// Make an unknown/not yet determined operation
    inline static ScoredOperations unknown(int score, int count) {
        return ScoredOperations(score, Operations::unknown(count));
    }

    /// Make a sentinel impossible value
    inline static ScoredOperations impossible() {
        return ScoredOperations(std::numeric_limits<int>::min(), Operations());
    }

    /// Allow conversion to an integer
    inline operator int() const {
        return score;
    }
};

/// Write a score and its operations to a stream
ostream& operator<<(ostream& out, const ScoredOperations& operations);

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
    /// Get the alignment score of the anchor (and the operations involved) 
    inline const ScoredOperations& score() const {
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
    inline ZipCodeDecoder* start_hint() const {
        return start_decoder;
    }

    /// Get the graph distance from wherever the start hint is positioned back
    /// to the actual start of the anchor.
    inline size_t start_hint_offset() const {
        return start_offset;
    }
    
    /// Get the distance-finding hint information (i.e. "zip code") for
    /// accelerating distance queries from the end of this anchor, or null if
    /// none is set.
    inline ZipCodeDecoder* end_hint() const {
        return end_decoder;
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

    // Construction
    
    /// Compose a read start position, graph start position, and match length into an Anchor.
    /// Can also bring along a distance hint and a seed number.
    inline Anchor(size_t read_start, const pos_t& graph_start, size_t length, size_t margin_before, size_t margin_after, const ScoredOperations& score, size_t seed_number = std::numeric_limits<size_t>::max(), ZipCodeDecoder* hint = nullptr, size_t hint_start = 0) : start(read_start), size(length), margin_before(margin_before), margin_after(margin_after), start_pos(graph_start), end_pos(advance(graph_start, length)), points(score), start_seed(seed_number), end_seed(seed_number), start_decoder(hint), end_decoder(hint), start_offset(hint_start), end_offset(length - hint_start), seed_length(margin_before + length + margin_after) {
        // Nothing to do!
    }
    
    /// Compose two Anchors into an Anchor that represents coming in through
    /// the first one and going out through the second, like a tunnel. Useful
    /// for representing chains as chainable items.
    inline Anchor(const Anchor& first, const Anchor& last, size_t extra_margin_before, size_t extra_margin_after, const ScoredOperations& score) : start(first.read_start()), size(last.read_end() - first.read_start()), margin_before(first.margin_before + extra_margin_before), margin_after(last.margin_after + extra_margin_after), start_pos(first.graph_start()), end_pos(last.graph_end()), points(score), start_seed(first.seed_start()), end_seed(last.seed_end()), start_decoder(first.start_hint()), end_decoder(last.end_hint()), start_offset(first.start_offset), end_offset(last.end_offset), seed_length((first.base_seed_length() + last.base_seed_length()) / 2) {
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
    ScoredOperations points;
    size_t start_seed;
    size_t end_seed;
    ZipCodeDecoder* start_decoder;
    ZipCodeDecoder* end_decoder;
    size_t start_offset;
    size_t end_offset;
    size_t seed_length;
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


    /// Construct a default, unset TracedScore
    inline TracedScore(): _score(ScoredOperations()), _source(nowhere()) {
        // Nothing to do!
    }

    /// Construct a TracedScore from a score and a source
    inline TracedScore(const ScoredOperations& score, size_t source): _score(score), _source(source) {
        // Nothing to do
    }

    // Make movable and copyable
    TracedScore(const TracedScore& other) = default;
    TracedScore(TracedScore&& other) = default;
    TracedScore& operator=(const TracedScore& other) = default;
    TracedScore& operator=(TracedScore&& other) = default;

    
    /// What's the default value for an empty table cell? Syntactic sugar to
    /// make it clearer when we mean an unset value.
    inline static TracedScore unset() {
        return TracedScore();
    }
    
    /// Max in a score from a DP table. If it wins, record provenance.
    void max_in(const vector<TracedScore>& options, size_t option_number);
    
    /// Get a score from a table of scores and record provenance in it.
    static TracedScore score_from(const vector<TracedScore>& options, size_t option_number);
    
    /// Add (or remove) points along a route to somewhere, as part of an operation. Return a modified copy.
    TracedScore add(const ScoredOperations& adjustment) const;
    
    /// Compare for equality.
    /// Only score and source matter for equality and comparison; the oprtation
    /// totals just ride along.
    inline bool operator==(const TracedScore& other) const {
        return score() == other.score() && source() == other.source();
    }
    
    /// Compare for inequality
    inline bool operator!=(const TracedScore& other) const {
        return !(*this == other);
    }
    
    /// Compare for less-than
    inline bool operator<(const TracedScore& other) const {
        return score() < other.score() || (score() == other.score() && source() < other.source());
    }
    
    /// Compare for greater-than
    inline bool operator>(const TracedScore& other) const {
        return score() > other.score() || (score() == other.score() && source() > other.source());
    }
    
    /// Subtraction to yield a difference in points and operations
    inline ScoredOperations operator-(const TracedScore& other) const {
        return score() - other.score();
    }
    
    /// Get the score value and associated operations
    inline const ScoredOperations& score() const {
        return _score;
    }

    /// Get the source index
    inline size_t source() const {
        return _source;
    }

    

private:

    /// Number of points and the operations they came from
    ScoredOperations _score;
    /// Index of source score among possibilities/traceback pointer
    size_t _source;
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
 * with their penalties from the optimal score (and the operation count
 * deltas). The best_past_ending_score_ever is *not* always the source of the
 * first traceback, if there is a tie.
 *
 * Tracebacks are constrained to be nonoverlapping by stopping each traceback
 * when the optimum place to come from has already been used. The second-best
 * place to come from is *not* considered. It might be possible that two
 * returned tracebacks could be pasted together to get a higher score, but it
 * won't be possible to recombine two tracebacks to get a higher score; no
 * edges followed between items will ever need to be cut.
 */
vector<pair<vector<size_t>, ScoredOperations>> chain_items_traceback(const vector<TracedScore>& chain_scores,
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
vector<pair<ScoredOperations, vector<size_t>>> find_best_chains(const VectorView<Anchor>& to_chain,
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
pair<ScoredOperations, vector<size_t>> find_best_chain(const VectorView<Anchor>& to_chain,
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
