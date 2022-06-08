#ifndef VG_GBWT_EXTENDER_HPP_INCLUDED
#define VG_GBWT_EXTENDER_HPP_INCLUDED

/** \file 
 * Haplotype-consistent seed extension in GBWTGraph.
 */

#include <functional>
#include <unordered_set>

#include "aligner.hpp"

#include <gbwtgraph/cached_gbwtgraph.h>

namespace vg {

//------------------------------------------------------------------------------

/**
 * A result of the gapless extension of a seed.
 * - The extension is a path starting from offset 'offset' of node path.front().
 * - Search state 'state' corresponds to the path.
 * - The extension covers semiopen interval [read_interval.first, read_interval.second)
 *   of the read.
 * - Vector 'mismatch_positions' contains the mismatching read positions in sorted order.
 * - 'score' is an alignment score (bigger is better).
 * - Flags 'left_full' and 'right_full' indicate whether the extension covers the
 *   start/end of the read.
 */
struct GaplessExtension
{
    typedef std::pair<handle_t, int64_t> seed_type; // (handle, read_offset - node_offset).

    // In the graph.
    std::vector<handle_t>     path;
    size_t                    offset;
    gbwt::BidirectionalState  state;

    // In the read.
    std::pair<size_t, size_t> read_interval;
    std::vector<size_t>       mismatch_positions;

    // Alignment properties.
    int32_t                   score;
    bool                      left_full, right_full;

    // For internal use.
    bool                      left_maximal, right_maximal;
    uint32_t                  internal_score; // Total number of mismatches.
    uint32_t                  old_score;      // Mismatches before the current flank.

    /// Length of the extension.
    size_t length() const { return this->read_interval.second - this->read_interval.first; }

    /// Is the extension empty?
    bool empty() const { return (this->length() == 0); }

    /// Is the extension a full-length alignment?
    bool full() const { return (this->left_full & this->right_full); }

    /// Is the extension an exact match?
    bool exact() const { return this->mismatch_positions.empty(); }

    /// Number of mismatches in the extension.
    size_t mismatches() const { return this->mismatch_positions.size(); }

    /// Does the extension contain the seed?
    bool contains(const HandleGraph& graph, seed_type seed) const;

    /// Return the starting position of the extension.
    Position starting_position(const HandleGraph& graph) const;

    /// Return the position after the extension.
    Position tail_position(const HandleGraph& graph) const;

    /// Return the node offset after the extension.
    size_t tail_offset(const HandleGraph& graph) const;

    /// Number of shared (read position, graph position) pairs in the extensions.
    size_t overlap(const HandleGraph& graph, const GaplessExtension& another) const;
    
    /// Take a prefix of the extension as a new GaplessExtension object.
    GaplessExtension prefix(size_t prefix_length) const;

    /// Convert the extension into a Path.
    Path to_path(const HandleGraph& graph, const std::string& sequence) const;

    /// For priority queues.
    bool operator<(const GaplessExtension& another) const {
        return (this->score < another.score);
    }

    /// Two extensions are equal if the same read interval matches the same search state
    /// with the same node offset.
    bool operator==(const GaplessExtension& another) const {
        return (this->read_interval == another.read_interval && this->state == another.state && this->offset == another.offset);
    }

    /// Two extensions are not equal if the state, the read interval, or the node offset is different.
    bool operator!=(const GaplessExtension& another) const {
        return !(this->operator==(another));
    }
};

//------------------------------------------------------------------------------

/**
 * An utility class that masks all characters except the specified ones with X,
 * which is assumed to not exist in the alignment target.
 */
class ReadMasker {
public:
    /// Creates a new ReadMasker that masks all characters except the specified
    /// ones.
    explicit ReadMasker(const std::string& valid_chars);

    /// Applies the mask to the given sequence.
    void operator()(std::string& sequence) const;

private:
    std::vector<char> mask;
};

//------------------------------------------------------------------------------

/**
 * A class that supports haplotype-consistent seed extension using GBWTGraph. Each seed
 * is a pair of matching read/graph positions and each extension is a gapless alignment
 * of an interval of the read to a haplotype.
 * A cluster is an unordered set of distinct seeds. Seeds in the same node with the same
 * (read_offset - node_offset) difference are considered equivalent.
 * GaplessExtender also needs an Aligner object for scoring the extension candidates.
 */
class GaplessExtender {
public:
    typedef GaplessExtension::seed_type seed_type;
    typedef pair_hash_set<seed_type>    cluster_type;

    /// The default value for the maximum number of mismatches.
    constexpr static size_t MAX_MISMATCHES = 4;

    /// Two full-length alignments are distinct, if the fraction of overlapping
    /// position pairs is at most this.
    constexpr static double OVERLAP_THRESHOLD = 0.8;

    /// Create an empty GaplessExtender.
    GaplessExtender();

    /// Create a GaplessExtender using the given GBWTGraph and Aligner objects.
    GaplessExtender(const gbwtgraph::GBWTGraph& graph, const Aligner& aligner);

    /// Convert (graph position, read offset) to a seed.
    static seed_type to_seed(pos_t pos, size_t read_offset) {
        return seed_type(gbwtgraph::GBWTGraph::node_to_handle(gbwt::Node::encode(id(pos), is_rev(pos))),
                         static_cast<int64_t>(read_offset) - static_cast<int64_t>(offset(pos)));
    }

    /// Get the graph position from a seed.
    static pos_t get_pos(seed_type seed) {
        gbwt::node_type node = gbwtgraph::GBWTGraph::handle_to_node(seed.first);
        return make_pos_t(gbwt::Node::id(node), gbwt::Node::is_reverse(node), get_node_offset(seed));
    }

    /// Get the handle from a seed.
    static handle_t get_handle(seed_type seed) {
        return seed.first;
    }

    /// Get the node offset from a seed.
    static size_t get_node_offset(seed_type seed) {
        return (seed.second < 0 ? -(seed.second) : 0);
    }

    /// Get the read offset from a seed.
    static size_t get_read_offset(seed_type seed) {
        return (seed.second < 0 ? 0 : seed.second);
    }

    /**
     * Find the highest-scoring extension for each seed in the cluster.
     * If there is a full-length extension with at most max_mismatches
     * mismatches, sort them in descending order by score and return the
     * best non-overlapping full-length extensions. Two extensions overlap
     * if the fraction of identical base mappings is greater than
     * overlap_threshold.
     * If there are no good enough full-length extensions, trim the
     * extensions to maximize the score and remove duplicates. In this
     * case, the extensions are sorted by read interval.
     * Use full_length_extensions() to determine the type of the returned
     * extension set.
     * The sequence that will be aligned is passed by value. All non-ACGT
     * characters are masked with character X, which should not match any
     * character in the graph.
     * Allow any number of mismatches in the initial node, at least
     * max_mismatches mismatches in the entire extension, and at least
     * max_mismatches / 2 mismatches on each flank.
     * Use the provided CachedGBWTGraph or allocate a new one.
     */
    std::vector<GaplessExtension> extend(cluster_type& cluster, std::string sequence, const gbwtgraph::CachedGBWTGraph* cache = nullptr, size_t max_mismatches = MAX_MISMATCHES, double overlap_threshold = OVERLAP_THRESHOLD) const;

    /**
     * Determine whether the extension set contains non-overlapping
     * full-length extensions sorted in descending order by score. Use
     * the same value of max_mismatches as in extend().
     */
    static bool full_length_extensions(const std::vector<GaplessExtension>& result, size_t max_mismatches = MAX_MISMATCHES);

    const gbwtgraph::GBWTGraph* graph;
    const Aligner*              aligner;
    ReadMasker                  mask;
};

//------------------------------------------------------------------------------

/**
 * An alignment found by WFAAligner, consisting of a sequence of nodes, a sequence
 * of edits, and a starting offset in the initial node. The alignment score is
 * computed using the match / mismatch / gap open / gap extend parameters in the
 * Aligner object given to the WFAAligner. Full-length bonuses are not used.
 *
 * Note: The aligner merges consecutive edit operations of the same type. Hence an
 * edit may span multiple nodes.
 *
 * This struct is an aggregate and can be aggregate-initialized with a
 * brace-enclosed initializer list.
 */
struct WFAAlignment {
    enum Edit { match, mismatch, insertion, deletion };
    
    /// Generate a WFAAlignment from a GaplessExtension. This can't be a
    /// constructor because then WFAAlignment wouldn't be able to be
    /// aggregate-initialized. 
    static WFAAlignment from_extension(const GaplessExtension& extension);

    /// Sequence of oriented nodes.
    std::vector<handle_t> path;

    /// Sequence of edit operations and their lengths.
    std::vector<std::pair<Edit, uint32_t>> edits;

    /// Starting offset in the initial node.
    uint32_t node_offset = 0;

    /// Starting offset in the sequence.
    uint32_t seq_offset = 0;

    /// Length of the alignment in the sequence.
    uint32_t length = 0;

    /// Alignment score.
    int32_t score = 0;

    /// Returns true if the alignment is empty.
    bool empty() const { return (this->length == 0); }

    /// Returns true if the alignment is an insertion without a location.
    bool unlocalized_insertion() const;

    /// Computes the node offset after the alignment in the final node.
    uint32_t final_offset(const gbwtgraph::GBWTGraph& graph) const;

    /// Transforms the alignment to the other strand.
    void flip(const gbwtgraph::GBWTGraph& graph, const std::string& sequence);

    /// Appends an edit operation, merging it with the latest edit if possible.
    /// Ignores empty edits.
    void append(Edit edit, uint32_t length);
    
    /// Concatenate another WFAAlignment onto this one, assuming a single
    /// shared match base. Either may be empty.
    void join_on_shared_match(const WFAAlignment& second, int match_score);
    
    /// Convert the WFAAlignment into a Path.
    Path to_path(const HandleGraph& graph, const std::string& sequence) const;

    /// Prints some debug information about the alignment.
    std::ostream& print(const gbwtgraph::GBWTGraph& graph, std::ostream& out) const;
};

//------------------------------------------------------------------------------

/**
 * A class that supports haplotype-consistent seed extension in a GBWTGraph using the
 * WFA algorithm:
 *
 *   Marco-Sola, Moure, Moreto, Espinosa: Fast gap-affine pairwise alignment using the
 *   wavefront algorithm. Bioinformatics, 2021.
 *
 * The algorithm either tries to connect two seeds or extends a seed to the start/end
 * of the read.
 *
 * WFAExtender also needs an Aligner object for scoring the extension candidates.
 * While VG wants to maximize a four-parameter alignment score, WFA minimizes a
 * three-parameter score. We use the conversion between the parameters from:
 *
 *   Eizenga, Paten: Improving the time and space complexity of the WFA algorithm
 *   and generalizing its scoring. bioRxiv, 2022.
 *
 * VG scores a gap of length `n` as `gap_open + (n - 1) * gap_extend`, while WFA
 * papers use `gap_open + n * gap_extend`. Hence we use `gap_open - gap_extend` as
 * the effective four-parameter gap open score inside this aligner.
 *
 * NOTE: Most internal arithmetic operations use 32-bit integers.
 */
class WFAExtender {
public:
    /// Create an empty WFAExtender.
    WFAExtender();

    /// Create a WFAExtender using the given GBWTGraph and Aligner objects.
    /// TODO: Provide an error model for specifying bounds for acceptable
    /// alignment scores.
    WFAExtender(const gbwtgraph::GBWTGraph& graph, const Aligner& aligner);

    /**
     * Align the sequence to a haplotype between the two graph positions.
     *
     * The sequence that will be aligned is passed by value. All non-ACGT
     * characters are masked with character X, which should not match any
     * character in the graph.
     *
     * Returns an empty alignment if there is no alignment with an acceptable
     * score.
     */
    WFAAlignment connect(std::string sequence, pos_t from, pos_t to) const;

    /**
     * A special case of connect() for aligning the sequence to a haplotype
     * starting at the given position. If there is no alignment for the
     * entire sequence with an acceptable score, returns the highest-scoring
     * partial alignment.
     *
     * NOTE: This creates a suffix of the alignment by aligning a prefix of
     * the sequence.
     * TODO: Should we use full-length bonuses?
     */
    WFAAlignment suffix(const std::string& sequence, pos_t from) const;

    /**
     * A special case of connect() for aligning the sequence to a haplotype
     * ending at the given position. If there is no alignment for the entire
     * sequence with an acceptable score, returns the highest-scoring partial
     * alignment.
     *
     * NOTE: This creates a prefix of the alignment by aligning a suffix of
     * the sequence.
     * TODO: Should we use full-length bonuses?
     */
    WFAAlignment prefix(const std::string& sequence, pos_t to) const;

    const gbwtgraph::GBWTGraph* graph;
    ReadMasker                  mask;
    const Aligner*              aligner;

    /// TODO: Remove when unnecessary.
    bool debug = false;
};

//------------------------------------------------------------------------------

} // namespace vg

#endif // VG_GBWT_EXTENDER_HPP_INCLUDED
