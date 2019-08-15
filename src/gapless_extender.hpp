#ifndef VG_GAPLESS_EXTENDER_HPP_INCLUDED
#define VG_GAPLESS_EXTENDER_HPP_INCLUDED

/** \file 
 * Haplotype-consistent gapless seed extension.
 */

#include <functional>

#include "aligner.hpp"
#include "gbwt_helper.hpp"
#include "hash_map.hpp"

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

    /// Return the starting position of the extension.
    Position starting_position(const GBWTGraph& graph) const;

    /// Return the position after the extension.
    Position tail_position(const GBWTGraph& graph) const;

    /// Return the node offset after the extension.
    size_t tail_offset(const GBWTGraph& graph) const;

    /// Convert the extension into a Path.
    Path to_path(const GBWTGraph& graph, const std::string& sequence) const;

    /// For priority queues.
    bool operator<(const GaplessExtension& another) const {
        return (this->score < another.score);
    }
};

//------------------------------------------------------------------------------

/**
 * A class that supports haplotype-consistent seed extension using GBWTGraph. Each seed
 * is a pair of matching read/graph positions and each extension is a gapless alignment
 * of an interval of the read to a haplotype.
 * A cluster is an unordered set of distinct seeds. Seeds in the same node with the same
 * (read_offset - node_offset) difference are considered equivalent. All seeds in a cluster
 * should correspond to the same alignment or positions near it.
 * GaplessExtender also needs an Aligner object for scoring the extension candidates.
 */
class GaplessExtender {
public:
    typedef std::pair<handle_t, int64_t> seed_type; // (handle, read_offset - node_offset).
    typedef pair_hash_set<seed_type>     cluster_type;

    /// The default value for the maximum number of mismatches.
    constexpr static size_t MAX_MISMATCHES = 4;

    /// Create an empty GaplessExtender.
    GaplessExtender();

    /// Create a GaplessExtender using the given GBWTGraph and Aligner objects.
    explicit GaplessExtender(const GBWTGraph& graph, const Aligner& aligner);

    /// Convert (graph position, read offset) to a seed.
    static seed_type to_seed(pos_t pos, size_t read_offset) {
        return seed_type(GBWTGraph::node_to_handle(gbwt::Node::encode(id(pos), is_rev(pos))),
                         static_cast<int64_t>(read_offset) - static_cast<int64_t>(offset(pos)));
    }

    /// Get the graph position from a seed.
    static pos_t get_pos(seed_type seed) {
        gbwt::node_type node = GBWTGraph::handle_to_node(seed.first);
        return make_pos_t(gbwt::Node::id(node), gbwt::Node::is_reverse(node), get_node_offset(seed));
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
     * Find the best full-length alignment for the sequence within the cluster with
     * at most max_mismatches mismatches.
     * If that is not possible, find the set of highest-scoring maximal extensions
     * of the seeds, allowing any number of mismatches in the seed node and
     * max_mismatches / 2 mismatches on each flank. Flanks may have more mismatches
     * if it does not bring the total beyond max_mismatches. Then call trim() if
     * trim_extensions i set.
     * The extensions are sorted by their coordinates in the sequence.
     */
    std::vector<GaplessExtension> extend(cluster_type& cluster, const std::string& sequence, size_t max_mismatches = MAX_MISMATCHES, bool trim_extensions = true) const;

    /**
     * Try to improve the score of each extension by trimming mismatches from the flanks.
     * Do not trim full-length alignments with <= max_mismatches mismatches.
     * Use the provided CachedGBWT or allocate a new one.
     * Note that extend() already calls this by default.
     */
    void trim(std::vector<GaplessExtension>& extensions, size_t max_mismatches = MAX_MISMATCHES, const gbwt::CachedGBWT* cache = nullptr) const;

    const GBWTGraph* graph;
    const Aligner*   aligner;
};

//------------------------------------------------------------------------------

} // namespace vg

#endif // VG_GAPLESS_EXTENDER_HPP_INCLUDED
