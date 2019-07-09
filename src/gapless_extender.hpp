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
 * - The extension covers semiopen interval 'read_interval' of the read.
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
    uint32_t                  flank_score;    // Mismatches in the current flank.

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
    seed_type to_seed(pos_t pos, size_t read_offset) const;

    /**
     * Find the best full-length alignment for the sequence within the cluster with
     * at most max_mismatches mismatches.
     * If that is not possible, find the set of highest-scoring maximal extensions
     * of the seeds, allowing any number of mismatches in the seed node and
     * max_mismatches / 2 mismatches on each flank. Flanks may have more mismatches
     * if it does not bring the total beyond max_mismatches. Then try to trim
     * mismatches from the flanks if it improves the score.
     * The extensions are sorted by their coordinates in the sequence.
     */
    std::vector<GaplessExtension> extend(cluster_type& cluster, const std::string& sequence, size_t max_mismatches = MAX_MISMATCHES) const;

    /**
     * Try to improve the score of each extension by trimming mismatches from the flanks.
     * Do not trim full-length alignments with <= max_mismatches mismatches.
     * Note that extend() already calls this.
     */
    void trim(std::vector<GaplessExtension>& extensions, size_t max_mismatches = MAX_MISMATCHES) const;

    const GBWTGraph* graph;
    const Aligner*   aligner;
};

//------------------------------------------------------------------------------

} // namespace vg

#endif // VG_GAPLESS_EXTENDER_HPP_INCLUDED
