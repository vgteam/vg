#ifndef VG_GAPLESS_EXTENDER_HPP_INCLUDED
#define VG_GAPLESS_EXTENDER_HPP_INCLUDED

/** \file 
 * Haplotype-consistent gapless seed extension.
 */

#include <functional>

#include "gbwt_helper.hpp"

namespace vg {

//------------------------------------------------------------------------------

/**
 * A result of the gapless extension of a seed.
 * - 'path' covers 'core_interval' of the read, starting from 'offset' in the initial
 *   node.
 * - The intervals are semi-open [first, second) intervals.
 * - 'state' is the search state corresponding to 'path'.
 * - If 'is_full_alignment' is set, the path is a full-length alignment and may contain
 *   mismatches. Otherwise it is a maximal unambiguous exact extension of a seed.
 * - If the path is not a full-length alignment, 'flanked_interval' is a read interval
 *   extending 'core_interval'. The extension may contain mismatches, but it starts and
 *   ends with a match.
 * - The mismatching read positions are stored in 'mismatch_positions' in sorted order.
 */
struct GaplessExtension
{
    std::vector<handle_t>     path;
    size_t                    offset;

    gbwt::BidirectionalState  state;
    std::pair<size_t, size_t> core_interval;
    bool                      is_full_alignment;

    std::pair<size_t, size_t> flanked_interval;
    std::vector<size_t>       mismatch_positions;

    /// Length of the core interval.
    size_t core_length() const { return this->core_interval.second - this->core_interval.first; }

    /// Length of the flanked interval.
    size_t flanked_length() const { return this->flanked_interval.second - this->flanked_interval.first; }

    /// Is the extension empty?
    bool empty() const { return (this->core_length() == 0); }

    /// Is the extension a full-length alignment?
    bool full() const { return this->is_full_alignment; }

    /// Is the extension an exact match?
    bool exact() const { return this->mismatch_positions.empty(); }

    /// Number of mismatches in the extension.
    size_t mismatches() const { return this->mismatch_positions.size(); }

    /// Return the starting position of the core interval.
    Position starting_position(const GBWTGraph& graph) const;

    /// Return the position after the core interval.
    Position tail_position(const GBWTGraph& graph) const;

    /// Return the node offset after the core interval
    size_t tail_offset(const GBWTGraph& graph) const;

    /// Convert the extension into a Path.
    Path to_path(const GBWTGraph& graph, const std::string& sequence) const;
};

//------------------------------------------------------------------------------

/**
 * A class that supports haplotype-consistent seed extension using GBWTGraph. Each seed
 * is a pair of matching read/graph positions and each extension is a gapless alignment
 * of an interval of the read to a haplotype.
 * A cluster is a vector of distinct seeds. The cluster may initially be in an arbitrary
 * order, but it will be sorted during the extension. All seeds in a cluster should
 * correspond to the same alignment or positions near it.
 */
class GaplessExtender {
public:
    typedef std::pair<size_t, pos_t> seed_type;
    typedef std::vector<seed_type>   cluster_type;

    /// The default value for the maximum number of mismatches.
    constexpr static size_t MAX_MISMATCHES = 4;

    /// Create an empty GaplessExtender.
    GaplessExtender();

    /// Create a GaplessExtender using the given GBWTGraph.
    explicit GaplessExtender(const GBWTGraph& graph);

    /**
     * 1. Call extend_seeds(). If there is a full-length alignment, return the result.
     * 2. Call maximal_extensions().
     * 3. Call extend_flanks() with max_mismatches / 2 mismatches.
     * The extensions are sorted by (core_interval.first, core_interval.second).
     */
    std::vector<GaplessExtension> extend(cluster_type& cluster, const std::string& sequence, size_t max_mismatches = MAX_MISMATCHES) const;

    /**
     * Find a full-length extension of the seeds with up to 'max_mismatches' mismatches.
     * Return an alignment with the fewest number of mismatches or an empty extension if
     * no full-length alignment exists.
     */
    GaplessExtension extend_seeds(cluster_type& cluster, const std::string& sequence, size_t max_mismatches = MAX_MISMATCHES, bool cluster_is_sorted = false) const;

    /**
     * Find the maximal unambiguous extension for each seed. Returns the set of distinct
     * extensions. A maximal unambiguous extension ends when further extensions either
     * contain mismatches or branch with the same character.
     */
    std::vector<GaplessExtension> maximal_extensions(cluster_type& cluster, const std::string& sequence, bool cluster_is_sorted = false) const;

    /**
     * Extends the flanks of each unambiguous extension with up to 'max_mismatches'
     * mismatches in each direction. Trims the mismatches if the flanks start/end with
     * them. Updates 'flanked_interval' to match the longest trimmed interval.
     * Note that the extensions must be exact, non-empty, and non-full.
     */
    void extend_flanks(std::vector<GaplessExtension>& extensions, const std::string& sequence, size_t max_mismatches = MAX_MISMATCHES / 2) const;

    const GBWTGraph* graph;
};

//------------------------------------------------------------------------------

} // namespace vg

#endif // VG_GAPLESS_EXTENDER_HPP_INCLUDED
