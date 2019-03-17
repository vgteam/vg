#ifndef VG_GAPLESS_EXTENDER_HPP_INCLUDED
#define VG_GAPLESS_EXTENDER_HPP_INCLUDED

/** \file 
 * Haplotype-consistent gapless seed extension.
 */

#include "gbwt_helper.hpp"

namespace vg {

//------------------------------------------------------------------------------

/**
 * A class that supports haplotype-consistent seed extension using GBWTGraph.
 * The extension finds, for a given cluster of matching sequence/graph positions,
 * the best haplotype-consistent gapless alignment starting from any pair of
 * positions in the cluster. The search can be constrained by giving a maximum
 * number of mismatches.
 */
class GaplessExtender {
public:
    /// The default value for the maximum number of mismatches.
    constexpr static size_t MAX_MISMATCHES = 4;

    /// Create an empty GaplessExtender.
    GaplessExtender();

    /// Create a GaplessExtender using the given GBWTGraph and set the maximum
    /// number of mismatches allowed.
    explicit GaplessExtender(const GBWTGraph& graph, size_t max_mismatches = MAX_MISMATCHES);

    /// Extend the exact match hits against the given input sequence into a Path representing
    /// the best alignment to the graph and a score value. The cluster can be in an arbitrary
    /// order, but it will be sorted during the call. Each pair in the cluster consists of
    /// matching sequence/graph positions. Some positions may occur multiple times, and the
    /// matches in the cluster may agree or conflict.
    std::pair<Path, size_t> extend_seeds(std::vector<std::pair<size_t, pos_t>>& cluster, const std::string& sequence);

    const GBWTGraph* graph;
    size_t           max_mismatches;
};

//------------------------------------------------------------------------------

} // namespace vg

#endif // VG_GAPLESS_EXTENDER_HPP_INCLUDED
