#ifndef VG_MOSAIC_WRITER_HPP_INCLUDED
#define VG_MOSAIC_WRITER_HPP_INCLUDED

#include <atomic>
#include <cstddef>

namespace vg {

/// Counters for the mosaic writer. A member of each VCFOutputCaller, so that runs count
/// separately; `mutable` there because the writing paths are const.
struct MosaicCounters {
    /// Runs with no position to walk from. Panel haplotypes are often fragments, so this is
    /// reported rather than expected to be zero.
    std::atomic<size_t> unwalkable{0};
    /// Of those, the ones that are only a head: the run could be walked from a later site, so the
    /// walkable rest is written separately.
    std::atomic<size_t> head_clipped{0};
    /// Run boundaries across which the first run's haplotype could be followed, and those across
    /// which it could not, which a reference fill or a new fragment has to cover.
    std::atomic<size_t> extended{0}, gap_left{0}, patched{0};
    /// Rows whose own haplotype does not span them, rewritten as a reference substitution.
    std::atomic<size_t> row_to_ref{0};
    /// Run boundaries between a parent and a child snarl, at the child's own boundary nodes.
    std::atomic<size_t> nested_enter{0}, nested_leave{0};
    /// Rows the current direction could not walk but the other direction could, at an inversion.
    std::atomic<size_t> direction_broken{0}, extended_left{0};
};

}

#endif
