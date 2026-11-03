#ifndef VG_HAPLOTYPE_DEPTH_HPP_INCLUDED
#define VG_HAPLOTYPE_DEPTH_HPP_INCLUDED

/**
 * \file haplotype_depth.hpp
 *
 * How many haplotypes read each base of a graph, read off along any path.
 */

#include "walk_interval.hpp"

#include <gbwtgraph/gbz.h>

#include <cstdint>
#include <functional>
#include <iosfwd>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace vg {

/**
 * A set of regions on the full paths of a graph, by path name.
 *
 * A path is named as for_each_path_start() names it, so a region on a
 * haplotype is in coordinates of its assembly contig. Each path's regions are
 * sorted and merged.
 */
class PathRegions {
public:
    /// A stretch [first, second) of a full path.
    using Region = std::pair<uint64_t, uint64_t>;

    /**
     * Read the regions in a BED file. Columns past the 3rd are ignored.
     *
     * @throws std::runtime_error on a malformed line.
     */
    explicit PathRegions(std::istream& bed);

    /// Return the regions of the named path, sorted and merged, or null if it has none.
    const std::vector<Region>* regions_of(const std::string& name) const;
    /// Return the names of the paths that have regions.
    std::vector<std::string> path_names() const;

private:
    std::unordered_map<std::string, std::vector<Region>> regions;
};

/**
 * A tally of how many haplotypes read each base of each node of a graph.
 *
 * It counts distinct haplotypes, by sample and haplotype number, or the visits
 * they make, which exceed the haplotypes where one reads a node more than
 * once. Haplotype paths count, and reference paths if asked. Gref paths and
 * generic paths never count. With regions, a path counts only its bases
 * inside its regions, so a path without regions counts nothing.
 *
 * Most nodes are read whole by every haplotype that reads them and hold one
 * count. A node that regions cut holds a count for each stretch of its bases.
 */
class HaplotypeDepth {
public:
    struct Options {
        /// Count visits instead of distinct haplotypes.
        bool count_visits = false;
        /// Count reference paths as haplotypes too.
        bool count_references = false;
        /// Count only the bases of each path inside its regions. Null counts every base.
        const PathRegions* regions = nullptr;
    };

    /// Count the haplotypes of the GBZ, with OpenMP threads.
    HaplotypeDepth(const gbwtgraph::GBZ& gbz, const Options& options);

    /**
     * Call visit(start, end, depth) for each maximal run of bases of equal
     * depth along the path, in order, including runs of depth 0.
     *
     * Coordinates are on the full path, which starts offset bases before the
     * path's first base.
     */
    void for_each_run(const PathSteps& path, uint64_t offset,
                      const std::function<void(uint64_t, uint64_t, uint32_t)>& visit) const;

    /// Return the bases that the counted haplotypes read on the node: its
    /// depth summed over its bases.
    uint64_t bases_read(handlegraph::nid_t id, uint32_t length) const;

    /// Return the number of haplotypes counted, with or without regions.
    size_t haplotype_count() const { return haplotypes; }
    /// Return the number of counted haplotypes that have no regions, and so add nothing.
    size_t haplotypes_without_regions() const { return without_regions; }
    /// Return the names of paths with regions that are not counted paths in the graph.
    const std::vector<std::string>& uncounted_region_paths() const { return unmatched_regions; }

private:
    /// Bases [from, to) of a node, in its forward orientation.
    struct NodeStretch {
        handlegraph::nid_t node;
        uint32_t from;
        uint32_t to;
    };
    /**
     * Add one to the count of each base in the stretches that one haplotype reads.
     *
     * A stretch of a whole node adds to its count here. A stretch of part of a
     * node is appended to cut_stretches. With merge, overlapping stretches of a
     * node count once, so the count is of haplotypes and not of visits.
     */
    void add_stretches(std::vector<NodeStretch> stretches, const gbwtgraph::GBWTGraph& graph, bool merge,
                       std::vector<NodeStretch>& cut_stretches);

    handlegraph::nid_t min_id = 0;
    /// The count of each node read whole, by node ID minus min_id.
    std::vector<uint32_t> whole;
    /// For each node that regions cut, where its count changes along its
    /// forward bases: the offset and the change, sorted by offset.
    std::unordered_map<handlegraph::nid_t, std::vector<std::pair<uint32_t, int32_t>>> cuts;

    size_t haplotypes = 0;
    size_t without_regions = 0;
    std::vector<std::string> unmatched_regions;
};

}

#endif
