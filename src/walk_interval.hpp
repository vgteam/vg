#ifndef VG_WALK_INTERVAL_HPP_INCLUDED
#define VG_WALK_INTERVAL_HPP_INCLUDED

/**
 * \file walk_interval.hpp
 *
 * Stretches of sequence along walks in a GBWTGraph, and the paths of its GBWT
 * as walks with coordinates.
 */

#include <gbwtgraph/gbwtgraph.h>

#include <cstdint>
#include <functional>
#include <iosfwd>
#include <string>
#include <string_view>
#include <vector>

namespace vg {

/**
 * A stretch of graph sequence.
 *
 * It holds the bases from offset to offset + length along the walk. The walk
 * visits only the nodes that hold those bases. So offset is less than the
 * length of the first node, and the last base is on the last node.
 */
struct WalkInterval {
    std::vector<handlegraph::handle_t> walk;
    /// Offset of the first base in the first node of the walk.
    size_t offset = 0;
    size_t length = 0;

    /// Return the bases, read along the walk.
    std::string sequence(const gbwtgraph::GBWTGraph& graph) const;
    /// Return the same bases read the other way, along the reversed walk with
    /// each handle flipped.
    WalkInterval reverse(const gbwtgraph::GBWTGraph& graph) const;
    /// Return bases [start, end) of this interval. Requires start < end <= length.
    WalkInterval subinterval(const gbwtgraph::GBWTGraph& graph, size_t start, size_t end) const;
    /// Return the total length of the nodes on the walk.
    size_t walk_length(const gbwtgraph::GBWTGraph& graph) const;
};

/**
 * Parse a walk written as in GAF, such as >1<2>3.
 *
 * Each visit is > (forward) or < (reverse), followed by a node ID.
 *
 * @throws std::runtime_error if the text is not such a walk, or if it visits a
 * node that is not in the graph.
 */
std::vector<handlegraph::handle_t> parse_walk(std::string_view text, const gbwtgraph::GBWTGraph& graph);

/// Return the walk written as parse_walk() reads it.
std::string walk_text(const std::vector<handlegraph::handle_t>& walk, const gbwtgraph::GBWTGraph& graph);

/// Write the interval as the 12 required columns of a GAF record. The query is
/// the interval's sequence under the given name, aligned exactly. The caller
/// writes any tags and the newline.
void write_gaf_columns(std::ostream& out, const gbwtgraph::GBWTGraph& graph, const std::string& name,
                       const WalkInterval& interval);

/**
 * Read a GAF line written by write_gaf_columns().
 *
 * @returns The interval, and stores the query name in name.
 * @throws std::runtime_error if the line is not a GAF record of a walk in the graph.
 */
WalkInterval parse_gaf_columns(const std::string& line, const gbwtgraph::GBWTGraph& graph, std::string& name);

/**
 * Call visit(path_id, name, start) for each path in the GBWT, in order.
 *
 * The name is the path's name without any subrange. The start is where the
 * path starts on the full path of that name, which is the start of its
 * subrange, or 0 without one.
 *
 * A haplotype path's start is its phase block, and its name leaves the phase
 * block out. Graphs built from assemblies, such as Minigraph-Cactus graphs,
 * store the offset of each haplotype fragment on its contig as the phase block.
 */
void for_each_path_start(const gbwt::GBWT& index,
                         const std::function<void(gbwt::size_type, const std::string&, uint64_t)>& visit);

/**
 * The steps of one GBWT path, with the path coordinate where each starts.
 *
 * It finds the walk interval of any stretch of the path.
 */
class PathSteps {
public:
    struct Step {
        /// Path coordinate of the first base on the step.
        uint64_t start;
        /// The handle the path visits.
        handlegraph::handle_t handle;
    };

    /// Extract the path with this GBWT path ID from the graph's GBWT.
    PathSteps(const gbwtgraph::GBWTGraph& graph, gbwt::size_type path_id);

    const std::vector<Step>& steps() const { return step_list; }
    /// Return the length of the path.
    uint64_t length() const { return path_length; }
    /// Return bases [start, end) of the path. Requires start < end <= length().
    WalkInterval interval(uint64_t start, uint64_t end) const;
    /// Return the base at the given path coordinate, as the path reads it.
    char base(const gbwtgraph::GBWTGraph& graph, uint64_t coord) const;

private:
    std::vector<Step> step_list;
    uint64_t path_length = 0;
};

}

#endif
