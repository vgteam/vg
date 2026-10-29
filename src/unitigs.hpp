#ifndef VG_UNITIGS_HPP_INCLUDED
#define VG_UNITIGS_HPP_INCLUDED

/**
 * \file unitigs.hpp
 *
 * Haplotype unitigs: a set of graph walks whose sequences contain every k-mer
 * spelled by the haplotypes in a GBZ, so that a linear k-mer scanner run over
 * them sees the same k-mers as one run over every haplotype.
 *
 * A k-mer here is a haplotype-consistent graph walk spelling k bases, so the
 * same sequence at two places in the graph is two k-mers. Two k-mers are
 * adjacent if some haplotype has them one base apart. A unitig is a maximal
 * chain of k-mers where each has exactly one successor, that successor has
 * exactly one predecessor, and no haplotype starts or ends in between. Every
 * k-mer is in exactly one unitig, and every unitig is a substring of every
 * haplotype that contains any of its k-mers, so its support is constant.
 *
 * With a gref (a set of node-disjoint paths in the graph, see gref.hpp), the
 * k-mers that lie along a single gref path are left out, because the gref
 * sequence already contains them. What remains are the junction tigs.
 *
 * The gref paths and the tigs together are a cover: sequences that a linear
 * tool can scan in place of the haplotypes. cover_lift.hpp translates what such
 * a tool finds back to the graph.
 *
 * Only the GBWTGraph interface is used here.
 */

#include <gbwtgraph/gbz.h>

#include <cstdint>
#include <functional>
#include <iosfwd>
#include <limits>
#include <string>
#include <unordered_map>
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

/// Where a walk interval lies on the gref. See GrefIndex::place().
struct GrefPlacement {
    /// Comma-separated PATH:START-END:STRAND:CIGAR segments (1-based,
    /// inclusive, CIGAR in path order), or *:<length>I if no base is on the gref.
    std::string text;
    /// The path of the first segment, or GrefIndex::NO_PATH if no base is on the gref.
    uint32_t path = std::numeric_limits<uint32_t>::max();
    /// The lowest coordinate of the first segment. Placements sort along the
    /// gref by path, then coord.
    uint64_t coord = 0;
    /// Number of bases on gref paths that the paths read forward.
    size_t forward_bases = 0;
    /// Number of bases on gref paths that the paths read in reverse.
    size_t reverse_bases = 0;
    /// True if the whole interval lies along one gref path, so that the gref
    /// path visits it.
    bool along_one_path = false;
};

/**
 * An index of where the gref paths lie in the graph.
 *
 * For each node, it holds the gref path that visits it and the path coordinate
 * of its first base. For each gref path, it holds the steps. With these, it
 * places walk intervals on the gref.
 */
class GrefIndex {
public:
    /// Index the GBWT paths whose name starts with the given prefix.
    GrefIndex(const gbwtgraph::GBZ& gbz, const std::string& name_prefix);

    static constexpr uint32_t NO_PATH = std::numeric_limits<uint32_t>::max();

    struct Position {
        uint32_t path = NO_PATH;
        /// The path reads the node in reverse.
        bool is_reverse = false;
        /// Path coordinate of the first base the path reads on the node.
        uint64_t offset = 0;
    };

    /// A stretch of a walk interval that lies contiguously along one gref path.
    struct Run {
        /// Bases [begin, end) of the interval.
        size_t begin, end;
        uint32_t path;
        /// True if path coordinates decrease along the interval.
        bool is_reverse;
        /// Path coordinate of the base at begin.
        uint64_t coord;

        uint64_t coord_at(size_t position) const {
            return is_reverse ? coord - (position - begin) : coord + (position - begin);
        }
        uint64_t last() const { return coord_at(end - 1); }
    };

    /// Where the node lies on the gref, or a Position with path == NO_PATH.
    const Position& find(handlegraph::nid_t id) const;

    size_t path_count() const { return paths.size(); }
    /// Path name without any subrange.
    const std::string& path_name(uint32_t path) const { return paths[path].name; }
    /// Amount to add to a path coordinate to get a coordinate on the full path.
    uint64_t path_shift(uint32_t path) const { return paths[path].shift; }
    /// Return the length of the path.
    uint64_t path_length(uint32_t path) const;
    /// Return the path with the given name, or NO_PATH. The name may have a
    /// subrange, as vg paths writes it. It may leave the subrange out if no
    /// other gref path has the same name.
    uint32_t find_path(const std::string& name) const;
    /// Return bases [start, end) of the path, in path coordinates. Requires
    /// start < end <= path_length(path).
    WalkInterval interval(uint32_t path, uint64_t start, uint64_t end) const;
    /// Base at the given path coordinate, as read by the path.
    char base(uint32_t path, uint64_t coord) const;

    /// Number of bases on gref paths.
    uint64_t total_length() const;
    /// Number of node visits that were ignored because the node was already on the gref.
    size_t repeated_visits() const { return repeats; }

    /// Return the runs of the interval, in order along it.
    std::vector<Run> runs(const WalkInterval& interval) const;
    /// Return where the interval lies on the gref. Bases between two runs on
    /// one path become edits. Bases between runs on different paths become an
    /// insertion at the end of the earlier segment.
    GrefPlacement place(const WalkInterval& interval) const;

private:
    struct Step {
        /// Path coordinate of the first base on the step.
        uint64_t start;
        /// The handle the path visits.
        handlegraph::handle_t handle;
    };
    /// A gref path.
    struct Path {
        /// Name without any subrange.
        std::string name;
        /// Where the path starts on the full path.
        uint64_t shift;
        std::vector<Step> steps;
    };

    const gbwtgraph::GBWTGraph& graph;
    handlegraph::nid_t min_id = 0;
    std::vector<Position> positions;
    std::vector<Path> paths;
    /// Index in paths of each gref path, by name with and without its subrange.
    std::unordered_map<std::string, uint32_t> path_ranks;
    size_t repeats = 0;
};

/// A walk interval in the cover, with its support and gref placement. Spell
/// its sequence with interval.sequence().
struct Tig {
    WalkInterval interval;
    /// Number of haplotype occurrences of each of its k-mers.
    size_t support = 0;
    /// Where it lies on the gref. Empty text without a gref.
    GrefPlacement placement;
};

struct UnitigStatistics {
    size_t unitigs = 0;
    size_t tigs = 0;
    size_t bases = 0;
    size_t kmers = 0;
    size_t longest = 0;
    size_t unanchored = 0;
};

/**
 * Find the haplotype unitigs of order k (k >= 2) in the graph, which must have
 * a bidirectional GBWT. With a gref, leave out the k-mers that lie along a gref
 * path, and cut the unitigs into the tigs between them. Each tig is reported
 * once, in the orientation with the most bases on forward gref strands. Uses
 * OpenMP threads.
 */
std::vector<Tig> haplotype_unitigs(const gbwtgraph::GBWTGraph& graph, size_t k,
                                   const GrefIndex* gref = nullptr,
                                   UnitigStatistics* statistics = nullptr);

/**
 * Find the gref paths (named with the prefix) that no other path in the GBZ
 * contains, as (name, length) pairs. Each gref path is cut from some other
 * path, so this only happens when the GBZ has lost that source path, for
 * example through haplotype sampling. The gref sequence then has k-mers that
 * no remaining haplotype has, and the unitigs, which count every path in the
 * GBWT as a haplotype, include them too.
 */
std::vector<std::pair<std::string, size_t>> unsupported_gref_paths(const gbwtgraph::GBZ& gbz,
                                                                   const std::string& name_prefix);

}

#endif
