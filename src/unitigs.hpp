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
 * Only the GBWTGraph interface is used here.
 */

#include <gbwtgraph/gbz.h>

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace vg {

/**
 * Where the gref paths lie in the graph: for each node, the gref path that
 * visits it and the path coordinate of its first base.
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

    /// Where the node lies on the gref, or a Position with path == NO_PATH.
    const Position& find(handlegraph::nid_t id) const;

    size_t path_count() const { return names.size(); }
    /// Path name without any subrange.
    const std::string& path_name(uint32_t path) const { return names[path]; }
    /// Amount to add to a path coordinate to get a coordinate on the full path.
    uint64_t path_shift(uint32_t path) const { return shifts[path]; }
    /// Base at the given path coordinate, as read by the path.
    char base(uint32_t path, uint64_t coord) const;

    /// Number of bases on gref paths.
    uint64_t total_length() const;
    /// Number of node visits that were ignored because the node was already on the gref.
    size_t repeated_visits() const { return repeats; }

private:
    const gbwtgraph::GBWTGraph& graph;
    handlegraph::nid_t min_id = 0;
    std::vector<Position> positions;
    std::vector<std::string> names;
    std::vector<uint64_t> shifts;
    /// For each path, the coordinate where each step starts, and the handle as read.
    std::vector<std::vector<uint64_t>> step_starts;
    std::vector<std::vector<handlegraph::handle_t>> step_handles;
    size_t repeats = 0;
};

/// A walk in the graph whose sequence is part of the cover.
struct Tig {
    std::vector<handlegraph::handle_t> walk;
    /// Offset of the first base in the first node of the walk.
    size_t offset = 0;
    std::string sequence;
    /// Number of haplotype occurrences of each of its k-mers.
    size_t support = 0;
    /// Placement on the gref, as comma-separated PATH:START-END:STRAND:CIGAR
    /// segments (1-based, inclusive, CIGAR in path order). Empty without a gref.
    std::string gref;

    /// Sort order: gref path and coordinate of the first segment, if any.
    uint32_t sort_path = GrefIndex::NO_PATH;
    uint64_t sort_coord = 0;
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

}

#endif
