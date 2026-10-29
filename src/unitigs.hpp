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
 * tool can scan in place of the haplotypes. feature_lift.hpp translates what such
 * a tool finds back to the graph.
 *
 * Only the GBWTGraph interface is used here.
 */

#include "gref.hpp"
#include "walk_interval.hpp"

#include <gbwtgraph/gbz.h>

#include <vector>

namespace vg {

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

}

#endif
