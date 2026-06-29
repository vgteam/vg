#ifndef VG_HAPLOTYPE_ASSIGNMENT_HPP_INCLUDED
#define VG_HAPLOTYPE_ASSIGNMENT_HPP_INCLUDED

/** \file
 * Post-processing utility that, for a final Giraffe Alignment, computes the
 * set of haplotypes that traverse the alignment's graph path, optionally
 * extended outward to the nearest enclosing snarl boundaries.
 *
 * Uses the already-loaded GBZ (GBWT + GBWTGraph) and SnarlDistanceIndex that
 * the engine owns. No new index files are required.
 *
 * Thread-safe after construction: assign() is const and does not mutate any
 * member state. The precomputed node_to_enclosing_snarl table is read-only.
 */

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include <gbwt/gbwt.h>
#include <gbwtgraph/gbz.h>

#include "snarl_distance_index.hpp"

namespace vg {

// Forward declaration to avoid including vg.pb.h from this header.
class Alignment;

/// Result of assigning haplotypes to a single alignment.
struct HaplotypeAssignmentResult {
    /// Compatible haplotype sequence IDs (GBWT sequence IDs), possibly capped.
    std::vector<gbwt::size_type> seq_ids;
    /// Number of haplotypes consistent with the final SearchState, before any
    /// capping by max_haplotypes_to_report.
    size_t candidate_count = 0;
    /// Left snarl-boundary node (in forward orientation) if the alignment was
    /// extended to a boundary. 0 if no extension was requested, or if
    /// extension walked off the graph, hit a branch, or was not performed.
    nid_t left_boundary_node = 0;
    /// Right snarl-boundary node, analogous to left_boundary_node.
    nid_t right_boundary_node = 0;
    /// True if seq_ids was truncated to max_haplotypes_to_report.
    bool truncated = false;
};

/**
 * Builds a dense "node id -> enclosing snarl id" table once from a
 * SnarlDistanceIndex, and uses it together with the GBZ to answer haplotype
 * membership queries for alignment paths.
 *
 * Usage:
 *   HaplotypeAssigner assigner(gbz, dist_index);
 *   auto result = assigner.assign(aln, /extend_to_snarls=/ true);
 */
class HaplotypeAssigner {
public:
    /// Build the enclosing-snarl table immediately. The GBZ and
    /// SnarlDistanceIndex must outlive this object. If verbose is true, logs
    /// a one-line summary (snarl count and table footprint) to std::cerr.
    HaplotypeAssigner(const gbwtgraph::GBZ& gbz,
                      const SnarlDistanceIndex& dist_index,
                      bool verbose = false);

    /// Compute haplotypes compatible with the alignment's path. Const and
    /// thread-safe with respect to concurrent callers.
    HaplotypeAssignmentResult assign(const Alignment& aln, bool extend_to_snarls) const;

    /// Cap on the number of haplotype IDs returned per alignment.
    size_t max_haplotypes_to_report = 1000;
    /// Cap on the number of graph nodes traversed in each direction when
    /// extend_to_snarls is true.
    size_t max_extension_nodes = 1000;

private:
    const gbwtgraph::GBZ* gbz_;
    const SnarlDistanceIndex* dist_index_;

    /// node_to_enclosing_snarl_[id] is the dense integer id of the snarl that
    /// immediately encloses the node with the given node id. 0 is a sentinel
    /// meaning "no enclosing snarl" (nodes directly under the root).
    std::vector<uint32_t> node_to_enclosing_snarl_;

    /// Return the enclosing-snarl id for the given node id, or 0 if the node
    /// is out of range.
    uint32_t enclosing_snarl(nid_t node_id) const;
};

} // namespace vg

#endif // VG_HAPLOTYPE_ASSIGNMENT_HPP_INCLUDED
