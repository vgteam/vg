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
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include <gbwt/gbwt.h>
#include <gbwtgraph/gbz.h>

#include "snarl_distance_index.hpp"

namespace vg {

// Forward declaration to avoid including vg.pb.h from this header.
class Alignment;

/// Result of assigning haplotypes to a single alignment.
struct HaplotypeAssignmentResult {
    /// One maximal run of the read's alignment path that is carried, as a single
    /// walk, by at least one haplotype. The read path is tiled greedily left to
    /// right into such segments: a boundary between two segments is a
    /// recombination point (an edge no single haplotype takes). When the whole
    /// read is haplotype-consistent there is exactly one segment spanning it.
    struct Segment {
        size_t start_mapping   = 0;   ///< first path-mapping index (inclusive)
        size_t end_mapping     = 0;   ///< one past the last path-mapping index
        size_t covered_bp      = 0;   ///< sum of node lengths over the segment
        size_t candidate_count = 0;   ///< haplotypes carrying the whole segment (pre-cap)
        /// Carrying-haplotype GBWT sequence ids for this segment (capped at
        /// max_haplotypes_to_report).
        std::vector<gbwt::size_type> seq_ids;
        /// Distinct haplotype/path names for this segment (filled by the caller).
        std::vector<std::string> haplotype_names;
        /// Representative haplotype name for this segment (reference-first, else
        /// lowest seq id; filled by the caller). Empty if no candidates.
        std::string representative_name;
    };

    /// Compatible haplotype sequence IDs (GBWT sequence IDs), possibly capped.
    /// These are the haplotypes covering the BEST (longest) segment — for a
    /// fully-covered read that is the whole path (identical to the old
    /// behavior); for a mosaic read it is the longest single consistent block.
    std::vector<gbwt::size_type> seq_ids;
    /// Distinct haplotype/path names (e.g. "HG00290#1#chr10") for `seq_ids`,
    /// deduplicated and sorted. Collapses the raw sequence ids (which double-
    /// count orientations, repeats, and fragments) into the actual set of
    /// haplotypes — far smaller and human-readable. Populated by the caller
    /// from GBWT metadata; empty if metadata is unavailable.
    std::vector<std::string> haplotype_names;
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
    /// A single representative haplotype name chosen to break the (often large)
    /// tie when many haplotypes carry the read's path identically. Policy:
    /// prefer a reference-sample haplotype; among the eligible class take the
    /// lowest GBWT sequence id (deterministic). Empty if no candidates.
    /// Populated by the caller. See compute in giraffe_engine.cpp.
    std::string representative_name;
    /// True if the representative was chosen because it is a reference-sample
    /// (or generic backbone) haplotype; false if it is just the lowest-id
    /// fallback (no reference among the candidates).
    bool representative_is_reference = false;

    /// Greedy left-to-right decomposition of the read path into maximal
    /// haplotype-consistent segments (in path order). Size 1 ⇔ fully_covered.
    /// A read that spans a recombination, tandem repeat, or inversion in a way
    /// no single haplotype realizes will have >1 segment — the mosaic.
    std::vector<Segment> segments;
    /// True iff some haplotype carries the ENTIRE read path as a single walk
    /// (segments.size() == 1 spanning the whole path). This is the old
    /// "clean" case; false means the read is a mosaic and no single haplotype
    /// is consistent end to end.
    bool fully_covered = false;
    /// Total bp of the read's node path (sum of node lengths).
    size_t path_length_bp = 0;
    /// bp covered by the best (longest) single segment — how much of the read
    /// path the reported haplotype(s) cover as one consistent block.
    size_t best_covered_bp = 0;
    /// Index into `segments` of the best (longest) segment, whose haplotypes are
    /// promoted to the top-level seq_ids/representative. SIZE_MAX if no segment.
    size_t best_segment = std::numeric_limits<size_t>::max();
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
