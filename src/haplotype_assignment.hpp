#ifndef VG_HAPLOTYPE_ASSIGNMENT_HPP_INCLUDED
#define VG_HAPLOTYPE_ASSIGNMENT_HPP_INCLUDED

/** \file
 * Find the haplotypes in a GBZ that carry an alignment's path.
 */

#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <gbwt/gbwt.h>
#include <gbwtgraph/gbz.h>

#include "snarl_distance_index.hpp"

namespace vg {

// Declared here so this header need not include vg.pb.h.
class Alignment;

/**
 * The haplotypes that carry an alignment's path.
 *
 * The path is split into segments, each one as long as some haplotype carries
 * it unbroken. A path that one haplotype carries end to end is one segment.
 * The longest segment is the best one, and the top-level fields describe it.
 */
struct HaplotypeAssignmentResult {
    /**
     * A run of the path's mappings that at least one haplotype carries.
     *
     * The run ends where none of the haplotypes that carry it go on to the
     * next node of the path.
     */
    struct Segment {
        size_t start_mapping   = 0;   ///< Index of the first mapping in the run.
        size_t end_mapping     = 0;   ///< One past the index of the last mapping.
        size_t covered_bp      = 0;   ///< Total length of the run's nodes.
        /// Number of times GBWT sequences carry the run. A sequence that
        /// carries it twice counts twice.
        size_t candidate_count = 0;
        /// The GBWT sequence IDs that carry the run, at most
        /// max_haplotypes_to_report of them.
        std::vector<gbwt::size_type> seq_ids;
        /// Distinct names of the haplotypes in seq_ids. The caller fills them in.
        std::vector<std::string> haplotype_names;
        /// One haplotype from seq_ids to report alone, preferably a reference.
        /// The caller fills it in.
        std::string representative_name;
    };

    /// The best segment's seq_ids.
    std::vector<gbwt::size_type> seq_ids;
    /// The best segment's haplotype names, as sorted PanSN names. The caller
    /// fills them in. Without GBWT metadata a sequence ID stands in for a name.
    std::vector<std::string> haplotype_names;
    /// The best segment's candidate_count. It is more than seq_ids.size() when
    /// the cap cut seq_ids short.
    size_t candidate_count = 0;
    /// The boundary node before the path, in the path's orientation, of the
    /// smallest snarl that contains the whole path. 0 unless extension was
    /// asked for, the whole path is one segment, some snarl contains it, and
    /// the walk got that far.
    nid_t left_boundary_node = 0;
    /// The boundary node after the path, found the same way as
    /// left_boundary_node.
    nid_t right_boundary_node = 0;
    /// The best segment's representative_name. The caller fills it in.
    std::string representative_name;
    /// True if the representative is from a reference sample or the generic
    /// sample. The caller fills it in.
    bool representative_is_reference = false;

    /// The segments, in path order.
    std::vector<Segment> segments;
    /// True if one segment spans the whole path.
    bool fully_covered = false;
    /// Total length of the path's nodes.
    size_t path_length_bp = 0;
    /// Total length of the best segment's nodes.
    size_t best_covered_bp = 0;
    /// Index of the best segment in segments, or SIZE_MAX if there are none.
    size_t best_segment = std::numeric_limits<size_t>::max();
};

/**
 * Finds the haplotypes in a GBZ that carry an alignment's path.
 *
 * It can also walk out from the path while those haplotypes agree, to find
 * where they leave the smallest snarl that contains the path. Thread-safe.
 */
class HaplotypeAssigner {
public:
    /// The GBZ and the distance index must outlive the assigner.
    HaplotypeAssigner(const gbwtgraph::GBZ& gbz,
                      const SnarlDistanceIndex& dist_index);

    /**
     * Find the haplotypes that carry the alignment's path.
     *
     * If extend_to_snarls is true, one segment covers the whole path, and a
     * snarl contains the whole path, also walk outward from each end while the
     * haplotypes all take the same node, out to the smallest such snarl's
     * boundary nodes. This sets the boundary nodes and keeps only the
     * haplotypes that make the whole walk.
     */
    HaplotypeAssignmentResult assign(const Alignment& aln, bool extend_to_snarls) const;

    /// The most sequence IDs to list per segment.
    size_t max_haplotypes_to_report = 1000;
    /// The most nodes to walk out from each end of the path when extending.
    size_t max_extension_nodes = 1000;

private:
    const gbwtgraph::GBZ* gbz_;
    const SnarlDistanceIndex* dist_index_;

    /// Return the snarls that contain the node, innermost first, as canonical
    /// net handles. A snarl's boundary nodes are in the snarl around it, not in
    /// the snarl itself. A node in no snarl gets an empty list.
    std::vector<net_handle_t> enclosing_snarls(nid_t node_id) const;

    /// Return true if the snarl, as a canonical net handle, contains the node.
    bool in_snarl(nid_t node_id, const net_handle_t& snarl) const;
};

/// Return the name of a GBWT path's haplotype, as sample#phase#contig. A
/// sample or contig that the metadata does not name appears as its number.
///
/// This is not the graph's name for the path, which also has the fragment.
/// The fragments of a haplotype share this name.
std::string haplotype_name(const gbwt::Metadata& metadata, gbwt::size_type path_id);

/// Return the distinct names, sorted, of the haplotypes with these GBWT
/// sequence IDs.
///
/// A haplotype can have many sequence IDs: one per orientation, per visit
/// and per fragment of its path. A sequence ID that GBWT metadata does not
/// name stands for itself, as a number.
std::vector<std::string> haplotype_names(const gbwtgraph::GBZ& gbz,
                                         const std::vector<gbwt::size_type>& seq_ids);

/// Choose one haplotype to stand for all the haplotypes with these GBWT
/// sequence IDs.
///
/// Prefer a reference sample or the generic path sample, and then the lowest
/// sequence ID, so that the choice is the same on every run. Return the
/// chosen haplotype's name, and whether it is a reference. Return an empty
/// name if there are no sequence IDs or no GBWT metadata.
std::pair<std::string, bool> choose_representative(const gbwtgraph::GBZ& gbz,
                                                   const std::vector<gbwt::size_type>& seq_ids);

} // namespace vg

#endif // VG_HAPLOTYPE_ASSIGNMENT_HPP_INCLUDED
