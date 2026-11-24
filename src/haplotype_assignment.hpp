#ifndef VG_HAPLOTYPE_ASSIGNMENT_HPP_INCLUDED
#define VG_HAPLOTYPE_ASSIGNMENT_HPP_INCLUDED

/** \file
 * Find the haplotypes in a GBZ that carry an alignment's path.
 */

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
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
 * A run of an alignment's path that at least one haplotype follows without a
 * break.
 *
 * The run ends where none of the haplotypes that follow it go on to the next
 * node of the path.
 */
struct HaplotypeSegment {
    /// Index of the first mapping in the run.
    size_t first_mapping = 0;
    /// One past the index of the last mapping.
    size_t past_last_mapping = 0;
    /// Total length of the run's nodes.
    size_t length_bp = 0;
    /// Number of times GBWT sequences follow the run. A sequence that
    /// follows it twice counts twice.
    size_t visit_count = 0;
    /// The GBWT sequence IDs that follow the run, at most
    /// max_haplotypes_to_report of them.
    std::vector<gbwt::size_type> seq_ids;
};

/// Return the segment that covers the most bases, the first one if several
/// do, or nullptr if there are no segments. Its haplotypes carry the most of
/// the path.
const HaplotypeSegment* longest_segment(const std::vector<HaplotypeSegment>& segments);

/**
 * The nodes where the haplotypes that follow a whole path leave the smallest
 * snarl that contains it. Left is before the path and right is after it, in
 * the path's orientation. A side is unset if its node was not found.
 */
struct SnarlBoundaries {
    std::optional<nid_t> left;
    std::optional<nid_t> right;
};

/**
 * Answers two questions about an alignment's path, from a GBZ: which
 * haplotypes follow it, and where those haplotypes leave the snarl around it.
 * Thread-safe.
 */
class HaplotypeAssigner {
public:
    /// The GBZ and the distance index must outlive the assigner.
    HaplotypeAssigner(const gbwtgraph::GBZ& gbz,
                      const SnarlDistanceIndex& dist_index);

    /**
     * Split the alignment's path into segments, in path order, each one as
     * long as some haplotype follows it. A path that one haplotype follows
     * end to end is one segment. A node that no haplotype visits is in no
     * segment.
     */
    std::vector<HaplotypeSegment> segments(const Alignment& aln) const;

    /**
     * Find the boundary nodes of the smallest snarl that contains the
     * alignment's whole path.
     *
     * Walk out from each end of the path while the haplotypes that follow the
     * whole path all take the same node, through child snarls, to the first
     * node outside that snarl. The walk to the right uses only the haplotypes
     * that were left at the end of the walk to the left.
     *
     * A side is unset if no haplotype follows the whole path, no snarl
     * contains it, the haplotypes part before they leave the snarl, or the
     * walk takes more than max_extension_nodes steps.
     */
    SnarlBoundaries snarl_boundaries(const Alignment& aln) const;

    /// The most sequence IDs to list per segment.
    size_t max_haplotypes_to_report = 1000;
    /// The most nodes to walk out from each end of the path.
    size_t max_extension_nodes = 1000;

private:
    const gbwtgraph::GBZ* gbz_;
    const SnarlDistanceIndex* dist_index_;

    /// Follow the nodes from nodes[start] with a GBWT search, as far as some
    /// haplotype goes. Return the search state, and set end to one past the
    /// last node reached. If no haplotype visits nodes[start], the state is
    /// empty and end is start.
    gbwt::BidirectionalState follow(const std::vector<gbwt::node_type>& nodes,
                                    size_t start, size_t& end) const;

    /// Return the snarls that contain the node, innermost first, as canonical
    /// net handles. A snarl's boundary nodes are in the snarl around it, not in
    /// the snarl itself. A node in no snarl gets an empty list.
    std::vector<net_handle_t> enclosing_snarls(nid_t node_id) const;

    /// Return true if the snarl, as a canonical net handle, contains the node.
    bool in_snarl(nid_t node_id, const net_handle_t& snarl) const;
};

/// Return the name of a GBWT path's haplotype, as sample#phase#contig. A
/// sample or contig that the metadata does not name appears as its number.
/// The fragments of a haplotype share this name, so it is shorter than the
/// graph's name for the path, which also has the fragment.
///
/// A generic path gets its contig alone, which is the graph's name for it.
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
