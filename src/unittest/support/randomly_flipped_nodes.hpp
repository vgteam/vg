#ifndef VG_UNITTEST_RANDOMLY_FLIPPED_NODES_HPP_INCLUDED
#define VG_UNITTEST_RANDOMLY_FLIPPED_NODES_HPP_INCLUDED

/**
 * \file randomly_flipped_nodes.hpp
 * Utility for creating a copy of a HandleGraph with a random subset of nodes
 * flipped in orientation.
 */

#include <random>
#include <bdsg/hash_graph.hpp>
#include "handle.hpp"

namespace vg {
namespace unittest {

/**
 * Return a copy of the given graph with approximately p_flip fraction of its
 * nodes reversed in their local forward orientation. When a node is flipped,
 * its sequence is reverse-complemented and all edges that connected to its
 * forward orientation now connect to its reverse orientation, and vice versa.
 *
 * The returned graph preserves node IDs.
 */
template<typename URNG>
bdsg::HashGraph randomly_flipped_nodes(const HandleGraph& source, double p_flip, URNG& generator) {
    bdsg::HashGraph result;

    std::uniform_real_distribution<double> dist(0.0, 1.0);

    // Track which nodes get flipped
    std::unordered_set<nid_t> flipped;

    // Copy all nodes, flipping some
    source.for_each_handle([&](const handle_t& handle) {
        nid_t id = source.get_id(handle);
        if (dist(generator) < p_flip) {
            // Flip this node: store its reverse complement sequence as forward
            result.create_handle(source.get_sequence(source.flip(handle)), id);
            flipped.insert(id);
        } else {
            // Keep this node as-is
            result.create_handle(source.get_sequence(handle), id);
        }
    });

    // Copy all edges, adjusting for flipped nodes.
    // An edge (left, right) means: leave left in its orientation, enter right
    // in its orientation. If we flipped a node, we need to toggle the
    // orientation on that side of the edge.
    source.for_each_edge([&](const edge_t& edge) {
        handle_t left = edge.first;
        handle_t right = edge.second;

        nid_t left_id = source.get_id(left);
        bool left_is_reverse = source.get_is_reverse(left);

        nid_t right_id = source.get_id(right);
        bool right_is_reverse = source.get_is_reverse(right);

        // If we flipped a node, toggle the orientation for that side
        if (flipped.count(left_id)) {
            left_is_reverse = !left_is_reverse;
        }
        if (flipped.count(right_id)) {
            right_is_reverse = !right_is_reverse;
        }

        result.create_edge(
            result.get_handle(left_id, left_is_reverse),
            result.get_handle(right_id, right_is_reverse)
        );

        return true;
    });

    return result;
}

} // namespace unittest
} // namespace vg

#endif
