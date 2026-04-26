#pragma once
#include <cstddef>
#include <functional>
#include <span>
#include <utility>
#include <bdsg/snarl_distance_index.hpp>
#include <handlegraph/handle_graph.hpp>

namespace vg {

// Read-only view of a snarl's child net-graph.
//
// Encapsulates the "get outgoing handle → follow_edges → map to snarl child"
// kernel shared by topo_sort_children and populate_distance_matrix_row.
// Both callers need the same three-step operation:
//   1. compute the handle pointing out of `child` in direction `go_left`
//   2. follow graph edges from that handle
//   3. for each landing node, resolve its snarl-level ancestor + direction
//
// Dijkstra priority-queue management and topo-sort BFS logic stay with the
// respective callers.
class SnarlChildGraph {
public:
    using temp_record_ref_t  = SnarlDistanceIndex::temp_record_ref_t;
    using TempIndex          = SnarlDistanceIndex::TemporaryDistanceIndex;

    // `children` may be empty if the caller doesn't need the children() accessor.
    // All other methods only require temp_index, snarl_index, and graph.
    SnarlChildGraph(TempIndex& temp_index,
                    temp_record_ref_t snarl_index,
                    std::span<const temp_record_ref_t> children,
                    const handlegraph::HandleGraph* graph);

    std::span<const temp_record_ref_t> children() const noexcept;

    // For each graph edge leaving `child` in direction `go_left`, invoke the
    // callback with:
    //   neighbor        — snarl-level ancestor of the landing node
    //   neighbor_rev    — true iff `neighbor` is entered from its right side
    //                     (chain traversed right-to-left, or node reversed)
    //   edge_distance   — traversal length of `neighbor`:
    //                     TEMP_NODE → sequence length of the landing node;
    //                     TEMP_CHAIN → chain.min_length (∞ if disconnected)
    //   arriving_node_id — graph node id of the immediate landing handle
    //                      (before ancestor resolution; needed by callers that
    //                       must preserve the original follow_edges semantics
    //                       for is_simple detection)
    //
    // Boundary nodes ARE included in the callback; callers filter as needed.
    void for_each_outgoing(
        temp_record_ref_t child,
        bool go_left,
        const std::function<void(temp_record_ref_t neighbor,
                                 bool neighbor_rev,
                                 size_t edge_distance,
                                 handlegraph::nid_t arriving_node_id)>&) const;

    // Returns {boundary_node_ref, start_node_rev} for start=true,
    //         {boundary_node_ref, end_node_rev}   for start=false.
    std::pair<temp_record_ref_t, bool> boundary(bool start) const;

private:
    TempIndex&                           temp_index_;
    temp_record_ref_t                    snarl_index_;
    std::span<const temp_record_ref_t>   children_;
    const handlegraph::HandleGraph*      graph_;
};

} // namespace vg
