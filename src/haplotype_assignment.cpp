#include "haplotype_assignment.hpp"

#include <algorithm>
#include <functional>
#include <iostream>

#include <vg/vg.pb.h>

namespace vg {

using namespace std;

// Recursively walk the snarl decomposition starting from `net` and stamp every
// node reached with `current_snarl_id`. When we descend into a non-root snarl,
// allocate a fresh dense id and pass that id to all nodes in its child chains.
// Boundary nodes of a snarl live in the parent chain, so they keep the parent
// snarl's id (which is the desired "enclosing snarl" semantics).
static void walk_decomposition(const SnarlDistanceIndex& dist_index,
                               const net_handle_t& net,
                               uint32_t current_snarl_id,
                               uint32_t& next_snarl_id,
                               vector<uint32_t>& node_to_enclosing_snarl) {
    if (dist_index.is_node(net)) {
        nid_t id = dist_index.node_id(net);
        if (id >= 0 && static_cast<size_t>(id) < node_to_enclosing_snarl.size()) {
            node_to_enclosing_snarl[id] = current_snarl_id;
        }
        return;
    }

    uint32_t child_snarl_id = current_snarl_id;
    if (dist_index.is_snarl(net) && !dist_index.is_root(net)) {
        child_snarl_id = next_snarl_id++;
    }

    dist_index.for_each_child(net, [&](const net_handle_t& child) {
        walk_decomposition(dist_index, child, child_snarl_id,
                           next_snarl_id, node_to_enclosing_snarl);
        return true;
    });
}

HaplotypeAssigner::HaplotypeAssigner(const gbwtgraph::GBZ& gbz,
                                     const SnarlDistanceIndex& dist_index,
                                     bool verbose)
    : gbz_(&gbz), dist_index_(&dist_index)
{
    const auto& graph = gbz_->graph;
    size_t table_size = static_cast<size_t>(graph.max_node_id()) + 1;
    node_to_enclosing_snarl_.assign(table_size, 0u);

    uint32_t next_snarl_id = 1; // 0 is reserved for "no enclosing snarl".
    walk_decomposition(*dist_index_, dist_index_->get_root(),
                       0u, next_snarl_id, node_to_enclosing_snarl_);

    if (verbose) {
        size_t bytes = node_to_enclosing_snarl_.size() * sizeof(uint32_t);
        cerr << "[HaplotypeAssigner] " << (next_snarl_id - 1)
             << " distinct snarls; node_to_enclosing_snarl table: "
             << table_size << " entries (" << bytes << " bytes)" << endl;
    }
}

uint32_t HaplotypeAssigner::enclosing_snarl(nid_t node_id) const {
    if (node_id < 0) {
        return 0;
    }
    size_t idx = static_cast<size_t>(node_id);
    if (idx >= node_to_enclosing_snarl_.size()) {
        return 0;
    }
    return node_to_enclosing_snarl_[idx];
}

HaplotypeAssignmentResult HaplotypeAssigner::assign(const Alignment& aln,
                                                    bool extend_to_snarls) const {
    HaplotypeAssignmentResult result;

    // Step A: extract the alignment path as GBWT-encoded nodes.
    const auto& path = aln.path();
    if (path.mapping_size() == 0) {
        return result;
    }
    vector<gbwt::node_type> path_nodes;
    path_nodes.reserve(path.mapping_size());
    for (int i = 0; i < path.mapping_size(); ++i) {
        const auto& pos = path.mapping(i).position();
        path_nodes.push_back(gbwt::Node::encode(pos.node_id(), pos.is_reverse()));
    }

    // Step B: walk the GBWT along the alignment path, keeping a bidirectional
    // state so we can also extend outward afterwards.
    gbwt::BidirectionalState state = gbz_->index.bdFind(path_nodes.front());
    if (state.empty()) {
        return result;
    }
    for (size_t i = 1; i < path_nodes.size(); ++i) {
        state = gbz_->index.bdExtendForward(state, path_nodes[i]);
        if (state.empty()) {
            return result;
        }
    }

    // Step C: optionally extend outward to the nearest enclosing snarl
    // boundary. We stop extending as soon as:
    //   - there are no GBWT-consistent neighbors (walked off the graph),
    //   - there is more than one neighbor with non-empty state (branching),
    //   - the single neighbor's enclosing-snarl id differs from the starting
    //     node's snarl id (boundary found: record that neighbor's node id),
    //   - max_extension_nodes steps have been taken on that side.
    // Extensions narrow the BidirectionalState to haplotypes consistent with
    // the extended region, which is the point of the option.
    if (extend_to_snarls) {
        const uint32_t start_snarl_id =
            enclosing_snarl(gbwt::Node::id(path_nodes.front()));

        // Left extension.
        for (size_t step = 0; step < max_extension_nodes; ++step) {
            vector<gbwt::BidirectionalState> nexts;
            gbz_->graph.follow_paths(state, /*backward=*/true,
                [&](const gbwt::BidirectionalState& next) -> bool {
                    if (!next.empty()) {
                        nexts.push_back(next);
                    }
                    return true;
                });
            if (nexts.empty() || nexts.size() > 1) {
                break;
            }
            // The newly added backward node is stored reversed; flip to get
            // its forward-orientation GBWT encoding.
            gbwt::node_type pred = gbwt::Node::reverse(nexts[0].backward.node);
            nid_t pred_id = gbwt::Node::id(pred);
            if (enclosing_snarl(pred_id) != start_snarl_id) {
                result.left_boundary_node = pred_id;
                break;
            }
            state = nexts[0];
        }

        // Right extension.
        for (size_t step = 0; step < max_extension_nodes; ++step) {
            vector<gbwt::BidirectionalState> nexts;
            gbz_->graph.follow_paths(state, /*backward=*/false,
                [&](const gbwt::BidirectionalState& next) -> bool {
                    if (!next.empty()) {
                        nexts.push_back(next);
                    }
                    return true;
                });
            if (nexts.empty() || nexts.size() > 1) {
                break;
            }
            gbwt::node_type succ = nexts[0].forward.node;
            nid_t succ_id = gbwt::Node::id(succ);
            if (enclosing_snarl(succ_id) != start_snarl_id) {
                result.right_boundary_node = succ_id;
                break;
            }
            state = nexts[0];
        }
    }

    // Step D: enumerate compatible haplotype sequence ids from the forward
    // SearchState of the (possibly extended) bidirectional state.
    result.candidate_count = state.size();
    size_t n = std::min(state.size(), max_haplotypes_to_report);
    result.seq_ids.reserve(n);
    const gbwt::SearchState& fwd = state.forward;
    for (size_t i = 0; i < n; ++i) {
        result.seq_ids.push_back(gbz_->index.locate(fwd.node, fwd.range.first + i));
    }
    result.truncated = (n < state.size());

    return result;
}

} // namespace vg
