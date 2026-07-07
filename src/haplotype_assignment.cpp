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

    // Step A: extract the alignment path as GBWT-encoded nodes, plus each node's
    // length in bases (for coverage accounting).
    const auto& path = aln.path();
    if (path.mapping_size() == 0) {
        return result;
    }
    vector<gbwt::node_type> path_nodes;
    vector<size_t> node_bp;
    path_nodes.reserve(path.mapping_size());
    node_bp.reserve(path.mapping_size());
    for (int i = 0; i < path.mapping_size(); ++i) {
        const auto& pos = path.mapping(i).position();
        path_nodes.push_back(gbwt::Node::encode(pos.node_id(), pos.is_reverse()));
        node_bp.push_back(gbz_->graph.get_length(
            gbz_->graph.get_handle(pos.node_id(), pos.is_reverse())));
    }
    const size_t n = path_nodes.size();
    for (size_t bp : node_bp) {
        result.path_length_bp += bp;
    }

    // Step B: greedily tile the read path into maximal haplotype-consistent
    // segments. From each start we walk the GBWT with a bidirectional state as
    // far forward as any haplotype allows; when no haplotype can cross the next
    // edge (a recombination point — common in tandem repeats, inversions, and
    // dense-variant regions) we close the segment and restart the search at that
    // node. A single haplotype that carries the whole read yields exactly one
    // segment; a mosaic read yields several. Because covering a walk [i..j]
    // requires covering its prefix [i..j-1], the run we stop at is genuinely
    // maximal: no haplotype carries one more node.
    size_t i = 0;
    while (i < n) {
        gbwt::BidirectionalState state = gbz_->index.bdFind(path_nodes[i]);
        if (state.empty()) {
            // Node not represented in the GBWT in this orientation; cannot
            // anchor here. Advance so we always make progress.
            ++i;
            continue;
        }
        size_t j = i + 1;
        while (j < n) {
            gbwt::BidirectionalState next =
                gbz_->index.bdExtendForward(state, path_nodes[j]);
            if (next.empty()) {
                break;
            }
            state = next;
            ++j;
        }

        // Step C: only for the clean case — when the very first segment already
        // spans the whole read — optionally extend outward to the nearest
        // enclosing snarl boundaries, narrowing the haplotype set to those
        // consistent with the flanking region (unchanged from the original
        // behavior). Mosaic segments are not extended: their ends are
        // recombination points where the consistent set changes by construction.
        const bool whole_path = (i == 0 && j == n);
        if (extend_to_snarls && whole_path) {
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

        // Enumerate the haplotypes carrying this whole segment.
        HaplotypeAssignmentResult::Segment seg;
        seg.start_mapping = i;
        seg.end_mapping   = j;
        for (size_t k = i; k < j; ++k) {
            seg.covered_bp += node_bp[k];
        }
        seg.candidate_count = state.size();
        const size_t cap = std::min(state.size(), max_haplotypes_to_report);
        const gbwt::SearchState& fwd = state.forward;
        seg.seq_ids.reserve(cap);
        for (size_t t = 0; t < cap; ++t) {
            seg.seq_ids.push_back(gbz_->index.locate(fwd.node, fwd.range.first + t));
        }
        result.segments.push_back(std::move(seg));

        i = j;  // restart at the recombination point (or end).
    }

    if (result.segments.empty()) {
        return result;  // Nothing in the path could be anchored.
    }

    // Step D: pick the best (longest) segment and promote it to the top level so
    // the existing naming / representative / surjection machinery operates on the
    // haplotype(s) covering the most of the read path.
    size_t best = 0;
    for (size_t s = 1; s < result.segments.size(); ++s) {
        if (result.segments[s].covered_bp > result.segments[best].covered_bp) {
            best = s;
        }
    }
    result.best_segment    = best;
    result.best_covered_bp = result.segments[best].covered_bp;
    result.fully_covered   = (result.segments.size() == 1 &&
                              result.segments.front().start_mapping == 0 &&
                              result.segments.front().end_mapping == n);

    result.seq_ids         = result.segments[best].seq_ids;
    result.candidate_count = result.segments[best].candidate_count;
    result.truncated       =
        (result.segments[best].seq_ids.size() < result.segments[best].candidate_count);

    return result;
}

} // namespace vg
