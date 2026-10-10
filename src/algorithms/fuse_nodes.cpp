/**
 * \file fuse_nodes.cpp
 *
 * Defines an algorithm to fuse one pair of adjacent nodes into one node.
 */

#include "fuse_nodes.hpp"

#include <unordered_set>
#include <vector>

namespace vg {
namespace algorithms {

using namespace std;

/// Check whether two nodes can be fused.
bool can_fuse_nodes(const HandleGraph& graph, handle_t left, handle_t right) {
    if (graph.get_id(left) == graph.get_id(right)) {
        // A node cannot be concatenated onto itself.
        return false;
    }
    bool clean = true;
    // left must have nothing attached on its right side.
    graph.follow_edges(left, false, [&](const handle_t&) {
        clean = false;
        return false;
    });
    if (!clean) {
        return false;
    }
    // right must have nothing attached on its left side.
    graph.follow_edges(right, true, [&](const handle_t&) {
        clean = false;
        return false;
    });
    return clean;
}

/// Fuse two nodes into one node.
handle_t fuse_nodes(handlegraph::MutablePathDeletableHandleGraph* graph,
                    handle_t left, handle_t right) {

    // Create the fused node.
    handle_t fused = graph->create_handle(graph->get_sequence(left)
                                         + graph->get_sequence(right));

    // Collect all neighbours of left and right.
    unordered_set<handle_t> predecessors;
    unordered_set<handle_t> successors;
    graph->follow_edges(left, true, [&](const handle_t& p) {
        predecessors.insert(p);
    });
    graph->follow_edges(right, false, [&](const handle_t& n) {
        successors.insert(n);
    });

    // Handle a neighbour that is left or right itself, so the edges built below
    // land on fused.
    auto translate = [&](const handle_t& h) {
        if (h == left || h == right) {
            return fused;
        }
        if (h == graph->flip(left) || h == graph->flip(right)) {
            return graph->flip(fused);
        }
        return h;
    };
    // Build the edges to the neighbours.
    for (const handle_t& p : predecessors) {
        handle_t from = translate(p);
        if (!graph->has_edge(from, fused)) {
            graph->create_edge(from, fused);
        }
    }
    for (const handle_t& n : successors) {
        handle_t to = translate(n);
        if (!graph->has_edge(fused, to)) {
            graph->create_edge(fused, to);
        }
    }

    // Collect every step that visits left or right.
    vector<pair<step_handle_t, bool>> to_rewrite;
    for (handle_t original : {left, right}) {
        graph->for_each_step_on_handle(original, [&](const step_handle_t& s) {
            bool flipped = graph->get_is_reverse(original)
                != graph->get_is_reverse(graph->get_handle_of_step(s));
            to_rewrite.emplace_back(s, flipped);
        });
    }
    // Rewrite each collected step to point at fused.
    for (const auto& step_and_flip : to_rewrite) {
        graph->rewrite_segment(step_and_flip.first,
                               graph->get_next_step(step_and_flip.first),
                               {step_and_flip.second ? graph->flip(fused) : fused});
    }

    // Tear down the old nodes and the edges on them.
    for (handle_t original : {left, right}) {
        unordered_set<edge_t> to_remove;
        graph->follow_edges(original, false, [&](const handle_t& h) {
            to_remove.insert(graph->edge_handle(original, h));
        });
        graph->follow_edges(original, true, [&](const handle_t& h) {
            to_remove.insert(graph->edge_handle(h, original));
        });
        for (const edge_t& e : to_remove) {
            graph->destroy_edge(e);
        }
        graph->destroy_handle(original);
    }

    return fused;
}

}
}
