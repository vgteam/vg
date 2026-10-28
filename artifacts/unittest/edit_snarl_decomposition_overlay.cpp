/// \file edit_snarl_decomposition_overlay.cpp
///
/// Unit tests for EditSnarlDecompositionOverlay, laid over both implementations of a decomposition.

#include <atomic>
#include <functional>
#include <vector>

#include "../../src/unittest/catch.hpp"
#include "support/decomposition_pair.hpp"
#include "../edit_snarl_decomposition_overlay.hpp"
#include "../../src/integrated_snarl_finder.hpp"
#include "support/nested_graph.hpp"
#include "../../src/vg.hpp"

namespace vg {
namespace unittest {

using namespace std;

namespace {

/// Get the children of `parent`, canonicalized.
vector<net_handle_t> children(const SnarlDecomposition& decomposition, const net_handle_t& parent) {
    vector<net_handle_t> out;
    decomposition.for_each_child(parent, [&](const net_handle_t& child) {
        out.push_back(decomposition.canonical(child));
    });
    return out;
}

/// Run `check` against the overlay laid over each implementation of the decomposition.
void on_both(NestedGraph& nested,
             const function<void(const SnarlDecomposition&, EditSnarlDecompositionOverlay&)>& check) {
    IntegratedSnarlFinder finder(nested.graph);
    DecompositionPair pair(nested.graph, finder);
    SECTION("over a SnarlManagerDecomposition") {
        EditSnarlDecompositionOverlay overlay(pair.adapter, nested.graph);
        check(pair.adapter, overlay);
    }
    SECTION("over a SnarlDistanceIndex") {
        EditSnarlDecompositionOverlay overlay(pair.distance_index, nested.graph);
        check(pair.distance_index, overlay);
    }
}

}

TEST_CASE("With no edits, the overlay shows the inner decomposition's tree",
          "[edit_snarl_decomposition_overlay]") {
    NestedGraph nested;
    on_both(nested, [&](const SnarlDecomposition& inner, EditSnarlDecompositionOverlay& overlay) {
        // Walk both trees together, comparing each element's kind and number of children.
        function<void(const net_handle_t&, const net_handle_t&)> compare =
            [&](const net_handle_t& inner_net, const net_handle_t& net) {
                REQUIRE(overlay.is_root(net) == inner.is_root(inner_net));
                REQUIRE(overlay.is_snarl(net) == inner.is_snarl(inner_net));
                REQUIRE(overlay.is_chain(net) == inner.is_chain(inner_net));
                REQUIRE(overlay.is_node(net) == inner.is_node(inner_net));
                REQUIRE_FALSE(overlay.is_edit(net));
                const vector<net_handle_t> inner_children = children(inner, inner_net);
                const vector<net_handle_t> overlay_children = children(overlay, net);
                REQUIRE(inner_children.size() == overlay_children.size());
                for (size_t i = 0; i < inner_children.size(); ++i) {
                    REQUIRE(overlay.canonical(overlay.get_parent(overlay_children[i]))
                            == overlay.canonical(net));
                    if (!inner.is_node(inner_children[i])) {
                        compare(inner_children[i], overlay_children[i]);
                    }
                }
            };
        compare(inner.get_root(), overlay.get_root());
    });
}

TEST_CASE("An edit inside one node is a leaf child of the node",
          "[edit_snarl_decomposition_overlay]") {
    NestedGraph nested;
    on_both(nested, [&](const SnarlDecomposition&, EditSnarlDecompositionOverlay& overlay) {
        const handle_t node = nested.handle(4);
        const net_handle_t edit = overlay.add_edit(node, 1, node, 3);
        const net_handle_t node_net = overlay.canonical(overlay.get_net(node, &nested.graph));

        REQUIRE(overlay.is_edit(edit));
        REQUIRE_FALSE(overlay.is_snarl(edit));
        REQUIRE_FALSE(overlay.is_chain(edit));
        REQUIRE_FALSE(overlay.is_node(edit));
        REQUIRE_FALSE(overlay.is_sentinel(edit));
        REQUIRE_FALSE(overlay.is_root(edit));

        REQUIRE(overlay.get_parent(edit) == node_net);
        REQUIRE(children(overlay, node_net) == vector<net_handle_t>{edit});
        REQUIRE(children(overlay, edit).empty());
        REQUIRE(overlay.follow_net_edges(edit, &nested.graph, false,
                                         [](const net_handle_t&) { return false; }));

        REQUIRE(overlay.get_site_start(edit) == node);
        REQUIRE(overlay.get_site_start_offset(edit) == 1);
        REQUIRE(overlay.get_site_end(edit) == node);
        REQUIRE(overlay.get_site_end_offset(edit) == 3);
        REQUIRE(overlay.starts_at(edit) == SnarlDecomposition::START);
        REQUIRE(overlay.ends_at(edit) == SnarlDecomposition::END);

        // Read the other way, the edit covers the same bases of the flipped node.
        const net_handle_t back = overlay.flip(edit);
        REQUIRE(back != edit);
        REQUIRE(overlay.canonical(back) == edit);
        REQUIRE(overlay.is_edit(back));
        REQUIRE(overlay.get_site_start(back) == nested.graph.flip(node));
        REQUIRE(overlay.get_site_start_offset(back) == 1);
        REQUIRE(overlay.get_site_end(back) == nested.graph.flip(node));
        REQUIRE(overlay.get_site_end_offset(back) == 3);
        REQUIRE(overlay.starts_at(back) == SnarlDecomposition::END);
        REQUIRE(overlay.ends_at(back) == SnarlDecomposition::START);
    });
}

TEST_CASE("An edit across nodes is a child of the lowest element holding them",
          "[edit_snarl_decomposition_overlay]") {
    NestedGraph nested;
    on_both(nested, [&](const SnarlDecomposition&, EditSnarlDecompositionOverlay& overlay) {
        const net_handle_t outer = find_snarl(overlay, nested.graph, 1, 8);
        const net_handle_t middle = find_snarl(overlay, nested.graph, 2, 7);

        SECTION("Two nodes of a top-level chain") {
            const net_handle_t edit = overlay.add_edit(nested.handle(1), 2, nested.handle(8), 1);
            REQUIRE(overlay.get_parent(edit) == overlay.canonical(overlay.get_parent(outer)));
        }
        SECTION("Two nodes of a child chain") {
            const net_handle_t edit = overlay.add_edit(nested.handle(2), 0, nested.handle(7), 1);
            REQUIRE(overlay.get_parent(edit) == overlay.canonical(overlay.get_parent(middle)));
            // The edit comes after the chain's own children.
            const vector<net_handle_t> listed = children(overlay, overlay.get_parent(edit));
            REQUIRE(listed.back() == edit);
            REQUIRE(listed.size() > 1);
        }
        SECTION("A node of a child chain and a node in a snarl beside it") {
            const net_handle_t edit = overlay.add_edit(nested.handle(3), 0, nested.handle(6), 1);
            REQUIRE(overlay.get_parent(edit) == overlay.canonical(middle));
        }
    });
}

TEST_CASE("A snarl's bounds are its boundary nodes' inward sides",
          "[edit_snarl_decomposition_overlay]") {
    NestedGraph nested;
    on_both(nested, [&](const SnarlDecomposition&, EditSnarlDecompositionOverlay& overlay) {
        const net_handle_t forward = find_snarl(overlay, nested.graph, 2, 7);
        REQUIRE(overlay.get_site_start(forward) == nested.handle(2));
        REQUIRE(overlay.get_site_start_offset(forward) == 1);
        REQUIRE(overlay.get_site_end(forward) == nested.handle(7));
        REQUIRE(overlay.get_site_end_offset(forward) == 0);

        const net_handle_t backward = overlay.flip(forward);
        REQUIRE(overlay.get_site_start(backward) == nested.handle(7, true));
        REQUIRE(overlay.get_site_start_offset(backward) == 1);
        REQUIRE(overlay.get_site_end(backward) == nested.handle(2, true));
        REQUIRE(overlay.get_site_end_offset(backward) == 0);
    });
}

TEST_CASE("Several threads can add edits and read them at once",
          "[edit_snarl_decomposition_overlay]") {
    NestedGraph nested;
    on_both(nested, [&](const SnarlDecomposition&, EditSnarlDecompositionOverlay& overlay) {
        const int threads = 8;
        const int per_thread = 500;
        vector<vector<net_handle_t>> added(threads);
        // Each thread adds edits inside node 4, CTGA, and reads each back as it goes.
        atomic<bool> all_read {true};
#pragma omp parallel for num_threads(threads)
        for (int t = 0; t < threads; ++t) {
            for (int i = 0; i < per_thread; ++i) {
                const size_t offset = (t + i) % 4;
                const net_handle_t edit =
                    overlay.add_edit(nested.handle(4), offset, nested.handle(4), offset + 1);
                added[t].push_back(edit);
                if (overlay.get_site_start_offset(edit) != offset) {
                    all_read = false;
                }
            }
        }
        REQUIRE(all_read);
        const net_handle_t node4 = overlay.canonical(overlay.get_net(nested.handle(4), &nested.graph));
        REQUIRE(children(overlay, node4).size() == threads * per_thread);
        for (int t = 0; t < threads; ++t) {
            for (int i = 0; i < per_thread; ++i) {
                REQUIRE(overlay.get_site_end_offset(added[t][i]) == (t + i) % 4 + 1);
            }
        }
    });
}

TEST_CASE("An edit with an offset past the end of its node is refused",
          "[edit_snarl_decomposition_overlay]") {
    NestedGraph nested;
    on_both(nested, [&](const SnarlDecomposition&, EditSnarlDecompositionOverlay& overlay) {
        REQUIRE_THROWS_AS(overlay.add_edit(nested.handle(4), 5, nested.handle(4), 5),
                          std::invalid_argument);
    });
}

}
}
