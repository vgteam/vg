/// \file edit_snarl_decomposition_overlay.cpp
///
/// Unit tests for EditSnarlDecompositionOverlay, laid over both implementations of a decomposition.

#include <functional>
#include <vector>

#include "catch.hpp"
#include "decomposition_pair.hpp"
#include "edit_snarl_decomposition_overlay.hpp"
#include "integrated_snarl_finder.hpp"
#include "vg.hpp"

namespace vg {
namespace unittest {

using namespace std;

namespace {

/// A graph of three nested snarls: 1 to 8 holds 2 to 7, which holds 3 to 5. Node 6 is the other
/// branch of 2 to 7, and node 4 the other branch of 3 to 5.
struct NestedGraph {
    VG graph;
    vector<Node*> nodes;

    NestedGraph() {
        for (const string& sequence : {"GCA", "T", "G", "CTGA", "GCA", "T", "G", "CTGA"}) {
            nodes.push_back(graph.create_node(sequence));
        }
        auto edge = [&](int from, int to) { graph.create_edge(nodes[from - 1], nodes[to - 1]); };
        edge(1, 2);
        edge(1, 8);
        edge(2, 3);
        edge(2, 6);
        edge(3, 4);
        edge(3, 5);
        edge(4, 5);
        edge(5, 7);
        edge(6, 7);
        edge(7, 8);
    }

    handle_t handle(nid_t id, bool reverse = false) {
        return graph.get_handle(id, reverse);
    }
};

/// Get the children of `parent`, canonicalized.
vector<net_handle_t> children(const SnarlDecomposition& decomposition, const net_handle_t& parent) {
    vector<net_handle_t> out;
    decomposition.for_each_child(parent, [&](const net_handle_t& child) {
        out.push_back(decomposition.canonical(child));
    });
    return out;
}

/// Find the snarl of `decomposition` bounded by nodes `a` and `b`, by searching down from the root.
///
/// The decomposition may orient the snarl either way. Return the traversal of it that starts at
/// `a`.
net_handle_t find_snarl(const SnarlDecomposition& decomposition, const HandleGraph& graph, nid_t a,
                        nid_t b) {
    net_handle_t found;
    bool have = false;
    function<void(const net_handle_t&)> search = [&](const net_handle_t& here) {
        decomposition.for_each_child(here, [&](const net_handle_t& child) {
            if (decomposition.is_snarl(child)) {
                const nid_t start = graph.get_id(decomposition.get_handle(
                    decomposition.get_bound(child, false, true), &graph));
                const nid_t end = graph.get_id(decomposition.get_handle(
                    decomposition.get_bound(child, true, false), &graph));
                if (start == a && end == b) {
                    found = decomposition.starts_at(child) == SnarlDecomposition::START
                                ? child : decomposition.flip(child);
                    have = true;
                } else if (start == b && end == a) {
                    found = decomposition.starts_at(child) == SnarlDecomposition::END
                                ? child : decomposition.flip(child);
                    have = true;
                }
            }
            if (!decomposition.is_node(child)) {
                search(child);
            }
        });
    };
    search(decomposition.get_root());
    REQUIRE(have);
    return found;
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

TEST_CASE("remove_edits removes the edits at and below an element, and keeps the others",
          "[edit_snarl_decomposition_overlay]") {
    NestedGraph nested;
    on_both(nested, [&](const SnarlDecomposition&, EditSnarlDecompositionOverlay& overlay) {
        const net_handle_t middle = find_snarl(overlay, nested.graph, 2, 7);
        const net_handle_t inside = overlay.add_edit(nested.handle(4), 0, nested.handle(4), 2);
        const net_handle_t outside = overlay.add_edit(nested.handle(1), 0, nested.handle(1), 1);
        const net_handle_t node4 = overlay.canonical(overlay.get_net(nested.handle(4), &nested.graph));
        const net_handle_t node1 = overlay.canonical(overlay.get_net(nested.handle(1), &nested.graph));

        overlay.remove_edits(middle);
        REQUIRE(children(overlay, node4).empty());
        REQUIRE(children(overlay, node1) == vector<net_handle_t>{outside});
        // The edit that was not removed still answers through its handle.
        REQUIRE(overlay.get_site_end_offset(outside) == 1);
        // A new edit may take the removed one's place, and the kept one is unchanged.
        const net_handle_t again = overlay.add_edit(nested.handle(4), 1, nested.handle(4), 3);
        REQUIRE(children(overlay, node4) == vector<net_handle_t>{again});
        REQUIRE(overlay.get_site_end_offset(again) == 3);
        REQUIRE(overlay.get_site_end_offset(outside) == 1);
        (void) inside;
    });
}

TEST_CASE("remove_edit removes one edit and keeps its siblings",
          "[edit_snarl_decomposition_overlay]") {
    NestedGraph nested;
    on_both(nested, [&](const SnarlDecomposition&, EditSnarlDecompositionOverlay& overlay) {
        const net_handle_t first = overlay.add_edit(nested.handle(4), 0, nested.handle(4), 1);
        const net_handle_t second = overlay.add_edit(nested.handle(4), 2, nested.handle(4), 3);
        const net_handle_t node4 = overlay.canonical(overlay.get_net(nested.handle(4), &nested.graph));
        overlay.remove_edit(first);
        REQUIRE(children(overlay, node4) == vector<net_handle_t>{second});
        REQUIRE(overlay.get_site_start_offset(second) == 2);
        overlay.remove_edit(second);
        REQUIRE(children(overlay, node4).empty());
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
