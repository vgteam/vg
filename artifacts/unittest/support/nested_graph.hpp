#ifndef VG_UNITTEST_NESTED_GRAPH_HPP_INCLUDED
#define VG_UNITTEST_NESTED_GRAPH_HPP_INCLUDED

/** \file
 * A small graph of three nested snarls, and a way to find a snarl of a decomposition by its
 * boundary nodes, for tests of code that works on the snarl tree.
 */

#include <functional>
#include <vector>

#include <handlegraph/snarl_decomposition.hpp>

#include "../../../src/unittest/catch.hpp"
#include "../../../src/vg.hpp"

namespace vg {
namespace unittest {

using namespace std;
using handlegraph::net_handle_t;
using handlegraph::SnarlDecomposition;

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

    handle_t handle(nid_t id, bool reverse = false) const {
        return graph.get_handle(id, reverse);
    }
};

/// Find the snarl of `decomposition` bounded by nodes `a` and `b`, by searching down from the root.
///
/// The decomposition may orient the snarl either way. Return the traversal of it that starts at
/// `a`.
inline net_handle_t find_snarl(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                               nid_t a, nid_t b) {
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

}
}

#endif
