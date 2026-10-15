/// \file follow_net_elements.cpp
///
/// Unit tests for follow_net_elements(), which steps along a chain one net
/// graph element at a time.

#include "catch.hpp"
#include "../follow_net_elements.hpp"
#include "../integrated_snarl_finder.hpp"
#include "../snarl_manager_decomposition.hpp"

#include <bdsg/hash_graph.hpp>

#include <vector>

namespace vg::unittest {

using handlegraph::SnarlDecomposition;

TEST_CASE("follow_net_elements() from a chain node gives the next snarl once as START_END",
          "[follow_net_elements]") {
    // The chain is node 1, the snarl from 1 to 4, and node 4. In the snarl,
    // an edge from the end of 2 into the end of 3 lets a walk from 1 turn
    // around, so follow_net_edges() gives the snarl START_END and
    // START_START.
    bdsg::HashGraph graph;
    handle_t n1 = graph.create_handle("A");
    handle_t n2 = graph.create_handle("C");
    handle_t n3 = graph.create_handle("G");
    handle_t n4 = graph.create_handle("T");
    graph.create_edge(n1, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n4);
    graph.create_edge(n3, n4);
    graph.create_edge(n2, graph.flip(n3));
    IntegratedSnarlFinder finder(graph);
    SnarlManager manager = finder.find_snarls();
    SnarlManagerDecomposition view(manager, graph);

    net_handle_t node1 = view.get_net(n1, &graph);
    std::vector<net_handle_t> edges;
    view.follow_net_edges(node1, &graph, false, [&](const net_handle_t& next) {
        edges.push_back(next);
    });
    REQUIRE(edges.size() == 2);

    std::vector<net_handle_t> elements;
    follow_net_elements(view, node1, &graph, false, [&](const net_handle_t& next) {
        elements.push_back(next);
    });
    REQUIRE(elements.size() == 1);
    REQUIRE(view.is_snarl(elements[0]));
    REQUIRE(view.starts_at(elements[0]) == SnarlDecomposition::START);
    REQUIRE(view.ends_at(elements[0]) == SnarlDecomposition::END);
}

}
