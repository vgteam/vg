/// \file decomposition_neighbors.cpp
///
/// Unit tests for for_each_neighbor(), which gives one traversal of each
/// neighboring net graph element.

#include "catch.hpp"
#include "../decomposition_neighbors.hpp"
#include "../integrated_snarl_finder.hpp"
#include "../snarl_manager_decomposition.hpp"

#include <bdsg/hash_graph.hpp>

#include <vector>

namespace vg::unittest {

using handlegraph::SnarlDecomposition;

TEST_CASE("for_each_neighbor() from a chain node gives the next snarl once as START_END",
          "[for_each_neighbor]") {
    // The chain is node 1, the snarl from 1 to 4, and node 4. In the snarl,
    // an edge from the end of 2 into the end of 3 lets a walk from 1 turn
    // around, so follow_net_edges() gives the snarl START_END and
    // START_START. for_each_neighbor() gives only the START_END traversal.
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
    std::vector<net_handle_t> adjacent;
    view.follow_net_edges(node1, &graph, false, [&](const net_handle_t& next) {
        adjacent.push_back(next);
    });
    REQUIRE(adjacent.size() == 2);

    std::vector<net_handle_t> neighbors;
    for_each_neighbor(view, node1, &graph, false, [&](const net_handle_t& next) {
        neighbors.push_back(next);
    });
    REQUIRE(neighbors.size() == 1);
    REQUIRE(view.is_snarl(neighbors[0]));
    REQUIRE(view.starts_at(neighbors[0]) == SnarlDecomposition::START);
    REQUIRE(view.ends_at(neighbors[0]) == SnarlDecomposition::END);
}

}
