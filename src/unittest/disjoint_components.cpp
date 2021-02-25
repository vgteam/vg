#include "catch.hpp"
#include "algorithms/disjoint_components.hpp"
#include "handle.hpp"

#include <bdsg/hash_graph.hpp>

namespace vg {
namespace unittest {

using namespace std;

TEST_CASE("disjoint_components produces correct results", "[handle]") {
    
    bdsg::HashGraph graph;
    handle_t h1 = graph.create_handle("A");
    handle_t h2 = graph.create_handle("C");
    handle_t h3 = graph.create_handle("T");
    handle_t h4 = graph.create_handle("G");
    
    graph.create_edge(h1, h2);
    graph.create_edge(h3, h4);
    
    bdsg::HashGraph comp1, comp2;
    handle_t c1 = comp1.create_handle(graph.get_sequence(h1), graph.get_id(h1));
    handle_t c2 = comp1.create_handle(graph.get_sequence(h2), graph.get_id(h2));
    handle_t c3 = comp2.create_handle(graph.get_sequence(h3), graph.get_id(h3));
    handle_t c4 = comp2.create_handle(graph.get_sequence(h4), graph.get_id(h4));
    
    comp1.create_edge(c1, c2);
    comp2.create_edge(c3, c4);
    
    auto disj_comps = algorithms::disjoint_components(graph);
    
    REQUIRE(disj_comps.size() == 2);
    
    if (handlealgs::are_equivalent(&comp1, &disj_comps.front())) {
        REQUIRE(handlealgs::are_equivalent(&comp1, &disj_comps.front()));
        REQUIRE(handlealgs::are_equivalent(&comp2, &disj_comps.back()));
    }
    else {
        REQUIRE(handlealgs::are_equivalent(&comp2, &disj_comps.front()));
        REQUIRE(handlealgs::are_equivalent(&comp1, &disj_comps.back()));
    }
}
}
}
