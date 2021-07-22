/// \file dagify.cpp
///  
/// unit tests for the handle based dagify algorithm
///

#include "../handle.hpp"
#include "../algorithms/locally_expand_graph.hpp"
#include "catch.hpp"

#include "bdsg/hash_graph.hpp"

namespace vg {
namespace unittest {

TEST_CASE("locally_expand_graph produces expected results", "[algorithms][handle][localexpand]") {
    
    bdsg::HashGraph parent;
    bdsg::HashGraph subgraph;
    
    handle_t p1 = parent.create_handle("GAT");
    handle_t p2 = parent.create_handle("TA");
    handle_t p3 = parent.create_handle("CA");
    handle_t p4 = parent.create_handle("T");
    handle_t p5 = parent.create_handle("GGG");
    handle_t p6 = parent.create_handle("A");
    
    parent.create_edge(p1, p2);
    parent.create_edge(p2, p3);
    parent.create_edge(p3, p4);
    parent.create_edge(p3, p5);
    parent.create_edge(p4, p5);
    parent.create_edge(p5, p6);
    
    handle_t s2 = subgraph.create_handle(parent.get_sequence(p2),
                                         parent.get_id(p2));
    
    SECTION("locally_expand_graph can expand forward") {

        algorithms::locally_expand_graph(parent, subgraph, s2, 1);

        bdsg::HashGraph expected;

        handle_t e2 = expected.create_handle(parent.get_sequence(p2),
                                             parent.get_id(p2));
        handle_t e3 = expected.create_handle(parent.get_sequence(p3),
                                             parent.get_id(p3));
        expected.create_edge(e2, e3);


        REQUIRE(handlealgs::are_equivalent(&subgraph, &expected));

    }

    SECTION("locally_expand_graph can expand backwards") {

        algorithms::locally_expand_graph(parent, subgraph, subgraph.flip(s2), 1);

        bdsg::HashGraph expected;

        handle_t e1 = expected.create_handle(parent.get_sequence(p1),
                                             parent.get_id(p1));
        handle_t e2 = expected.create_handle(parent.get_sequence(p2),
                                             parent.get_id(p2));
        expected.create_edge(e1, e2);



        REQUIRE(handlealgs::are_equivalent(&subgraph, &expected));

    }
    
    SECTION("locally_expand_graph expands along shortest paths") {

        algorithms::locally_expand_graph(parent, subgraph, s2, 3);

        bdsg::HashGraph expected;

        handle_t e2 = expected.create_handle(parent.get_sequence(p2),
                                             parent.get_id(p2));
        handle_t e3 = expected.create_handle(parent.get_sequence(p3),
                                             parent.get_id(p3));
        handle_t e4 = expected.create_handle(parent.get_sequence(p4),
                                             parent.get_id(p4));
        handle_t e5 = expected.create_handle(parent.get_sequence(p5),
                                             parent.get_id(p5));

        expected.create_edge(e2, e3);
        expected.create_edge(e3, e4);
        expected.create_edge(e3, e5);

//        cerr << "expanded" << endl;
//        subgraph.for_each_handle([&](const handle_t& h) {
//            cerr << subgraph.get_id(h) << ": " << subgraph.get_sequence(h) << endl;
//            subgraph.follow_edges(h, true, [&](const handle_t& n) {
//                cerr << "\t" << subgraph.get_id(n) << (subgraph.get_is_reverse(n) ? "-" : "+") << " <-" << endl;
//            });
//            subgraph.follow_edges(h, false, [&](const handle_t& n) {
//                cerr << "\t-> " << subgraph.get_id(n) << (subgraph.get_is_reverse(n) ? "-" : "+") << endl;
//            });
//        });
//
//        cerr << "expected" << endl;
//        expected.for_each_handle([&](const handle_t& h) {
//            cerr << expected.get_id(h) << ": " << expected.get_sequence(h) << endl;
//            expected.follow_edges(h, true, [&](const handle_t& n) {
//                cerr << "\t" << expected.get_id(n) << (expected.get_is_reverse(n) ? "-" : "+") << " <-" << endl;
//            });
//            expected.follow_edges(h, false, [&](const handle_t& n) {
//                cerr << "\t-> " << expected.get_id(n) << (expected.get_is_reverse(n) ? "-" : "+") << endl;
//            });
//        });
        
        REQUIRE(handlealgs::are_equivalent(&subgraph, &expected));

    }
}
}
}
