/// \file incremental_subgraph.cpp
///  
/// Unit tests for the IncrementalSubgraph
///

#include <iostream>

#include "../incremental_subgraph.hpp"
#include "catch.hpp"

#include <bdsg/hash_graph.hpp>

namespace vg {
namespace unittest {
using namespace std;

using bdsg::HashGraph;

bool verify_incremental_max_distance(const IncrementalSubgraph& subgraph) {
    
    vector<int64_t> dists(subgraph.get_node_count(), numeric_limits<int64_t>::min());
    dists[0] = subgraph.max_distance_from_start(subgraph.handle_at_order(0));
    for (int i = 0; i < subgraph.get_node_count(); ++i) {
        handle_t h = subgraph.handle_at_order(i);
        int64_t dist_thru = dists[i] + subgraph.get_length(h);
        subgraph.follow_edges(h, subgraph.extracting_left(), [&](const handle_t& next) {
            size_t j = subgraph.order_of(next);
            dists[j] = max(dists[j], dist_thru);
        });
    }
    for (int i = 0; i < subgraph.get_node_count(); ++i) {
        if (subgraph.max_distance_from_start(subgraph.handle_at_order(i)) != dists[i]) {
            return false;
        }
    }
    return true;
}

TEST_CASE("IncrementalSubgraph behaves correctly on a trivial example",
          "[incremental]") {
    
    HashGraph graph;
    handle_t h0 = graph.create_handle("ACC");
    handle_t h1 = graph.create_handle("GCT");
    handle_t h2 = graph.create_handle("GAAC");
    
    graph.create_edge(h0, h1);
    graph.create_edge(h1, h2);
    
    pos_t start = make_pos_t(graph.get_id(h0), false, 2);
    
    IncrementalSubgraph subgraph(graph, start, false);
    
    // basic queries work
    REQUIRE(subgraph.get_node_count() == 1);
    handle_t s0 = subgraph.handle_at_order(0);
    REQUIRE(subgraph.get_sequence(s0) == graph.get_sequence(h0));
    REQUIRE(subgraph.get_underlying_handle(s0) == h0);
    REQUIRE(subgraph.get_underlying_handle(subgraph.flip(s0)) == graph.flip(h0));
    REQUIRE(subgraph.min_distance_from_start(s0) == -2);
    REQUIRE(subgraph.max_distance_from_start(s0) == -2);
    REQUIRE(subgraph.has_node(subgraph.get_id(s0)));
    REQUIRE(subgraph.get_handle(subgraph.get_id(s0)) == s0);
    REQUIRE(subgraph.order_of(s0) == 0);
    REQUIRE(verify_incremental_max_distance(subgraph));
    
    // we shouldn't have any edges yet
    subgraph.follow_edges(s0, false, [&](const handle_t& next) {
        REQUIRE(false);
    });
    subgraph.follow_edges(s0, true, [&](const handle_t& next) {
        REQUIRE(false);
    });
    subgraph.follow_edges(subgraph.flip(s0), false, [&](const handle_t& next) {
        REQUIRE(false);
    });
    subgraph.follow_edges(subgraph.flip(s0), true, [&](const handle_t& next) {
        REQUIRE(false);
    });
    
    // we can extend the graph
    REQUIRE(subgraph.is_extendable());
    handle_t s1 = subgraph.extend();
    REQUIRE(subgraph.get_underlying_handle(s1) == h1);
    REQUIRE(subgraph.min_distance_from_start(s1) == 1);
    REQUIRE(subgraph.max_distance_from_start(s1) == 1);
    REQUIRE(subgraph.handle_at_order(1) == s1);
    REQUIRE(subgraph.order_of(s1) == 1);
    REQUIRE(verify_incremental_max_distance(subgraph));
    
    // edges are added correctly and can traverse in all directions
    int count = 0;
    subgraph.follow_edges(s0, false, [&](const handle_t& next) {
        REQUIRE(next == s1);
        count++;
    });
    REQUIRE(count == 1);
    count = 0;
    subgraph.follow_edges(s0, true, [&](const handle_t& next) {
        REQUIRE(false);
    });
    subgraph.follow_edges(subgraph.flip(s0), false, [&](const handle_t& next) {
        REQUIRE(false);
    });
    subgraph.follow_edges(subgraph.flip(s0), true, [&](const handle_t& next) {
        REQUIRE(next == subgraph.flip(s1));
        count++;
    });
    REQUIRE(count == 1);
    count = 0;
    subgraph.follow_edges(s1, false, [&](const handle_t& next) {
        REQUIRE(false);
    });
    subgraph.follow_edges(s1, true, [&](const handle_t& next) {
        REQUIRE(next == s0);
        count++;
    });
    REQUIRE(count == 1);
    count = 0;
    subgraph.follow_edges(subgraph.flip(s1), false, [&](const handle_t& next) {
        REQUIRE(next == subgraph.flip(s0));
        count++;
    });
    REQUIRE(count == 1);
    count = 0;
    subgraph.follow_edges(subgraph.flip(s1), true, [&](const handle_t& next) {
        REQUIRE(false);
    });
    
    // we can extend until hitting the end
    
    REQUIRE(subgraph.is_extendable());
    handle_t s2 = subgraph.extend();
    REQUIRE(subgraph.get_underlying_handle(s2) == h2);
    REQUIRE(subgraph.min_distance_from_start(s2) == 4);
    REQUIRE(subgraph.max_distance_from_start(s2) == 4);
    REQUIRE(subgraph.handle_at_order(2) == s2);
    REQUIRE(subgraph.order_of(s2) == 2);
    REQUIRE(verify_incremental_max_distance(subgraph));
    
    REQUIRE(!subgraph.is_extendable());
}

TEST_CASE("IncrementalSubgraph can extract in either direction on either strand",
          "[incremental]") {
    
    HashGraph graph;
    handle_t h0 = graph.create_handle("AC");
    handle_t h1 = graph.create_handle("GCT");
    
    graph.create_edge(h0, h1);
    
    handle_t s0, s1;
    int count = 0;
    
    {
        // extract to left on forward strand
        pos_t start = make_pos_t(graph.get_id(h1), false, 3);
        IncrementalSubgraph subgraph(graph, start, true);
        
        s0 = subgraph.handle_at_order(0);
        s1 = subgraph.extend();
        
        REQUIRE(subgraph.get_underlying_handle(s0) == h1);
        REQUIRE(subgraph.get_underlying_handle(s1) == h0);
        REQUIRE(verify_incremental_max_distance(subgraph));
        
        count = 0;
        subgraph.follow_edges(s0, false, [&](const handle_t& next) {
            count++;
        });
        REQUIRE(count == 0);
        count = 0;
        subgraph.follow_edges(s0, true, [&](const handle_t& next) {
            REQUIRE(next == s1);
            count++;
        });
        REQUIRE(count == 1);
        count = 0;
        subgraph.follow_edges(subgraph.flip(s0), false, [&](const handle_t& next) {
            REQUIRE(next == subgraph.flip(s1));
            count++;
        });
        REQUIRE(count == 1);
        count = 0;
        subgraph.follow_edges(subgraph.flip(s0), true, [&](const handle_t& next) {
            count++;
        });
        REQUIRE(count == 0);
        count = 0;
        subgraph.follow_edges(s1, false, [&](const handle_t& next) {
            REQUIRE(next == s0);
            count++;
        });
        REQUIRE(count == 1);
        count = 0;
        subgraph.follow_edges(s1, true, [&](const handle_t& next) {
            count++;
        });
        REQUIRE(count == 0);
        count = 0;
        subgraph.follow_edges(subgraph.flip(s1), false, [&](const handle_t& next) {
            count++;
        });
        REQUIRE(count == 0);
        count = 0;
        subgraph.follow_edges(subgraph.flip(s1), true, [&](const handle_t& next) {
            REQUIRE(next == subgraph.flip(s0));
            count++;
        });
        REQUIRE(count == 1);
        count = 0;
    }
    
    {
        // extract to right on reverse strand
        pos_t start = make_pos_t(graph.get_id(h1), true, 1);
        IncrementalSubgraph subgraph(graph, start, false);
        
        s0 = subgraph.handle_at_order(0);
        s1 = subgraph.extend();
        
        REQUIRE(subgraph.get_underlying_handle(s0) == graph.flip(h1));
        REQUIRE(subgraph.get_underlying_handle(s1) == graph.flip(h0));
        REQUIRE(verify_incremental_max_distance(subgraph));
        
        count = 0;
        subgraph.follow_edges(s0, false, [&](const handle_t& next) {
            REQUIRE(next == s1);
            count++;
        });
        REQUIRE(count == 1);
        count = 0;
        subgraph.follow_edges(s0, true, [&](const handle_t& next) {
            count++;
        });
        REQUIRE(count == 0);
        count = 0;
        subgraph.follow_edges(subgraph.flip(s0), false, [&](const handle_t& next) {
            count++;
        });
        REQUIRE(count == 0);
        count = 0;
        subgraph.follow_edges(subgraph.flip(s0), true, [&](const handle_t& next) {
            REQUIRE(next == subgraph.flip(s1));
            count++;
        });
        REQUIRE(count == 1);
        count = 0;
        subgraph.follow_edges(s1, false, [&](const handle_t& next) {
            count++;
        });
        REQUIRE(count == 0);
        count = 0;
        subgraph.follow_edges(s1, true, [&](const handle_t& next) {
            REQUIRE(next == s0);
            count++;
        });
        REQUIRE(count == 1);
        count = 0;
        subgraph.follow_edges(subgraph.flip(s1), false, [&](const handle_t& next) {
            REQUIRE(next == subgraph.flip(s0));
            count++;
        });
        REQUIRE(count == 1);
        count = 0;
        subgraph.follow_edges(subgraph.flip(s1), true, [&](const handle_t& next) {
            count++;
        });
        REQUIRE(count == 0);
        count = 0;
    }
    
    {
        // extract to left on reverse strand
        pos_t start = make_pos_t(graph.get_id(h0), true, 1);
        IncrementalSubgraph subgraph(graph, start, true);
        
        s0 = subgraph.handle_at_order(0);
        s1 = subgraph.extend();
        
        REQUIRE(subgraph.get_underlying_handle(s0) == graph.flip(h0));
        REQUIRE(subgraph.get_underlying_handle(s1) == graph.flip(h1));
        REQUIRE(verify_incremental_max_distance(subgraph));
        
        count = 0;
        subgraph.follow_edges(s0, false, [&](const handle_t& next) {
            count++;
        });
        REQUIRE(count == 0);
        count = 0;
        subgraph.follow_edges(s0, true, [&](const handle_t& next) {
            REQUIRE(next == s1);
            count++;
        });
        REQUIRE(count == 1);
        count = 0;
        subgraph.follow_edges(subgraph.flip(s0), false, [&](const handle_t& next) {
            REQUIRE(next == subgraph.flip(s1));
            count++;
        });
        REQUIRE(count == 1);
        count = 0;
        subgraph.follow_edges(subgraph.flip(s0), true, [&](const handle_t& next) {
            count++;
        });
        REQUIRE(count == 0);
        count = 0;
        subgraph.follow_edges(s1, false, [&](const handle_t& next) {
            REQUIRE(next == s0);
            count++;
        });
        REQUIRE(count == 1);
        count = 0;
        subgraph.follow_edges(s1, true, [&](const handle_t& next) {
            count++;
        });
        REQUIRE(count == 0);
        count = 0;
        subgraph.follow_edges(subgraph.flip(s1), false, [&](const handle_t& next) {
            count++;
        });
        REQUIRE(count == 0);
        count = 0;
        subgraph.follow_edges(subgraph.flip(s1), true, [&](const handle_t& next) {
            REQUIRE(next == subgraph.flip(s0));
            count++;
        });
        REQUIRE(count == 1);
        count = 0;
    }
}

TEST_CASE("IncrementalSubgraph can handle branching paths and bubbles correctly",
          "[incremental]") {
    
    SECTION("Two deletions") {

        HashGraph graph;
        handle_t h0 = graph.create_handle("AC");
        handle_t h1 = graph.create_handle("A");
        handle_t h2 = graph.create_handle("T");
        handle_t h3 = graph.create_handle("C");
        handle_t h4 = graph.create_handle("G");

        graph.create_edge(h0, h1);
        graph.create_edge(h0, h2);
        graph.create_edge(h1, h2);
        graph.create_edge(h2, h4);
        graph.create_edge(h2, h3);
        graph.create_edge(h3, h4);

        pos_t start = make_pos_t(graph.get_id(h0), false, 1);

        IncrementalSubgraph subgraph(graph, start, false);
        handle_t s0 = subgraph.handle_at_order(0);
        handle_t s1 = subgraph.extend();
        handle_t s2 = subgraph.extend();
        handle_t s3 = subgraph.extend();
        handle_t s4 = subgraph.extend();

        REQUIRE(!subgraph.is_extendable());

        // did we get everything in topological order
        REQUIRE(subgraph.get_underlying_handle(s0) == h0);
        REQUIRE(subgraph.get_underlying_handle(s1) == h1);
        REQUIRE(subgraph.get_underlying_handle(s2) == h2);
        REQUIRE(subgraph.get_underlying_handle(s3) == h3);
        REQUIRE(subgraph.get_underlying_handle(s4) == h4);
    }

    SECTION("Unreachable branch") {

        HashGraph graph;
        handle_t h0 = graph.create_handle("C");
        handle_t h1 = graph.create_handle("A");
        handle_t h2 = graph.create_handle("T");

        graph.create_edge(h0, h2);
        graph.create_edge(h1, h2);

        pos_t start = make_pos_t(graph.get_id(h0), false, 0);

        IncrementalSubgraph subgraph(graph, start, false);
        handle_t s0 = subgraph.handle_at_order(0);
        handle_t s1 = subgraph.extend();

        REQUIRE(!subgraph.is_extendable());

        REQUIRE(subgraph.get_underlying_handle(s0) == h0);
        REQUIRE(subgraph.get_underlying_handle(s1) == h2);
        REQUIRE(verify_incremental_max_distance(subgraph));
    }

    SECTION("Respects max distance") {

        HashGraph graph;
        handle_t h0 = graph.create_handle("ACC");
        handle_t h1 = graph.create_handle("AA");
        handle_t h2 = graph.create_handle("GT");

        graph.create_edge(h0, h1);
        graph.create_edge(h1, h2);

        pos_t start = make_pos_t(graph.get_id(h0), false, 1);

        IncrementalSubgraph subgraph(graph, start, false, 4);
        handle_t s0 = subgraph.handle_at_order(0);
        handle_t s1 = subgraph.extend();

        REQUIRE(!subgraph.is_extendable());

        REQUIRE(subgraph.get_underlying_handle(s0) == h0);
        REQUIRE(subgraph.get_underlying_handle(s1) == h1);
        REQUIRE(verify_incremental_max_distance(subgraph));

    }
    
    SECTION("Respects max distance along different paths") {
        
        HashGraph graph;
        handle_t h0 = graph.create_handle("ACC");
        handle_t h1 = graph.create_handle("G");
        handle_t h2 = graph.create_handle("G");
        handle_t h3 = graph.create_handle("AGGT");
        handle_t h4 = graph.create_handle("TTCA");
        handle_t h5 = graph.create_handle("CG");
        
        graph.create_edge(h0, h2);
        graph.create_edge(h0, h3);
        graph.create_edge(h1, h2);
        graph.create_edge(h2, h5);
        graph.create_edge(h3, h4);
        graph.create_edge(h4, h5);
        
        pos_t start = make_pos_t(graph.get_id(h0), false, 1);
        
        IncrementalSubgraph subgraph(graph, start, false);
        
        handle_t s0 = subgraph.handle_at_order(0);
        handle_t s1 = subgraph.extend();
        handle_t s2 = subgraph.extend();
        handle_t s3 = subgraph.extend();
        handle_t s4 = subgraph.extend();
        
        REQUIRE(!subgraph.is_extendable());
        
        // did we assign the final node the right distance, even
        // though we would find the longer path first?
        REQUIRE(subgraph.get_underlying_handle(s4) == h5);
        REQUIRE(subgraph.min_distance_from_start(s4) == 3);
        REQUIRE(verify_incremental_max_distance(subgraph));
        
    }
}

TEST_CASE("IncrementalSubgraph remains acyclic when extracting from within a self-loop",
          "[incremental]") {
    
    HashGraph graph;
    handle_t h0 = graph.create_handle("A");
    handle_t h1 = graph.create_handle("AAAAAA");
    
    graph.create_edge(h0, h0);
    graph.create_edge(h0, h1);
    
    pos_t start = make_pos_t(graph.get_id(h0), false, 0);
    
    IncrementalSubgraph subgraph(graph, start, false, 10);
    
    REQUIRE(handlealgs::is_directed_acyclic(&subgraph));
    while (subgraph.is_extendable()) {
        subgraph.extend();
        REQUIRE(handlealgs::is_directed_acyclic(&subgraph));
    }
    
//    subgraph.for_each_handle([&](const handle_t& h) {
//        cerr << subgraph.get_id(h) << " " << subgraph.get_sequence(h) << endl;
//        subgraph.follow_edges(h, true, [&](const handle_t& n) {
//            cerr << "\t" << subgraph.get_id(n) << " <-" << endl;
//        });
//        subgraph.follow_edges(h, false, [&](const handle_t& n) {
//            cerr << "\t-> " << subgraph.get_id(n) << endl;
//        });
//    });
}

}
}
        
