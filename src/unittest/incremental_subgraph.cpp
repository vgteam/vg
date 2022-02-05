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

TEST_CASE("IncrementalSubgraph does not crash on a particular gnarly graph",
          "[incremental]") {
    
    // encountered this one in the wild
    
    bdsg::HashGraph graph;
    
    graph.create_handle("GA", 134336223);
    graph.create_handle("T", 134336281);
    graph.create_handle("T", 134336248);
    graph.create_handle("T", 134336291);
    graph.create_handle("G", 134336221);
    graph.create_handle("A", 134336259);
    graph.create_handle("C", 134336230);
    graph.create_handle("T", 134336346);
    graph.create_handle("C", 134336181);
    graph.create_handle("A", 134336289);
    graph.create_handle("G", 134336257);
    graph.create_handle("TAAGAC", 134336256);
    graph.create_handle("T", 134336254);
    graph.create_handle("G", 134336208);
    graph.create_handle("C", 134336284);
    graph.create_handle("GA", 134336307);
    graph.create_handle("G", 134336233);
    graph.create_handle("GCGCGTGTCGTGGGGCCGGCGGGCGGCGGGGA", 134336195);
    graph.create_handle("C", 134336323);
    graph.create_handle("C", 134336260);
    graph.create_handle("C", 134336327);
    graph.create_handle("A", 134336179);
    graph.create_handle("T", 134336358);
    graph.create_handle("T", 134336229);
    graph.create_handle("C", 134336240);
    graph.create_handle("C", 134336340);
    graph.create_handle("CACAGT", 134336232);
    graph.create_handle("A", 134336219);
    graph.create_handle("C", 134336237);
    graph.create_handle("C", 134336186);
    graph.create_handle("C", 134336308);
    graph.create_handle("A", 134336274);
    graph.create_handle("A", 134336255);
    graph.create_handle("A", 134336330);
    graph.create_handle("G", 134336177);
    graph.create_handle("T", 134336350);
    graph.create_handle("T", 134336296);
    graph.create_handle("A", 134336214);
    graph.create_handle("T", 134336316);
    graph.create_handle("C", 134336220);
    graph.create_handle("GTT", 134336213);
    graph.create_handle("C", 134336302);
    graph.create_handle("A", 134336188);
    graph.create_handle("A", 134336329);
    graph.create_handle("C", 134336258);
    graph.create_handle("G", 134336222);
    graph.create_handle("TGTCGGCGGGCGCGGGGGCGGTTCTCGGCGGC", 134336197);
    graph.create_handle("G", 134336319);
    graph.create_handle("C", 134336278);
    graph.create_handle("T", 134336324);
    graph.create_handle("A", 134336343);
    graph.create_handle("C", 134336242);
    graph.create_handle("CAT", 134336317);
    graph.create_handle("C", 134336348);
    graph.create_handle("T", 134336180);
    graph.create_handle("A", 134336269);
    graph.create_handle("C", 134336351);
    graph.create_handle("G", 134336293);
    graph.create_handle("G", 134336294);
    graph.create_handle("GAGAAAGAGAAAGAAGGGCGTGTCGTTGGTGT", 134336194);
    graph.create_handle("C", 134336311);
    graph.create_handle("C", 134336314);
    graph.create_handle("G", 134336353);
    graph.create_handle("GG", 134336309);
    graph.create_handle("T", 134336357);
    graph.create_handle("G", 134336344);
    graph.create_handle("T", 134336306);
    graph.create_handle("G", 134336342);
    graph.create_handle("C", 134336320);
    graph.create_handle("C", 134336175);
    graph.create_handle("C", 134336290);
    graph.create_handle("TC", 134336206);
    graph.create_handle("TC", 134336272);
    graph.create_handle("C", 134336321);
    graph.create_handle("T", 134336253);
    graph.create_handle("C", 134336250);
    graph.create_handle("C", 134336271);
    graph.create_handle("T", 134336246);
    graph.create_handle("T", 134336241);
    graph.create_handle("C", 134336273);
    graph.create_handle("AG", 134336322);
    graph.create_handle("A", 134336224);
    graph.create_handle("T", 134336204);
    graph.create_handle("C", 134336345);
    graph.create_handle("C", 134336185);
    graph.create_handle("A", 134336265);
    graph.create_handle("C", 134336182);
    graph.create_handle("T", 134336202);
    graph.create_handle("GGCCGCGAGAGCCGGAGAACTCGGGAGGGAGA", 134336192);
    graph.create_handle("C", 134336249);
    graph.create_handle("GTCGCGGCGGGTCTGGGGG", 134336198);
    graph.create_handle("G", 134336189);
    graph.create_handle("G", 134336207);
    graph.create_handle("C", 134336292);
    graph.create_handle("A", 134336199);
    graph.create_handle("G", 134336336);
    graph.create_handle("CGGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA", 134336193);
    graph.create_handle("C", 134336212);
    graph.create_handle("T", 134336287);
    graph.create_handle("GG", 134336245);
    graph.create_handle("C", 134336354);
    graph.create_handle("C", 134336339);
    graph.create_handle("A", 134336216);
    graph.create_handle("CAGACAGACAGACAGACAGACA", 134336295);
    graph.create_handle("G", 134336215);
    graph.create_handle("T", 134336337);
    graph.create_handle("TTG", 134336201);
    graph.create_handle("T", 134336347);
    graph.create_handle("T", 134336298);
    graph.create_handle("C", 134336190);
    graph.create_handle("T", 134336310);
    graph.create_handle("G", 134336332);
    graph.create_handle("C", 134336263);
    graph.create_handle("G", 134336266);
    graph.create_handle("C", 134336187);
    graph.create_handle("C", 134336217);
    graph.create_handle("T", 134336267);
    graph.create_handle("T", 134336178);
    graph.create_handle("C", 134336239);
    graph.create_handle("G", 134336300);
    graph.create_handle("G", 134336328);
    graph.create_handle("C", 134336333);
    graph.create_handle("A", 134336184);
    graph.create_handle("CT", 134336279);
    graph.create_handle("C", 134336297);
    graph.create_handle("T", 134336275);
    graph.create_handle("A", 134336334);
    graph.create_handle("C", 134336228);
    graph.create_handle("G", 134336205);
    graph.create_handle("G", 134336252);
    graph.create_handle("C", 134336282);
    graph.create_handle("CT", 134336335);
    graph.create_handle("C", 134336326);
    graph.create_handle("GCGGTCCCCGGCCGCGGCCCCGACGGCGTGGG", 134336196);
    graph.create_handle("C", 134336268);
    graph.create_handle("G", 134336325);
    graph.create_handle("C", 134336276);
    graph.create_handle("C", 134336191);
    graph.create_handle("G", 134336211);
    graph.create_handle("C", 134336299);
    graph.create_handle("C", 134336238);
    graph.create_handle("C", 134336277);
    graph.create_handle("C", 134336352);
    graph.create_handle("G", 134336359);
    graph.create_handle("G", 134336210);
    graph.create_handle("G", 134336244);
    graph.create_handle("A", 134336356);
    graph.create_handle("A", 134336176);
    graph.create_handle("G", 134336262);
    graph.create_handle("G", 134336173);
    graph.create_handle("T", 134336341);
    graph.create_handle("G", 134336338);
    graph.create_handle("G", 134336174);
    graph.create_handle("C", 134336203);
    graph.create_handle("A", 134336209);
    graph.create_handle("CG", 134336270);
    graph.create_handle("C", 134336251);
    graph.create_handle("T", 134336331);
    graph.create_handle("C", 134336315);
    graph.create_handle("A", 134336172);
    graph.create_handle("A", 134336183);
    graph.create_handle("T", 134336349);
    graph.create_handle("A", 134336243);
    graph.create_handle("C", 134336200);
    graph.create_handle("G", 134336355);
    graph.create_handle("T", 134336234);
    graph.create_handle("T", 134336286);
    graph.create_handle("C", 134336305);
    graph.create_handle("C", 134336283);
    graph.create_handle("C", 134336288);
    graph.create_handle("C", 134336313);
    graph.create_handle("G", 134336285);
    graph.create_handle("C", 134336264);
    graph.create_handle("A", 134336218);
    graph.create_handle("T", 134336301);
    graph.create_handle("G", 134336236);
    graph.create_handle("C", 134336225);
    graph.create_handle("A", 134336318);
    graph.create_handle("T", 134336261);
    graph.create_handle("TC", 134336304);
    graph.create_handle("C", 134336280);
    graph.create_handle("C", 134336164);
    graph.create_handle("G", 134336231);
    graph.create_handle("G", 134336226);
    graph.create_handle("T", 134336312);
    graph.create_handle("C", 134336247);
    graph.create_handle("A", 134336227);
    graph.create_handle("C", 134336235);
    graph.create_handle("A", 134336165);
    graph.create_handle("G", 134336167);
    graph.create_handle("T", 134336303);
    
    graph.create_edge(graph.get_handle(134336223, true), graph.get_handle(134336225, false));
    graph.create_edge(graph.get_handle(134336228, false), graph.get_handle(134336223, false));
    graph.create_edge(graph.get_handle(134336224, false), graph.get_handle(134336223, false));
    graph.create_edge(graph.get_handle(134336281, false), graph.get_handle(134336284, false));
    graph.create_edge(graph.get_handle(134336248, false), graph.get_handle(134336250, false));
    graph.create_edge(graph.get_handle(134336248, false), graph.get_handle(134336254, false));
    graph.create_edge(graph.get_handle(134336291, false), graph.get_handle(134336292, false));
    graph.create_edge(graph.get_handle(134336223, false), graph.get_handle(134336221, false));
    graph.create_edge(graph.get_handle(134336259, true), graph.get_handle(134336260, false));
    graph.create_edge(graph.get_handle(134336230, false), graph.get_handle(134336231, true));
    graph.create_edge(graph.get_handle(134336230, false), graph.get_handle(134336233, false));
    graph.create_edge(graph.get_handle(134336346, true), graph.get_handle(134336348, false));
    graph.create_edge(graph.get_handle(134336181, false), graph.get_handle(134336182, false));
    graph.create_edge(graph.get_handle(134336181, false), graph.get_handle(134336184, false));
    graph.create_edge(graph.get_handle(134336289, true), graph.get_handle(134336290, false));
    graph.create_edge(graph.get_handle(134336257, true), graph.get_handle(134336258, false));
    graph.create_edge(graph.get_handle(134336259, false), graph.get_handle(134336257, false));
    graph.create_edge(graph.get_handle(134336257, false), graph.get_handle(134336256, false));
    graph.create_edge(graph.get_handle(134336254, false), graph.get_handle(134336255, true));
    graph.create_edge(graph.get_handle(134336208, false), graph.get_handle(134336209, true));
    graph.create_edge(graph.get_handle(134336284, false), graph.get_handle(134336287, false));
    graph.create_edge(graph.get_handle(134336284, false), graph.get_handle(134336285, true));
    graph.create_edge(graph.get_handle(134336307, true), graph.get_handle(134336310, false));
    graph.create_edge(graph.get_handle(134336307, true), graph.get_handle(134336308, false));
    graph.create_edge(graph.get_handle(134336233, false), graph.get_handle(134336234, false));
    graph.create_edge(graph.get_handle(134336233, false), graph.get_handle(134336235, false));
    graph.create_edge(graph.get_handle(134336195, false), graph.get_handle(134336196, false));
    graph.create_edge(graph.get_handle(134336323, false), graph.get_handle(134336324, false));
    graph.create_edge(graph.get_handle(134336323, false), graph.get_handle(134336325, false));
    graph.create_edge(graph.get_handle(134336260, false), graph.get_handle(134336261, false));
    graph.create_edge(graph.get_handle(134336260, false), graph.get_handle(134336262, true));
    graph.create_edge(graph.get_handle(134336327, false), graph.get_handle(134336331, false));
    graph.create_edge(graph.get_handle(134336179, true), graph.get_handle(134336182, false));
    graph.create_edge(graph.get_handle(134336229, false), graph.get_handle(134336230, false));
    graph.create_edge(graph.get_handle(134336240, true), graph.get_handle(134336245, false));
    graph.create_edge(graph.get_handle(134336240, true), graph.get_handle(134336241, false));
    graph.create_edge(graph.get_handle(134336340, false), graph.get_handle(134336341, false));
    graph.create_edge(graph.get_handle(134336340, false), graph.get_handle(134336344, true));
    graph.create_edge(graph.get_handle(134336232, true), graph.get_handle(134336234, false));
    graph.create_edge(graph.get_handle(134336221, false), graph.get_handle(134336219, false));
    graph.create_edge(graph.get_handle(134336237, false), graph.get_handle(134336238, false));
    graph.create_edge(graph.get_handle(134336237, false), graph.get_handle(134336240, true));
    graph.create_edge(graph.get_handle(134336186, false), graph.get_handle(134336187, false));
    graph.create_edge(graph.get_handle(134336186, false), graph.get_handle(134336188, true));
    graph.create_edge(graph.get_handle(134336186, false), graph.get_handle(134336189, false));
    graph.create_edge(graph.get_handle(134336308, false), graph.get_handle(134336311, false));
    graph.create_edge(graph.get_handle(134336274, true), graph.get_handle(134336275, false));
    graph.create_edge(graph.get_handle(134336276, false), graph.get_handle(134336274, false));
    graph.create_edge(graph.get_handle(134336257, false), graph.get_handle(134336255, false));
    graph.create_edge(graph.get_handle(134336256, false), graph.get_handle(134336255, false));
    graph.create_edge(graph.get_handle(134336330, true), graph.get_handle(134336331, false));
    graph.create_edge(graph.get_handle(134336332, false), graph.get_handle(134336330, false));
    graph.create_edge(graph.get_handle(134336177, false), graph.get_handle(134336181, false));
    graph.create_edge(graph.get_handle(134336350, true), graph.get_handle(134336351, false));
    graph.create_edge(graph.get_handle(134336296, false), graph.get_handle(134336297, false));
    graph.create_edge(graph.get_handle(134336215, false), graph.get_handle(134336214, false));
    graph.create_edge(graph.get_handle(134336316, false), graph.get_handle(134336320, false));
    graph.create_edge(graph.get_handle(134336316, false), graph.get_handle(134336323, false));
    graph.create_edge(graph.get_handle(134336220, false), graph.get_handle(134336221, true));
    graph.create_edge(graph.get_handle(134336220, false), graph.get_handle(134336222, false));
    graph.create_edge(graph.get_handle(134336213, false), graph.get_handle(134336241, false));
    graph.create_edge(graph.get_handle(134336302, false), graph.get_handle(134336304, false));
    graph.create_edge(graph.get_handle(134336188, true), graph.get_handle(134336190, false));
    graph.create_edge(graph.get_handle(134336330, false), graph.get_handle(134336329, false));
    graph.create_edge(graph.get_handle(134336258, false), graph.get_handle(134336260, false));
    graph.create_edge(graph.get_handle(134336222, false), graph.get_handle(134336223, true));
    graph.create_edge(graph.get_handle(134336197, false), graph.get_handle(134336198, false));
    graph.create_edge(graph.get_handle(134336319, false), graph.get_handle(134336320, false));
    graph.create_edge(graph.get_handle(134336278, false), graph.get_handle(134336279, false));
    graph.create_edge(graph.get_handle(134336278, false), graph.get_handle(134336283, false));
    graph.create_edge(graph.get_handle(134336324, false), graph.get_handle(134336326, false));
    graph.create_edge(graph.get_handle(134336344, false), graph.get_handle(134336343, false));
    graph.create_edge(graph.get_handle(134336242, false), graph.get_handle(134336243, true));
    graph.create_edge(graph.get_handle(134336242, false), graph.get_handle(134336244, false));
    graph.create_edge(graph.get_handle(134336242, false), graph.get_handle(134336247, false));
    graph.create_edge(graph.get_handle(134336318, false), graph.get_handle(134336317, false));
    graph.create_edge(graph.get_handle(134336348, false), graph.get_handle(134336349, false));
    graph.create_edge(graph.get_handle(134336348, false), graph.get_handle(134336350, true));
    graph.create_edge(graph.get_handle(134336348, false), graph.get_handle(134336352, false));
    graph.create_edge(graph.get_handle(134336180, true), graph.get_handle(134336181, false));
    graph.create_edge(graph.get_handle(134336269, true), graph.get_handle(134336271, false));
    graph.create_edge(graph.get_handle(134336351, false), graph.get_handle(134336354, false));
    graph.create_edge(graph.get_handle(134336351, false), graph.get_handle(134336355, false));
    graph.create_edge(graph.get_handle(134336351, false), graph.get_handle(134336357, false));
    graph.create_edge(graph.get_handle(134336293, false), graph.get_handle(134336298, false));
    graph.create_edge(graph.get_handle(134336293, false), graph.get_handle(134336294, true));
    graph.create_edge(graph.get_handle(134336294, true), graph.get_handle(134336301, false));
    graph.create_edge(graph.get_handle(134336194, false), graph.get_handle(134336195, false));
    graph.create_edge(graph.get_handle(134336311, false), graph.get_handle(134336312, false));
    graph.create_edge(graph.get_handle(134336311, false), graph.get_handle(134336313, false));
    graph.create_edge(graph.get_handle(134336314, false), graph.get_handle(134336315, false));
    graph.create_edge(graph.get_handle(134336314, false), graph.get_handle(134336318, true));
    graph.create_edge(graph.get_handle(134336314, false), graph.get_handle(134336321, true));
    graph.create_edge(graph.get_handle(134336353, false), graph.get_handle(134336357, false));
    graph.create_edge(graph.get_handle(134336309, false), graph.get_handle(134336311, false));
    graph.create_edge(graph.get_handle(134336357, false), graph.get_handle(134336358, false));
    graph.create_edge(graph.get_handle(134336357, false), graph.get_handle(134336359, true));
    graph.create_edge(graph.get_handle(134336344, true), graph.get_handle(134336345, false));
    graph.create_edge(graph.get_handle(134336346, false), graph.get_handle(134336344, false));
    graph.create_edge(graph.get_handle(134336344, true), graph.get_handle(134336347, false));
    graph.create_edge(graph.get_handle(134336306, false), graph.get_handle(134336308, false));
    graph.create_edge(graph.get_handle(134336342, false), graph.get_handle(134336343, true));
    graph.create_edge(graph.get_handle(134336320, false), graph.get_handle(134336323, false));
    graph.create_edge(graph.get_handle(134336175, false), graph.get_handle(134336176, false));
    graph.create_edge(graph.get_handle(134336175, false), graph.get_handle(134336179, true));
    graph.create_edge(graph.get_handle(134336175, false), graph.get_handle(134336181, false));
    graph.create_edge(graph.get_handle(134336290, false), graph.get_handle(134336291, false));
    graph.create_edge(graph.get_handle(134336206, false), graph.get_handle(134336210, false));
    graph.create_edge(graph.get_handle(134336206, false), graph.get_handle(134336212, false));
    graph.create_edge(graph.get_handle(134336272, false), graph.get_handle(134336273, false));
    graph.create_edge(graph.get_handle(134336272, false), graph.get_handle(134336274, true));
    graph.create_edge(graph.get_handle(134336322, false), graph.get_handle(134336321, false));
    graph.create_edge(graph.get_handle(134336253, false), graph.get_handle(134336258, false));
    graph.create_edge(graph.get_handle(134336250, false), graph.get_handle(134336251, false));
    graph.create_edge(graph.get_handle(134336250, false), graph.get_handle(134336252, false));
    graph.create_edge(graph.get_handle(134336250, false), graph.get_handle(134336255, true));
    graph.create_edge(graph.get_handle(134336250, false), graph.get_handle(134336261, false));
    graph.create_edge(graph.get_handle(134336271, false), graph.get_handle(134336272, false));
    graph.create_edge(graph.get_handle(134336246, false), graph.get_handle(134336247, false));
    graph.create_edge(graph.get_handle(134336241, false), graph.get_handle(134336242, false));
    graph.create_edge(graph.get_handle(134336241, false), graph.get_handle(134336246, false));
    graph.create_edge(graph.get_handle(134336273, false), graph.get_handle(134336275, false));
    graph.create_edge(graph.get_handle(134336273, false), graph.get_handle(134336277, false));
    graph.create_edge(graph.get_handle(134336322, true), graph.get_handle(134336323, false));
    graph.create_edge(graph.get_handle(134336226, false), graph.get_handle(134336224, false));
    graph.create_edge(graph.get_handle(134336204, false), graph.get_handle(134336205, false));
    graph.create_edge(graph.get_handle(134336345, false), graph.get_handle(134336348, false));
    graph.create_edge(graph.get_handle(134336185, false), graph.get_handle(134336186, false));
    graph.create_edge(graph.get_handle(134336265, true), graph.get_handle(134336267, false));
    graph.create_edge(graph.get_handle(134336270, false), graph.get_handle(134336265, false));
    graph.create_edge(graph.get_handle(134336182, false), graph.get_handle(134336183, true));
    graph.create_edge(graph.get_handle(134336182, false), graph.get_handle(134336185, false));
    graph.create_edge(graph.get_handle(134336202, false), graph.get_handle(134336207, true));
    graph.create_edge(graph.get_handle(134336202, false), graph.get_handle(134336203, false));
    graph.create_edge(graph.get_handle(134336202, false), graph.get_handle(134336204, false));
    graph.create_edge(graph.get_handle(134336192, false), graph.get_handle(134336193, false));
    graph.create_edge(graph.get_handle(134336249, false), graph.get_handle(134336250, false));
    graph.create_edge(graph.get_handle(134336198, false), graph.get_handle(134336199, true));
    graph.create_edge(graph.get_handle(134336189, false), graph.get_handle(134336190, false));
    graph.create_edge(graph.get_handle(134336209, false), graph.get_handle(134336207, false));
    graph.create_edge(graph.get_handle(134336292, false), graph.get_handle(134336293, false));
    graph.create_edge(graph.get_handle(134336292, false), graph.get_handle(134336297, false));
    graph.create_edge(graph.get_handle(134336292, false), graph.get_handle(134336295, true));
    graph.create_edge(graph.get_handle(134336199, true), graph.get_handle(134336200, false));
    graph.create_edge(graph.get_handle(134336336, false), graph.get_handle(134336338, false));
    graph.create_edge(graph.get_handle(134336193, false), graph.get_handle(134336194, false));
    graph.create_edge(graph.get_handle(134336212, false), graph.get_handle(134336248, false));
    graph.create_edge(graph.get_handle(134336212, false), graph.get_handle(134336213, false));
    graph.create_edge(graph.get_handle(134336287, false), graph.get_handle(134336288, false));
    graph.create_edge(graph.get_handle(134336287, false), graph.get_handle(134336289, true));
    graph.create_edge(graph.get_handle(134336245, false), graph.get_handle(134336247, false));
    graph.create_edge(graph.get_handle(134336354, false), graph.get_handle(134336358, false));
    graph.create_edge(graph.get_handle(134336339, false), graph.get_handle(134336340, false));
    graph.create_edge(graph.get_handle(134336339, false), graph.get_handle(134336343, true));
    graph.create_edge(graph.get_handle(134336218, false), graph.get_handle(134336216, false));
    graph.create_edge(graph.get_handle(134336295, true), graph.get_handle(134336296, false));
    graph.create_edge(graph.get_handle(134336216, false), graph.get_handle(134336215, false));
    graph.create_edge(graph.get_handle(134336215, true), graph.get_handle(134336217, false));
    graph.create_edge(graph.get_handle(134336337, false), graph.get_handle(134336338, false));
    graph.create_edge(graph.get_handle(134336201, false), graph.get_handle(134336202, false));
    graph.create_edge(graph.get_handle(134336347, false), graph.get_handle(134336348, false));
    graph.create_edge(graph.get_handle(134336298, false), graph.get_handle(134336302, false));
    graph.create_edge(graph.get_handle(134336298, false), graph.get_handle(134336299, false));
    graph.create_edge(graph.get_handle(134336298, false), graph.get_handle(134336301, false));
    graph.create_edge(graph.get_handle(134336190, false), graph.get_handle(134336191, false));
    graph.create_edge(graph.get_handle(134336190, false), graph.get_handle(134336199, true));
    graph.create_edge(graph.get_handle(134336310, false), graph.get_handle(134336311, false));
    graph.create_edge(graph.get_handle(134336332, true), graph.get_handle(134336333, false));
    graph.create_edge(graph.get_handle(134336334, false), graph.get_handle(134336332, false));
    graph.create_edge(graph.get_handle(134336263, false), graph.get_handle(134336264, false));
    graph.create_edge(graph.get_handle(134336263, false), graph.get_handle(134336265, true));
    graph.create_edge(graph.get_handle(134336263, false), graph.get_handle(134336266, false));
    graph.create_edge(graph.get_handle(134336266, false), graph.get_handle(134336267, false));
    graph.create_edge(graph.get_handle(134336266, false), graph.get_handle(134336291, false));
    graph.create_edge(graph.get_handle(134336266, false), graph.get_handle(134336290, false));
    graph.create_edge(graph.get_handle(134336187, false), graph.get_handle(134336190, false));
    graph.create_edge(graph.get_handle(134336217, false), graph.get_handle(134336218, true));
    graph.create_edge(graph.get_handle(134336267, false), graph.get_handle(134336268, false));
    graph.create_edge(graph.get_handle(134336267, false), graph.get_handle(134336269, true));
    graph.create_edge(graph.get_handle(134336178, false), graph.get_handle(134336179, true));
    graph.create_edge(graph.get_handle(134336240, false), graph.get_handle(134336239, false));
    graph.create_edge(graph.get_handle(134336300, false), graph.get_handle(134336304, false));
    graph.create_edge(graph.get_handle(134336328, false), graph.get_handle(134336332, true));
    graph.create_edge(graph.get_handle(134336333, false), graph.get_handle(134336335, false));
    graph.create_edge(graph.get_handle(134336333, false), graph.get_handle(134336336, false));
    graph.create_edge(graph.get_handle(134336184, false), graph.get_handle(134336185, false));
    graph.create_edge(graph.get_handle(134336279, false), graph.get_handle(134336283, false));
    graph.create_edge(graph.get_handle(134336279, false), graph.get_handle(134336280, true));
    graph.create_edge(graph.get_handle(134336297, false), graph.get_handle(134336298, false));
    graph.create_edge(graph.get_handle(134336275, false), graph.get_handle(134336278, false));
    graph.create_edge(graph.get_handle(134336275, false), graph.get_handle(134336281, false));
    graph.create_edge(graph.get_handle(134336275, false), graph.get_handle(134336282, true));
    graph.create_edge(graph.get_handle(134336275, false), graph.get_handle(134336283, false));
    graph.create_edge(graph.get_handle(134336334, true), graph.get_handle(134336335, false));
    graph.create_edge(graph.get_handle(134336334, true), graph.get_handle(134336337, false));
    graph.create_edge(graph.get_handle(134336228, true), graph.get_handle(134336229, false));
    graph.create_edge(graph.get_handle(134336205, false), graph.get_handle(134336206, false));
    graph.create_edge(graph.get_handle(134336205, false), graph.get_handle(134336207, true));
    graph.create_edge(graph.get_handle(134336205, false), graph.get_handle(134336208, false));
    graph.create_edge(graph.get_handle(134336252, false), graph.get_handle(134336298, false));
    graph.create_edge(graph.get_handle(134336252, false), graph.get_handle(134336253, false));
    graph.create_edge(graph.get_handle(134336252, false), graph.get_handle(134336261, false));
    graph.create_edge(graph.get_handle(134336252, false), graph.get_handle(134336272, false));
    graph.create_edge(graph.get_handle(134336282, true), graph.get_handle(134336283, false));
    graph.create_edge(graph.get_handle(134336335, false), graph.get_handle(134336339, false));
    graph.create_edge(graph.get_handle(134336335, false), graph.get_handle(134336342, false));
    graph.create_edge(graph.get_handle(134336326, false), graph.get_handle(134336327, false));
    graph.create_edge(graph.get_handle(134336326, false), graph.get_handle(134336328, false));
    graph.create_edge(graph.get_handle(134336326, false), graph.get_handle(134336330, true));
    graph.create_edge(graph.get_handle(134336196, false), graph.get_handle(134336197, false));
    graph.create_edge(graph.get_handle(134336268, false), graph.get_handle(134336271, false));
    graph.create_edge(graph.get_handle(134336325, false), graph.get_handle(134336326, false));
    graph.create_edge(graph.get_handle(134336325, false), graph.get_handle(134336329, true));
    graph.create_edge(graph.get_handle(134336276, true), graph.get_handle(134336278, false));
    graph.create_edge(graph.get_handle(134336191, false), graph.get_handle(134336200, false));
    graph.create_edge(graph.get_handle(134336191, false), graph.get_handle(134336192, false));
    graph.create_edge(graph.get_handle(134336215, false), graph.get_handle(134336211, false));
    graph.create_edge(graph.get_handle(134336299, false), graph.get_handle(134336300, false));
    graph.create_edge(graph.get_handle(134336299, false), graph.get_handle(134336302, false));
    graph.create_edge(graph.get_handle(134336238, false), graph.get_handle(134336241, false));
    graph.create_edge(graph.get_handle(134336277, false), graph.get_handle(134336278, false));
    graph.create_edge(graph.get_handle(134336352, false), graph.get_handle(134336353, false));
    graph.create_edge(graph.get_handle(134336352, false), graph.get_handle(134336356, true));
    graph.create_edge(graph.get_handle(134336210, false), graph.get_handle(134336248, false));
    graph.create_edge(graph.get_handle(134336210, false), graph.get_handle(134336249, false));
    graph.create_edge(graph.get_handle(134336210, false), graph.get_handle(134336211, true));
    graph.create_edge(graph.get_handle(134336210, false), graph.get_handle(134336214, true));
    graph.create_edge(graph.get_handle(134336210, false), graph.get_handle(134336218, true));
    graph.create_edge(graph.get_handle(134336210, false), graph.get_handle(134336229, false));
    graph.create_edge(graph.get_handle(134336210, false), graph.get_handle(134336234, false));
    graph.create_edge(graph.get_handle(134336244, false), graph.get_handle(134336248, false));
    graph.create_edge(graph.get_handle(134336356, true), graph.get_handle(134336357, false));
    graph.create_edge(graph.get_handle(134336176, false), graph.get_handle(134336185, false));
    graph.create_edge(graph.get_handle(134336176, false), graph.get_handle(134336182, false));
    graph.create_edge(graph.get_handle(134336262, true), graph.get_handle(134336263, false));
    graph.create_edge(graph.get_handle(134336173, true), graph.get_handle(134336175, false));
    graph.create_edge(graph.get_handle(134336341, false), graph.get_handle(134336345, false));
    graph.create_edge(graph.get_handle(134336341, false), graph.get_handle(134336348, false));
    graph.create_edge(graph.get_handle(134336338, false), graph.get_handle(134336339, false));
    graph.create_edge(graph.get_handle(134336174, false), graph.get_handle(134336175, false));
    graph.create_edge(graph.get_handle(134336174, false), graph.get_handle(134336177, false));
    graph.create_edge(graph.get_handle(134336203, false), graph.get_handle(134336205, false));
    graph.create_edge(graph.get_handle(134336209, true), graph.get_handle(134336210, false));
    graph.create_edge(graph.get_handle(134336209, true), graph.get_handle(134336212, false));
    graph.create_edge(graph.get_handle(134336270, true), graph.get_handle(134336271, false));
    graph.create_edge(graph.get_handle(134336251, false), graph.get_handle(134336298, false));
    graph.create_edge(graph.get_handle(134336251, false), graph.get_handle(134336253, false));
    graph.create_edge(graph.get_handle(134336251, false), graph.get_handle(134336257, true));
    graph.create_edge(graph.get_handle(134336251, false), graph.get_handle(134336272, false));
    graph.create_edge(graph.get_handle(134336331, false), graph.get_handle(134336333, false));
    graph.create_edge(graph.get_handle(134336331, false), graph.get_handle(134336334, true));
    graph.create_edge(graph.get_handle(134336315, false), graph.get_handle(134336316, false));
    graph.create_edge(graph.get_handle(134336315, false), graph.get_handle(134336317, true));
    graph.create_edge(graph.get_handle(134336315, false), graph.get_handle(134336319, false));
    graph.create_edge(graph.get_handle(134336172, true), graph.get_handle(134336175, false));
    graph.create_edge(graph.get_handle(134336172, true), graph.get_handle(134336178, false));
    graph.create_edge(graph.get_handle(134336183, true), graph.get_handle(134336186, false));
    graph.create_edge(graph.get_handle(134336349, false), graph.get_handle(134336351, false));
    graph.create_edge(graph.get_handle(134336243, true), graph.get_handle(134336248, false));
    graph.create_edge(graph.get_handle(134336200, false), graph.get_handle(134336202, false));
    graph.create_edge(graph.get_handle(134336200, false), graph.get_handle(134336201, false));
    graph.create_edge(graph.get_handle(134336355, false), graph.get_handle(134336359, true));
    graph.create_edge(graph.get_handle(134336234, false), graph.get_handle(134336236, true));
    graph.create_edge(graph.get_handle(134336234, false), graph.get_handle(134336237, false));
    graph.create_edge(graph.get_handle(134336286, false), graph.get_handle(134336287, false));
    graph.create_edge(graph.get_handle(134336305, false), graph.get_handle(134336309, false));
    graph.create_edge(graph.get_handle(134336305, false), graph.get_handle(134336306, false));
    graph.create_edge(graph.get_handle(134336283, false), graph.get_handle(134336284, false));
    graph.create_edge(graph.get_handle(134336283, false), graph.get_handle(134336286, false));
    graph.create_edge(graph.get_handle(134336288, false), graph.get_handle(134336292, false));
    graph.create_edge(graph.get_handle(134336288, false), graph.get_handle(134336296, false));
    graph.create_edge(graph.get_handle(134336288, false), graph.get_handle(134336290, false));
    graph.create_edge(graph.get_handle(134336313, false), graph.get_handle(134336314, false));
    graph.create_edge(graph.get_handle(134336289, false), graph.get_handle(134336285, false));
    graph.create_edge(graph.get_handle(134336264, false), graph.get_handle(134336267, false));
    graph.create_edge(graph.get_handle(134336219, false), graph.get_handle(134336218, false));
    graph.create_edge(graph.get_handle(134336218, true), graph.get_handle(134336220, false));
    graph.create_edge(graph.get_handle(134336301, false), graph.get_handle(134336302, false));
    graph.create_edge(graph.get_handle(134336301, false), graph.get_handle(134336303, false));
    graph.create_edge(graph.get_handle(134336236, true), graph.get_handle(134336237, false));
    graph.create_edge(graph.get_handle(134336239, false), graph.get_handle(134336236, false));
    graph.create_edge(graph.get_handle(134336225, false), graph.get_handle(134336229, false));
    graph.create_edge(graph.get_handle(134336225, false), graph.get_handle(134336226, true));
    graph.create_edge(graph.get_handle(134336322, false), graph.get_handle(134336318, false));
    graph.create_edge(graph.get_handle(134336261, false), graph.get_handle(134336263, false));
    graph.create_edge(graph.get_handle(134336261, false), graph.get_handle(134336267, false));
    graph.create_edge(graph.get_handle(134336304, false), graph.get_handle(134336307, true));
    graph.create_edge(graph.get_handle(134336304, false), graph.get_handle(134336305, false));
    graph.create_edge(graph.get_handle(134336280, true), graph.get_handle(134336284, false));
    graph.create_edge(graph.get_handle(134336164, false), graph.get_handle(134336175, false));
    graph.create_edge(graph.get_handle(134336164, false), graph.get_handle(134336180, true));
    graph.create_edge(graph.get_handle(134336164, false), graph.get_handle(134336167, true));
    graph.create_edge(graph.get_handle(134336231, true), graph.get_handle(134336234, false));
    graph.create_edge(graph.get_handle(134336232, false), graph.get_handle(134336231, false));
    graph.create_edge(graph.get_handle(134336226, true), graph.get_handle(134336230, false));
    graph.create_edge(graph.get_handle(134336227, false), graph.get_handle(134336226, false));
    graph.create_edge(graph.get_handle(134336312, false), graph.get_handle(134336314, false));
    graph.create_edge(graph.get_handle(134336247, false), graph.get_handle(134336248, false));
    graph.create_edge(graph.get_handle(134336247, false), graph.get_handle(134336249, false));
    graph.create_edge(graph.get_handle(134336231, false), graph.get_handle(134336227, false));
    graph.create_edge(graph.get_handle(134336235, false), graph.get_handle(134336236, true));
    graph.create_edge(graph.get_handle(134336165, true), graph.get_handle(134336175, false));
    graph.create_edge(graph.get_handle(134336167, false), graph.get_handle(134336165, false));
    graph.create_edge(graph.get_handle(134336167, true), graph.get_handle(134336175, false));
    graph.create_edge(graph.get_handle(134336303, false), graph.get_handle(134336304, false));
    
    
    IncrementalSubgraph subgraph(graph, make_pos_t(134336263, true, 1), true, 34);
    
    // goal is just to not crash
    while (subgraph.is_extendable()) {
        subgraph.extend();
    }
}

}
}
        
