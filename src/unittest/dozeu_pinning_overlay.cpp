/// \file dozeru_pinning_overlay.cpp
///  
/// unit tests for the overlay that sockets handle graphs that may have
/// nodes with no sequence into the XdropAligner for pinning
///

#include <iostream>
#include "../dozeu_pinning_overlay.hpp"
#include "catch.hpp"

#include "bdsg/hash_graph.hpp"

namespace vg {
namespace unittest {

    TEST_CASE("DozeuPinningOverlay produces expected topology for small test graph", "[xdrop][pinned][overlay]") {
        
        bdsg::HashGraph graph;
        
        handle_t h1 = graph.create_handle("", 1);
        handle_t h2 = graph.create_handle("CGGTG", 2);
        handle_t h3 = graph.create_handle("AGAA", 3);
        handle_t h4 = graph.create_handle("TTG", 4);
        
        graph.create_edge(h1, h2);
        graph.create_edge(h1, h3);
        graph.create_edge(h2, h3);
        graph.create_edge(h3, h4);
        
        DozeuPinningOverlay overlay(&graph, false);
        
        REQUIRE(overlay.performed_duplications());
        
        REQUIRE(overlay.get_node_count() == 4);
        size_t count = 0;
        bool found1 = false, found2 = false, found3 = false, found4 = false;
        handle_t o1, o2, o3, o4;
        overlay.for_each_handle([&](const handle_t& h) {
            if (overlay.get_sequence(h) == graph.get_sequence(h3)) {
                id_t node_id = overlay.get_id(h);
                if (node_id == 3) {
                    found1 = true;
                    o1 = h;
                }
                else if (node_id != 1 && node_id != 2 && node_id != 4) {
                    found2 = true;
                    o2 = h;
                }
            }
            else if (overlay.get_sequence(h) == graph.get_sequence(h2)) {
                REQUIRE(overlay.get_id(h) == 2);
                found3 = true;
                o3 = h;
            }
            else if (overlay.get_sequence(h) == graph.get_sequence(h4)) {
                REQUIRE(overlay.get_id(h) == 4);
                found4 = true;
                o4 = h;
            }
            else {
                REQUIRE(false);
            }
            
            ++count;
        });
        REQUIRE(count == 4);
        REQUIRE(found1);
        REQUIRE(found2);
        REQUIRE(found3);
        REQUIRE(found4);
        found1 = found2 = found3 = found4 = false;
        count = 0;
        
        REQUIRE(overlay.get_underlying_handle(o1) == h3);
        REQUIRE(overlay.get_underlying_handle(o2) == h3);
        REQUIRE(overlay.get_underlying_handle(o3) == h2);
        REQUIRE(overlay.get_underlying_handle(o4) == h4);
        REQUIRE(overlay.get_underlying_handle(overlay.flip(o1)) == graph.flip(h3));
        REQUIRE(overlay.get_underlying_handle(overlay.flip(o2)) == graph.flip(h3));
        REQUIRE(overlay.get_underlying_handle(overlay.flip(o3)) == graph.flip(h2));
        REQUIRE(overlay.get_underlying_handle(overlay.flip(o4)) == graph.flip(h4));
        
        for (handle_t h : {o1, o2, o3, o4}) {
            REQUIRE(overlay.has_node(overlay.get_id(h)));
        }
        
        unordered_set<id_t> ids;
        ids.insert(overlay.get_id(o1));
        ids.insert(overlay.get_id(o2));
        ids.insert(overlay.get_id(o3));
        ids.insert(overlay.get_id(o4));
        
        REQUIRE(ids.size() == 4);
        
        for (id_t i : ids) {
            REQUIRE(i <= overlay.max_node_id());
            REQUIRE(i >= overlay.min_node_id());
        }
        
        for (id_t i = 1; i < 50; ++i) {
            if (ids.count(i)) {
                REQUIRE(overlay.has_node(i));
            }
            else {
                REQUIRE(!overlay.has_node(i));
            }
        }
        
        for (handle_t h : {o1, o2, o3, o4}) {
            REQUIRE(h == overlay.flip(overlay.flip(h)));
        }
        
        // from o1
        
        overlay.follow_edges(o1, false, [&](const handle_t& h) {
            REQUIRE(h == o4);
            ++count;
        });
        REQUIRE(count == 1);
        count = 0;
        
        overlay.follow_edges(o1, true, [&](const handle_t& h) {
            REQUIRE(h == o3);
            ++count;
        });
        REQUIRE(count == 1);
        count = 0;
        
        overlay.follow_edges(overlay.flip(o1), false, [&](const handle_t& h) {
            REQUIRE(h == overlay.flip(o3));
            ++count;
        });
        REQUIRE(count == 1);
        count = 0;
        
        overlay.follow_edges(overlay.flip(o1), true, [&](const handle_t& h) {
            REQUIRE(h == overlay.flip(o4));
            ++count;
        });
        REQUIRE(count == 1);
        count = 0;
        
        // from o2
        
        overlay.follow_edges(o2, false, [&](const handle_t& h) {
            REQUIRE(h == o4);
            ++count;
        });
        REQUIRE(count == 1);
        count = 0;
        
        overlay.follow_edges(o2, true, [&](const handle_t& h) {
            ++count;
        });
        REQUIRE(count == 0);
        count = 0;
        
        overlay.follow_edges(overlay.flip(o2), false, [&](const handle_t& h) {
            ++count;
        });
        REQUIRE(count == 0);
        count = 0;
        
        overlay.follow_edges(overlay.flip(o2), true, [&](const handle_t& h) {
            REQUIRE(h == overlay.flip(o4));
            ++count;
        });
        REQUIRE(count == 1);
        count = 0;
        
        // from o3
        
        overlay.follow_edges(o3, false, [&](const handle_t& h) {
            REQUIRE(h == o1);
            ++count;
        });
        REQUIRE(count == 1);
        count = 0;
        
        overlay.follow_edges(o3, true, [&](const handle_t& h) {
            ++count;
        });
        REQUIRE(count == 0);
        count = 0;
        
        overlay.follow_edges(overlay.flip(o3), false, [&](const handle_t& h) {
            ++count;
        });
        REQUIRE(count == 0);
        count = 0;
        
        overlay.follow_edges(overlay.flip(o3), true, [&](const handle_t& h) {
            REQUIRE(h == overlay.flip(o1));
            ++count;
        });
        REQUIRE(count == 1);
        count = 0;
        
        // from o4
        
        overlay.follow_edges(o4, false, [&](const handle_t& h) {
            ++count;
        });
        REQUIRE(count == 0);
        count = 0;
        
        overlay.follow_edges(o4, true, [&](const handle_t& h) {
            if (h == o1) {
                found1 = true;
            }
            if (h == o2) {
                found2 = true;
            }
            ++count;
        });
        REQUIRE(count == 2);
        REQUIRE(found1);
        REQUIRE(found2);
        found1 = found2 = false;
        count = 0;
        
        overlay.follow_edges(overlay.flip(o4), false, [&](const handle_t& h) {
            if (h == overlay.flip(o1)) {
                found1 = true;
            }
            if (h == overlay.flip(o2)) {
                found2 = true;
            }
            ++count;
        });
        REQUIRE(count == 2);
        REQUIRE(found1);
        REQUIRE(found2);
        found1 = found2 = false;
        count = 0;
        
        overlay.follow_edges(overlay.flip(o4), true, [&](const handle_t& h) {
            ++count;
        });
        REQUIRE(count == 0);
        count = 0;
    }
}
}
