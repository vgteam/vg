/// \file packed_graph.cpp
///  
/// Unit tests for the PackedGraph class.
///

#include <iostream>

#include "../json2pb.h"
#include "../packed_graph.hpp"
#include "catch.hpp"

namespace vg {
namespace unittest {
using namespace std;

    TEST_CASE("PackedGraph's defragmentation preserves topology", "[packed][handle]") {
        
        PackedGraph graph;
        
        auto check_path = [&](const path_handle_t& p, const vector<handle_t>& occs) {
            
            occurrence_handle_t occ;
            for (int i = 0; i < occs.size(); i++){
                if (i == 0) {
                    occ = graph.get_first_occurrence(p);
                }
                
                REQUIRE(graph.get_path_handle_of_occurrence(occ) == p);
                REQUIRE(graph.get_occurrence(occ) == occs[i]);
                REQUIRE(graph.has_previous_occurrence(occ) == (i > 0));
                REQUIRE(graph.has_next_occurrence(occ) == (i < occs.size() - 1));
                
                if (i != occs.size() - 1) {
                    occ = graph.get_next_occurrence(occ);
                }
            }
            
            for (int i = occs.size() - 1; i >= 0; i--){
                if (i == occs.size() - 1) {
                    occ = graph.get_last_occurrence(p);
                }
                
                REQUIRE(graph.get_path_handle_of_occurrence(occ) == p);
                REQUIRE(graph.get_occurrence(occ) == occs[i]);
                REQUIRE(graph.has_previous_occurrence(occ) == (i > 0));
                REQUIRE(graph.has_next_occurrence(occ) == (i < occs.size() - 1));
                
                if (i != 0) {
                    occ = graph.get_previous_occurrence(occ);
                }
            }
        };
        
        auto check_flips = [&](const path_handle_t& p, const vector<handle_t>& occs) {
            auto flipped = occs;
            for (size_t i = 0; i < occs.size(); i++) {
                
                graph.apply_orientation(graph.flip(graph.forward(flipped[i])));
                flipped[i] = graph.flip(flipped[i]);
                check_path(p, flipped);
                
                graph.apply_orientation(graph.flip(graph.forward(flipped[i])));
                flipped[i] = graph.flip(flipped[i]);
                check_path(p, flipped);
            }
        };
        
        handle_t h1 = graph.create_handle("A");
        handle_t h2 = graph.create_handle("A");
        handle_t h3 = graph.create_handle("C");
        handle_t h4 = graph.create_handle("A");
        handle_t h5 = graph.create_handle("G");
        
        graph.create_edge(h1, h2);
        graph.create_edge(h1, h3);
        graph.create_edge(h2, h3);
        graph.create_edge(h3, h5);
        graph.create_edge(h3, h4);
        graph.create_edge(h4, h5);
        
        path_handle_t p0 = graph.create_path_handle("0");
        path_handle_t p1 = graph.create_path_handle("1");
        path_handle_t p2 = graph.create_path_handle("2");
        
        
        graph.append_occurrence(p0, h3);
        graph.append_occurrence(p0, h4);
        graph.append_occurrence(p0, h5);
        
        graph.append_occurrence(p1, h1);
        graph.append_occurrence(p1, h3);
        graph.append_occurrence(p1, h5);
        
        graph.append_occurrence(p2, h1);
        graph.append_occurrence(p2, h2);
        graph.append_occurrence(p2, h3);
        graph.append_occurrence(p2, h4);
        graph.append_occurrence(p2, h5);
        
        
        // delete enough nodes/edges/path memberships to trigger defrag
        
        graph.destroy_path(p0);
        graph.destroy_path(p2);
        graph.destroy_handle(h2);
        graph.destroy_handle(h4);
        
        REQUIRE(graph.get_sequence(h1) == "A");
        REQUIRE(graph.get_sequence(h3) == "C");
        REQUIRE(graph.get_sequence(h5) == "G");
        
        bool found = false;
        graph.follow_edges(h1, false, [&](const handle_t& next) {
            if (next == h3) {
                found = true;
            }
            else {
                REQUIRE(false);
            }
            return true;
        });
        REQUIRE(found);
        
        found = false;
        graph.follow_edges(h3, false, [&](const handle_t& next) {
            if (next == h5) {
                found = true;
            }
            else {
                REQUIRE(false);
            }
            return true;
        });
        REQUIRE(found);
        
        check_flips(p1, {h1, h3, h5});
    }
}
}
        
