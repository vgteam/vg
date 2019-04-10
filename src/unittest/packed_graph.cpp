/// \file packed_graph.cpp
///  
/// Unit tests for the PackedGraph class.
///

#include <iostream>
#include <sstream>

#include "../json2pb.h"
#include "../packed_graph.hpp"
#include "random_graph.hpp"

#include "algorithms/are_equivalent.hpp"

#include "catch.hpp"


namespace vg {
namespace unittest {
using namespace std;

    TEST_CASE("PackedGraph's reallocations do not change topology", "[packed][handle]") {
        
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
        
        handle_t h1 = graph.create_handle("ATGTAG");
        handle_t h2 = graph.create_handle("ACCCC");
        handle_t h3 = graph.create_handle("C");
        handle_t h4 = graph.create_handle("ATT");
        handle_t h5 = graph.create_handle("GGCA");
        
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
        
        SECTION("Defragmentation during dynamic deletes does not change topology") {
            
            // delete enough nodes/edges/path memberships to trigger defrag
            
            graph.destroy_path(p0);
            graph.destroy_path(p2);
            graph.destroy_handle(h2);
            graph.destroy_handle(h4);
            
            REQUIRE(graph.get_sequence(h1) == "ATGTAG");
            REQUIRE(graph.get_sequence(h3) == "C");
            REQUIRE(graph.get_sequence(h5) == "GGCA");
            
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
        
        SECTION("Compactification does not change topology") {
            
            // delete some things, but not enough to trigger defragmentation
            graph.destroy_path(p2);
            graph.destroy_handle(h2);
            
            // reallocate and compress down to the smaller size
            graph.compactify();
            
            REQUIRE(graph.get_sequence(h1) == "ATGTAG");
            REQUIRE(graph.get_sequence(h3) == "C");
            REQUIRE(graph.get_sequence(h4) == "ATT");
            REQUIRE(graph.get_sequence(h5) == "GGCA");
            
            int count = 0;
            bool found1 = false, found2 = false;
            graph.follow_edges(h1, false, [&](const handle_t& h) {
                if (h == h3) {
                    found1 = true;
                }
                count++;
            });
            REQUIRE(found1);
            REQUIRE(count == 1);
            
            count = 0;
            found1 = false, found2 = false;
            graph.follow_edges(h1, true, [&](const handle_t& h) {
                count++;
            });
            REQUIRE(count == 0);
            
            count = 0;
            found1 = false, found2 = false;
            graph.follow_edges(h3, false, [&](const handle_t& h) {
                if (h == h4) {
                    found1 = true;
                }
                if (h == h5) {
                    found2 = true;
                }
                count++;
            });
            REQUIRE(found1);
            REQUIRE(found2);
            REQUIRE(count == 2);
            
            count = 0;
            found1 = false, found2 = false;
            graph.follow_edges(h3, true, [&](const handle_t& h) {
                if (h == h1) {
                    found1 = true;
                }
                count++;
            });
            REQUIRE(found1);
            REQUIRE(count == 1);
            
            count = 0;
            found1 = false, found2 = false;
            graph.follow_edges(h4, false, [&](const handle_t& h) {
                if (h == h5) {
                    found1 = true;
                }
                count++;
            });
            REQUIRE(found1);
            REQUIRE(count == 1);
            
            count = 0;
            found1 = false, found2 = false;
            graph.follow_edges(h4, true, [&](const handle_t& h) {
                if (h == h3) {
                    found1 = true;
                }
                count++;
            });
            REQUIRE(found1);
            REQUIRE(count == 1);
            
            count = 0;
            found1 = false, found2 = false;
            graph.follow_edges(h5, false, [&](const handle_t& h) {
                count++;
            });
            REQUIRE(count == 0);
            
            count = 0;
            found1 = false, found2 = false;
            graph.follow_edges(h5, true, [&](const handle_t& h) {
                if (h == h3) {
                    found1 = true;
                }
                else if (h == h4) {
                    found2 = true;
                }
                count++;
            });
            REQUIRE(found1);
            REQUIRE(found2);
            REQUIRE(count == 2);
            
            check_flips(p0, {h3, h4, h5});
            check_flips(p1, {h1, h3, h5});
        }
    }
    
    TEST_CASE("PackedGraph serialization works on randomized graphs", "[packed]") {
        
        int num_graphs = 100;
        int seq_length = 200;
        int num_variants = 30;
        int long_var_length = 10;
        
        for (int i = 0; i < num_graphs; i++) {
            
            PackedGraph graph;
            random_graph(seq_length, long_var_length, num_variants, &graph);
            
            stringstream strm;
            graph.serialize(strm);
            strm.seekg(0);
            PackedGraph loaded(strm);
            
            REQUIRE(algorithms::are_equivalent_with_paths(&graph, &loaded));
        }
    }
}
}
        
