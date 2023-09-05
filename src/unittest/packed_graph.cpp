/// \file packed_graph.cpp
///  
/// Unit tests for the PackedGraph class.
///

#include <iostream>
#include <sstream>
#include <set>

#include "vg/io/json2pb.h"
#include "random_graph.hpp"
#include "randomness.hpp"

#include "bdsg/packed_graph.hpp"

#include "catch.hpp"


namespace vg {
namespace unittest {
using namespace std;

    TEST_CASE("PackedGraph's reallocations do not change topology", "[packed][handle]") {
        
        bdsg::PackedGraph graph;
        
        auto check_path = [&](const path_handle_t& p, const vector<handle_t>& steps) {
            
            REQUIRE(graph.get_step_count(p) == steps.size());
            
            step_handle_t step = graph.path_begin(p);
            for (int i = 0; i < steps.size(); i++) {
                
                REQUIRE(graph.get_path_handle_of_step(step) == p);
                REQUIRE(graph.get_handle_of_step(step) == steps[i]);
                
                if (graph.get_is_circular(p)) {
                    REQUIRE(graph.has_next_step(step));
                    REQUIRE(graph.has_previous_step(step));
                }
                else {
                    REQUIRE(graph.has_next_step(step) == i + 1 < steps.size());
                    REQUIRE(graph.has_previous_step(step) == i > 0);
                }
                
                step = graph.get_next_step(step);
            }
            
            if (graph.get_is_circular(p) && !graph.is_empty(p)) {
                REQUIRE(step == graph.path_begin(p));
            }
            else {
                REQUIRE(step == graph.path_end(p));
            }
            
            step = graph.path_back(p);
            
            for (int i = steps.size() - 1; i >= 0; i--) {
                
                REQUIRE(graph.get_path_handle_of_step(step) == p);
                REQUIRE(graph.get_handle_of_step(step) == steps[i]);
                
                if (graph.get_is_circular(p)) {
                    REQUIRE(graph.has_next_step(step));
                    REQUIRE(graph.has_previous_step(step));
                }
                else {
                    REQUIRE(graph.has_next_step(step) == i + 1 < steps.size());
                    REQUIRE(graph.has_previous_step(step) == i > 0);
                }
                
                step = graph.get_previous_step(step);
            }
            
            if (graph.get_is_circular(p) && !graph.is_empty(p)) {
                REQUIRE(step == graph.path_back(p));
            }
            else {
                REQUIRE(step == graph.path_front_end(p));
            }
        };
        
        auto check_flips = [&](const path_handle_t& p, const vector<handle_t>& steps) {
            auto flipped = steps;
            for (size_t i = 0; i < steps.size(); i++) {
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

        
        graph.append_step(p0, h3);
        graph.append_step(p0, h4);
        graph.append_step(p0, h5);
        
        graph.append_step(p1, h1);
        graph.append_step(p1, h3);
        graph.append_step(p1, h5);
        
        graph.append_step(p2, h1);
        graph.append_step(p2, h2);
        graph.append_step(p2, h3);
        graph.append_step(p2, h4);
        graph.append_step(p2, h5);
        
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
        
        SECTION("Tighten vector allocations does not change topology") {
            
            // delete some things, but not enough to trigger defragmentation
            graph.destroy_path(p2);
            graph.destroy_handle(h2);
            // reallocate and compress down to the smaller size
            graph.optimize(false);
            
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
        
        SECTION("Optimizing with ID reassignment does not change topology") {
            
            // delete some things, but not enough to trigger defragmentation
            graph.destroy_path(p2);
            graph.destroy_handle(h2);
            // reallocate and compress down to the smaller size, reassigning IDs
            graph.optimize(true);
            
            set<id_t> seen_ids;
            
            int count = 0;
            bool found1 = false, found2 = false, found3 = false, found4 = false;
            graph.for_each_handle([&](const handle_t& handle) {
                if (graph.get_sequence(handle) == "ATGTAG") {
                    h1 = handle;
                    found1 = true;
                }
                else if (graph.get_sequence(handle) == "C") {
                    h3 = handle;
                    found2 = true;
                }
                else if (graph.get_sequence(handle) == "ATT") {
                    h4 = handle;
                    found3 = true;
                }
                else if (graph.get_sequence(handle) == "GGCA") {
                    h5 = handle;
                    found4 = true;
                }
                else {
                    REQUIRE(false);
                }
                count++;
                
                seen_ids.insert(graph.get_id(handle));
                
                REQUIRE(graph.get_id(handle) >= 1);
                REQUIRE(graph.get_id(handle) <= 4);
            });
            
            REQUIRE(found1);
            REQUIRE(found2);
            REQUIRE(found3);
            REQUIRE(found4);
            REQUIRE(count == 4);
            REQUIRE(seen_ids.size() == 4);
            
            count = 0;
            found1 = found2 = found3 = found4 = false;
            
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
        
        default_random_engine gen(test_seed_source());
        uniform_int_distribution<int> circ_distr(0, 1);
        
        for (int i = 0; i < num_graphs; i++) {
            
            bdsg::PackedGraph graph;
            random_graph(seq_length, long_var_length, num_variants, &graph);
            
            graph.for_each_path_handle([&](const path_handle_t& path) {
                graph.set_circularity(path, circ_distr(gen));
            });
            
            stringstream strm;
            graph.serialize(strm);
            strm.seekg(0);
            bdsg::PackedGraph loaded;
            loaded.deserialize(strm);
            
            REQUIRE(handlealgs::are_equivalent_with_paths(&graph, &loaded));
        }
    }
}
}
        
