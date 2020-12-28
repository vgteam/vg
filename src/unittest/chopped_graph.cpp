#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <set>

#include "chopped_graph.hpp"
#include "catch.hpp"

#include "bdsg/hash_graph.hpp"
#include "bdsg/overlays/path_position_overlays.hpp"

namespace vg {
namespace unittest {

TEST_CASE("ChoppedGraph produces expected results for small test cases", "[chop][overlay]") {
    
    
    bdsg::HashGraph graph;
    
    handle_t h1 = graph.create_handle("TCGGCAG");
    handle_t h2 = graph.create_handle("GA");
    handle_t h3 = graph.create_handle("T");
    handle_t h4 = graph.create_handle("TTCC");
    
    graph.create_edge(h1, h2);
    graph.create_edge(h1, h3);
    graph.create_edge(h2, h4);
    graph.create_edge(h3, h4);
    graph.create_edge(h4, h1);
    
    path_handle_t p1 = graph.create_path_handle("p1");
    path_handle_t p2 = graph.create_path_handle("p2");
    graph.append_step(p1, h1);
    graph.append_step(p1, h2);
    graph.append_step(p1, h4);
    graph.append_step(p2, graph.flip(h4));
    graph.append_step(p2, graph.flip(h3));
    graph.append_step(p2, graph.flip(h1));
    graph.set_circularity(p2, true);
    
    bdsg::PositionOverlay pos_overlay(&graph);
    
    //ChoppedPathPositionGraph<2> chopped(pos_overlay);
    //vector<string> chopped_seqs{"TC", "GG", "CA", "G", "GA", "T", "TT", "CC"};
    
    ChoppedPathPositionGraph<3> chopped(pos_overlay);
    vector<string> chopped_seqs{"TCG", "GCA", "G", "GA", "T", "TTC", "C"};
    
    vector<handle_t> chopped_handles(chopped_seqs.size());
    vector<bool> chopped_found(chopped_seqs.size(), false);
    
    chopped.for_each_handle([&](const handle_t& c) {
        size_t i;
        for (i = 0; i < chopped_seqs.size(); ++i) {
            if (chopped.get_sequence(c) == chopped_seqs[i]) {
                chopped_found[i] = true;
                chopped_handles[i] = c;
                break;
            }
        }
        REQUIRE(i != chopped_seqs.size());
    });
    
    
    SECTION("ChoppedGraph methods work correctly") {
        
        REQUIRE(chopped.get_node_count() == chopped_handles.size());
        REQUIRE(chopped.get_total_length() == graph.get_total_length());
        
        nid_t min_id = numeric_limits<nid_t>::max();
        nid_t max_id = numeric_limits<nid_t>::min();
        size_t total_degree = 0;
        for (size_t i = 0; i < chopped_handles.size(); ++i) {
            
            REQUIRE(chopped_found[i]);
            handle_t c = chopped_handles[i];
            
            min_id = min(chopped.get_id(c), min_id);
            max_id = max(chopped.get_id(c), max_id);
            
            REQUIRE(c == chopped.flip(chopped.flip(c)));
            REQUIRE(chopped.get_sequence(chopped.flip(c)) == reverse_complement(chopped.get_sequence(c)));
            REQUIRE(chopped.get_length(c) == chopped.get_sequence(c).size());
            REQUIRE(c == chopped.get_handle(chopped.get_id(c)));
            
            for (auto d : {c, chopped.flip(c)}) {
                for (size_t j = 0; j < chopped.get_length(d); ++j) {
                    REQUIRE(chopped.get_base(d, j) == chopped.get_sequence(d)[j]);
                    for (size_t l = 0; l <= chopped.get_length(d); ++l) {
                        REQUIRE(chopped.get_subsequence(d, j, l) == chopped.get_sequence(d).substr(j, l));
                    }
                }
            }
            
            for (bool left : {true, false}) {
                set<handle_t> by_left, by_rev;
                chopped.follow_edges(c, left, [&](const handle_t& n) {
                    by_left.insert(n);
                });
                chopped.follow_edges(chopped.flip(c), !left, [&](const handle_t& n) {
                    by_rev.insert(chopped.flip(n));
                });
                REQUIRE(by_left == by_rev);
                REQUIRE(chopped.get_degree(c, left) == by_left.size());
                REQUIRE(chopped.get_degree(chopped.flip(c), !left) == by_left.size());
                total_degree += by_left.size();
            }
            
        }
        REQUIRE(chopped.get_edge_count() * 2 == total_degree);
        REQUIRE(chopped.max_node_id() == max_id);
        REQUIRE(chopped.min_node_id() == min_id);
        
        for (auto c : chopped_handles) {
            for (auto d : {c, chopped.flip(c)}) {
                for (auto e : chopped_handles) {
                    for (auto f : {e, chopped.flip(e)}) {
                        bool found_edge = false;
                        chopped.follow_edges(d, false, [&](const handle_t& n) {
                            found_edge = (n == f);
                            return !found_edge;
                        });
                        
                        REQUIRE(chopped.has_edge(d, f) == found_edge);
                    }
                }
            }
        }
    }
    
    vector<path_handle_t> unchopped_paths{p1, p2};
    vector<path_handle_t> paths;
    
    for (auto p : unchopped_paths) {
        REQUIRE(chopped.has_path(graph.get_path_name(p)));
        paths.push_back(chopped.get_path_handle(graph.get_path_name(p)));
    }
    
    SECTION("ChoppedPathGraph methods work correctly") {
        
        REQUIRE(!chopped.has_path("fake"));
        
        REQUIRE(chopped.get_path_count() == unchopped_paths.size());
        
        vector<bool> found_path(paths.size(), false);
        chopped.for_each_path_handle([&](const path_handle_t& q) {
            int count_found = 0;
            for (int i = 0; i < paths.size(); ++i) {
                if (paths[i] == q) {
                    found_path[i] = true;
                    ++count_found;
                }
            }
            REQUIRE(count_found == 1);
        });
        for (size_t i = 0; i < found_path.size(); ++i) {
            REQUIRE(found_path[i]);
        }
        
        for (size_t i = 0; i < unchopped_paths.size(); ++i) {
            auto p = unchopped_paths[i];
            auto q = paths[i];
            
            REQUIRE(chopped.get_path_name(q) == graph.get_path_name(p));
            REQUIRE(chopped.get_is_circular(q) == graph.get_is_circular(p));
            REQUIRE(chopped.get_step_count(q) >= graph.get_step_count(p));
            
            string path_seq, chopped_path_seq, chopped_rev_path_seq;
            for (handle_t h : graph.scan_path(p)) {
                path_seq += graph.get_sequence(h);
            }
            bool did_first_iter = false;
            size_t step_count1 = 0, step_count2 = 0;;
            for (auto s = chopped.path_begin(q);
                 s != chopped.path_end(q) && !(did_first_iter && s == chopped.path_begin(q));
                 s = chopped.get_next_step(s)) {
                
                chopped_path_seq += chopped.get_sequence(chopped.get_handle_of_step(s));
                REQUIRE(chopped.get_path_handle_of_step(s) == q);
                REQUIRE(chopped.has_next_step(s) == (chopped.get_is_circular(q) ||
                                                     chopped.get_next_step(s) != chopped.path_end(q)));
                REQUIRE(chopped.has_previous_step(s) == (chopped.get_is_circular(q) ||
                                                         chopped.get_previous_step(s) != chopped.path_front_end(q)));
                ++step_count1;
                did_first_iter = true;
            }
            did_first_iter = false;
            for (auto s = chopped.path_back(q);
                 s != chopped.path_front_end(q) && !(did_first_iter && s == chopped.path_back(q));
                 s = chopped.get_previous_step(s)) {
                
                chopped_rev_path_seq += chopped.get_sequence(chopped.flip(chopped.get_handle_of_step(s)));
                REQUIRE(chopped.get_path_handle_of_step(s) == q);
                ++step_count2;
                did_first_iter = true;
            }
            REQUIRE(chopped.get_step_count(p) == step_count1);
            REQUIRE(step_count1 == step_count2);
            REQUIRE(chopped_path_seq == path_seq);
            REQUIRE(chopped_rev_path_seq == reverse_complement(path_seq));
        }
        
        for (size_t i = 0; i < chopped_handles.size(); ++i) {
            handle_t c = chopped_handles[i];
            size_t correct_count;
            //if (i == 4 || i == 5) {
            if (i == 3 || i == 4) {
                correct_count = 1;
            }
            else {
                correct_count = 2;
            }
            REQUIRE(chopped.get_step_count(c) == correct_count);
            size_t direct_count = 0;
            chopped.for_each_step_on_handle(c, [&](const step_handle_t& s) {
                REQUIRE(c == chopped.forward(chopped.get_handle_of_step(s)));
                ++direct_count;
            });
            REQUIRE(direct_count == correct_count);
        }
    }
    
    SECTION("ChoppedPathPositionGraph methods work correctly") {
        
        for (size_t i = 0; i < paths.size(); ++i) {
            auto p = unchopped_paths[i];
            auto q = paths[i];
                        
            REQUIRE(chopped.get_path_length(q) == pos_overlay.get_path_length(p));
            
            size_t observed_pos = 0;
            chopped.for_each_step_in_path(q, [&](const step_handle_t& s) {
                REQUIRE(chopped.get_position_of_step(s) == observed_pos);
                size_t len = chopped.get_length(chopped.get_handle_of_step(s));
                for (size_t j = observed_pos; j < observed_pos + len; ++j) {
                    REQUIRE(chopped.get_step_at_position(q, j) == s);
                }
                observed_pos += len;
            });
            REQUIRE(chopped.get_path_length(q) == observed_pos);
        }
        
    }
}

}
}
