/**
 * \file 
 * unittest/source_sink_overlay.cpp: test cases for the source and sink node adding overlay.
 */

#include "catch.hpp"

#include "random_graph.hpp"

#include "../source_sink_overlay.hpp"
#include "../kmer.hpp"
#include "../vg.hpp"
#include "vg/io/json2pb.h"

#include <iostream>
#include <vector>
#include <unordered_set>

namespace vg {
namespace unittest {

using namespace std;

TEST_CASE("SourceSinkOverlay adds a source and a sink to a 1-node graph", "[overlay]") {

    // Make a vg graph
    VG vg;
    Node* n1 = vg.create_node("GATTACA");
    
    // Make an overlay
    SourceSinkOverlay overlay(&vg, 10);
    
    SECTION("we see the right graph") {
    
        bool found_node = false;
        bool found_source = false;
        bool found_sink = false;
        
        // Make sure the nodes are as expected
        overlay.for_each_handle([&](const handle_t& handle) {
            
#ifdef debug
            cerr << "Observed node " << overlay.get_id(handle) << " orientation " << overlay.get_is_reverse(handle)
                << " with sequence " << overlay.get_sequence(handle) << endl;
#endif
        
            if (overlay.get_id(handle) == n1->id()) {
                // This is the real node
                REQUIRE(!found_node);
                found_node = true;
                
                REQUIRE(overlay.get_sequence(handle) == "GATTACA");
                REQUIRE(overlay.get_degree(handle, false) == 1);
                REQUIRE(overlay.get_degree(handle, true) == 1);
                
                unordered_set<pair<id_t, bool>> overlay_neighbors;
                overlay.follow_edges(handle, false, [&](const handle_t& neighbor) {
                    overlay_neighbors.emplace(overlay.get_id(neighbor), overlay.get_is_reverse(neighbor));
                });
                
                REQUIRE(overlay_neighbors.size() == 1);
                REQUIRE(overlay_neighbors.count(make_pair(overlay.get_id(overlay.get_sink_handle()), false)));
                
                overlay_neighbors.clear();
                overlay.follow_edges(handle, true, [&](const handle_t& neighbor) {
                    overlay_neighbors.emplace(overlay.get_id(neighbor), overlay.get_is_reverse(neighbor));
                });
                
                REQUIRE(overlay_neighbors.size() == 1);
                REQUIRE(overlay_neighbors.count(make_pair(overlay.get_id(overlay.get_source_handle()), false)));
                
            } else if (handle == overlay.get_source_handle()) {
                // This is the fake source
                REQUIRE(!found_source);
                found_source = true;
                
                REQUIRE(overlay.get_sequence(handle) == "##########");
                REQUIRE(overlay.get_degree(handle, false) == 1);
                REQUIRE(overlay.get_degree(handle, true) == 0);
                
                unordered_set<pair<id_t, bool>> overlay_neighbors;
                overlay.follow_edges(handle, false, [&](const handle_t& neighbor) {
                    overlay_neighbors.emplace(overlay.get_id(neighbor), overlay.get_is_reverse(neighbor));
                });
                
                REQUIRE(overlay_neighbors.size() == 1);
                REQUIRE(overlay_neighbors.count(make_pair(n1->id(), false)));
                
                overlay_neighbors.clear();
                overlay.follow_edges(handle, true, [&](const handle_t& neighbor) {
                    overlay_neighbors.emplace(overlay.get_id(neighbor), overlay.get_is_reverse(neighbor));
                });
                
                REQUIRE(overlay_neighbors.empty());
                
            } else if (handle == overlay.get_sink_handle()) {
                // This is the fake sink
                REQUIRE(!found_sink);
                found_sink = true;
                
                REQUIRE(overlay.get_sequence(handle) == "$$$$$$$$$$");
                REQUIRE(overlay.get_degree(handle, false) == 0);
                REQUIRE(overlay.get_degree(handle, true) == 1);
                
                unordered_set<pair<id_t, bool>> overlay_neighbors;
                overlay.follow_edges(handle, false, [&](const handle_t& neighbor) {
                    overlay_neighbors.emplace(overlay.get_id(neighbor), overlay.get_is_reverse(neighbor));
                });
                
                REQUIRE(overlay_neighbors.empty());
                
                overlay_neighbors.clear();
                overlay.follow_edges(handle, true, [&](const handle_t& neighbor) {
                    overlay_neighbors.emplace(overlay.get_id(neighbor), overlay.get_is_reverse(neighbor));
                });
                
                REQUIRE(overlay_neighbors.size() == 1);
                REQUIRE(overlay_neighbors.count(make_pair(n1->id(), false)));
                
            } else {
                // This is an unrecognized node!
                REQUIRE(false);
            }
        });
        
        REQUIRE(found_node);
        REQUIRE(found_source);
        REQUIRE(found_sink);
    }
    
}

TEST_CASE("SourceSinkOverlay agrees with VG::add_start_end_markers in a tiny graph", "[overlay]") {
    const string graph_json = R"({"node":[{"sequence":"CAAATAAG","id":"1"},{"sequence":"A","id":"2"},{"sequence":"G","id":"3"},{"sequence":"T","id":"4"},{"sequence":"C","id":"5"},{"sequence":"TTG","id":"6"},{"sequence":"A","id":"7"},{"sequence":"G","id":"8"},{"sequence":"AAATTTTCTGGAGTTCTAT","id":"9"},{"sequence":"A","id":"10"},{"sequence":"T","id":"11"},{"sequence":"ATAT","id":"12"},{"sequence":"A","id":"13"},{"sequence":"T","id":"14"},{"sequence":"CCAACTCTCTG","id":"15"}],"edge":[{"from":"1","to":"2"},{"from":"1","to":"3"},{"from":"2","to":"4"},{"from":"2","to":"5"},{"from":"3","to":"4"},{"from":"3","to":"5"},{"from":"4","to":"6"},{"from":"5","to":"6"},{"from":"6","to":"7"},{"from":"6","to":"8"},{"from":"7","to":"9"},{"from":"8","to":"9"},{"from":"9","to":"10"},{"from":"9","to":"11"},{"from":"10","to":"12"},{"from":"11","to":"12"},{"from":"12","to":"13"},{"from":"12","to":"14"},{"from":"13","to":"15"},{"from":"14","to":"15"}],"path":[{"name":"x","mapping":[{"position":{"node_id":"1"},"edit":[{"from_length":8,"to_length":8}],"rank":"1"},{"position":{"node_id":"3"},"edit":[{"from_length":1,"to_length":1}],"rank":"2"},{"position":{"node_id":"5"},"edit":[{"from_length":1,"to_length":1}],"rank":"3"},{"position":{"node_id":"6"},"edit":[{"from_length":3,"to_length":3}],"rank":"4"},{"position":{"node_id":"8"},"edit":[{"from_length":1,"to_length":1}],"rank":"5"},{"position":{"node_id":"9"},"edit":[{"from_length":19,"to_length":19}],"rank":"6"},{"position":{"node_id":"11"},"edit":[{"from_length":1,"to_length":1}],"rank":"7"},{"position":{"node_id":"12"},"edit":[{"from_length":4,"to_length":4}],"rank":"8"},{"position":{"node_id":"14"},"edit":[{"from_length":1,"to_length":1}],"rank":"9"},{"position":{"node_id":"15"},"edit":[{"from_length":11,"to_length":11}],"rank":"10"}]}]})";
    
    Graph graph;
    json2pb(graph, graph_json);
    
    VG produced(graph);
    
    id_t highest_id = produced.max_node_id();
    id_t start_id = highest_id + 1;
    id_t end_id = start_id + 1;
    
    Node* start_node = nullptr;
    Node* end_node = nullptr;
    
    VG copy = produced;
    copy.add_start_end_markers(10, '#', '$', start_node, end_node, start_id, end_id);
    
    SourceSinkOverlay overlay(&produced, 10, start_id, end_id);
    
    // We know this graph has no tipless components.
    
    SECTION("graphs agree") {
    
        // Get handles for all the nodes
        unordered_map<id_t, handle_t> copy_handles;
        copy.for_each_handle([&](const handle_t& handle) {
            copy_handles[copy.get_id(handle)] = handle;
        });
        
        unordered_map<id_t, handle_t> overlay_handles;
        overlay.for_each_handle([&](const handle_t& handle) {
            overlay_handles[overlay.get_id(handle)] = handle;
        });
        
        REQUIRE(copy_handles.size() == overlay_handles.size());
        
        for (auto& kv : copy_handles) {
            // Unpack the copy handle
            auto& id = kv.first;
            auto& copy_handle = kv.second;
            
            // Find the corresponding overlay handle
            REQUIRE(overlay_handles.count(id));
            auto& overlay_handle = overlay_handles.at(id);
            
            // Both should be forward
            REQUIRE(!copy.get_is_reverse(copy_handle));
            REQUIRE(!overlay.get_is_reverse(overlay_handle));
            
            // Both should have the same sequence
            REQUIRE(copy.get_sequence(copy_handle) == overlay.get_sequence(overlay_handle));
            
            for (bool go_left : {false, true}) {
                // Both should have the same neighbors
                unordered_set<pair<id_t, bool>> copy_neighbors;
                unordered_set<pair<id_t, bool>> overlay_neighbors;
                
#ifdef debug
                cerr << "Look " << (go_left ? "left" : "right") << " from " << id << endl;
#endif
                
                copy.follow_edges(copy_handle, go_left, [&](const handle_t& neighbor) {
#ifdef debug
                    cerr << "\tIn copy find node " << copy.get_id(neighbor) << " orientation " << copy.get_is_reverse(neighbor) << endl;
#endif
                    
                    copy_neighbors.emplace(copy.get_id(neighbor), copy.get_is_reverse(neighbor));
                });
                
                overlay.follow_edges(overlay_handle, go_left, [&](const handle_t& neighbor) {
#ifdef debug
                    cerr << "\tIn overlay find node " << overlay.get_id(neighbor) << " orientation " << overlay.get_is_reverse(neighbor) << endl;
#endif
                    
                    overlay_neighbors.emplace(overlay.get_id(neighbor), overlay.get_is_reverse(neighbor));
                });
                
                REQUIRE(copy_neighbors.size() == overlay_neighbors.size());
                for (auto& item : copy_neighbors) {
                    REQUIRE(overlay_neighbors.count(item));
                }
            }
        }
        
    }
        
    
    SECTION("kmer generation works") {
        // Now test kmer generation
        size_t size_limit = 10000;
        temp_file::remove(write_gcsa_kmers_to_tmpfile(copy, 10, size_limit, start_id, end_id));
        size_limit = 10000;
        temp_file::remove(write_gcsa_kmers_to_tmpfile(overlay, 10, size_limit, start_id, end_id));
    }
    
    SECTION("for_each_handle works in parallel mode") {
    
        size_t found = 0;
    
        overlay.for_each_handle([&](const handle_t& handle) {
            
            #pragma omp critical
            {
            
#ifdef debug
                cerr << "Observed node " << overlay.get_id(handle) << " orientation " << overlay.get_is_reverse(handle)
                    << " with sequence " << overlay.get_sequence(handle) << endl;
#endif
        
            
                found++;
            }
        }, true);
        
        REQUIRE(found == produced.get_node_count() + 2);
    }
}

TEST_CASE("SourceSinkOverlay agrees with VG::add_start_end_markers in a random graph", "[overlay]") {

    for (size_t trial = 0; trial < 1000; trial++) {
        
        VG random;
        random_graph(100, 3, 30, &random);
        
#ifdef debug
        cerr << "Trial " << trial << ": " << pb2json(random.graph) << endl;
#endif
        
        id_t highest_id = random.max_node_id();
        id_t start_id = highest_id + 1;
        id_t end_id = start_id + 1;
        
        Node* start_node = nullptr;
        Node* end_node = nullptr;
        
        VG copy = random;
        copy.add_start_end_markers(10, '#', '$', start_node, end_node, start_id, end_id);
        
        SourceSinkOverlay overlay(&random, 10, start_id, end_id);
        
        // Now we compare copy and overlay. They ought to match. Except that if
        // we have to break into a component with no real tips, we might break
        // in at a different place. So work out which components to not worry about so much.
        unordered_set<id_t> in_tipless_component;
        vector<pair<unordered_set<id_t>, vector<handle_t>>> components = handlealgs::weakly_connected_components_with_tips(&random);
        for (auto& component : components) {
            if (component.second.empty()) {
                for (auto& id : component.first) {
                    // Mark all the nodes in this tipless component as in a tipless component
                    in_tipless_component.insert(id);
                }
            }
        }
        
        // Get handles for all the nodes
        unordered_map<id_t, handle_t> copy_handles;
        copy.for_each_handle([&](const handle_t& handle) {
            copy_handles[copy.get_id(handle)] = handle;
        });
        
        unordered_map<id_t, handle_t> overlay_handles;
        overlay.for_each_handle([&](const handle_t& handle) {
            overlay_handles[overlay.get_id(handle)] = handle;
        });
        
        REQUIRE(copy_handles.size() == overlay_handles.size());
        
        for (auto& kv : copy_handles) {
            // Unpack the copy handle
            auto& id = kv.first;
            auto& copy_handle = kv.second;
            
            // Find the corresponding overlay handle
            REQUIRE(overlay_handles.count(id));
            auto& overlay_handle = overlay_handles.at(id);
            
            // Both should be forward
            REQUIRE(!copy.get_is_reverse(copy_handle));
            REQUIRE(!overlay.get_is_reverse(overlay_handle));
            
            // Both should have the same sequence
            REQUIRE(copy.get_sequence(copy_handle) == overlay.get_sequence(overlay_handle));
            
            for (bool go_left : {false, true}) {
                // Both should have the same neighbors
                unordered_set<pair<id_t, bool>> copy_neighbors;
                unordered_set<pair<id_t, bool>> overlay_neighbors;
                
#ifdef debug
                cerr << "Look " << (go_left ? "left" : "right") << " from " << id << endl;
#endif
                
                copy.follow_edges(copy_handle, go_left, [&](const handle_t& neighbor) {
#ifdef debug
                    cerr << "\tIn copy find node " << copy.get_id(neighbor) << " orientation " << copy.get_is_reverse(neighbor) << endl;
#endif
                    
                    if ((in_tipless_component.count(id) && (copy.get_id(neighbor) == start_id || copy.get_id(neighbor) == end_id)) ||
                        (in_tipless_component.count(copy.get_id(neighbor)) && (id == start_id || id == end_id))) {
                        // This edge constitutes breaking into a tipless component. Ignore it.
#ifdef debug
                        cerr << "\t\tDrop entry to tipless component" << endl;
#endif
                        return;
                    }
                    
                    copy_neighbors.emplace(copy.get_id(neighbor), copy.get_is_reverse(neighbor));
                });
                
                overlay.follow_edges(overlay_handle, go_left, [&](const handle_t& neighbor) {
#ifdef debug
                    cerr << "\tIn overlay find node " << overlay.get_id(neighbor) << " orientation " << overlay.get_is_reverse(neighbor) << endl;
#endif
                    
                    if ((in_tipless_component.count(id) && (overlay.get_id(neighbor) == start_id || overlay.get_id(neighbor) == end_id)) ||
                        (in_tipless_component.count(overlay.get_id(neighbor)) && (id == start_id || id == end_id))) {
                        // This edge constitutes breaking into a tipless component. Ignore it.
#ifdef debug
                        cerr << "\t\tDrop entry to tipless component" << endl;
#endif
                        return;
                    }
                     
                    overlay_neighbors.emplace(overlay.get_id(neighbor), overlay.get_is_reverse(neighbor));
                });
                
                REQUIRE(copy_neighbors.size() == overlay_neighbors.size());
                for (auto& item : copy_neighbors) {
                    REQUIRE(overlay_neighbors.count(item));
                }
            }
        }
        
    
    }

}

}
}
