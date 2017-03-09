//
//  xg.cpp
//  
// Tests for xg algorithms
//

#include "catch.hpp"
#include "vg.hpp"
#include "xg.hpp"
#include <stdio.h>

namespace vg {
    namespace unittest {
        using namespace std;

        TEST_CASE("Path-based distance approximation in XG produces expected results", "[xg][mapping]") {
            
            VG vg;
            
            Node* n0 = vg.create_node("CGA");
            Node* n1 = vg.create_node("TTGG");
            Node* n2 = vg.create_node("CCGT");
            Node* n3 = vg.create_node("C");
            Node* n4 = vg.create_node("GT");
            Node* n5 = vg.create_node("GATAA");
            Node* n6 = vg.create_node("CGG");
            Node* n7 = vg.create_node("ACA");
            Node* n8 = vg.create_node("GCCG");
            Node* n9 = vg.create_node("A");
            Node* n10 = vg.create_node("C");
            Node* n11 = vg.create_node("G");
            Node* n12 = vg.create_node("T");
            Node* n13 = vg.create_node("A");
            Node* n14 = vg.create_node("C");
            Node* n15 = vg.create_node("C");
            
            vg.create_edge(n0, n1);
            vg.create_edge(n2, n0, true, true);
            vg.create_edge(n1, n3);
            vg.create_edge(n2, n3);
            vg.create_edge(n3, n4, false, true);
            vg.create_edge(n4, n5, true, false);
            vg.create_edge(n5, n6);
            vg.create_edge(n8, n6, false, true);
            vg.create_edge(n6, n7, false, true);
            vg.create_edge(n7, n9, true, true);
            vg.create_edge(n9, n10, true, false);
            vg.create_edge(n10, n11, false, false);
            vg.create_edge(n12, n11, false, true);
            vg.create_edge(n13, n12, false, false);
            vg.create_edge(n14, n13, true, false);
            vg.create_edge(n15, n14, true, true);
            
            Graph graph = vg.graph;
            
            Path* path = graph.add_path();
            path->set_name("path");
            Mapping* mapping = path->add_mapping();
            mapping->mutable_position()->set_node_id(n0->id());
            mapping->set_rank(1);
            mapping = path->add_mapping();
            mapping->mutable_position()->set_node_id(n2->id());
            mapping->set_rank(2);
            mapping = path->add_mapping();
            mapping->mutable_position()->set_node_id(n3->id());
            mapping->set_rank(3);
            mapping = path->add_mapping();
            mapping->mutable_position()->set_node_id(n4->id());
            mapping->mutable_position()->set_is_reverse(true);
            mapping->set_rank(4);
            mapping = path->add_mapping();
            mapping->mutable_position()->set_node_id(n5->id());
            mapping->set_rank(5);
            mapping = path->add_mapping();
            mapping->mutable_position()->set_node_id(n6->id());
            mapping->set_rank(6);
            mapping = path->add_mapping();
            mapping->mutable_position()->set_node_id(n8->id());
            mapping->mutable_position()->set_is_reverse(true);
            mapping->set_rank(7);
            
            xg::XG xg_index(graph);
            
            SECTION("Distance approxmation produces exactly correct path distances when positions are on path") {
                
                int64_t dist = xg_index.closest_shared_path_oriented_distance(n2->id(), 0, false,
                                                                              n5->id(), 0, false, 10);
                REQUIRE(dist == 7);
            }
            
            SECTION("Distance approxmation produces negative distance when order is reversed when positions are on path") {
                
                int64_t dist = xg_index.closest_shared_path_oriented_distance(n5->id(), 0, false,
                                                                              n2->id(), 0, false, 10);
                REQUIRE(dist == -7);
                
            }
            
            SECTION("Distance approxmation produces correctly signed distances on reverse strand when positions are on path") {
                
                int64_t dist = xg_index.closest_shared_path_oriented_distance(n5->id(),
                                                                              n5->sequence().size(),
                                                                              true,
                                                                              n2->id(),
                                                                              n2->sequence().size(),
                                                                              true,
                                                                              10);
                REQUIRE(dist == 7);
                
                dist = xg_index.closest_shared_path_oriented_distance(n2->id(),
                                                                      n2->sequence().size(),
                                                                      true,
                                                                      n5->id(),
                                                                      n5->sequence().size(),
                                                                      true,
                                                                      10);
                REQUIRE(dist == -7);
                
            }
            
            SECTION("Distance approxmation produces expected distances when positions are not on path") {
                
                int64_t dist = xg_index.closest_shared_path_oriented_distance(n1->id(),
                                                                              3,
                                                                              false,
                                                                              n7->id(),
                                                                              3,
                                                                              true,
                                                                              10);
                REQUIRE(dist == 15);
                
                dist = xg_index.closest_shared_path_oriented_distance(n7->id(),
                                                                      3,
                                                                      true,
                                                                      n1->id(),
                                                                      3,
                                                                      false ,
                                                                      10);
                REQUIRE(dist == -15);
                
            }
            
            SECTION("Distance approxmation produces expected output when search to path must traverse every type of edge") {
                
                int64_t dist = xg_index.closest_shared_path_oriented_distance(n1->id(),
                                                                              3,
                                                                              false,
                                                                              n15->id(),
                                                                              1,
                                                                              false,
                                                                              20);
                REQUIRE(dist == 22);
                
            }
            
            SECTION("Distance approxmation produces expected output when positions are on opposite strands") {
                int64_t dist = xg_index.closest_shared_path_oriented_distance(n1->id(),
                                                                              3,
                                                                              true,
                                                                              n7->id(),
                                                                              3,
                                                                              true,
                                                                              10);
                REQUIRE(dist == std::numeric_limits<int64_t>::max());
            }
            
            SECTION("Distance approxmation produces expected output when positions are in disconnected components") {
                int64_t dist = xg_index.closest_shared_path_oriented_distance(n1->id(),
                                                                              3,
                                                                              false,
                                                                              n7->id(),
                                                                              3,
                                                                              true,
                                                                              0);
                REQUIRE(dist == std::numeric_limits<int64_t>::max());
            }
            
        }
    }
}
