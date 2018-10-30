//
//  xg.cpp
//  
// Tests for xg algorithms
//

#include "catch.hpp"
#include "vg.hpp"
#include "xg.hpp"
#include "graph.hpp"
#include <stdio.h>

namespace vg {
    namespace unittest {
        using namespace std;

TEST_CASE("We can build an xg index on a nice graph", "[xg]") {

    string graph_json = R"(
    {"node":[{"id":1,"sequence":"GATT"},
    {"id":2,"sequence":"ACA"}],
    "edge":[{"to":2,"from":1}]}
    )";
    
    // Load the JSON
    Graph proto_graph;
    json2pb(proto_graph, graph_json.c_str(), graph_json.size());
    
    // Build the xg index
    xg::XG xg_index(proto_graph);

    Graph graph = xg_index.graph_context_id(make_pos_t(1, false, 0), 100);
    sort_by_id_dedup_and_clean(graph);

    REQUIRE(graph.node_size() == 2);
    REQUIRE(graph.edge_size() == 1);

}

TEST_CASE("We can build an xg index on a nasty graph", "[xg]") {

    string graph_json = R"(
    {"node":[{"id":1,"sequence":"GATT"},
    {"id":2,"sequence":"ACA"},
    {"id":9999,"sequence":"AAA"}],
    "edge":[{"to":2,"from":1}]}
    )";
    
    // Load the JSON
    Graph proto_graph;
    json2pb(proto_graph, graph_json.c_str(), graph_json.size());
    
    // Build the xg index
    xg::XG xg_index(proto_graph);

    Graph graph = xg_index.graph_context_id(make_pos_t(1, false, 0), 100);
    sort_by_id_dedup_and_clean(graph);

    REQUIRE(graph.node_size() == 2);
    REQUIRE(graph.edge_size() == 1);

}

TEST_CASE("We can build an xg index on a very nasty graph", "[xg]") {
        // We have a node 9999 in here to bust some MEM we don't want, to trigger the condition we are trying to test
    string graph_json = R"(
    {"node":[{"id":1430,"sequence":"GAGATCGTGCTACCGCACTCCATGCACTCTAG"},
    {"id":1428,"sequence":"C"},
    {"id":1429,"sequence":"T"},
    {"id":1431,"sequence":"CCTGGGCAACAGAACGAGATGCTGTC"},
    {"id":1427,"sequence":"AGGTTGGAGTGAGC"},
    {"id":1432,"sequence":"ACA"},
    {"id":1433,"sequence":"ACAACAACAACAACAA"},
    {"id":1426,"sequence":"GAGGCAGGAAAATCACTTGAACCGGGAGGCGG"},
    {"id":1434,"sequence":"T"},
    {"id":1435,"sequence":"C"},
    {"id":1424,"sequence":"C"},
    {"id":1425,"sequence":"T"},
    {"id":1436,"sequence":"AACAACAACAA"},
    {"id":1423,"sequence":"TCGGGAGGC"},
    {"id":1437,"sequence":"T"},
    {"id":1438,"sequence":"C"},
    {"id":1421,"sequence":"T"},
    {"id":1422,"sequence":"C"},
    {"id":1439,"sequence":"AACAACAACAACAA"},
    {"id":1420,"sequence":"TA"},
    {"id":1440,"sequence":"A"},
    {"id":1441,"sequence":"C"},
    {"id":1418,"sequence":"A"},
    {"id":1419,"sequence":"C"},
    {"id":1442,"sequence":"AA"},
    {"id":1417,"sequence":"TGCCTGTAATCCCAG"},
    {"id":1443,"sequence":"C"},
    {"id":1444,"sequence":"A"},
    {"id":1416,"sequence":"AAATACAAGTATTAGCCAGGCATTGTGGCAGG"},
    {"id":1445,"sequence":"TTCTCACATCTAAAACAGAGTTCCTGGTTCCA"},
    {"id":9999,"sequence":"TGTTGTTGTTGTTGTGTTNCTCATTTCGTTGCCAAT"}
    ],"edge":[{"to":1431,"from":1430},
    {"to":1430,"from":1428},
    {"to":1430,"from":1429},
    {"to":1432,"from":1431},
    {"to":1433,"from":1431},
    {"to":1428,"from":1427},
    {"to":1429,"from":1427},
    {"to":1433,"from":1432},
    {"to":1434,"from":1433},
    {"to":1435,"from":1433},
    {"to":1427,"from":1426},
    {"to":1436,"from":1434},
    {"to":1436,"from":1435},
    {"to":1426,"from":1424},
    {"to":1426,"from":1425},
    {"to":1437,"from":1436},
    {"to":1438,"from":1436},
    {"to":1424,"from":1423},
    {"to":1425,"from":1423},
    {"to":1439,"from":1437},
    {"to":1439,"from":1438},
    {"to":1423,"from":1421},
    {"to":1423,"from":1422},
    {"to":1440,"from":1439},
    {"to":1441,"from":1439},
    {"to":1421,"from":1420},
    {"to":1422,"from":1420},
    {"to":1442,"from":1440},
    {"to":1442,"from":1441},
    {"to":1420,"from":1418},
    {"to":1420,"from":1419},
    {"to":1443,"from":1442},
    {"to":1444,"from":1442},
    {"to":1418,"from":1417},
    {"to":1419,"from":1417},
    {"to":1445,"from":1443},
    {"to":1445,"from":1444},
    {"to":1417,"from":1416}],"path":[{"name":"17","mapping":[{"position":{"node_id":1416},"rank":1040},
    {"position":{"node_id":1417},"rank":1041},
    {"position":{"node_id":1419},"rank":1042},
    {"position":{"node_id":1420},"rank":1043},
    {"position":{"node_id":1422},"rank":1044},
    {"position":{"node_id":1423},"rank":1045},
    {"position":{"node_id":1425},"rank":1046},
    {"position":{"node_id":1426},"rank":1047},
    {"position":{"node_id":1427},"rank":1048},
    {"position":{"node_id":1429},"rank":1049},
    {"position":{"node_id":1430},"rank":1050},
    {"position":{"node_id":1431},"rank":1051},
    {"position":{"node_id":1433},"rank":1052},
    {"position":{"node_id":1435},"rank":1053},
    {"position":{"node_id":1436},"rank":1054},
    {"position":{"node_id":1438},"rank":1055},
    {"position":{"node_id":1439},"rank":1056},
    {"position":{"node_id":1441},"rank":1057},
    {"position":{"node_id":1442},"rank":1058},
    {"position":{"node_id":1444},"rank":1059},
    {"position":{"node_id":1445},"rank":1060}]}]}
    )";
    
    // Load the JSON
    Graph proto_graph;
    json2pb(proto_graph, graph_json.c_str(), graph_json.size());

    sort_by_id_dedup_and_clean(proto_graph);
    // Build the xg index
    xg::XG xg_index(proto_graph);

    SECTION("Context extraction gets something") {
        Graph graph = xg_index.graph_context_id(make_pos_t(1420, false, 0), 30);
        REQUIRE(graph.node_size() > 0);
    }
    
    SECTION("Extracting path ranges works as expected") {
        Graph graph;
        
        SECTION("We can extract within a single node") {
            xg_index.get_path_range("17", 5, 15, graph, false);
            
            // We should just get node 1416
            REQUIRE(graph.node_size() == 1);
            REQUIRE(graph.node(0).id() == 1416);
            // And the one dangling edge
            REQUIRE(graph.edge_size() == 1);
            for (size_t i = 0; i < graph.edge_size(); i++) {
                // All the edges should be normal
                REQUIRE(graph.edge(i).from_start() == false);
                REQUIRE(graph.edge(i).to_end() == false);
            }
        }
        
        SECTION("We can extract across two nodes") {
            xg_index.get_path_range("17", 5, 40, graph, false);
            
            // We should just get node 1416 and 1417
            REQUIRE(graph.node_size() == 2);
            REQUIRE(graph.node(0).id() >= 1416);
            REQUIRE(graph.node(0).id() <= 1417);
            REQUIRE(graph.node(1).id() >= 1416);
            REQUIRE(graph.node(1).id() <= 1417);
            REQUIRE(graph.node(0).id() != graph.node(1).id());
            // And the one real and 2 dangling edges
            REQUIRE(graph.edge_size() == 3);
            for (size_t i = 0; i < graph.edge_size(); i++) {
                // All the edges should be normal
                REQUIRE(graph.edge(i).from_start() == false);
                REQUIRE(graph.edge(i).to_end() == false);
            }
        }
        
        
    }

}

TEST_CASE("We can build the xg index on a small graph with discontinuous node ids that don't start at 1", "[xg]") {

    string graph_json = R"(
    {"node":[{"id":10,"sequence":"GATT"},
    {"id":20,"sequence":"ACA"}],
    "edge":[{"to":20,"from":10}]}
    )";

    // Load the JSON
    Graph proto_graph;
    json2pb(proto_graph, graph_json.c_str(), graph_json.size());

    sort_by_id_dedup_and_clean(proto_graph);
    // Build the xg index
    xg::XG xg_index(proto_graph);

    Graph graph = xg_index.graph_context_id(make_pos_t(10, false, 0), 100);
    sort_by_id_dedup_and_clean(graph);

    REQUIRE(graph.node_size() == 2);
    REQUIRE(graph.edge_size() == 1);


}

TEST_CASE("Target to alignment extraction", "[xg-target-to-aln]") {

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

    SECTION("Subpath getting gives us the expected 1bp alignment") {
        Alignment target = xg_index.target_alignment("path", 1, 2, "feature", false);
        REQUIRE(alignment_from_length(target) == 2 - 1);
    }

    SECTION("Subpath getting gives us the expected 10bp alignment") {
        Alignment target = xg_index.target_alignment("path", 10, 20, "feature", false);
        REQUIRE(alignment_from_length(target) == 20 - 10);
    }

    SECTION("Subpath getting gives us the expected 14bp alignment") {
        Alignment target = xg_index.target_alignment("path", 0, 14, "feature", false);
        REQUIRE(alignment_from_length(target) == 14);
    }

    SECTION("Subpath getting gives us the expected 21bp alignment") {
        Alignment target = xg_index.target_alignment("path", 0, 21, "feature", false);
        REQUIRE(alignment_from_length(target) == 21);
    }

    SECTION("Subpath getting gives us the expected inverted 7bp alignment") {
        Alignment target = xg_index.target_alignment("path", 0, 7, "feature", true);
        REQUIRE(alignment_from_length(target) == 7);
        REQUIRE(target.path().mapping(0).position().node_id() == n2->id());
        REQUIRE(target.path().mapping(1).position().node_id() == n0->id());
        REQUIRE(target.path().mapping(0).position().is_reverse() == true);
        REQUIRE(target.path().mapping(1).position().is_reverse() == true);
    }

}

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
                                                                      n5->id(), 0, false, false, 10);
        REQUIRE(dist == 7);
    }
    
    SECTION("Distance approxmation produces negative distance when order is reversed when positions are on path") {
        
        int64_t dist = xg_index.closest_shared_path_oriented_distance(n5->id(), 0, false,
                                                                      n2->id(), 0, false, false, 10);
        REQUIRE(dist == -7);
        
    }
    
    SECTION("Distance approxmation produces correctly signed distances on reverse strand when positions are on path") {
        
        int64_t dist = xg_index.closest_shared_path_oriented_distance(n5->id(),
                                                                      n5->sequence().size(),
                                                                      true,
                                                                      n2->id(),
                                                                      n2->sequence().size(),
                                                                      true,
                                                                      false,
                                                                      10);
        REQUIRE(dist == 7);
        
        dist = xg_index.closest_shared_path_oriented_distance(n2->id(),
                                                              n2->sequence().size(),
                                                              true,
                                                              n5->id(),
                                                              n5->sequence().size(),
                                                              true,
                                                              false,
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
                                                                      false,
                                                                      10);
        REQUIRE(dist == 15);
        
        dist = xg_index.closest_shared_path_oriented_distance(n7->id(),
                                                              3,
                                                              true,
                                                              n1->id(),
                                                              3,
                                                              false ,
                                                              false,
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
                                                                      false,
                                                                      10);
        REQUIRE(dist == std::numeric_limits<int64_t>::max());
    }
    
    SECTION("Distance approxmation produces expected output when cannot reach each other within the maximum distance") {
        int64_t dist = xg_index.closest_shared_path_oriented_distance(n1->id(),
                                                                      3,
                                                                      false,
                                                                      n7->id(),
                                                                      3,
                                                                      true,
                                                                      false,
                                                                      0);
        REQUIRE(dist == std::numeric_limits<int64_t>::max());
    }
    
    SECTION("Distance jumping produces expected result when start position and jump position are on path") {
        vector<tuple<int64_t, bool, size_t>> jump_pos = xg_index.jump_along_closest_path(n0->id(),
                                                                                         false,
                                                                                         1,
                                                                                         8,
                                                                                         10);
        
        REQUIRE(jump_pos.size() == 1);
        REQUIRE(get<0>(jump_pos[0]) == n4->id());
        REQUIRE(get<1>(jump_pos[0]) == true);
        REQUIRE(get<2>(jump_pos[0]) == 1);
        
        jump_pos = xg_index.jump_along_closest_path(n0->id(),
                                                    true,
                                                    2,
                                                    -8,
                                                    10);
        
        
        REQUIRE(jump_pos.size() == 1);
        REQUIRE(get<0>(jump_pos[0]) == n4->id());
        REQUIRE(get<1>(jump_pos[0]) == false);
        REQUIRE(get<2>(jump_pos[0]) == 1);
        
        jump_pos = xg_index.jump_along_closest_path(n4->id(),
                                                    true,
                                                    1,
                                                    -8,
                                                    10);
        
        
        REQUIRE(jump_pos.size() == 1);
        REQUIRE(get<0>(jump_pos[0]) == n0->id());
        REQUIRE(get<1>(jump_pos[0]) == false);
        REQUIRE(get<2>(jump_pos[0]) == 1);
        
        jump_pos = xg_index.jump_along_closest_path(n4->id(),
                                                    false,
                                                    1,
                                                    8,
                                                    10);
        
        
        REQUIRE(jump_pos.size() == 1);
        REQUIRE(get<0>(jump_pos[0]) == n0->id());
        REQUIRE(get<1>(jump_pos[0]) == true);
        REQUIRE(get<2>(jump_pos[0]) == 2);
    }
    
    SECTION("Distance jumping doesn't go past the end of a path") {
        vector<tuple<int64_t, bool, size_t>> jump_pos = xg_index.jump_along_closest_path(n0->id(),
                                                                                         false,
                                                                                         1,
                                                                                         200,
                                                                                         10);
        
        REQUIRE(jump_pos.empty());
    }
    
    SECTION("Distance jumping produces expected results when it needs to traverse an edge to the path") {
        vector<tuple<int64_t, bool, size_t>> jump_pos = xg_index.jump_along_closest_path(n1->id(),
                                                                                         false,
                                                                                         1,
                                                                                         9,
                                                                                         10);
        
        REQUIRE(jump_pos.size() == 1);
        REQUIRE(get<0>(jump_pos[0]) == n5->id());
        REQUIRE(get<1>(jump_pos[0]) == false);
        REQUIRE(get<2>(jump_pos[0]) == 3);
        
        jump_pos = xg_index.jump_along_closest_path(n1->id(),
                                                    false,
                                                    3,
                                                    7,
                                                    10);
        
        REQUIRE(jump_pos.size() == 1);
        REQUIRE(get<0>(jump_pos[0]) == n5->id());
        REQUIRE(get<1>(jump_pos[0]) == false);
        REQUIRE(get<2>(jump_pos[0]) == 3);
        
        jump_pos = xg_index.jump_along_closest_path(n1->id(),
                                                    true,
                                                    3,
                                                    -9,
                                                    10);
        
        REQUIRE(jump_pos.size() == 1);
        REQUIRE(get<0>(jump_pos[0]) == n5->id());
        REQUIRE(get<1>(jump_pos[0]) == true);
        REQUIRE(get<2>(jump_pos[0]) == 2);
        
        jump_pos = xg_index.jump_along_closest_path(n1->id(),
                                                    true,
                                                    1,
                                                    -7,
                                                    10);
        
        REQUIRE(jump_pos.size() == 1);
        REQUIRE(get<0>(jump_pos[0]) == n5->id());
        REQUIRE(get<1>(jump_pos[0]) == true);
        REQUIRE(get<2>(jump_pos[0]) == 2);
    }
    
    SECTION("Distance jumping doesn't search past the maximum search length") {
        vector<tuple<int64_t, bool, size_t>> jump_pos = xg_index.jump_along_closest_path(n1->id(),
                                                                                         true,
                                                                                         2,
                                                                                         9,
                                                                                         1);
        
        REQUIRE(jump_pos.empty());
    }
    
    
    SECTION("Distance jumping can traverse every edge type") {
        vector<tuple<int64_t, bool, size_t>> jump_pos = xg_index.jump_along_closest_path(n15->id(),
                                                                                         false,
                                                                                         1,
                                                                                         -12,
                                                                                         15);
        
        REQUIRE(jump_pos.size() == 1);
        REQUIRE(get<0>(jump_pos[0]) == n6->id());
        REQUIRE(get<1>(jump_pos[0]) == false);
        REQUIRE(get<2>(jump_pos[0]) == 1);
    }
}
        
TEST_CASE("Path component memoization produces expected results", "[xg]") {
    
    string graph_json = R"({"node": [{"sequence": "AAACCC", "id": 1}, {"sequence": "CACACA", "id": 2}, {"sequence": "CACACA", "id": 3}, {"sequence": "TTTTGG", "id": 4}, {"sequence": "ACGTAC", "id": 5}], "path": [{"name": "one", "mapping": [{"position": {"node_id": 1}, "rank": 1}, {"position": {"node_id": 2}, "rank": 2}]}, {"name": "three", "mapping": [{"position": {"node_id": 2}, "rank": 1}, {"position": {"node_id": 3}, "rank": 2}]}, {"name": "two", "mapping": [{"position": {"node_id": 4}, "rank": 1}, {"position": {"node_id": 5}, "rank": 2}]}], "edge": [{"from": 1, "to": 2}, {"from": 2, "to": 3}, {"from": 4, "to": 5}]})";
    
    // Load the JSON
    Graph proto_graph;
    json2pb(proto_graph, graph_json.c_str(), graph_json.size());
    
    // Build the xg index
    xg::XG xg_index(proto_graph);
    
    for (int node_id : {1, 2, 3}) {
        for (size_t path_rank : xg_index.paths_of_node(node_id)){
            for (int other_node_id : {1, 2, 3}) {
                for (size_t other_path_rank : xg_index.paths_of_node(other_node_id)) {
                    REQUIRE(xg_index.paths_on_same_component(path_rank, other_path_rank));
                }
            }
            
            for (int other_node_id : {4, 5}) {
                for (size_t other_path_rank : xg_index.paths_of_node(other_node_id)) {
                    REQUIRE(!xg_index.paths_on_same_component(path_rank, other_path_rank));
                }
            }
        }
    }
    
    for (int node_id : {4, 5}) {
        for (size_t path_rank : xg_index.paths_of_node(node_id)){
            for (int other_node_id : {1, 2, 3}) {
                for (size_t other_path_rank : xg_index.paths_of_node(other_node_id)) {
                    REQUIRE(!xg_index.paths_on_same_component(path_rank, other_path_rank));
                }
            }
            
            for (int other_node_id : {4, 5}) {
                for (size_t other_path_rank : xg_index.paths_of_node(other_node_id)) {
                    REQUIRE(xg_index.paths_on_same_component(path_rank, other_path_rank));
                }
            }
        }
    }
}

TEST_CASE("Looping over XG handles in parallel works", "[xg]") {

    string graph_json = R"(
    {"node":[{"id":1,"sequence":"GATT"},
    {"id":2,"sequence":"ACA"}],
    "edge":[{"to":2,"from":1}]}
    )";
    
    // Load the JSON
    Graph proto_graph;
    json2pb(proto_graph, graph_json.c_str(), graph_json.size());
    
    // Build the xg index
    xg::XG xg_index(proto_graph);

    size_t count = 0;

    xg_index.for_each_handle([&](const handle_t& got) {
        #pragma omp critical
        count++;
    }, true);
    
    REQUIRE(count == 2);

}

}
}
