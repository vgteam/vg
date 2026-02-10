//
//  xg.cpp
//  
// Tests for xg algorithms
//

#include "catch.hpp"
#include "vg.hpp"
#include "xg.hpp"
#include "graph.hpp"
#include "algorithms/subgraph.hpp"
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
    xg::XG xg_index;
    xg_index.from_path_handle_graph(VG(proto_graph));

    VG vg_graph;
    algorithms::extract_context(xg_index, vg_graph, xg_index.get_handle(1), 0, 100);
    Graph& graph = vg_graph.graph;
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
    xg::XG xg_index;
    xg_index.from_path_handle_graph(VG(proto_graph));

    VG vg_graph;
    algorithms::extract_context(xg_index, vg_graph, xg_index.get_handle(1), 0, 100);
    Graph& graph = vg_graph.graph;
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
    xg::XG xg_index;
    xg_index.from_path_handle_graph(VG(proto_graph));

    SECTION("Context extraction gets something") {
        VG graph;
        algorithms::extract_context(xg_index, graph, xg_index.get_handle(1420), 0, 30);
        REQUIRE(graph.get_node_count() > 0);
    }
    
    SECTION("Extracting path ranges works as expected") {
        VG graph;
        
        SECTION("We can extract within a single node") {
            algorithms::extract_path_range(xg_index, xg_index.get_path_handle("17"), 5, 15, graph);
            
            // We should just get node 1416
            REQUIRE(graph.graph.node_size() == 1);
            REQUIRE(graph.graph.node(0).id() == 1416);
            // And no dangling edges
            REQUIRE(graph.graph.edge_size() == 0);
            for (size_t i = 0; i < graph.graph.edge_size(); i++) {
                // All the edges should be normal
                REQUIRE(graph.graph.edge(i).from_start() == false);
                REQUIRE(graph.graph.edge(i).to_end() == false);
            }
        }
        
        SECTION("We can extract across two nodes") {
            algorithms::extract_path_range(xg_index, xg_index.get_path_handle("17"), 5, 40, graph);
            
            // We should just get node 1416 and 1417
            REQUIRE(graph.graph.node_size() == 2);
            REQUIRE(graph.graph.node(0).id() >= 1416);
            REQUIRE(graph.graph.node(0).id() <= 1417);
            REQUIRE(graph.graph.node(1).id() >= 1416);
            REQUIRE(graph.graph.node(1).id() <= 1417);
            REQUIRE(graph.graph.node(0).id() != graph.graph.node(1).id());
            // And the one real and 0 dangling edges
            REQUIRE(graph.graph.edge_size() == 1);
            for (size_t i = 0; i < graph.graph.edge_size(); i++) {
                // All the edges should be normal
                REQUIRE(graph.graph.edge(i).from_start() == false);
                REQUIRE(graph.graph.edge(i).to_end() == false);
            }
        }
    }
    
    SECTION("Looking up steps by positions works") {
        auto path = xg_index.get_path_handle("17");
        
        step_handle_t found;
        size_t first_found = 0;
        
        for (size_t i = 0; i < 150; i++) {
            // Scan a bunch of the path
            step_handle_t found_here = xg_index.get_step_at_position(path, i);
            if (i == 0 || found_here != found) {
                if (i > 0) {
                    // How long was the last thing we saw?
                    size_t length = i - first_found;
                    // Make sure it is right
                    REQUIRE(length == xg_index.get_length(xg_index.get_handle_of_step(found)));
                }
                
                // Remember what we saw here
                found = found_here;
                first_found = i;
            }
        }
    }

}

TEST_CASE("We can build and scan an XG index for a problematic graph", "[xg]") {
    string graph_json = R"(
    {"node":[
      {"id":1,"sequence":"AAA"},
      {"id":2,"sequence":"C"},
      {"id":3,"sequence":"G"},
      {"id":5,"sequence":"T"},
      {"id":4,"sequence":"AAA"}
    ],
    "edge":[
      {"to":2,"from":1},
      {"to":3,"from":1},
      {"to":4,"from":2},
      {"to":5,"from":3},
      {"to":4,"from":5}
    ],
    "path":[
      {"name":"reference","mapping":[
        {"position":{"node_id":1},"rank":1},
        {"position":{"node_id":2},"rank":2},
        {"position":{"node_id":4},"rank":3}
      ]}
    ]}
    )";
    
    // Load the JSON
    Graph proto_graph;
    json2pb(proto_graph, graph_json.c_str(), graph_json.size());

    // Build the xg index (without any sorting)
    xg::XG xg_index;
    xg_index.from_path_handle_graph(VG(proto_graph));

    REQUIRE(xg_index.get_node_count() == 5);
    
    // We will cross every edge from every node, so this should come to 2 * edge count.
    size_t edge_obs_count = 0;
    
    xg_index.for_each_handle([&](const handle_t& here) {
        xg_index.follow_edges(here, false, [&](const handle_t& there) {
            edge_obs_count++;
        });
        xg_index.follow_edges(here, true, [&](const handle_t& there) {
            edge_obs_count++;
        });
    });
    
    REQUIRE(edge_obs_count == 10);
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
    xg::XG xg_index;
    xg_index.from_path_handle_graph(VG(proto_graph));

    VG vg_graph;
    algorithms::extract_context(xg_index, vg_graph, xg_index.get_handle(10), 0, 100);
    Graph& graph = vg_graph.graph;
    sort_by_id_dedup_and_clean(graph);

    REQUIRE(graph.node_size() == 2);
    REQUIRE(graph.edge_size() == 1);


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
    xg::XG xg_index;
    xg_index.from_path_handle_graph(VG(proto_graph));

    size_t count = 0;

    xg_index.for_each_handle([&](const handle_t& got) {
        #pragma omp critical
        count++;
    }, true);
    
    REQUIRE(count == 2);

}


TEST_CASE("Vectorization of xg works correctly", "[xg]") {
    string graph_json = R"(
    {"edge": [
        {"from": "5", "to": "6"},
        {"from": "7", "to": "9"},
        {"from": "12", "to": "13"},
        {"from": "12", "to": "14"},
        {"from": "8", "to": "9"},
        {"from": "1", "to": "2"},
        {"from": "1", "to": "3"},
        {"from": "4", "to": "6"},
        {"from": "6", "to": "7"},
        {"from": "6", "to": "8"},
        {"from": "2", "to": "4"},
        {"from": "2", "to": "5"},
        {"from": "10", "to": "12"},
        {"from": "9", "to": "10"},
        {"from": "9", "to": "11"},
        {"from": "11", "to": "12"},
        {"from": "13", "to": "15"},
        {"from": "14", "to": "15"},
        {"from": "3", "to": "4"},
        {"from": "3", "to": "5"}
    ], "node": [
        {"id": "5", "sequence": "C"},
        {"id": "7", "sequence": "A"},
        {"id": "12", "sequence": "ATAT"},
        {"id": "8", "sequence": "G"},
        {"id": "1", "sequence": "CAAATAAG"},
        {"id": "4", "sequence": "T"},
        {"id": "6", "sequence": "TTG"},
        {"id": "15", "sequence": "CCAACTCTCTG"},
        {"id": "2", "sequence": "A"},
        {"id": "10", "sequence": "A"},
        {"id": "9", "sequence": "AAATTTTCTGGAGTTCTAT"},
        {"id": "11", "sequence": "T"},
        {"id": "13", "sequence": "A"},
        {"id": "14", "sequence": "T"},
        {"id": "3", "sequence": "G"}
    ], "path": [
        {"mapping": [
            {"edit": [{"from_length": 8, "to_length": 8}], "position": {"node_id": "1"}, "rank": "1"},
            {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "3"}, "rank": "2"},
            {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "5"}, "rank": "3"},
            {"edit": [{"from_length": 3, "to_length": 3}], "position": {"node_id": "6"}, "rank": "4"},
            {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "8"}, "rank": "5"},
            {"edit": [{"from_length": 19, "to_length": 19}], "position": {"node_id": "9"}, "rank": "6"},
            {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "11"}, "rank": "7"},
            {"edit": [{"from_length": 4, "to_length": 4}], "position": {"node_id": "12"}, "rank": "8"},
            {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "14"}, "rank": "9"},
            {"edit": [{"from_length": 11, "to_length": 11}], "position": {"node_id": "15"}, "rank": "10"}
        ], "name": "x"}]}
    )";
    
    // Load the JSON
    Graph proto_graph;
    json2pb(proto_graph, graph_json.c_str(), graph_json.size());

    // Build the xg index (without any sorting)
    xg::XG xg_index;
    xg_index.from_path_handle_graph(VG(proto_graph));

    REQUIRE(xg_index.get_node_count() == 15);
    
    SECTION("edge ranks are unique") {
        // Collect all the unique edge ranks we observe for all the edges in the XG.
        unordered_set<size_t> unique_edge_ranks;
        xg_index.for_each_edge([&](const edge_t& edge) {
            size_t edge_index = xg_index.edge_index(edge);
#ifdef debug
            cerr << "Edge " << xg_index.get_id(edge.first) << (xg_index.get_is_reverse(edge.first) ? '-' : '+') << " -> "
                << xg_index.get_id(edge.second) << (xg_index.get_is_reverse(edge.second) ? '-' : '+')
                << " has index " << edge_index << endl;
#endif
            unique_edge_ranks.insert(edge_index);
        });
        
        REQUIRE(unique_edge_ranks.size() == xg_index.get_edge_count());
    }
    
    SECTION("edge ranks are defined for all ways of articulating edges that exist") {
        vector<handle_t> forwards;
        vector<handle_t> reverse;
        xg_index.for_each_handle([&](const handle_t& h) {
            forwards.push_back(h);
            reverse.push_back(xg_index.flip(h));
        });
        
        for (size_t i = 0; i < forwards.size(); i++) {
            for (size_t j = 0; j < forwards.size(); j++) {
                for (bool flip1 : {false, true}) {
                    for (bool flip2 : {false, true}) {
                        // For all possible combinations of handles and orientations
                        handle_t from = (flip1 ? reverse : forwards)[i];
                        handle_t to = (flip2 ? reverse : forwards)[j];
                        if (xg_index.has_edge(from, to)) {
                            // If the edge exists
                            
                            // Make sure it exists the other way
                            REQUIRE(xg_index.has_edge(xg_index.flip(to), xg_index.flip(from)));
                            
                            // Make sure that both directions have the same index, even if someone didn't use edge_handle.
                            REQUIRE(xg_index.edge_index(std::make_pair(from, to)) ==
                                xg_index.edge_index(std::make_pair(xg_index.flip(to), xg_index.flip(from))));
                        }
                    }
                }
            }
        }
        
    }
    
    
    SECTION("node offsets are unique and map back to the right nodes") {
        // Do the same for the node vector offsets, except mapping offset to node ID. 
        map<size_t, nid_t> id_at_offset;
        xg_index.for_each_handle([&](const handle_t& h) {
            nid_t node = xg_index.get_id(h);
            // TODO: We're assuming this is 0-based. See https://github.com/vgteam/libhandlegraph/issues/41.
            size_t offset = xg_index.node_vector_offset(node);
            
#ifdef debug
            if (id_at_offset.count(offset)) {
                cerr << "Found offset " << offset << " for node " << node << " occupied by " << id_at_offset.at(offset) << endl;
            } else {
                cerr << "Found free offset " << offset << " for node " << node << endl;
            }
            // TODO: We're providing 1-based indexes here. See https://github.com/vgteam/libhandlegraph/issues/41.
            cerr << "Mapping back offset " << offset << " gives node " << xg_index.node_at_vector_offset(offset + 1) << endl;
#endif
            
            id_at_offset[offset] = node;
        });
        
        // Make sure all the nodes have unique offsets.
        REQUIRE(id_at_offset.size() == xg_index.get_node_count());
        
        for (auto& kv : id_at_offset) {
            // Make sure each offset maps back to the right node.
            // TODO: We're providing 1-based indexes here. See https://github.com/vgteam/libhandlegraph/issues/41.
            REQUIRE(xg_index.node_at_vector_offset(kv.first + 1) == kv.second);
        }
        
        // Now make sure all intermediate offsets map properly.
        auto it = id_at_offset.begin();
        auto prev = it;
        ++it;
        while (it != id_at_offset.end()) {
            // Go through adjacent node starts in sequence order
            for (size_t i = prev->first; i < it->first; i++) {
                // Every base before the first base of the next node should belong to the previous node.
                // TODO: We're providing 1-based indexes here. See https://github.com/vgteam/libhandlegraph/issues/41.
                nid_t node_at_i = xg_index.node_at_vector_offset(i + 1);
                
#ifdef debug
                cerr << "Base at index " << i << " maps to node " << node_at_i << endl;
#endif
                
                REQUIRE(node_at_i == prev->second);
            }
            prev = it;
            ++it;
        }
        
    }
}




}
}
