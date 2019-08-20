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

}
}
