#include <stdio.h>
#include <iostream>
#include <set>
#include "json2pb.h"
#include "vg.hpp"
#include "catch.hpp"
#include "snarls.hpp"
#include "distance.hpp"
#include "testTraversal.cpp"
#include "genotypekit.hpp"
#include <fstream>
namespace vg {
    namespace unittest {
        TEST_CASE( "Test distance index",
                   "[dist]" ) {
        VG graph;
        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("CTGA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n8);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n6);
        Edge* e5 = graph.create_edge(n3, n4);
        Edge* e6 = graph.create_edge(n3, n5);
        Edge* e7 = graph.create_edge(n4, n5);
        Edge* e8 = graph.create_edge(n5, n7);
        Edge* e9 = graph.create_edge(n6, n7);
        Edge* e10 = graph.create_edge(n7, n8);

        Snarl top_snarl;
        top_snarl.mutable_start()->set_node_id(n1->id());
        top_snarl.mutable_end()->set_node_id(n8->id());
          
        vector<vector<Snarl>> top_chains;
        top_chains.emplace_back();
        auto& top_chain1 = top_chains.back();
        top_chain1.emplace_back();
        auto& nested_snarl1 = top_chain1.back();

        nested_snarl1.mutable_start()->set_node_id(n2->id());
        nested_snarl1.mutable_end()->set_node_id(n7->id());
        nested_snarl1.set_type(ULTRABUBBLE);
        *nested_snarl1.mutable_parent() = top_snarl;
        nested_snarl1.set_directed_acyclic_net_graph(true);
        nested_snarl1.set_start_self_reachable(false);
        nested_snarl1.set_end_self_reachable(false);
        nested_snarl1.set_start_end_reachable(true);
        
        vector<Snarl> top_unary_snarls;

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls();
        DistanceIndex di = makeDistanceIndex(&graph, snarl_manager); 
        
        SECTION( "TRAVERSE GRAPH" ) {
            REQUIRE(true);
        }
    }
    }    
}
