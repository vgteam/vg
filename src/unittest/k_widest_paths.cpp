/// \file yen.cpp
///  
/// Unit tests for the Yen's K Widest Paths algorithm
///

#include <iostream>
#include <string>
#include "../algorithms/k_widest_paths.hpp"
#include "../handle.hpp"
#include "../json2pb.h"
#include "../proto_handle_graph.hpp"
#include "catch.hpp"

#include <vg/vg.pb.h>

#include <bdsg/hash_graph.hpp>


namespace vg {
namespace unittest {
using namespace std;

using bdsg::HashGraph;

TEST_CASE("K Widest Paths works correctly on tiny example", "[k_widest_paths][algorithms]") {
    
    HashGraph graph;
    
    handle_t start = graph.create_handle("GAT");
    handle_t middle = graph.create_handle("TA");
    handle_t snp1 = graph.create_handle("C");
    handle_t snp2 = graph.create_handle("T");
    handle_t end = graph.create_handle("A");
    
    graph.create_edge(start, middle);
    graph.create_edge(middle, snp1);
    graph.create_edge(middle, snp2);
    graph.create_edge(snp1, end);
    graph.create_edge(snp2, end);

    unordered_map<handle_t, double> node_weights;
    node_weights[start] = 1.;
    node_weights[middle] = 0.75;
    node_weights[snp1] = 1.;
    node_weights[snp2] = 0.5;
    node_weights[end] = 2.;

    unordered_map<edge_t, double> edge_weights;
    edge_weights[graph.edge_handle(start, middle)] = 1.;
    edge_weights[graph.edge_handle(middle, snp1)] = 1.;
    edge_weights[graph.edge_handle(middle, snp2)] = 1.;
    edge_weights[graph.edge_handle(snp1, end)] = 1.;
    edge_weights[graph.edge_handle(snp2, end)] = 1.;
        
    function<double(handle_t)> node_weight_callback = [&](handle_t h) {
        if (!node_weights.count(h)) {
            // canonicalize to forward strand so we can use callback in both directions
            h = graph.flip(h);
        }
        return node_weights[h];
    };

    function<double(edge_t)> edge_weight_callback = [&](edge_t e) {
        return edge_weights[e];
    };
    
    // Track what we reach and at what distance
    unordered_map<handle_t, size_t> seen;
    
    SECTION("Does widest_dijkstra find the widest path (through snp2)") {

        double width;
        vector<handle_t> path;
        tie(width, path) = algorithms::widest_dijkstra(&graph, start, end, node_weight_callback, edge_weight_callback,
                                                       [&](handle_t) { return false; },
                                                       [&](edge_t) { return false; });


        REQUIRE(width == 0.75);
        REQUIRE(path == vector<handle_t>({start, middle, snp1, end}));        
    }

    SECTION("Does widest_dijkstra find the widest path (through snp2) in reverse direction") {

        double width;
        vector<handle_t> path;
        tie(width, path) = algorithms::widest_dijkstra(&graph, graph.flip(end), graph.flip(start), node_weight_callback, edge_weight_callback,
                                                       [&](handle_t) { return false; },
                                                       [&](edge_t) { return false; });


        REQUIRE(width == 0.75);
        REQUIRE(path == vector<handle_t>({graph.flip(end), graph.flip(snp1), graph.flip(middle), graph.flip(start)}));
    }


    SECTION("Does yens algorithm give same path when K=1") {

        vector<pair<double, vector<handle_t>>> widest_paths;
        widest_paths = algorithms::yens_k_widest_paths(&graph, start, end, 1, node_weight_callback, edge_weight_callback);

        REQUIRE(widest_paths.size() == 1);
        REQUIRE(widest_paths[0].first == 0.75);
        REQUIRE(widest_paths[0].second == vector<handle_t>({start, middle, snp1, end}));
        
    }

    SECTION("Does yens algorithm find a second path when K=2") {

        vector<pair<double, vector<handle_t>>> widest_paths;
        widest_paths = algorithms::yens_k_widest_paths(&graph, start, end, 2, node_weight_callback, edge_weight_callback);

        REQUIRE(widest_paths.size() == 2);
        REQUIRE(widest_paths[0].first == 0.75);
        REQUIRE(widest_paths[0].second == vector<handle_t>({start, middle, snp1, end}));
        REQUIRE(widest_paths[1].first == 0.5);
        REQUIRE(widest_paths[1].second == vector<handle_t>({start, middle, snp2, end}));
        
        
    }

    SECTION("Does yens algorithm find a second path when K>2") {

        vector<pair<double, vector<handle_t>>> widest_paths;
        widest_paths = algorithms::yens_k_widest_paths(&graph, start, end, 5, node_weight_callback, edge_weight_callback);

        REQUIRE(widest_paths.size() == 2);
        REQUIRE(widest_paths[0].first == 0.75);
        REQUIRE(widest_paths[0].second == vector<handle_t>({start, middle, snp1, end}));
        REQUIRE(widest_paths[1].first == 0.5);
        REQUIRE(widest_paths[1].second == vector<handle_t>({start, middle, snp2, end}));
        
    }

    SECTION("Does yesn algorithm find the equivalent paths when looking from sink to source?") {
        vector<pair<double, vector<handle_t>>> widest_paths;
        widest_paths = algorithms::yens_k_widest_paths(&graph, graph.flip(end), graph.flip(start), 5, node_weight_callback, edge_weight_callback);

        REQUIRE(widest_paths.size() == 2);
        REQUIRE(widest_paths[0].first == 0.75);
        REQUIRE(widest_paths[0].second == vector<handle_t>({graph.flip(end), graph.flip(snp1), graph.flip(middle), graph.flip(start)}));
        REQUIRE(widest_paths[1].first == 0.5);
        REQUIRE(widest_paths[1].second == vector<handle_t>({graph.flip(end), graph.flip(snp2), graph.flip(middle), graph.flip(start)}));
        
    }

}

TEST_CASE("K Widest Paths works correctly on small Wikipedia example", "[k_widest_paths][algorithms]") {
    
    HashGraph graph;

    handle_t C = graph.create_handle("C");
    handle_t D = graph.create_handle("A");
    handle_t E = graph.create_handle("C");
    handle_t F = graph.create_handle("A");
    handle_t G = graph.create_handle("C");
    handle_t H = graph.create_handle("A");
    
    graph.create_edge(C, D);
    graph.create_edge(C, E);
    graph.create_edge(D, E);
    graph.create_edge(D, F);
    graph.create_edge(E, F);
    graph.create_edge(E, G);
    graph.create_edge(F, G);
    graph.create_edge(F, H);
    graph.create_edge(G, H);

    unordered_map<edge_t, double> edge_weights;
    edge_weights[graph.edge_handle(C, D)] = 3.;
    edge_weights[graph.edge_handle(C, E)] = 2.;
    edge_weights[graph.edge_handle(D, E)] = 1.;
    edge_weights[graph.edge_handle(D, F)] = 4.;
    edge_weights[graph.edge_handle(E, F)] = 2.;
    edge_weights[graph.edge_handle(E, G)] = 3.;
    edge_weights[graph.edge_handle(F, G)] = 2.;
    edge_weights[graph.edge_handle(F, H)] = 1.;
    edge_weights[graph.edge_handle(G, H)] = 2.;

    function<double(handle_t)> node_weight_callback = [&](handle_t h) {
        return 10;
    };

    function<double(edge_t)> edge_weight_callback = [&](edge_t e) {
        return edge_weights[e];
    };
    
    // Track what we reach and at what distance
    unordered_map<handle_t, size_t> seen;
    
    SECTION("Does yens algorithm give the correct paths as verified by hand") {

        vector<pair<double, vector<handle_t>>> widest_paths;
        widest_paths = algorithms::yens_k_widest_paths(&graph, C, H, 100, node_weight_callback, edge_weight_callback);

        // correct number of paths
        REQUIRE(widest_paths.size() == 8);
        // paths listed in decreasing order based on score
        for (int i = 1; i < widest_paths.size(); ++i) {
            REQUIRE(widest_paths[i].first <= widest_paths[i-1].first);
        }
        // check each path as verified by hand
        map<vector<handle_t>, double> path_to_score;
        for (auto sp : widest_paths) {
            path_to_score[sp.second] = sp.first;
        }
        REQUIRE(path_to_score.at(vector<handle_t>({C, D, F, H})) == 1.);
        REQUIRE(path_to_score.at(vector<handle_t>({C, D, F, G, H})) == 2.);
        REQUIRE(path_to_score.at(vector<handle_t>({C, D, E, G, H})) == 1.);
        REQUIRE(path_to_score.at(vector<handle_t>({C, D, E, F, H})) == 1.);
        REQUIRE(path_to_score.at(vector<handle_t>({C, D, E, F, G, H})) == 1.);
        REQUIRE(path_to_score.at(vector<handle_t>({C, E, G, H})) == 2.);
        REQUIRE(path_to_score.at(vector<handle_t>({C, E, F, H})) == 1.);
        REQUIRE(path_to_score.at(vector<handle_t>({C, E, F, G, H})) == 2.);
        
    }



}


}
}
        
