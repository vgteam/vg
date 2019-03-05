#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <set>
#include "json2pb.h"
#include "vg.hpp"
#include "catch.hpp"
#include "snarls.hpp"
#include "distance.hpp"
#include "genotypekit.hpp"
#include "random_graph.hpp"
#include "seed_clusterer.hpp"
#include <fstream>
#include <random>
#include <time.h>

//#define print

namespace vg {
namespace unittest {
    TEST_CASE( "multiple clusters in a chain",
                   "[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("CTGA");
        Node* n9 = graph.create_node("GCA");
        Node* n10 = graph.create_node("T");
        Node* n11 = graph.create_node("G");
        Node* n12 = graph.create_node("CTGA");
        Node* n13 = graph.create_node("GCA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n9);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n4);
        Edge* e5 = graph.create_edge(n3, n4);
        Edge* e6 = graph.create_edge(n4, n5);
        Edge* e7 = graph.create_edge(n4, n5, false, true);
        Edge* e8 = graph.create_edge(n5, n6);
        Edge* e9 = graph.create_edge(n5, n6, true, false);
        Edge* e10 = graph.create_edge(n6, n7);
        Edge* e11 = graph.create_edge(n6, n8);
        Edge* e12 = graph.create_edge(n7, n8);
        Edge* e13 = graph.create_edge(n8, n10);
        Edge* e14 = graph.create_edge(n9, n10);
        Edge* e15 = graph.create_edge(n10, n11);
        Edge* e16 = graph.create_edge(n10, n12);
        Edge* e17 = graph.create_edge(n11, n13);
        Edge* e18 = graph.create_edge(n12, n13);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls();

        const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
        const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
        const Snarl* snarl3 = snarl_manager.into_which_snarl(3, false);

        DistanceIndex dist_index (&graph, &snarl_manager, 20);
        SnarlSeedClusterer clusterer;

        SECTION( "One cluster" ) {
 
            id_t seed_nodes[] = {2, 3, 4, 7, 8, 9, 11};
            //all are in the same cluster
            vector<pos_t> seeds;
            for (id_t n : seed_nodes) {
                seeds.push_back(make_pos_t(n, false, 0));
            }

            vector<hash_set<size_t>> clusters = clusterer.cluster_seeds(seeds, 10, 
                           snarl_manager, dist_index); 
            REQUIRE(clusters.size() == 1); 

        }
        SECTION( "Two clusters" ) {
 
            vector<id_t> seed_nodes( {2, 3, 4, 7, 8, 10, 11});
            //Clusters should be {2, 3, 4}, {7, 8, 10, 11}
            //Distance from pos on 4 to pos on 7 is 9, including positions
            vector<pos_t> seeds;
            for (id_t n : seed_nodes) {
                seeds.push_back(make_pos_t(n, false, 0));
            }

            vector<hash_set<size_t>> clusters = clusterer.cluster_seeds(
                                         seeds, 9,  snarl_manager, dist_index); 
            for (hash_set<size_t> c : clusters) { 
                for (size_t s : c) {cerr << seed_nodes[s] << " ";}
                cerr << endl; 
            }
            REQUIRE( clusters.size() == 2);
            REQUIRE (( (clusters[0].count(0) == 1 &&
                       clusters[0].count(1) == 1 &&
                       clusters[0].count(2) == 1 &&
                       clusters[1].count(3) == 1 &&
                       clusters[1].count(4) == 1 &&
                       clusters[1].count(5) == 1 &&
                       clusters[1].count(6) == 1  ) ||

                     ( clusters[1].count(0) == 1 &&
                       clusters[1].count(1) == 1 &&
                       clusters[1].count(2) == 1 &&
                       clusters[0].count(3) == 1 &&
                       clusters[0].count(4) == 1 &&
                       clusters[0].count(5) == 1 &&
                       clusters[0].count(6) == 1  )));

        }
    }//End test case
TEST_CASE( "Clusters in snarl","[dist]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GGGGGGGGGGGG");//12 Gs
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("G");
        Node* n9 = graph.create_node("AA");
        Node* n10 = graph.create_node("G");
        Node* n11 = graph.create_node("G");
        Node* n12 = graph.create_node("G");
        Node* n13 = graph.create_node("GA");
        Node* n14 = graph.create_node("G");
        Node* n15 = graph.create_node("G");
        Node* n16 = graph.create_node("G");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n13);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n16);
        Edge* e27 = graph.create_edge(n16, n9);
        Edge* e5 = graph.create_edge(n3, n4);
        Edge* e6 = graph.create_edge(n3, n5);
        Edge* e7 = graph.create_edge(n4, n6);
        Edge* e8 = graph.create_edge(n5, n6);
        Edge* e9 = graph.create_edge(n6, n7);
        Edge* e10 = graph.create_edge(n6, n8);
        Edge* e11 = graph.create_edge(n7, n8);
        Edge* e12 = graph.create_edge(n8, n9);
        Edge* e13 = graph.create_edge(n9, n10);
        Edge* e14 = graph.create_edge(n9, n11);
        Edge* e15 = graph.create_edge(n10, n11);
        Edge* e16 = graph.create_edge(n11, n12);
        Edge* e17 = graph.create_edge(n11, n2);
        Edge* e18 = graph.create_edge(n12, n1);
        Edge* e19 = graph.create_edge(n13, n14);
        Edge* e20 = graph.create_edge(n13, n15);
        Edge* e21 = graph.create_edge(n14, n15);
        Edge* e22 = graph.create_edge(n15, n12);
        Edge* e23 = graph.create_edge(n2, n2, true, false);
        Edge* e24 = graph.create_edge(n11, n11, false, true);
        Edge* e25 = graph.create_edge(n1, n1, true, false);
        Edge* e26 = graph.create_edge(n12, n12, false, true);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls();
        DistanceIndex dist_index (&graph, &snarl_manager, 20);

        SnarlSeedClusterer clusterer;

        SECTION( "One cluster" ) {
        }
    }//end test case
}
}
