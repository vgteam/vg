#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <set>
#include "json2pb.h"
#include "vg.hpp"
#include "catch.hpp"
#include "snarls.hpp"
#include "distance.hpp"
#include "genotypekit.hpp"
#include "random_graph.hpp"
#include "seed_clusterer.hpp"
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
    TEST_CASE( "Revese in chain","[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("G");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("G");
        Node* n9 = graph.create_node("AA");
        Node* n10 = graph.create_node("G");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n10);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n4);
        Edge* e5 = graph.create_edge(n3, n5);
        Edge* e6 = graph.create_edge(n4, n5);
        Edge* e7 = graph.create_edge(n5, n6);
        Edge* e8 = graph.create_edge(n5, n7);
        Edge* e9 = graph.create_edge(n6, n7);
        Edge* e10 = graph.create_edge(n7, n7, false, true);
        Edge* e11 = graph.create_edge(n7, n8);
        Edge* e12 = graph.create_edge(n7, n9);
        Edge* e13 = graph.create_edge(n8, n9);
        Edge* e14 = graph.create_edge(n9, n10);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls();
        DistanceIndex dist_index (&graph, &snarl_manager, 20);

        SnarlSeedClusterer clusterer;

        SECTION( "One cluster" ) {
            vector<pos_t> seeds;
            seeds.push_back(make_pos_t(3, false, 0));
            seeds.push_back(make_pos_t(4, false, 0));

            vector<hash_set<size_t>> clusters = clusterer.cluster_seeds(
                                        seeds, 10,  snarl_manager, dist_index); 
            for (hash_set<size_t> c : clusters) { 
                for (size_t s : c) {cerr << s << " ";}
                cerr << endl; 
            }
            REQUIRE( clusters.size() == 1);
        }
    }//end test case


    TEST_CASE( "Clusters in snarl","[cluster]" ) {
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

        SECTION( "Two clusters in a chain and loop of snarl boundary" ) {
            vector<pos_t> seeds;
            seeds.push_back(make_pos_t(2, false, 0));
            seeds.push_back(make_pos_t(3, false, 0));
            seeds.push_back(make_pos_t(5, false, 0));
            seeds.push_back(make_pos_t(16, false, 0));
            //New cluster
            seeds.push_back(make_pos_t(5, false, 10));
            seeds.push_back(make_pos_t(6, false, 0));
            seeds.push_back(make_pos_t(8, false, 0));

            vector<hash_set<size_t>> clusters = clusterer.cluster_seeds(
                                         seeds, 4,  snarl_manager, dist_index); 
            for (hash_set<size_t> c : clusters) { 
                for (size_t s : c) {cerr << s << " ";}
                cerr << endl; 
            }
            REQUIRE( clusters.size() == 2);
            REQUIRE (( (clusters[0].count(0) == 1 &&
                       clusters[0].count(1) == 1 &&
                       clusters[0].count(2) == 1 &&
                       clusters[0].count(3) == 1 &&
                       clusters[1].count(4) == 1 &&
                       clusters[1].count(5) == 1 &&
                       clusters[1].count(6) == 1)  ||

                     ( clusters[1].count(0) == 1 &&
                       clusters[1].count(1) == 1 &&
                       clusters[1].count(2) == 1 &&
                       clusters[0].count(3) == 1 &&
                       clusters[0].count(4) == 1 &&
                       clusters[0].count(5) == 1 &&
                       clusters[0].count(6) == 1  )));
        }
        SECTION( "Five clusters" ) {
            vector<pos_t> seeds;
            seeds.push_back(make_pos_t(3, false, 0));
            seeds.push_back(make_pos_t(5, false, 0));
            //New cluster
            seeds.push_back(make_pos_t(5, false, 10));
            //New cluster
            seeds.push_back(make_pos_t(6, false, 0));
            seeds.push_back(make_pos_t(8, false, 0));
            //New cluster
            seeds.push_back(make_pos_t(16, false, 0));
            //New cluster
            seeds.push_back(make_pos_t(13, false, 1));
            seeds.push_back(make_pos_t(14, false, 0));
            seeds.push_back(make_pos_t(15, false, 0));

            vector<hash_set<size_t>> clusters = clusterer.cluster_seeds(
                                         seeds, 3,  snarl_manager, dist_index); 
            for (hash_set<size_t> c : clusters) { 
                for (size_t s : c) {cerr << s << " ";}
                cerr << endl; 
            }
            REQUIRE( clusters.size() == 5);
        }
    }//end test case

    TEST_CASE("Random graphs", "[cluster]"){

        for (int i = 0; i < 1; i++) {
            // For each random graph
            VG graph;
            random_graph(1000, 20, 100, &graph);


            CactusSnarlFinder bubble_finder(graph);
            SnarlManager snarl_manager = bubble_finder.find_snarls();
            DistanceIndex dist_index (&graph, &snarl_manager, 20);

            SnarlSeedClusterer clusterer;

            vector<const Snarl*> allSnarls;
            auto addSnarl = [&] (const Snarl* s) {
                allSnarls.push_back(s);
            };
            snarl_manager.for_each_snarl_preorder(addSnarl);

            uniform_int_distribution<int> randSnarlIndex(0, allSnarls.size()-1);
            default_random_engine generator(time(NULL));
            for (size_t k = 0; k < 100 ; k++) {
                vector<pos_t> seeds;
                for (int j = 0; j < 100; j++) {
                    //Check clusters of 100 random positions 
                    const Snarl* snarl1 = allSnarls[randSnarlIndex(generator)];

                    pair<unordered_set<Node*>, unordered_set<Edge*>> contents1 =
                           snarl_manager.shallow_contents(snarl1, graph, true);
  
                    vector<Node*> nodes1 (contents1.first.begin(), contents1.first.end());


                    uniform_int_distribution<int> randNodeIndex1(0,nodes1.size()-1);
 
                    Node* node1 = nodes1[randNodeIndex1(generator)];
                    id_t nodeID1 = node1->id();
 
                    off_t offset1 = uniform_int_distribution<int>(0,node1->sequence().size() - 1)(generator);

                    pos_t pos = make_pos_t(nodeID1,
                        uniform_int_distribution<int>(0,1)(generator) == 0,offset1 );
                    seeds.push_back(pos);

                }
                int64_t lim = 20;// Distance between clusters
                vector<hash_set<size_t>> clusters = clusterer.cluster_seeds(
                                      seeds, lim, snarl_manager, dist_index); 

                for (size_t a = 0; a < clusters.size(); a++) {
                    // For each cluster -cluster this cluster to ensure that 
                    // there is only one
                    hash_set<size_t> clust = clusters[a];
                    
                    //vector of verified clusters. should be only one after 
                    //running through all seeds in cluster a 
                    vector<hash_set<size_t>> checked_clusters;
                    for (size_t i1 : clust) {
                        // For each clustered position
                        hash_set<size_t> assignments;
                        for (size_t b = 0 ; b < checked_clusters.size() ; b++) {
                            //For each new cluster that we're making

                            for (size_t i2 : checked_clusters[b]) {
                                //Ever seed is close to at least one seed in the same cluster
                                pos_t pos1 = seeds[i1];
                                pos_t pos2 = seeds[i2];
                                size_t len1 = graph.get_length(graph.get_handle(get_id(pos1), false));
                                pos_t rev1 = make_pos_t(get_id(pos1), 
                                                    !is_rev(pos1),
                                                    len1 - get_offset(pos1)); 
                                size_t len2 = graph.get_length(graph.get_handle(get_id(pos2), false));
                                pos_t rev2 = make_pos_t(get_id(pos2), 
                                                    !is_rev(pos2),
                                                    len2 - get_offset(pos2)); 
                                int64_t dist1 = dist_index.minDistance(pos1, pos2);
                                int64_t dist2 = dist_index.minDistance(pos1, rev2);
                                int64_t dist3 = dist_index.minDistance(rev1, pos2);
                                int64_t dist4 = dist_index.minDistance(rev1, rev2);
                                int64_t dist = DistanceIndex::minPos({dist1, 
                                                   dist2, dist3, dist4});
                                if ( dist != -1 && dist <= lim+2) {
                                    assignments.insert(b);
                                }

                            }
                        }
                        vector<hash_set<size_t>> new_checked_clusters;
                        hash_set<size_t> curr_cluster;
                        curr_cluster.insert(i1);
                        for (size_t b = 0 ; b < checked_clusters.size() ; b++) {
                            if (assignments.count(b) > 0) {
                                curr_cluster.insert(checked_clusters[b].begin(),
                                                    checked_clusters[b].end());
                            } else {
                                new_checked_clusters.push_back(checked_clusters[b]);
                            } 
                        }
                        new_checked_clusters.push_back(curr_cluster);
                        checked_clusters = move(new_checked_clusters);
                        for ( size_t b = 0; b < a ; b ++) {
                            // For each other cluster
                            hash_set<size_t> clust2 = clusters[b];
                            for (size_t i2 : clust2) {
                                // And each position in *that* cluster
                                pos_t pos1 = seeds[i1];
                                pos_t pos2 = seeds[i2];
                                size_t len1 = graph.get_length(graph.get_handle(get_id(pos1), false));
                                pos_t rev1 = make_pos_t(get_id(pos1), 
                                                    !is_rev(pos1),
                                                    len1 - get_offset(pos1)); 
                                size_t len2 = graph.get_length(graph.get_handle(get_id(pos2), false));
                                pos_t rev2 = make_pos_t(get_id(pos2), 
                                                    !is_rev(pos2),
                                                    len2 - get_offset(pos2)); 
                                int64_t dist1 = dist_index.minDistance(pos1, pos2);
                                int64_t dist2 = dist_index.minDistance(pos1, rev2);
                                int64_t dist3 = dist_index.minDistance(rev1, pos2);
                                int64_t dist4 = dist_index.minDistance(rev1, rev2);
                                if (!(    (dist1 == -1 || dist1 >= lim-2) 
                                      && (dist2 == -1 ||  dist2 >= lim-2) 
                                      && (dist3 == -1 ||  dist3 >= lim-2)  
                                      && (dist4 == -1 || dist4 >= lim-2))){
                                    //graph.serialize_to_file("testGraph");
                                    cerr << "These should have been in the same cluster: ";
                                    cerr << pos1 << " " << pos2 << endl;
                                    cerr << dist1 << " " << dist2 << " " << dist3<< " " << dist4 << endl;
                                    dist_index.printSelf();
                                REQUIRE(false);
                                };
                            }
                        }
                    }
                    if (checked_clusters.size() != 1) {
                        //graph.serialize_to_file("testGraph");
                        cerr << "These should be different clusters: " << endl;
                        for (hash_set<size_t> c : checked_clusters) {
                            cerr << "cluster: " ; 
                            for (size_t i1 : c) {
                                cerr << seeds[i1] << " ";
                            }
                            cerr << endl;
                        }
                    }
                    REQUIRE(checked_clusters.size() == 1);
                }
;
            }
        }
    } //end test case
}
}
