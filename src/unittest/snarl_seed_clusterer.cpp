#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <set>
#include "vg/io/json2pb.h"
#include "../vg.hpp"
#include "bdsg/hash_graph.hpp"
#include "catch.hpp"
#include "../snarls.hpp"
#include "../cactus_snarl_finder.hpp"
#include "../integrated_snarl_finder.hpp"
#include "../genotypekit.hpp"
#include "random_graph.hpp"
#include "../snarl_seed_clusterer.hpp"
#include <random>
#include <time.h>
#include <structures/union_find.hpp>

//#define print

namespace vg {
namespace unittest {

    TEST_CASE( "cluster one node",
                   "[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");


        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);
        SnarlDistanceIndexClusterer clusterer(dist_index, &graph);
        
        //graph.to_dot(cerr);

        SECTION( "One cluster" ) {
 
            id_t seed_nodes[] = {1, 1};
            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (bool use_minimizers : {true, false} ) {
                for (id_t n : seed_nodes) {
                    pos_t pos = make_pos_t(n, false, 0);
                    auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 10); 
                REQUIRE(clusters.size() == 1); 
            }


        }
    }
    TEST_CASE( "Looping chain", "[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("ACACGTTGC");
        Node* n2 = graph.create_node("TCTCCACCGGCAAGTTTCACTTCACTT");
        Node* n3 = graph.create_node("A");
        Node* n4 = graph.create_node("AT");
        Node* n5 = graph.create_node("CGTGGGG");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n5);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n4);
        Edge* e5 = graph.create_edge(n3, n4);
        Edge* e6 = graph.create_edge(n4, n5);


        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);
        SnarlDistanceIndexClusterer clusterer(dist_index, &graph);
        
        //graph.to_dot(cerr);

        SECTION( "Two cluster" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(make_pos_t(2, false, 1));
            positions.emplace_back(make_pos_t(2, true, 7));
            //all are in the same cluster
            for (bool use_minimizers : {true, false} ) {
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (auto& pos : positions) {
                    auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 15); 
                REQUIRE(clusters.size() == 2); 
            }


        }
    }
    TEST_CASE( "cluster one node with loop",
                   "[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCAATGGACA");

        Edge* e1 = graph.create_edge(n1, n1);


        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);
        SnarlDistanceIndexClusterer clusterer(dist_index, &graph);
        
        //graph.to_dot(cerr);

        SECTION( "One cluster" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(make_pos_t(1, false, 0));
            positions.emplace_back(make_pos_t(1, true, 0));
            //all are in the same cluster
            for (bool use_minimizers : {true, false} ) {
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (auto& pos : positions) {
                    auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0,chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 5); 
                REQUIRE(clusters.size() == 1); 
            }


        }
    }
    TEST_CASE( "two tips", "[cluster][bug]" ) {
        VG graph;

        Node* n1 = graph.create_node("AGGGAAGATGTCGTGAAG");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("GA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);


        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);
        SnarlDistanceIndexClusterer clusterer(dist_index, &graph);
        
        //graph.to_dot(cerr);

        SECTION( "One cluster on the same node" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(make_pos_t(2, false, 0));
            positions.emplace_back(make_pos_t(1, false, 5));
            //all are in the same cluster
            for (bool use_minimizers : {true, false} ) {
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (pos_t pos : positions) {
                    auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 15); 
                REQUIRE(clusters.size() == 1); 
            }
        }
    }


    TEST_CASE( "cluster simple chain",
                   "[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("T");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n2, n4);
        Edge* e4 = graph.create_edge(n3, n4);
        Edge* e5 = graph.create_edge(n4, n5);
        Edge* e6 = graph.create_edge(n4, n6);
        Edge* e7 = graph.create_edge(n5, n7);
        Edge* e8 = graph.create_edge(n6, n7);


        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);
        SnarlDistanceIndexClusterer clusterer(dist_index, &graph);
        
        //graph.to_dot(cerr);

        SECTION( "One cluster on the same node" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(make_pos_t(4, false, 0));
            positions.emplace_back(make_pos_t(4, false, 1));
            positions.emplace_back(make_pos_t(4, false, 3));
            //all are in the same cluster
            for (bool use_minimizers : {true, false} ) {
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (pos_t pos : positions) {
                    auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 2); 
                REQUIRE(clusters.size() == 1); 
            }


        }
        SECTION( "One cluster on opposite sides of a snp" ) {
 
            id_t seed_nodes[] = {2, 3, 5};
            //all are in the same cluster
            for (bool use_minimizers : {true, false} ) {
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (id_t n : seed_nodes) {
                    pos_t pos = make_pos_t(n, false, 0);
                    auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 10); 
                REQUIRE(clusters.size() == 1); 
            }


        }
        SECTION( "Three clusters on opposite sides of a snp" ) {
 
            id_t seed_nodes[] = {2, 3, 5};
            //all are in the same cluster
            for (bool use_minimizers : {true, false} ) {
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (id_t n : seed_nodes) {
                    pos_t pos = make_pos_t(n, false, 0);
                    auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0,chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 4); 
                REQUIRE(clusters.size() == 3); 
            }


        }
        SECTION( "Three clusters on opposite sides of a snp in two reads" ) {
 
            id_t seed_nodes[] = {2, 4};
            //all are in the same cluster
            vector<vector<SnarlDistanceIndexClusterer::Seed>> seeds (2);

            pos_t pos = make_pos_t(2, false, 0);
            seeds[0].push_back({ pos, 0});
            pos = make_pos_t(3, false, 0);
            seeds[0].push_back({ pos, 0});

            pos = make_pos_t(5, false, 0);
            seeds[1].push_back({ pos, 0});

            vector<vector<SnarlDistanceIndexClusterer::Cluster>> clusters = clusterer.cluster_seeds(seeds,  5, 5); 
            REQUIRE(clusters.size() == 2); 
            REQUIRE(clusters[0][0].fragment == clusters[1][0].fragment);
            REQUIRE(clusters[0].size() == 2);
            REQUIRE(clusters[1].size() == 1);
        }
        SECTION( "Three clusters on opposite sides of a snp in two reads going in the other direction in the chain" ) {
 
            //all are in the same cluster
            vector<vector<SnarlDistanceIndexClusterer::Seed>> seeds (2);

            pos_t pos = make_pos_t(5, false, 0);
            seeds[0].push_back({ pos, 0});
            pos = make_pos_t(6, false, 0);
            seeds[0].push_back({ pos, 0});

            pos = make_pos_t(1, false, 0);
            seeds[1].push_back({ pos, 0});

            vector<vector<SnarlDistanceIndexClusterer::Cluster>> clusters = clusterer.cluster_seeds(seeds, 10, 10); 
            REQUIRE(clusters.size() == 2); 
            REQUIRE(clusters[0][0].fragment == clusters[1][0].fragment);
            REQUIRE(clusters[0].size() == 2);
            REQUIRE(clusters[1].size() == 1);
        }
        SECTION( "Distances are correct" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(make_pos_t(1, false, 1));
            positions.emplace_back(make_pos_t(2, false, 0));
            positions.emplace_back(make_pos_t(4, false, 3));
            positions.emplace_back(make_pos_t(7, false, 0));
            //all are in the same cluster
            for (bool use_minimizers : {true, false} ) {
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (pos_t pos : positions) {
                    auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                REQUIRE(clusterer.distance_between_seeds(seeds[0], seeds[1], false) == 2); 
                REQUIRE(clusterer.distance_between_seeds(seeds[0], seeds[2], false) == 6); 
                REQUIRE(clusterer.distance_between_seeds(seeds[0], seeds[3], false) == 8); 
                REQUIRE(clusterer.distance_between_seeds(seeds[1], seeds[3], false) == 6); 
                REQUIRE(clusterer.distance_between_seeds(seeds[1], seeds[0], false) == 2); 
                REQUIRE(clusterer.distance_between_seeds(seeds[2], seeds[0], false) == 6); 
                REQUIRE(clusterer.distance_between_seeds(seeds[3], seeds[0], false) == 8); 
                REQUIRE(clusterer.distance_between_seeds(seeds[3], seeds[1], false) == 6); 
            }


        }
    }

    TEST_CASE( "cluster simple chain with multiple connected components",
                   "[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("T");
        Node* n8 = graph.create_node("TTTTTTTTT");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n2, n4);
        Edge* e4 = graph.create_edge(n3, n4);
        Edge* e5 = graph.create_edge(n4, n5);
        Edge* e6 = graph.create_edge(n4, n6);
        Edge* e7 = graph.create_edge(n5, n7);
        Edge* e8 = graph.create_edge(n6, n7);


        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);
        SnarlDistanceIndexClusterer clusterer(dist_index, &graph);
        
        //graph.to_dot(cerr);

        SECTION( "One cluster on the same node plus extra node" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(make_pos_t(4, false, 0));
            positions.emplace_back(make_pos_t(4, false, 1));
            positions.emplace_back(make_pos_t(4, false, 3));
            positions.emplace_back(make_pos_t(8, false, 3));
            //all are in the same cluster
            for (bool use_minimizers : {true, false} ) {
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (pos_t pos : positions) {
                    auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 2); 
                REQUIRE(clusters.size() == 2); 
            }


        }
        SECTION( "One cluster on opposite sides of a snp" ) {
 
            id_t seed_nodes[] = {2, 3, 5, 8};
            //all are in the same cluster
            for (bool use_minimizers : {true, false} ) {
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (id_t n : seed_nodes) {
                    pos_t pos = make_pos_t(n, false, 0);
                    auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds,  10); 
                REQUIRE(clusters.size() == 2); 
            }


        }
        SECTION( "Three clusters on opposite sides of a snp" ) {
 
            id_t seed_nodes[] = {2, 3, 5, 8};
            //all are in the same cluster
            for (bool use_minimizers : {true, false} ) {
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (id_t n : seed_nodes) {
                    pos_t pos = make_pos_t(n, false, 0);
                    auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 4); 
                REQUIRE(clusters.size() == 4); 
            }


        }
        SECTION( "Three clusters on opposite sides of a snp in two reads" ) {
 
            id_t seed_nodes[] = {2, 4};
            //all are in the same cluster
            vector<vector<SnarlDistanceIndexClusterer::Seed>> seeds (2);

            pos_t pos = make_pos_t(2, false, 0);
            seeds[0].push_back({ pos, 0});
            pos = make_pos_t(3, false, 0);
            seeds[0].push_back({ pos, 0});

            pos = make_pos_t(5, false, 0);
            seeds[1].push_back({ pos, 0});

            vector<vector<SnarlDistanceIndexClusterer::Cluster>> clusters = clusterer.cluster_seeds(seeds, 5, 5); 
            REQUIRE(clusters.size() == 2); 
            REQUIRE(clusters[0][0].fragment == clusters[1][0].fragment);
            REQUIRE(clusters[0].size() == 2);
            REQUIRE(clusters[1].size() == 1);
        }
        SECTION( "Three clusters on opposite sides of a snp in two reads going in the other direction in the chain" ) {
 
            //all are in the same cluster
            vector<vector<SnarlDistanceIndexClusterer::Seed>> seeds (2);

            pos_t pos = make_pos_t(5, false, 0);
            seeds[0].push_back({ pos, 0});
            pos = make_pos_t(6, false, 0);
            seeds[0].push_back({ pos, 0});

            pos = make_pos_t(1, false, 0);
            seeds[1].push_back({ pos, 0});

            vector<vector<SnarlDistanceIndexClusterer::Cluster>> clusters = clusterer.cluster_seeds(seeds, 10, 10); 
            REQUIRE(clusters.size() == 2); 
            REQUIRE(clusters[0][0].fragment == clusters[1][0].fragment);
            REQUIRE(clusters[0].size() == 2);
            REQUIRE(clusters[1].size() == 1);
        }
    }

    TEST_CASE( "cluster long snarl in chain",
                   "[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("GGC");
        Node* n2 = graph.create_node("GCA");
        Node* n3 = graph.create_node("GCAGCACATGCACATC"); //16
        Node* n4 = graph.create_node("GCA");
        Node* n5 = graph.create_node("GCAAGCACATGCACATCCA");
        Node* n6 = graph.create_node("GCA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n2, n3);
        Edge* e3 = graph.create_edge(n2, n4);
        Edge* e4 = graph.create_edge(n3, n5);
        Edge* e5 = graph.create_edge(n4, n5);
        Edge* e6 = graph.create_edge(n1, n6);
        Edge* e7 = graph.create_edge(n6, n2);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);
        SnarlDistanceIndexClusterer clusterer(dist_index, &graph);
        
        //graph.to_dot(cerr);

        SECTION( "Two clusters including snarl" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(make_pos_t(2, true, 0));
            positions.emplace_back(make_pos_t(3, false, 8));
            positions.emplace_back(make_pos_t(5, false, 0));
            //all are in the same cluster
            for (bool use_minimizers : {true, false} ) {
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (pos_t pos : positions) {
                    auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 5); 
                REQUIRE(clusters.size() == 2); 
            }


        }
        SECTION( "Three clusters not including snarl" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(make_pos_t(2, true, 0));
            positions.emplace_back(make_pos_t(3, false, 8));
            positions.emplace_back(make_pos_t(5, false, 0));
            //all are in the same cluster
            for (bool use_minimizers : {true, false} ) {
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (pos_t pos : positions) {
                    auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 2); 
                REQUIRE(clusters.size() == 3); 
            }


        }
    }

    TEST_CASE("Use path through big snarl", "[cluster]") {
        //Chain: 1 - (snarl 2-7) - 8
        
        VG graph;
        
        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("C");
        Node* n3 = graph.create_node("A");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("AGTA");
        Node* n9 = graph.create_node("AGTAAGTA");
        Node* n10 = graph.create_node("A");
        Node* n11 = graph.create_node("AGTAAAA");
        Node* n12 = graph.create_node("AG");
        Node* n13 = graph.create_node("AGT");
        Node* n14 = graph.create_node("AG");
        Node* n15 = graph.create_node("AGTA");
        
        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e3 = graph.create_edge(n2, n4);
        Edge* e4 = graph.create_edge(n3, n4);
        Edge* e5 = graph.create_edge(n4, n5);
        Edge* e6 = graph.create_edge(n4, n6);
        Edge* e7 = graph.create_edge(n5, n6);
        Edge* e8 = graph.create_edge(n6, n2, false, true);
        Edge* e9 = graph.create_edge(n6, n7);
        Edge* e10 = graph.create_edge(n7, n8);
        Edge* e11 = graph.create_edge(n4, n9);
        Edge* e12 = graph.create_edge(n9, n7);
        Edge* e13 = graph.create_edge(n8, n11);
        Edge* e14 = graph.create_edge(n8, n10);
        Edge* e15 = graph.create_edge(n10, n12);
        Edge* e16 = graph.create_edge(n10, n13);
        Edge* e17 = graph.create_edge(n11, n12);
        Edge* e18 = graph.create_edge(n11, n15);
        Edge* e19 = graph.create_edge(n12, n14);
        Edge* e20 = graph.create_edge(n14, n15);
        Edge* e21 = graph.create_edge(n11, n14);
        
        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);
        SnarlDistanceIndexClusterer clusterer(distance_index, &graph);
        SECTION("one cluster in same snarl") {
            vector<pos_t> positions;
            positions.emplace_back(make_pos_t(10, false, 0));
            positions.emplace_back(make_pos_t(12, false, 1));
            //all are in the same cluster
            for (bool use_minimizers : {false, true} ) {
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (pos_t pos : positions) {
                    auto chain_info = MIPayload::encode(get_minimizer_distances(distance_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 2); 
                REQUIRE(clusters.size() == 1); 
            }
        }
        SECTION("two clusters in same snarl") {
            vector<pos_t> positions;
            positions.emplace_back(make_pos_t(10, false, 0));
            positions.emplace_back(make_pos_t(12, false, 1));
            //all are in the same cluster
            for (bool use_minimizers : {false, true} ) {
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (pos_t pos : positions) {
                    auto chain_info = MIPayload::encode(get_minimizer_distances(distance_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 1); 
                REQUIRE(clusters.size() == 2); 
            }
        }
        SECTION("one cluster in same snarl separated by one node") {
            vector<pos_t> positions;
            positions.emplace_back(make_pos_t(10, false, 0));
            positions.emplace_back(make_pos_t(14, false, 0));
            //all are in the same cluster
            for (bool use_minimizers : {false, true} ) {
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (pos_t pos : positions) {
                    auto chain_info = MIPayload::encode(get_minimizer_distances(distance_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 3); 
                REQUIRE(clusters.size() == 1); 
            }
        }
        SECTION("two clusters in same snarl separated by one node") {
            vector<pos_t> positions;
            positions.emplace_back(make_pos_t(10, false, 0));
            positions.emplace_back(make_pos_t(14, false, 0));
            //all are in the same cluster
            for (bool use_minimizers : {false, true} ) {
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (pos_t pos : positions) {
                    auto chain_info = MIPayload::encode(get_minimizer_distances(distance_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 2); 
                REQUIRE(clusters.size() == 2); 
            }
        }
        SECTION("two clusters using path in different snarl") {
            vector<pos_t> positions;
            positions.emplace_back(make_pos_t(5, false, 0));
            positions.emplace_back(make_pos_t(12, false, 0));
            //all are in the same cluster
            for (bool use_minimizers : {false, true} ) {
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (pos_t pos : positions) {
                    auto chain_info = MIPayload::encode(get_minimizer_distances(distance_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 9); 
                REQUIRE(clusters.size() == 2); 
            }
        }
        SECTION("one cluster using path in different snarl") {
            vector<pos_t> positions;
            positions.emplace_back(make_pos_t(5, false, 0));
            positions.emplace_back(make_pos_t(12, false, 0));
            //all are in the same cluster
            for (bool use_minimizers : {false, true} ) {
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (pos_t pos : positions) {
                    auto chain_info = MIPayload::encode(get_minimizer_distances(distance_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 10); 
                REQUIRE(clusters.size() == 1); 
            }
        }
        SECTION("one cluster") {
            vector<pos_t> positions;
            positions.emplace_back(make_pos_t(2, false, 0));
            positions.emplace_back(make_pos_t(4, false, 0));
            positions.emplace_back(make_pos_t(9, true, 2));
            positions.emplace_back(make_pos_t(7, false, 0));
            //all are in the same cluster
            for (bool use_minimizers : {true, false} ) {
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (pos_t pos : positions) {
                    auto chain_info = MIPayload::encode(get_minimizer_distances(distance_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 8); 
                REQUIRE(clusters.size() == 1); 
            }
        }
        SECTION("two clusters") {
            vector<pos_t> positions;
            positions.emplace_back(make_pos_t(12, false, 0));
            positions.emplace_back(make_pos_t(7, false, 0));
            //all are in the same cluster
            for (bool use_minimizers : {true, false} ) {
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (pos_t pos : positions) {
                    auto chain_info = MIPayload::encode(get_minimizer_distances(distance_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 4); 
                REQUIRE(clusters.size() == 2); 
            }
        }
        SECTION("distance") {
            vector<pos_t> positions;
            positions.emplace_back(make_pos_t(12, false, 0));
            positions.emplace_back(make_pos_t(7, false, 0));
            positions.emplace_back(make_pos_t(1, false, 0));
            positions.emplace_back(make_pos_t(5, false, 0));
            //all are in the same cluster
            for (bool use_minimizers : {true, false} ) {
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (pos_t pos : positions) {
                    auto chain_info = MIPayload::encode(get_minimizer_distances(distance_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                REQUIRE(clusterer.distance_between_seeds(seeds[0], seeds[1], false) == 6); 
                REQUIRE(clusterer.distance_between_seeds(seeds[3], seeds[2], false) == 7); 
                REQUIRE(clusterer.distance_between_seeds(seeds[2], seeds[3], false) == 7); 
            }
        }
    }

    TEST_CASE( "Weird loop with three components of the root",
                   "[cluster]" ) {
        //THis is a symmetrical graph with two weird loopy things on the ends of a chain from 4 to 15
        VG graph;

        Node* n1 = graph.create_node("G");
        Node* n2 = graph.create_node("G");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("G");
        Node* n5 = graph.create_node("G");
        Node* n6 = graph.create_node("G");
        Node* n7 = graph.create_node("AACAT"); //5
        Node* n8 = graph.create_node("GACAT");
        Node* n9 = graph.create_node("CACAT");
        Node* n10 = graph.create_node("CACAT");
        Node* n11 = graph.create_node("A");
        Node* n12 = graph.create_node("A");
        Node* n13 = graph.create_node("A");
        Node* n14 = graph.create_node("A");
        Node* n15 = graph.create_node("C");
        Node* n16 = graph.create_node("G");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n2, n1);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n3, n1);
        Edge* e5 = graph.create_edge(n1, n4);
        Edge* e6 = graph.create_edge(n4, n5);
        Edge* e7 = graph.create_edge(n4, n6);
        Edge* e8 = graph.create_edge(n5, n6);
        Edge* e9 = graph.create_edge(n6, n7);
        Edge* e10 = graph.create_edge(n6, n7);
        Edge* e11 = graph.create_edge(n6, n10);
        Edge* e26 = graph.create_edge(n7, n10);
        Edge* e12 = graph.create_edge(n7, n8);
        Edge* e13 = graph.create_edge(n7, n9);
        Edge* e14 = graph.create_edge(n8, n9);
        Edge* e15 = graph.create_edge(n9, n11);
        Edge* e16 = graph.create_edge(n10, n9);
        Edge* e17 = graph.create_edge(n10, n11);
        Edge* e18 = graph.create_edge(n11, n12);
        Edge* e19 = graph.create_edge(n11, n13);
        Edge* e20 = graph.create_edge(n12, n13);
        Edge* e21 = graph.create_edge(n13, n14);
        Edge* e22 = graph.create_edge(n14, n15);
        Edge* e23 = graph.create_edge(n14, n16);
        Edge* e24 = graph.create_edge(n16, n15);
        Edge* e25 = graph.create_edge(n15, n14);
        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);
        SnarlDistanceIndexClusterer clusterer(dist_index, &graph);
        
        //graph.to_dot(cerr);

        SECTION( "Three clusters going across snarl" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(make_pos_t(2, false, 0));
            positions.emplace_back(make_pos_t(11, false, 0));
            positions.emplace_back(make_pos_t(8, false, 2));
            //all are in the same cluster
            for (bool use_minimizers : {true, false} ) {
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (pos_t pos : positions) {
                    auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 4); 
                REQUIRE(clusters.size() == 3); 
            }


        }
        SECTION( "A bunch of nodes in the snarl" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(make_pos_t(6, true, 0));
            positions.emplace_back(make_pos_t(8, false, 0));
            positions.emplace_back(make_pos_t(8, false, 2));
            positions.emplace_back(make_pos_t(10, false, 0));
            positions.emplace_back(make_pos_t(10, false, 2));
            positions.emplace_back(make_pos_t(8, false, 2));
            positions.emplace_back(make_pos_t(7, false, 2));
            positions.emplace_back(make_pos_t(9, false, 0));
            positions.emplace_back(make_pos_t(13, false, 0));
            positions.emplace_back(make_pos_t(7, false, 0));
            //all are in the same cluster
            for (bool use_minimizers : {true, false} ) {
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (pos_t pos : positions) {
                    auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 3); 
                REQUIRE(clusters.size() == 2); 
            }
        }
        SECTION( "A bunch of nodes in the snarl on the other side" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(make_pos_t(6, true, 0));
            positions.emplace_back(make_pos_t(9, false, 0));
            positions.emplace_back(make_pos_t(9, false, 2));
            positions.emplace_back(make_pos_t(8, false, 0));
            positions.emplace_back(make_pos_t(8, false, 2));
            positions.emplace_back(make_pos_t(8, false, 2));
            positions.emplace_back(make_pos_t(10, false, 2));
            positions.emplace_back(make_pos_t(13, false, 0));
            //all are in the same cluster
            for (bool use_minimizers : {true, false} ) {
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (pos_t pos : positions) {
                    auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 3); 
                REQUIRE(clusters.size() == 2); 
            }
        }
    }
    TEST_CASE( "Cluster looping, multicomponent",
                   "[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("TGCAGT");
        Node* n7 = graph.create_node("T");
        Node* n8 = graph.create_node("CTGA");
        Node* n9 = graph.create_node("GCA");
        Node* n10 = graph.create_node("T");
        Node* n11 = graph.create_node("T");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n10);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n4);
        Edge* e5 = graph.create_edge(n3, n4);
        Edge* e6 = graph.create_edge(n4, n5);
        Edge* e7 = graph.create_edge(n4, n6);
        Edge* e8 = graph.create_edge(n5, n6);
        Edge* e9 = graph.create_edge(n6, n7);
        Edge* e10 = graph.create_edge(n6, n8);
        Edge* e11 = graph.create_edge(n7, n8);
        Edge* e12 = graph.create_edge(n8, n9);
        Edge* e13 = graph.create_edge(n8, n10);
        Edge* e14 = graph.create_edge(n9, n11);
        Edge* e15 = graph.create_edge(n10, n11);
        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);
        SnarlDistanceIndexClusterer clusterer(dist_index, &graph);

    
        
        //graph.to_dot(cerr);

        SECTION( "Test distance values" ) {
            net_handle_t node1 = dist_index.get_parent(dist_index.get_node_net_handle(n1->id()));
            net_handle_t snarl82 = dist_index.get_parent(node1); 

            if (dist_index.node_id(dist_index.get_bound(snarl82, false, false)) == n2->id()) {
                //If the snarl is from 2rev to 8rev
                REQUIRE(dist_index.distance_to_parent_bound(snarl82, false, node1) == std::numeric_limits<size_t>::max());
                REQUIRE(dist_index.distance_to_parent_bound(snarl82, false, dist_index.flip(node1)) == std::numeric_limits<size_t>::max());
                REQUIRE(dist_index.distance_to_parent_bound(snarl82, true, node1) == 0);
                REQUIRE(dist_index.distance_to_parent_bound(snarl82, true, dist_index.flip(node1)) == std::numeric_limits<size_t>::max());
            } else {
                REQUIRE(dist_index.distance_to_parent_bound(snarl82, false, node1) == 0);
                REQUIRE(dist_index.distance_to_parent_bound(snarl82, false, dist_index.flip(node1)) == std::numeric_limits<size_t>::max());
                REQUIRE(dist_index.distance_to_parent_bound(snarl82, true, node1) == std::numeric_limits<size_t>::max());
                REQUIRE(dist_index.distance_to_parent_bound(snarl82, true, dist_index.flip(node1)) == std::numeric_limits<size_t>::max());
            }


            net_handle_t node3 = dist_index.get_parent(dist_index.get_node_net_handle(n3->id()));
            net_handle_t snarl24 = dist_index.get_parent(node3); 

            if (dist_index.node_id(dist_index.get_bound(snarl24, false, false)) == n2->id()) {
                REQUIRE(dist_index.distance_to_parent_bound(snarl24, true, dist_index.flip(node3)) == 0);
                REQUIRE(dist_index.distance_to_parent_bound(snarl24, true, node3) == std::numeric_limits<size_t>::max());
                REQUIRE(dist_index.distance_to_parent_bound(snarl24, false, node3) == 0);
                REQUIRE(dist_index.distance_to_parent_bound(snarl24, false, dist_index.flip(node3)) == std::numeric_limits<size_t>::max());
            } else {
                REQUIRE(dist_index.distance_to_parent_bound(snarl24, true, dist_index.flip(node3)) == std::numeric_limits<size_t>::max());
                REQUIRE(dist_index.distance_to_parent_bound(snarl24, true, node3) == 0);
                REQUIRE(dist_index.distance_to_parent_bound(snarl24, false, node3) == std::numeric_limits<size_t>::max());
                REQUIRE(dist_index.distance_to_parent_bound(snarl24, false, dist_index.flip(node3)) == 0);
            }

            net_handle_t node6 = dist_index.get_node_net_handle(n6->id());
            net_handle_t chain66 = dist_index.get_parent(node6);
            net_handle_t node5 = dist_index.get_parent(dist_index.get_node_net_handle(n5->id()));
            net_handle_t snarl46 = dist_index.get_parent(node5);
            if (dist_index.node_id(dist_index.get_bound(snarl46, false, false)) == n6->id()) {
                snarl46 = dist_index.flip(snarl46);
            }
            REQUIRE(dist_index.distance_in_parent(chain66, snarl46, dist_index.flip(node6)) == 0);
        }
        SECTION( "Two clusters" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(make_pos_t(1, false, 0));
            positions.emplace_back(make_pos_t(3, false, 0));
            positions.emplace_back(make_pos_t(9, false, 0));
            positions.emplace_back(make_pos_t(10, false, 0));
            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (bool use_minimizers : {true, false} ) {
                for (pos_t pos : positions) {
                    auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 5); 
                REQUIRE(clusters.size() == 2); 
            }


        }
        SECTION( "Two clusters" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(make_pos_t(2, false, 0));
            positions.emplace_back(make_pos_t(8, false, 0));
            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (bool use_minimizers : {true, false} ) {
                for (pos_t pos : positions) {
                    auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 5); 
                REQUIRE(clusters.size() == 2); 
            }


        }
        SECTION( "One cluster" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(make_pos_t(5, false, 0));
            positions.emplace_back(make_pos_t(7, false, 0));
            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (bool use_minimizers : {true, false} ) {
                seeds.clear();
                for (pos_t pos : positions) {
                    auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 9); 
                REQUIRE(clusters.size() == 1); 
            }
        }
        SECTION( "Two clusters" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(make_pos_t(3, false, 0));
            positions.emplace_back(make_pos_t(7, false, 0));
            positions.emplace_back(make_pos_t(11, false, 0));
            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (bool use_minimizers : {true, false} ) {
                seeds.clear();
                for (pos_t pos : positions) {
                    auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 10); 
                REQUIRE(clusters.size() == 2); 
            }
        }
        SECTION( "One cluster" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(make_pos_t(3, false, 0));
            positions.emplace_back(make_pos_t(7, false, 0));
            positions.emplace_back(make_pos_t(11, false, 0));
            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (bool use_minimizers : {true, false} ) {
                seeds.clear();
                for (pos_t pos : positions) {
                    auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 11); 
                REQUIRE(clusters.size() == 1); 
            }
        }

    }
    TEST_CASE( "looping chain of nested unary snarls",
                   "[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("T");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n2, n4);
        Edge* e4 = graph.create_edge(n3, n4);
        Edge* e5 = graph.create_edge(n4, n5);
        Edge* e6 = graph.create_edge(n4, n6);
        Edge* e7 = graph.create_edge(n5, n6);
        Edge* e8 = graph.create_edge(n6, n6, false, true);
        Edge* e9 = graph.create_edge(n1, n1, true, false);


        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);
        SnarlDistanceIndexClusterer clusterer(dist_index, &graph);
        
        //graph.to_dot(cerr);

        SECTION( "One cluster taking loop" ) {
 
            for (bool use_minimizers : {true, false} ) {
                id_t seed_nodes[] = {1, 4};
                //all are in the same cluster
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (id_t n : seed_nodes) {
                    pos_t pos = make_pos_t(n, false, 0);
                    auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }

                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 6); 
                REQUIRE(clusters.size() == 1); 
            }

        }
        SECTION( "One cluster on boundary" ) {
 
            for (bool use_minimizers : {true, false} ) {
                id_t seed_nodes[] = {2, 4};
                //all are in the same cluster
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (id_t n : seed_nodes) {
                    pos_t pos = make_pos_t(n, false, 0);
                    auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }

                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 3); 
                REQUIRE(clusters.size() == 1); 
            }

        }
        SECTION( "One fragment cluster on boundary" ) {
 
            id_t seed_nodes[] = {2, 4};
            //all are in the same cluster
            vector<vector<SnarlDistanceIndexClusterer::Seed>> seeds (2);

            pos_t pos = make_pos_t(2, false, 0);
            seeds[0].push_back({ pos, 0});

            pos = make_pos_t(4, false, 0);
            seeds[1].push_back({ pos, 0});

            vector<vector<SnarlDistanceIndexClusterer::Cluster>> clusters = clusterer.cluster_seeds(seeds, 3, 3); 
            REQUIRE(clusters.size() == 2); 
            REQUIRE(clusters[0][0].fragment == clusters[1][0].fragment);

        }
        SECTION( "One cluster on boundary" ) {
 
            for (bool use_minimizers : {true, false} ) {
                id_t seed_nodes[] = {3, 4};
                //all are in the same cluster
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (id_t n : seed_nodes) {
                    pos_t pos = make_pos_t(n, false, 0);
                    auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }

                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 3); 
                REQUIRE(clusters.size() == 1); 

            }
        }
    }
    TEST_CASE( "chain with loop",
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

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n4);
        Edge* e5 = graph.create_edge(n3, n4);
        Edge* e6 = graph.create_edge(n3, n5);
        Edge* e7 = graph.create_edge(n4, n5);
        Edge* e8 = graph.create_edge(n4, n6);
        Edge* e9 = graph.create_edge(n5, n6);
        Edge* e10 = graph.create_edge(n6, n7);
        Edge* e11 = graph.create_edge(n6, n7, false, true);
        Edge* e12 = graph.create_edge(n6, n8);
        Edge* e13 = graph.create_edge(n7, n8);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);
        SnarlDistanceIndexClusterer clusterer(dist_index, &graph);
        
        //graph.to_dot(cerr);

        SECTION( "One cluster taking loop" ) {
 
            for (bool use_minimizers : {true, false} ) {
                id_t seed_nodes[] = {4, 5};
                //all are in the same cluster
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (id_t n : seed_nodes) {
                    pos_t pos = make_pos_t(n, false, 0);
                    auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }

                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 11); 
                REQUIRE(clusters.size() == 1); 
            }

        }
        SECTION( "One cluster not taking loop" ) {
 
            for (bool use_minimizers : {true, false} ) {
                id_t seed_nodes[] = {4, 5, 3};
                //all are in the same cluster
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (id_t n : seed_nodes) {
                    pos_t pos = make_pos_t(n, false, 0);
                    auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }

                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 3); 
                REQUIRE(clusters.size() == 1); 
            }
        }
        SECTION( "One cluster not taking loop" ) {
 
            for (bool use_minimizers : {true, false} ) {
                id_t seed_nodes[] = {4, 5, 6};
                //all are in the same cluster
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (id_t n : seed_nodes) {
                    pos_t pos = make_pos_t(n, false, 0);
                    auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }

                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 8); 
                REQUIRE(clusters.size() == 1); 
            }

        }
        SECTION( "Two clusters" ) {
 
            for (bool use_minimizers : {true, false} ) {
                id_t seed_nodes[] = {4, 5, 1};
                //all are in the same cluster
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (id_t n : seed_nodes) {
                    pos_t pos = make_pos_t(n, false, 0);
                    auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }

                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 3); 
                REQUIRE(clusters.size() == 3); 
            }

        }
    }
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

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);
        SnarlDistanceIndexClusterer clusterer(dist_index, &graph);
        
        //graph.to_dot(cerr);

        SECTION( "One cluster with seed struct" ) {
 
            for (bool use_minimizers : {true, false} ) {
                id_t seed_nodes[] = {2, 3, 4, 7, 8, 9, 11};
                //all are in the same cluster
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (id_t n : seed_nodes) {
                    pos_t pos = make_pos_t(n, false, 0);
                    auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                    if (use_minimizers) {
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }

                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 10); 
                REQUIRE(clusters.size() == 1); 
            }
        }
        SECTION( "Two clusters" ) {
            for (bool use_minimizers : {true, false} ) {
 
                vector<id_t> seed_nodes( {2, 3, 4, 7, 8, 10, 11});
                //Clusters should be {2, 3, 4}, {7, 8, 10, 11}
                //Distance from pos on 4 to pos on 7 is 8, including one position
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (id_t n : seed_nodes) {
                    pos_t pos = make_pos_t(n, false, 0);
                    if (use_minimizers) {
                        auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }


                vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 7); 
                vector<hash_set<size_t>> cluster_sets;
                for (auto& c : clusters) {
                    hash_set<size_t> h;
                    for (size_t s : c.seeds) {
                        h.insert(s);
                    }
                    cluster_sets.push_back(h);
                }
                REQUIRE( clusters.size() == 2);
                REQUIRE (( (cluster_sets[0].count(0) == 1 &&
                           cluster_sets[0].count(1) == 1 &&
                           cluster_sets[0].count(2) == 1 &&
                           cluster_sets[1].count(3) == 1 &&
                           cluster_sets[1].count(4) == 1 &&
                           cluster_sets[1].count(5) == 1 &&
                           cluster_sets[1].count(6) == 1  ) ||

                         ( cluster_sets[1].count(0) == 1 &&
                           cluster_sets[1].count(1) == 1 &&
                           cluster_sets[1].count(2) == 1 &&
                           cluster_sets[0].count(3) == 1 &&
                           cluster_sets[0].count(4) == 1 &&
                           cluster_sets[0].count(5) == 1 &&
                           cluster_sets[0].count(6) == 1  )));

            }
        }
        SECTION( "One fragment cluster of the same node" ) {
 
            vector<id_t> seed_nodes( {2, 3});
            vector<id_t> seed_nodes1({2, 7, 8, 10, 11});
            //Clusters should be {2, 3, 4}, {2}, {7, 8, 10, 11}
            //One fragment cluster
            //Distance from pos on 4 to pos on 7 is 8, including one position
            //
            for (bool use_minimizers : {true, false} ) {
                vector<SnarlDistanceIndexClusterer::Seed> seeds ;
                for (id_t n : seed_nodes) {
                    pos_t pos = make_pos_t(n, false, 0);
                    if (use_minimizers) {
                        auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                vector<SnarlDistanceIndexClusterer::Seed> seeds1;
                for (id_t n : seed_nodes1) {
                    pos_t pos = make_pos_t(n, false, 0);
                    if (use_minimizers) {
                        auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                        seeds1.push_back({ pos, 0, chain_info});
                    } else {
                        seeds1.push_back({ pos, 0});
                    }
                }
                vector<vector<SnarlDistanceIndexClusterer::Seed>> all_seeds;
                all_seeds.push_back(seeds);
                all_seeds.push_back(seeds1);


                vector<vector<SnarlDistanceIndexClusterer::Cluster>> paired_clusters = clusterer.cluster_seeds(all_seeds, 7, 15); 
                //Should be [[<[0,1,2], 0>],[<[3,4,5,6], 0>]] 
                REQUIRE( paired_clusters.size() == 2);
                REQUIRE( paired_clusters[0].size() == 1);
                REQUIRE( paired_clusters[1].size() == 2);
                REQUIRE( paired_clusters[0][0].fragment == paired_clusters[1][0].fragment);
                REQUIRE( paired_clusters[1][0].fragment == paired_clusters[1][1].fragment);
            }
        }
        SECTION( "One fragment cluster" ) {
            for (bool use_minimizers : {true, false}) {
 
                vector<id_t> seed_nodes( {2, 3, 4});
                vector<id_t> seed_nodes1({7, 8, 10, 11});
                //Clusters should be {2, 3, 4}, {7, 8, 10, 11}
                //One fragment cluster
                //Distance from pos on 4 to pos on 7 is 8, including one position
                vector<SnarlDistanceIndexClusterer::Seed> seeds ;
                for (id_t n : seed_nodes) {
                    pos_t pos = make_pos_t(n, false, 0);
                    if (use_minimizers) {
                        auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                vector<SnarlDistanceIndexClusterer::Seed> seeds1;
                for (id_t n : seed_nodes1) {
                    pos_t pos = make_pos_t(n, false, 0);
                    if (use_minimizers) {
                        auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                        seeds1.push_back({ pos, 0, chain_info});
                    } else {
                        seeds1.push_back({ pos, 0});
                    }
                }
                vector<vector<SnarlDistanceIndexClusterer::Seed>> all_seeds;
                all_seeds.push_back(seeds);
                all_seeds.push_back(seeds1);


                vector<vector<SnarlDistanceIndexClusterer::Cluster>> paired_clusters = clusterer.cluster_seeds(all_seeds, 7, 15); 
                //Should be [[<[0,1,2], 0>],[<[3,4,5,6], 0>]] 
                REQUIRE( paired_clusters.size() == 2);
                REQUIRE( paired_clusters[0].size() == 1);
                REQUIRE( paired_clusters[1].size() == 1);
                REQUIRE( paired_clusters[0][0].seeds.size() == 3);
                REQUIRE( paired_clusters[1][0].seeds.size() == 4);
                REQUIRE( paired_clusters[0][0].fragment == paired_clusters[1][0].fragment);
            }
        }
        SECTION( "Two fragment clusters with seed structs" ) {
 
            vector<id_t> seed_nodes( {2, 3, 4});
            vector<id_t> seed_nodes1({7, 8, 10, 11});
            //Fragment clusters should be {2, 3, 4}, {7, 8, 10, 11}
            //Distance from pos on 4 to pos on 7 is 8, including one position
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : seed_nodes) {
                pos_t pos = make_pos_t(n, false, 0);
                auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                seeds.push_back({ pos, 0, chain_info});
            }
            vector<SnarlDistanceIndexClusterer::Seed> seeds1;
            for (id_t n : seed_nodes1) {
                pos_t pos = make_pos_t(n, false, 0);
                auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                seeds1.push_back({ pos, 0, chain_info});
            }
            vector<vector<SnarlDistanceIndexClusterer::Seed>> all_seeds;
            all_seeds.push_back(seeds);
            all_seeds.push_back(seeds1);


            vector<vector<SnarlDistanceIndexClusterer::Cluster>> paired_clusters = clusterer.cluster_seeds(all_seeds, 2, 7); 
            // read_clusters = [ [[0,1,2]],[[3,4],[5,6]] ]
            // fragment_clusters = [ [0,1,2], [3,4,5,6] ]
            REQUIRE( paired_clusters.size() == 2) ;
            REQUIRE( paired_clusters[0].size() == 1);
            REQUIRE( paired_clusters[1].size() == 2);
            REQUIRE( paired_clusters[0][0].fragment != paired_clusters[1][0].fragment);
            REQUIRE( paired_clusters[0][0].fragment != paired_clusters[1][1].fragment);
            REQUIRE( paired_clusters[1][0].fragment == paired_clusters[1][1].fragment);

        }
        SECTION( "Two fragment clusters" ) {
 
            vector<id_t> seed_nodes( {2, 3, 4});
            vector<id_t> seed_nodes1({7, 8, 10, 11});
            //Fragment clusters should be {2, 3, 4}, {7, 8, 10, 11}
            //Distance from pos on 4 to pos on 7 is 8, including one position
            vector<SnarlDistanceIndexClusterer::Seed> seeds ;
            for (id_t n : seed_nodes) {
                pos_t pos = make_pos_t(n, false, 0);
                auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                seeds.push_back({ pos, 0, chain_info});
            }
            vector<SnarlDistanceIndexClusterer::Seed> seeds1;
            for (id_t n : seed_nodes1) {
                pos_t pos = make_pos_t(n, false, 0);
                auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                seeds1.push_back({ pos, 0, chain_info});
            }
            vector<vector<SnarlDistanceIndexClusterer::Seed>> all_seeds;
            all_seeds.push_back(seeds);
            all_seeds.push_back(seeds1);


            vector<vector<SnarlDistanceIndexClusterer::Cluster>> paired_clusters = clusterer.cluster_seeds(all_seeds, 2, 7); 
            // read_clusters = [ [[0,1,2]],[[3,4],[5,6]] ]
            // fragment_clusters = [ [0,1,2], [3,4,5,6] ]
            REQUIRE( paired_clusters.size() == 2) ;
            REQUIRE( paired_clusters[0].size() == 1);
            REQUIRE( paired_clusters[1].size() == 2);
            REQUIRE( paired_clusters[0][0].fragment != paired_clusters[1][0].fragment);
            REQUIRE( paired_clusters[0][0].fragment != paired_clusters[1][1].fragment);
            REQUIRE( paired_clusters[1][0].fragment == paired_clusters[1][1].fragment);

        }
    }//End test case

    TEST_CASE( "Reverse in chain right","[cluster]" ) {
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
        Node* n11 = graph.create_node("GGGGGGGGGG");//10

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n10);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n4);
        Edge* e5 = graph.create_edge(n3, n5);
        Edge* e6 = graph.create_edge(n4, n5);
        Edge* e7 = graph.create_edge(n5, n6);
        Edge* e8 = graph.create_edge(n5, n11);
        Edge* e9 = graph.create_edge(n11, n7);
        Edge* e10 = graph.create_edge(n6, n7);
        Edge* e11 = graph.create_edge(n8, n8, false, true);
        Edge* e12 = graph.create_edge(n7, n8);
        Edge* e13 = graph.create_edge(n7, n9);
        Edge* e14 = graph.create_edge(n8, n9);
        Edge* e15 = graph.create_edge(n9, n10);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);


        SnarlDistanceIndexClusterer clusterer(dist_index, &graph);

        SECTION( "Same snarl" ) {
            vector<id_t> seed_nodes ({3, 4});
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : seed_nodes) {
                pos_t pos = make_pos_t(n, false, 0);
                auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                seeds.push_back({ pos, 0, chain_info});
            }


            vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 13); 

            REQUIRE( clusters.size() == 1);
        }
        SECTION( "Different snarl" ) {
            vector<SnarlDistanceIndexClusterer::Seed> seeds;

            vector<pos_t> pos_ts;
            pos_ts.emplace_back(3, false, 0);
            pos_ts.emplace_back(11, false, 9);
            for (pos_t pos : pos_ts) {
                seeds.push_back({ pos, 0});
            }



            vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 8); 


            REQUIRE( clusters.size() == 1);
        }
    }//end test case
    TEST_CASE( "Reverse in chain left","[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("G");
        Node* n6 = graph.create_node("TGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("G");
        Node* n9 = graph.create_node("AA");
        Node* n10 = graph.create_node("G");
        Node* n11 = graph.create_node("G");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n10);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n4);
        Edge* e5 = graph.create_edge(n3, n5);
        Edge* e6 = graph.create_edge(n4, n5);
        Edge* e7 = graph.create_edge(n5, n6);
        Edge* e8 = graph.create_edge(n5, n7);
        Edge* e9 = graph.create_edge(n6, n7);
        Edge* e10 = graph.create_edge(n7, n8);
        Edge* e11 = graph.create_edge(n7, n9);
        Edge* e12 = graph.create_edge(n8, n9);
        Edge* e13 = graph.create_edge(n9, n10);
        Edge* e14 = graph.create_edge(n11, n5);
        Edge* e15 = graph.create_edge(n11, n5, true, false);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);

        SnarlDistanceIndexClusterer clusterer(dist_index, &graph);

        SECTION( "One cluster" ) {
            vector<id_t> seed_nodes ({7, 7, 6});
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : seed_nodes) {
                pos_t pos = make_pos_t(n, false, 0);
                seeds.push_back({ pos, 0});
            }

            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 20); 


            REQUIRE( clusters.size() == 1);
        }
        SECTION( "two clusters" ) {
            vector<id_t> seed_nodes ({2, 6});
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : seed_nodes) {
                pos_t pos = make_pos_t(n, false, 0);
                auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                seeds.push_back({ pos, 0, chain_info});
            }

            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds,  20); 


        }
        SECTION( "different snarl" ) {
            vector<id_t> seed_nodes ({8, 6});
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : seed_nodes) {
                pos_t pos = make_pos_t(n, false, 0);
                auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                seeds.push_back({ pos, 0, chain_info});
            }

            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 20); 


            REQUIRE( clusters.size() == 1);
        }
    }//end test case


    TEST_CASE( "Loop on node","[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GGGGGGGGGGGG");//12 Gs
        Node* n6 = graph.create_node("T");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n2, n4);
        Edge* e4 = graph.create_edge(n3, n4);
        Edge* e27 = graph.create_edge(n4, n5);
        Edge* e5 = graph.create_edge(n4, n6);
        Edge* e6 = graph.create_edge(n5, n6);
        Edge* e7 = graph.create_edge(n5, n5);

        
        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);

        SnarlDistanceIndexClusterer clusterer(dist_index, &graph);

        SECTION( "One cluster taking node loop" ) {
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            vector<pos_t> pos_ts;
            pos_ts.emplace_back(5, false, 0);
            pos_ts.emplace_back(5, true, 0);

            for (pos_t pos : pos_ts){

                auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                seeds.push_back({ pos, 0, chain_info});
            }
            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 3); 

            REQUIRE( clusters.size() == 1);
        }
    }
    TEST_CASE( "Loop on first node in a top-level chain","[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GGGGGGGGGGGG");//12 Gs
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("CTGA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n1);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n2);
        Edge* e5 = graph.create_edge(n3, n4);
        Edge* e6 = graph.create_edge(n3, n5);
        Edge* e7 = graph.create_edge(n4, n5);
        Edge* e8 = graph.create_edge(n5, n6);
        Edge* e9 = graph.create_edge(n5, n7);
        Edge* e10 = graph.create_edge(n6, n7);
        Edge* e11 = graph.create_edge(n7, n8);

        
        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);

        SnarlDistanceIndexClusterer clusterer(dist_index, &graph);

        SECTION( "One cluster across top-level snarl" ) {
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            vector<pos_t> pos_ts;
            pos_ts.emplace_back(1, false, 0);
            pos_ts.emplace_back(4, true, 0);

            for (pos_t pos : pos_ts){

                auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                seeds.push_back({ pos, 0, chain_info});
            }
            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 10); 

            REQUIRE( clusters.size() == 1);
        }
        SECTION( "Two clusters across top-level snarl" ) {
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            vector<pos_t> pos_ts;
            pos_ts.emplace_back(1, false, 0);
            pos_ts.emplace_back(4, true, 0);

            for (pos_t pos : pos_ts){

                auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                seeds.push_back({ pos, 0, chain_info});
            }
            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 5); 

            REQUIRE( clusters.size() == 2);
        }
    }
    TEST_CASE( "Chain connected to node in top-level snarl","[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GGGGGGGGGGGG");//12 Gs
        Node* n6 = graph.create_node("T");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n1, n3, false, true);
        Edge* e4 = graph.create_edge(n2, n4);
        Edge* e5 = graph.create_edge(n2, n5);
        Edge* e6 = graph.create_edge(n3, n5);
        Edge* e7 = graph.create_edge(n4, n5);
        Edge* e8 = graph.create_edge(n5, n6);
        Edge* e9 = graph.create_edge(n5, n6, false, true);

        
        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);

        SnarlDistanceIndexClusterer clusterer(dist_index, &graph);

        SECTION( "One cluster across top-level snarl" ) {
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            vector<pos_t> pos_ts;
            pos_ts.emplace_back(2, false, 0);
            pos_ts.emplace_back(6, true, 0);

            for (pos_t pos : pos_ts){

                auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                seeds.push_back({ pos, 0, chain_info});
            }
            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 20); 

            REQUIRE( clusters.size() == 1);
        }
        SECTION( "Two clusters across top-level snarl" ) {
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            vector<pos_t> pos_ts;
            pos_ts.emplace_back(2, false, 0);
            pos_ts.emplace_back(6, true, 0);

            for (pos_t pos : pos_ts){

                auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                seeds.push_back({ pos, 0, chain_info});
            }
            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 5); 

            REQUIRE( clusters.size() == 2);
        }
    }
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
 
        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);

        SnarlDistanceIndexClusterer clusterer(dist_index, &graph);


        SECTION( "Two clusters in a chain and loop of snarl boundary" ) {
            vector<pos_t> pos_ts;
            pos_ts.emplace_back(2, false, 0);
            pos_ts.emplace_back(3, false, 0);
            pos_ts.emplace_back(5, false, 0);
            pos_ts.emplace_back(16, false, 0);
            //New cluster
            pos_ts.emplace_back(5, false, 10);
            pos_ts.emplace_back(6, false, 0);
            pos_ts.emplace_back(8, false, 0);

            for (bool use_minimizers : {true, false}) {
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (pos_t pos : pos_ts){

                    if (use_minimizers) {
                        auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                        seeds.push_back({ pos, 0,chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 3); 

                REQUIRE( clusters.size() == 2);

                vector<hash_set<size_t>> cluster_sets;
                for (auto& c : clusters) {
                    hash_set<size_t> h;
                    for (size_t s : c.seeds) {
                        h.insert(s);
                    }
                    cluster_sets.push_back(h);
                }
                REQUIRE (( (cluster_sets[0].count(0) == 1 &&
                           cluster_sets[0].count(1) == 1 &&
                           cluster_sets[0].count(2) == 1 &&
                           cluster_sets[0].count(3) == 1 &&
                           cluster_sets[1].count(4) == 1 &&
                           cluster_sets[1].count(5) == 1 &&
                           cluster_sets[1].count(6) == 1)  ||

                         ( cluster_sets[1].count(0) == 1 &&
                           cluster_sets[1].count(1) == 1 &&
                           cluster_sets[1].count(2) == 1 &&
                           cluster_sets[1].count(3) == 1 &&
                           cluster_sets[0].count(4) == 1 &&
                           cluster_sets[0].count(5) == 1 &&
                           cluster_sets[0].count(6) == 1  )));
            }
        }
        SECTION( "Four clusters" ) {
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            vector<pos_t> pos_ts;
            pos_ts.emplace_back(3, false, 0);
            pos_ts.emplace_back(5, false, 0);
            pos_ts.emplace_back(16, false, 0);
            //New cluster
            pos_ts.emplace_back(5, false, 8);
            //new_cluster
            pos_ts.emplace_back(6, false, 0);
            pos_ts.emplace_back(8, false, 0);
            //New_cluster
            pos_ts.emplace_back(13, false, 1);
            pos_ts.emplace_back(14, false, 0);
            pos_ts.emplace_back(15, false, 0);

            for (pos_t pos : pos_ts){
                seeds.push_back({ pos, 0});
            }
            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 3); 

            REQUIRE( clusters.size() == 4);

            vector<vector<SnarlDistanceIndexClusterer::Seed>> all_seeds;

            all_seeds.push_back(seeds);

            vector<vector<SnarlDistanceIndexClusterer::Cluster>> paired_clusters = clusterer.cluster_seeds(all_seeds, 3, 3); 

            REQUIRE( paired_clusters.size() == 1);
            REQUIRE( paired_clusters[0].size() == 4);
            REQUIRE( paired_clusters[0][0].fragment != paired_clusters[0][1].fragment);
            REQUIRE( paired_clusters[0][0].fragment != paired_clusters[0][2].fragment);
            REQUIRE( paired_clusters[0][0].fragment != paired_clusters[0][3].fragment);
            REQUIRE( paired_clusters[0][1].fragment != paired_clusters[0][2].fragment);
            REQUIRE( paired_clusters[0][1].fragment != paired_clusters[0][3].fragment);
            REQUIRE( paired_clusters[0][2].fragment != paired_clusters[0][3].fragment);

            //New fragment clusters
        } SECTION ("Four fragment clusters") {
            vector<vector<SnarlDistanceIndexClusterer::Seed>> all_seeds;
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            vector<pos_t>pos_ts;
            pos_ts.emplace_back(3, false, 0);
            pos_ts.emplace_back(5, false, 0);
            pos_ts.emplace_back(16, false, 0);
            //New cluster
            pos_ts.emplace_back(6, false, 0);
            pos_ts.emplace_back(8, false, 0);
            for (pos_t pos : pos_ts){
                seeds.push_back({ pos, 0});
            }
            all_seeds.push_back(seeds);
            seeds.clear();
            pos_ts.clear();
            //New cluster
            pos_ts.emplace_back(5, false, 8);
            //New cluster
            pos_ts.emplace_back(13, false, 1);
            pos_ts.emplace_back(14, false, 0);
            pos_ts.emplace_back(15, false, 0);
            for (pos_t pos : pos_ts){
                seeds.push_back({ pos, 0});
            }
            all_seeds.push_back(seeds);

            vector<vector<SnarlDistanceIndexClusterer::Cluster>> paired_clusters = clusterer.cluster_seeds(all_seeds, 3, 3);

            REQUIRE( paired_clusters.size() == 2);
            REQUIRE( paired_clusters[0].size() == 2);
            REQUIRE( paired_clusters[1].size() == 2);
            REQUIRE( paired_clusters[0][0].fragment != paired_clusters[0][1].fragment);
            REQUIRE( paired_clusters[0][0].fragment != paired_clusters[1][0].fragment);
            REQUIRE( paired_clusters[0][0].fragment != paired_clusters[1][1].fragment);
            REQUIRE( paired_clusters[0][1].fragment != paired_clusters[1][0].fragment);
            REQUIRE( paired_clusters[0][1].fragment != paired_clusters[1][1].fragment);

            //New fragment clusters

            paired_clusters = clusterer.cluster_seeds(all_seeds, 3, 5);

            REQUIRE( paired_clusters.size() == 2);
            REQUIRE( paired_clusters[0].size() == 2);
            REQUIRE( paired_clusters[1].size() == 2);
            REQUIRE( paired_clusters[0][0].fragment == paired_clusters[0][1].fragment);
            REQUIRE( paired_clusters[0][0].fragment == paired_clusters[1][0].fragment);
            REQUIRE( paired_clusters[0][0].fragment != paired_clusters[1][1].fragment);
            REQUIRE( paired_clusters[0][1].fragment == paired_clusters[1][0].fragment);
            REQUIRE( paired_clusters[0][1].fragment != paired_clusters[1][1].fragment);
            REQUIRE( paired_clusters[1][0].fragment != paired_clusters[1][1].fragment);
        }
        SECTION( "Same node, same cluster" ) {
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            vector<pos_t> pos_ts;
            pos_ts.emplace_back(5, false, 0);
            pos_ts.emplace_back(5, false, 11);
            pos_ts.emplace_back(5, false, 5);

            for (pos_t pos : pos_ts){
                seeds.push_back({ pos, 0});
            }
            vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 7); 


            REQUIRE( clusters.size() == 1);
        }
    }//end test case
    TEST_CASE("Top level root", "[cluster]") {
        VG graph;

        Node* n1 = graph.create_node("GTGCACA");//8
        Node* n2 = graph.create_node("GTGCACA");
        Node* n3 = graph.create_node("GT");
        Node* n4 = graph.create_node("GATTCTTATAG");//11

        Edge* e1 = graph.create_edge(n1, n3);
        Edge* e2 = graph.create_edge(n1, n4);
        Edge* e3 = graph.create_edge(n3, n2);
        Edge* e4 = graph.create_edge(n3, n4, false, true);
        Edge* e5 = graph.create_edge(n2, n4);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);


        SnarlDistanceIndexClusterer clusterer(dist_index, &graph);


        SECTION("One cluster") {
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            vector<pos_t> pos_ts;
            pos_ts.emplace_back(2, false, 0);
            pos_ts.emplace_back(1, false, 7);
            pos_ts.emplace_back(1, false, 2);
            pos_ts.emplace_back(1, true, 5);
            pos_ts.emplace_back(3, false, 3);

            for (pos_t pos : pos_ts){
                seeds.push_back({ pos, 0});
            }
            vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 10); 


            REQUIRE( clusters.size() == 1);
        }

    }
    TEST_CASE("Top level unary snarl", "[cluster]") {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n2, n7);
        Edge* e4 = graph.create_edge(n3, n4);
        Edge* e5 = graph.create_edge(n3, n5);
        Edge* e6 = graph.create_edge(n4, n6);
        Edge* e7 = graph.create_edge(n5, n6);
        Edge* e8 = graph.create_edge(n6, n7);
        Edge* e9 = graph.create_edge(n1, n1, true, false);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);


        SnarlDistanceIndexClusterer clusterer(dist_index, &graph);



        // We end up with a big unary snarl of 7 rev -> 7 rev
        // Inside that we have a chain of two normal snarls 2 rev -> 3 fwd,      and 3 fwd -> 6 fwd
        // And inside 2 rev -> 3 fwd, we get 1 rev -> 1 rev as another unar     y snarl.

        // We name the snarls for the distance index by their start nodes.
        SECTION("Distances in root") {
            net_handle_t root = dist_index.get_root();
            net_handle_t chain = dist_index.get_parent(dist_index.get_node_net_handle(1));
            REQUIRE(dist_index.get_parent(chain) == root);
        }

        SECTION("Top level cluster") {
            vector<id_t> ids({1, 2, 7});
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                seeds.push_back({ pos, 0});
            }

            vector<SnarlDistanceIndexClusterer::Cluster> clusters= clusterer.cluster_seeds(seeds, 10); 


            REQUIRE( clusters.size() == 1);
        }
        SECTION("One cluster") {
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            vector<pos_t> pos_ts;
            pos_ts.emplace_back(1, false, 0);
            pos_ts.emplace_back(2, false, 0);
            pos_ts.emplace_back(7, false, 0);
            pos_ts.emplace_back(4, false, 0);

            for (pos_t pos : pos_ts){
                seeds.push_back({ pos, 0});
            }
            vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 10); 


            REQUIRE( clusters.size() == 1);
        }
        SECTION("One cluster") {
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            vector<pos_t> pos_ts;
            pos_ts.emplace_back(2, false, 0);
            pos_ts.emplace_back(4, false, 0);

            for (pos_t pos : pos_ts){
                seeds.push_back({ pos, 0});
            }
            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 10); 



            REQUIRE( clusters.size() == 1);
        }
        SECTION("Two clusters") {
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            vector<pos_t> pos_ts;
            pos_ts.emplace_back(2, false, 0);
            pos_ts.emplace_back(4, false, 1);
            pos_ts.emplace_back(6, false, 0);

            for (pos_t pos : pos_ts){
                seeds.push_back({ pos, 0});
            }
            vector<SnarlDistanceIndexClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 5); 


            REQUIRE( clusters.size() == 2);
        }
        SECTION("No clusters") {
            vector<SnarlDistanceIndexClusterer::Seed> seeds;

            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 5); 


            REQUIRE( clusters.size() == 0);
        }
    }
    TEST_CASE( "Long chain",
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
        Node* n9 = graph.create_node("TTA");
        Node* n10 = graph.create_node("G");
        Node* n11 = graph.create_node("CTGA");
        Node* n12 = graph.create_node("G");
        Node* n13 = graph.create_node("CTGA");
        Node* n14 = graph.create_node("CTGA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n6);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n4);
        Edge* e5 = graph.create_edge(n3, n5);
        Edge* e6 = graph.create_edge(n4, n5);
        Edge* e7 = graph.create_edge(n5, n6);
        Edge* e8 = graph.create_edge(n6, n7);
        Edge* e9 = graph.create_edge(n6, n8);
        Edge* e10 = graph.create_edge(n7, n8);
        Edge* e11 = graph.create_edge(n8, n9);
        Edge* e12 = graph.create_edge(n8, n12);
        Edge* e13 = graph.create_edge(n9, n10);
        Edge* e14 = graph.create_edge(n9, n11);
        Edge* e15 = graph.create_edge(n10, n11);
        Edge* e16 = graph.create_edge(n11, n12);
        Edge* e17 = graph.create_edge(n12, n13);
        Edge* e18 = graph.create_edge(n12, n14);
        Edge* e19 = graph.create_edge(n13, n14);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);


        SnarlDistanceIndexClusterer clusterer(dist_index, &graph);

        SECTION("Snarl then seed") {

            vector<id_t> ids({3, 5, 6, 11});
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                seeds.push_back({ pos, 0});
            }

            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 5); 


            REQUIRE( clusters.size() == 2);

        }
        SECTION("Seed then snarl") {

            vector<id_t> ids({1, 2, 3, 5, 6, 11, 10});
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                seeds.push_back({ pos, 0});
            }

            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 5); 


            REQUIRE( clusters.size() == 2);

        }
        SECTION("Only seeds") {

            vector<id_t> ids({1, 6, 14});
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                seeds.push_back({ pos, 0});
            }

            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 4); 


            REQUIRE( clusters.size() == 2);

        }
        SECTION("Only seeds two reads") {

            vector<id_t> ids({1, 6, 14});
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                seeds.push_back({ pos, 0});
            }
            vector<id_t> ids1({8, 12});
            vector<SnarlDistanceIndexClusterer::Seed> seeds1;
            for (id_t n : ids1) {
                pos_t pos = make_pos_t(n, false, 0);
                seeds1.push_back({ pos, 0});
            }

            vector<vector<SnarlDistanceIndexClusterer::Seed>> all_seeds;
            all_seeds.emplace_back(seeds);
            all_seeds.emplace_back(seeds1);

            vector<vector<SnarlDistanceIndexClusterer::Cluster>> clusters =  clusterer.cluster_seeds(all_seeds, 4, 5); 


            REQUIRE( clusters.size() == 2);
            REQUIRE( clusters[0].size() == 2);
            REQUIRE( clusters[1].size() == 1);
            REQUIRE( clusters[0][0].fragment == clusters[0][1].fragment);
            REQUIRE( clusters[0][0].fragment == clusters[1][0].fragment);

        }
        SECTION("Only snarls") {

            vector<id_t> ids({4, 5, 9});
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                seeds.push_back({ pos, 0});
            }

            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 9); 


            REQUIRE( clusters.size() == 1);

        }
        SECTION("Skip snarl") {

            vector<id_t> ids({7, 10, 13});
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                seeds.push_back({ pos, 0});
            }

            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 6); 

            REQUIRE( clusters.size() == 1);
        }
    }

    TEST_CASE( "Disconnected graph",
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
//Disconnected
        Node* n9 = graph.create_node("T");
        Node* n10 = graph.create_node("G");
        Node* n11 = graph.create_node("CTGA");
        Node* n12 = graph.create_node("G");
        Node* n13 = graph.create_node("CTGA");

        Node* n14 = graph.create_node("AGCCGTGTGC");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n3, n4);
        Edge* e5 = graph.create_edge(n3, n5);
        Edge* e6 = graph.create_edge(n4, n5);
        Edge* e7 = graph.create_edge(n5, n6);
        Edge* e8 = graph.create_edge(n5, n7);
        Edge* e9 = graph.create_edge(n6, n8);
        Edge* e10 = graph.create_edge(n7, n8);

        Edge* e11 = graph.create_edge(n9, n10);
        Edge* e12 = graph.create_edge(n9, n11);
        Edge* e13 = graph.create_edge(n10, n11);
        Edge* e14 = graph.create_edge(n11, n12);
        Edge* e15 = graph.create_edge(n11, n13);
        Edge* e16 = graph.create_edge(n12, n13);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);


        SnarlDistanceIndexClusterer clusterer(dist_index, &graph);

        SECTION("Two clusters") {
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            vector<pos_t> pos_ts;
            pos_ts.emplace_back(2, false, 0);
            pos_ts.emplace_back(3, false, 0);
            pos_ts.emplace_back(9, false, 0);

            for (pos_t pos : pos_ts){
                seeds.push_back({ pos, 0});
            }
            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 5); 


            REQUIRE( clusters.size() == 2);

        }
        SECTION("Two clusters with seed structs") {

            vector<id_t> ids({2, 3, 9});
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                seeds.push_back({ pos, 0});
            }

            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 5); 


            REQUIRE( clusters.size() == 2);

        }
        SECTION("Two clusters with seed structs") {

            vector<id_t> ids({2, 3, 5, 9, 10});
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                seeds.push_back({ pos, 0});
            }

            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 5); 


            REQUIRE( clusters.size() == 2);

        }
        SECTION("Two top level clusters") {

            vector<id_t> ids({1, 3, 11});
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                seeds.push_back({ pos, 0});
            }
            vector<id_t> ids1({5, 13});
            vector<SnarlDistanceIndexClusterer::Seed> seeds1;
            for (id_t n : ids1) {
                pos_t pos = make_pos_t(n, false, 0);
                seeds1.push_back({ pos, 0});
            }
            //Clusters are 
            //Read 1: {1, 3} in a fragment cluster with Read 2: {5}
            //Read 1: {11} in a fragment cluster with Read 2: {13}
            vector<vector<SnarlDistanceIndexClusterer::Seed>> all_seeds;
            all_seeds.emplace_back(seeds);
            all_seeds.emplace_back(seeds1);


            vector<vector<SnarlDistanceIndexClusterer::Cluster>> clusters =  clusterer.cluster_seeds(all_seeds, 5, 10); 


            REQUIRE( clusters.size() == 2);
            REQUIRE( clusters[0].size() == 2);
            REQUIRE( clusters[1].size() == 2);
            REQUIRE( clusters[0][0].fragment != clusters[0][1].fragment);
            REQUIRE( clusters[1][0].fragment != clusters[1][1].fragment);

            REQUIRE(( clusters[0][0].fragment == clusters[1][0].fragment || clusters[0][0].fragment == clusters[1][1].fragment));
            REQUIRE(( clusters[0][1].fragment == clusters[1][0].fragment || clusters[0][1].fragment == clusters[1][1].fragment));


        }
        SECTION("Disconnected node") {

            vector<id_t> ids({1, 3, 11, 14, 14});
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                seeds.push_back({ pos, 0});
            }
            vector<id_t> ids1({5, 13});
            vector<SnarlDistanceIndexClusterer::Seed> seeds1;
            for (id_t n : ids1) {
                pos_t pos = make_pos_t(n, false, 0);
                seeds1.push_back({ pos, 0});
            }
            //Clusters are 
            //Read 1: {1, 3} in a fragment cluster with Read 2: {5}
            //Read 1: {11} in a fragment cluster with Read 2: {13}
            //Read 1 : {14, 14}
            vector<vector<SnarlDistanceIndexClusterer::Seed>> all_seeds;
            all_seeds.emplace_back(seeds);
            all_seeds.emplace_back(seeds1);


            vector<vector<SnarlDistanceIndexClusterer::Cluster>> clusters =  clusterer.cluster_seeds(all_seeds, 5, 10); 


            REQUIRE( clusters.size() == 2);
            REQUIRE( clusters[0].size() == 3);
            REQUIRE( clusters[1].size() == 2);
            REQUIRE( clusters[0][0].fragment != clusters[0][1].fragment);
            REQUIRE( clusters[1][0].fragment != clusters[1][1].fragment);

            REQUIRE(( clusters[0][0].fragment == clusters[1][0].fragment || clusters[0][0].fragment == clusters[1][1].fragment));
            REQUIRE(( clusters[0][1].fragment == clusters[1][0].fragment || clusters[0][1].fragment == clusters[1][1].fragment));


        }
    }
    TEST_CASE("Simple nested chain", "[cluster]") {
        VG graph;
 
        Node* n1 = graph.create_node("GAC");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("GCAA");
        Node* n9 = graph.create_node("GTGACTAAGA");//10
        Node* n10 = graph.create_node("GTGACTAAGA");//10
 
        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n10);
        Edge* e3 = graph.create_edge(n10, n8);
        Edge* e4 = graph.create_edge(n2, n3);
        Edge* e5 = graph.create_edge(n2, n4);
        Edge* e6 = graph.create_edge(n3, n4);
        Edge* e7 = graph.create_edge(n4, n5);
        Edge* e8 = graph.create_edge(n4, n6);
        Edge* e9 = graph.create_edge(n5, n6);
        Edge* e10 = graph.create_edge(n6, n7);
        Edge* e11 = graph.create_edge(n7, n8);
        Edge* e12 = graph.create_edge(n8, n9);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);


        SnarlDistanceIndexClusterer clusterer(dist_index, &graph);

        SECTION("Only seeds on nodes in inner chain one cluster") {
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            vector<pos_t> pos_ts;
            pos_ts.emplace_back(2, false, 0);
            pos_ts.emplace_back(4, false, 0);
            pos_ts.emplace_back(7, false, 0);

            for (pos_t pos : pos_ts){
                seeds.push_back({ pos, 0});
            }
            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 7); 


            REQUIRE( clusters.size() == 1);

        }
        SECTION("Only seeds on nodes in inner chain two clusters") {
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            vector<pos_t> pos_ts;
            pos_ts.emplace_back(2, false, 0);
            pos_ts.emplace_back(4, false, 0);
            pos_ts.emplace_back(7, false, 0);

            for (pos_t pos : pos_ts){
                seeds.push_back({ pos, 0});
            }
            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 4); 


            REQUIRE( clusters.size() == 2);

        }
        SECTION("Only seeds on nodes in inner chain two clusters with outer nodes") {
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            vector<pos_t> pos_ts;
            pos_ts.emplace_back(1, false, 0);
            pos_ts.emplace_back(2, false, 0);
            pos_ts.emplace_back(4, false, 0);
            pos_ts.emplace_back(7, false, 0);
            pos_ts.emplace_back(8, true, 0);

            for (pos_t pos : pos_ts){
                seeds.push_back({ pos, 0});
            }
            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 4); 


            REQUIRE( clusters.size() == 2);
            REQUIRE((clusters[0].seeds.size() == 3 || clusters[0].seeds.size() == 2));
            REQUIRE((clusters[1].seeds.size() == 3 || clusters[1].seeds.size() == 2));

        }
        SECTION("One fragment cluster") {
            vector<vector<pos_t>> pos_ts;
            pos_ts.emplace_back();
            pos_ts.emplace_back();
            pos_ts[0].emplace_back(1, false, 0);
            pos_ts[0].emplace_back(2, false, 0);
            pos_ts[0].emplace_back(4, false, 0);
            pos_ts[1].emplace_back(7, false, 0);
            pos_ts[1].emplace_back(8, true, 0);

            for (bool use_minimizers : {true, false}) {
                vector<vector<SnarlDistanceIndexClusterer::Seed>> seeds(2);
                for (size_t read_num = 0 ; read_num < pos_ts.size() ; read_num ++) {
                    for (pos_t pos : pos_ts[read_num]){
                        if (use_minimizers) {
                            auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                            seeds[read_num].push_back({ pos, 0, chain_info});
                        } else {
                            seeds[read_num].push_back({ pos, 0});
                        }
                    }
                }
                
                vector<vector<SnarlDistanceIndexClusterer::Cluster>> clusters =  clusterer.cluster_seeds(seeds, 4, 10); 

                REQUIRE( clusters.size() == 2);
                REQUIRE(clusters[0].size() == 1);
                REQUIRE(clusters[1].size() == 1);
                REQUIRE(clusters[0][0].fragment == clusters[1][0].fragment);
            }


        }
        SECTION("One fragment cluster") {
            vector<pos_t> pos_ts;
            pos_ts.emplace_back(1, false, 0);
            pos_ts.emplace_back(2, false, 0);
            pos_ts.emplace_back(3, false, 0);
            pos_ts.emplace_back(4, false, 0);
            pos_ts.emplace_back(5, false, 0);
            pos_ts.emplace_back(7, false, 0);
            pos_ts.emplace_back(8, true, 0);

            for (bool use_minimizers : {true, false}) {
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (pos_t pos : pos_ts){
                    if (use_minimizers) {
                        auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                        seeds.push_back({ pos, 0, chain_info});
                    } else {
                        seeds.push_back({ pos, 0});
                    }
                }
                REQUIRE(clusterer.distance_between_seeds(seeds[0], seeds[1],false) == 3);
                REQUIRE(clusterer.distance_between_seeds(seeds[1], seeds[0],false) == 3);
                REQUIRE(clusterer.distance_between_seeds(seeds[0], seeds[2],false) == 4);
                REQUIRE(clusterer.distance_between_seeds(seeds[0], seeds[3],false) == 4);
                REQUIRE(clusterer.distance_between_seeds(seeds[2], seeds[4],false) == 5);
            }


        }
    }//End test case

    TEST_CASE("Top level loop creates looping chain", "[cluster]") {
        VG graph;
 
        Node* n1 = graph.create_node("G");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGAAAAAAAAAAAA"); //15
        Node* n5 = graph.create_node("GCAA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("A");
        Node* n9 = graph.create_node("T");
        Node* n10 = graph.create_node("G");
        Node* n11 = graph.create_node("GGGGG");
 
        Edge* e1 = graph.create_edge(n9, n1);
        Edge* e2 = graph.create_edge(n9, n11);
        Edge* e3 = graph.create_edge(n1, n2);
        Edge* e4 = graph.create_edge(n1, n8);
        Edge* e5 = graph.create_edge(n2, n3);
        Edge* e6 = graph.create_edge(n2, n4);
        Edge* e7 = graph.create_edge(n3, n5);
        Edge* e8 = graph.create_edge(n4, n5);
        Edge* e9 = graph.create_edge(n5, n6);
        Edge* e10 = graph.create_edge(n5, n7);
        Edge* e11 = graph.create_edge(n6, n7);
        Edge* e12 = graph.create_edge(n7, n8);
        Edge* e13 = graph.create_edge(n8, n10);
        Edge* e16 = graph.create_edge(n10, n9);
        Edge* e17 = graph.create_edge(n2, n2, true, false);
        Edge* e18 = graph.create_edge(n11, n10);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);


        SnarlDistanceIndexClusterer clusterer(dist_index, &graph);

        SECTION("Two clusters") {
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            vector<pos_t> pos_ts;
            pos_ts.emplace_back(2, false, 0);
            pos_ts.emplace_back(3, false, 0);
            pos_ts.emplace_back(8, false, 0);

            for (pos_t pos : pos_ts){
                seeds.push_back({ pos, 0});
            }
            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 3); 


            REQUIRE( clusters.size() == 2);

        }
        SECTION("One cluster") {
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            vector<pos_t> pos_ts;
            pos_ts.emplace_back(1, false, 0);
            pos_ts.emplace_back(2, false, 0);
            pos_ts.emplace_back(7, false, 0);

            for (pos_t pos : pos_ts){
                seeds.push_back({ pos, 0});
            }
            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 6); 


            REQUIRE( clusters.size() == 1);

        }
        SECTION("One cluster taking chain loop") {
            vector<id_t> ids({8, 9, 10});
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                seeds.push_back({ pos, 0});
            }

            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 3); 


            REQUIRE( clusters.size() == 1);

        }
    }//End test case


    TEST_CASE( "Nested unary snarls","[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("G");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("G");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n2, n4);
        Edge* e4 = graph.create_edge(n3, n4);
        Edge* e5 = graph.create_edge(n4, n5);
        Edge* e6 = graph.create_edge(n4, n6);
        Edge* e7 = graph.create_edge(n5, n6);
        Edge* e8 = graph.create_edge(n6, n7);
        Edge* e9 = graph.create_edge(n6, n8);
        Edge* e10 = graph.create_edge(n7, n8);
        Edge* e11 = graph.create_edge(n8, n8, false, true);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);



        SnarlDistanceIndexClusterer clusterer(dist_index, &graph);
        //Unary snarl at 8 nested in unary snarl at 6 nested in 
        //unary snarl at  4 nested in regular snarl at 2 (ending at 3)
        //nested in unary snarl at 1

        SECTION( "One cluster" ) {
            vector<id_t> ids({4, 3});
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                seeds.push_back({ pos, 0});
            }

            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 10); 


            REQUIRE( clusters.size() == 1);
        }
        SECTION( "One cluster nested" ) {
            vector<id_t> ids({5, 3});
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                seeds.push_back({ pos, 0});
            }

            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 10); 


            REQUIRE( clusters.size() == 1);
        }
        SECTION( "Three clusters" ) {
            vector<id_t> ids({2, 3, 8});
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                seeds.push_back({ pos, 0});
            }

            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 3); 



            REQUIRE( clusters.size() == 3);
        }
        SECTION( "One cluster taking loop" ) {
            vector<id_t> ids({2, 3});
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                seeds.push_back({ pos, 0});
            }

            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 15); 


            REQUIRE( clusters.size() == 1);
        }
    }//end test case
    TEST_CASE( "Top level snarl","[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n2, n5);
        Edge* e4 = graph.create_edge(n3, n4);
        Edge* e5 = graph.create_edge(n3, n5);
        Edge* e6 = graph.create_edge(n4, n5);


        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);

        SnarlDistanceIndexClusterer clusterer(dist_index, &graph);

        SECTION( "Top level seeds" ) {
            vector<id_t> ids({1, 2, 4});
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                seeds.push_back({ pos, 0});
            }

            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 3); 


            REQUIRE( clusters.size() == 2);
        }
    }
    TEST_CASE( "Two tip right","[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("GACCT");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("G");
        Node* n6 = graph.create_node("CTGA");
        Node* n7 = graph.create_node("G");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n3, n4);
        Edge* e5 = graph.create_edge(n3, n5);
        Edge* e6 = graph.create_edge(n6, n1);
        Edge* e7 = graph.create_edge(n6, n7);
        Edge* e8 = graph.create_edge(n7, n1);
        Edge* e9 = graph.create_edge(n1, n1, true, false);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);



        SnarlDistanceIndexClusterer clusterer(dist_index, &graph);

        SECTION( "Two cluster" ) {
            vector<id_t> ids({4, 5});
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                seeds.push_back({ pos, 0});
            }

            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 5); 


            REQUIRE( clusters.size() == 2);
        }

        SECTION( "One clusters" ) {
            vector<id_t> ids({4, 5, 3});
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                seeds.push_back({ pos, 0});
            }

            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 10); 


            REQUIRE( clusters.size() == 1);
        }

        SECTION( "One cluster loop" ) {
            vector<id_t> ids({4, 5});
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                seeds.push_back({ pos, 0});
            }

            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 18); 


            REQUIRE( clusters.size() == 1);
        }
    }
    TEST_CASE( "Two tips","[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("CATCCTCCTCGATT");//14
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("GA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n1, n1);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);

        SnarlDistanceIndexClusterer clusterer(dist_index, &graph);
        net_handle_t snarl = dist_index.get_parent(dist_index.get_parent(dist_index.get_node_net_handle(n1->id())));
        REQUIRE(!dist_index.is_simple_snarl(snarl));

        SECTION( "One cluster" ) {
            vector<pos_t> positions;
            positions.emplace_back(make_pos_t(1, true, 8));
            positions.emplace_back(make_pos_t(3, false, 1));
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (auto pos : positions) {
                seeds.push_back({pos, 0});
            }

            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 10); 


            REQUIRE( clusters.size() == 1);
        }

    }
    TEST_CASE( "Two tip left","[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("G");
        Node* n6 = graph.create_node("G");
        Node* n7 = graph.create_node("CTGA");

        Edge* e1 = graph.create_edge(n1, n3);
        Edge* e2 = graph.create_edge(n2, n3);
        Edge* e3 = graph.create_edge(n3, n4);
        Edge* e4 = graph.create_edge(n3, n5);
        Edge* e5 = graph.create_edge(n4, n5);
        Edge* e6 = graph.create_edge(n5, n6);
        Edge* e7 = graph.create_edge(n5, n7);
        Edge* e8 = graph.create_edge(n6, n7);
        Edge* e9 = graph.create_edge(n5, n5, false, true);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);

        SnarlDistanceIndexClusterer clusterer(dist_index, &graph);

        SECTION( "One cluster" ) {
            vector<id_t> ids({1, 2});
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                seeds.push_back({  pos, 0});
            }

            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 5); 


            REQUIRE( clusters.size() == 2);
        }

        SECTION( "Two clusters" ) {
            vector<id_t> ids({1, 2, 3});
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                seeds.push_back({ pos, 0});
            }

            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 5); 


            REQUIRE( clusters.size() == 1);
        }
        SECTION( "Two clusters with snarl" ) {
            vector<id_t> ids({1, 2, 4});
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                seeds.push_back({ pos, 0});
            }

            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 5); 


            REQUIRE( clusters.size() == 1);
        }
        SECTION( "One cluster with loop" ) {
            vector<id_t> ids({1, 2});
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                seeds.push_back({ pos, 0});
            }

            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 7); 


            REQUIRE( clusters.size() == 1);
        }
    }
    TEST_CASE( "trivial snarls on the ends of a chain","[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("G");
        Node* n7 = graph.create_node("C");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n2, n3);
        Edge* e3 = graph.create_edge(n3, n4);
        Edge* e4 = graph.create_edge(n3, n5);
        Edge* e5 = graph.create_edge(n4, n5);
        Edge* e6 = graph.create_edge(n5, n6);
        Edge* e7 = graph.create_edge(n6, n7);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);

        SnarlDistanceIndexClusterer clusterer(dist_index, &graph);

        SECTION( "One cluster" ) {
            vector<id_t> ids({1, 2});
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                seeds.push_back({  pos, 0});
            }

            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 5); 


            REQUIRE( clusters.size() == 1);
        }

        SECTION( "One cluster across snarl" ) {
            vector<id_t> ids({2, 6});
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                seeds.push_back({ pos, 0});
            }

            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 7); 


            REQUIRE( clusters.size() == 1);
        }
        SECTION( "Two clusters " ) {
            vector<id_t> ids({1, 6});
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                seeds.push_back({ pos, 0});
            }

            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 5); 


            REQUIRE( clusters.size() == 2);
        }
        SECTION( "One cluster with snarl" ) {
            vector<id_t> ids({1, 2, 4, 6});
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                seeds.push_back({ pos, 0});
            }

            vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 7); 


            REQUIRE( clusters.size() == 1);
        }
    }



    //TEST_CASE("Load graph", "[cluster][load_cluster]"){

    //    ifstream vg_stream("testGraph.hg");
    //    HashGraph graph(vg_stream);
    //    vg_stream.close();
    //    IntegratedSnarlFinder snarl_finder(graph);
    //    SnarlDistanceIndex dist_index;
    //    fill_in_distance_index(&dist_index, &graph, &snarl_finder);


    //    SnarlDistanceIndexClusterer clusterer(dist_index, &graph);

    //    size_t read_lim = 10;// Distance between read clusters
    //    size_t fragment_lim = 15;// Distance between fragment clusters



    //    vector<SnarlDistanceIndexClusterer::Seed> seeds;
    //    vector<pos_t> pos_ts;
    //    pos_ts.emplace_back(6, false, 4);
    //    pos_ts.emplace_back(8, false, 0);
    //    pos_ts.emplace_back(9, false, 0);

    //    for (pos_t pos : pos_ts) {
    //        auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
    //        seeds.push_back({ pos, 0, chain_info});
    //    }
    //    vector<SnarlDistanceIndexClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, read_lim); 
    //    REQUIRE(clusters.size() == 1);
    //}//end test case

    /*
    TEST_CASE("Failed graph", "[failed_cluster]"){

        HashGraph graph;
        graph.deserialize("testGraph.hg");
        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);


        dist_index.print_self();

        SnarlDistanceIndexClusterer clusterer(dist_index, &graph);



        vector<vector<pos_t>> pos_ts(2);
        pos_ts[0].emplace_back(30, false, 0);
        pos_ts[0].emplace_back(22, false, 0);
        pos_t pos1 = pos_ts[0][0];
        pos_t pos2 = pos_ts[0][1];
        net_handle_t node31 = dist_index.get_node_net_handle(30);
        
       size_t dist = dist_index.minimum_distance(get_id(pos1), get_is_rev(pos1), get_offset(pos1), get_id(pos2), get_is_rev(pos2), get_offset(pos2), true, &graph);
       cerr << "DISTANCE BETWEEN " << pos1 << " and " << pos2 << " = " << dist << endl;

        //for (bool use_minimizers : {true, false}) {

        //    vector<vector<SnarlDistanceIndexClusterer::Seed>> seeds(2);
        //    for (size_t read_num = 0 ; read_num < pos_ts.size() ; read_num++) {
        //            for (pos_t pos : pos_ts[read_num]) {

        //                if (use_minimizers) {
        //                    auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
        //                    seeds[read_num].push_back({ pos, 0, chain_info});
        //                } else {
        //                    seeds[read_num].push_back({ pos, 0});
        //                }
        //            }
        //    }

        //    vector<vector<SnarlDistanceIndexClusterer::Cluster>> clusters =  clusterer.cluster_seeds(seeds, 15, 35); 

        //    REQUIRE(clusters.size() == 1);
        //}
        REQUIRE(false);
    }
    */
    TEST_CASE("Random graphs", "[cluster_random]"){


        for (int i = 0; i < 0; i++) {
            // For each random graph
            
            default_random_engine generator(time(NULL));
            uniform_int_distribution<int> variant_count(1, 70);
            uniform_int_distribution<int> chrom_len(10, 200);

            //Make a random graph with three chromosomes of random lengths
            HashGraph graph;
            random_graph({chrom_len(generator),chrom_len(generator),chrom_len(generator)}, 30, variant_count(generator), &graph);
            graph.serialize("testGraph.hg");


            IntegratedSnarlFinder snarl_finder(graph);
            SnarlDistanceIndex dist_index;
            fill_in_distance_index(&dist_index, &graph, &snarl_finder, 5);
    


            SnarlDistanceIndexClusterer clusterer(dist_index, &graph);


            vector<id_t> all_nodes;
            graph.for_each_handle([&](const handle_t& h)->bool{
                id_t id = graph.get_id(h);
                all_nodes.push_back(id);
                return true;
            });


            uniform_int_distribution<int> randPosIndex(0, all_nodes.size()-1);
            for (bool use_minimizers : {true, false}) {

                for (size_t k = 0; k < 10 ; k++) {

                    vector<vector<SnarlDistanceIndexClusterer::Seed>> all_seeds(2);
                    size_t read_lim = 15;// Distance between read clusters
                    size_t fragment_lim = 35;// Distance between fragment clusters
                    for (size_t read = 0 ; read < 2 ; read ++) {
                        uniform_int_distribution<int> randPosCount(3, 70);
                        for (int j = 0; j < randPosCount(generator); j++) {
                            //Check clusters of j random positions 
 
                            id_t nodeID1 = all_nodes[randPosIndex(generator)];
                            handle_t node1 = graph.get_handle(nodeID1);
 
                            offset_t offset1 = uniform_int_distribution<int>(0,graph.get_length(node1) - 1)(generator);

                            pos_t pos = make_pos_t(nodeID1,
                                uniform_int_distribution<int>(0,1)(generator) == 0,offset1 );


                        
                        if (use_minimizers) {
                            auto chain_info = MIPayload::encode(get_minimizer_distances(dist_index, pos));
                            all_seeds[read].push_back({ pos, 0, chain_info});
                        } else {
                            all_seeds[read].push_back({ pos, 0});
                        }

                        }
                    }
                    vector<vector<SnarlDistanceIndexClusterer::Cluster>> paired_clusters = clusterer.cluster_seeds(all_seeds, read_lim, fragment_lim); 
                   
                    vector<vector<pos_t>> fragment_clusters;

                    for (size_t read_num = 0 ; read_num < 2 ; read_num ++) {
                        auto& one_read_clusters = paired_clusters[read_num];
                        if (one_read_clusters.size() > 0) {
                            for (size_t a = 0; a < one_read_clusters.size(); a++) {
                                // For each cluster -cluster this cluster to ensure that 
                                // there is only one
                                vector<size_t> clust = one_read_clusters[a].seeds;
                                size_t fragment_cluster = one_read_clusters[a].fragment;
                                if (fragment_cluster >= fragment_clusters.size()) {
                                    fragment_clusters.resize(fragment_cluster+1);
                                }
                                
                                structures::UnionFind new_clusters (clust.size(), false);

                                for (size_t i1 = 0 ; i1 < clust.size() ; i1++) {
                                    pos_t pos1 = all_seeds[read_num][clust[i1]].pos;
                                    fragment_clusters[fragment_cluster].emplace_back(pos1);
                                    size_t len1 = dist_index.minimum_length(dist_index.get_node_net_handle(get_id(pos1)));;
                                    pos_t rev1 = make_pos_t(get_id(pos1), !is_rev(pos1),len1 - get_offset(pos1)-1); 

                                    for (size_t b = 0 ; b < one_read_clusters.size() ; b++) {
                                        if (b != a) {
                                            //For each other cluster
                                            vector<size_t> clust2 = one_read_clusters[b].seeds;
                                            for (size_t i2 = 0 ; i2 < clust2.size() ; i2++) {
                                                //And each position in each other cluster,
                                                //make sure that this position is far away from i1
                                                pos_t pos2 = all_seeds[read_num][clust2[i2]].pos;
                                                size_t len2 = dist_index.minimum_length(dist_index.get_node_net_handle(get_id(pos2)));
                                                pos_t rev2 = make_pos_t(get_id(pos2), 
                                                                 !is_rev(pos2),
                                                                 len2 - get_offset(pos2)-1); 

                                                size_t dist1 = dist_index.minimum_distance(get_id(pos1), get_is_rev(pos1), get_offset(pos1), get_id(pos2), get_is_rev(pos2), get_offset(pos2), false, &graph);
                                                size_t dist2 = dist_index.minimum_distance(get_id(pos1), get_is_rev(pos1), get_offset(pos1), get_id(rev2), get_is_rev(rev2), get_offset(rev2), false, &graph);
                                                size_t dist3 = dist_index.minimum_distance(get_id(rev1), get_is_rev(rev1), get_offset(rev1), get_id(pos2), get_is_rev(pos2), get_offset(pos2), false, &graph);
                                                size_t dist4 = dist_index.minimum_distance(get_id(rev1), get_is_rev(rev1), get_offset(rev1), get_id(rev2), get_is_rev(rev2), get_offset(rev2), false, &graph);
                                                size_t dist = std::min(std::min(dist1, 
                                                                   dist2), std::min( dist3, dist4));
                                                if (dist != clusterer.distance_between_seeds(all_seeds[read_num][clust[i1]],
                                                                                                 all_seeds[read_num][clust2[i2]],

                                                                                                 false)) {
                                                    graph.serialize("testGraph.hg");
                                                    cerr << "Distance between " << pos1 << " and " << pos2 << " should be " << dist << endl;

                                                }
                                                REQUIRE(dist == clusterer.distance_between_seeds(all_seeds[read_num][clust[i1]],
                                                                                                 all_seeds[read_num][clust2[i2]],

                                                                                                 false));
                                                if ( dist != -1 && dist <= read_lim) {
                                                    dist_index.print_self();
                                                    graph.serialize("testGraph.hg");
                                                    cerr << "These should have been in the same read cluster: " ;
                                                    cerr << pos1 << " and " << pos2 << endl;
                                                    cerr << dist1 << " " << dist2 << " " << dist3 << " " << dist4 << endl;
                                                    REQUIRE(false);
                                                }
                                                
                                            }
                                        }
                                    }
                                    for (size_t i2 = 0 ; i2 < clust.size() ; i2++) {
                                        //For each position in the same cluster
                                        pos_t pos2 = all_seeds[read_num][clust[i2]].pos;
                                        size_t len2 = dist_index.minimum_length(dist_index.get_node_net_handle(get_id(pos2)));
                                        pos_t rev2 = make_pos_t(get_id(pos2), 
                                                             !is_rev(pos2),
                                                             len2 - get_offset(pos2)-1); 
                                        size_t dist = dist_index.minimum_distance(get_id(pos1), get_is_rev(pos1), get_offset(pos1), get_id(pos2), get_is_rev(pos2), get_offset(pos2), true, &graph);
                                        if ( dist != -1 && dist <= read_lim) {
                                            new_clusters.union_groups(i1, i2);
                                        }

                                    }
                                }
                                auto actual_clusters = new_clusters.all_groups();
                                if (actual_clusters.size() != 1) {
                                    dist_index.print_self();
                                    graph.serialize("testGraph.hg");
                                    cerr << "These should be different read clusters: " << endl;
                                    for (auto c : actual_clusters) {
                                        cerr << "cluster: " ; 
                                        for (size_t i1 : c) {
                                            cerr << all_seeds[read_num][clust[i1]].pos << " ";
                                        }
                                        cerr << endl;
                                    }
                                }
                                REQUIRE(actual_clusters.size() == 1);
                            }
                        }
                    }
                    for (size_t a = 0; a < fragment_clusters.size(); a++) {
                        // For each cluster -cluster this cluster to ensure that 
                        // there is only one
                        vector<pos_t> clust = fragment_clusters[a];
                        
                        structures::UnionFind new_clusters (clust.size(), false);

                        for (size_t i1 = 0 ; i1 < clust.size() ; i1++) {
                            pos_t pos1 = clust[i1];
                            size_t len1 = graph.get_length(graph.get_handle(get_id(pos1), false));
                            pos_t rev1 = make_pos_t(get_id(pos1), 
                                                !is_rev(pos1),
                                                len1 - get_offset(pos1)-1); 

                            for (size_t b = 0 ; b < fragment_clusters.size() ; b++) {
                                if (b != a) {
                                    //For each other cluster
                                    vector<pos_t> clust2 = fragment_clusters[b];
                                    for (size_t i2 = 0 ; i2 < clust2.size() ; i2++) {
                                        //And each position in each other cluster,
                                        //make sure that this position is far away from i1
                                        pos_t pos2 = clust2[i2];
                                        size_t len2 = graph.get_length(graph.get_handle(get_id(pos2), false));
                                        pos_t rev2 = make_pos_t(get_id(pos2), 
                                                         !is_rev(pos2),
                                                         len2 - get_offset(pos2)-1); 

                                        size_t dist1 = dist_index.minimum_distance(get_id(pos1), get_is_rev(pos1), get_offset(pos1), get_id(pos2), get_is_rev(pos2), get_offset(pos2), false, &graph);
                                        size_t dist2 = dist_index.minimum_distance(get_id(pos1), get_is_rev(pos1), get_offset(pos1), get_id(rev2), get_is_rev(rev2), get_offset(rev2), false, &graph);
                                        size_t dist3 = dist_index.minimum_distance(get_id(rev1), get_is_rev(rev1), get_offset(rev1), get_id(pos2), get_is_rev(pos2), get_offset(pos2), false, &graph);
                                        size_t dist4 = dist_index.minimum_distance(get_id(rev1), get_is_rev(rev1), get_offset(rev1), get_id(rev2), get_is_rev(rev2), get_offset(rev2), false, &graph);
                                        size_t dist = std::min(std::min(dist1, dist2), std::min( dist3, dist4));
                                        if ( dist != -1 && dist <= fragment_lim) {
                                            dist_index.print_self();
                                            graph.serialize("testGraph.hg");
                                            cerr << "These should have been in the same fragment cluster: " ;
                                            cerr << pos1 << " and " << pos2 << endl;
                                            cerr << dist1 << " " << dist2 << " " << dist3 << " " << dist4 << endl;
                                            REQUIRE(false);
                                        }
                                        
                                    }
                                }
                            }
                            for (size_t i2 = 0 ; i2 < clust.size() ; i2++) {
                                //For each position in the same cluster
                                pos_t pos2 = clust[i2];
                                size_t len2 = graph.get_length(graph.get_handle(get_id(pos2), false));
                                pos_t rev2 = make_pos_t(get_id(pos2), 
                                                     !is_rev(pos2),
                                                     len2 - get_offset(pos2)-1); 
                                size_t dist1 = dist_index.minimum_distance(get_id(pos1), get_is_rev(pos1), get_offset(pos1), get_id(pos2), get_is_rev(pos2), get_offset(pos2), false, &graph);
                                size_t dist2 = dist_index.minimum_distance(get_id(pos1), get_is_rev(pos1), get_offset(pos1), get_id(rev2), get_is_rev(rev2), get_offset(rev2), false, &graph);
                                size_t dist3 = dist_index.minimum_distance(get_id(rev1), get_is_rev(rev1), get_offset(rev1), get_id(pos2), get_is_rev(pos2), get_offset(pos2), false, &graph);
                                size_t dist4 = dist_index.minimum_distance(get_id(rev1), get_is_rev(rev1), get_offset(rev1), get_id(rev2), get_is_rev(rev2), get_offset(rev2), false, &graph);
                                size_t dist = std::min(std::min(dist1, 
                                                   dist2), std::min( dist3, dist4));
                                if ( dist != -1 && dist <= fragment_lim) {
                                    new_clusters.union_groups(i1, i2);
                                }

                            }
                        }
                        auto actual_clusters = new_clusters.all_groups();
                        if (actual_clusters.size() != 1) {
                                            dist_index.print_self();
                            graph.serialize("testGraph.hg");
                            cerr << "These should be different fragment clusters: " << endl;
                            for (auto c : actual_clusters) {
                                cerr << "cluster: " ; 
                                for (size_t i1 : c) {
                                    cerr << clust[i1] << " ";
                                }
                                cerr << endl;
                            }
                        }
                        REQUIRE(actual_clusters.size() == 1);
                    }
                }
            }
        }
    } //end test case
}
}
