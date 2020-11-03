//
//  phased_genome.cpp
//
//  Unit tests for PhasedGenome object
//

#include <stdio.h>
#include <iostream>
#include <list>
#include "vg/io/json2pb.h"
#include <vg/vg.pb.h>
#include "catch.hpp"
#include "../utility.hpp"
#include "../phased_genome.hpp"
#include "../multipath_alignment.hpp"
#include "../nodetraversal.hpp"
#include "../vg.hpp"
#include "../genotypekit.hpp"
#include "../cactus_snarl_finder.hpp"

namespace vg {
    namespace unittest {
        
        TEST_CASE( "PhasedGenome can be built correctly",
                  "[phasing][mcmc]" ) {
            
            SECTION( "PhasedGenome constructor runs") {
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("G");
                Node* n4 = graph.create_node("CTGA");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n4);
                graph.create_edge(n3, n4);
                
                CactusSnarlFinder bubble_finder(graph);
                SnarlManager snarl_manager = bubble_finder.find_snarls();
                
                PhasedGenome genome(snarl_manager);
                
            }
            
            SECTION( "PhasedGenome can add haplotypes") {
                
                // construct graph
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("G");
                Node* n4 = graph.create_node("CTGA");
                Node* n5 = graph.create_node("TTG");
                Node* n6 = graph.create_node("CGGATA");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n4);
                graph.create_edge(n3, n4);
                graph.create_edge(n4, n5);
                graph.create_edge(n5, n6);
                // segregating inversion
                graph.create_edge(n4, n5, false, true);
                graph.create_edge(n5, n6, true, false);
                
                CactusSnarlFinder bubble_finder(graph);
                SnarlManager snarl_manager = bubble_finder.find_snarls();
                
                PhasedGenome genome(snarl_manager);
                
                // construct haplotypes
                
                list<NodeTraversal> haplotype_1;
                list<NodeTraversal> haplotype_2;
                
                haplotype_1.push_back(NodeTraversal(n1));
                haplotype_1.push_back(NodeTraversal(n2));
                haplotype_1.push_back(NodeTraversal(n4));
                haplotype_1.push_back(NodeTraversal(n5));
                haplotype_1.push_back(NodeTraversal(n6));
                
                haplotype_2.push_back(NodeTraversal(n1));
                haplotype_2.push_back(NodeTraversal(n3));
                haplotype_2.push_back(NodeTraversal(n4));
                haplotype_2.push_back(NodeTraversal(n5, true)); // reversed
                haplotype_2.push_back(NodeTraversal(n6));
                
                genome.add_haplotype(haplotype_1.begin(), haplotype_1.end());
                genome.add_haplotype(haplotype_2.begin(), haplotype_2.end());
                
            }
            
            SECTION( "PhasedGenome can iterate through haplotypes") {
                
                // construct graph
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("G");
                Node* n4 = graph.create_node("CTGA");
                Node* n5 = graph.create_node("TTG");
                Node* n6 = graph.create_node("CGGATA");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n4);
                graph.create_edge(n3, n4);
                graph.create_edge(n4, n5);
                graph.create_edge(n5, n6);
                // segregating inversion
                graph.create_edge(n4, n5, false, true);
                graph.create_edge(n5, n6, true, false);
                
                CactusSnarlFinder bubble_finder(graph);
                SnarlManager snarl_manager = bubble_finder.find_snarls();
                
                PhasedGenome genome(snarl_manager);
                
                // construct haplotypes
                
                list<NodeTraversal> haplotype_1;
                list<NodeTraversal> haplotype_2;
                
                haplotype_1.push_back(NodeTraversal(n1));
                haplotype_1.push_back(NodeTraversal(n2));
                haplotype_1.push_back(NodeTraversal(n4));
                haplotype_1.push_back(NodeTraversal(n5));
                haplotype_1.push_back(NodeTraversal(n6));
                
                haplotype_2.push_back(NodeTraversal(n1));
                haplotype_2.push_back(NodeTraversal(n3));
                haplotype_2.push_back(NodeTraversal(n4));
                haplotype_2.push_back(NodeTraversal(n5, true)); // reversed
                haplotype_2.push_back(NodeTraversal(n6));
                
                genome.add_haplotype(haplotype_1.begin(), haplotype_1.end());
                genome.add_haplotype(haplotype_2.begin(), haplotype_2.end());
                
                auto iter = genome.begin(0);
                REQUIRE((*iter) == NodeTraversal(n1, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n2, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n4, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n5, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n6, false));
                iter++;
                REQUIRE(iter == genome.end(0));
                
                iter = genome.begin(1);
                REQUIRE((*iter) == NodeTraversal(n1, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n3, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n4, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n5, true));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n6, false));
                iter++;
                REQUIRE(iter == genome.end(1));
                
            }
            
        }
        
        // TODO: temporarily disabling these tests because they're failing, but they're for an
        // experimental feature that isn't included in the VG paper
        
        TEST_CASE("PhasedGenome can be indexed and edited correctly",
                  "[phasing][mcmc]" ) {
        
            SECTION( "PhasedGenome can build site indices with simple sites") {
                
                // construct graph
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("G");
                Node* n4 = graph.create_node("CTGA");
                Node* n5 = graph.create_node("TTG");
                Node* n6 = graph.create_node("CGGATA");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n4);
                graph.create_edge(n3, n4);
                graph.create_edge(n4, n5);
                graph.create_edge(n5, n6);
                // segregating inversion
                graph.create_edge(n4, n5, false, true);
                graph.create_edge(n5, n6, true, false);
                
                CactusSnarlFinder bubble_finder(graph);
                SnarlManager snarl_manager = bubble_finder.find_snarls();
                
                PhasedGenome genome(snarl_manager);
                
                // construct haplotypes
                
                list<NodeTraversal> haplotype_1;
                list<NodeTraversal> haplotype_2;
                
                haplotype_1.push_back(NodeTraversal(n1));
                haplotype_1.push_back(NodeTraversal(n2));
                haplotype_1.push_back(NodeTraversal(n4));
                haplotype_1.push_back(NodeTraversal(n5));
                haplotype_1.push_back(NodeTraversal(n6));
                
                haplotype_2.push_back(NodeTraversal(n1));
                haplotype_2.push_back(NodeTraversal(n3));
                haplotype_2.push_back(NodeTraversal(n4));
                haplotype_2.push_back(NodeTraversal(n5, true)); // reversed
                haplotype_2.push_back(NodeTraversal(n6));
                
                genome.add_haplotype(haplotype_1.begin(), haplotype_1.end());
                genome.add_haplotype(haplotype_2.begin(), haplotype_2.end());
                
                // index sites
                
                genome.build_indices();
                
            }
            
            SECTION( "PhasedGenome can build indices with nested sites") {
                
                // construct graph
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("G");
                Node* n4 = graph.create_node("CTGA");
                Node* n5 = graph.create_node("TTG");
                Node* n6 = graph.create_node("CGGATA");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n6);
                graph.create_edge(n2, n3);
                graph.create_edge(n2, n4);
                graph.create_edge(n3, n5);
                graph.create_edge(n4, n5);
                graph.create_edge(n5, n6);
                
                CactusSnarlFinder bubble_finder(graph);
                SnarlManager snarl_manager = bubble_finder.find_snarls();
                
                PhasedGenome genome(snarl_manager);
                
                // construct haplotypes
                
                list<NodeTraversal> haplotype_1;
                list<NodeTraversal> haplotype_2;
                
                haplotype_1.push_back(NodeTraversal(n1));
                haplotype_1.push_back(NodeTraversal(n6));
                
                haplotype_2.push_back(NodeTraversal(n1));
                haplotype_2.push_back(NodeTraversal(n2));
                haplotype_2.push_back(NodeTraversal(n3));
                haplotype_2.push_back(NodeTraversal(n5));
                haplotype_2.push_back(NodeTraversal(n6));
                
                genome.add_haplotype(haplotype_1.begin(), haplotype_1.end());
                genome.add_haplotype(haplotype_2.begin(), haplotype_2.end());
                
                // index sites
                
                genome.build_indices();
                
            }
            
            SECTION( "PhasedGenome can set a simple allele") {
                
                // construct graph
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("G");
                Node* n4 = graph.create_node("CTGA");
                Node* n5 = graph.create_node("TTG");
                Node* n6 = graph.create_node("CGGATA");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n4);
                graph.create_edge(n3, n4);
                graph.create_edge(n4, n5);
                graph.create_edge(n5, n6);
                // segregating inversion
                graph.create_edge(n4, n5, false, true);
                graph.create_edge(n5, n6, true, false);
                
                CactusSnarlFinder bubble_finder(graph);
                SnarlManager snarl_manager = bubble_finder.find_snarls();
                
                PhasedGenome genome(snarl_manager);
                
                // construct haplotypes
                
                list<NodeTraversal> haplotype_1;
                list<NodeTraversal> haplotype_2;
                
                haplotype_1.push_back(NodeTraversal(n1));
                haplotype_1.push_back(NodeTraversal(n2));
                haplotype_1.push_back(NodeTraversal(n4));
                haplotype_1.push_back(NodeTraversal(n5));
                haplotype_1.push_back(NodeTraversal(n6));
                
                haplotype_2.push_back(NodeTraversal(n1));
                haplotype_2.push_back(NodeTraversal(n3));
                haplotype_2.push_back(NodeTraversal(n4));
                haplotype_2.push_back(NodeTraversal(n5, true)); // reversed
                haplotype_2.push_back(NodeTraversal(n6));
                
                genome.add_haplotype(haplotype_1.begin(), haplotype_1.end());
                genome.add_haplotype(haplotype_2.begin(), haplotype_2.end());
                
                // index sites
                
                genome.build_indices();
                
                list<NodeTraversal> allele;
                allele.push_back(NodeTraversal(n2));
                
                auto iter = genome.begin(1);
                REQUIRE((*iter) == NodeTraversal(n1, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n3, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n4, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n5, true));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n6, false));
                iter++;
                REQUIRE(iter == genome.end(1));
                
                // find the site we want to switch
                for (const Snarl* site : snarl_manager.top_level_snarls()) {
                    if (site->start().node_id() == n1->id() || site->end().node_id() == n1->id()) {
                        // set the allele
                        genome.set_allele(*site, allele.begin(), allele.end(), 1);
                        break;
                    }
                }
                
                iter = genome.begin(1);
                REQUIRE((*iter) == NodeTraversal(n1, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n2, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n4, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n5, true));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n6, false));
                iter++;
                REQUIRE(iter == genome.end(1));
            }
            
            SECTION( "PhasedGenome can retrieve an allele" ) {
                // construct graph
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("G");
                Node* n4 = graph.create_node("CTGA");
                Node* n5 = graph.create_node("TTG");
                Node* n6 = graph.create_node("CGGATA");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n4);
                graph.create_edge(n3, n4);
                graph.create_edge(n4, n5);
                graph.create_edge(n5, n6);
                // segregating inversion
                graph.create_edge(n4, n5, false, true);
                graph.create_edge(n5, n6, true, false);
                
                CactusSnarlFinder bubble_finder(graph);
                SnarlManager snarl_manager = bubble_finder.find_snarls();
                
                PhasedGenome genome(snarl_manager);
                
                // construct haplotypes
                
                list<NodeTraversal> haplotype_1;
                list<NodeTraversal> haplotype_2;
                
                haplotype_1.push_back(NodeTraversal(n1));
                haplotype_1.push_back(NodeTraversal(n2));
                haplotype_1.push_back(NodeTraversal(n4));
                haplotype_1.push_back(NodeTraversal(n5));
                haplotype_1.push_back(NodeTraversal(n6));
                
                haplotype_2.push_back(NodeTraversal(n1));
                haplotype_2.push_back(NodeTraversal(n3));
                haplotype_2.push_back(NodeTraversal(n4));
                haplotype_2.push_back(NodeTraversal(n5, true)); // reversed
                haplotype_2.push_back(NodeTraversal(n6));
                
                int haplo_id_1 = genome.add_haplotype(haplotype_1.begin(), haplotype_1.end());
                int haplo_id_2 = genome.add_haplotype(haplotype_2.begin(), haplotype_2.end());
                
                // index sites
                
                genome.build_indices();
                
                // find the site we want to retrieve
                for (const Snarl* site : snarl_manager.top_level_snarls()) {
                    if (site->start().node_id() == n1->id()) {
                        vector<NodeTraversal> allele_1 = genome.get_allele(*site, haplo_id_1);
                        vector<NodeTraversal> allele_2 = genome.get_allele(*site, haplo_id_2);
                        
                        REQUIRE(allele_1.size() == 1);
                        REQUIRE(allele_2.size() == 1);
                        
                        REQUIRE(allele_1[0] == NodeTraversal(n2, false));
                        REQUIRE(allele_2[0] == NodeTraversal(n3, false));
                    }
                    else if (site->end().node_id() == n1->id()) {
                        vector<NodeTraversal> allele_1 = genome.get_allele(*site, haplo_id_1);
                        vector<NodeTraversal> allele_2 = genome.get_allele(*site, haplo_id_2);
                        
                        REQUIRE(allele_1.size() == 1);
                        REQUIRE(allele_2.size() == 1);
                        
                        
                        REQUIRE(allele_1[0] == NodeTraversal(n2, true));
                        REQUIRE(allele_2[0] == NodeTraversal(n3, true));
                    }
                }
            }
            
            SECTION( "PhasedGenome can swap a simple allele between chromosomes") {
                
                // construct graph
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("G");
                Node* n4 = graph.create_node("CTGA");
                Node* n5 = graph.create_node("TTG");
                Node* n6 = graph.create_node("CGGATA");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n4);
                graph.create_edge(n3, n4);
                graph.create_edge(n4, n5);
                graph.create_edge(n5, n6);
                // segregating inversion
                graph.create_edge(n4, n5, false, true);
                graph.create_edge(n5, n6, true, false);
                
                CactusSnarlFinder bubble_finder(graph);
                SnarlManager snarl_manager = bubble_finder.find_snarls();
                
                PhasedGenome genome(snarl_manager);
                
                // construct haplotypes
                
                list<NodeTraversal> haplotype_1;
                list<NodeTraversal> haplotype_2;
                
                haplotype_1.push_back(NodeTraversal(n1));
                haplotype_1.push_back(NodeTraversal(n2));
                haplotype_1.push_back(NodeTraversal(n4));
                haplotype_1.push_back(NodeTraversal(n5));
                haplotype_1.push_back(NodeTraversal(n6));
                
                haplotype_2.push_back(NodeTraversal(n1));
                haplotype_2.push_back(NodeTraversal(n3));
                haplotype_2.push_back(NodeTraversal(n4));
                haplotype_2.push_back(NodeTraversal(n5, true)); // reversed
                haplotype_2.push_back(NodeTraversal(n6));
                
                genome.add_haplotype(haplotype_1.begin(), haplotype_1.end());
                genome.add_haplotype(haplotype_2.begin(), haplotype_2.end());
                
                // index sites
                
                genome.build_indices();
                
                auto iter = genome.begin(0);
                REQUIRE((*iter) == NodeTraversal(n1, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n2, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n4, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n5, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n6, false));
                iter++;
                REQUIRE(iter == genome.end(0));
                
                iter = genome.begin(1);
                REQUIRE((*iter) == NodeTraversal(n1, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n3, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n4, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n5, true));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n6, false));
                iter++;
                REQUIRE(iter == genome.end(1));
                
                // find the site we want to switch
                for (const Snarl* site : snarl_manager.top_level_snarls()) {
                    if (site->start().node_id() == n1->id() || site->end().node_id() == n1->id()) {
                        // swap the allele
                        genome.swap_alleles(*site, 0, 1);
                        break;
                    }
                }
                
                iter = genome.begin(0);
                REQUIRE((*iter) == NodeTraversal(n1, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n3, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n4, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n5, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n6, false));
                iter++;
                REQUIRE(iter == genome.end(0));
                
                iter = genome.begin(1);
                REQUIRE((*iter) == NodeTraversal(n1, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n2, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n4, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n5, true));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n6, false));
                iter++;
                REQUIRE(iter == genome.end(1));
            }
            
            SECTION( "PhasedGenome can perform successive allele swaps") {
                
                // construct graph
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("G");
                Node* n4 = graph.create_node("CTGA");
                Node* n5 = graph.create_node("TTG");
                Node* n6 = graph.create_node("CGGATA");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n4);
                graph.create_edge(n3, n4);
                graph.create_edge(n4, n5);
                graph.create_edge(n5, n6);
                // segregating inversion
                graph.create_edge(n4, n5, false, true);
                graph.create_edge(n5, n6, true, false);
                
                CactusSnarlFinder bubble_finder(graph);
                SnarlManager snarl_manager = bubble_finder.find_snarls();
                
                PhasedGenome genome(snarl_manager);
                
                // add haplotypes
                
                list<NodeTraversal> haplotype_1;
                list<NodeTraversal> haplotype_2;
                
                haplotype_1.push_back(NodeTraversal(n1));
                haplotype_1.push_back(NodeTraversal(n2));
                haplotype_1.push_back(NodeTraversal(n4));
                haplotype_1.push_back(NodeTraversal(n5));
                haplotype_1.push_back(NodeTraversal(n6));
                
                haplotype_2.push_back(NodeTraversal(n1));
                haplotype_2.push_back(NodeTraversal(n3));
                haplotype_2.push_back(NodeTraversal(n4));
                haplotype_2.push_back(NodeTraversal(n5, true)); // reversed
                haplotype_2.push_back(NodeTraversal(n6));
                
                genome.add_haplotype(haplotype_1.begin(), haplotype_1.end());
                genome.add_haplotype(haplotype_2.begin(), haplotype_2.end());
                
                // index sites
                
                genome.build_indices();
                
                // find the site we want to switch
                for (const Snarl* site : snarl_manager.top_level_snarls()) {
                    if (site->start().node_id() == n1->id() || site->end().node_id() == n1->id()) {
                        // swap the allele twice
                        genome.swap_alleles(*site, 0, 1);
                        genome.swap_alleles(*site, 0, 1);
                        break;
                    }
                }
                
                auto iter = genome.begin(0);
                REQUIRE((*iter) == NodeTraversal(n1, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n2, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n4, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n5, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n6, false));
                iter++;
                REQUIRE(iter == genome.end(0));
                
                iter = genome.begin(1);
                REQUIRE((*iter) == NodeTraversal(n1, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n3, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n4, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n5, true));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n6, false));
                iter++;
                REQUIRE(iter == genome.end(1));

            }
            
            SECTION( "PhasedGenome can perform successive allele swaps involving a deletion") {
                
                // construct graph
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("G");
                Node* n4 = graph.create_node("CTGA");
                Node* n5 = graph.create_node("TTG");
                Node* n6 = graph.create_node("CGGATA");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n1, n4);
                graph.create_edge(n2, n4);
                graph.create_edge(n3, n4);
                graph.create_edge(n4, n5);
                graph.create_edge(n5, n6);
                // segregating inversion
                graph.create_edge(n4, n5, false, true);
                graph.create_edge(n5, n6, true, false);
                
                CactusSnarlFinder bubble_finder(graph);
                SnarlManager snarl_manager = bubble_finder.find_snarls();
                
                PhasedGenome genome(snarl_manager);
                
                list<NodeTraversal> haplotype_1;
                list<NodeTraversal> haplotype_2;
                
                haplotype_1.push_back(NodeTraversal(n1));
                haplotype_1.push_back(NodeTraversal(n4));
                haplotype_1.push_back(NodeTraversal(n5));
                haplotype_1.push_back(NodeTraversal(n6));
                
                haplotype_2.push_back(NodeTraversal(n1));
                haplotype_2.push_back(NodeTraversal(n3));
                haplotype_2.push_back(NodeTraversal(n4));
                haplotype_2.push_back(NodeTraversal(n5, true)); // reversed
                haplotype_2.push_back(NodeTraversal(n6));
                
                genome.add_haplotype(haplotype_1.begin(), haplotype_1.end());
                genome.add_haplotype(haplotype_2.begin(), haplotype_2.end());
                
                // index sites
                
                genome.build_indices();
                
                // find the site we want to switch
                for (const Snarl* site : snarl_manager.top_level_snarls()) {
                    if (site->start().node_id() == n1->id() || site->end().node_id() == n1->id()) {
                        // swap the allele twice
                        genome.swap_alleles(*site, 0, 1);
                        genome.swap_alleles(*site, 0, 1);
                        break;
                    }
                }
                
                auto iter = genome.begin(0);
                REQUIRE((*iter) == NodeTraversal(n1, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n4, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n5, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n6, false));
                iter++;
                REQUIRE(iter == genome.end(0));
                
                iter = genome.begin(1);
                REQUIRE((*iter) == NodeTraversal(n1, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n3, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n4, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n5, true));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n6, false));
                iter++;
                REQUIRE(iter == genome.end(1));
                
            }
            
            SECTION( "PhasedGenome can perform successive allele edits in a nested allele") {
                
                // construct graph
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("G");
                Node* n4 = graph.create_node("CTGA");
                Node* n5 = graph.create_node("TTG");
                Node* n6 = graph.create_node("CGGATA");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n6);
                graph.create_edge(n2, n3);
                graph.create_edge(n2, n4);
                graph.create_edge(n3, n5);
                graph.create_edge(n4, n5);
                graph.create_edge(n5, n6);
                
                CactusSnarlFinder bubble_finder(graph);
                SnarlManager snarl_manager = bubble_finder.find_snarls();
                
                PhasedGenome genome(snarl_manager);
                
                // add haplotypes
                
                list<NodeTraversal> haplotype_1;
                list<NodeTraversal> haplotype_2;
                
                haplotype_1.push_back(NodeTraversal(n1));
                haplotype_1.push_back(NodeTraversal(n2));
                haplotype_1.push_back(NodeTraversal(n3));
                haplotype_1.push_back(NodeTraversal(n5));
                haplotype_1.push_back(NodeTraversal(n6));
                
                haplotype_2.push_back(NodeTraversal(n1));
                haplotype_2.push_back(NodeTraversal(n6));
                
                genome.add_haplotype(haplotype_1.begin(), haplotype_1.end());
                genome.add_haplotype(haplotype_2.begin(), haplotype_2.end());
                
                // index sites
                
                genome.build_indices();
                
                list<NodeTraversal> allele;
                allele.push_back(NodeTraversal(n4));
                
                // find the site we want to switch
                REQUIRE(snarl_manager.top_level_snarls().size() > 0);
                const Snarl* site = snarl_manager.top_level_snarls()[0];
                REQUIRE(snarl_manager.children_of(site).size() > 0);
                const Snarl* subsite = snarl_manager.children_of(site)[0];
                
                genome.swap_alleles(*site, 0, 1);
                genome.set_allele(*subsite, allele.begin(), allele.end(), 1);
                
                auto iter = genome.begin(0);
                REQUIRE((*iter) == NodeTraversal(n1, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n6, false));
                iter++;
                REQUIRE(iter == genome.end(0));
                
                iter = genome.begin(1);
                REQUIRE((*iter) == NodeTraversal(n1, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n2, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n4, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n5, false));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n6, false));
                iter++;
                REQUIRE(iter == genome.end(1));
            }
        }
        
        
        TEST_CASE("PhasedGenome can perform computations on multipath alignments",
                  "[phasing][mcmc]" ) {
            
            SECTION( "PhasedGenome can compute the score of a restricted multipath alignment on one haplotype") {
                
                // construct graph
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("CTGT");
                Node* n3 = graph.create_node("G");
                Node* n4 = graph.create_node("A");
                Node* n5 = graph.create_node("TTG");
                Node* n6 = graph.create_node("CGGATA");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n6);
                graph.create_edge(n2, n3);
                graph.create_edge(n2, n4);
                graph.create_edge(n3, n5);
                graph.create_edge(n4, n5);
                graph.create_edge(n5, n6);
                
                CactusSnarlFinder bubble_finder(graph);
                SnarlManager snarl_manager = bubble_finder.find_snarls();
                
                PhasedGenome genome(snarl_manager);
                
                REQUIRE(snarl_manager.top_level_snarls().size() > 0);
                const Snarl* site = snarl_manager.top_level_snarls()[0];
                REQUIRE(snarl_manager.children_of(site).size() > 0);
                const Snarl* subsite = snarl_manager.children_of(site)[0];
                
                // add haplotypes
                
                list<NodeTraversal> haplotype_1;
                
                haplotype_1.push_back(NodeTraversal(n1));
                haplotype_1.push_back(NodeTraversal(n2));
                haplotype_1.push_back(NodeTraversal(n3));
                haplotype_1.push_back(NodeTraversal(n5));
                haplotype_1.push_back(NodeTraversal(n6));
                
                genome.add_haplotype(haplotype_1.begin(), haplotype_1.end());
                
                // index sites
                
                genome.build_indices();
                
                // construct a multipath alignment
                
                multipath_alignment_t multipath_aln;
                multipath_aln.set_sequence("CACTGTATTGCGGATA");
                
                // add subpaths
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                subpath_t* subpath0 = multipath_aln.mutable_subpath(0);
                subpath_t* subpath1 = multipath_aln.mutable_subpath(1);
                subpath_t* subpath2 = multipath_aln.mutable_subpath(2);
                subpath_t* subpath3 = multipath_aln.mutable_subpath(3);
                
                // set edges between subpaths
                subpath0->add_next(1);
                subpath0->add_next(2);
                subpath1->add_next(3);
                subpath2->add_next(3);
                
                // designate mappings and scores
                path_mapping_t* mapping0 = subpath0->mutable_path()->add_mapping();
                mapping0->mutable_position()->set_node_id(1);
                mapping0->mutable_position()->set_offset(1);
                edit_t* edit00 = mapping0->add_edit();
                edit00->set_from_length(2);
                edit00->set_to_length(2);
                
                path_mapping_t* mapping1 = subpath0->mutable_path()->add_mapping();
                mapping1->mutable_position()->set_node_id(2);
                edit_t* edit10 = mapping1->add_edit();
                edit10->set_from_length(4);
                edit10->set_to_length(4);
                
                subpath0->set_score(6);
                
                path_mapping_t* mapping2 = subpath1->mutable_path()->add_mapping();
                mapping2->mutable_position()->set_node_id(3);
                edit_t* edit20 = mapping2->add_edit();
                edit20->set_from_length(1);
                edit20->set_to_length(1);
                edit20->set_sequence("A");
                
                subpath1->set_score(-4);
                
                path_mapping_t* mapping3 = subpath2->mutable_path()->add_mapping();
                mapping3->mutable_position()->set_node_id(4);
                edit_t* edit30 = mapping3->add_edit();
                edit30->set_from_length(1);
                edit30->set_to_length(1);
                
                subpath2->set_score(1);
                
                path_mapping_t* mapping4 = subpath3->mutable_path()->add_mapping();
                mapping4->mutable_position()->set_node_id(5);
                edit_t* edit40 = mapping4->add_edit();
                edit40->set_from_length(3);
                edit40->set_to_length(3);
                
                path_mapping_t* mapping5 = subpath3->mutable_path()->add_mapping();
                mapping5->mutable_position()->set_node_id(6);
                edit_t* edit50 = mapping5->add_edit();
                edit50->set_from_length(6);
                edit50->set_to_length(6);
                
                subpath3->set_score(9);
                
                identify_start_subpaths(multipath_aln);
                
                REQUIRE( genome.optimal_score_on_genome(multipath_aln, graph) == 6 - 4 + 9 );
                
                list<NodeTraversal> allele;
                allele.push_back(NodeTraversal(n4));
                
                genome.set_allele(*subsite, allele.begin(), allele.end(), 0);
                
                REQUIRE( genome.optimal_score_on_genome(multipath_aln, graph) == 6 + 1 + 9 );
            }
            
            
            
            SECTION( "PhasedGenome can compute the score of a restricted multipath alignment on two haplotype") {
                
                // construct graph
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("CTGT");
                Node* n3 = graph.create_node("G");
                Node* n4 = graph.create_node("A");
                Node* n5 = graph.create_node("TTG");
                Node* n6 = graph.create_node("CGGATA");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n6);
                graph.create_edge(n2, n3);
                graph.create_edge(n2, n4);
                graph.create_edge(n3, n5);
                graph.create_edge(n4, n5);
                graph.create_edge(n5, n6);
                
                // construct phased genome
                
                CactusSnarlFinder bubble_finder(graph);
                SnarlManager snarl_manager = bubble_finder.find_snarls();
                
                PhasedGenome genome(snarl_manager);
                
                REQUIRE(snarl_manager.top_level_snarls().size() > 0);
                const Snarl* site = snarl_manager.top_level_snarls()[0];
                REQUIRE(snarl_manager.children_of(site).size() > 0);
                const Snarl* subsite = snarl_manager.children_of(site)[0];
                
                list<NodeTraversal> haplotype_1;
                list<NodeTraversal> haplotype_2;
                
                haplotype_1.push_back(NodeTraversal(n1));
                haplotype_1.push_back(NodeTraversal(n2));
                haplotype_1.push_back(NodeTraversal(n3));
                haplotype_1.push_back(NodeTraversal(n5));
                haplotype_1.push_back(NodeTraversal(n6));
                
                haplotype_2.push_back(NodeTraversal(n1));
                haplotype_2.push_back(NodeTraversal(n2));
                haplotype_2.push_back(NodeTraversal(n4));
                haplotype_2.push_back(NodeTraversal(n5));
                haplotype_2.push_back(NodeTraversal(n6));
                
                genome.add_haplotype(haplotype_1.begin(), haplotype_1.end());
                genome.add_haplotype(haplotype_2.begin(), haplotype_2.end());
                
                // index sites
                
                genome.build_indices();
                
                // construct a multipath alignment
                
                multipath_alignment_t multipath_aln;
                multipath_aln.set_sequence("CACTGTATTGCGGATA");
                
                // add subpaths
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                subpath_t* subpath0 = multipath_aln.mutable_subpath(0);
                subpath_t* subpath1 = multipath_aln.mutable_subpath(1);
                subpath_t* subpath2 = multipath_aln.mutable_subpath(2);
                subpath_t* subpath3 = multipath_aln.mutable_subpath(3);
                
                // set edges between subpaths
                subpath0->add_next(1);
                subpath0->add_next(2);
                subpath1->add_next(3);
                subpath2->add_next(3);
                
                // designate mappings and scores
                path_mapping_t* mapping0 = subpath0->mutable_path()->add_mapping();
                mapping0->mutable_position()->set_node_id(1);
                mapping0->mutable_position()->set_offset(1);
                edit_t* edit00 = mapping0->add_edit();
                edit00->set_from_length(2);
                edit00->set_to_length(2);
                
                path_mapping_t* mapping1 = subpath0->mutable_path()->add_mapping();
                mapping1->mutable_position()->set_node_id(2);
                edit_t* edit10 = mapping1->add_edit();
                edit10->set_from_length(4);
                edit10->set_to_length(4);
                
                subpath0->set_score(6);
                
                path_mapping_t* mapping2 = subpath1->mutable_path()->add_mapping();
                mapping2->mutable_position()->set_node_id(3);
                edit_t* edit20 = mapping2->add_edit();
                edit20->set_from_length(1);
                edit20->set_to_length(1);
                edit20->set_sequence("A");
                
                subpath1->set_score(-4);
                
                path_mapping_t* mapping3 = subpath2->mutable_path()->add_mapping();
                mapping3->mutable_position()->set_node_id(4);
                edit_t* edit30 = mapping3->add_edit();
                edit30->set_from_length(1);
                edit30->set_to_length(1);
                
                subpath2->set_score(1);
                
                path_mapping_t* mapping4 = subpath3->mutable_path()->add_mapping();
                mapping4->mutable_position()->set_node_id(5);
                edit_t* edit40 = mapping4->add_edit();
                edit40->set_from_length(3);
                edit40->set_to_length(3);
                
                path_mapping_t* mapping5 = subpath3->mutable_path()->add_mapping();
                mapping5->mutable_position()->set_node_id(6);
                edit_t* edit50 = mapping5->add_edit();
                edit50->set_from_length(6);
                edit50->set_to_length(6);
                
                subpath3->set_score(9);
                
                identify_start_subpaths(multipath_aln);
                
                REQUIRE( genome.optimal_score_on_genome(multipath_aln, graph) == 6 + 1 + 9 );
            }
            
            SECTION( "PhasedGenome can compute the score of a restricted multipath alignment on the reverse strand") {
                
                // construct graph
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("CTGT");
                Node* n3 = graph.create_node("G");
                Node* n4 = graph.create_node("A");
                Node* n5 = graph.create_node("TTG");
                Node* n6 = graph.create_node("CGGATA");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n5, false, true);
                graph.create_edge(n2, n3);
                graph.create_edge(n2, n4);
                graph.create_edge(n2, n6, true, false);
                graph.create_edge(n3, n5);
                graph.create_edge(n4, n5);
                graph.create_edge(n5, n6);
                
                CactusSnarlFinder bubble_finder(graph);
                SnarlManager snarl_manager = bubble_finder.find_snarls();
                
                PhasedGenome genome(snarl_manager);
                
                // constuct haplotypes
                
                list<NodeTraversal> haplotype_1;
                
                haplotype_1.push_back(NodeTraversal(n1));
                haplotype_1.push_back(NodeTraversal(n5, true));
                haplotype_1.push_back(NodeTraversal(n3, true));
                haplotype_1.push_back(NodeTraversal(n2, true));
                haplotype_1.push_back(NodeTraversal(n6));
                
                genome.add_haplotype(haplotype_1.begin(), haplotype_1.end());
                
                // index sites
                
                genome.build_indices();
                
                // construct a multipath alignment
                
                multipath_alignment_t multipath_aln;
                multipath_aln.set_sequence("CCGCTGTATTGTGC");
                
                // add subpaths
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                subpath_t* subpath0 = multipath_aln.mutable_subpath(0);
                subpath_t* subpath1 = multipath_aln.mutable_subpath(1);
                subpath_t* subpath2 = multipath_aln.mutable_subpath(2);
                subpath_t* subpath3 = multipath_aln.mutable_subpath(3);
                
                // set edges between subpaths
                subpath0->add_next(1);
                subpath0->add_next(2);
                subpath1->add_next(3);
                subpath2->add_next(3);
                
                // designate mappings and scores
                path_mapping_t* mapping0 = subpath0->mutable_path()->add_mapping();
                mapping0->mutable_position()->set_node_id(6);
                mapping0->mutable_position()->set_offset(3);
                mapping0->mutable_position()->set_is_reverse(true);
                edit_t* edit00 = mapping0->add_edit();
                edit00->set_from_length(3);
                edit00->set_to_length(3);
                
                path_mapping_t* mapping1 = subpath0->mutable_path()->add_mapping();
                mapping1->mutable_position()->set_node_id(2);
                edit_t* edit10 = mapping1->add_edit();
                edit10->set_from_length(4);
                edit10->set_to_length(4);
                
                subpath0->set_score(7);
                
                path_mapping_t* mapping2 = subpath1->mutable_path()->add_mapping();
                mapping2->mutable_position()->set_node_id(3);
                edit_t* edit20 = mapping2->add_edit();
                edit20->set_from_length(1);
                edit20->set_to_length(1);
                edit20->set_sequence("A");
                
                subpath1->set_score(-4);
                
                path_mapping_t* mapping3 = subpath2->mutable_path()->add_mapping();
                mapping3->mutable_position()->set_node_id(4);
                edit_t* edit30 = mapping3->add_edit();
                edit30->set_from_length(1);
                edit30->set_to_length(1);
                
                subpath2->set_score(1);
                
                path_mapping_t* mapping4 = subpath3->mutable_path()->add_mapping();
                mapping4->mutable_position()->set_node_id(5);
                edit_t* edit40 = mapping4->add_edit();
                edit40->set_from_length(3);
                edit40->set_to_length(3);
                
                path_mapping_t* mapping5 = subpath3->mutable_path()->add_mapping();
                mapping5->mutable_position()->set_node_id(1);
                mapping5->mutable_position()->set_is_reverse(true);
                edit_t* edit50 = mapping5->add_edit();
                edit50->set_from_length(3);
                edit50->set_to_length(3);
                
                subpath3->set_score(6);
                
                identify_start_subpaths(multipath_aln);
                
                REQUIRE( genome.optimal_score_on_genome(multipath_aln, graph) == 7 - 4 + 6 );
                
            }
        }
    
    
        TEST_CASE("PhasedGenome can compute alignment likelihoods",
                  "[phasing][mcmc]" ) {
            
            VG graph;
            
            Node* n0 = graph.create_node("CCT");
            Node* n1 = graph.create_node("GCA");
            Node* n2 = graph.create_node("A");
            Node* n3 = graph.create_node("T");
            Node* n4 = graph.create_node("GAA");
            Node* n5 = graph.create_node("TAT");
            
            graph.create_edge(n0, n1);
            graph.create_edge(n1, n2);
            graph.create_edge(n1, n3);
            graph.create_edge(n2, n4);
            graph.create_edge(n3, n4);
            graph.create_edge(n4, n5);
            
            CactusSnarlFinder bubble_finder(graph);
            SnarlManager snarl_manager = bubble_finder.find_snarls();
            
            PhasedGenome genome(snarl_manager);
            
            SECTION("The likelihood is correct for a simple alignment") {
                
                
                // constuct haplotypes
                
                list<NodeTraversal> haplotype_1;
                
                haplotype_1.push_back(NodeTraversal(n1));
                haplotype_1.push_back(NodeTraversal(n2));
                haplotype_1.push_back(NodeTraversal(n4));
                
                genome.add_haplotype(haplotype_1.begin(), haplotype_1.end());
                
                // index sites
                
                genome.build_indices();
                
                multipath_alignment_t multipath_aln;
                multipath_aln.set_sequence("GCA");
                multipath_aln.set_mapping_quality(60);
                
                subpath_t* subpath = multipath_aln.add_subpath();
                path_t* path = subpath->mutable_path();
                path_mapping_t* mapping = path->add_mapping();
                position_t* position = mapping->mutable_position();
                position->set_node_id(n1->id());
                position->set_is_reverse(false);
                edit_t* edit = mapping->add_edit();
                edit->set_from_length(3);
                edit->set_to_length(3);
                subpath->set_score(3);
                
                identify_start_subpaths(multipath_aln);
                
                //
                REQUIRE(abs(genome.read_log_likelihood(multipath_aln, 1.0) - 3.0) < .001);
                REQUIRE(abs(genome.read_log_likelihood(multipath_aln, 2.0) - 6.0) < .001);
            }
            
            SECTION("The likelihood doubles if it is present on two haplotypes") {
                
                
                // constuct haplotypes
                
                list<NodeTraversal> haplotype_1;
                
                haplotype_1.push_back(NodeTraversal(n1));
                haplotype_1.push_back(NodeTraversal(n2));
                haplotype_1.push_back(NodeTraversal(n4));
                
                genome.add_haplotype(haplotype_1.begin(), haplotype_1.end());
                
                list<NodeTraversal> haplotype_2;
                
                haplotype_2.push_back(NodeTraversal(n1));
                haplotype_2.push_back(NodeTraversal(n3));
                haplotype_2.push_back(NodeTraversal(n4));
                
                genome.add_haplotype(haplotype_2.begin(), haplotype_2.end());
                
                // index sites
                
                genome.build_indices();
                
                multipath_alignment_t multipath_aln;
                multipath_aln.set_sequence("GCA");
                multipath_aln.set_mapping_quality(60);
                
                subpath_t* subpath = multipath_aln.add_subpath();
                path_t* path = subpath->mutable_path();
                path_mapping_t* mapping = path->add_mapping();
                position_t* position = mapping->mutable_position();
                position->set_node_id(n1->id());
                position->set_is_reverse(false);
                edit_t* edit = mapping->add_edit();
                edit->set_from_length(3);
                edit->set_to_length(3);
                subpath->set_score(3);
                
                identify_start_subpaths(multipath_aln);
                
                REQUIRE(abs(genome.read_log_likelihood(multipath_aln, 1.0) - (3.0 + log(2.0))) < .001);
            }
            
            SECTION("Likelihood can be calculated with a branching multipath alignment") {
                
                // constuct haplotypes
                
                list<NodeTraversal> haplotype_1;
                
                haplotype_1.push_back(NodeTraversal(n1));
                haplotype_1.push_back(NodeTraversal(n2));
                haplotype_1.push_back(NodeTraversal(n4));
                
                genome.add_haplotype(haplotype_1.begin(), haplotype_1.end());
                
                // index sites
                
                genome.build_indices();
                
                multipath_alignment_t multipath_aln;
                multipath_aln.set_sequence("GCAAGAA");
                multipath_aln.set_mapping_quality(60);
                
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                subpath_t* subpath0 = multipath_aln.mutable_subpath(0);
                subpath_t* subpath1 = multipath_aln.mutable_subpath(1);
                subpath_t* subpath2 = multipath_aln.mutable_subpath(2);
                subpath_t* subpath3 = multipath_aln.mutable_subpath(3);
                
                path_t* path0 = subpath0->mutable_path();
                path_mapping_t* mapping0 = path0->add_mapping();
                position_t* position0 = mapping0->mutable_position();
                position0->set_node_id(n1->id());
                position0->set_is_reverse(false);
                edit_t* edit0 = mapping0->add_edit();
                edit0->set_from_length(3);
                edit0->set_to_length(3);
                subpath0->set_score(3);
                
                path_t* path1 = subpath1->mutable_path();
                path_mapping_t* mapping1 = path1->add_mapping();
                position_t* position1 = mapping1->mutable_position();
                position1->set_node_id(n2->id());
                position1->set_is_reverse(false);
                edit_t* edit1 = mapping1->add_edit();
                edit1->set_from_length(1);
                edit1->set_to_length(1);
                subpath1->set_score(1);
                
                path_t* path2 = subpath2->mutable_path();
                path_mapping_t* mapping2 = path2->add_mapping();
                position_t* position2 = mapping2->mutable_position();
                position2->set_node_id(n3->id());
                position2->set_is_reverse(false);
                edit_t* edit2 = mapping2->add_edit();
                edit2->set_from_length(1);
                edit2->set_to_length(1);
                subpath2->set_score(-2);
                
                path_t* path3 = subpath3->mutable_path();
                path_mapping_t* mapping3 = path3->add_mapping();
                position_t* position3 = mapping3->mutable_position();
                position3->set_node_id(n4->id());
                position3->set_is_reverse(false);
                edit_t* edit3 = mapping3->add_edit();
                edit3->set_from_length(3);
                edit3->set_to_length(3);
                subpath3->set_score(3);
                
                subpath0->add_next(1);
                subpath0->add_next(2);
                subpath1->add_next(3);
                subpath2->add_next(3);
                
                identify_start_subpaths(multipath_aln);
                
                REQUIRE(abs(genome.read_log_likelihood(multipath_aln, 1.0) - 7.0) < .001);
            }
            
            SECTION("Log likeihood can be calculated for a read on the reverse strand") {
                // constuct haplotypes
                
                list<NodeTraversal> haplotype_1;
                
                haplotype_1.push_back(NodeTraversal(n1));
                haplotype_1.push_back(NodeTraversal(n3));
                haplotype_1.push_back(NodeTraversal(n4));
                
                genome.add_haplotype(haplotype_1.begin(), haplotype_1.end());
                
                // index sites
                
                genome.build_indices();
                
                multipath_alignment_t multipath_aln;
                multipath_aln.set_sequence("TTCTTGC");
                multipath_aln.set_mapping_quality(60);
                
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                subpath_t* subpath0 = multipath_aln.mutable_subpath(0);
                subpath_t* subpath1 = multipath_aln.mutable_subpath(1);
                subpath_t* subpath2 = multipath_aln.mutable_subpath(2);
                subpath_t* subpath3 = multipath_aln.mutable_subpath(3);
                
                path_t* path0 = subpath0->mutable_path();
                path_mapping_t* mapping0 = path0->add_mapping();
                position_t* position0 = mapping0->mutable_position();
                position0->set_node_id(n4->id());
                position0->set_is_reverse(true);
                edit_t* edit0 = mapping0->add_edit();
                edit0->set_from_length(3);
                edit0->set_to_length(3);
                subpath0->set_score(3);
                
                path_t* path1 = subpath1->mutable_path();
                path_mapping_t* mapping1 = path1->add_mapping();
                position_t* position1 = mapping1->mutable_position();
                position1->set_node_id(n2->id());
                position1->set_is_reverse(true);
                edit_t* edit1 = mapping1->add_edit();
                edit1->set_from_length(1);
                edit1->set_to_length(1);
                subpath1->set_score(1);
                
                path_t* path2 = subpath2->mutable_path();
                path_mapping_t* mapping2 = path2->add_mapping();
                position_t* position2 = mapping2->mutable_position();
                position2->set_node_id(n3->id());
                position2->set_is_reverse(true);
                edit_t* edit2 = mapping2->add_edit();
                edit2->set_from_length(1);
                edit2->set_to_length(1);
                subpath2->set_score(-2);
                
                path_t* path3 = subpath3->mutable_path();
                path_mapping_t* mapping3 = path3->add_mapping();
                position_t* position3 = mapping3->mutable_position();
                position3->set_node_id(n1->id());
                position3->set_is_reverse(true);
                edit_t* edit3 = mapping3->add_edit();
                edit3->set_from_length(3);
                edit3->set_to_length(3);
                subpath3->set_score(3);
                
                subpath0->add_next(1);
                subpath0->add_next(2);
                subpath1->add_next(3);
                subpath2->add_next(3);
                
                identify_start_subpaths(multipath_aln);
                
                REQUIRE(abs(genome.read_log_likelihood(multipath_aln, 1.0) - 4.0) < .001);
            }
            
            SECTION("Log likeihood can be calculated for a read with mid-node subpath adjacencies") {
                
                // constuct haplotypes
                
                list<NodeTraversal> haplotype_1;
                
                haplotype_1.push_back(NodeTraversal(n1));
                haplotype_1.push_back(NodeTraversal(n3));
                haplotype_1.push_back(NodeTraversal(n4));
                
                genome.add_haplotype(haplotype_1.begin(), haplotype_1.end());
                
                // index sites
                
                genome.build_indices();
                
                multipath_alignment_t multipath_aln;
                multipath_aln.set_sequence("TTC");
                multipath_aln.set_mapping_quality(60);
                
                
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                subpath_t* subpath0 = multipath_aln.mutable_subpath(0);
                subpath_t* subpath1 = multipath_aln.mutable_subpath(1);
                
                path_t* path0 = subpath0->mutable_path();
                path_mapping_t* mapping0 = path0->add_mapping();
                position_t* position0 = mapping0->mutable_position();
                position0->set_node_id(n4->id());
                position0->set_is_reverse(true);
                edit_t* edit0 = mapping0->add_edit();
                edit0->set_from_length(2);
                edit0->set_to_length(2);
                subpath0->set_score(2);
                
                path_t* path1 = subpath1->mutable_path();
                path_mapping_t* mapping1 = path1->add_mapping();
                position_t* position1 = mapping1->mutable_position();
                position1->set_node_id(n4->id());
                position1->set_is_reverse(true);
                position1->set_offset(2);
                edit_t* edit1 = mapping1->add_edit();
                edit1->set_from_length(1);
                edit1->set_to_length(1);
                subpath1->set_score(1);
                
                subpath0->add_next(1);
                                
                identify_start_subpaths(multipath_aln);
                
                REQUIRE(abs(genome.read_log_likelihood(multipath_aln, 1.0) - 3.0) < .001);
            }
            
            SECTION("Log likeihood can be calculated for a read multiple with mappings per subpath") {
                
                // constuct haplotypes
                
                list<NodeTraversal> haplotype_1;
                
                haplotype_1.push_back(NodeTraversal(n0));
                haplotype_1.push_back(NodeTraversal(n1));
                haplotype_1.push_back(NodeTraversal(n2));
                haplotype_1.push_back(NodeTraversal(n4));
                haplotype_1.push_back(NodeTraversal(n5));
                
                genome.add_haplotype(haplotype_1.begin(), haplotype_1.end());
                
                // index sites
                
                genome.build_indices();
                
                multipath_alignment_t multipath_aln;
                multipath_aln.set_sequence("CCTGCAAGAATAT");
                multipath_aln.set_mapping_quality(60);
                
                
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                subpath_t* subpath0 = multipath_aln.mutable_subpath(0);
                subpath_t* subpath1 = multipath_aln.mutable_subpath(1);
                subpath_t* subpath2 = multipath_aln.mutable_subpath(2);
                subpath_t* subpath3 = multipath_aln.mutable_subpath(3);
                
                path_t* path0 = subpath0->mutable_path();
                path0->add_mapping();
                path0->add_mapping();
                path_mapping_t* mapping0 = path0->mutable_mapping(0);
                path_mapping_t* mapping4 = path0->mutable_mapping(1);
                
                position_t* position0 = mapping0->mutable_position();
                position0->set_node_id(n0->id());
                position0->set_is_reverse(false);
                edit_t* edit0 = mapping0->add_edit();
                edit0->set_from_length(3);
                edit0->set_to_length(3);
                position_t* position4 = mapping4->mutable_position();
                position4->set_node_id(n1->id());
                position4->set_is_reverse(false);
                edit_t* edit4 = mapping4->add_edit();
                edit4->set_from_length(3);
                edit4->set_to_length(3);
                subpath0->set_score(6);
                
                path_t* path1 = subpath1->mutable_path();
                path_mapping_t* mapping1 = path1->add_mapping();
                position_t* position1 = mapping1->mutable_position();
                position1->set_node_id(n2->id());
                position1->set_is_reverse(false);
                edit_t* edit1 = mapping1->add_edit();
                edit1->set_from_length(1);
                edit1->set_to_length(1);
                subpath1->set_score(1);
                
                path_t* path2 = subpath2->mutable_path();
                path_mapping_t* mapping2 = path2->add_mapping();
                position_t* position2 = mapping2->mutable_position();
                position2->set_node_id(n3->id());
                position2->set_is_reverse(false);
                edit_t* edit2 = mapping2->add_edit();
                edit2->set_from_length(1);
                edit2->set_to_length(1);
                subpath2->set_score(-2);
                
                
                path_t* path3 = subpath3->mutable_path();
                path3->add_mapping();
                path3->add_mapping();
                path_mapping_t* mapping3 = path3->mutable_mapping(0);
                path_mapping_t* mapping5 = path3->mutable_mapping(1);
                
                position_t* position3 = mapping3->mutable_position();
                position3->set_node_id(n4->id());
                position3->set_is_reverse(false);
                edit_t* edit3 = mapping3->add_edit();
                edit3->set_from_length(3);
                edit3->set_to_length(3);
                position_t* position5 = mapping5->mutable_position();
                position5->set_node_id(n5->id());
                position5->set_is_reverse(false);
                edit_t* edit5 = mapping5->add_edit();
                edit5->set_from_length(3);
                edit5->set_to_length(3);
                subpath3->set_score(6);
                
                subpath0->add_next(1);
                subpath0->add_next(2);
                subpath1->add_next(3);
                subpath2->add_next(3);
                
                identify_start_subpaths(multipath_aln);
                
                REQUIRE(abs(genome.read_log_likelihood(multipath_aln, 1.0) - 13.0) < .001);
                
            }
            
            
            SECTION("Log likeihood is split up between partial alignments") {
                
                // constuct haplotypes
                
                list<NodeTraversal> haplotype_1;
                
                haplotype_1.push_back(NodeTraversal(n0));
                haplotype_1.push_back(NodeTraversal(n1));
                haplotype_1.push_back(NodeTraversal(n2));
                haplotype_1.push_back(NodeTraversal(n4));
                haplotype_1.push_back(NodeTraversal(n5));
                
                genome.add_haplotype(haplotype_1.begin(), haplotype_1.end());
                
                // index sites
                
                genome.build_indices();
                
                multipath_alignment_t multipath_aln;
                multipath_aln.set_sequence("CCTGCAAGAATAT");
                multipath_aln.set_mapping_quality(60);
                
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                subpath_t* subpath0 = multipath_aln.mutable_subpath(0);
                subpath_t* subpath2 = multipath_aln.mutable_subpath(1);
                subpath_t* subpath3 = multipath_aln.mutable_subpath(2);
                
                path_t* path0 = subpath0->mutable_path();
                path0->add_mapping();
                path0->add_mapping();
                path_mapping_t* mapping0 = path0->mutable_mapping(0);
                path_mapping_t* mapping4 = path0->mutable_mapping(1);
                
                position_t* position0 = mapping0->mutable_position();
                position0->set_node_id(n0->id());
                position0->set_is_reverse(false);
                edit_t* edit0 = mapping0->add_edit();
                edit0->set_from_length(3);
                edit0->set_to_length(3);
                position_t* position4 = mapping4->mutable_position();
                position4->set_node_id(n1->id());
                position4->set_is_reverse(false);
                edit_t* edit4 = mapping4->add_edit();
                edit4->set_from_length(3);
                edit4->set_to_length(3);
                subpath0->set_score(6);
                
                path_t* path2 = subpath2->mutable_path();
                path_mapping_t* mapping2 = path2->add_mapping();
                position_t* position2 = mapping2->mutable_position();
                position2->set_node_id(n3->id());
                position2->set_is_reverse(false);
                edit_t* edit2 = mapping2->add_edit();
                edit2->set_from_length(1);
                edit2->set_to_length(1);
                subpath2->set_score(-2);
                
                path_t* path3 = subpath3->mutable_path();
                path3->add_mapping();
                path3->add_mapping();
                path_mapping_t* mapping3 = path3->mutable_mapping(0);
                path_mapping_t* mapping5 = path3->mutable_mapping(1);
                
                position_t* position3 = mapping3->mutable_position();
                position3->set_node_id(n4->id());
                position3->set_is_reverse(false);
                edit_t* edit3 = mapping3->add_edit();
                edit3->set_from_length(3);
                edit3->set_to_length(3);
                position_t* position5 = mapping5->mutable_position();
                position5->set_node_id(n5->id());
                position5->set_is_reverse(false);
                edit_t* edit5 = mapping5->add_edit();
                edit5->set_from_length(3);
                edit5->set_to_length(3);
                subpath3->set_score(6);
                
                subpath0->add_next(1);
                subpath2->add_next(2);
                
                identify_start_subpaths(multipath_aln);
                
                //cerr << genome.read_log_likelihood(multipath_aln, 1.0) << endl;
                REQUIRE(abs(genome.read_log_likelihood(multipath_aln, 1.0) - (6.0 + log(2.0))) < .001);
                
            }
        }
    }
}
