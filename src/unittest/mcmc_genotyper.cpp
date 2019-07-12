/// \file multipath_alignment.cpp
///  
/// unit tests for mcmc genotyper construction and utility functions
///

#include <stdio.h>
#include <iostream>

#include <gbwt/dynamic_gbwt.h>

#include "../haplotypes.hpp"
#include "../json2pb.h"
#include <vg/vg.pb.h>
#include "../vg.hpp"
#include "../multipath_alignment.hpp"
#include "../utility.hpp"
#include "../mcmc_genotyper.hpp"
#include "../snarls.hpp"
#include "../genotypekit.hpp"

#include "catch.hpp"

namespace vg {
    namespace unittest {
        
        TEST_CASE( "Test1") {
            
            SECTION( "An empty Graph and empty multipath alignment") {
                VG graph;
				
				CactusSnarlFinder bubble_finder(graph);
                SnarlManager snarl_manager = bubble_finder.find_snarls();
                
                string read = string(""); //empty read 
                MultipathAlignment multipath_aln;
                multipath_aln.set_sequence(read);

                // no subpaths or mappings on an empty graph and empty alignment                
				
                // check requirements 
                // create MCMC_Genotyper object 
                // create and instantiate vector of mp alignments
                // call run_genotyper and use genome object to test if it what I expect 
                /** (genome.num_haplotREQUIREypes == 0);*/

            

            }
 
        }
        TEST_CASE( "Test2") {
            
            SECTION( "A Graph with one node") {
                VG graph;
				
                Node* n1 = graph.create_node("GCA");
                
                graph.create_edge(n1, n1); //if the edge cannot be created returns null

				CactusSnarlFinder bubble_finder(graph);
                SnarlManager snarl_manager = bubble_finder.find_snarls();

                
                string read = string("GCA");
                MultipathAlignment multipath_aln;
                multipath_aln.set_sequence(read);
                
                 // add subpaths
                Subpath* subpath0 = multipath_aln.add_subpath();
                 
                // designate mappings
                Mapping* mapping0 = subpath0->mutable_path()->add_mapping();
                mapping0->mutable_position()->set_node_id(1);
                Edit* edit0 = mapping0->add_edit();
                edit0->set_from_length(3); //a match 
                edit0->set_to_length(3);

                subpath0->set_score(1);

                multipath_aln.add_start(0);

                // check requirements
                // TODO: check that my output is what I expect it to be  
                
                
            

            }
 
        }
        TEST_CASE("Test3"){
                  SECTION( "A Graph with one snarl") {   
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
                
                string read = string("GCATCTGA");
                MultipathAlignment multipath_aln;
                multipath_aln.set_sequence(read);
                
                // add subpaths with same topology as graph
                Subpath* subpath0 = multipath_aln.add_subpath();
                Subpath* subpath1 = multipath_aln.add_subpath();
                Subpath* subpath2 = multipath_aln.add_subpath();
                Subpath* subpath3 = multipath_aln.add_subpath();
                
                // set edges between subpaths
                subpath0->add_next(1);
                subpath0->add_next(2);
                subpath1->add_next(3);
                subpath2->add_next(3);
                
                // set scores 
                subpath0->set_score(1);
                subpath1->set_score(1);
                subpath2->set_score(-4);
                subpath3->set_score(1);

                // designate mappings
                Mapping* mapping0 = subpath0->mutable_path()->add_mapping();
                mapping0->mutable_position()->set_node_id(1);
                Edit* edit0 = mapping0->add_edit();
                edit0->set_from_length(3); //a match 
                edit0->set_to_length(3);
                
                Mapping* mapping1 = subpath1->mutable_path()->add_mapping();
                mapping1->mutable_position()->set_node_id(2);
                Edit* edit1 = mapping1->add_edit();
                edit1->set_from_length(1); //a match 
                edit1->set_to_length(1);
                
                Mapping* mapping2 = subpath2->mutable_path()->add_mapping();
                mapping2->mutable_position()->set_node_id(3);
                Edit* edit2 = mapping2->add_edit();
                edit2->set_from_length(1); //a snip
                edit2->set_to_length(1);
                edit2->set_sequence("T"); 

                Mapping* mapping3 = subpath3->mutable_path()->add_mapping();
                mapping3->mutable_position()->set_node_id(4);
                Edit* edit3 = mapping3->add_edit();
                edit3->set_from_length(4); //a match 
                edit3->set_to_length(4);
                
                multipath_aln.add_start(0);

            
                MCMCGenotyper mcmc_genotyper = MCMCGenotyper(snarl_manager);    //object
                vector<MultipathAlignment> multipath_aln_vector = vector<MultipathAlignment>({multipath_aln}); 
                //instantiating an empty vector: 
                //vector<a> v;
                //v.push_back(multipath_align)
                PhasedGenome genome = mcmc_genotyper.run_genotype(multipath_aln_vector);

                // check requirements
                REQUIRE(genome.num_haplotypes() == 2);
                // haplotype 1 
                auto iter = genome.begin(0); //Haplotype ID
                REQUIRE((*iter) == NodeTraversal(n1)); 
                iter++;
                REQUIRE((*iter) == NodeTraversal(n2));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n4));        
                iter++;
                REQUIRE(iter == genome.end(0));

                // haplotype 2
                iter = genome.begin(1);
                REQUIRE((*iter) == NodeTraversal(n1)); 
                iter++;
                REQUIRE((*iter) == NodeTraversal(n2));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n4));        
                iter++;
                REQUIRE(iter == genome.end(1));
                
                
                
            }
        
        }
        TEST_CASE("Test4"){
                SECTION( "A Graph with two snarls") {   
                VG graph;
				
                
                Node* n1 = graph.create_node("GCA"); //gets ID # in incremintal order starting at 1, used in mapping
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("G");
                Node* n4 = graph.create_node("CTGA"); //this node is part of both snarls
                Node* n5 = graph.create_node("A");
                Node* n6 = graph.create_node("G");
                Node* n7 = graph.create_node("CCC");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n4);
                graph.create_edge(n3, n4);
                graph.create_edge(n3, n4);
                graph.create_edge(n4, n5);
                graph.create_edge(n4, n6);
                graph.create_edge(n5, n7);
                graph.create_edge(n6, n7);
				
				CactusSnarlFinder bubble_finder(graph);
                SnarlManager snarl_manager = bubble_finder.find_snarls();
                
                string read = string("GCATCTGAGCCC");
                MultipathAlignment multipath_aln;
                multipath_aln.set_sequence(read);
                
                // add subpaths with same topology as graph
                Subpath* subpath0 = multipath_aln.add_subpath();
                Subpath* subpath1 = multipath_aln.add_subpath();
                Subpath* subpath2 = multipath_aln.add_subpath();
                Subpath* subpath3 = multipath_aln.add_subpath();
                Subpath* subpath4 = multipath_aln.add_subpath();
                Subpath* subpath5 = multipath_aln.add_subpath();
                Subpath* subpath6 = multipath_aln.add_subpath();
                
                // set edges between subpaths
                subpath0->add_next(1);
                subpath0->add_next(2);
                subpath1->add_next(3);
                subpath2->add_next(3);
                subpath3->add_next(4);
                subpath3->add_next(5);
                subpath4->add_next(6);
                subpath5->add_next(6);
                
                // set scores
                subpath0->set_score(1);
                subpath1->set_score(1);
                subpath2->set_score(-4);
                subpath3->set_score(1);
                subpath4->set_score(-4);
                subpath5->set_score(1);
                subpath6->set_score(1); 
                
                // designate mappings
                Mapping* mapping0 = subpath0->mutable_path()->add_mapping();
                mapping0->mutable_position()->set_node_id(1);
                Edit* edit0 = mapping0->add_edit();
                edit0->set_from_length(3); //a match 
                edit0->set_to_length(3);
                
                Mapping* mapping1 = subpath1->mutable_path()->add_mapping();
                mapping1->mutable_position()->set_node_id(2);
                Edit* edit1 = mapping1->add_edit();
                edit1->set_from_length(1); //a match 
                edit1->set_to_length(1);
                
                Mapping* mapping2 = subpath2->mutable_path()->add_mapping();
                mapping2->mutable_position()->set_node_id(3);
                Edit* edit2 = mapping2->add_edit();
                edit2->set_from_length(1); //a mismatch , snp
                edit2->set_to_length(1);
                edit2->set_sequence("T");

                Mapping* mapping3 = subpath3->mutable_path()->add_mapping();
                mapping3->mutable_position()->set_node_id(4);
                Edit* edit3 = mapping3->add_edit();
                edit3->set_from_length(3); //a match 
                edit3->set_to_length(3);

                Mapping* mapping4 = subpath4->mutable_path()->add_mapping();
                mapping4->mutable_position()->set_node_id(5);
                Edit* edit4 = mapping4->add_edit();
                edit4->set_from_length(1); //a mismatch, snp
                edit4->set_to_length(1);
                edit4->set_sequence("G");

                Mapping* mapping5 = subpath5->mutable_path()->add_mapping();
                mapping5->mutable_position()->set_node_id(6);
                Edit* edit5 = mapping5->add_edit();
                edit5->set_from_length(1); //a match 
                edit5->set_to_length(1);

                Mapping* mapping6 = subpath6->mutable_path()->add_mapping();
                mapping6->mutable_position()->set_node_id(7);
                Edit* edit6 = mapping6->add_edit();
                edit6->set_from_length(3); //a match 
                edit6->set_to_length(3);
                
                multipath_aln.add_start(0); //integer given is the subpath #
				
                //check requirements
            }
        }
        TEST_CASE("Test5"){
            SECTION( "A Graph with nested snarls") {   
                VG graph;
				
                
                Node* n1 = graph.create_node("GGG"); //gets ID in incremintal order 
                Node* n2 = graph.create_node("CCC");
                Node* n3 = graph.create_node("A");
                Node* n4 = graph.create_node("T");
                Node* n5 = graph.create_node("G");
                Node* n6 = graph.create_node("CTTG");
                Node* n7 = graph.create_node("TAC");
                Node* n8 = graph.create_node("C");

                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n7);
                graph.create_edge(n2, n3);
                graph.create_edge(n2, n4);
                graph.create_edge(n3, n5);
                graph.create_edge(n4, n5);
                graph.create_edge(n5, n6);
                graph.create_edge(n7, n8);
                graph.create_edge(n8, n6);
				
				CactusSnarlFinder bubble_finder(graph);
                SnarlManager snarl_manager = bubble_finder.find_snarls();
                PhasedGenome genome = PhasedGenome(snarl_manager);
                
                string read = string("GGGCCCAGCTGG");
                MultipathAlignment multipath_aln;
                multipath_aln.set_sequence(read);
                
                // add subpaths with same topology as graph
                Subpath* subpath0 = multipath_aln.add_subpath();
                Subpath* subpath1 = multipath_aln.add_subpath();
                Subpath* subpath2 = multipath_aln.add_subpath();
                Subpath* subpath3 = multipath_aln.add_subpath();
                Subpath* subpath4 = multipath_aln.add_subpath();
                Subpath* subpath5 = multipath_aln.add_subpath();
                Subpath* subpath6 = multipath_aln.add_subpath();
                Subpath* subpath7 = multipath_aln.add_subpath();
                
                // set edges between subpaths
                subpath0->add_next(1);
                subpath0->add_next(6);
                subpath1->add_next(2);
                subpath1->add_next(3);
                subpath2->add_next(4);
                subpath3->add_next(4);
                subpath4->add_next(5);
                subpath6->add_next(7);
                subpath7->add_next(5);

                // set scores
                subpath0->set_score(1);
                subpath1->set_score(1);
                subpath2->set_score(1);
                subpath3->set_score(-4);
                subpath4->set_score(1);
                subpath5->set_score(1);
                subpath6->set_score(-3);
                subpath7->set_score(-10);
                
                // designate mappings
                Mapping* mapping0 = subpath0->mutable_path()->add_mapping();
                mapping0->mutable_position()->set_node_id(1);
                Edit* edit0 = mapping0->add_edit();
                edit0->set_from_length(3); //a match 
                edit0->set_to_length(3);
                
                Mapping* mapping1 = subpath0->mutable_path()->add_mapping();
                mapping1->mutable_position()->set_node_id(2);
                Edit* edit1 = mapping1->add_edit();
                edit1->set_from_length(3); //a match 
                edit1->set_to_length(3);
                
                Mapping* mapping2 = subpath1->mutable_path()->add_mapping();
                mapping2->mutable_position()->set_node_id(3);
                Edit* edit2 = mapping2->add_edit();
                edit2->set_from_length(1); //a match 
                edit2->set_to_length(1);
                
                Mapping* mapping3 = subpath3->mutable_path()->add_mapping();
                mapping3->mutable_position()->set_node_id(4);
                Edit* edit3 = mapping3->add_edit();
                edit3->set_from_length(1); //a mismatch, snp
                edit3->set_to_length(1);
                edit3->set_sequence("A");

                Mapping* mapping4 = subpath4->mutable_path()->add_mapping();
                mapping4->mutable_position()->set_node_id(5);
                Edit* edit4 = mapping4->add_edit();
                edit4->set_from_length(1); //a match 
                edit4->set_to_length(1);

                Mapping* mapping5 = subpath5->mutable_path()->add_mapping();
                mapping5->mutable_position()->set_node_id(6);
                Edit* edit5 = mapping5->add_edit();
                edit5->set_from_length(3); //a match 
                edit5->set_to_length(3);

                Mapping* mapping6 = subpath6->mutable_path()->add_mapping();
                mapping6->mutable_position()->set_node_id(7);
                Edit* edit6 = mapping6->add_edit();
                edit6->set_from_length(2); //a mismatch 
                edit6->set_to_length(2); // will offset 2 positions in read and graph 
                edit6->set_sequence("TA");
                Edit* edit6b = mapping6->add_edit();
                edit6b->set_from_length(1); // a match
                edit6b->set_to_length(1);
                

                Mapping* mapping7 = subpath7->mutable_path()->add_mapping();
                mapping7->mutable_position()->set_node_id(8);
                Edit* edit7 = mapping7->add_edit();
                edit7->set_from_length(1); //a mismatch 
                edit7->set_to_length(1);
                edit7->set_sequence("T");
                Edit* edit7b = mapping7->add_edit(); // a deletion , gap open 
                edit7b->set_from_length(1); // denotes the length of the gap observed in alt seq
                edit7b->set_to_length(0); // length in alt sequence is zero because it is a gap
                edit7b->set_sequence("C");

                multipath_aln.add_start(0); //integer given is the subpath #
				
                // check requirements
               
            }    

        }
    }
}


