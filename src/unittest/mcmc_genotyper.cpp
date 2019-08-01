/// \file multipath_alignment.cpp
///  
/// unit tests for mcmc genotyper construction and utility functions
///

#include <stdio.h>
#include <iostream>
#include "../multipath_mapper.hpp"
#include <gbwt/dynamic_gbwt.h>
#include "../build_index.hpp"
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
        // seed for random engine;
        const int seed =0;
        const int n_iterations = 10;

        /**
         * Moves an object from one vector and pushes it to the target vector
         * vect is *assumed* to have only one MultipathAlignment object
         */
        void accumulate_alns(vector<MultipathAlignment> vect, vector<MultipathAlignment> target){
            //moves the first object only 
            move(vect.begin(), vect.end(), back_inserter(target)); 
            // erase vect because it is not in an indeterminate state 
            vect.erase(vect.begin(), vect.end()); 
        }

        
        TEST_CASE("Test1") {
            
            SECTION("Tests run_genotyper() that takes a graph with one node and vector of one short read") {
                VG graph;
				
                Node* n1 = graph.create_node("GCA");
                cerr << "this is the node address"<< &n1 << endl;

                // name a path 
                path_handle_t path_handle = graph.create_path_handle("x");
                // append a visit to a node (via handle) to the given path 
                graph.append_step(path_handle, graph.get_handle(n1->id()));
                
                //Cactus does not currently support finding snarls in graph of single-node connected components
                SnarlManager snarl_manager; 
            
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

                MCMCGenotyper mcmc_genotyper = MCMCGenotyper(snarl_manager, graph, n_iterations, seed);   
                vector<MultipathAlignment> multipath_aln_vector = vector<MultipathAlignment>({multipath_aln}); 
                double log_base = gssw_dna_recover_log_base(1,4,.5,1e-12);
                unique_ptr<PhasedGenome> genome = mcmc_genotyper.run_genotype(multipath_aln_vector, log_base);

                // check requirements
                // here the optimal solution is that both haplotypes match to the read
                REQUIRE(genome->num_haplotypes() == 2);
                REQUIRE(multipath_aln.start_size() > 0);
                
                // haplotype 1 
                auto iter = genome->begin(0); //Haplotype ID
                REQUIRE(iter != genome->end(0));
                // check that iter is pointing to "valid" node and not empty
                cerr << "this is the iter node pointer"<< (*iter).node << endl;
                REQUIRE((*iter).node != nullptr); 
                REQUIRE((*iter).node == n1);
                REQUIRE((*iter).backward == false);
                REQUIRE((*iter) == NodeTraversal(n1)); 
                iter++;
                REQUIRE(iter == genome->end(0));

                // haplotype 2
                iter = genome->begin(1);
                // check that iter is pointing to "valid" node and not empty
                REQUIRE(iter != genome->end(1));
                REQUIRE((*iter).node != nullptr);
                REQUIRE((*iter).node == n1);
                REQUIRE((*iter).backward == false);
                REQUIRE((*iter) == NodeTraversal(n1)); 
                iter++;
                REQUIRE(iter == genome->end(1));


            }
 
        }
        TEST_CASE("Test2"){
            
            SECTION("Tests run_genotyper() with a simple graph containing one snarl and a vector with only one read") {   
                VG graph;
				
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("G");
                Node* n4 = graph.create_node("CTGA");

                // picked the suboptimal path so vg can have something to do
                path_handle_t path_handle = graph.create_path_handle("x");
                graph.append_step(path_handle, graph.get_handle(n1->id()));
                graph.append_step(path_handle, graph.get_handle(n3->id()));
                graph.append_step(path_handle, graph.get_handle(n4->id()));
                
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

            
                MCMCGenotyper mcmc_genotyper = MCMCGenotyper(snarl_manager, graph, n_iterations, seed);
                vector<MultipathAlignment> multipath_aln_vector = vector<MultipathAlignment>({multipath_aln}); 
                double log_base = gssw_dna_recover_log_base(1,4,.5,1e-12);
                unique_ptr<PhasedGenome> genome = mcmc_genotyper.run_genotype(multipath_aln_vector, log_base);

                // check requirements
                REQUIRE(genome->num_haplotypes() == 2);
                REQUIRE(multipath_aln.start_size() > 0);
                
                // haplotype 1 
                auto iter = genome->begin(0); //Haplotype ID
                // check that iter is pointing to "valid" node and not empty
                REQUIRE(iter != genome->end(0));
                REQUIRE((*iter).node != nullptr);
                REQUIRE((*iter) == NodeTraversal(n1)); 
                iter++;
                REQUIRE((*iter) == NodeTraversal(n2));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n4));        
                iter++;
                REQUIRE(iter == genome->end(0));

                // haplotype 2
                iter = genome->begin(1);
                // check that iter is pointing to "valid" node and not empty
                REQUIRE(iter != genome->end(1));
                REQUIRE((*iter).node != nullptr);

                REQUIRE((*iter) == NodeTraversal(n1)); 
                iter++;
                REQUIRE((*iter) == NodeTraversal(n2));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n4));        
                iter++;
                REQUIRE(iter == genome->end(1));
 
            } 
        }                  
        TEST_CASE("Test3"){
            
            SECTION( "Tests run_genotyper() with a graph containing two snarls and a vector with one read") {   
                VG graph;

                Node* n1 = graph.create_node("GCA"); //gets ID # in incremintal order starting at 1, used in mapping
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("G");
                Node* n4 = graph.create_node("CTGA"); //this node is part of both snarls
                Node* n5 = graph.create_node("A");
                Node* n6 = graph.create_node("G");
                Node* n7 = graph.create_node("CCC");

                path_handle_t path_handle = graph.create_path_handle("x");
                graph.append_step(path_handle, graph.get_handle(n1->id()));
                graph.append_step(path_handle, graph.get_handle(n3->id()));
                graph.append_step(path_handle, graph.get_handle(n4->id()));
                graph.append_step(path_handle, graph.get_handle(n5->id()));
                graph.append_step(path_handle, graph.get_handle(n7->id()));
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n4);
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
                
                multipath_aln.add_start(0); //integer given is the subpath number


               MCMCGenotyper mcmc_genotyper = MCMCGenotyper(snarl_manager, graph, n_iterations, seed);
                vector<MultipathAlignment> multipath_aln_vector = vector<MultipathAlignment>({multipath_aln}); 
                //instantiating an empty vector: 
                //vector<a> v;
                //v.push_back(multipath_align)

                double log_base = gssw_dna_recover_log_base(1,4,.5,1e-12);
                unique_ptr<PhasedGenome> genome = mcmc_genotyper.run_genotype(multipath_aln_vector, log_base);
				
                //check requirements
                REQUIRE(genome->num_haplotypes() == 2);
                REQUIRE(multipath_aln.start_size() > 0);
                // haplotype 1 
                auto iter = genome->begin(0); //Haplotype ID

                // check that iter is pointing to "valid" node and not empty
                REQUIRE(iter != genome->end(0));
                REQUIRE((*iter).node != nullptr);

                REQUIRE((*iter) == NodeTraversal(n1)); 
                iter++;
                REQUIRE((*iter) == NodeTraversal(n2));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n4));        
                iter++;
                REQUIRE((*iter) == NodeTraversal(n6));        
                iter++;
                REQUIRE((*iter) == NodeTraversal(n7));        
                iter++;
                REQUIRE(iter == genome->end(0));

                // haplotype 2
                iter = genome->begin(1);

                // check that iter is pointing to "valid" node and not empty
                REQUIRE(iter != genome->end(1));
                REQUIRE((*iter).node != nullptr);

                REQUIRE((*iter) == NodeTraversal(n1)); 
                iter++;
                REQUIRE((*iter) == NodeTraversal(n2));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n4));        
                iter++;
                REQUIRE((*iter) == NodeTraversal(n6));        
                iter++;
                REQUIRE((*iter) == NodeTraversal(n7));        
                iter++;
                REQUIRE(iter == genome->end(1));


            }
        }
        TEST_CASE("Test4"){
            
            SECTION( "Tests run_genotyper() using a graph containing nested snarls and a vector with one read") {   
                VG graph;
				
                
                Node* n1 = graph.create_node("GGG"); //gets ID in incremintal order 
                Node* n2 = graph.create_node("CCC");
                Node* n3 = graph.create_node("A");
                Node* n4 = graph.create_node("T");
                Node* n5 = graph.create_node("G");
                Node* n6 = graph.create_node("CTTG");
                Node* n7 = graph.create_node("TAC");
                Node* n8 = graph.create_node("C");

                path_handle_t path_handle = graph.create_path_handle("x");
                graph.append_step(path_handle, graph.get_handle(n1->id()));
                graph.append_step(path_handle, graph.get_handle(n2->id()));
                graph.append_step(path_handle, graph.get_handle(n4->id()));
                graph.append_step(path_handle, graph.get_handle(n5->id()));
                graph.append_step(path_handle, graph.get_handle(n6->id()));
                graph.append_step(path_handle, graph.get_handle(n7->id()));
                graph.append_step(path_handle, graph.get_handle(n8->id()));

                
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

                multipath_aln.add_start(0); //integer given is the subpath number
				
                // match and mismatch represent the point system used for penalty
                double log_base = gssw_dna_recover_log_base(1,4,.5,1e-12);

                MCMCGenotyper mcmc_genotyper = MCMCGenotyper(snarl_manager, graph, n_iterations, seed);
                vector<MultipathAlignment> multipath_aln_vector = vector<MultipathAlignment>({multipath_aln}); 
                unique_ptr<PhasedGenome> genome = mcmc_genotyper.run_genotype(multipath_aln_vector, log_base);

                //check requirements
                REQUIRE(genome->num_haplotypes() == 2);
                REQUIRE(multipath_aln.start_size() > 0);

                // haplotype 1 
                auto iter = genome->begin(0); //Haplotype ID
                REQUIRE((*iter) == NodeTraversal(n1)); 
                iter++;
                REQUIRE((*iter) == NodeTraversal(n2));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n3));        
                iter++;
                REQUIRE((*iter) == NodeTraversal(n5));        
                iter++;
                REQUIRE((*iter) == NodeTraversal(n6));        
                iter++;
                REQUIRE(iter == genome->end(0));

                // haplotype 2
                iter = genome->begin(1);
                REQUIRE((*iter) == NodeTraversal(n1)); 
                iter++;
                REQUIRE((*iter) == NodeTraversal(n2));
                iter++;
                REQUIRE((*iter) == NodeTraversal(n3));        
                iter++;
                REQUIRE((*iter) == NodeTraversal(n5));        
                iter++;
                REQUIRE((*iter) == NodeTraversal(n6));        
                iter++;
                REQUIRE(iter == genome->end(1));
               
            }    

        }
        TEST_CASE("Test5"){
            SECTION("Tests run_genotyper using graph with two snarls and several reads"){
                VG graph;

                Node* n1 = graph.create_node("GCA"); //gets ID # in incremintal order starting at 1, used in mapping
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("G");
                Node* n4 = graph.create_node("CTGA"); //this node is part of both snarls
                Node* n5 = graph.create_node("A");
                Node* n6 = graph.create_node("G");
                Node* n7 = graph.create_node("CCC");

                path_handle_t path_handle = graph.create_path_handle("x");
                graph.append_step(path_handle, graph.get_handle(n1->id()));
                graph.append_step(path_handle, graph.get_handle(n2->id()));
                graph.append_step(path_handle, graph.get_handle(n3->id()));
                graph.append_step(path_handle, graph.get_handle(n4->id()));
                graph.append_step(path_handle, graph.get_handle(n5->id()));
                graph.append_step(path_handle, graph.get_handle(n6->id()));
                graph.append_step(path_handle, graph.get_handle(n7->id()));
                
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n4);
                graph.create_edge(n3, n4);
                graph.create_edge(n4, n5);
                graph.create_edge(n4, n6);
                graph.create_edge(n5, n7);
                graph.create_edge(n6, n7);
				
				CactusSnarlFinder bubble_finder(graph);
                SnarlManager snarl_manager = bubble_finder.find_snarls();

                // Configure GCSA temp directory to the system temp directory
                gcsa::TempFile::setDirectory(temp_file::get_dir());
                // And make it quiet
                gcsa::Verbosity::set(gcsa::Verbosity::SILENT);
                
                // Make pointers to fill in
                gcsa::GCSA* gcsaidx = nullptr;
                gcsa::LCPArray* lcpidx = nullptr;
                
                // Build the GCSA index
                build_gcsa_lcp(graph, gcsaidx, lcpidx, 16, 3); 
            
                
                // Build the xg index
                //defining an XG and a variable called xg_index and calling the constructor
                 XG xg_index(graph.graph); //VG uses a Graph as internal structure 
                
                // Make a multipath mapper to map against the graph.
                MultipathMapper multipath_mapper(&xg_index, gcsaidx, lcpidx); 

                string read1 = string("GCATCTGAGCCC"); 
                string read2 = string("GCATCTGAGCCC");
                string read3 = string("GCAGCTGAACCC");
                string read4 = string("GCAGCTGAACCC");
                
                Alignment aln1, aln2, aln3, aln4;
                aln1.set_sequence(read1);
                aln2.set_sequence(read2);
                aln3.set_sequence(read3);
                aln4.set_sequence(read4);

                vector<MultipathAlignment> multipath_alns_out1, multipath_alns_out2, multipath_alns_out3, multipath_alns_out4 ;
                MCMCGenotyper mcmc_genotyper = MCMCGenotyper(snarl_manager, graph, n_iterations, seed);
                vector<MultipathAlignment> multipath_aln_vector = vector<MultipathAlignment>(); 

                // map read in alignment to graph and make multipath alignments  
                multipath_mapper.multipath_map(aln1, multipath_alns_out1, 1);
                multipath_mapper.multipath_map(aln2, multipath_alns_out2, 1);
                multipath_mapper.multipath_map(aln3, multipath_alns_out3, 1);
                multipath_mapper.multipath_map(aln4, multipath_alns_out4, 1);
                
                
                //accumulate MultipathAlignment objects 
                accumulate_alns(multipath_alns_out1, multipath_aln_vector);
                accumulate_alns(multipath_alns_out2, multipath_aln_vector);
                accumulate_alns(multipath_alns_out3, multipath_aln_vector);
                accumulate_alns(multipath_alns_out4, multipath_aln_vector);

                double log_base = gssw_dna_recover_log_base(1,4,.5,1e-12);
                //pass vector with accumulated MultipathAlignment objects to run_genotype()
                unique_ptr<PhasedGenome> genome = mcmc_genotyper.run_genotype(multipath_aln_vector, log_base); 
                
                // create a set of 2 possible solutions
                vector<NodeTraversal> soln1, soln2;
                soln1 = {NodeTraversal(n1), NodeTraversal(n2), NodeTraversal(n4), NodeTraversal(n5), NodeTraversal(n7)};
                soln2 = {NodeTraversal(n1), NodeTraversal(n3), NodeTraversal(n4), NodeTraversal(n6), NodeTraversal(n7)};

                set<vector<NodeTraversal>> solns_set;
                solns_set.insert(soln1);
                solns_set.insert(soln2);

                // check requirements 
                REQUIRE(genome->num_haplotypes() == 2);

                // move the genome haplotype into a vector
                vector<NodeTraversal> haplotype1, haplotype2;
                copy(genome->begin(0), genome->end(0), back_inserter(haplotype1));
                copy(genome->begin(1), genome->end(1), back_inserter(haplotype2));
                
                // check haplotypes are indeed in optimal solution set   
                REQUIRE(solns_set.count(haplotype1));
                REQUIRE(solns_set.count(haplotype2));

                // Clean up the GCSA/LCP index
                delete gcsaidx;
                delete lcpidx;
            }
        }
    }

}


