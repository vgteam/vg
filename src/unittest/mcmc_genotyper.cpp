/// \file multipath_alignment.cpp
///  
/// unit tests for mcmc genotyper construction and utility functions
///
#include <xg.hpp>
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
        
        const int seed = 0;
        const int n_iterations = 100;

        
        TEST_CASE("Test1") {
            //Returns optimal phased genome on a 1-node graph with 1 short read 
            SECTION("Test1: Requires haplotype pair to match truth set") {
                VG graph;
				
                Node* n1 = graph.create_node("GCA");
                //cerr << "this is the node address"<< &n1 << endl;

                // name a path 
                path_handle_t path_handle = graph.create_path_handle("x");
                // append a visit to a node (via Tests run_genotyper() that takes handle) to the given path 
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
                
            
                // create a set of 2 possible solutions
                vector<NodeTraversal> soln1;
                soln1 = {NodeTraversal(n1)};
                

                set<vector<NodeTraversal>> solns_set;
                solns_set.insert(soln1);

                 // move the genome haplotype into a vector
                vector<NodeTraversal> haplotype1, haplotype2;
                copy(genome->begin(0), genome->end(0), back_inserter(haplotype1));
                copy(genome->begin(1), genome->end(1), back_inserter(haplotype2));
                
                // check haplotypes are indeed in optimal solution set
                bool pass = false;
                if(solns_set.count(haplotype1)&& solns_set.count(haplotype2)){
                    pass = true;
                }
                // requires both haplotypes to match the one read
                REQUIRE(pass);   
            }
 
        }
        TEST_CASE("Test2"){
            //Returns optimal phased genome on a 4-node graph, snarl with 1 short read
            SECTION("Test2: Requires haplotype pair to match truth set") {   
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
                //cerr << "**************" <<endl;
                //genome->print_phased_genome();
                // check requirements
                REQUIRE(genome->num_haplotypes() == 2);
                REQUIRE(multipath_aln.start_size() > 0);
                
            
                // create a set of 2 possible solutions
                vector<NodeTraversal> soln1;
                soln1 = {NodeTraversal(n1), NodeTraversal(n2), NodeTraversal(n4)};

                set<vector<NodeTraversal>> solns_set;
                solns_set.insert(soln1);

                 // move the genome haplotype into a vector
                vector<NodeTraversal> haplotype1, haplotype2;
                copy(genome->begin(0), genome->end(0), back_inserter(haplotype1));
                copy(genome->begin(1), genome->end(1), back_inserter(haplotype2));
                
                // check haplotypes are indeed in optimal solution set  
                bool pass = false;
                if(solns_set.count(haplotype1) || solns_set.count(haplotype2)){
                    pass = true;
                }
                // requires at least one of the haplotypes to match truth set 
                // only one match is required because there is only one read 
                REQUIRE(pass);
                
            } 
        }                  
        TEST_CASE("Test3"){
            //Returns optimal phased genome on a 7-node graph two connected snarls with 1 short read
            SECTION("Test3: Requires haplotype pair to match truth set") {   
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
                xg::XG xg_index; 
                xg_index.from_path_handle_graph(graph);              

                // Make a multipath mapper to map against the graph.
                MultipathMapper multipath_mapper(&xg_index, gcsaidx, lcpidx); 
                
                vector<string> reads = {"GCATCTGAGCCC"};
                vector<Alignment> alns = {reads.size(), Alignment()};

                // set alignment sequence
                for(int i = 0; i< reads.size(); i++){
                    alns[i].set_sequence(reads[i]);
                }

               MCMCGenotyper mcmc_genotyper = MCMCGenotyper(snarl_manager, graph, n_iterations, seed);
               vector<MultipathAlignment> multipath_aln_vector = vector<MultipathAlignment>(); 

               vector<vector<MultipathAlignment>> vect = {reads.size(),vector<MultipathAlignment>() };
                    
                    
                // map read in alignment to graph and make multipath alignments 
                for(int i = 0; i< reads.size(); i++){
                    multipath_mapper.multipath_map(alns[i], vect[i], 1);
                }

                // accumulate the mapped reads in one vector
                for(int i = 0; i< reads.size(); i++){
                    move(vect[i].begin(), vect[i].end(), back_inserter(multipath_aln_vector)); 
                }
                    
                double log_base = gssw_dna_recover_log_base(1,4,.5,1e-12);
                unique_ptr<PhasedGenome> genome = mcmc_genotyper.run_genotype(multipath_aln_vector, log_base);

                 // create a set of 2 possible solutions
                vector<NodeTraversal> soln1;
                soln1 = {NodeTraversal(n1), NodeTraversal(n2), NodeTraversal(n4), NodeTraversal(n6), NodeTraversal(n7)};

                set<vector<NodeTraversal>> solns_set;
                solns_set.insert(soln1);

                // check requirements 
                REQUIRE(genome->num_haplotypes() == 2);

                // move the genome haplotype into a vector
                vector<NodeTraversal> haplotype1, haplotype2;
                copy(genome->begin(0), genome->end(0), back_inserter(haplotype1));
                copy(genome->begin(1), genome->end(1), back_inserter(haplotype2));
                
                bool pass = false;
                if(solns_set.count(haplotype1) || solns_set.count(haplotype2)){
                    pass = true;
                }
                // requires at least one of the haplotypes to match truth set 
                // only one match is required because there is only one read
                REQUIRE(pass);
                
                

            }
        }
        TEST_CASE("Test4"){
            //Returns optimal phased genome on a 8-node graph 2 nested snarls with 1 short read
            SECTION("Test4: Requires haplotype pair to match truth set") {   
                VG graph;
				
                
                Node* n1 = graph.create_node("GGG"); //gets ID in incremintal order 
                Node* n2 = graph.create_node("CCC");
                Node* n3 = graph.create_node("A");
                Node* n4 = graph.create_node("T");
                Node* n5 = graph.create_node("G");
                Node* n6 = graph.create_node("CTGG");
                Node* n7 = graph.create_node("TAC");
                Node* n8 = graph.create_node("C");

                path_handle_t path_handle = graph.create_path_handle("x");
                graph.append_step(path_handle, graph.get_handle(n1->id()));
                graph.append_step(path_handle, graph.get_handle(n7->id()));
                graph.append_step(path_handle, graph.get_handle(n8->id()));
                graph.append_step(path_handle, graph.get_handle(n6->id()));

                
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
                
                Mapping* mapping1 = subpath1->mutable_path()->add_mapping();
                mapping1->mutable_position()->set_node_id(2);
                Edit* edit1 = mapping1->add_edit();
                edit1->set_from_length(3); //a match 
                edit1->set_to_length(3);
                
                Mapping* mapping2 = subpath2->mutable_path()->add_mapping();
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
                edit6->set_sequence("CC");
                Edit* edit6b = mapping6->add_edit();
                edit6b->set_from_length(1); // a match
                edit6b->set_to_length(1);
                

                Mapping* mapping7 = subpath7->mutable_path()->add_mapping();
                mapping7->mutable_position()->set_node_id(8);
                Edit* edit7 = mapping7->add_edit();
                edit7->set_from_length(1); //a mismatch 
                edit7->set_to_length(1);
                edit7->set_sequence("A");
                Edit* edit7b = mapping7->add_edit(); // a deletion , gap open 
                edit7b->set_from_length(1); // denotes the length of the gap observed in alt seq
                edit7b->set_to_length(0); // length in alt sequence is zero because it is a gap
                edit7b->set_sequence("G");

                multipath_aln.add_start(0); //integer given is the subpath number
				
                // match and mismatch represent the point system used for penalty
                double log_base = gssw_dna_recover_log_base(1,4,.5,1e-12);

                MCMCGenotyper mcmc_genotyper = MCMCGenotyper(snarl_manager, graph, n_iterations, seed);
                vector<MultipathAlignment> multipath_aln_vector = vector<MultipathAlignment>({multipath_aln}); 
                unique_ptr<PhasedGenome> genome = mcmc_genotyper.run_genotype(multipath_aln_vector, log_base);

                 // create a set of possible solutions
                vector<NodeTraversal> soln1;
                soln1 = {NodeTraversal(n1), NodeTraversal(n2), NodeTraversal(n3), NodeTraversal(n5), NodeTraversal(n6)};

                set<vector<NodeTraversal>> solns_set;
                solns_set.insert(soln1); 

                // check requirements 
                REQUIRE(genome->num_haplotypes() == 2);

                // move the genome haplotype into a vector
                vector<NodeTraversal> haplotype1, haplotype2;
                copy(genome->begin(0), genome->end(0), back_inserter(haplotype1));
                copy(genome->begin(1), genome->end(1), back_inserter(haplotype2));
                
                vector<vector<NodeTraversal>> haplos = {haplotype1, haplotype2};
                
                bool pass = false;
                if(solns_set.count(haplotype1) || solns_set.count(haplotype2)){
                    pass = true;
                }
                // requires at least one of the haplotypes to match truth set 
                // only one match is required because there is only one read
                REQUIRE(pass);
            
            }    

        }
        TEST_CASE("Test5"){
            //Returns optimal phased genome on a 7-node graph containing 2 connected snarls, and 4 mapped reads
            SECTION("Test5: Requires haplotype pair to match truth set"){
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
                xg::XG xg_index; 
                xg_index.from_path_handle_graph(graph);              

                // Make a multipath mapper to map against the graph.
                MultipathMapper multipath_mapper(&xg_index, gcsaidx, lcpidx); 
                
                //vector<string> reads = {"GCATCTGAGCCC","GCATCTGAGCCC","GCAGCTGAACCC","GCAGCTGAACCC"}; //TGGA
                vector<string> reads = {"GCATCTGAGCCC","GCATCTGAGCCC","GCATCTGAACCC", "GCATCTGAACCC"};//TGTA
                vector<Alignment> alns = {reads.size(), Alignment()};

                // set alignment sequence
                for(int i = 0; i< reads.size(); i++){
                    alns[i].set_sequence(reads[i]);
                }
                
                
                MCMCGenotyper mcmc_genotyper = MCMCGenotyper(snarl_manager, graph, n_iterations, seed);
                vector<MultipathAlignment> multipath_aln_vector = vector<MultipathAlignment>(); 

                vector<vector<MultipathAlignment>> vect = {reads.size(),vector<MultipathAlignment>() };
                    
                    
                // map read in alignment to graph and make multipath alignments 
                for(int i = 0; i< reads.size(); i++){
                    multipath_mapper.multipath_map(alns[i], vect[i], 1);
                }
                    
                // accumulate the mapped reads in one vector
                for(int i = 0; i< reads.size(); i++){
                    move(vect[i].begin(), vect[i].end(), back_inserter(multipath_aln_vector)); 
                }

                double log_base = gssw_dna_recover_log_base(1,4,.5,1e-12);
                
                //pass vector with accumulated MultipathAlignment objects to run_genotype()
                unique_ptr<PhasedGenome> genome = mcmc_genotyper.run_genotype(multipath_aln_vector, log_base); 
                
                // // create a set of 2 possible solutions
                // // ex 2: TG GA
                // vector<NodeTraversal> soln1, soln2;
                // soln1 = {NodeTraversal(n1), NodeTraversal(n2), NodeTraversal(n4), NodeTraversal(n6), NodeTraversal(n7)};
                // soln2 = {NodeTraversal(n1), NodeTraversal(n3), NodeTraversal(n4), NodeTraversal(n5), NodeTraversal(n7)};

                // set<vector<NodeTraversal>> solns_set;
                // solns_set.insert(soln1);
                // solns_set.insert(soln2);

                //ex 1: TG, TA
                vector<NodeTraversal> soln1, soln2;
                soln1 = {NodeTraversal(n1), NodeTraversal(n2), NodeTraversal(n4), NodeTraversal(n6), NodeTraversal(n7)};
                soln2 = {NodeTraversal(n1), NodeTraversal(n2), NodeTraversal(n4), NodeTraversal(n5), NodeTraversal(n7)};

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
                bool pass = false; 
                if(solns_set.count(haplotype1) && solns_set.count(haplotype2)){
                    pass = true;
                }
                REQUIRE(pass);
                

                // Clean up the GCSA/LCP index
                delete gcsaidx;
                delete lcpidx;
            }
        }
        TEST_CASE("Test6"){
            //Returns optimal phased genome on a 7-node graph, two connected snarls with 8 short read
            SECTION("Test6: Requires haplotype pair to match truth set"){

                double count_correct = 0.0;
                double count_incorrect = 0.0;
                double count_half_correct = 0.0;
                double count_minus = 0.0;
                vector<double> results = vector<double>();
                
                int num_iterations = 1000;
                int max = 30;
                
                for(int seed_i = 0; seed_i < max; seed_i++){
                    
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
                    xg::XG xg_index; //VG uses a Graph as internal structure 
                    xg_index.from_path_handle_graph(graph);               //xg::XG xg_index();
                    
                    // Make a multipath mapper to map against the graph.
                    MultipathMapper multipath_mapper(&xg_index, gcsaidx, lcpidx); 

                    
                    vector<string> reads = {"GCATCTGAGCCC", "GCATCTGAGCCC", "GCAGCTGAACCC", "GCAGCTGAACCC","GCAGCTGAACCC", "GCAGCTGAACCC", "GCAGCTGAGCCC", "GCAGCTGAGCCC"};
                    vector<Alignment> alns = {reads.size(), Alignment()};
                    
                    // set alignment sequence
                    for(int i = 0; i< reads.size(); i++){
                        alns[i].set_sequence(reads[i]);
                    }
                    
                    MCMCGenotyper mcmc_genotyper = MCMCGenotyper(snarl_manager, graph, num_iterations, seed_i);
                    vector<MultipathAlignment> multipath_aln_vector = vector<MultipathAlignment>(); 

                    vector<vector<MultipathAlignment>> vect = {reads.size(),vector<MultipathAlignment>() };
                    
                    
                    // map read in alignment to graph and make multipath alignments 
                    for(int i = 0; i< reads.size(); i++){
                        multipath_mapper.multipath_map(alns[i], vect[i], 1);
                    }
                    
                    // accumulate the mapped reads in one vector
                    for(int i = 0; i< reads.size(); i++){
                        move(vect[i].begin(), vect[i].end(), back_inserter(multipath_aln_vector)); 
                    }
                    

                    double log_base = gssw_dna_recover_log_base(1,4,.5,1e-12);
                    //pass vector with accumulated MultipathAlignment objects to run_genotype()
                    unique_ptr<PhasedGenome> genome = mcmc_genotyper.run_genotype(multipath_aln_vector, log_base); 
                    

                    // create a set of 2 possible solutions
                    vector<NodeTraversal> soln1, soln2;
                    soln1 = {NodeTraversal(n1), NodeTraversal(n2), NodeTraversal(n4), NodeTraversal(n6), NodeTraversal(n7)};
                    soln2 = {NodeTraversal(n1), NodeTraversal(n3), NodeTraversal(n4), NodeTraversal(n5), NodeTraversal(n7)};

                    set<vector<NodeTraversal>> solns_set;
                    solns_set.insert(soln1);
                    solns_set.insert(soln2);

                    REQUIRE(genome->num_haplotypes() == 2);

                    // move the genome haplotype into a vector
                    vector<NodeTraversal> haplotype1, haplotype2;
                    copy(genome->begin(0), genome->end(0), back_inserter(haplotype1));
                    copy(genome->begin(1), genome->end(1), back_inserter(haplotype2));
                    
                    
                    if(genome->num_haplotypes() == 2){
                        if(solns_set.count(haplotype1) && solns_set.count(haplotype2)){
                            count_correct++;
                        }
                        else if(solns_set.count(haplotype1) || solns_set.count(haplotype2)){
                            count_half_correct++;
                        }else{
                            count_incorrect++;

                        }
                    }
                    
                    // Clean up the GCSA/LCP index
                    delete gcsaidx;
                    delete lcpidx;
                }

                // cerr <<"****************************DONE TESTING****************************" << endl;
                // cerr << count_incorrect << " tests "<< "out of " << max <<" did not match any haplotypes from haplotype pair " << endl;
                // cerr << count_correct << " tests" << " out of " << max <<" matched both haplotypes from the haplotype pair " <<endl;
                // cerr << count_half_correct << " tests " << "out of " << max << " matched one haplotype from haplotype pair  " <<endl;
                // int percent_half_correct = (count_half_correct/max)*100;
                // int percent_correct = (count_correct/max)*100;  
                // int percent_incorrect = (count_incorrect/max)*100;
                // cerr << endl;
                // cerr << percent_incorrect << "% percent with zero haplotypes matched from haplotype pair" <<endl;
                // cerr << percent_correct << "% percent with two matched haplotypes" <<endl;
                // cerr << percent_half_correct << "% percent with one haplotype mathced from haplotype pair" <<endl;
                
                
            }
            
        }

        TEST_CASE("Test7"){
            //Returns optimal phased genome on a 8-node graph 2 nested snarls with 1 short read
            SECTION("Test7:Requires haplotype pair to match truth set") {   
                VG graph;
				
                
                Node* n1 = graph.create_node("GGG"); //gets ID in incremintal order 
                Node* n2 = graph.create_node("CCC");
                Node* n3 = graph.create_node("A");
                Node* n4 = graph.create_node("T");
                Node* n5 = graph.create_node("G");
                Node* n6 = graph.create_node("CTGG");
                Node* n7 = graph.create_node("TAC");
                Node* n8 = graph.create_node("C");

                path_handle_t path_handle = graph.create_path_handle("x");
                graph.append_step(path_handle, graph.get_handle(n1->id()));
                graph.append_step(path_handle, graph.get_handle(n7->id()));
                graph.append_step(path_handle, graph.get_handle(n8->id()));
                graph.append_step(path_handle, graph.get_handle(n6->id()));

                
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
                xg::XG xg_index; 
                xg_index.from_path_handle_graph(graph);              

                // Make a multipath mapper to map against the graph.
                MultipathMapper multipath_mapper(&xg_index, gcsaidx, lcpidx); 
                
                vector<string> reads = {"GGGCCCAGCTGG"};
                vector<Alignment> alns = {reads.size(), Alignment()};

                // set alignment sequence
                for(int i = 0; i< reads.size(); i++){
                    alns[i].set_sequence(reads[i]);
                }

               MCMCGenotyper mcmc_genotyper = MCMCGenotyper(snarl_manager, graph, n_iterations, seed);
               vector<MultipathAlignment> multipath_aln_vector = vector<MultipathAlignment>(); 

               vector<vector<MultipathAlignment>> vect = {reads.size(),vector<MultipathAlignment>() };
                    
                    
                // map read in alignment to graph and make multipath alignments 
                for(int i = 0; i< reads.size(); i++){
                    multipath_mapper.multipath_map(alns[i], vect[i], 1);
                }


                // accumulate the mapped reads in one vector
                for(int i = 0; i< reads.size(); i++){
                    move(vect[i].begin(), vect[i].end(), back_inserter(multipath_aln_vector)); 
                }
                    
                double log_base = gssw_dna_recover_log_base(1,4,.5,1e-12);
                unique_ptr<PhasedGenome> genome = mcmc_genotyper.run_genotype(multipath_aln_vector, log_base);

                 // create a set of possible solutions
                vector<NodeTraversal> soln1;
                soln1 = {NodeTraversal(n1), NodeTraversal(n2), NodeTraversal(n3), NodeTraversal(n5), NodeTraversal(n6)};

                set<vector<NodeTraversal>> solns_set;
                solns_set.insert(soln1); 

                // check requirements 
                REQUIRE(genome->num_haplotypes() == 2);

                // move the genome haplotype into a vector
                vector<NodeTraversal> haplotype1, haplotype2;
                copy(genome->begin(0), genome->end(0), back_inserter(haplotype1));
                copy(genome->begin(1), genome->end(1), back_inserter(haplotype2));
                
                vector<vector<NodeTraversal>> haplos = {haplotype1, haplotype2};
                
                bool pass = false;
                if(solns_set.count(haplotype1) || solns_set.count(haplotype2)){
                    pass = true;
                }
                // requires at least one of the haplotypes to match truth set 
                // only one match is required because there is only one read
                REQUIRE(pass);
            
            }
        }

    }

}


