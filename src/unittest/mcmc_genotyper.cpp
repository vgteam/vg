/// \file mcmc_genotyper.cpp
///  
/// unit tests for mcmc genotyper construction and utility functions
///
#include <list>
#include <xg.hpp>
#include <stdio.h>
#include <iostream>
#include "../multipath_mapper.hpp"
#include <gbwt/dynamic_gbwt.h>
#include "../build_index.hpp"
#include "../haplotypes.hpp"
#include "vg/io/json2pb.h"
#include <vg/vg.pb.h>
#include "../vg.hpp"
#include "../multipath_alignment.hpp"
#include "../utility.hpp"
#include "../mcmc_genotyper.hpp"
#include "../snarls.hpp"
#include "../integrated_snarl_finder.hpp"
#include "../genotypekit.hpp"
#include "algorithms/min_cut_graph.hpp"
#include "../cactus_snarl_finder.hpp"

#include "catch.hpp"

// #define debug

// #define debug_snarl_graph

namespace vg {
    namespace unittest {
        
        const int seed = 0;
        const int n_iterations = 100;
        const int burn_in = n_iterations/2;
        const int gamma_freq = 50;

        
        TEST_CASE("Returns optimal phased genome on a 1-node graph with 1 short read ") {
            
            SECTION("Test1: Requires haplotype pair to match truth set") {
                VG graph;
				
                Node* n1 = graph.create_node("GCA");
#ifdef debug
                //cerr << "this is the node address"<< &n1 << endl;
#endif
                

                // name a path 
                path_handle_t path_handle = graph.create_path_handle("x");
                // append a visit to a node (via Tests run_genotyper() that takes handle) to the given path 
                graph.append_step(path_handle, graph.get_handle(n1->id()));
                
                //Cactus does not currently support finding snarls in graph of single-node connected components
                SnarlManager snarl_manager; 
            
                string read = string("GCA");
                multipath_alignment_t multipath_aln;
                multipath_aln.set_sequence(read);
                
                 // add subpaths
                subpath_t* subpath0 = multipath_aln.add_subpath();
                 
                // designate mappings
                path_mapping_t* mapping0 = subpath0->mutable_path()->add_mapping();
                mapping0->mutable_position()->set_node_id(1);
                edit_t* edit0 = mapping0->add_edit();
                edit0->set_from_length(3); //a match 
                edit0->set_to_length(3);

                subpath0->set_score(1);

                multipath_aln.add_start(0);        

                MCMCGenotyper mcmc_genotyper = MCMCGenotyper(snarl_manager, graph, n_iterations, seed, burn_in, gamma_freq);   
                vector<multipath_alignment_t> multipath_aln_vector = vector<multipath_alignment_t>({multipath_aln});
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
        TEST_CASE("Returns optimal phased genome on a 4-node graph and snarl with 1 short read"){
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
				
				IntegratedSnarlFinder bubble_finder(graph);
                SnarlManager snarl_manager = bubble_finder.find_snarls_parallel();
                
                string read = string("GCATCTGA");
                multipath_alignment_t multipath_aln;
                multipath_aln.set_sequence(read);
                
                // add subpaths with same topology as graph
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
                
                // set scores 
                subpath0->set_score(1);
                subpath1->set_score(1);
                subpath2->set_score(-4);
                subpath3->set_score(1);

                // designate mappings
                path_mapping_t* mapping0 = subpath0->mutable_path()->add_mapping();
                mapping0->mutable_position()->set_node_id(1);
                edit_t* edit0 = mapping0->add_edit();
                edit0->set_from_length(3); //a match 
                edit0->set_to_length(3);
                
                path_mapping_t* mapping1 = subpath1->mutable_path()->add_mapping();
                mapping1->mutable_position()->set_node_id(2);
                edit_t* edit1 = mapping1->add_edit();
                edit1->set_from_length(1); //a match 
                edit1->set_to_length(1);
                
                path_mapping_t* mapping2 = subpath2->mutable_path()->add_mapping();
                mapping2->mutable_position()->set_node_id(3);
                edit_t* edit2 = mapping2->add_edit();
                edit2->set_from_length(1); //a snip
                edit2->set_to_length(1);
                edit2->set_sequence("T"); 

                path_mapping_t* mapping3 = subpath3->mutable_path()->add_mapping();
                mapping3->mutable_position()->set_node_id(4);
                edit_t* edit3 = mapping3->add_edit();
                edit3->set_from_length(4); //a match 
                edit3->set_to_length(4);
                
                multipath_aln.add_start(0);

            
                MCMCGenotyper mcmc_genotyper = MCMCGenotyper(snarl_manager, graph, n_iterations, seed, burn_in, gamma_freq);
                vector<multipath_alignment_t> multipath_aln_vector = vector<multipath_alignment_t>({multipath_aln});
                double log_base = gssw_dna_recover_log_base(1,4,.5,1e-12);
                unique_ptr<PhasedGenome> genome = mcmc_genotyper.run_genotype(multipath_aln_vector, log_base);
#ifdef debug
                //cerr << "**************" <<endl;
                //genome->print_phased_genome();
#endif
                
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
        TEST_CASE("Returns optimal phased genome on a 7-node graph two connected snarls with 1 short read"){
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
				
                IntegratedSnarlFinder bubble_finder(graph);
                SnarlManager snarl_manager = bubble_finder.find_snarls_parallel();

                // Make GCSA quiet
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

               MCMCGenotyper mcmc_genotyper = MCMCGenotyper(snarl_manager, graph, n_iterations, seed,burn_in, gamma_freq);
               vector<multipath_alignment_t> multipath_aln_vector = vector<multipath_alignment_t>();

               vector<vector<multipath_alignment_t>> vect = {reads.size(),vector<multipath_alignment_t>() };
                    
                    
                // map read in alignment to graph and make multipath alignments 
                for(int i = 0; i< reads.size(); i++){
                    multipath_mapper.multipath_map(alns[i], vect[i]);
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
        TEST_CASE("Returns optimal phased genome on a 8-node graph 2 nested snarls with 1 read"){
            SECTION("Test4:Requires haplotype pair to match truth set") {   
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
                graph.append_step(path_handle, graph.get_handle(n2->id()));
                graph.append_step(path_handle, graph.get_handle(n4->id()));
                graph.append_step(path_handle, graph.get_handle(n5->id()));
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
				
                // IntegratedSnarlFinder bubble_finder(graph);
                // SnarlManager snarl_manager = bubble_finder.find_snarls_parallel();

                CactusSnarlFinder bubble_finder(graph);
                SnarlManager snarl_manager = bubble_finder.find_snarls();

                // Make GCSA quiet
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
               

               MCMCGenotyper mcmc_genotyper = MCMCGenotyper(snarl_manager, graph, n_iterations, seed, burn_in, gamma_freq);
               vector<multipath_alignment_t> multipath_aln_vector = vector<multipath_alignment_t>(); 

               vector<vector<multipath_alignment_t>> vect = {reads.size(),vector<multipath_alignment_t>() };
                    
                    
                // map read in alignment to graph and make multipath alignments 
                for(int i = 0; i< reads.size(); i++){
                    multipath_mapper.multipath_map(alns[i], vect[i]);
                }


                // accumulate the mapped reads in one vector
                for(int i = 0; i< reads.size(); i++){
                    move(vect[i].begin(), vect[i].end(), back_inserter(multipath_aln_vector)); 
                }
#ifdef debug

                for (const auto& subpath : multipath_aln_vector[0].subpath()) {
                    if(subpath.has_path()){
                        auto& path = subpath.path();
                        //for every mapping in the path
                        for(size_t i = 0; i < path.mapping_size(); i++){
                            auto& mapping = path.mapping(i);
                            int64_t node_id = mapping.position().node_id();
                            bool is_reverse = mapping.position().is_reverse(); 
                            cerr << node_id << "->";
                        }
                        cerr << endl;
                    }
                }        
#endif 
                double log_base = gssw_dna_recover_log_base(1,4,.5,1e-12);
                unique_ptr<PhasedGenome> genome = mcmc_genotyper.run_genotype(multipath_aln_vector, log_base);
                
#ifdef debug
                genome->print_phased_genome();
#endif
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
        TEST_CASE("Returns optimal phased genome on a 7-node graph containing 2 connected snarls and 4 mapped reads"){
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
				
				IntegratedSnarlFinder bubble_finder(graph);
                SnarlManager snarl_manager = bubble_finder.find_snarls_parallel();

                // Make GCSA quiet
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
                
                
                MCMCGenotyper mcmc_genotyper = MCMCGenotyper(snarl_manager, graph, n_iterations, seed, burn_in, gamma_freq);
                vector<multipath_alignment_t> multipath_aln_vector = vector<multipath_alignment_t>();

                vector<vector<multipath_alignment_t>> vect = {reads.size(),vector<multipath_alignment_t>() };
                    
                    
                // map read in alignment to graph and make multipath alignments 
                for(int i = 0; i< reads.size(); i++){
                    multipath_mapper.multipath_map(alns[i], vect[i]);
                }
                    
                // accumulate the mapped reads in one vector
                for(int i = 0; i< reads.size(); i++){
                    move(vect[i].begin(), vect[i].end(), back_inserter(multipath_aln_vector)); 
                }

                double log_base = gssw_dna_recover_log_base(1,4,.5,1e-12);
                
                //pass vector with accumulated multipath_alignment_t objects to run_genotype()
                unique_ptr<PhasedGenome> genome = mcmc_genotyper.run_genotype(multipath_aln_vector, log_base); 

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
                
                set<vector<NodeTraversal>> phased_genome = {haplotype1, haplotype2};
                set<vector<NodeTraversal>> rev_phased_genome = {haplotype2, haplotype1};



                // check haplotypes are indeed in optimal solution set  
                bool pass = false; 
                if(phased_genome == solns_set || rev_phased_genome == solns_set){
                    pass = true;
                }
                REQUIRE(pass);
                

                // Clean up the GCSA/LCP index
                delete gcsaidx;
                delete lcpidx;
            }
        }
        TEST_CASE("Returns optimal phased genome on a 7-node graph and two connected snarls with 8 short read"){
            SECTION("Test6: Requires haplotype pair to match truth set"){
                
                vector<int> failed_seeds;
                double count_correct = 0.0;
                double count_incorrect = 0.0;
                double count_minus = 0.0;
                vector<double> results = vector<double>();
                
                int num_iterations = 1000;
                int max = 30;
                int burn = 500;
                int freq = 500;

                
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
                    
                    IntegratedSnarlFinder bubble_finder(graph);
                    SnarlManager snarl_manager = bubble_finder.find_snarls_parallel();

                    // Make GCSA quiet
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
                    
                    MCMCGenotyper mcmc_genotyper = MCMCGenotyper(snarl_manager, graph, num_iterations, seed_i, burn, freq);
                    vector<multipath_alignment_t> multipath_aln_vector = vector<multipath_alignment_t>();

                    vector<vector<multipath_alignment_t>> vect = {reads.size(),vector<multipath_alignment_t>() };
                    
                    
                    // map read in alignment to graph and make multipath alignments 
                    for(int i = 0; i< reads.size(); i++){
                        multipath_mapper.multipath_map(alns[i], vect[i]);
                    }
                    
                    // accumulate the mapped reads in one vector
                    for(int i = 0; i< reads.size(); i++){
                        move(vect[i].begin(), vect[i].end(), back_inserter(multipath_aln_vector)); 
                    }
                    

                    double log_base = gssw_dna_recover_log_base(1,4,.5,1e-12);
                    //pass vector with accumulated multipath_alignment_t objects to run_genotype()
                    unique_ptr<PhasedGenome> genome = mcmc_genotyper.run_genotype(multipath_aln_vector, log_base); 
                    

                        // create a set of 2 possible solutions
                    vector<NodeTraversal> soln1, soln2, soln3;
                    soln1 = {NodeTraversal(n1), NodeTraversal(n2), NodeTraversal(n4), NodeTraversal(n6), NodeTraversal(n7)};
                    soln2 = {NodeTraversal(n1), NodeTraversal(n3), NodeTraversal(n4), NodeTraversal(n5), NodeTraversal(n7)};
                    soln3 = {NodeTraversal(n1),NodeTraversal(n3), NodeTraversal(n4), NodeTraversal(n6), NodeTraversal(n7)};

                    set<vector<NodeTraversal>> solns_set1, solns_set2;
                    solns_set1.insert(soln1);
                    solns_set1.insert(soln2);
                    
                    solns_set2.insert(soln2);
                    solns_set2.insert(soln3);


                    REQUIRE(genome->num_haplotypes() == 2);

                    // move the genome haplotype into a vector
                    vector<NodeTraversal> haplotype1, haplotype2;
                    copy(genome->begin(0), genome->end(0), back_inserter(haplotype1));
                    copy(genome->begin(1), genome->end(1), back_inserter(haplotype2));
                
                    set<vector<NodeTraversal>> phased_genome = {haplotype1, haplotype2};
                    set<vector<NodeTraversal>> rev_phased_genome = {haplotype2, haplotype1};
                     
                    if(genome->num_haplotypes() == 2){
                        if((phased_genome == solns_set1 || rev_phased_genome == solns_set1) || (phased_genome == solns_set2 || rev_phased_genome == solns_set2)){
                            // both haplotypes in phased genome are in the soln set
                            // prevents false positives
                            count_correct++;
                        }
                        else{
                            // one of the haplotypes from the phased genome are in the solution set
                            // or if haplo1=haplo2 then both haplotypes are in the solution set 
                            count_incorrect++;
                            failed_seeds.push_back(seed_i);
                        }
                    }
                    
                    // Clean up the GCSA/LCP index
                    delete gcsaidx;
                    delete lcpidx;
                }
#ifdef debug
                // cerr <<"****************************DONE TESTING****************************" << endl;
                // cerr << count_correct << " tests" << " out of " << max <<" matched both haplotypes from the haplotype pair " <<endl;
                // cerr << count_incorrect << " tests " << "out of " << max << " matched one haplotype from haplotype pair  " <<endl;
                // int percent_incorrect = (count_incorrect/max)*100;
                // int percent_correct = (count_correct/max)*100;  

                // cerr << endl;
            
                // cerr << percent_correct << "% percent with two matched haplotypes" <<endl;
                // cerr << percent_incorrect << "% percent with one haplotype mathced from haplotype pair" <<endl;
                // cerr <<endl;
                // cerr << failed_seeds.size() <<" seeds failed testing : "<<endl;
                // for(int i = 0; i < failed_seeds.size(); i++ ){
                //     cerr<< failed_seeds[i] << ", ";
                // }
                // cerr << endl;
#endif
    
            }
            
        }
        TEST_CASE("Returns optimal phased genome on a 7-node graph and two connected snarls with 9 short read"){
            SECTION("Test7: Requires haplotype pair to match truth set"){
                
                vector<int> failed_seeds;
                double count_correct = 0.0;
                double count_incorrect = 0.0;
                double count_half_correct = 0.0;
                double count_minus = 0.0;
                vector<double> results = vector<double>();
                
                int num_iterations = 500;
                int seed_i = std::chrono::system_clock::now().time_since_epoch().count();
                int burn = 250;
                int freq = 250;
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
                
                IntegratedSnarlFinder bubble_finder(graph);
                SnarlManager snarl_manager = bubble_finder.find_snarls_parallel();

                // Make GCSA quiet
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

                
                vector<string> reads = {"GCATCTGAGCCC","GCATCTGAGCCC", "GCATCTGAGCCC", "GCAGCTGAACCC", "GCAGCTGAACCC","GCAGCTGAACCC", "GCAGCTGAACCC", "GCAGCTGAGCCC", "GCAGCTGAGCCC","GCATCTGAACCC" };
                vector<Alignment> alns = {reads.size(), Alignment()};
                
                // set alignment sequence
                for(int i = 0; i< reads.size(); i++){
                    alns[i].set_sequence(reads[i]);
                }
                
                MCMCGenotyper mcmc_genotyper = MCMCGenotyper(snarl_manager, graph, num_iterations, seed_i, burn, freq);
                vector<multipath_alignment_t> multipath_aln_vector = vector<multipath_alignment_t>(); 

                vector<vector<multipath_alignment_t>> vect = {reads.size(),vector<multipath_alignment_t>() };
                
                
                // map read in alignment to graph and make multipath alignments 
                for(int i = 0; i< reads.size(); i++){
                    multipath_mapper.multipath_map(alns[i], vect[i]);
                }
                
                // accumulate the mapped reads in one vector
                for(int i = 0; i< reads.size(); i++){
                    move(vect[i].begin(), vect[i].end(), back_inserter(multipath_aln_vector)); 
                }
                
                double log_base = gssw_dna_recover_log_base(1,4,.5,1e-12);
                //pass vector with accumulated MultipathAlignment objects to run_genotype()
                unique_ptr<PhasedGenome> genome = mcmc_genotyper.run_genotype(multipath_aln_vector, log_base); 
                

                // create a set of possible solutions
                vector<NodeTraversal> soln1, soln2, soln3;
                soln1 = {NodeTraversal(n1), NodeTraversal(n2), NodeTraversal(n4), NodeTraversal(n6), NodeTraversal(n7)};
                soln2 = {NodeTraversal(n1), NodeTraversal(n3), NodeTraversal(n4), NodeTraversal(n5), NodeTraversal(n7)};
                soln3 = {NodeTraversal(n1),NodeTraversal(n3), NodeTraversal(n4), NodeTraversal(n6), NodeTraversal(n7)};

                set<vector<NodeTraversal>> solns_set1, solns_set2;
                solns_set1.insert(soln1);
                solns_set1.insert(soln2);
                
                solns_set2.insert(soln2);
                solns_set2.insert(soln3);


                REQUIRE(genome->num_haplotypes() == 2);

                // move the genome haplotype into a vector
                vector<NodeTraversal> haplotype1, haplotype2;
                copy(genome->begin(0), genome->end(0), back_inserter(haplotype1));
                copy(genome->begin(1), genome->end(1), back_inserter(haplotype2));
                
                set<vector<NodeTraversal>> phased_genome = {haplotype1, haplotype2};
                set<vector<NodeTraversal>> rev_phased_genome = {haplotype2, haplotype1};

                bool pass = false;
                if( (phased_genome == solns_set1 || rev_phased_genome == solns_set1) || (phased_genome == solns_set2 || rev_phased_genome == solns_set2) ){
                    pass = true;
                }
                // check if the haplotype set is equal to solution set 
                // take into account the haplotypes do not have inherent ordering 
                REQUIRE(pass); 
                 
                // Clean up the GCSA/LCP index
                delete gcsaidx;
                delete lcpidx;    
                
            }
            
        }
        TEST_CASE("snarl_map and snarl_graph work on 7 node graph") {

                VG graph;

                Node* n1 = graph.create_node("GCA"); 
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("G");
                Node* n4 = graph.create_node("CTGA"); 
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
				
				IntegratedSnarlFinder bubble_finder(graph);
                SnarlManager snarl_manager = bubble_finder.find_snarls_parallel();

                // Make GCSA quiet
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
                
                MCMCGenotyper mcmc_genotyper = MCMCGenotyper(snarl_manager, graph, n_iterations, seed, burn_in, gamma_freq);

                vector<multipath_alignment_t> multipath_aln_vector = vector<multipath_alignment_t>();

                vector<vector<multipath_alignment_t>> vect = {reads.size(),vector<multipath_alignment_t>() };
                
                
                // map read in alignment to graph and make multipath alignments 
                for(int i = 0; i< reads.size(); i++){
                    multipath_mapper.multipath_map(alns[i], vect[i]);
                }
                
                // accumulate the mapped reads in one vector
                for(int i = 0; i< reads.size(); i++){
                    move(vect[i].begin(), vect[i].end(), back_inserter(multipath_aln_vector)); 
                }
                
                
                // generate phased genome with haplotypes that have diff allele variants at snarl sites
                // unique_ptr<PhasedGenome> phased_genome = mcmc_genotyper.generate_initial_guess();
                unique_ptr<PhasedGenome> phased_genome(new PhasedGenome(snarl_manager));

                vector<NodeTraversal> haplotype_1;
                vector<NodeTraversal> haplotype_2;
                
                haplotype_1.push_back(NodeTraversal(n1));
                haplotype_1.push_back(NodeTraversal(n3));
                haplotype_1.push_back(NodeTraversal(n4));
                haplotype_1.push_back(NodeTraversal(n5));
                haplotype_1.push_back(NodeTraversal(n7));
                
                haplotype_2.push_back(NodeTraversal(n1));
                haplotype_2.push_back(NodeTraversal(n2));
                haplotype_2.push_back(NodeTraversal(n4));
                haplotype_2.push_back(NodeTraversal(n6));
                haplotype_2.push_back(NodeTraversal(n7));

                // construct haplotypes
                phased_genome->add_haplotype(haplotype_1.begin(), haplotype_1.end());
                phased_genome->add_haplotype(haplotype_2.begin(), haplotype_2.end());

                // index sites
                phased_genome->build_indices();
                
                unordered_map<pair<const Snarl*, const Snarl*>, int32_t> snarl_map = mcmc_genotyper.make_snarl_map(multipath_aln_vector ,  *phased_genome);
#ifdef debug_snarl_graph
                unordered_map<pair<const Snarl*, const Snarl*>, int32_t>::iterator it = snarl_map.begin();
                while(it != snarl_map.end())
                {
                    std::cout<<"weight: " << it->second<<std::endl;
                    it++;
                }
#endif
                algorithms::Graph snarl_graph = mcmc_genotyper.make_snarl_graph(snarl_map);

#ifdef debug_snarl_graph
                cout << "graph" <<endl;        
                
                for(size_t i =0; i < snarl_graph.get_node_ids().size(); i++){
                    cout << "node size " <<snarl_graph.get_node_ids().size() << endl; 
                    vector<size_t> node_ids = snarl_graph.get_node_ids();
                    
                    for (size_t k =0; k < node_ids.size(); k++){
                        cout << "node id : " <<node_ids[k] << endl; 
                    }
                    for(size_t j =0; j < snarl_graph.get_node_by_id(node_ids[i]).edges.size(); j++){
                        cout << "node "<< node_ids[i] <<" size of edges " <<snarl_graph.get_node_ids().size() << endl; 
                        algorithms::Node& node =  snarl_graph.get_node_by_id(node_ids[i]);
                        cout << "node "<<node_ids[i] <<"->" << node.edges[j].other <<endl;
                        cout << "edge weight: " << node.edges[j].weight <<endl;

                    }
                    
                }
#endif




        }
        TEST_CASE("snarl_map and snarl_graph work on 14 node graph") {
            VG graph;
				
            Node* n1 = graph.create_node("GGG");  
            Node* n2 = graph.create_node("CCC");
            Node* n3 = graph.create_node("A");
            Node* n4 = graph.create_node("T");
            Node* n5 = graph.create_node("G");
            Node* n6 = graph.create_node("CTGG");
            Node* n7 = graph.create_node("TAC");
            Node* n8 = graph.create_node("C");
            Node* n9 = graph.create_node("T");
            Node* n10 = graph.create_node("G");
            Node* n11 = graph.create_node("CTGA");
            Node* n12 = graph.create_node("A");
            Node* n13 = graph.create_node("G");
            Node* n14 = graph.create_node("CCC");

            path_handle_t path_handle = graph.create_path_handle("x");
            graph.append_step(path_handle, graph.get_handle(n1->id()));
            graph.append_step(path_handle, graph.get_handle(n7->id()));
            graph.append_step(path_handle, graph.get_handle(n8->id()));
            graph.append_step(path_handle, graph.get_handle(n6->id()));
            graph.append_step(path_handle, graph.get_handle(n9->id()));
            graph.append_step(path_handle, graph.get_handle(n11->id()));
            graph.append_step(path_handle, graph.get_handle(n12->id()));
            graph.append_step(path_handle, graph.get_handle(n14->id()));

            
            graph.create_edge(n1, n2);
            graph.create_edge(n1, n7);
            graph.create_edge(n2, n3);
            graph.create_edge(n2, n4);
            graph.create_edge(n3, n5);
            graph.create_edge(n4, n5);
            graph.create_edge(n5, n6);
            graph.create_edge(n7, n8);
            graph.create_edge(n8, n6);
            graph.create_edge(n6, n9);
            graph.create_edge(n6, n10);
            graph.create_edge(n9, n11);
            graph.create_edge(n10, n11);
            graph.create_edge(n11, n12);
            graph.create_edge(n11, n13);
            graph.create_edge(n12, n14);
            graph.create_edge(n13, n14);
            
            IntegratedSnarlFinder bubble_finder(graph);
            SnarlManager snarl_manager = bubble_finder.find_snarls_parallel();

            // Make GCSA quiet
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
            
            vector<string> reads = {"GGGCCCAGCTGG", "GGGCCCAGCTGGTCTGAGCCC", "GGGCCCAGCTGGTCTGAGCCC", 
                                    "GGGTACCCTGGTCTGAGCCC", "CTGGTCTGAGCCC", "GGGCCCTGCTGGGCTGAGCCC", 
                                    "GGGCCCTGCTGGGCTGAGCCC", "GGGCCCTGCTGGGCTGAGCCC", "GGGCCCAGCTGGTCTGAACCC",
                                     "GGGCCCAGCTGGTCTGAACCC"};
            vector<Alignment> alns = {reads.size(), Alignment()};

            // set alignment sequence
            for(int i = 0; i< reads.size(); i++){
                alns[i].set_sequence(reads[i]);
            }

            MCMCGenotyper mcmc_genotyper = MCMCGenotyper(snarl_manager, graph, n_iterations, seed, burn_in, gamma_freq);
            vector<multipath_alignment_t> multipath_aln_vector = vector<multipath_alignment_t>(); 

            vector<vector<multipath_alignment_t>> vect = {reads.size(),vector<multipath_alignment_t>() };
                
                
            // map read in alignment to graph and make multipath alignments 
            for(int i = 0; i< reads.size(); i++){
                multipath_mapper.multipath_map(alns[i], vect[i]);
            }


            // accumulate the mapped reads in one vector
            for(int i = 0; i< reads.size(); i++){
                move(vect[i].begin(), vect[i].end(), back_inserter(multipath_aln_vector)); 
            }

            // generate phased genome with haplotypes that have diff allele variants at snarl sites
            // unique_ptr<PhasedGenome> phased_genome = mcmc_genotyper.generate_initial_guess();
            unique_ptr<PhasedGenome> phased_genome(new PhasedGenome(snarl_manager));

            vector<NodeTraversal> haplotype_1;
            vector<NodeTraversal> haplotype_2;
            
            // haplotype_1.push_back(NodeTraversal(n1));
            // haplotype_1.push_back(NodeTraversal(n7));
            // haplotype_1.push_back(NodeTraversal(n8));
            // haplotype_1.push_back(NodeTraversal(n6));
            // haplotype_1.push_back(NodeTraversal(n10));
            // haplotype_1.push_back(NodeTraversal(n11));
            // haplotype_1.push_back(NodeTraversal(n13));
            // haplotype_1.push_back(NodeTraversal(n14));
            haplotype_1.push_back(NodeTraversal(n1));
            haplotype_1.push_back(NodeTraversal(n2));
            haplotype_1.push_back(NodeTraversal(n4));
            haplotype_1.push_back(NodeTraversal(n5));
            haplotype_1.push_back(NodeTraversal(n6));
            haplotype_1.push_back(NodeTraversal(n10));
            haplotype_1.push_back(NodeTraversal(n11));
            haplotype_1.push_back(NodeTraversal(n13));
            haplotype_1.push_back(NodeTraversal(n14));

            haplotype_2.push_back(NodeTraversal(n1));
            haplotype_2.push_back(NodeTraversal(n2));
            haplotype_2.push_back(NodeTraversal(n3));
            haplotype_2.push_back(NodeTraversal(n5));
            haplotype_2.push_back(NodeTraversal(n6));
            haplotype_2.push_back(NodeTraversal(n9));
            haplotype_2.push_back(NodeTraversal(n11));
            haplotype_2.push_back(NodeTraversal(n12));
            haplotype_2.push_back(NodeTraversal(n14));

            // construct haplotypes
            phased_genome->add_haplotype(haplotype_1.begin(), haplotype_1.end());
            phased_genome->add_haplotype(haplotype_2.begin(), haplotype_2.end());

            // index sites
            phased_genome->build_indices();

            unordered_map<pair<const Snarl*, const Snarl*>, int32_t> snarl_map = mcmc_genotyper.make_snarl_map(multipath_aln_vector ,  *phased_genome);
#ifdef debug_snarl_graph
            cout << "map" <<endl;        
            cout << "******************************************************"<<endl;
            unordered_map<pair<const Snarl*, const Snarl*>, int32_t>::iterator it = snarl_map.begin();
            while(it != snarl_map.end())
            {
                std::cout<<"weight: " << it->second<<std::endl;
                it++;
            }
            cout << "******************************************************"<<endl;
#endif
            algorithms::Graph snarl_graph = mcmc_genotyper.make_snarl_graph(snarl_map);

#ifdef debug_snarl_graph
            cout << "******************************************************"<<endl;
            cout << "graph" <<endl; 
            
            vector<size_t> node_ids = snarl_graph.get_node_ids();
            size_t num_nodes = node_ids.size();
            cout << "node size " <<num_nodes << endl; 
            for (size_t k =0; k < node_ids.size(); k++){
                    cout << "node id : " <<node_ids[k] << endl; 
            }       
            
            for(size_t i =0; i < num_nodes; i++){
                size_t id = node_ids[i];
                algorithms::Node& node_i = snarl_graph.get_node_by_id(id);
                for(size_t j =0; j < node_i.edges.size(); j++){
                    cout << "node "<< id <<" size of edges " <<node_i.edges.size() << endl; 
                    cout << "node "<< id <<"->" << node_i.edges[j].other <<endl;
                    cout << "edge weight: " << node_i.edges[j].weight <<endl;

                }
                
            }
            cout << "******************************************************"<<endl;
#endif

        }
        TEST_CASE("run mcmc_genotyper with karger-stein min cut and alt. proposal") {
            VG graph;
				
            Node* n1 = graph.create_node("GGG");  
            Node* n2 = graph.create_node("CCC");
            Node* n3 = graph.create_node("A");
            Node* n4 = graph.create_node("T");
            Node* n5 = graph.create_node("G");
            Node* n6 = graph.create_node("CTGG");
            Node* n7 = graph.create_node("TAC");
            Node* n8 = graph.create_node("C");
            Node* n9 = graph.create_node("T");
            Node* n10 = graph.create_node("G");
            Node* n11 = graph.create_node("CTGA");
            Node* n12 = graph.create_node("A");
            Node* n13 = graph.create_node("G");
            Node* n14 = graph.create_node("CCC");

            path_handle_t path_handle = graph.create_path_handle("x");
            graph.append_step(path_handle, graph.get_handle(n1->id()));
            graph.append_step(path_handle, graph.get_handle(n7->id()));
            graph.append_step(path_handle, graph.get_handle(n8->id()));
            graph.append_step(path_handle, graph.get_handle(n6->id()));
            graph.append_step(path_handle, graph.get_handle(n9->id()));
            graph.append_step(path_handle, graph.get_handle(n11->id()));
            graph.append_step(path_handle, graph.get_handle(n12->id()));
            graph.append_step(path_handle, graph.get_handle(n14->id()));

            
            graph.create_edge(n1, n2);
            graph.create_edge(n1, n7);
            graph.create_edge(n2, n3);
            graph.create_edge(n2, n4);
            graph.create_edge(n3, n5);
            graph.create_edge(n4, n5);
            graph.create_edge(n5, n6);
            graph.create_edge(n7, n8);
            graph.create_edge(n8, n6);
            graph.create_edge(n6, n9);
            graph.create_edge(n6, n10);
            graph.create_edge(n9, n11);
            graph.create_edge(n10, n11);
            graph.create_edge(n11, n12);
            graph.create_edge(n11, n13);
            graph.create_edge(n12, n14);
            graph.create_edge(n13, n14);
            
            IntegratedSnarlFinder bubble_finder(graph);
            SnarlManager snarl_manager = bubble_finder.find_snarls_parallel();

            // Make GCSA quiet
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
            
            vector<string> reads = {"GGGCCCAGCTGG", "GGGCCCAGCTGGTCTGAGCCC", "GGGCCCAGCTGGTCTGAGCCC", 
                                    "GGGTACCCTGGTCTGAGCCC", "CTGGTCTGAGCCC", "GGGCCCTGCTGGGCTGAGCCC", 
                                    "GGGCCCTGCTGGGCTGAGCCC", "GGGCCCTGCTGGGCTGAGCCC", "GGGCCCAGCTGGTCTGAACCC",
                                     "GGGCCCAGCTGGTCTGAACCC"};
            vector<Alignment> alns = {reads.size(), Alignment()};

            // set alignment sequence
            for(int i = 0; i< reads.size(); i++){
                alns[i].set_sequence(reads[i]);
            }

            MCMCGenotyper mcmc_genotyper = MCMCGenotyper(snarl_manager, graph, n_iterations, seed, burn_in, gamma_freq);
            vector<multipath_alignment_t> multipath_aln_vector = vector<multipath_alignment_t>(); 

            vector<vector<multipath_alignment_t>> vect = {reads.size(),vector<multipath_alignment_t>() };
                
                
            // map read in alignment to graph and make multipath alignments 
            for(int i = 0; i< reads.size(); i++){
                multipath_mapper.multipath_map(alns[i], vect[i]);
            }


            // accumulate the mapped reads in one vector
            for(int i = 0; i< reads.size(); i++){
                move(vect[i].begin(), vect[i].end(), back_inserter(multipath_aln_vector)); 
            }

            // generate phased genome with haplotypes that have diff allele variants at snarl sites
            unique_ptr<PhasedGenome> phased_genome(new PhasedGenome(snarl_manager));

            vector<NodeTraversal> haplotype_1;
            vector<NodeTraversal> haplotype_2;
            
            haplotype_1.push_back(NodeTraversal(n1));
            haplotype_1.push_back(NodeTraversal(n2));
            haplotype_1.push_back(NodeTraversal(n4));
            haplotype_1.push_back(NodeTraversal(n5));
            haplotype_1.push_back(NodeTraversal(n6));
            haplotype_1.push_back(NodeTraversal(n10));
            haplotype_1.push_back(NodeTraversal(n11));
            haplotype_1.push_back(NodeTraversal(n13));
            haplotype_1.push_back(NodeTraversal(n14));

            haplotype_2.push_back(NodeTraversal(n1));
            haplotype_2.push_back(NodeTraversal(n2));
            haplotype_2.push_back(NodeTraversal(n3));
            haplotype_2.push_back(NodeTraversal(n5));
            haplotype_2.push_back(NodeTraversal(n6));
            haplotype_2.push_back(NodeTraversal(n9));
            haplotype_2.push_back(NodeTraversal(n11));
            haplotype_2.push_back(NodeTraversal(n12));
            haplotype_2.push_back(NodeTraversal(n14));

            // construct haplotypes
            phased_genome->add_haplotype(haplotype_1.begin(), haplotype_1.end());
            phased_genome->add_haplotype(haplotype_2.begin(), haplotype_2.end());

            // index sites
            phased_genome->build_indices();
            unique_ptr<vector<unordered_set<size_t>>> gamma(nullptr);
            //stores the sites we swapped alleles at - returned by alt_proposal
            unique_ptr<unordered_set<size_t>> to_swap_back(nullptr);
            // build gamma set 
            gamma.reset(new vector<unordered_set<size_t>> (mcmc_genotyper.karger_stein(multipath_aln_vector, *phased_genome))); 

            to_swap_back.reset(new unordered_set<size_t> (mcmc_genotyper.alt_proposal_sample(*gamma, *phased_genome)));
            REQUIRE(!to_swap_back->empty());
            REQUIRE(to_swap_back->size() >1 );
            
            for(auto iter = to_swap_back->begin(); iter != to_swap_back->end(); ++iter){
                // cerr << "snarl num to swap backs" <<*iter << endl; 
            
            }


        }

    }

}


