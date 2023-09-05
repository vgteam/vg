#include "mcmc_caller.hpp"  
#include "graph_caller.hpp"
#include "algorithms/expand_context.hpp"
#include "memoizing_graph.hpp"
#include "phased_genome.hpp"
// #define debug

namespace vg {

    /**
     * MCMCCaller : Inherits from VCFOutputCaller    
     */         
    MCMCCaller::MCMCCaller(const PathPositionHandleGraph* path_position_handle_graph,
                           PhasedGenome& genome,
                           SnarlManager& snarl_manager,
                           const string& sample_name,
                           const vector<string>& ref_paths,
                           const vector<size_t>& ref_path_offsets,
                           const vector<size_t>& ref_path_lengths,
                           ostream& out_stream) :
        path_position_handle_graph(path_position_handle_graph), genome(genome), snarl_manager(snarl_manager), VCFOutputCaller(sample_name), 
        sample_name(sample_name), ref_paths(ref_paths), ref_path_offsets(ref_path_offsets),
        ref_path_lengths(ref_path_lengths), out_stream(out_stream) {
        
        if(path_position_handle_graph == nullptr){
            cerr << "graph is empty" <<endl;
            //exit not succesful
            exit(1);
        }
        
        /// keep track of offsets in the reference paths
        // used to offset variant reference position
        for (int i = 0; i < ref_paths.size(); ++i) {
            ref_offsets[ref_paths[i]] = i < ref_path_offsets.size() ? ref_path_offsets[i] : 0;
        }

    }
    
    MCMCCaller:: ~MCMCCaller(){

    }
    string MCMCCaller::vcf_header(const PathPositionHandleGraph& pph_graph, const vector<string>& ref_paths,
                                const vector<size_t>& contig_length_overrides) const {
      
        string header = VCFOutputCaller::vcf_header(pph_graph, ref_paths, contig_length_overrides);
        header += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
        header += "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
        header += "##SAMPLE=<ID=" + sample_name + ">\n";
        header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_name;
        assert(output_vcf.openForOutput(header));
        header += "\n";
        return header;
    }


    void MCMCCaller::call_top_level_snarls(bool recurse_on_fail) {

    // Used to recurse on children of parents that can't be called
    vector<const Snarl*> snarl_queue;

    // Run the snarl caller on a snarl, and queue up the children if it fails
    auto process_snarl = [&](const Snarl* snarl) {
        // cerr << "before call_snarl"<<endl;
        // cerr << pb2json(*snarl) <<endl;
        bool was_called = call_snarl(*snarl);
        // cerr << "after call_snarl"<<endl;
        if (!was_called && recurse_on_fail) {
            const vector<const Snarl*>& children = snarl_manager.children_of(snarl);
#pragma omp critical (snarl_queue)
            {
                snarl_queue.insert(snarl_queue.end(), children.begin(), children.end());
            }
        }
    };

    // Start with the top level snarls
    // snarl_manager.for_each_top_level_snarl_parallel(process_snarl);
    snarl_manager.for_each_top_level_snarl(process_snarl);

    // Then recurse on any children the snarl caller failed to handle
    while (!snarl_queue.empty()) {
        vector<const Snarl*> cur_queue;
        std::swap(snarl_queue, cur_queue);
#pragma omp parallel for
        for (int i = 0; i < cur_queue.size(); ++i) {
            process_snarl(cur_queue[i]);
        }
    }
  
    }

    bool MCMCCaller::call_snarl(const Snarl& snarl){
         // if we can't handle the snarl, then the GraphCaller framework will recurse on its children
        if (!is_traversable(snarl)) {
            cerr<< "snarl is not traversable" <<endl;
            return false;
        }
        
        PathTraversalFinder trav_finder(*path_position_handle_graph, snarl_manager, ref_paths);
        auto trav_results = trav_finder.find_path_traversals(snarl);
        vector<SnarlTraversal> ref_path = trav_results.first;


        
        //If it can't find any traversals, you can't output the snarl In VCF.
        if(ref_path.empty()){
            // continue the loop of snarl without printing VCF file
            cerr << "ref path empty" <<endl;
            return 0;
        }else{
            //In practice there should be only one reference path, 
            //so you can just choose the first one returned
            SnarlTraversal ref_trav = ref_path[0];
                  
            vector<pair<step_handle_t, step_handle_t> > steps = trav_results.second;
            //<start_step, end_step>
            pair<step_handle_t, step_handle_t> start_and_end_pair = steps[0];
            step_handle_t first_start_step = start_and_end_pair.first;

            
            function<string(const SnarlTraversal&)> trav_string = [&](const SnarlTraversal& trav) {
            string seq;
            for (int i = 0; i < trav.visit_size(); ++i) {
                seq += path_position_handle_graph->get_sequence(path_position_handle_graph->get_handle(trav.visit(i).node_id(), trav.visit(i).backward()));
            }
            return seq;
            };
            

            // set ref_path name and seq 
            path_handle_t path_handle = path_position_handle_graph->get_path_handle_of_step(first_start_step);
            string ref_path_name  = path_position_handle_graph->get_path_name(path_handle);
            // cerr << "ref_path_name " << ref_path_name << endl;
            //TODO: check that ref_path name is equal to snarl traversal name 
            string ref_path_seq = trav_string(ref_trav);
            
            
            //get haplotypes that pass snarl
            vector<id_t> haplos_pass_snarl = genome.get_haplotypes_with_snarl(&snarl);
            
            assert(!haplos_pass_snarl.empty());

#ifdef debug
    cerr <<endl;
    cerr << "REF " << ref_path_seq << endl;
    cerr << "POS "<<snarl.start().node_id() << "_"<< snarl.end().node_id() <<endl;
    if(haplos_pass_snarl.size()!=2){cerr << "less than two haplos pass this snarl" <<endl;};
#endif          
            //traversals of each haplotype through the snarl
            // vector<SnarlTraversal> haplo_travs = {haplos_pass_snarl.size(), SnarlTraversal()};
            vector<SnarlTraversal> haplo_travs;
            vector<vector<NodeTraversal>> haplo_nodes;
            vector<int> genotype;
            int match_ref = 0;
            int not_match_ref = 1;
            string prev_trav_seq = "";

            //build SnarlTraversal obj for each haplotype that passes through snarl 
            for(int i = 0; i <haplos_pass_snarl.size(); i++ ){

                SnarlTraversal fresh_trav;
                // get_allele() does not add the start of the snarl so we have to add it manually
                *fresh_trav.add_visit() = snarl.start();

                //push back NodeTraversal for haplotype id in ith position of array 
                haplo_nodes.push_back( genome.get_allele(snarl,haplos_pass_snarl[i]) );

                // iterate through NodeTraversal nodes
                for (auto iter = haplo_nodes[i].begin(); iter != haplo_nodes[i].end(); iter++){
                  
                    int64_t n_id = iter->node->id();
                    bool backward = iter->backward;
                    Visit* v = fresh_trav.add_visit(); 
                    v->set_node_id(n_id);
                    v->set_backward(backward);
#ifdef debug
   cerr << "node id " << n_id <<endl; 
#endif 
                    // *haplo_travs[i].add_visit() = *v;
                } 
                // get_allele() does not add the end of the snarl so we have to add it manually
                *fresh_trav.add_visit() = snarl.end();
                haplo_travs.push_back(fresh_trav);

                assert(!haplo_travs.empty());
               
                //get the sequence for the SnarlTraversal
                string trav_seq = trav_string(haplo_travs[i]);
#ifdef debug
   cerr << "haplo " << i  << " " << trav_seq << endl; 
#endif                  
                //compare the sequences 
                if(trav_seq.compare(ref_path_seq) != 0){
                    //strings don't match ref path
                    // check if h1=h0, if so index will be equal
                    if(i == 1  && prev_trav_seq.compare(trav_seq) == 0 ){
                        // if H1==H2, index will be equal
                        //cerr << "match 1|1" <<endl;
                        genotype.push_back(not_match_ref-1);
#ifdef debug
   cerr << "GT " << not_match_ref-1  << endl; 
#endif                         
                        
                    }else{
                        genotype.push_back(not_match_ref);
#ifdef debug
   cerr << "GT " << not_match_ref << endl; 
#endif   
                        not_match_ref++;
                        
                    }                        
                }else if(trav_seq.compare(ref_path_seq) == 0){
                    //strings match ref path 
                    genotype.push_back(match_ref); 
#ifdef debug
   cerr << "GT " << match_ref << endl; 
#endif 
                }
                prev_trav_seq = trav_seq; 
                
            }
            assert(!genotype.empty());

            //emit variant
            emit_variant(snarl, genotype, ref_trav, ref_path_name, haplo_travs);
        
        }
        return true;
    
    }

    void MCMCCaller::emit_variant(const Snarl& snarl, const vector<int>& genotype, SnarlTraversal ref_trav,
                                 const string& ref_path_name, const vector<SnarlTraversal>& haplo_travs) const{
        
        // convert traversal to string
        // function that converst SnarlTraversals to strings
        // usage: trav_string(SnarlTraversal)
        
        function<string(const SnarlTraversal&)> trav_string = [&](const SnarlTraversal& trav) {
            string seq;
            for (int i = 0; i < trav.visit_size(); ++i) {
                seq += path_position_handle_graph->get_sequence(path_position_handle_graph->get_handle(trav.visit(i).node_id(), trav.visit(i).backward()));
            }
            return seq;
        };

    
        vcflib::Variant out_variant;

        // when calling alt/alt, the reference traversal doesn't end up in called_traversals.
        // this should get changed, but in the meantime we add it back here (as we need it for
        // the VCF output)
        // udpate: the reference traversal will be there when re-genotyping, but we can leave this logic
        // in case we want to ever add an option to toggle this.
        vector<SnarlTraversal> site_traversals;
        vector<int> site_genotype;
        for (int i = 0; i < genotype.size(); ++i) {
            if (genotype[i] == 0) {
                // if haplo traversal matches the ref, add to container
                site_traversals.push_back(haplo_travs[i]);
                break;
            }
        }
        if(site_traversals.empty()){
            //if none of the haplotypes matched, get reference SnarlTraversal 
            // and convert to string
            site_traversals.push_back(ref_trav);
        }

        out_variant.ref = trav_string(ref_trav);
        
        // deduplicate alleles and compute the site traversals and genotype
        map<string, int> allele_to_gt;    
        allele_to_gt[out_variant.ref] = 0;    
        for (int i = 0; i < genotype.size(); ++i) {
            if (genotype[i] == 0) {
                site_genotype.push_back(0);
            } else {
                string allele_string = trav_string(haplo_travs[i]);
                if (allele_to_gt.count(allele_string)) {
                    site_genotype.push_back(allele_to_gt[allele_string]);
                } else {
                    site_traversals.push_back(haplo_travs[i]);
                    site_genotype.push_back(allele_to_gt.size());
                    allele_to_gt[allele_string] = site_genotype.back();
                }
            }
        }

        out_variant.alt.resize(allele_to_gt.size() - 1);
        out_variant.alleles.resize(allele_to_gt.size());
        for (auto& allele_gt : allele_to_gt) {
            if (allele_gt.second > 0) {
                out_variant.alt[allele_gt.second - 1] = allele_gt.first;
            }
            out_variant.alleles[allele_gt.second] = allele_gt.first;
        }

        // fill out the rest of the variant
        out_variant.sequenceName = ref_path_name;
        // +1 to convert to 1-based VCF
        out_variant.position = get_ref_position(snarl, ref_path_name).first + ref_offsets.find(ref_path_name)->second + 1;
        out_variant.id = std::to_string(snarl.start().node_id()) + "_" + std::to_string(snarl.end().node_id());
        out_variant.filter = "PASS";
        out_variant.updateAlleleIndexes();

        // add the genotype
        out_variant.format.push_back("GT");
        auto& genotype_vector = out_variant.samples[sample_name]["GT"];
        
        stringstream vcf_gt;
        for (int i = 0; i < site_genotype.size(); ++i) {
            vcf_gt << site_genotype[i];
            if (i != site_genotype.size() - 1) {
                vcf_gt << "/";
            }
        }
        genotype_vector.push_back(vcf_gt.str());
        
        // clean up the alleles to not have so many common prefixes
        flatten_common_allele_ends(out_variant, true);
        flatten_common_allele_ends(out_variant, false);

        // add variant to list 
        if (!out_variant.alt.empty()) {
            add_variant(out_variant);
        }
        

    }


    void MCMCCaller::update_vcf_info(const Snarl& snarl,
                                 const vector<SnarlTraversal>& traversals,
                                 const vector<int>& genotype,
                                 const string& sample_name,
                                 vcflib::Variant& variant) const{
        

    }

    bool MCMCCaller::is_traversable(const Snarl& snarl) {

        bool ret;
        bool is_reachable = snarl.start_end_reachable();
        bool is_DAG = snarl.directed_acyclic_net_graph();
        const Visit& start = snarl.start();
        int64_t node_id = start.node_id();
        bool has_node_start = path_position_handle_graph->has_node(node_id);
        bool has_node_end = path_position_handle_graph->has_node(snarl.end().node_id());

        // we need this to be true all the way down to use the RepresentativeTraversalFinder on our snarl.
        if(is_reachable && is_DAG  && has_node_start && has_node_end ){
            ret = true;
        }else{
            ret = false;
        }
        
        if (ret == true) {
            const vector<const Snarl*>& children = snarl_manager.children_of(&snarl);
            for (int i = 0; i < children.size() && ret; ++i) {
                ret = is_traversable(*children[i]);
            }
        }
        return ret;
    }


    pair<size_t, bool> MCMCCaller::get_ref_position(const Snarl& snarl, const string& ref_path_name) const {
    path_handle_t path_handle = path_position_handle_graph->get_path_handle(ref_path_name);

    handle_t start_handle = path_position_handle_graph->get_handle(snarl.start().node_id(), snarl.start().backward());
    map<size_t, step_handle_t> start_steps;
    path_position_handle_graph->for_each_step_on_handle(start_handle, [&](step_handle_t step) {
            if (path_position_handle_graph->get_path_handle_of_step(step) == path_handle) {
                start_steps[path_position_handle_graph->get_position_of_step(step)] = step;
            }
        });

    handle_t end_handle = path_position_handle_graph->get_handle(snarl.end().node_id(), snarl.end().backward());
    map<size_t, step_handle_t> end_steps;
    path_position_handle_graph->for_each_step_on_handle(end_handle, [&](step_handle_t step) {
            if (path_position_handle_graph->get_path_handle_of_step(step) == path_handle) {
                end_steps[path_position_handle_graph->get_position_of_step(step)] = step;
            }
        });

    assert(start_steps.size() > 0 && end_steps.size() > 0);
    step_handle_t start_step = start_steps.begin()->second;
    step_handle_t end_step = end_steps.begin()->second;
    bool scan_backward = path_position_handle_graph->get_is_reverse(path_position_handle_graph->get_handle_of_step(start_step));

    // if we're on a cycle, we keep our start step and find the end step by scanning the path
    if (start_steps.size() > 1 || end_steps.size() > 1) {
        bool found_end = false;
        if (scan_backward) {
            for (step_handle_t cur_step = start_step; path_position_handle_graph->has_previous_step(end_step) && !found_end;
                 cur_step = path_position_handle_graph->get_previous_step(cur_step)) {
                if (path_position_handle_graph->get_handle_of_step(cur_step) == end_handle) {
                    end_step = cur_step;
                    found_end = true;
                }
            }
            assert(found_end);
        } else {
            for (step_handle_t cur_step = start_step; path_position_handle_graph->has_next_step(end_step) && !found_end;
                 cur_step = path_position_handle_graph->get_next_step(cur_step)) {
                if (path_position_handle_graph->get_handle_of_step(cur_step) == end_handle) {
                    end_step = cur_step;
                    found_end = true;
                }
            }
            assert(found_end);
        }
    }
    
    size_t start_position = start_steps.begin()->first;
    size_t end_position = end_step == end_steps.begin()->second ? end_steps.begin()->first : path_position_handle_graph->get_position_of_step(end_step);
    bool backward = end_position < start_position;

    return make_pair(backward ? end_position : start_position, backward);
    }



    void MCMCCaller::flatten_common_allele_ends(vcflib::Variant& variant, bool backward) const {
        if (variant.alt.size() == 0) {
            return;
        }
        size_t min_len = variant.alleles[0].length();
        for (int i = 1; i < variant.alleles.size(); ++i) {
            min_len = std::min(min_len, variant.alleles[i].length());
        }
        // want to leave at least one in the reference position
        if (min_len > 0) {
            --min_len;
        }

        bool match = true;
        int shared_prefix_len = 0;
        for (int i = 0; i < min_len && match; ++i) {
            char c1 = std::toupper(variant.alleles[0][!backward ? i : variant.alleles[0].length() - 1 - i]);
            for (int j = 1; j < variant.alleles.size() && match; ++j) {
                char c2 = std::toupper(variant.alleles[j][!backward ? i : variant.alleles[j].length() - 1 - i]);
                match = c1 == c2;
            }
            if (match) {
                ++shared_prefix_len;
            }
        }

        if (!backward) {
            variant.position += shared_prefix_len;
        }
        for (int i = 0; i < variant.alleles.size(); ++i) {
            if (!backward) {
                variant.alleles[i] = variant.alleles[i].substr(shared_prefix_len);
            } else {
                variant.alleles[i] = variant.alleles[i].substr(0, variant.alleles[i].length() - shared_prefix_len);
            }
            if (i == 0) {
                variant.ref = variant.alleles[i];
            } else {
                variant.alt[i - 1] = variant.alleles[i];
            }
        }
    }

}  

// #undef debug
