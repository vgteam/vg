#include "graph_caller.hpp"
#include "algorithms/expand_context.hpp"

namespace vg {

GraphCaller::GraphCaller(SnarlCaller& snarl_caller,
                         SnarlManager& snarl_manager,
                         ostream& out_stream) :
    snarl_caller(snarl_caller), snarl_manager(snarl_manager), out_stream(out_stream) {
}

GraphCaller::~GraphCaller() {
}

void GraphCaller::call_top_level_snarls(bool recurse_on_fail) {

    // Used to recurse on children of parents that can't be called
    vector<const Snarl*> snarl_queue;

    // Run the snarl caller on a snarl, and queue up the children if it fails
    auto process_snarl = [&](const Snarl* snarl) {
        bool was_called = call_snarl(*snarl);
        if (!was_called && recurse_on_fail) {
            const vector<const Snarl*>& children = snarl_manager.children_of(snarl);
#pragma omp critical (snarl_queue)
            {
                snarl_queue.insert(snarl_queue.end(), children.begin(), children.end());
            }
        }
    };

    // Start with the top level snarls
    snarl_manager.for_each_top_level_snarl_parallel(process_snarl);

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

string GraphCaller::vcf_header(const PathHandleGraph& graph, const vector<string>& contigs) const {
    stringstream ss;
    ss << "##fileformat=VCFv4.2" << endl;
    for (auto contig : contigs) {
        size_t length = 0;
        for (handle_t handle : graph.scan_path(graph.get_path_handle(contig))) {
            length += graph.get_length(handle);
        }
        ss << "##contig=<ID=" << contig << ",length=" << length << ">" << endl;
    }
    return ss.str();
}

VCFGenotyper::VCFGenotyper(const PathHandleGraph& graph,
                           SnarlCaller& snarl_caller,
                           SnarlManager& snarl_manager,
                           vcflib::VariantCallFile& variant_file,
                           const string& sample_name,
                           const vector<string>& ref_paths,
                           FastaReference* ref_fasta,
                           FastaReference* ins_fasta,
                           ostream& out_stream) :
    GraphCaller(snarl_caller, snarl_manager, out_stream),
    graph(graph),
    input_vcf(variant_file),
    sample_name(sample_name),
    traversal_finder(graph, snarl_manager, variant_file, ref_paths, ref_fasta, ins_fasta) {
    
}

VCFGenotyper::~VCFGenotyper() {

}

bool VCFGenotyper::call_snarl(const Snarl& snarl) {

    // get our traversals out of the finder
    vector<pair<SnarlTraversal, vector<int>>> alleles;
    vector<vcflib::Variant*> variants;
    std::tie(alleles, variants) = traversal_finder.find_allele_traversals(snarl);

    if (!alleles.empty()) {

        // hmm, maybe find a way not to copy?
        vector<SnarlTraversal> travs;
        travs.reserve(alleles.size());
        for (const auto& ta : alleles) {
            travs.push_back(ta.first);
        }

        // find the reference traversal
        // todo: is it the reference always first?
        int ref_trav_idx = -1;
        for (int i = 0; i < alleles.size() && ref_trav_idx < 0; ++i) {
            if (std::all_of(alleles[i].second.begin(), alleles[i].second.end(), [](int x) {return x == 0;})) {
                ref_trav_idx = i;
            }
        }

        // use our support caller to choose our genotype (int traversal coordinates)
        vector<int> trav_genotype = snarl_caller.genotype(snarl, travs, ref_trav_idx, 2);

        assert(trav_genotype.empty() || trav_genotype.size() == 2);

        // map our genotype back to the vcf
        for (int i = 0; i < variants.size(); ++i) {
            vector<int> vcf_alleles;
            string vcf_genotype;
            vector<SnarlTraversal> vcf_traversals;
            
            if (trav_genotype.empty()) {
                vcf_genotype = "./.";
            } else {
                // map our traversal genotype to a vcf variant genotype
                // using the information out of the traversal finder
                for (int j = 0; j < trav_genotype.size(); ++j) {
                    int trav_allele = trav_genotype[j];
                    int vcf_allele = alleles[trav_allele].second[i];
                    vcf_genotype += std::to_string(vcf_allele);
                    if (j < trav_genotype.size() - 1) {
                        vcf_genotype += "/";
                    }
                    vcf_alleles.push_back(vcf_allele);
                    vcf_traversals.push_back(travs[trav_allele]);
                }
            }
            // create an output variant from the input one
            vcflib::Variant out_variant;
            out_variant.sequenceName = variants[i]->sequenceName;
            out_variant.quality = 23;
            out_variant.position = variants[i]->position;
            out_variant.id = variants[i]->id;
            out_variant.ref = variants[i]->ref;
            out_variant.alt = variants[i]->alt;
            out_variant.alleles = variants[i]->alleles;
            out_variant.setVariantCallFile(output_vcf);
            out_variant.filter = "PASS";
            out_variant.updateAlleleIndexes();

            // add the genotype
            out_variant.format.push_back("GT");
            auto& genotype_vector = out_variant.samples[sample_name]["GT"];
            genotype_vector.push_back(vcf_genotype);

            // add some info
            snarl_caller.update_vcf_info(snarl, vcf_traversals, vcf_alleles, sample_name, out_variant);

            // print the variant
            out_stream << out_variant << endl;
        }
        return true;
    }
    
    return false;

}

string VCFGenotyper::vcf_header(const PathHandleGraph& graph, const vector<string>& ref_paths) const {
    string header = GraphCaller::vcf_header(graph, ref_paths);
    header += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    header += "##FILTER=<ID=PASS,Description=\"All filters passed\">";
    header += "##SAMPLE=<ID=" + sample_name + ">\n";
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_name;
    assert(output_vcf.openForOutput(header));
    header += "\n";
    return header;
}


LegacyCaller::LegacyCaller(const PathHandleGraph& graph,
                           SupportBasedSnarlCaller& snarl_caller,
                           SnarlManager& snarl_manager,
                           const string& sample_name,
                           const vector<string>& ref_paths) :
    GraphCaller(snarl_caller, snarl_manager, out_stream),
    graph(graph),
    sample_name(sample_name),
    ref_paths(ref_paths) {
    
    is_vg = dynamic_cast<const VG*>(&graph) != nullptr;
    if (is_vg) {
        // our graph is in vg format.  we index the paths and make a traversal finder just
        // like in the old call code
        for (auto ref_path : ref_paths) {
            path_indexes.push_back(new PathIndex(graph, ref_path));
        }
        // map snarl to the first reference path that spans it
        function<PathIndex*(const Snarl&)> get_path_index = [&](const Snarl& site) -> PathIndex* {
            return find_index(site, path_indexes).second;
        };
        // initialize our traversal finder
        traversal_finder = new RepresentativeTraversalFinder(graph, snarl_manager,
                                                             max_search_depth,
                                                             max_search_width,
                                                             max_bubble_paths,
                                                             min_total_support_for_call,
                                                             min_total_support_for_call,
                                                             get_path_index,
                                                             [&](id_t id) { return snarl_caller.get_min_node_support(id);},
                                                             [&](edge_t edge) { return snarl_caller.get_edge_support(edge);});

    } else {
        // our graph is not in vg format.  we will make graphs for each site as needed and work with those
        traversal_finder = nullptr;
    }
}

LegacyCaller::~LegacyCaller() {
    delete traversal_finder;
    for (PathIndex* path_index : path_indexes) {
        delete path_index;
    }
}

bool LegacyCaller::call_snarl(const Snarl& snarl) {

    // if we can't handle the snarl, then the GraphCaller framework will recurse on its children
    if (!is_traversable(snarl)) {
        return false;
    }
           
    RepresentativeTraversalFinder* rep_trav_finder;
    vector<PathIndex*> site_path_indexes;
    function<PathIndex*(const Snarl&)> get_path_index;
    VG vg_graph;
    SupportBasedSnarlCaller& support_caller = dynamic_cast<SupportBasedSnarlCaller&>(snarl_caller);
    
    if (is_vg) {
        // our graph is in VG format, so we've sorted this out in the constructor
        rep_trav_finder = traversal_finder;
        get_path_index = [&](const Snarl& site) {
            return find_index(site, path_indexes).second;
        };
        
    } else {
        // our graph isn't in VG format.  we are using a (hopefully temporary) workaround
        // of converting the subgraph into VG.
        pair<unordered_set<id_t>, unordered_set<edge_t> > contents = snarl_manager.deep_contents(&snarl, graph, true);
        for (auto node_id : contents.first) {
            vg_graph.create_handle(graph.get_sequence(graph.get_handle(node_id)), node_id);
        }
        for (auto edge : contents.second) {
          vg_graph.create_edge(vg_graph.get_handle(graph.get_id(edge.first), vg_graph.get_is_reverse(edge.first)),
                               vg_graph.get_handle(graph.get_id(edge.second), vg_graph.get_is_reverse(edge.second)));
        }
        // add the paths to the subgraph
        algorithms::expand_context_with_paths(&graph, &vg_graph, 1);
        // and index them
        for (auto& ref_path : ref_paths) {
            site_path_indexes.push_back(new PathIndex(vg_graph, ref_path));
        }
        get_path_index = [&](const Snarl& site) -> PathIndex* {
            return find_index(site, site_path_indexes).second;
        };
        rep_trav_finder = new RepresentativeTraversalFinder(vg_graph, snarl_manager,
                                                            max_search_depth,
                                                            max_search_width,
                                                            max_bubble_paths,
                                                            min_total_support_for_call,
                                                            min_total_support_for_call,
                                                            get_path_index,
                                                            [&](id_t id) { return support_caller.get_min_node_support(id);},
                                                            // note: because our traversal finder and support caller have
                                                            // different graphs, they can't share edge handles
                                                            [&](edge_t edge) { return support_caller.get_edge_support(
                                                                    vg_graph.get_id(edge.first), vg_graph.get_is_reverse(edge.first),
                                                                    vg_graph.get_id(edge.second), vg_graph.get_is_reverse(edge.second));});
                                                            
    }

    PathIndex* path_index = get_path_index(snarl);
    if (path_index != nullptr) {
        // recursively genotype the site beginning here at the top level snarl
        vector<SnarlTraversal> called_traversals = top_down_genotype(snarl, *rep_trav_finder, 2);
    
        // emit our vcf variant
        if (!called_traversals.empty()) {
            string path_name = find_index(snarl, is_vg ? path_indexes : site_path_indexes).first;
            emit_variant(snarl, called_traversals, path_name);
        }
    }        
    if (!is_vg) {
        // delete the temporary vg subgraph and traversal finder we created for this snarl
        delete rep_trav_finder;
        for (PathIndex* path_index : site_path_indexes) {
            delete path_index;
        }
    }

    return true;
}

string LegacyCaller::vcf_header(const PathHandleGraph& graph, const vector<string>& ref_paths) const {
    string header = GraphCaller::vcf_header(graph, ref_paths);
    header += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    header += "##FILTER=<ID=PASS,Description=\"All filters passed\">";
    header += "##SAMPLE=<ID=" + sample_name + ">\n";
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_name;
    assert(output_vcf.openForOutput(header));
    header += "\n";
    return header;
}

vector<SnarlTraversal> LegacyCaller::top_down_genotype(const Snarl& snarl, TraversalFinder& trav_finder, int ploidy) const {
    vector<SnarlTraversal> genotype;

    // get the traversals through the site
    vector<SnarlTraversal> traversals = trav_finder.find_traversals(snarl);

    cerr << "our representative traversals" << endl;
    for (auto xx : traversals) {
        cerr << "  " << pb2json(xx) << endl;
    }

    // use our support caller to choose our genotype
    vector<int> trav_genotype = snarl_caller.genotype(snarl, traversals, 0, ploidy);

    assert(trav_genotype.empty() || trav_genotype.size() == ploidy);
    // todo: handle non-call
    
    genotype.resize(ploidy);

    // do we have two paths going through a given traversal?  This is handled
    // as a special case below
    bool hom = trav_genotype.size() == 2 && trav_genotype[0] == trav_genotype[1];
    
    for (int i = 0; i < trav_genotype.size() && (!hom || i < 1); ++i) {
        int allele = trav_genotype[i];
        const SnarlTraversal& traversal = traversals[allele];
        for (int j = 0; j < traversal.visit_size(); ++j) {
            if (traversal.visit(j).node_id() > 0) {
                *genotype[i].add_visit() = traversal.visit(j);
                if (hom && i == 0) {
                    *genotype[1].add_visit() = traversal.visit(j);
                }
            } else {
                // recursively determine the traversal    
                vector<SnarlTraversal> child_genotype = top_down_genotype(*snarl_manager.into_which_snarl(traversal.visit(j)),
                                                                          trav_finder, hom ? 2: 1);
                bool back_to_back = j > 0 && traversal.visit(j - 1).node_id() == 0;
                assert(!back_to_back || traversal.visit(j - 1).snarl().end() == traversal.visit(j).snarl().start());

                // todo: handle empty case
                for (int k = back_to_back ? 1 : 0; k < child_genotype[0].visit_size(); ++k) {
                    *genotype[i].add_visit() = child_genotype[0].visit(k);
                }
                if (hom) {
                    assert(child_genotype.size() == 2 && i == 0);
                    for (int k = back_to_back ? 1 : 0; k < child_genotype[1].visit_size(); ++k) {
                        *genotype[1].add_visit() = child_genotype[1].visit(k);
                    }
                }
            }
        }
    }

    return genotype;
}

void LegacyCaller::emit_variant(const Snarl& snarl, const vector<SnarlTraversal>& called_traversals, const string& ref_path_name) const {

    // converet the traversals to strings and map them to a genotype
    map<string, int> allele_to_gt;
    for (int i = 0; i < called_traversals.size(); ++i) {
        const SnarlTraversal& trav = called_traversals[i];
        string allele;
        for (int i = 0; i < trav.visit_size(); ++i) {
            allele += graph.get_sequence(graph.get_handle(trav.visit(i).node_id(), trav.visit(i).backward()));
        }
        if (!allele_to_gt.count(allele)) {
            allele_to_gt[allele] = allele_to_gt.size();
        }
    }

    vcflib::Variant out_variant;
    out_variant.sequenceName = ref_path_name;
    out_variant.quality = 23;
    out_variant.position = 0;
    out_variant.id = ".";
    for (auto& allele_gt : allele_to_gt) {
        if (allele_gt.second == 0) {
            out_variant.ref = allele_gt.first;
        } else {
            out_variant.alt.push_back(allele_gt.first);
        }
        out_variant.alleles.push_back(allele_gt.first);
    }
    out_variant.setVariantCallFile(output_vcf);
    out_variant.filter = "PASS";
    out_variant.updateAlleleIndexes();

    // add the genotype
    out_variant.format.push_back("GT");
    auto& genotype_vector = out_variant.samples[sample_name]["GT"];
    genotype_vector.push_back("0/0");

    // add some info
    snarl_caller.update_vcf_info(snarl, called_traversals, {0,0}, sample_name, out_variant);

    cout << out_variant << endl;
}

bool LegacyCaller::is_traversable(const Snarl& snarl) {
    // we need this to be true all the way down to use the RepresentativeTraversalFinder on our snarl.
    bool ret = snarl.start_end_reachable() && snarl.directed_acyclic_net_graph();
    if (ret == true) {
        const vector<const Snarl*>& children = snarl_manager.children_of(&snarl);
        for (int i = 0; i < children.size() && ret; ++i) {
            ret = is_traversable(*children[i]);
        }
    }
    return ret;
}

pair<string, PathIndex*> LegacyCaller::find_index(const Snarl& snarl, const vector<PathIndex*> path_indexes) const {
    assert(path_indexes.size() == ref_paths.size());
    for (int i = 0; i < path_indexes.size(); ++i) {
        PathIndex* path_index = path_indexes[i];
        if (path_index->by_id.count(snarl.start().node_id()) &&
            path_index->by_id.count(snarl.end().node_id())) {
            // This path threads through this site
            return make_pair(ref_paths[i], path_index);
        }
    }
    return make_pair("", nullptr);
}

}

