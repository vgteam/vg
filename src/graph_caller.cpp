#include "graph_caller.hpp"
#include "algorithms/expand_context.hpp"
#include "annotation.hpp"

//#define debug

namespace vg {

GraphCaller::GraphCaller(SnarlCaller& snarl_caller,
                         SnarlManager& snarl_manager) :
    snarl_caller(snarl_caller), snarl_manager(snarl_manager) {
}

GraphCaller::~GraphCaller() {
}

void GraphCaller::call_top_level_snarls(const HandleGraph& graph, RecurseType recurse_type) {

    // Used to recurse on children of parents that can't be called
    size_t thread_count = get_thread_count();
    vector<vector<const Snarl*>> snarl_queue(thread_count);

    // Run the snarl caller on a snarl, and queue up the children if it fails
    auto process_snarl = [&](const Snarl* snarl) {

        if (!snarl_manager.is_trivial(snarl, graph)) {

#ifdef debug
            cerr << "GraphCaller running call_snarl on " << pb2json(*snarl) << endl;
#endif

            bool was_called = call_snarl(*snarl);
            if (recurse_type == RecurseAlways || (!was_called && recurse_type == RecurseOnFail)) {
                const vector<const Snarl*>& children = snarl_manager.children_of(snarl);
                vector<const Snarl*>& thread_queue = snarl_queue[omp_get_thread_num()];
                thread_queue.insert(thread_queue.end(), children.begin(), children.end());
            }
        }
    };

    // Start with the top level snarls
    snarl_manager.for_each_top_level_snarl_parallel(process_snarl);

    // Then recurse on any children the snarl caller failed to handle
    while (!std::all_of(snarl_queue.begin(), snarl_queue.end(),
                        [](const vector<const Snarl*>& snarl_vec) {return snarl_vec.empty();})) {
        vector<const Snarl*> cur_queue;
        for (vector<const Snarl*>& thread_queue : snarl_queue) {
            cur_queue.reserve(cur_queue.size() + thread_queue.size());
            std::move(thread_queue.begin(), thread_queue.end(), std::back_inserter(cur_queue));
            thread_queue.clear();
        }

#pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < cur_queue.size(); ++i) {
            process_snarl(cur_queue[i]);
        }
    }
  
}

static void flip_snarl(Snarl& snarl) {
    Visit v = snarl.start();
    *snarl.mutable_start() = reverse(snarl.end());
    *snarl.mutable_end() = reverse(v);
}

void GraphCaller::call_top_level_chains(const HandleGraph& graph, size_t max_edges, size_t max_trivial, RecurseType recurse_type) {
    // Used to recurse on children of parents that can't be called
    size_t thread_count = get_thread_count();
    vector<vector<Chain>> chain_queue(thread_count);

    // Run the snarl caller on a chain. queue up the children if it fails
    auto process_chain = [&](const Chain* chain) {

#ifdef debug
        cerr << "calling top level chain ";
        for (const auto& i : *chain) {
            cerr << pb2json(*i.first) << "," << i.second << ",";
        }
        cerr << endl;
#endif
        // Break up the chain
        vector<Chain> chain_pieces = break_chain(graph, *chain, max_edges, max_trivial);

        for (Chain& chain_piece : chain_pieces) {
            // Make a fake snarl spanning the chain
            // It is important to remember that along with not actually being a snarl,
            // it's not managed by the snarl manager so functions looking into its nesting
            // structure will not work
            Snarl fake_snarl;
            *fake_snarl.mutable_start() = chain_piece.front().second == true ? reverse(chain_piece.front().first->end()) :
                chain_piece.front().first->start();
            *fake_snarl.mutable_end() = chain_piece.back().second == true ? reverse(chain_piece.back().first->start()) :
                chain_piece.back().first->end();

#ifdef debug
            cerr << "calling fake snarl " << pb2json(fake_snarl) << endl;
#endif
            
            bool was_called = call_snarl(fake_snarl);
            if (recurse_type == RecurseAlways || (!was_called && recurse_type == RecurseOnFail)) {
                vector<Chain>& thread_queue = chain_queue[omp_get_thread_num()];                
                for (pair<const Snarl*, bool> chain_link : chain_piece) {
                    const deque<Chain>& child_chains = snarl_manager.chains_of(chain_link.first);
                    thread_queue.insert(thread_queue.end(), child_chains.begin(), child_chains.end());
                }
            }
        }
    };

    // Start with the top level snarls
    snarl_manager.for_each_top_level_chain_parallel(process_chain);

    // Then recurse on any children the snarl caller failed to handle
    while (!std::all_of(chain_queue.begin(), chain_queue.end(),
                        [](const vector<Chain>& chain_vec) {return chain_vec.empty();})) {
        vector<Chain> cur_queue;
        for (vector<Chain>& thread_queue : chain_queue) {
            cur_queue.reserve(cur_queue.size() + thread_queue.size());
            std::move(thread_queue.begin(), thread_queue.end(), std::back_inserter(cur_queue));
            thread_queue.clear();
        }

#pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < cur_queue.size(); ++i) {
            process_chain(&cur_queue[i]);
        }
    
    }
}

vector<Chain> GraphCaller::break_chain(const HandleGraph& graph, const Chain& chain, size_t max_edges, size_t max_trivial) {
    
    vector<Chain> chain_frags;

    // keep track of the current fragment and add it to chain_frags as soon as it gets too big
    Chain frag;
    size_t frag_edge_count = 0;
    size_t frag_triv_count = 0;
    
    for (const pair<const Snarl*, bool>& link : chain) {
        // todo: we're getting the contents here as well as within the caller.
        auto contents = snarl_manager.deep_contents(link.first, graph, false);

        // todo: use annotation from snarl itself?
        bool trivial = contents.second.empty();

        if ((trivial && frag_triv_count > max_trivial) ||
            (contents.second.size() + frag_edge_count > max_edges)) {
            // adding anything more to the chain would make it too long, so we
            // add it to the output and clear the current fragment
            if (!frag.empty() && frag_triv_count < frag.size()) {
                chain_frags.push_back(frag);
            }
            frag.clear();
            frag_edge_count = 0;
            frag_triv_count = 0;
        }

        if (!trivial || (frag_triv_count < max_trivial)) {
            // we start a new fragment or add to an existing fragment
            frag.push_back(link);
            frag_edge_count += contents.second.size();
            if (trivial) {
                ++frag_triv_count;
            }
        }
    }

    // and the last one
    if (!frag.empty()) {
        chain_frags.push_back(frag);
    }

    return chain_frags;
}
    
VCFOutputCaller::VCFOutputCaller(const string& sample_name) : sample_name(sample_name), translation(nullptr), include_nested(false)
{
    output_variants.resize(get_thread_count());
}

VCFOutputCaller::~VCFOutputCaller() {
}

string VCFOutputCaller::vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                                   const vector<size_t>& contig_length_overrides) const {
    stringstream ss;
    ss << "##fileformat=VCFv4.2" << endl;    
    for (int i = 0; i < contigs.size(); ++i) {
        const string& contig = contigs[i];
        size_t length;
        if (i < contig_length_overrides.size()) {
            // length override provided
            length = contig_length_overrides[i];
        } else {
            length = 0;
            for (handle_t handle : graph.scan_path(graph.get_path_handle(contig))) {
                length += graph.get_length(handle);
            }
        }
        ss << "##contig=<ID=" << contig << ",length=" << length << ">" << endl;
    }
    if (include_nested) {
        ss << "##INFO=<ID=LV,Number=1,Type=Integer,Description=\"Level in the snarl tree (0=top level)\">" << endl;
        ss << "##INFO=<ID=PS,Number=1,Type=String,Description=\"ID of variant corresponding to parent snarl\">" << endl;
    }
    ss << "##INFO=<ID=AT,Number=R,Type=String,Description=\"Allele Traversal as path in graph\">" << endl;
    return ss.str();
}

void VCFOutputCaller::add_variant(vcflib::Variant& var) const {
    var.setVariantCallFile(output_vcf);
    stringstream ss;
    ss << var;
    string dest;
    zstdutil::CompressString(ss.str(), dest);
    // the Variant object is too big to keep in memory when there are many genotypes, so we
    // store it in a zstd-compressed string
    output_variants[omp_get_thread_num()].push_back(make_pair(make_pair(var.sequenceName, var.position), dest));
}

void VCFOutputCaller::write_variants(ostream& out_stream, const SnarlManager* snarl_manager) {
    assert(include_nested == false || snarl_manager != nullptr);
    if (include_nested) {
        update_nesting_info_tags(snarl_manager);
    }
    vector<pair<pair<string, size_t>, string>> all_variants;
    for (const auto& buf : output_variants) {
        all_variants.reserve(all_variants.size() + buf.size());
        std::move(buf.begin(), buf.end(), std::back_inserter(all_variants));
    }
    std::sort(all_variants.begin(), all_variants.end(), [](const pair<pair<string, size_t>, string>& v1,
                                                           const pair<pair<string, size_t>, string>& v2) {
            return v1.first.first < v2.first.first || (v1.first.first == v2.first.first && v1.first.second < v2.first.second);
        });
    for (auto v : all_variants) {
        string dest;
        zstdutil::DecompressString(v.second, dest);
        out_stream << dest << endl;
    }
}

static int countAlts(vcflib::Variant& var, int alleleIndex) {
    int alts = 0;
    for (map<string, map<string, vector<string> > >::iterator s = var.samples.begin(); s != var.samples.end(); ++s) {
        map<string, vector<string> >& sample = s->second;
        map<string, vector<string> >::iterator gt = sample.find("GT");
        if (gt != sample.end()) {
            map<int, int> genotype = vcflib::decomposeGenotype(gt->second.front());
            for (map<int, int>::iterator g = genotype.begin(); g != genotype.end(); ++g) {
                if (g->first == alleleIndex) {
                    alts += g->second;
                }
            }
        }
    }
    return alts;
}

static int countAlleles(vcflib::Variant& var) {
    int alleles = 0;
    for (map<string, map<string, vector<string> > >::iterator s = var.samples.begin(); s != var.samples.end(); ++s) {
        map<string, vector<string> >& sample = s->second;
        map<string, vector<string> >::iterator gt = sample.find("GT");
        if (gt != sample.end()) {
            map<int, int> genotype = vcflib::decomposeGenotype(gt->second.front());
            for (map<int, int>::iterator g = genotype.begin(); g != genotype.end(); ++g) {
		if (g->first != vcflib::NULL_ALLELE) {
		    alleles += g->second;
		}
            }
        }
    }
    return alleles;
}

// this isn't from vcflib, but seems to make more sense than just returning the number of samples in the file again and again
static int countSamplesWithData(vcflib::Variant& var) {
    int samples_with_data = 0;
    for (map<string, map<string, vector<string> > >::iterator s = var.samples.begin(); s != var.samples.end(); ++s) {
        map<string, vector<string> >& sample = s->second;
        map<string, vector<string> >::iterator gt = sample.find("GT");
        bool has_data = false;
        if (gt != sample.end()) {
            map<int, int> genotype = vcflib::decomposeGenotype(gt->second.front());
            for (map<int, int>::iterator g = genotype.begin(); g != genotype.end(); ++g) {
		if (g->first != vcflib::NULL_ALLELE) {
                    has_data = true;
                    break;
		}
            }
        }
        if (has_data) {
            ++samples_with_data;
        }
    }
    return samples_with_data;
}

void VCFOutputCaller::vcf_fixup(vcflib::Variant& var) const {
    // copied from https://github.com/vgteam/vcflib/blob/master/src/vcffixup.cpp
    
    stringstream ns;
    ns << countSamplesWithData(var);
    var.info["NS"].clear();
    var.info["NS"].push_back(ns.str());

    var.info["AC"].clear();
    var.info["AF"].clear();
    var.info["AN"].clear();

    int allelecount = countAlleles(var);
    stringstream an;
    an << allelecount;
    var.info["AN"].push_back(an.str());

    for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
        string& allele = *a;
        int altcount = countAlts(var, var.getAltAlleleIndex(allele) + 1);
        stringstream ac;
        ac << altcount;
        var.info["AC"].push_back(ac.str());
        stringstream af;
        double faf = (double) altcount / (double) allelecount;
        if(faf != faf) faf = 0;
        af << faf;
        var.info["AF"].push_back(af.str());
    }
}

void VCFOutputCaller::set_translation(const unordered_map<nid_t, pair<nid_t, size_t>>* translation) {
    translation = translation;
}

void VCFOutputCaller::set_nested(bool nested) {
    include_nested = nested;
}

void VCFOutputCaller::add_allele_path_to_info(vcflib::Variant& v, int allele, const SnarlTraversal& trav,
                                              bool reversed, bool one_based) const {
    auto& trav_info = v.info["AT"];
    assert(allele < trav_info.size());

    vector<int> nodes;
    nodes.reserve(trav.visit_size());
    const Visit* prev_visit = nullptr;
    unordered_map<nid_t, pair<nid_t, size_t>>::const_iterator prev_trans;
    
    for (size_t i = 0; i < trav.visit_size(); ++i) {
        size_t j = !reversed ? i : trav.visit_size() - 1 - i;
        const Visit& visit = trav.visit(j);
        nid_t node_id = visit.node_id();
        bool skip = false;
        // todo: check one_based? (we kind of ignore that when writing the snarl name, so maybe not pertienent)
        if (translation) {
            auto i = translation->find(node_id);
            if (i == translation->end()) {
                throw runtime_error("Error [vg deconstruct]: Unable to find node " + std::to_string(node_id) + " in translation file");
            }
            if (prev_visit) {
                nid_t prev_node_id = prev_visit->node_id();
                if (prev_trans->second.first == i->second.first && node_id != prev_node_id) {
                    // here is a case where we have two consecutive nodes that map back to
                    // the same source node.
                    // todo: sanity check! (could verify if translation node properly covered)
                    skip = true;
                }
            }
            node_id = i->second.first;
            prev_trans = i;
        }

        if (!skip) {
            bool vrev = visit.backward() != reversed;
            trav_info[allele] += (vrev ? "<" : ">");
            trav_info[allele] += std::to_string(node_id);
        }
        prev_visit = &visit;
    }
}

string VCFOutputCaller::trav_string(const HandleGraph& graph, const SnarlTraversal& trav) const {
    string seq;
    for (int i = 0; i < trav.visit_size(); ++i) {
        const Visit& visit = trav.visit(i);
        if (visit.node_id() > 0) {
            seq += graph.get_sequence(graph.get_handle(visit.node_id(), visit.backward()));
        } else {
            seq += print_snarl(visit.snarl(), true);
        }
    }
    return seq;    
}

void VCFOutputCaller::emit_variant(const PathPositionHandleGraph& graph, SnarlCaller& snarl_caller,
                                   const Snarl& snarl, const vector<SnarlTraversal>& called_traversals,
                                   const vector<int>& genotype, int ref_trav_idx, const unique_ptr<SnarlCaller::CallInfo>& call_info,
                                   const string& ref_path_name, int ref_offset, bool genotype_snarls, int ploidy,
                                   function<string(const vector<SnarlTraversal>&, const vector<int>&, int, int, int)> trav_to_string) {
    
#ifdef debug
    cerr << "emitting variant for " << pb2json(snarl) << endl;
    for (int i = 0; i < called_traversals.size(); ++i) {
        if (i == ref_trav_idx) {
            cerr << "*";
        }
        cerr << "ct[" << i << "]=" << pb2json(called_traversals[i]) << endl;
    }
    for (int i = 0; i < genotype.size(); ++i) {
        cerr << "gt[" << i << "]=" << genotype[i] << endl;
    }
#endif

    if (trav_to_string == nullptr) {
        trav_to_string = [&](const vector<SnarlTraversal>& travs, const vector<int>& travs_genotype, int trav_allele, int genotype_allele, int ref_trav_idx) {
            return trav_string(graph, travs[trav_allele]);    
        };
    }

    vcflib::Variant out_variant;

    vector<SnarlTraversal> site_traversals = {called_traversals[ref_trav_idx]};
    vector<int> site_genotype;
    auto ref_gt_it = std::find(genotype.begin(), genotype.end(), ref_trav_idx);
    out_variant.ref = trav_to_string(called_traversals, genotype, ref_trav_idx,
                                     ref_gt_it != genotype.end() ? ref_gt_it - genotype.begin() : 0,
                                     ref_trav_idx);
    
    // deduplicate alleles and compute the site traversals and genotype
    map<string, int> allele_to_gt;    
    allele_to_gt[out_variant.ref] = 0;    
    for (int i = 0; i < genotype.size(); ++i) {
        if (genotype[i] == ref_trav_idx) {
            site_genotype.push_back(0);
        } else {
            string allele_string = trav_to_string(called_traversals, genotype, genotype[i], i, ref_trav_idx);
            if (allele_to_gt.count(allele_string)) {
                site_genotype.push_back(allele_to_gt[allele_string]);
            } else {
                site_traversals.push_back(called_traversals[genotype[i]]);
                site_genotype.push_back(allele_to_gt.size());
                allele_to_gt[allele_string] = site_genotype.back();
            }
        }
    }

    // add on fixed number of uncalled traversals if we're making a ref-call
    // with genotype_snarls set to true
    if (genotype_snarls && site_traversals.size() <= 1) {
        // note: we're adding all the strings here and sorting to make this deterministic
        // at the cost of speed
        map<string, const SnarlTraversal*> allele_map;
        for (int i = 0; i < called_traversals.size(); ++i) {
            // todo: verify index below.  it's for uncalled traversals so not important tho
            string allele_string = trav_to_string(called_traversals, genotype, i, max(0, (int)genotype.size() - 1), ref_trav_idx);
            if (!allele_map.count(allele_string)) {
                allele_map[allele_string] = &called_traversals[i];
            }
        }
        // pick out the first "max_uncalled_alleles" traversals to add
        int i = 0;
        for (auto ai = allele_map.begin(); i < max_uncalled_alleles && ai != allele_map.end(); ++i, ++ai) {
            if (!allele_to_gt.count(ai->first)) {
                allele_to_gt[ai->first] = allele_to_gt.size();
                site_traversals.push_back(*ai->second);
            }
        }
    }

    out_variant.alt.resize(allele_to_gt.size() - 1);
    out_variant.alleles.resize(allele_to_gt.size());
    
    // init the traversal info
    out_variant.info["AT"].resize(allele_to_gt.size());

    for (auto& allele_gt : allele_to_gt) {
#ifdef debug
        cerr << "allele " << allele_gt.first << " -> gt " << allele_gt.second << endl;
#endif
        if (allele_gt.second > 0) {
            out_variant.alt[allele_gt.second - 1] = allele_gt.first;
        }
        out_variant.alleles[allele_gt.second] = allele_gt.first;

        // update the traversal info
        add_allele_path_to_info(out_variant, allele_gt.second, site_traversals.at(allele_gt.second), false, false); 
    }

    // resolve subpath naming
    string basepath_name = ref_path_name;
    size_t basepath_offset = 0;
    auto subpath_info = Paths::parse_subpath_name(ref_path_name);
    if (get<0>(subpath_info)) {
        basepath_name = get<1>(subpath_info);
        basepath_offset = get<2>(subpath_info);
    }
    // fill out the rest of the variant
    out_variant.sequenceName = basepath_name;
    // +1 to convert to 1-based VCF
    out_variant.position = get<0>(get_ref_interval(graph, snarl, ref_path_name)) + ref_offset + 1 + basepath_offset;
    out_variant.id = print_snarl(snarl, false);
    out_variant.filter = "PASS";
    out_variant.updateAlleleIndexes();

    // add the genotype
    out_variant.format.push_back("GT");
    auto& genotype_vector = out_variant.samples[sample_name]["GT"];
    
    stringstream vcf_gt;
    if (!genotype.empty()) {
        for (int i = 0; i < site_genotype.size(); ++i) {
            vcf_gt << site_genotype[i];
            if (i != site_genotype.size() - 1) {
                vcf_gt << "/";
            }
        }
    } else {
        for (int i = 0; i < ploidy; ++i) {
            vcf_gt << ".";
            if (i != ploidy - 1) {
                vcf_gt << "/";
            }
        }
    }
                    
    genotype_vector.push_back(vcf_gt.str());

    // add some support info
    snarl_caller.update_vcf_info(snarl, site_traversals, site_genotype, call_info, sample_name, out_variant);

    // if genotype_snarls, then we only flatten up to the snarl endpoints
    // (this is when we are in genotyping mode and want consistent calls regardless of the sample)
    int64_t flatten_len_s = 0;
    int64_t flatten_len_e = 0;
    if (genotype_snarls) {
        flatten_len_s = graph.get_length(graph.get_handle(snarl.start().node_id()));
        assert(flatten_len_s >= 0);
        flatten_len_e = graph.get_length(graph.get_handle(snarl.end().node_id()));
    }
    // clean up the alleles to not have so man common prefixes
    flatten_common_allele_ends(out_variant, true, flatten_len_e);
    flatten_common_allele_ends(out_variant, false, flatten_len_s);
#ifdef debug
    for (int i = 0; i < site_traversals.size(); ++i) {
        cerr << " site trav[" << i << "]=" << pb2json(site_traversals[i]) << endl;
    }
    for (int i = 0; i < site_genotype.size(); ++i) {
        cerr << " site geno[" << i << "]=" << site_genotype[i] << endl;
    }
#endif

    if (genotype_snarls || !out_variant.alt.empty()) {
        add_variant(out_variant);
    }
}

tuple<int64_t, int64_t, bool, step_handle_t, step_handle_t> VCFOutputCaller::get_ref_interval(
    const PathPositionHandleGraph& graph, const Snarl& snarl, const string& ref_path_name) const {
    path_handle_t path_handle = graph.get_path_handle(ref_path_name);

    handle_t start_handle = graph.get_handle(snarl.start().node_id(), snarl.start().backward());
    map<size_t, step_handle_t> start_steps;
    graph.for_each_step_on_handle(start_handle, [&](step_handle_t step) {
            if (graph.get_path_handle_of_step(step) == path_handle) {
                start_steps[graph.get_position_of_step(step)] = step;
            }
        });

    handle_t end_handle = graph.get_handle(snarl.end().node_id(), snarl.end().backward());
    map<size_t, step_handle_t> end_steps;
    graph.for_each_step_on_handle(end_handle, [&](step_handle_t step) {
            if (graph.get_path_handle_of_step(step) == path_handle) {
                end_steps[graph.get_position_of_step(step)] = step;
            }
        });

    assert(start_steps.size() > 0 && end_steps.size() > 0);
    step_handle_t start_step = start_steps.begin()->second;
    step_handle_t end_step = end_steps.begin()->second;
    // just because we found a pair of steps on our path that correspond to the snarl ends, doesn't
    // mean the path threads the snarl.  verify that we can actaully walk, either forwards or backwards
    // along the path from the start node and hit then end node in the right orientation. 
    bool start_rev = graph.get_is_reverse(graph.get_handle_of_step(start_step)) != snarl.start().backward();
    bool end_rev = graph.get_is_reverse(graph.get_handle_of_step(end_step)) != snarl.end().backward();
    bool found_end = start_rev == end_rev && start_rev == start_steps.begin()->first > end_steps.begin()->first;
        
    // if we're on a cycle, we keep our start step and find the end step by scanning the path
    if (start_steps.size() > 1 || end_steps.size() > 1) {
        found_end = false;
        // try each start step
        for (auto i = start_steps.begin(); i != start_steps.end() && !found_end; ++i) {
            start_step = i->second;
            bool scan_backward = graph.get_is_reverse(graph.get_handle_of_step(start_step)) != snarl.start().backward();
            if (scan_backward) {
                // if we're going backward, we expect to reach the end backward
                end_handle = graph.get_handle(snarl.end().node_id(), !snarl.end().backward());
            }            
            if (scan_backward) {
                for (step_handle_t cur_step = start_step; graph.has_previous_step(cur_step) && !found_end;
                     cur_step = graph.get_previous_step(cur_step)) {
                    if (graph.get_handle_of_step(cur_step) == end_handle) {
                        end_step = cur_step;
                        found_end = true;
                    }
                }
            } else {
                for (step_handle_t cur_step = start_step; graph.has_next_step(cur_step) && !found_end;
                     cur_step = graph.get_next_step(cur_step)) {
                    if (graph.get_handle_of_step(cur_step) == end_handle) {
                        end_step = cur_step;
                        found_end = true;
                    }
                }
            }
        }
    }
    int64_t start_position = start_steps.begin()->first;
    step_handle_t out_start_step = start_step;
    int64_t end_position = end_step == end_steps.begin()->second ? end_steps.begin()->first : graph.get_position_of_step(end_step);
    step_handle_t out_end_step = end_step == end_steps.begin()->second ? end_steps.begin()->second : end_step;
    bool backward = end_position < start_position;
    

    if (!found_end) {
        // oops, once of the above checks failed.  we tell caller we coudlnt find by hacking in a -1 coordinate.
        start_position = -1;
        end_position = -1;
    }

    if (backward) {
        return make_tuple(end_position, start_position, backward, out_end_step, out_start_step);
    } else {
        return make_tuple(start_position, end_position, backward, out_start_step, out_end_step);
    }
}

pair<string, int64_t> VCFOutputCaller::get_ref_position(const PathPositionHandleGraph& graph, const Snarl& snarl, const string& ref_path_name,
                                                        int64_t ref_path_offset) const {
    string basepath_name = ref_path_name;
    size_t basepath_offset = 0;
    auto subpath_info = Paths::parse_subpath_name(ref_path_name);
    if (get<0>(subpath_info)) {
        basepath_name = get<1>(subpath_info);
        basepath_offset = get<2>(subpath_info);
    }
    // +1 to convert to 1-based VCF
    int64_t position = get<0>(get_ref_interval(graph, snarl, ref_path_name)) + ref_path_offset + 1 + basepath_offset;
    return make_pair(basepath_name, position);
}

void VCFOutputCaller::flatten_common_allele_ends(vcflib::Variant& variant, bool backward, size_t len_override) const {
    if (variant.alt.size() == 0) {
        return;
    }

    // find the minimum allele length to make sure we don't delete an entire allele
    size_t min_allele_len = variant.alleles[0].length();
    for (int i = 1; i < variant.alleles.size(); ++i) {
        min_allele_len = std::min(min_allele_len, variant.alleles[i].length());
    }

    // the maximum number of bases we want ot zip up, applying override if provided
    size_t max_flatten_len = len_override > 0 ? len_override : min_allele_len;
    
    // want to leave at least one in the reference position
    if (max_flatten_len == min_allele_len) {
        --max_flatten_len;
    }
    
    bool match = true;
    int shared_prefix_len = 0;
    for (int i = 0; i < max_flatten_len && match; ++i) {
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

string VCFOutputCaller::print_snarl(const Snarl& snarl, bool in_brackets) const {
    // todo, should we canonicalize here by putting lexicographic lowest node first?
    nid_t start_node = snarl.start().node_id();
    nid_t end_node = snarl.end().node_id();
    if (translation) {
        auto i = translation->find(start_node);
        if (i == translation->end()) {
            throw runtime_error("Error [VCFOutputCaller]: Unable to find node " + std::to_string(start_node) + " in translation file");
        }
        start_node = i->second.first;
        i = translation->find(end_node);
        if (i == translation->end()) {
            throw runtime_error("Error [VCFOutputCaller]: Unable to find node " + std::to_string(end_node) + " in translation file");
        }
        end_node = i->second.first;
    }
    stringstream ss;
    if (in_brackets) {
        ss << "(";
    }
    ss << (snarl.start().backward() ? "<" : ">") << start_node << (snarl.end().backward() ? "<" : ">") << end_node;
    if (in_brackets) {
        ss << ")";
    }
    return ss.str();
}

void VCFOutputCaller::scan_snarl(const string& allele_string, function<void(const string&, Snarl&)> callback) const {
    int left = -1;
    int last = 0;
    Snarl snarl;
    string frag;
    for (int i = 0; i < allele_string.length(); ++i) {
        if (allele_string[i] == '(') {
            assert(left == -1);
            if (last < i) {
                frag = allele_string.substr(last, i-last);
                callback(frag, snarl);
            }
            left = i;
        } else if (allele_string[i] == ')') {
            assert(left >= 0 && i > left + 3);
            frag = allele_string.substr(left + 1, i - left - 1);
            auto toks = split_delims(frag, "><");
            assert(toks.size() == 2);
            assert(frag[0] == '<' || frag[0] == '>');
            int64_t start = std::stoi(toks[0]);
            snarl.mutable_start()->set_node_id(start);
            snarl.mutable_start()->set_backward(frag[0] == '<');
            assert(frag[toks[0].size() + 1] == '<' || frag[toks[0].size() + 1] == '<');
            int64_t end = std::stoi(toks[1]);
            snarl.mutable_end()->set_node_id(abs(end));
            snarl.mutable_end()->set_backward(frag[toks[0].size() + 1] == '<');
            callback("", snarl);
            left = -1;
            last = i + 1;
        }
    }
    if (last == 0) {
        callback(allele_string, snarl);
    } else {
        frag = allele_string.substr(last);
        callback(frag, snarl);
    }
}

GAFOutputCaller::GAFOutputCaller(AlignmentEmitter* emitter, const string& sample_name, const vector<string>& ref_paths,
                                 size_t trav_padding) :
    emitter(emitter),
    gaf_sample_name(sample_name),
    ref_paths(ref_paths.begin(), ref_paths.end()),
    trav_padding(trav_padding) {
    
}

GAFOutputCaller::~GAFOutputCaller() {
}

void GAFOutputCaller::emit_gaf_traversals(const PathHandleGraph& graph, const string& snarl_name,
                                          const vector<SnarlTraversal>& travs,
                                          int64_t ref_trav_idx,
                                          const string& ref_path_name, int64_t ref_path_position,
                                          const TraversalSupportFinder* support_finder) {
    assert(emitter != nullptr);
    vector<Alignment> aln_batch;
    aln_batch.reserve(travs.size());

    stringstream ss;
    if (!ref_path_name.empty()) {
        ss << ref_path_name << "#" << ref_path_position << "#";
    }
    ss << snarl_name << "#" << gaf_sample_name;
    string variant_id = ss.str();

    // create allele ordering where reference is 0
    vector<int> alleles;
    if (ref_trav_idx >= 0) {
        alleles.push_back(ref_trav_idx);
    }
    for (int i = 0; i < travs.size(); ++i) {
        if (i != ref_trav_idx) {
            alleles.push_back(i);
        }
    }
    // make an alignment for each traversal
    for (int i = 0; i < alleles.size(); ++i) {
        const SnarlTraversal& trav = travs[alleles[i]];
        Alignment trav_aln;
        if (trav_padding > 0) {
            trav_aln = to_alignment(pad_traversal(graph, trav), graph);
        } else {
            trav_aln = to_alignment(trav, graph);
        }
        trav_aln.set_name(variant_id + "#" + std::to_string(i));
        if (support_finder) {
            int64_t support = support_finder->support_val(support_finder->get_traversal_support(trav));
            set_annotation(trav_aln, "support", std::to_string(support));
        }        
        aln_batch.push_back(trav_aln);
    }
    emitter->emit_singles(std::move(aln_batch)); 
}

void GAFOutputCaller::emit_gaf_variant(const PathHandleGraph& graph, const string& snarl_name,
                                       const vector<SnarlTraversal>& travs,
                                       const vector<int>& genotype,
                                       int64_t ref_trav_idx,
                                       const string& ref_path_name, int64_t ref_path_position,
                                       const TraversalSupportFinder* support_finder) {
    assert(emitter != nullptr);

    // pretty bare bones for now, just output the genotype as a pair of traversals
    // todo: we could embed some basic information (likelihood, ploidy, sample etc) in the gaf
    vector<SnarlTraversal> gt_travs;
    for (int allele : genotype) {
        gt_travs.push_back(travs[allele]);
    }
    emit_gaf_traversals(graph, snarl_name, gt_travs, ref_trav_idx, ref_path_name, ref_path_position, support_finder);
}

SnarlTraversal GAFOutputCaller::pad_traversal(const PathHandleGraph& graph, const SnarlTraversal& trav) const {

    assert(trav.visit_size() >= 2);

    SnarlTraversal out_trav;

    // traversal endpoints
    handle_t start_handle = graph.get_handle(trav.visit(0).node_id(), trav.visit(0).backward());
    handle_t end_handle = graph.get_handle(trav.visit(trav.visit_size() - 1).node_id(), trav.visit(trav.visit_size() - 1).backward());

    // find a reference path that touches the start node
    // todo: we could be more clever by finding the longest one or something
    path_handle_t reference_path;
    step_handle_t reference_step;
    bool found = false;
    size_t padding = 0;
    graph.for_each_step_on_handle(start_handle, [&](step_handle_t step_handle) {
            reference_path = graph.get_path_handle_of_step(step_handle);
            string name = graph.get_path_name(reference_path);
            if (!Paths::is_alt(name) && (ref_paths.empty() || ref_paths.count(name))) {
                reference_step = step_handle;
                found = true;
            }
            return !found;
        });

    // add left padding
    if (found) {
        deque<Visit> left_padding;

        if (graph.get_is_reverse(start_handle) == graph.get_is_reverse(graph.get_handle_of_step(reference_step))) {
            // path and handle oriented the same, we can just backtrack along the path to get previous stuff
            for (step_handle_t step = graph.get_previous_step(reference_step);
                 step != graph.path_front_end(reference_path) && padding < trav_padding;
                 step = graph.get_previous_step(step)) {
                left_padding.push_front(to_visit(graph, graph.get_handle_of_step(step)));
                padding += graph.get_length(graph.get_handle_of_step(step));
            }
        } else {
            // path and handle oriented differently, we go forward in the path, flipping each step
            for (step_handle_t step = graph.get_next_step(reference_step);
                 step != graph.path_end(reference_path) && padding < trav_padding;
                 step = graph.get_next_step(step)) {
                left_padding.push_front(to_visit(graph, graph.get_handle_of_step(step)));
                padding += graph.get_length(graph.get_handle_of_step(step));
            }
        }

        for (const Visit& visit : left_padding) {
            *out_trav.add_visit() = visit;
        }
    }

    // copy over center
    for (int i = 0; i < trav.visit_size(); ++i) {
        *out_trav.add_visit() = trav.visit(i);
    }

    // go through the whole thing again with the end
    found = false;
    padding = 0;
    graph.for_each_step_on_handle(end_handle, [&](step_handle_t step_handle) {
            reference_path = graph.get_path_handle_of_step(step_handle);
            string name = graph.get_path_name(reference_path);
            if (!Paths::is_alt(name) && (ref_paths.empty() || ref_paths.count(name))) {
                reference_step = step_handle;
                found = true;
            }
            return !found;
        });

    // add right padding
    if (found) {
        if (graph.get_is_reverse(end_handle) == graph.get_is_reverse(graph.get_handle_of_step(reference_step))) {
            // path and handle oriented the same, we can just continue along the path to get next stuff
            for (step_handle_t step = graph.get_next_step(reference_step);
                 step != graph.path_end(reference_path) && padding < trav_padding;
                 step = graph.get_next_step(step)) {
                Visit* visit = out_trav.add_visit();
                *visit = to_visit(graph, graph.get_handle_of_step(step));
                padding += graph.get_length(graph.get_handle_of_step(step));
            }
        } else {
            // path and handle oriented differently, we go backward in the path, flipping each step
            for (step_handle_t step = graph.get_previous_step(reference_step);
                 step != graph.path_front_end(reference_path) && padding < trav_padding;
                 step = graph.get_previous_step(step)) {
                Visit* visit = out_trav.add_visit();
                *visit = to_visit(graph, graph.flip(graph.get_handle_of_step(step)));
                padding += graph.get_length(graph.get_handle_of_step(step));
            }
        }
    }
        
    return out_trav;
}

void VCFOutputCaller::update_nesting_info_tags(const SnarlManager* snarl_manager) {

    // index the snarl tree by name
    unordered_map<string, const Snarl*> name_to_snarl;
    Snarl flipped_snarl;
    snarl_manager->for_each_snarl_preorder([&](const Snarl* snarl) {
            name_to_snarl[print_snarl(*snarl)] = snarl;
            // also add a map from the flipped snarl (as call sometimes messes with orientation)
            flipped_snarl.mutable_start()->set_node_id(snarl->end().node_id());
            flipped_snarl.mutable_start()->set_backward(!snarl->end().backward());
            flipped_snarl.mutable_end()->set_node_id(snarl->start().node_id());
            flipped_snarl.mutable_end()->set_backward(!snarl->start().backward());
            name_to_snarl[print_snarl(flipped_snarl)] = snarl;
        });

    // pass 1) index sites in vcf
    // (todo: this could be done more quickly upstream)
    unordered_set<string> names_in_vcf;
    for (auto& thread_buf : output_variants) {
        for (auto& output_variant_record : thread_buf) {
            string output_variant_string;
            zstdutil::DecompressString(output_variant_record.second, output_variant_string);
            vector<string> toks = split_delims(output_variant_string, "\t", 4);
            names_in_vcf.insert(toks[2]);
        }
    }

    // determine the tags from the index
    function<pair<size_t, string>(const string&)> get_lv_ps_tags = [&](const string& name) {
        string parent_name;
        size_t ancestor_count = 0;
        const Snarl* snarl = name_to_snarl.at(name);
        assert(snarl != nullptr);
        // walk up the snarl tree
        while (snarl = snarl_manager->parent_of(snarl)) {
            string cur_name = print_snarl(*snarl);
            if (names_in_vcf.count(cur_name)) {
                // only count snarls that are in the vcf
                ++ancestor_count;
                if (parent_name.empty()) {
                    // remembert the first parent
                    parent_name = cur_name;
                }
            }
        }
        return make_pair(ancestor_count, parent_name);
    };
    
    // pass 2) add the LV and PS tags
#pragma omp parallel for
    for (uint64_t i = 0; i < output_variants.size(); ++i) {
        auto& thread_buf = output_variants[i];
        for (auto& output_variant_record : thread_buf) {
            string output_variant_string;
            zstdutil::DecompressString(output_variant_record.second, output_variant_string);
            //string& output_variant_string = output_variant_record.second;
            vector<string> toks = split_delims(output_variant_string, "\t", 9);
            const string& name = toks[2];
            
            pair<size_t, string> lv_ps = get_lv_ps_tags(name);
            string nesting_tags = ";LV=" + std::to_string(lv_ps.first);
            if (lv_ps.first != 0) {
                assert(!lv_ps.second.empty());
                nesting_tags += ";PS=" + lv_ps.second;
            }

            // rewrite the output string using the updated info toks
            output_variant_string.clear();
            for (size_t i = 0; i < toks.size(); ++i) {
                output_variant_string += toks[i];
                if (i == 7) {
                    output_variant_string += nesting_tags;
                }
                if (i != toks.size() - 1) {
                    output_variant_string += "\t";
                }
            }
            output_variant_record.second.clear();
            zstdutil::CompressString(output_variant_string, output_variant_record.second);
        }
    }
}

VCFGenotyper::VCFGenotyper(const PathHandleGraph& graph,
                           SnarlCaller& snarl_caller,
                           SnarlManager& snarl_manager,
                           vcflib::VariantCallFile& variant_file,
                           const string& sample_name,
                           const vector<string>& ref_paths,
                           const vector<int>& ref_path_ploidies,
                           FastaReference* ref_fasta,
                           FastaReference* ins_fasta,
                           AlignmentEmitter* aln_emitter,
                           bool traversals_only,
                           bool gaf_output,
                           size_t trav_padding) :
    GraphCaller(snarl_caller, snarl_manager),
    VCFOutputCaller(sample_name),
    GAFOutputCaller(aln_emitter, sample_name, ref_paths, trav_padding),
    graph(graph),
    input_vcf(variant_file),
    traversal_finder(graph, snarl_manager, variant_file, ref_paths, ref_fasta, ins_fasta, snarl_caller.get_skip_allele_fn()),
    traversals_only(traversals_only),
    gaf_output(gaf_output) {

    scan_contig_lengths();

    assert(ref_paths.size() == ref_path_ploidies.size());
    for (int i = 0; i < ref_paths.size(); ++i) {
        path_to_ploidy[ref_paths[i]] = ref_path_ploidies[i];
    }
}

VCFGenotyper::~VCFGenotyper() {

}

bool VCFGenotyper::call_snarl(const Snarl& snarl) {

    // could be that our graph is a subgraph of the graph the snarls were computed from
    // so bypass snarls we can't process
    if (!graph.has_node(snarl.start().node_id()) || !graph.has_node(snarl.end().node_id())) {
        return false;
    }

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

        // find a path range corresponding to our snarl by way of the VCF variants.
        tuple<string, size_t, size_t> ref_positions = get_ref_positions(variants);

        // just print the traversals if requested
        if (traversals_only) {
            assert(gaf_output);
            // todo: can't get ref position here without pathposition graph
            emit_gaf_traversals(graph, print_snarl(snarl), travs, ref_trav_idx, "", -1);
            return true;
        }
        
        // use our support caller to choose our genotype (int traversal coordinates)
        vector<int> trav_genotype;
        unique_ptr<SnarlCaller::CallInfo> trav_call_info;
        std::tie(trav_genotype, trav_call_info) = snarl_caller.genotype(snarl, travs, ref_trav_idx, path_to_ploidy[get<0>(ref_positions)],
                                                                        get<0>(ref_positions),  make_pair(get<1>(ref_positions), get<2>(ref_positions)));

        assert(trav_genotype.size() <= 2);

        if (gaf_output) {
            // todo: can't get ref position here without pathposition graph
            emit_gaf_variant(graph, print_snarl(snarl), travs, trav_genotype, ref_trav_idx, "", -1);
            return true;
        }

        // map our genotype back to the vcf
        for (int i = 0; i < variants.size(); ++i) {
            vector<int> vcf_alleles;
            set<int> used_vcf_alleles;
            string vcf_genotype;
            vector<SnarlTraversal> vcf_traversals(variants[i]->alleles.size());            
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
                    used_vcf_alleles.insert(vcf_allele);
                    vcf_traversals[vcf_allele] = travs[trav_allele];
                }
                // add traversals that correspond to vcf genotypes that are not
                // present in the traversal_genotypes
                for (int j = 0; j < travs.size(); ++j) {
                    int vcf_allele = alleles[j].second[i];
                    if (!used_vcf_alleles.count(vcf_allele)) {
                        vcf_traversals[vcf_allele] = travs[j];
                        used_vcf_alleles.insert(vcf_allele);
                    }
                }
            }
            // create an output variant from the input one
            vcflib::Variant out_variant;
            out_variant.sequenceName = variants[i]->sequenceName;
            out_variant.position = variants[i]->position;
            out_variant.id = variants[i]->id;
            out_variant.ref = variants[i]->ref;
            out_variant.alt = variants[i]->alt;
            out_variant.alleles = variants[i]->alleles;
            out_variant.filter = "PASS";
            out_variant.updateAlleleIndexes();

            // add the genotype
            out_variant.format.push_back("GT");
            auto& genotype_vector = out_variant.samples[sample_name]["GT"];
            genotype_vector.push_back(vcf_genotype);

            // add some info
            snarl_caller.update_vcf_info(snarl, vcf_traversals, vcf_alleles, trav_call_info, sample_name, out_variant);

            // print the variant
            add_variant(out_variant);
        }
        return true;
    }
    
    return false;

}

string VCFGenotyper::vcf_header(const PathHandleGraph& graph, const vector<string>& ref_paths,
                                const vector<size_t>& contig_length_overrides) const {
    assert(contig_length_overrides.empty()); // using this override makes no sense

    // get the contig length overrides from the VCF
    vector<size_t> vcf_contig_lengths;
    auto length_map = scan_contig_lengths();
    for (int i = 0; i < ref_paths.size(); ++i) {
        vcf_contig_lengths.push_back(length_map[ref_paths[i]]);
    }
    
    string header = VCFOutputCaller::vcf_header(graph, ref_paths, vcf_contig_lengths);
    header += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    snarl_caller.update_vcf_header(header);
    header += "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
    header += "##SAMPLE=<ID=" + sample_name + ">\n";
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_name;
    assert(output_vcf.openForOutput(header));
    header += "\n";
    return header;
}

tuple<string, size_t, size_t> VCFGenotyper::get_ref_positions(const vector<vcflib::Variant*>& variants) const {
    // if there is more than one path in our snarl (unlikely for most graphs we'll vcf-genoetype)
    // then we return the one with the biggest interval
    map<string, pair<size_t, size_t>> path_offsets;
    for (const vcflib::Variant* var : variants) {
        if (path_offsets.count(var->sequenceName)) {
            pair<size_t, size_t>& record = path_offsets[var->sequenceName];
            record.first = std::min((size_t)var->position, record.first);
            record.second = std::max((size_t)var->position + var->ref.length(), record.second);
        } else {
            path_offsets[var->sequenceName] = make_pair(var->position, var->position + var->ref.length());
        }
    }

    string ref_path;
    size_t ref_range_size = 0;
    pair<size_t, size_t> ref_range;
    for (auto& path_offset : path_offsets) {
        size_t len = path_offset.second.second - path_offset.second.first;
        if (len > ref_range_size) {
            ref_range_size = len;
            ref_path = path_offset.first;
            ref_range = path_offset.second;
        }
    }

    return make_tuple(ref_path, ref_range.first, ref_range.second);
}

unordered_map<string, size_t> VCFGenotyper::scan_contig_lengths() const {

    unordered_map<string, size_t> ref_lengths;
    
    // copied from dumpContigsFromHeader.cpp in vcflib
    vector<string> headerLines = split(input_vcf.header, "\n");
    for(vector<string>::iterator it = headerLines.begin(); it != headerLines.end(); it++) {
        if((*it).substr(0,8) == "##contig"){
            string contigInfo = (*it).substr(10, (*it).length() -11);
            vector<string> info = split(contigInfo, ",");
            string id;
            int64_t length = -1;
            for(vector<string>::iterator sub = info.begin(); sub != info.end(); sub++) {
                vector<string> subfield = split((*sub), "=");
                if(subfield[0] == "ID"){
                    id = subfield[1];
                }
                if(subfield[0] == "length"){
                    length = parse<int>(subfield[1]);
                }
            }
            if (!id.empty() && length >= 0) {
                ref_lengths[id] = length;
            }
        }
    }

    return ref_lengths;
}


LegacyCaller::LegacyCaller(const PathPositionHandleGraph& graph,
                           SupportBasedSnarlCaller& snarl_caller,
                           SnarlManager& snarl_manager,
                           const string& sample_name,
                           const vector<string>& ref_paths,
                           const vector<size_t>& ref_path_offsets,
                           const vector<int>& ref_path_ploidies) :
    GraphCaller(snarl_caller, snarl_manager),
    VCFOutputCaller(sample_name),
    graph(graph),
    ref_paths(ref_paths) {

    for (int i = 0; i < ref_paths.size(); ++i) {
        ref_offsets[ref_paths[i]] = i < ref_path_offsets.size() ? ref_path_offsets[i] : 0;
        ref_ploidies[ref_paths[i]] = i < ref_path_ploidies.size() ? ref_path_ploidies[i] : 2;
    }
    
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
                                                             0,
                                                             0,
                                                             get_path_index,
                                                             [&](id_t id) { return snarl_caller.get_support_finder().get_min_node_support(id);},
                                                             [&](edge_t edge) { return snarl_caller.get_support_finder().get_edge_support(edge);});

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
    bool was_called = false;
    
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
        size_t total_snarl_length = 0;
        for (auto node_id : contents.first) {
            handle_t new_handle = vg_graph.create_handle(graph.get_sequence(graph.get_handle(node_id)), node_id);
            if (node_id != snarl.start().node_id() && node_id != snarl.end().node_id()) {
                total_snarl_length += vg_graph.get_length(new_handle);
            }
        }
        for (auto edge : contents.second) {
            vg_graph.create_edge(vg_graph.get_handle(graph.get_id(edge.first), vg_graph.get_is_reverse(edge.first)),
                                 vg_graph.get_handle(graph.get_id(edge.second), vg_graph.get_is_reverse(edge.second)));
            total_snarl_length += 1;
        }
        // add the paths to the subgraph
        algorithms::expand_context_with_paths(&graph, &vg_graph, 1);
        // and index them
        for (auto& ref_path : ref_paths) {
            if (vg_graph.has_path(ref_path)) {
                site_path_indexes.push_back(new PathIndex(vg_graph, ref_path));
            } else {
                site_path_indexes.push_back(nullptr);
            }
        }
        get_path_index = [&](const Snarl& site) -> PathIndex* {
            return find_index(site, site_path_indexes).second;
        };
        // determine the support threshold for the traversal finder.  if we're using average
        // support, then we don't use any (set to 0), other wise, use the minimum support for a call
        SupportBasedSnarlCaller& support_caller = dynamic_cast<SupportBasedSnarlCaller&>(snarl_caller);
        size_t threshold = support_caller.get_support_finder().get_average_traversal_support_switch_threshold();
        double support_cutoff = total_snarl_length <= threshold ? support_caller.get_min_total_support_for_call() : 0;
        rep_trav_finder = new RepresentativeTraversalFinder(vg_graph, snarl_manager,
                                                            max_search_depth,
                                                            max_search_width,
                                                            max_bubble_paths,
                                                            support_cutoff,
                                                            support_cutoff,
                                                            get_path_index,
                                                            [&](id_t id) { return support_caller.get_support_finder().get_min_node_support(id);},
                                                            // note: because our traversal finder and support caller have
                                                            // different graphs, they can't share edge handles
                                                            [&](edge_t edge) { return support_caller.get_support_finder().get_edge_support(
                                                                    vg_graph.get_id(edge.first), vg_graph.get_is_reverse(edge.first),
                                                                    vg_graph.get_id(edge.second), vg_graph.get_is_reverse(edge.second));});
                                                            
    }

    PathIndex* path_index = get_path_index(snarl);
    if (path_index != nullptr) {
        string path_name = find_index(snarl, is_vg ? path_indexes : site_path_indexes).first;

        // orient the snarl along the reference path
        tuple<size_t, size_t, bool, step_handle_t, step_handle_t> ref_interval = get_ref_interval(graph, snarl, path_name);
        if (get<2>(ref_interval) == true) {
            snarl_manager.flip(&snarl);
        }

        // recursively genotype the site beginning here at the top level snarl
        vector<SnarlTraversal> called_traversals;
        // these integers map the called traversals to their positions in the list of all traversals
        // of the top level snarl.  
        vector<int> genotype;
        int ploidy = ref_ploidies[path_name];
        std::tie(called_traversals, genotype) = top_down_genotype(snarl, *rep_trav_finder, ploidy,
                                                                  path_name, make_pair(get<0>(ref_interval), get<1>(ref_interval)));
    
        if (!called_traversals.empty()) {
            // regenotype our top-level traversals now that we know they aren't nested, and we have a
            // good idea of all the sizes
            unique_ptr<SnarlCaller::CallInfo> call_info;
            std::tie(called_traversals, genotype, call_info) = re_genotype(snarl, *rep_trav_finder, called_traversals, genotype, ploidy,
                                                                           path_name, make_pair(get<0>(ref_interval), get<1>(ref_interval)));

            // emit our vcf variant
            emit_variant(graph, snarl_caller, snarl, called_traversals, genotype, 0, call_info, path_name, ref_offsets.find(path_name)->second, false,
                         ploidy);

            was_called = true;
        }
    }        
    if (!is_vg) {
        // delete the temporary vg subgraph and traversal finder we created for this snarl
        delete rep_trav_finder;
        for (PathIndex* path_index : site_path_indexes) {
            delete path_index;
        }
    }

    return was_called;
}

string LegacyCaller::vcf_header(const PathHandleGraph& graph, const vector<string>& ref_paths,
                                const vector<size_t>& contig_length_overrides) const {
    string header = VCFOutputCaller::vcf_header(graph, ref_paths, contig_length_overrides);
    header += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    snarl_caller.update_vcf_header(header);
    header += "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
    header += "##SAMPLE=<ID=" + sample_name + ">\n";
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_name;
    assert(output_vcf.openForOutput(header));
    header += "\n";
    return header;
}

pair<vector<SnarlTraversal>, vector<int>> LegacyCaller::top_down_genotype(const Snarl& snarl, TraversalFinder& trav_finder, int ploidy,
                                                                          const string& ref_path_name, pair<size_t, size_t> ref_interval) const {

    // get the traversals through the site
    vector<SnarlTraversal> traversals = trav_finder.find_traversals(snarl);

    // use our support caller to choose our genotype
    vector<int> trav_genotype;
    unique_ptr<SnarlCaller::CallInfo> trav_call_info;
    std::tie(trav_genotype, trav_call_info) = snarl_caller.genotype(snarl, traversals, 0, ploidy, ref_path_name, ref_interval);
    if (trav_genotype.empty()) {
        return make_pair(vector<SnarlTraversal>(), vector<int>());
    }

    assert(trav_genotype.size() == ploidy);

    vector<SnarlTraversal> called_travs(ploidy);

    // do we have two paths going through a given traversal?  This is handled
    // as a special case below
    bool hom = trav_genotype.size() == 2 && trav_genotype[0] == trav_genotype[1];
    
    for (int i = 0; i < trav_genotype.size() && (!hom || i < 1); ++i) {
        int allele = trav_genotype[i];
        const SnarlTraversal& traversal = traversals[allele];
        Visit prev_end;
        for (int j = 0; j < traversal.visit_size(); ++j) {
            if (traversal.visit(j).node_id() > 0) {
                *called_travs[i].add_visit() = traversal.visit(j);
                if (hom && i == 0) {
                    *called_travs[1].add_visit() = traversal.visit(j);
                }
            } else {
                // recursively determine the traversal
                const Snarl* into_snarl = snarl_manager.into_which_snarl(traversal.visit(j));
                bool flipped = traversal.visit(j).backward();
                if (flipped) {
                    // we're always processing our snarl from start to end, so make sure
                    // it lines up with the parent (note that we've oriented the root along the ref path)
                    snarl_manager.flip(into_snarl);
                }
                vector<SnarlTraversal> child_genotype = top_down_genotype(*into_snarl,
                                                                          trav_finder, hom ? 2: 1, ref_path_name, ref_interval).first;                
                if (child_genotype.empty()) {
                    return make_pair(vector<SnarlTraversal>(), vector<int>());
                }
                bool back_to_back = j > 0 && traversal.visit(j - 1).node_id() == 0 && prev_end == into_snarl->start();

                for (int k = back_to_back ? 1 : 0; k < child_genotype[0].visit_size(); ++k) {
                    *called_travs[i].add_visit() = child_genotype[0].visit(k);
                }
                if (hom) {
                    assert(child_genotype.size() == 2 && i == 0);
                    for (int k = back_to_back ? 1 : 0; k < child_genotype[1].visit_size(); ++k) {
                        *called_travs[1].add_visit() = child_genotype[1].visit(k);
                    }
                }
                prev_end = into_snarl->end();
                if (flipped) {
                    // leave our snarl like we found it
                    snarl_manager.flip(into_snarl);
                }
            }
        }
    }

    return make_pair(called_travs, trav_genotype);
}

SnarlTraversal LegacyCaller::get_reference_traversal(const Snarl& snarl, TraversalFinder& trav_finder) const {

    // get the ref traversal through the site
    // todo: don't avoid so many traversal recomputations
    SnarlTraversal traversal = trav_finder.find_traversals(snarl)[0];
    SnarlTraversal out_traversal;

    Visit prev_end;
    for (int i = 0; i < traversal.visit_size(); ++i) {
        const Visit& visit = traversal.visit(i);
        if (visit.node_id() != 0) {
            *out_traversal.add_visit() = visit;
        } else {
            const Snarl* into_snarl = snarl_manager.into_which_snarl(visit);
            if (visit.backward()) {
                snarl_manager.flip(into_snarl);
            }
            bool back_to_back = i > 0 && traversal.visit(i - 1).node_id() == 0 && prev_end == into_snarl->start();

            SnarlTraversal child_ref = get_reference_traversal(*into_snarl, trav_finder);
            for (int j = back_to_back ? 1 : 0; j < child_ref.visit_size(); ++j) {
                *out_traversal.add_visit() = child_ref.visit(j);
            }
            prev_end = into_snarl->end();
            if (visit.backward()) {
                // leave our snarl like we found it
                snarl_manager.flip(into_snarl);
            }
        }
    }
    return out_traversal;    
}

tuple<vector<SnarlTraversal>, vector<int>, unique_ptr<SnarlCaller::CallInfo>>
LegacyCaller::re_genotype(const Snarl& snarl, TraversalFinder& trav_finder,
                          const vector<SnarlTraversal>& in_traversals,
                          const vector<int>& in_genotype,
                          int ploidy,
                          const string& ref_path_name,
                          pair<size_t, size_t> ref_interval) const {
    
    assert(in_traversals.size() == in_genotype.size());
    
    // create a set of unique traversal candidates that must include the reference first
    vector<SnarlTraversal> rg_traversals;
    // add our reference traversal to the front
    for (int i = 0; i < in_traversals.size() && !rg_traversals.empty(); ++i) {
        if (in_genotype[i] == 0) {
            rg_traversals.push_back(in_traversals[i]);
        }
    }
    if (rg_traversals.empty()) {
        rg_traversals.push_back(get_reference_traversal(snarl, trav_finder));
    }
    set<int> gt_set = {0};
    for (int i = 0; i < in_traversals.size(); ++i) {
        if (!gt_set.count(in_genotype[i])) {
            rg_traversals.push_back(in_traversals[i]);
            gt_set.insert(in_genotype[i]);
        }
    }
    
    // re-genotype the candidates
    vector<int> rg_genotype;
    unique_ptr<SnarlCaller::CallInfo> rg_call_info;
    std::tie(rg_genotype, rg_call_info) = snarl_caller.genotype(snarl, rg_traversals, 0, ploidy, ref_path_name, ref_interval);

    return make_tuple(rg_traversals, rg_genotype, std::move(rg_call_info));
}

bool LegacyCaller::is_traversable(const Snarl& snarl) {
    // we need this to be true all the way down to use the RepresentativeTraversalFinder on our snarl.
    bool ret = snarl.start_end_reachable() && snarl.directed_acyclic_net_graph() &&
       graph.has_node(snarl.start().node_id()) && graph.has_node(snarl.end().node_id());
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
        if (path_index != nullptr &&
            path_index->by_id.count(snarl.start().node_id()) &&
            path_index->by_id.count(snarl.end().node_id())) {
            // This path threads through this site
            return make_pair(ref_paths[i], path_index);
        }
    }
    return make_pair("", nullptr);
}

FlowCaller::FlowCaller(const PathPositionHandleGraph& graph,
                       SupportBasedSnarlCaller& snarl_caller,
                       SnarlManager& snarl_manager,
                       const string& sample_name,
                       TraversalFinder& traversal_finder,
                       const vector<string>& ref_paths,
                       const vector<size_t>& ref_path_offsets,
                       const vector<int>& ref_path_ploidies,
                       AlignmentEmitter* aln_emitter,
                       bool traversals_only,
                       bool gaf_output,
                       size_t trav_padding,
                       bool genotype_snarls,
                       const pair<int64_t, int64_t>& ref_allele_length_range) :
    GraphCaller(snarl_caller, snarl_manager),
    VCFOutputCaller(sample_name),
    GAFOutputCaller(aln_emitter, sample_name, ref_paths, trav_padding),
    graph(graph),
    traversal_finder(traversal_finder),
    ref_paths(ref_paths),
    traversals_only(traversals_only),
    gaf_output(gaf_output),
    genotype_snarls(genotype_snarls),
    ref_allele_length_range(ref_allele_length_range)
{
    for (int i = 0; i < ref_paths.size(); ++i) {
        ref_offsets[ref_paths[i]] = i < ref_path_offsets.size() ? ref_path_offsets[i] : 0;
        ref_path_set.insert(ref_paths[i]);
        ref_ploidies[ref_paths[i]] = i < ref_path_ploidies.size() ? ref_path_ploidies[i] : 2;
    }

}
   
FlowCaller::~FlowCaller() {

}

bool FlowCaller::call_snarl(const Snarl& managed_snarl) {

    // todo: In order to experiment with merging consecutive snarls to make longer traversals,
    // I am experimenting with sending "fake" snarls through this code.  So make a local
    // copy to work on to do things like flip -- calling any snarl_manager code that
    // wants a pointer will crash. 
    Snarl snarl = managed_snarl;

    if (snarl.start().node_id() == snarl.end().node_id() ||
        !graph.has_node(snarl.start().node_id()) || !graph.has_node(snarl.end().node_id())) {
        // can't call one-node or out-of graph snarls.
        return false;
    }
    // toggle average flow / flow width based on snarl length.  this is a bit inconsistent with
    // downstream which uses the longest traversal length, but it's a bit chicken and egg
    // todo: maybe use snarl length for everything?
    const auto& support_finder = dynamic_cast<SupportBasedSnarlCaller&>(snarl_caller).get_support_finder();    
    bool greedy_avg_flow = false;
    {
        auto snarl_contents = snarl_manager.deep_contents(&snarl, graph, false);
        if (snarl_contents.second.size() > max_snarl_edges) {
            // size cap needed as FlowCaller doesn't have nesting support yet
            return false;
        }        
        size_t len_threshold = support_finder.get_average_traversal_support_switch_threshold();
        size_t length = 0;
        for (auto i = snarl_contents.first.begin(); i != snarl_contents.first.end() && length < len_threshold; ++i) {
            length += graph.get_length(graph.get_handle(*i));
        }
        greedy_avg_flow = length > len_threshold;
    }
    
    handle_t start_handle = graph.get_handle(snarl.start().node_id(), snarl.start().backward());
    handle_t end_handle = graph.get_handle(snarl.end().node_id(), snarl.end().backward());

    // as we're writing to VCF, we need a reference path through the snarl.  we
    // look it up directly from the graph, and abort if we can't find one
    set<string> start_path_names;
    graph.for_each_step_on_handle(start_handle, [&](step_handle_t step_handle) {
            string name = graph.get_path_name(graph.get_path_handle_of_step(step_handle));
            if (!Paths::is_alt(name) && (ref_path_set.empty() || ref_path_set.count(name))) {
                start_path_names.insert(name);
            }
            return true;
        });
    
    set<string> end_path_names;
    if (!start_path_names.empty()) {
        graph.for_each_step_on_handle(end_handle, [&](step_handle_t step_handle) {
                string name = graph.get_path_name(graph.get_path_handle_of_step(step_handle));
                if (!Paths::is_alt(name) && (ref_path_set.empty() || ref_path_set.count(name))) {                
                    end_path_names.insert(name);
                }
                return true;
            });
    }
    
    // we do the full intersection (instead of more quickly finding the first common path)
    // so that we always take the lexicographically lowest path, rather than depending
    // on the order of iteration which could change between implementations / runs.
    vector<string> common_names;
    std::set_intersection(start_path_names.begin(), start_path_names.end(),
                          end_path_names.begin(), end_path_names.end(),
                          std::back_inserter(common_names));

    if (common_names.empty()) {
        return false;
    }

    string& ref_path_name = common_names.front();

    // find the reference traversal and coordinates using the path position graph interface
    tuple<int64_t, int64_t, bool, step_handle_t, step_handle_t> ref_interval = get_ref_interval(graph, snarl, ref_path_name);
    if (get<0>(ref_interval) == -1) {
        // could not find reference path interval consisten with snarl due to orientation conflict
        return false;
    }
    if (get<2>(ref_interval) == true) {
        // calling code assumes snarl forward on reference
        flip_snarl(snarl);
        ref_interval = get_ref_interval(graph, snarl, ref_path_name);
    }

    step_handle_t cur_step = get<3>(ref_interval);
    step_handle_t last_step = get<4>(ref_interval);
    if (get<2>(ref_interval)) {
        std::swap(cur_step, last_step);
    }
    bool start_backwards = snarl.start().backward() != graph.get_is_reverse(graph.get_handle_of_step(cur_step));
    
    SnarlTraversal ref_trav;
    while (true) {
        handle_t cur_handle = graph.get_handle_of_step(cur_step);
        Visit* visit = ref_trav.add_visit();
        visit->set_node_id(graph.get_id(cur_handle));
        visit->set_backward(start_backwards ? !graph.get_is_reverse(cur_handle) : graph.get_is_reverse(cur_handle));
        if (graph.get_id(cur_handle) == snarl.end().node_id()) {
            break;
        } else if (get<2>(ref_interval) == true) {
            if (!graph.has_previous_step(cur_step)) {
                cerr << "Warning [vg call]: Unable, due to bug or corrupt path information, to trace reference path through snarl " << pb2json(managed_snarl) << endl;
                return false;
            }
            cur_step = graph.get_previous_step(cur_step);
        } else {
            if (!graph.has_next_step(cur_step)) {
                cerr << "Warning [vg call]: Unable, due to bug or corrupt path information, to trace reference path through snarl " << pb2json(managed_snarl) << endl;
                return false;
            }
            cur_step = graph.get_next_step(cur_step);
        }
        // todo: we can compute flow at the same time
    }
    assert(ref_trav.visit(0) == snarl.start() && ref_trav.visit(ref_trav.visit_size() - 1) == snarl.end());

    // optional reference length clamp can, ex, avoid trying to resolve a giant snarl
    if (ref_trav.visit_size() > 1 && ref_allele_length_range.first > 0 || ref_allele_length_range.second < numeric_limits<int64_t>::max()) {
        size_t ref_trav_len = 0;
        for (size_t j = 1; j < ref_trav.visit_size() - 1; ++j) {
            ref_trav_len += graph.get_length(graph.get_handle(ref_trav.visit(j).node_id()));
        }
        if (ref_trav_len < ref_allele_length_range.first || ref_trav_len > ref_allele_length_range.second) {
            return false;
        }        
    }

    vector<SnarlTraversal> travs;
    FlowTraversalFinder* flow_trav_finder = dynamic_cast<FlowTraversalFinder*>(&traversal_finder);
    if (flow_trav_finder != nullptr) {
        // find the max flow traversals using specialized interface that accepts avg heurstic toggle
        pair<vector<SnarlTraversal>, vector<double>> weighted_travs = flow_trav_finder->find_weighted_traversals(snarl, greedy_avg_flow);
        travs = std::move(weighted_travs.first);
    } else {
        // find the traversals using the generic interface
        travs = traversal_finder.find_traversals(snarl);
    }

    if (travs.empty()) {
        cerr << "Warning [vg call]: Unable, due to bug or corrupt graph, to search for any traversals through snarl " << pb2json(managed_snarl) << endl;
        return false;
    }

    // find the reference traversal in the list of results from the traversal finder
    int ref_trav_idx = -1;
    for (int i = 0; i < travs.size() && ref_trav_idx < 0; ++i) {
        // todo: is there a way to speed this up?
        if (travs[i] == ref_trav) {
            ref_trav_idx = i;
        }
    }

    if (ref_trav_idx == -1) {
        ref_trav_idx = travs.size();
        // we didn't get the reference traversal from the finder, so we add it here
        travs.push_back(ref_trav);
    }

    bool ret_val = true;

    if (traversals_only) {
        assert(gaf_output);
        pair<string, int64_t> pos_info = get_ref_position(graph, snarl, ref_path_name, ref_offsets[ref_path_name]);
        emit_gaf_traversals(graph, print_snarl(snarl), travs, ref_trav_idx, pos_info.first, pos_info.second, &support_finder);
    } else {
        // use our support caller to choose our genotype
        vector<int> trav_genotype;
        unique_ptr<SnarlCaller::CallInfo> trav_call_info;
        int ploidy = ref_ploidies[ref_path_name];
        std::tie(trav_genotype, trav_call_info) = snarl_caller.genotype(snarl, travs, ref_trav_idx, ploidy, ref_path_name,
                                                                        make_pair(get<0>(ref_interval), get<1>(ref_interval)));

        assert(trav_genotype.empty() || trav_genotype.size() == ploidy);

        if (!gaf_output) {
            emit_variant(graph, snarl_caller, snarl, travs, trav_genotype, ref_trav_idx, trav_call_info, ref_path_name,
                         ref_offsets[ref_path_name], genotype_snarls, ploidy);
        } else {
            pair<string, int64_t> pos_info = get_ref_position(graph, snarl, ref_path_name, ref_offsets[ref_path_name]);
            emit_gaf_variant(graph, print_snarl(snarl), travs, trav_genotype, ref_trav_idx, pos_info.first, pos_info.second, &support_finder);
        }

        ret_val = trav_genotype.size() == ploidy;
    }
        
    return ret_val;
}

string FlowCaller::vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                              const vector<size_t>& contig_length_overrides) const {
    string header = VCFOutputCaller::vcf_header(graph, contigs, contig_length_overrides);
    header += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    snarl_caller.update_vcf_header(header);
    header += "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
    header += "##SAMPLE=<ID=" + sample_name + ">\n";
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_name;
    assert(output_vcf.openForOutput(header));
    header += "\n";
    return header;
}


NestedFlowCaller::NestedFlowCaller(const PathPositionHandleGraph& graph,
                                   SupportBasedSnarlCaller& snarl_caller,
                                   SnarlManager& snarl_manager,
                                   const string& sample_name,
                                   TraversalFinder& traversal_finder,
                                   const vector<string>& ref_paths,
                                   const vector<size_t>& ref_path_offsets,
                                   const vector<int>& ref_path_ploidies,
                                   AlignmentEmitter* aln_emitter,
                                   bool traversals_only,
                                   bool gaf_output,
                                   size_t trav_padding,
                                   bool genotype_snarls) :
    GraphCaller(snarl_caller, snarl_manager),
    VCFOutputCaller(sample_name),
    GAFOutputCaller(aln_emitter, sample_name, ref_paths, trav_padding),
    graph(graph),
    traversal_finder(traversal_finder),
    ref_paths(ref_paths),
    traversals_only(traversals_only),
    gaf_output(gaf_output),
    genotype_snarls(genotype_snarls),
    nested_support_finder(dynamic_cast<NestedCachedPackedTraversalSupportFinder&>(snarl_caller.get_support_finder())){
  
    for (int i = 0; i < ref_paths.size(); ++i) {
        ref_offsets[ref_paths[i]] = i < ref_path_offsets.size() ? ref_path_offsets[i] : 0;
        ref_path_set.insert(ref_paths[i]);
        ref_ploidies[ref_paths[i]] = i < ref_path_ploidies.size() ? ref_path_ploidies[i] : 2;
    }

}
   
NestedFlowCaller::~NestedFlowCaller() {

}

bool NestedFlowCaller::call_snarl(const Snarl& managed_snarl) {
    
    // remember the calls for each child snarl in this table
    CallTable call_table;

    bool called = call_snarl_recursive(managed_snarl, -1, "", make_pair(0, 0), call_table);

    if (called) { 
        emit_snarl_recursive(managed_snarl, -1, call_table);
    }

    return called;
}

bool NestedFlowCaller::call_snarl_recursive(const Snarl& managed_snarl, int max_ploidy,
                                            const string& parent_ref_path_name, pair<size_t, size_t> parent_ref_interval,
                                            CallTable& call_table) {

    // todo: In order to experiment with merging consecutive snarls to make longer traversals,
    // I am experimenting with sending "fake" snarls through this code.  So make a local
    // copy to work on to do things like flip -- calling any snarl_manager code that
    // wants a pointer will crash. 
    Snarl snarl = managed_snarl;

    // hook into our table entry
    CallRecord& record = call_table[managed_snarl];
    
    // get some reference information if possible
    // todo: make a function
    
    handle_t start_handle = graph.get_handle(snarl.start().node_id(), snarl.start().backward());
    handle_t end_handle = graph.get_handle(snarl.end().node_id(), snarl.end().backward());

    // as we're writing to VCF, we need a reference path through the snarl.  we
    // look it up directly from the graph, and abort if we can't find one
    set<string> start_path_names;
    graph.for_each_step_on_handle(start_handle, [&](step_handle_t step_handle) {
            string name = graph.get_path_name(graph.get_path_handle_of_step(step_handle));
            if (!Paths::is_alt(name) && (ref_path_set.empty() || ref_path_set.count(name))) {
                start_path_names.insert(name);
            }
            return true;
        });
    
    set<string> end_path_names;
    if (!start_path_names.empty()) {
        graph.for_each_step_on_handle(end_handle, [&](step_handle_t step_handle) {
                string name = graph.get_path_name(graph.get_path_handle_of_step(step_handle));
                if (!Paths::is_alt(name) && (ref_path_set.empty() || ref_path_set.count(name))) {                
                    end_path_names.insert(name);
                }
                return true;
            });
    }
    
    // we do the full intersection (instead of more quickly finding the first common path)
    // so that we always take the lexicographically lowest path, rather than depending
    // on the order of iteration which could change between implementations / runs.
    vector<string> common_names;
    std::set_intersection(start_path_names.begin(), start_path_names.end(),
                          end_path_names.begin(), end_path_names.end(),
                          std::back_inserter(common_names));

    string ref_path_name;
    SnarlTraversal ref_trav;
    int ref_trav_idx = -1;
    tuple<int64_t, int64_t, bool, step_handle_t, step_handle_t> ref_interval;
    string gt_ref_path_name;
    pair<size_t, size_t> gt_ref_interval;
    
    if (!common_names.empty()) {
        ref_path_name = common_names.front();

        // find the reference traversal and coordinates using the path position graph interface
        ref_interval = get_ref_interval(graph, snarl, ref_path_name);
        if (get<0>(ref_interval) == -1) {
            // no reference path found due to orientation conflict
            return false;
        }
        if (get<2>(ref_interval) == true) {
            // calling code assumes snarl forward on reference
            flip_snarl(snarl);
            ref_interval = get_ref_interval(graph, snarl, ref_path_name);
        }

        step_handle_t cur_step = get<3>(ref_interval);
        step_handle_t last_step = get<4>(ref_interval);
        if (get<2>(ref_interval)) {
            std::swap(cur_step, last_step);
        }
        bool start_backwards = snarl.start().backward() != graph.get_is_reverse(graph.get_handle_of_step(cur_step));
    
        while (true) {
            handle_t cur_handle = graph.get_handle_of_step(cur_step);
            Visit* visit = ref_trav.add_visit();
            visit->set_node_id(graph.get_id(cur_handle));
            visit->set_backward(start_backwards ? !graph.get_is_reverse(cur_handle) : graph.get_is_reverse(cur_handle));
            if (graph.get_id(cur_handle) == snarl.end().node_id()) {
                break;
            } else if (get<2>(ref_interval) == true) {
                if (!graph.has_previous_step(cur_step)) {
                    cerr << "Warning [vg call]: Unable, due to bug or corrupt path information, to trace reference path through snarl " << pb2json(snarl) << endl;
                    return false;
                }
                cur_step = graph.get_previous_step(cur_step);
            } else {
                if (!graph.has_next_step(cur_step)) {
                    cerr << "Warning [vg call]: Unable, due to bug or corrupt path information, to trace reference path through snarl " << pb2json(snarl) << endl;
                    return false;
                }
                cur_step = graph.get_next_step(cur_step);
            }
            // todo: we can compute flow at the same time
        }
        assert(ref_trav.visit(0) == snarl.start() && ref_trav.visit(ref_trav.visit_size() - 1) == snarl.end());

        gt_ref_path_name = ref_path_name;
        gt_ref_interval = make_pair(get<0>(ref_interval), get<1>(ref_interval));
        if (max_ploidy == -1) {
            max_ploidy = ref_ploidies[ref_path_name];
        }        
    } else {
        // if we have no reference infromation, try to get it from the parent snarl
        gt_ref_path_name = parent_ref_path_name;
        gt_ref_interval = parent_ref_interval;
        if (gt_ref_path_name.empty()) {
            // there's just no reference path through this snarl
            return false;
        }
        assert(max_ploidy >= 0);
    }
    
    // recurse on the children
    // todo: do we need to make this iterative for deep snarl trees? 
    const vector<const Snarl*>& children = snarl_manager.children_of(&managed_snarl);

    for (const Snarl* child : children) {
        if (!snarl_manager.is_trivial(child, graph)) {
            bool called = call_snarl_recursive(*child, max_ploidy, gt_ref_path_name, gt_ref_interval, call_table);
            if (!called) {
                return false;
            }
        }
    }

#ifdef debug
    cerr << "recursively calling " << pb2json(managed_snarl) << " with " << children.size() << " children"
         << " and ref_path " << gt_ref_path_name << " and parent ref_path " << parent_ref_path_name << endl << endl;
#endif

    // abstract away the child snarls in the graph.  traversals will bypass them via
    // "virtual" edges
    SnarlGraph snarl_graph(&graph, snarl_manager, children);

    if (snarl.start().node_id() == snarl.end().node_id() ||
        !graph.has_node(snarl.start().node_id()) || !graph.has_node(snarl.end().node_id())) {
        // can't call one-node or out-of graph snarls.
        return false;
    }
    // toggle average flow / flow width based on snarl length.  this is a bit inconsistent with
    // downstream which uses the longest traversal length, but it's a bit chicken and egg
    // todo: maybe use snarl length for everything?
    const auto& support_finder = dynamic_cast<SupportBasedSnarlCaller&>(snarl_caller).get_support_finder();
    
    bool greedy_avg_flow = false;
    {
        auto snarl_contents = snarl_manager.shallow_contents(&snarl, graph, false);
        if (max(snarl_contents.first.size(), snarl_contents.second.size()) > max_snarl_shallow_size) {
            return false;
        }
        size_t len_threshold = support_finder.get_average_traversal_support_switch_threshold();
        size_t length = 0;
        for (auto i = snarl_contents.first.begin(); i != snarl_contents.first.end() && length < len_threshold; ++i) {
            length += graph.get_length(graph.get_handle(*i));
        }
        greedy_avg_flow = length > len_threshold;
    }

    vector<SnarlTraversal> travs;
    FlowTraversalFinder* flow_trav_finder = dynamic_cast<FlowTraversalFinder*>(&traversal_finder);
    if (flow_trav_finder != nullptr) {
        // find the max flow traversals using specialized interface that accepts avg heurstic toggle and overlay
        pair<vector<SnarlTraversal>, vector<double>> weighted_travs = flow_trav_finder->find_weighted_traversals(snarl, greedy_avg_flow, &snarl_graph);
        travs = std::move(weighted_travs.first);
           
    } else {
        // find the traversals using the generic interface
        assert(false);
        travs = traversal_finder.find_traversals(snarl);
    }

    // todo: we need to make reference traversal nesting aware
#ifdef debug
    for (int i = 0; i < travs.size(); ++i) {
        cerr << "[" << i << "]: " << pb2json(travs[i]) << endl;
    }
#endif
    
    // find the reference traversal in the list of results from the traversal finder
    if (!ref_path_name.empty()) {
        for (int i = 0; i < travs.size() && ref_trav_idx < 0; ++i) {
            // todo: is there a way to speed this up?
            if (travs[i] == ref_trav) {
                ref_trav_idx = i;
            }
        }

        if (ref_trav_idx == -1) {
            ref_trav_idx = travs.size();
            // we didn't get the reference traversal from the finder, so we add it here
            travs.push_back(ref_trav);
#ifdef debug
            cerr << "[ref]: " << pb2json(ref_trav) << endl;
#endif
        }
    }
    // store the reference traversal information, which could be empty
    record.ref_path_name = ref_path_name;
    record.ref_trav_idx = ref_trav_idx;

    // in the snarl graph, snarls a represented by a snarl end point and that's it.  here we fix up the traversals
    // to actually embed the snarls
    // todo: should be able to avoid copy here!
    vector<SnarlTraversal> embedded_travs = travs;
    for (int i = 0; i < embedded_travs.size(); ++i) {
        SnarlTraversal& traversal = embedded_travs[i];
        if (i != ref_trav_idx) {
            snarl_graph.embed_snarls(traversal);
        } else {
            snarl_graph.embed_ref_path_snarls(traversal);
        }
    }

    bool ret_val = true;

    if (traversals_only) {
        assert(gaf_output);
        for (SnarlTraversal& traversal : travs) {
            snarl_graph.embed_snarls(traversal);
        }
        pair<string, int64_t> pos_info = get_ref_position(graph, snarl, ref_path_name, 0);
        emit_gaf_traversals(graph, print_snarl(snarl), travs, ref_trav_idx, pos_info.first, pos_info.second, &support_finder);
    } else {
        // use our support caller to choose our genotype
        for (int ploidy = 1; ploidy <= max_ploidy; ++ploidy) {
            vector<int> trav_genotype;
            unique_ptr<SnarlCaller::CallInfo> trav_call_info;
            std::tie(trav_genotype, trav_call_info) = snarl_caller.genotype(snarl, travs, ref_trav_idx, ploidy, gt_ref_path_name,  gt_ref_interval);
            assert(trav_genotype.empty() || trav_genotype.size() == ploidy);

            // update the traversal finder with summary support statistics from this call
            // todo: be smarted about ploidy here
            NestedCachedPackedTraversalSupportFinder::SupportMap& child_support_map = nested_support_finder.child_support_map;
            // todo: re-use information that was produced in genotype!!
            int max_trav_size = 0;
            vector<Support> genotype_supports = nested_support_finder.get_traversal_genotype_support(embedded_travs, trav_genotype, {}, ref_trav_idx, &max_trav_size);
            Support total_site_support = std::accumulate(genotype_supports.begin(), genotype_supports.end(), Support());
            // todo: do we want to use max_trav_size, or something derived from the genotype? 
            child_support_map[snarl] = make_tuple(total_site_support, total_site_support, max_trav_size);

            // and now we need to update our own table with the genotype            
            if (record.genotype_by_ploidy.size() < ploidy) {
                record.genotype_by_ploidy.resize(ploidy);
            }
            record.genotype_by_ploidy[ploidy-1].first = trav_genotype;
            record.genotype_by_ploidy[ploidy-1].second.reset(trav_call_info.release());
            record.travs = embedded_travs;
        
            ret_val = trav_genotype.size() == ploidy;
        }
    }

    return ret_val;
}

bool NestedFlowCaller::emit_snarl_recursive(const Snarl& managed_snarl, int ploidy, CallTable& call_table) {
    // fetch the current snarl from the table
    CallRecord& record = call_table[managed_snarl];

    // only emit snarl with reference backbone:
    // todo: emit when no call (at least optionally)
    if (record.ref_trav_idx >= 0 && !record.genotype_by_ploidy.empty() && ploidy != 0) {

        if (ploidy < 0) {
            ploidy = ref_ploidies[record.ref_path_name];
        }
        
        pair<vector<int>, unique_ptr<SnarlCaller::CallInfo>>& genotype = record.genotype_by_ploidy[ploidy - 1];

        // compute count how many times a nested snarl appears in the genotype.  this will be the ploidy
        // it gets emitted with
        // todo: feed into flatten_alt_allele!
        map<Snarl, int, NestedCachedPackedTraversalSupportFinder::snarl_less> nested_ploidy;
        for (int allele : genotype.first) {
            const SnarlTraversal& allele_trav = record.travs[allele];
            for (size_t i = 0; i < allele_trav.visit_size(); ++i) {
                const Visit& visit = allele_trav.visit(i);
                if (visit.node_id() == 0) {
                    ++nested_ploidy[visit.snarl()];
                }
            }
        }

        // recurse on the children
        // todo: do we need to make this iterative for deep snarl trees? 
        const vector<const Snarl*>& children = snarl_manager.children_of(&managed_snarl);
        
        for (const Snarl* child : children) {
            if (!snarl_manager.is_trivial(child, graph)) {
                emit_snarl_recursive(*child, nested_ploidy[*child], call_table);
            }
        }

#ifdef debug
        cerr << "Recursively emitting " << pb2json(managed_snarl) << "with ploidy " << ploidy << endl;
#endif
        function<string(const vector<SnarlTraversal>&, const vector<int>&, int, int, int)> trav_to_flat_string =
            [&](const vector<SnarlTraversal>& travs, const vector<int>& travs_genotype, int trav_allele, int genotype_allele, int ref_trav_idx) {

            string allele_string = trav_string(graph, travs[trav_allele]);
            if (trav_allele == ref_trav_idx) {
                return flatten_reference_allele(allele_string, call_table);
            } else {
                int allele_ploidy = std::max((int)std::count(travs_genotype.begin(), travs_genotype.end(), trav_allele), 1);
                return flatten_alt_allele(allele_string, std::min(allele_ploidy-1, genotype_allele), allele_ploidy, call_table);
            }
        };
                
        if (!gaf_output) {
            emit_variant(graph, snarl_caller, managed_snarl, record.travs, genotype.first, record.ref_trav_idx, genotype.second, record.ref_path_name,
                         ref_offsets[record.ref_path_name], genotype_snarls, ploidy, trav_to_flat_string);
        } else {
            // todo:
            //    emit_gaf_variant(graph, snarl, travs, trav_genotype);
        }
    }

    return true;
}

string NestedFlowCaller::flatten_reference_allele(const string& nested_allele, const CallTable& call_table) const {

    string flat_allele;

    scan_snarl(nested_allele, [&](const string& fragment, Snarl& snarl) {
            if (!fragment.empty()) {
                flat_allele += fragment;
            } else {
                const CallRecord& record = call_table.at(snarl);
                assert(record.ref_trav_idx >= 0);
                if (record.travs.empty()) {
                    flat_allele += "<***>";
                    assert(false);
                } else{
                    const SnarlTraversal& traversal = record.travs[record.ref_trav_idx];
                    string nested_snarl_allele = trav_string(graph, traversal);
                    flat_allele += flatten_reference_allele(nested_snarl_allele, call_table);
                }
            }       
        });
    
    return flat_allele;
}

string NestedFlowCaller::flatten_alt_allele(const string& nested_allele, int allele, int ploidy, const CallTable& call_table) const {

    string flat_allele;
#ifdef debug
    cerr << "Flattening " << nested_allele << " at allele " << allele << endl;
#endif
    scan_snarl(nested_allele, [&](const string& fragment, Snarl& snarl) {
            if (!fragment.empty()) {
                flat_allele += fragment;
            } else {
                const CallRecord& record = call_table.at(snarl);
#ifdef debug
                cerr << "got record with " << record.travs.size() << " travs and " << record.genotype_by_ploidy.size() << " gts" << endl;
#endif
                int fallback_allele = -1;
                if (record.genotype_by_ploidy[ploidy-1].first.empty()) {
                    // there's no call here. but we really want to emit something, so try picking the reference
                    // or first allele
                    if (record.ref_trav_idx >= 0) {
                        fallback_allele = record.ref_trav_idx;
                    } else if (!record.travs.empty()) {
                        fallback_allele = 0;
                    }
                }
                if (fallback_allele >= (int)record.travs.size()) {
                    flat_allele += "<...>";
                } else {
                    // todo: passing in a single ploidy simplisitic, would need to derive from the calls when
                    // reucrising
                    // in practice, the results will nearly the same but still needs fixing
                    // we try to get the the allele from the genotype if possible, but fallback on the fallback_allele
                    int trav_allele = fallback_allele >= 0 ? fallback_allele : record.genotype_by_ploidy[ploidy-1].first[allele];
                    const SnarlTraversal& traversal = record.travs[trav_allele];
                    string nested_snarl_allele = trav_string(graph, traversal);
                    flat_allele += flatten_alt_allele(nested_snarl_allele, allele, ploidy, call_table);
                }                
            }  
        });
    
    return flat_allele;
}



string NestedFlowCaller::vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                              const vector<size_t>& contig_length_overrides) const {
    string header = VCFOutputCaller::vcf_header(graph, contigs, contig_length_overrides);
    header += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    snarl_caller.update_vcf_header(header);
    header += "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
    header += "##SAMPLE=<ID=" + sample_name + ">\n";
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_name;
    assert(output_vcf.openForOutput(header));
    header += "\n";
    return header;
}


SnarlGraph::SnarlGraph(const HandleGraph* backing_graph, SnarlManager& snarl_manager, vector<const Snarl*> snarls) :
    backing_graph(backing_graph),
    snarl_manager(snarl_manager) {
    for (const Snarl* snarl : snarls) {
        if (!snarl_manager.is_trivial(snarl, *backing_graph)) {
            this->snarls[backing_graph->get_handle(snarl->start().node_id(), snarl->start().backward())] =
                make_pair(backing_graph->get_handle(snarl->end().node_id(), snarl->end().backward()), true);
            this->snarls[backing_graph->get_handle(snarl->end().node_id(), !snarl->end().backward())] =
                make_pair(backing_graph->get_handle(snarl->start().node_id(), !snarl->start().backward()), false);
        }
    }
}

pair<bool, handle_t> SnarlGraph::node_to_snarl(handle_t handle) const {
    auto i = snarls.find(handle);
    if (i != snarls.end()) {
        return make_pair(true, i->second.first);
    } else {
        return make_pair(false, handle);
    }
}

tuple<bool, handle_t, edge_t> SnarlGraph::edge_to_snarl_edge(edge_t edge) const {
    auto i = snarls.find(edge.first);
    edge_t out_edge;
    handle_t out_node;
    bool out_found = false;
    if (i != snarls.end()) {
        // edge is from snarl start to after snarl end
        out_edge.first = i->second.first;
        out_edge.second = edge.second;
        out_node = edge.first;
        out_found = true;
    } else {
        // reverse of above
        i = snarls.find(backing_graph->flip(edge.second));
        if (i != snarls.end()) {
            out_edge.first = edge.first;
            out_edge.second = backing_graph->flip(i->second.first);
            out_node = edge.second;
            out_found = true;
        }
    }
    // note that we only have those two cases since our internal map contains
    // both orientations of the snarl.

    return make_tuple(out_found, out_node, out_edge);
}

void SnarlGraph::embed_snarl(Visit& visit) {
    handle_t handle = backing_graph->get_handle(visit.node_id(), visit.backward());
    auto it = snarls.find(handle);
    if (it != snarls.end()) {
        // edit the Visit in place to replace id, with the full snarl
        Snarl* snarl = visit.mutable_snarl();
        snarl->mutable_start()->set_node_id(visit.node_id());
        snarl->mutable_start()->set_backward(visit.backward());
        handle_t other = it->second.first;
        snarl->mutable_end()->set_node_id(backing_graph->get_id(other));
        snarl->mutable_end()->set_backward(backing_graph->get_is_reverse(other));
        if (it->second.second == false) {
            // put the snarl in an orientation consisten with other indexes
            swap(*snarl->mutable_start(), *snarl->mutable_end());
            snarl->mutable_start()->set_backward(!snarl->start().backward());
            snarl->mutable_end()->set_backward(!snarl->end().backward());
        }
        visit.set_node_id(0);
    }    
}

void SnarlGraph::embed_snarls(SnarlTraversal& traversal) {
    for (size_t i = 0; i < traversal.visit_size(); ++i) {
        Visit& visit = *traversal.mutable_visit(i);
        if (visit.node_id() > 0) { 
            embed_snarl(visit);
        }
    }
}

void SnarlGraph::embed_ref_path_snarls(SnarlTraversal& traversal) {
    vector<Visit> out_trav;
    size_t snarl_count = 0;
    bool in_snarl = false;
    handle_t snarl_end;
    for (size_t i = 0; i < traversal.visit_size(); ++i) {
        Visit& visit = *traversal.mutable_visit(i);
        handle_t handle = backing_graph->get_handle(visit.node_id(), visit.backward());
        if (in_snarl) {
            // nothing to do if we're in a snarl except check for the end and come out
            if (handle == snarl_end) {
                in_snarl = false;
            } 
        } else {
            // if we're not in a snarl, check for a new one
            auto it = snarls.find(handle);
            if (it != snarls.end()) {
                embed_snarl(visit);
                snarl_end = it->second.first;
                in_snarl = true;
                ++snarl_count;
            }
            out_trav.push_back(visit);
        }
    }

    // switch in the updated traversal
    if (snarl_count > 0) {
        traversal.clear_visit();
        for (Visit& visit : out_trav) {
            *traversal.add_visit() = visit;
        }
    }
}

bool SnarlGraph::follow_edges_impl(const handle_t& handle, bool go_left, const std::function<bool(const handle_t&)>& iteratee) const {
    if (!go_left) {        
        auto i = snarls.find(handle);
        if (i == snarls.end()) {
            return backing_graph->follow_edges(handle, go_left, iteratee);
        } else {
            return backing_graph->follow_edges(i->second.first, go_left, iteratee);
        }
    } else {
        return this->follow_edges_impl(backing_graph->flip(handle), !go_left, iteratee);
    }
}

// a lot of these don't strictly make sense.  ex, we would want has_node to
// hide stuff inside snarls.  but... we don't want to pay the cost of maintining
// structures for functions that aren't used..
bool SnarlGraph::has_node(nid_t node_id) const {
    return backing_graph->has_node(node_id);
}
handle_t SnarlGraph::get_handle(const nid_t& node_id, bool is_reverse) const {
    return backing_graph->get_handle(node_id, is_reverse);
}
nid_t SnarlGraph::get_id(const handle_t& handle) const {
    return backing_graph->get_id(handle);
}
bool SnarlGraph::get_is_reverse(const handle_t& handle) const {
    return backing_graph->get_is_reverse(handle);
}
handle_t SnarlGraph::flip(const handle_t& handle) const {
    return backing_graph->flip(handle);
}
size_t SnarlGraph::get_length(const handle_t& handle) const {
    return backing_graph->get_length(handle);
}
std::string SnarlGraph::get_sequence(const handle_t& handle) const {
    return backing_graph->get_sequence(handle);
}
size_t SnarlGraph::get_node_count() const {
    return backing_graph->get_node_count();
}
nid_t SnarlGraph::min_node_id() const {
    return backing_graph->min_node_id();
}
nid_t SnarlGraph::max_node_id() const {
    return backing_graph->max_node_id();
}
bool SnarlGraph::for_each_handle_impl(const std::function<bool(const handle_t&)>& iteratee, bool parallel) const {
    return backing_graph->for_each_handle(iteratee, parallel);
}

}

