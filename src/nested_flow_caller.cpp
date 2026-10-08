#include <atomic>
#include <charconv>
#include <chrono>
#include <cstdio>
#include <limits>

#include <omp.h>

#include "nested_flow_caller.hpp"
#include "symbolic_allele.hpp"
#include "read_likelihood_caller.hpp"
#include "algorithms/expand_context.hpp"
#include "annotation.hpp"
#include "gref.hpp"
#include "traversal_clusters.hpp"
#include "utility.hpp"

//#define debug

namespace vg {

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
        // Prefer base reference paths over derived gref paths.  Test the whole gref
        // namespace, not just the fragment suffix: a gref copy of the reference sorts
        // before the path it was copied from (gref_x < x).
        ref_path_name = common_names.front();  // default to first (lexicographically smallest)
        for (const string& name : common_names) {
            if (!GrefCover::is_gref_derived(name)) {
                ref_path_name = name;
                break;
            }
        }

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
            max_ploidy = ploidy_regions.ploidy_at(ref_path_name, get<0>(ref_interval),
                                                  ref_offset_of(ref_offsets, ref_path_name),
                                                  ref_ploidy_of(ref_ploidies, ref_path_name));
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
        // find the max flow traversals using specialized interface that accepts avg heurstic toggle
        // and overlay
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
    // The interval the ploidy above was decided over. emit_snarl_recursive derives the ploidy
    // again from it, and indexes genotype_by_ploidy by that ploidy.
    record.ref_path_interval = make_pair((int64_t)gt_ref_interval.first,
                                         (int64_t)gt_ref_interval.second);
    record.ref_trav_idx = ref_trav_idx;

    // in the snarl graph, snarls a represented by a snarl end point and that's it.  here we fix up
    // the traversals to actually embed the snarls todo: should be able to avoid copy here!
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
            // Must agree with the ploidy the genotype was decided at, since genotype_by_ploidy is
            // indexed by it, so the record's own interval is used.
            ploidy = ploidy_regions.ploidy_at(record.ref_path_name, record.ref_path_interval.first,
                                              ref_offset_of(ref_offsets, record.ref_path_name),
                                              ref_ploidy_of(ref_ploidies, record.ref_path_name));
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
            bool added = emit_variant(graph, snarl_caller, managed_snarl, record.travs, genotype.first, record.ref_trav_idx, genotype.second, record.ref_path_name,
                                 ref_offset_of(ref_offsets, record.ref_path_name), genotype_snarls, ploidy, trav_to_flat_string);
            if (!added) {
                return false;
            }
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
                    // there's no call here. but we really want to emit something, so try picking
                    // the reference or first allele
                    if (record.ref_trav_idx >= 0) {
                        fallback_allele = record.ref_trav_idx;
                    } else if (!record.travs.empty()) {
                        fallback_allele = 0;
                    }
                }
                if (fallback_allele >= (int)record.travs.size()) {
                    flat_allele += "<...>";
                } else {
                    // todo: passing in a single ploidy simplisitic, would need to derive from the
                    // calls when reucrising in practice, the results will nearly the same but still
                    // needs fixing we try to get the allele from the genotype if possible, but
                    // fallback on the fallback_allele
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
}

