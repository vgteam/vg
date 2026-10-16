#include <atomic>
#include <charconv>
#include <chrono>
#include <cstdio>
#include <limits>

#include <omp.h>

#include "gaf_output_caller.hpp"
#include "symbolic_allele.hpp"
#include "read_likelihood_caller.hpp"
#include "algorithms/expand_context.hpp"
#include "annotation.hpp"
#include "gref.hpp"
#include "traversal_clusters.hpp"
#include "utility.hpp"

//#define debug

namespace vg {

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
    if (ref_trav_idx >= 0 && ref_trav_idx < (int64_t)travs.size()) {
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
    // gt_travs is a NEW vector holding one entry per called allele, so ref_trav_idx -- an index into
    // travs -- does not address it.  Passing it through unremapped made emit_gaf_traversals index a
    // 2-element vector with an index into the full traversal list.  The genotype can also carry the
    // negative star/missing markers, which address nothing at all.
    vector<SnarlTraversal> gt_travs;
    int64_t gt_ref_trav_idx = -1;
    for (int allele : genotype) {
        if (allele < 0 || allele >= (int)travs.size()) {
            continue;
        }
        if (allele == ref_trav_idx && gt_ref_trav_idx < 0) {
            gt_ref_trav_idx = (int64_t)gt_travs.size();
        }
        gt_travs.push_back(travs[allele]);
    }
    if (gt_travs.empty()) {
        return;
    }
    emit_gaf_traversals(graph, snarl_name, gt_travs, gt_ref_trav_idx, ref_path_name, ref_path_position, support_finder);
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
            // path and handle oriented the same, we can just backtrack along the path to get
            // previous stuff
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
}

