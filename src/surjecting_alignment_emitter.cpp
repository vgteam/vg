/**
 * \file surjecting_alignment_emitter.cpp
 * Implementation for SurjectingAlignmentEmitter
 */


#include "surjecting_alignment_emitter.hpp"
#include "hts_alignment_emitter.hpp"

#include <map>

namespace vg {

using namespace std;

unique_ptr<AlignmentEmitter> get_alignment_emitter_with_surjection(const string& filename, const string& format, 
                                                                   const vector<path_handle_t> paths_in_dict_order, size_t max_threads,
                                                                   const HandleGraph* graph) {
    
    // Build a path length map by name
    // TODO: be able to pass along dict order!
    map<string, int64_t> path_length;
    if (format == "SAM" || format == "BAM" || format == "CRAM") {
        // Make sure we actually have a PathPositionalHandleGraph
        const PathPositionHandleGraph* path_graph = dynamic_cast<const PathPositionHandleGraph*>(graph);
        if (path_graph == nullptr) {
            cerr << "error[vg::get_alignment_emitter_with_surjection]: Graph does not contain paths to surject into." << endl;
            exit(1);
        }
    
        // We will actually use it, so fill it in.
        for (auto& path_handle : paths_in_dict_order) {
            path_length[path_graph->get_path_name(path_handle)] = path_graph->get_path_length(path_handle);
        }
    }
    
    // Get the non-surjecting emitter
    unique_ptr<AlignmentEmitter> emitter = get_alignment_emitter(filename, format, path_length, max_threads, graph);
    
    if (format == "SAM" || format == "BAM" || format == "CRAM") {
        // Need to surject
        
        // Make sure (again) we actually have a PathPositionalHandleGraph
        const PathPositionHandleGraph* path_graph = dynamic_cast<const PathPositionHandleGraph*>(graph);
        assert(path_graph != nullptr);
        
        // Make a set of the path handles to surject into
        unordered_set<path_handle_t> target_paths(paths_in_dict_order.begin(), paths_in_dict_order.end());
        // Interpose a surjecting AlignmentEmitter
        emitter = make_unique<SurjectingAlignmentEmitter>(path_graph, target_paths, std::move(emitter));
    }
    
    return emitter;
}

SurjectingAlignmentEmitter::SurjectingAlignmentEmitter(const PathPositionHandleGraph* graph, unordered_set<path_handle_t> paths,
    unique_ptr<AlignmentEmitter>&& backing) : surjector(graph), paths(paths), backing(std::move(backing)) {
    
    // Nothing to do!
    
}

void SurjectingAlignmentEmitter::surject_alignments_in_place(vector<Alignment>& alns) const {
    for (auto& aln : alns) {
        // Surject each alignment and annotate with surjected path position
        aln = surjector.surject(aln, paths, surject_subpath_global);
    }
}

void SurjectingAlignmentEmitter::emit_singles(vector<Alignment>&& aln_batch) {
    // Intercept the batch on its way
    vector<Alignment> aln_batch_caught(aln_batch);
    // Surject it in place
    surject_alignments_in_place(aln_batch_caught);
    // Forward it along
    backing->emit_singles(std::move(aln_batch_caught));
}

void SurjectingAlignmentEmitter::emit_mapped_singles(vector<vector<Alignment>>&& alns_batch) {
    // Intercept the batch on its way
    vector<vector<Alignment>> alns_batch_caught(alns_batch);
    for (auto& mappings : alns_batch_caught) {
        // Surject all mappings in place
        surject_alignments_in_place(mappings);
    }
    // Forward it along
    backing->emit_mapped_singles(std::move(alns_batch_caught));
}

void SurjectingAlignmentEmitter::emit_pairs(vector<Alignment>&& aln1_batch, vector<Alignment>&& aln2_batch, vector<int64_t>&& tlen_limit_batch) {
    // Intercept the batch on its way
    vector<Alignment> aln1_batch_caught(aln1_batch);
    vector<Alignment> aln2_batch_caught(aln2_batch);
    // Surject it in place
    surject_alignments_in_place(aln1_batch_caught);
    surject_alignments_in_place(aln2_batch_caught);
    // Forward it along
    backing->emit_pairs(std::move(aln1_batch_caught), std::move(aln2_batch_caught), std::move(tlen_limit_batch));
}

void SurjectingAlignmentEmitter::emit_mapped_pairs(vector<vector<Alignment>>&& alns1_batch, vector<vector<Alignment>>&& alns2_batch, vector<int64_t>&& tlen_limit_batch) {
    // Intercept the batch on its way
    vector<vector<Alignment>> alns1_batch_caught(alns1_batch);
    vector<vector<Alignment>> alns2_batch_caught(alns2_batch);
    for (auto& mappings : alns1_batch_caught) {
        // Surject all mappings in place
        surject_alignments_in_place(mappings);
    }
    for (auto& mappings : alns2_batch_caught) {
        // Surject all mappings in place
        surject_alignments_in_place(mappings);
    }
    // Forward it along
    backing->emit_mapped_pairs(std::move(alns1_batch_caught), std::move(alns2_batch_caught), std::move(tlen_limit_batch));
}

}
