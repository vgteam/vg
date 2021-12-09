/**
 * \file surjecting_alignment_emitter.cpp
 * Implementation for SurjectingAlignmentEmitter
 */


#include "surjecting_alignment_emitter.hpp"
#include "hts_alignment_emitter.hpp"

#include <map>

namespace vg {

using namespace std;

SurjectingAlignmentEmitter::SurjectingAlignmentEmitter(const PathPositionHandleGraph* graph, unordered_set<path_handle_t> paths,
    unique_ptr<AlignmentEmitter>&& backing, bool prune_suspicious_anchors) : surjector(graph), paths(paths), backing(std::move(backing)) {
    
    // Configure the surjector
    surjector.prune_suspicious_anchors = prune_suspicious_anchors;
    
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
