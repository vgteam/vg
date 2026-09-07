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
    unique_ptr<AlignmentEmitter>&& backing, bool prune_suspicious_anchors, bool add_graph_alignment_tag, bool report_supplementary,
    bool add_off_ref_position_tag, bool left_align) : surjector(graph), paths(paths), backing(std::move(backing)) {
    
    // Configure the surjector
    surjector.prune_suspicious_anchors = prune_suspicious_anchors;
    surjector.annotate_with_graph_alignment = add_graph_alignment_tag;
    surjector.report_supplementary = report_supplementary;
    surjector.annotate_off_reference_pos = add_off_ref_position_tag;
    surjector.left_align = left_align;
}

void SurjectingAlignmentEmitter::emit_singles(vector<Alignment>&& aln_batch) {
    // Intercept the batch on its way
    vector<Alignment> aln_batch_caught(aln_batch);
    // Surject it in place
    surjector.surject_in_place(aln_batch_caught, paths, surject_subpath_global);
    // Forward it along
    backing->emit_singles(std::move(aln_batch_caught));
}

void SurjectingAlignmentEmitter::emit_mapped_singles(vector<vector<Alignment>>&& alns_batch) {
    // Intercept the batch on its way
    vector<vector<Alignment>> alns_batch_caught(alns_batch);
    for (auto& mappings : alns_batch_caught) {
        // Surject all mappings in place
        surjector.surject_in_place(mappings, paths, surject_subpath_global);
    }
    // Forward it along
    backing->emit_mapped_singles(std::move(alns_batch_caught));
}

void SurjectingAlignmentEmitter::emit_pairs(vector<Alignment>&& aln1_batch, vector<Alignment>&& aln2_batch, vector<int64_t>&& tlen_limit_batch) {
    // Intercept the batch on its way
    vector<Alignment> aln1_batch_caught(aln1_batch);
    vector<Alignment> aln2_batch_caught(aln2_batch);
    // Surject non-supplementary in place and gather supplementary
    vector<Alignment> supplementary_alns;
    surjector.surject_paired_in_place(aln1_batch_caught, aln2_batch_caught, supplementary_alns, supplementary_alns, paths, surject_subpath_global);
    // Forward them along
    backing->emit_pairs(std::move(aln1_batch_caught), std::move(aln2_batch_caught), std::move(tlen_limit_batch));
    backing->emit_singles(std::move(supplementary_alns));
}

void SurjectingAlignmentEmitter::emit_mapped_pairs(vector<vector<Alignment>>&& alns1_batch, vector<vector<Alignment>>&& alns2_batch, vector<int64_t>&& tlen_limit_batch) {
    // Intercept the batch on its way
    vector<vector<Alignment>> alns1_batch_caught(alns1_batch);
    vector<vector<Alignment>> alns2_batch_caught(alns2_batch);
    vector<vector<Alignment>> supplementary_batch;
    for (size_t i = 0; i < alns1_batch_caught.size(); ++i) {
        supplementary_batch.emplace_back();
        surjector.surject_paired_in_place(alns1_batch_caught[i], alns2_batch_caught[i], supplementary_batch.back(), supplementary_batch.back(), paths, surject_subpath_global);
    }
    // Forward it along
    backing->emit_mapped_pairs(std::move(alns1_batch_caught), std::move(alns2_batch_caught), std::move(tlen_limit_batch));
    backing->emit_mapped_singles(std::move(supplementary_batch));
}

void SurjectingAlignmentEmitter::emit_extra_message(const std::string& tag, std::string&& data) {
    backing->emit_extra_message(tag, std::move(data));
}

}
