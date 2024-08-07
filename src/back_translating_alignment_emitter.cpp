/**
 * \file back_translating_alignment_emitter.cpp
 * Implementation for BackTranslatingAlignmentEmitter
 */


#include "back_translating_alignment_emitter.hpp"
#include "algorithms/back_translate.hpp"

namespace vg {

using namespace std;

BackTranslatingAlignmentEmitter::BackTranslatingAlignmentEmitter(const NamedNodeBackTranslation* translation,
    unique_ptr<AlignmentEmitter>&& backing) : translation(translation), backing(std::move(backing)) {
    // Nothing to do!
}

void BackTranslatingAlignmentEmitter::back_translate_alignments_in_place(vector<Alignment>& alns) const {
    for (auto& aln : alns) {
        algorithms::back_translate_in_place(translation, *aln.mutable_path()); 
    }
}

void BackTranslatingAlignmentEmitter::emit_singles(vector<Alignment>&& aln_batch) {
    // Intercept the batch on its way
    vector<Alignment> aln_batch_caught(aln_batch);
    // Process it in place
    back_translate_alignments_in_place(aln_batch_caught);
    // Forward it along
    backing->emit_singles(std::move(aln_batch_caught));
}

void BackTranslatingAlignmentEmitter::emit_mapped_singles(vector<vector<Alignment>>&& alns_batch) {
    // Intercept the batch on its way
    vector<vector<Alignment>> alns_batch_caught(alns_batch);
    for (auto& mappings : alns_batch_caught) {
        // Surject all mappings in place
        back_translate_alignments_in_place(mappings);
    }
    // Forward it along
    backing->emit_mapped_singles(std::move(alns_batch_caught));
}

void BackTranslatingAlignmentEmitter::emit_pairs(vector<Alignment>&& aln1_batch, vector<Alignment>&& aln2_batch, vector<int64_t>&& tlen_limit_batch) {
    // Intercept the batch on its way
    vector<Alignment> aln1_batch_caught(aln1_batch);
    vector<Alignment> aln2_batch_caught(aln2_batch);
    // Process it in place
    back_translate_alignments_in_place(aln1_batch_caught);
    back_translate_alignments_in_place(aln2_batch_caught);
    // Forward it along
    backing->emit_pairs(std::move(aln1_batch_caught), std::move(aln2_batch_caught), std::move(tlen_limit_batch));
}

void BackTranslatingAlignmentEmitter::emit_mapped_pairs(vector<vector<Alignment>>&& alns1_batch, vector<vector<Alignment>>&& alns2_batch, vector<int64_t>&& tlen_limit_batch) {
    // Intercept the batch on its way
    vector<vector<Alignment>> alns1_batch_caught(alns1_batch);
    vector<vector<Alignment>> alns2_batch_caught(alns2_batch);
    for (auto& mappings : alns1_batch_caught) {
        // Process all mappings in place
        back_translate_alignments_in_place(mappings);
    }
    for (auto& mappings : alns2_batch_caught) {
        // Process all mappings in place
        back_translate_alignments_in_place(mappings);
    }
    // Forward it along
    backing->emit_mapped_pairs(std::move(alns1_batch_caught), std::move(alns2_batch_caught), std::move(tlen_limit_batch));
}

void BackTranslatingAlignmentEmitter::emit_extra_message(const std::string& tag, std::string&& data) {
    backing->emit_extra_message(tag, std::move(data));
}

}
