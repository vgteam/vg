/**
 * \file tee_alignment_emitter.cpp
 * Implementation for TeeAlignmentEmitter
 */

#include "tee_alignment_emitter.hpp"

#include <stdexcept>

namespace vg {

TeeAlignmentEmitter::TeeAlignmentEmitter(std::vector<std::unique_ptr<vg::io::AlignmentEmitter>>&& destinations) :
    destinations(std::move(destinations)) {
    // TODO: Is it worth having code in the constructor *just* to check the invariant?
    for (auto& destination : this->destinations) {
        if (!destination) {
            throw std::invalid_argument("Downstream alignment emitters cannot be null");
        }
    }
}

void TeeAlignmentEmitter::emit_singles(std::vector<Alignment>&& aln_batch) {
    if (destinations.empty()) {
        return;
    }
    for (size_t i = 0; i + 1 < destinations.size(); i++) {
        // Emit copies to all but the last destination
        destinations[i]->emit_singles(std::vector<Alignment>(aln_batch));
    }
    // And hand the original off to the last destination
    destinations.back()->emit_singles(std::move(aln_batch));
}

// TODO: If we were in a fancy language with macros we could say "just do that
// to all the functions". But we aren't.

void TeeAlignmentEmitter::emit_mapped_singles(std::vector<std::vector<Alignment>>&& alns_batch) {
    if (destinations.empty()) {
        return;
    }
    for (size_t i = 0; i + 1 < destinations.size(); i++) {
        destinations[i]->emit_mapped_singles(std::vector<std::vector<Alignment>>(alns_batch));
    }
    destinations.back()->emit_mapped_singles(std::move(alns_batch));
}

void TeeAlignmentEmitter::emit_pairs(std::vector<Alignment>&& aln1_batch, std::vector<Alignment>&& aln2_batch,
                                     std::vector<int64_t>&& tlen_limit_batch) {
    if (destinations.empty()) {
        return;
    }
    for (size_t i = 0; i + 1 < destinations.size(); i++) {
        destinations[i]->emit_pairs(std::vector<Alignment>(aln1_batch),
                                    std::vector<Alignment>(aln2_batch),
                                    std::vector<int64_t>(tlen_limit_batch));
    }
    destinations.back()->emit_pairs(std::move(aln1_batch), std::move(aln2_batch), std::move(tlen_limit_batch));
}

void TeeAlignmentEmitter::emit_mapped_pairs(std::vector<std::vector<Alignment>>&& alns1_batch,
                                            std::vector<std::vector<Alignment>>&& alns2_batch,
                                            std::vector<int64_t>&& tlen_limit_batch) {
    if (destinations.empty()) {
        return;
    }
    for (size_t i = 0; i + 1 < destinations.size(); i++) {
        destinations[i]->emit_mapped_pairs(std::vector<std::vector<Alignment>>(alns1_batch),
                                           std::vector<std::vector<Alignment>>(alns2_batch),
                                           std::vector<int64_t>(tlen_limit_batch));
    }
    destinations.back()->emit_mapped_pairs(std::move(alns1_batch), std::move(alns2_batch),
                                           std::move(tlen_limit_batch));
}

void TeeAlignmentEmitter::emit_extra_message(const std::string& tag, std::string&& data) {
    if (destinations.empty()) {
        return;
    }
    for (size_t i = 0; i + 1 < destinations.size(); i++) {
        destinations[i]->emit_extra_message(tag, std::string(data));
    }
    destinations.back()->emit_extra_message(tag, std::move(data));
}

}
