/**
 * \file tee_alignment_emitter.cpp
 * Implementation for TeeAlignmentEmitter
 */

#include "tee_alignment_emitter.hpp"

#include <stdexcept>

namespace vg {

TeeAlignmentEmitter::TeeAlignmentEmitter(std::vector<std::unique_ptr<vg::io::AlignmentEmitter>>&& destinations) :
    destinations(std::move(destinations)) {
    for (auto& destination : this->destinations) {
        if (!destination) {
            throw std::invalid_argument("TeeAlignmentEmitter cannot emit through a null emitter");
        }
    }
}

void TeeAlignmentEmitter::emit_singles(std::vector<Alignment>&& aln_batch) {
    if (destinations.empty()) {
        return;
    }
    for (size_t i = 0; i + 1 < destinations.size(); i++) {
        destinations[i]->emit_singles(std::vector<Alignment>(aln_batch));
    }
    destinations.back()->emit_singles(std::move(aln_batch));
}

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
