/**
 * \file tee_alignment_emitter.cpp
 * Implementation for TeeAlignmentEmitter
 */

#include "tee_alignment_emitter.hpp"

#include <stdexcept>

namespace vg {

TeeAlignmentEmitter::TeeAlignmentEmitter(std::unique_ptr<vg::io::AlignmentEmitter>&& primary,
                                         std::unique_ptr<vg::io::AlignmentEmitter>&& copy) :
    primary(std::move(primary)), copy(std::move(copy)) {
    if (!this->primary || !this->copy) {
        throw std::invalid_argument("TeeAlignmentEmitter needs two emitters");
    }
}

void TeeAlignmentEmitter::emit_singles(std::vector<Alignment>&& aln_batch) {
    primary->emit_singles(std::vector<Alignment>(aln_batch));
    copy->emit_singles(std::move(aln_batch));
}

void TeeAlignmentEmitter::emit_mapped_singles(std::vector<std::vector<Alignment>>&& alns_batch) {
    primary->emit_mapped_singles(std::vector<std::vector<Alignment>>(alns_batch));
    copy->emit_mapped_singles(std::move(alns_batch));
}

void TeeAlignmentEmitter::emit_pairs(std::vector<Alignment>&& aln1_batch, std::vector<Alignment>&& aln2_batch,
                                     std::vector<int64_t>&& tlen_limit_batch) {
    primary->emit_pairs(std::vector<Alignment>(aln1_batch), std::vector<Alignment>(aln2_batch),
                        std::vector<int64_t>(tlen_limit_batch));
    copy->emit_pairs(std::move(aln1_batch), std::move(aln2_batch), std::move(tlen_limit_batch));
}

void TeeAlignmentEmitter::emit_mapped_pairs(std::vector<std::vector<Alignment>>&& alns1_batch,
                                            std::vector<std::vector<Alignment>>&& alns2_batch,
                                            std::vector<int64_t>&& tlen_limit_batch) {
    primary->emit_mapped_pairs(std::vector<std::vector<Alignment>>(alns1_batch),
                               std::vector<std::vector<Alignment>>(alns2_batch),
                               std::vector<int64_t>(tlen_limit_batch));
    copy->emit_mapped_pairs(std::move(alns1_batch), std::move(alns2_batch), std::move(tlen_limit_batch));
}

void TeeAlignmentEmitter::emit_extra_message(const std::string& tag, std::string&& data) {
    primary->emit_extra_message(tag, std::string(data));
    copy->emit_extra_message(tag, std::move(data));
}

}
