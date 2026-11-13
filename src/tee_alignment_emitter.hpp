#ifndef VG_TEE_ALIGNMENT_EMITTER_HPP_INCLUDED
#define VG_TEE_ALIGNMENT_EMITTER_HPP_INCLUDED

/** \file
 *
 * Holds an AlignmentEmitter that emits through several other AlignmentEmitters.
 */

#include "vg/io/alignment_emitter.hpp"

#include <memory>
#include <string>
#include <vector>

namespace vg {

/**
 * An AlignmentEmitter that emits every alignment through each of several other
 * AlignmentEmitters, which it owns.
 *
 * Each destination gets the same batches in the same order, and writes what it
 * would write if it were used alone. Destinations may move from or change the
 * batches they get, so every destination but the last gets its own copy.
 *
 * It is thread safe if all its destinations are.
 */
class TeeAlignmentEmitter : public vg::io::AlignmentEmitter {
public:

    /// Make a TeeAlignmentEmitter that emits through the given destinations,
    /// in order. With no destinations, it discards everything.
    ///
    /// Throws std::invalid_argument if any destination is null.
    explicit TeeAlignmentEmitter(std::vector<std::unique_ptr<vg::io::AlignmentEmitter>>&& destinations);

    /// Emit a batch of Alignments through every destination.
    void emit_singles(std::vector<Alignment>&& aln_batch) override;
    /// Emit a batch of Alignments with secondaries through every destination.
    void emit_mapped_singles(std::vector<std::vector<Alignment>>&& alns_batch) override;
    /// Emit a batch of pairs of Alignments, with the maximum pairing distance
    /// of each pair, through every destination.
    void emit_pairs(std::vector<Alignment>&& aln1_batch, std::vector<Alignment>&& aln2_batch,
                    std::vector<int64_t>&& tlen_limit_batch) override;
    /// Emit the mappings of a batch of pairs of Alignments, with the maximum
    /// pairing distance of each pair, through every destination.
    void emit_mapped_pairs(std::vector<std::vector<Alignment>>&& alns1_batch,
                           std::vector<std::vector<Alignment>>&& alns2_batch,
                           std::vector<int64_t>&& tlen_limit_batch) override;
    /// Emit some extra type-tagged data through every destination whose format
    /// supports it.
    void emit_extra_message(const std::string& tag, std::string&& data) override;

private:
    /// The emitters that every batch goes to, in order. None is null.
    std::vector<std::unique_ptr<vg::io::AlignmentEmitter>> destinations;
};

}

#endif
