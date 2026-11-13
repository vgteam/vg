#ifndef VG_TEE_ALIGNMENT_EMITTER_HPP_INCLUDED
#define VG_TEE_ALIGNMENT_EMITTER_HPP_INCLUDED

/** \file
 *
 * Holds an AlignmentEmitter that writes every alignment to two outputs.
 */

#include "vg/io/alignment_emitter.hpp"

#include <memory>
#include <string>
#include <vector>

namespace vg {

/**
 * A TeeAlignmentEmitter is an AlignmentEmitter that hands every batch of
 * alignments to two emitters it owns: a primary one and a copy one.
 *
 * Both emitters receive the same alignments in the same batches, so each one
 * writes what it would have written if it had been the only emitter. This lets
 * one mapping run write, for example, surjected SAM/BAM/CRAM through the
 * primary emitter and the graph alignments as GAM or GAF through the copy
 * emitter, instead of mapping the reads twice.
 *
 * The primary emitter gets a copy of each batch, and the copy emitter gets the
 * batch itself, since an emitter may move from or change what it is given.
 * Like any AlignmentEmitter, a TeeAlignmentEmitter may be called from several
 * threads at once, and relies on the two emitters it owns to be safe for that.
 */
class TeeAlignmentEmitter : public vg::io::AlignmentEmitter {
public:

    /// Make an emitter that sends every batch to both of the given emitters,
    /// and takes ownership of them. Neither may be null.
    TeeAlignmentEmitter(std::unique_ptr<vg::io::AlignmentEmitter>&& primary,
                        std::unique_ptr<vg::io::AlignmentEmitter>&& copy);

    /// Emit a batch of Alignments to both emitters.
    void emit_singles(std::vector<Alignment>&& aln_batch) override;
    /// Emit a batch of Alignments with secondaries to both emitters.
    void emit_mapped_singles(std::vector<std::vector<Alignment>>&& alns_batch) override;
    /// Emit a batch of pairs of Alignments, with the maximum pairing distance
    /// of each pair, to both emitters.
    void emit_pairs(std::vector<Alignment>&& aln1_batch, std::vector<Alignment>&& aln2_batch,
                    std::vector<int64_t>&& tlen_limit_batch) override;
    /// Emit the mappings of a batch of pairs of Alignments, with the maximum
    /// pairing distance of each pair, to both emitters.
    void emit_mapped_pairs(std::vector<std::vector<Alignment>>&& alns1_batch,
                           std::vector<std::vector<Alignment>>&& alns2_batch,
                           std::vector<int64_t>&& tlen_limit_batch) override;
    /// Emit some extra type-tagged data to both emitters, each of which keeps
    /// it if its format supports it.
    void emit_extra_message(const std::string& tag, std::string&& data) override;

private:
    /// The emitter for the main output, which gets a copy of each batch.
    std::unique_ptr<vg::io::AlignmentEmitter> primary;
    /// The emitter for the second output, which gets each batch itself.
    std::unique_ptr<vg::io::AlignmentEmitter> copy;
};

}

#endif
