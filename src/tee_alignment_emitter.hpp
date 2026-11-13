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
 * An AlignmentEmitter that duplicates alignments out to several other emitters.
 *
 * Note that, because different threads can call into the emit methods, the
 * different downstream emitters can see batches in different orders.
 */
class TeeAlignmentEmitter : public vg::io::AlignmentEmitter {
public:

    /// Make a TeeAlignmentEmitter that emits to the given destinations.
    ///
    /// An empty list of destinations is allowed.
    ///
    /// Throws std::invalid_argument if any destination is null.
    explicit TeeAlignmentEmitter(std::vector<std::unique_ptr<vg::io::AlignmentEmitter>>&& destinations);

    // AlignmentEmitter implementation methods
    
    void emit_singles(std::vector<Alignment>&& aln_batch) override;
    void emit_mapped_singles(std::vector<std::vector<Alignment>>&& alns_batch) override;
    void emit_pairs(std::vector<Alignment>&& aln1_batch, std::vector<Alignment>&& aln2_batch,
                    std::vector<int64_t>&& tlen_limit_batch) override;
    void emit_mapped_pairs(std::vector<std::vector<Alignment>>&& alns1_batch,
                           std::vector<std::vector<Alignment>>&& alns2_batch,
                           std::vector<int64_t>&& tlen_limit_batch) override;
    /// Pass some extra type-tagged data to every destination. Each one emits it
    /// only if its format supports it.
    void emit_extra_message(const std::string& tag, std::string&& data) override;

private:
    /// The downstream alignment emitters.
    std::vector<std::unique_ptr<vg::io::AlignmentEmitter>> destinations;
};

}

#endif
