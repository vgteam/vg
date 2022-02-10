#ifndef VG_BACK_TRANSLATING_ALIGNMENT_EMITTER_HPP_INCLUDED
#define VG_BACK_TRANSLATING_ALIGNMENT_EMITTER_HPP_INCLUDED

/** \file
 *
 * Holds a back-translating wrapper AlignmentEmitter.
 */


#include "vg/io/alignment_emitter.hpp"
#include "handle.hpp"

#include <unordered_set>
#include <string>
#include <vector>

namespace vg {

using namespace std;
                                                                           
/**
 * An AlignmentEmitter implementation that translates alignments into
 * named-segment space coordinates before emitting them via a backing
 * AlignmentEmitter, which it owns.
 */
class BackTranslatingAlignmentEmitter : public vg::io::AlignmentEmitter {
public:
    
    /**
     * Make an alignment emitter that translates alignments using the given
     * translation and emits them to the given backing AlignmentEmitter.
     * Takes ownership of the AlignmentEmitter.
     */
    BackTranslatingAlignmentEmitter(const NamedNodeBackTranslation* translation,
        unique_ptr<AlignmentEmitter>&& backing);
   
    /// Emit a batch of Alignments
    virtual void emit_singles(vector<Alignment>&& aln_batch);
    /// Emit batch of Alignments with secondaries. All secondaries must have is_secondary set already.
    virtual void emit_mapped_singles(vector<vector<Alignment>>&& alns_batch);
    /// Emit a batch of pairs of Alignments. The tlen_limit_batch, if
    /// specified, is the maximum pairing distance for ewch pair to flag
    /// properly paired, if the output format cares about such things. TODO:
    /// Move to a properly paired annotation that runs with the Alignment.
    virtual void emit_pairs(vector<Alignment>&& aln1_batch, vector<Alignment>&& aln2_batch,
        vector<int64_t>&& tlen_limit_batch);
    /// Emit the mappings of a batch of pairs of Alignments. All secondaries
    /// must have is_secondary set already. The tlen_limit_batch, if specified,
    /// is the maximum pairing distance for each pair to flag properly paired,
    /// if the output format cares about such things. TODO: Move to a properly
    /// paired annotation that runs with the Alignment.
    ///
    /// Both ends of each pair must have the same number of mappings.
    virtual void emit_mapped_pairs(vector<vector<Alignment>>&& alns1_batch,
        vector<vector<Alignment>>&& alns2_batch, vector<int64_t>&& tlen_limit_batch);
    
protected:
    /// Translation to use to translate node IDs to pieces of named segments.
    const NamedNodeBackTranslation* translation;

    /// AlignmentEmitter to emit to once done
    unique_ptr<AlignmentEmitter> backing;
    
    /// Back-translate alignments in place.
    void back_translate_alignments_in_place(vector<Alignment>& alns) const;
};

}

#endif
