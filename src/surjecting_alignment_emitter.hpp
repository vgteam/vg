#ifndef VG_SURJECTING_ALIGNMENT_EMITTER_HPP_INCLUDED
#define VG_SURJECTING_ALIGNMENT_EMITTER_HPP_INCLUDED

/** \file
 *
 * Holds a surjecting wrapper AlignmentEmitter.
 */


#include "surjector.hpp"
#include "vg/io/alignment_emitter.hpp"
#include "handle.hpp"

#include <unordered_set>
#include <string>
#include <vector>

namespace vg {

using namespace std;
                                                                           
/**
 * An AlignmentEmitter implementation that surjects alignments before emitting them via a backing AlignmentEmitter, which it owns.
 */
class SurjectingAlignmentEmitter : public vg::io::AlignmentEmitter {
public:
    
    /**
     * Surject alignments using the given graph, into the given paths, and send them to the given AlignmentEmitter.
     * Takes ownership of the AlignmentEmitter.
     * Copies the set of paths.
     *
     * If prune_suspicious_anchors is set, prunes out repetitive-looking
     * anchors when surjecting and lets those parts of reads be realigned.
     */
    SurjectingAlignmentEmitter(const PathPositionHandleGraph* graph,
        unordered_set<path_handle_t> paths, unique_ptr<AlignmentEmitter>&& backing,
        bool prune_suspicious_anchors = false);
   
    ///  Force full length alignment in surjection resolution 
    bool surject_subpath_global = true;
    
    
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

    /// Emit some extra type-tagged data, if the backing format supports it.
    virtual void emit_extra_message(const std::string& tag, std::string&& data);
    
protected:
    /// Surjector used to do the surjection
    Surjector surjector;

    /// Paths to surject into
    unordered_set<path_handle_t> paths;
    
    /// AlignmentEmitter to emit to once done
    unique_ptr<AlignmentEmitter> backing;
    
    /// Surject alignments in place.
    void surject_alignments_in_place(vector<Alignment>& alns) const;
    
    
    
    
};

}

#endif
