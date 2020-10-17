#ifndef VG_SURJECTING_ALIGNMENT_EMITTER_HPP_INCLUDED
#define VG_SURJECTING_ALIGNMENT_EMITTER_HPP_INCLUDED

/** \file
 *
 * Holds a surjecting wrapper AlignmentEmitter, and an entry point for getting AlignmentEmitters that can use it.
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
 * Get an AlignmentEmitter that emits in the given format.
 * Supported formats are: GAF, GAM, SAM, BAM, CRAM, JSON.
 * When emitting in linear formats (SAM/BAM/CRAM), automatically surjects, and uses spliced alignment output.
 * When surjecting, surjects into the paths in the given vector, with the order
 * used to define the order of names in the SAM/BAM/CRAM header.
 * When surjecting, graph must be a PahtPositionalHandleGraph.
 */
unique_ptr<vg::io::AlignmentEmitter> get_alignment_emitter_with_surjection(const string& filename, const string& format, 
                                                                           const vector<path_handle_t> paths_in_dict_order, size_t max_threads,
                                                                           const HandleGraph* graph);
                                                                           
/**
 * Produce a list of path handles in a fixed order, suitable for use with get_alignment_emitter_with_surjection(), by parsing a file.
 * The file may be an HTSlib-style "sequence dictionary" (consisting of SAM @SQ header lines), or a plain list of sequence names (which do not start with "@SQ").
 * If the file is not openable or contains no entries, reports an error and quits.
 * If the filename is itself an empty string, all non-alt-allele paths from the graph will be collected in arbitrary order.
 */
vector<path_handle_t> get_sequence_dictionary_handles(const string& filename, const PathPositionHandleGraph& graph);

/**
 * An AlignmentEmitter implementation that surjects alignments before emitting them via a backing AlignmentEmitter, which it owns.
 */
class SurjectingAlignmentEmitter : public vg::io::AlignmentEmitter {
public:
    
    /**
     * Surject alignments using the given graph, into the given paths, and send them to the given AlignmentEmitter.
     * Takes ownership of the AlignmentEmitter.
     * Copies the set of paths.
     */
    SurjectingAlignmentEmitter(const PathPositionHandleGraph* graph, unordered_set<path_handle_t> paths, unique_ptr<AlignmentEmitter>&& backing);
   
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
