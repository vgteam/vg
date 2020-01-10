#ifndef VG_ALIGNMENT_EMITTER_HPP_INCLUDED
#define VG_ALIGNMENT_EMITTER_HPP_INCLUDED

/**
 * \file alignment_emitter.hpp
 *
 * Defines a system for emitting alignments and groups of alignments in multiple formats.
 */

#include <mutex>
#include <thread>
#include <vector>
#include <deque>

#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

#include <vg/vg.pb.h>
#include <vg/io/protobuf_emitter.hpp>
#include <vg/io/stream_multiplexer.hpp>

#include "handle.hpp"

namespace vg {
using namespace std;

/**
 * Base class for a sink that takes alignments, possibly with pairing/secondary
 * relationships, and writes them out somewhere.
 *
 * All implementations must be thread safe.
 *
 * All implementations assume OMP threading.
 */
class AlignmentEmitter {
public:
    
    // These batched methods are the ones you need to implement. We batch for
    // efficiency. If there are any locks necessary to ensure thread safety,
    // you must lock once per batch instead of locking for every read.

    /// Emit a batch of Alignments
    virtual void emit_singles(vector<Alignment>&& aln_batch) = 0;
    /// Emit batch of Alignments with secondaries. All secondaries must have is_secondary set already.
    virtual void emit_mapped_singles(vector<vector<Alignment>>&& alns_batch) = 0;
    /// Emit a batch of pairs of Alignments. The tlen_limit_batch, if
    /// specified, is the maximum pairing distance for ewch pair to flag
    /// properly paired, if the output format cares about such things. TODO:
    /// Move to a properly paired annotation that runs with the Alignment.
    virtual void emit_pairs(vector<Alignment>&& aln1_batch, vector<Alignment>&& aln2_batch,
        vector<int64_t>&& tlen_limit_batch) = 0;
    /// Emit the mappings of a batch of pairs of Alignments. All secondaries
    /// must have is_secondary set already. The tlen_limit_batch, if specified,
    /// is the maximum pairing distance for each pair to flag properly paired,
    /// if the output format cares about such things. TODO: Move to a properly
    /// paired annotation that runs with the Alignment.
    ///
    /// Both ends of each pair must have the same number of mappings.
    virtual void emit_mapped_pairs(vector<vector<Alignment>>&& alns1_batch,
        vector<vector<Alignment>>&& alns2_batch, vector<int64_t>&& tlen_limit_batch) = 0;
    
    
    // These single-read methods have default implementations.
    
    /// Emit a single Alignment
    virtual void emit_single(Alignment&& aln);
    /// Emit a single Alignment with secondaries. All secondaries must have is_secondary set already.
    virtual void emit_mapped_single(vector<Alignment>&& alns);
    /// Emit a pair of Alignments. The tlen_limit, if specified, is the maximum
    /// pairing distance to flag properly paired, if the output format cares
    /// about such things. TODO: Move to a properly paired annotation that runs
    /// with the Alignment.
    virtual void emit_pair(Alignment&& aln1, Alignment&& aln2, int64_t tlen_limit = 0);
    /// Emit the mappings of a pair of Alignments. All secondaries must have
    /// is_secondary set already. The tlen_limit, if specified, is the maximum
    /// pairing distance to flag properly paired, if the output format cares
    /// about such things. TODO: Move to a properly paired annotation that runs
    /// with the Alignment.
    ///
    /// Both ends of the pair must have the same number of mappings.
    virtual void emit_mapped_pair(vector<Alignment>&& alns1, vector<Alignment>&& alns2,
        int64_t tlen_limit = 0);
    
    /// Allow destruction through base class pointer.
    virtual ~AlignmentEmitter() = default;
};

/// Get an AlignmentEmitter that can emit to the given file (or "-") in the
/// given format. A table of contig lengths is required for HTSlib formats.
/// Automatically applies per-thread buffering, but needs to know how many OMP
/// threads will be in use.
/// If the alignments are spliced at known splice sites (i.e. edges in the graph),
/// the graph can be provided in order to form spliced CIGAR strings for the HTSlib
/// formats. Other output formats ignore the graph.
unique_ptr<AlignmentEmitter> get_alignment_emitter(const string& filename, const string& format, 
    const map<string, int64_t>& path_length, size_t max_threads,
    const PathPositionHandleGraph* splicing_graph = nullptr);

/**
 * Discards all alignments.
 */
class NullAlignmentEmitter : public AlignmentEmitter {
public:
    inline virtual void emit_singles(vector<Alignment>&& aln_batch) {}
    inline virtual void emit_mapped_singles(vector<vector<Alignment>>&& alns_batch) {}
    inline virtual void emit_pairs(vector<Alignment>&& aln1_batch, vector<Alignment>&& aln2_batch,
        vector<int64_t>&& tlen_limit_batch) {}
    inline virtual void emit_mapped_pairs(vector<vector<Alignment>>&& alns1_batch,
        vector<vector<Alignment>>&& alns2_batch, vector<int64_t>&& tlen_limit_batch) {}
        
};

/**
 * Emit a TSV table describing alignments.
 */
class TSVAlignmentEmitter : public AlignmentEmitter {
public:
    
    /// Create a TSVAlignmentEmitter writing to the given file (or "-")
    TSVAlignmentEmitter(const string& filename, size_t max_threads);

    /// The default destructor should clean up the open file, if any.
    ~TSVAlignmentEmitter() = default;
    
    /// Emit a batch of Alignments.
    virtual void emit_singles(vector<Alignment>&& aln_batch);
    /// Emit a batch of Alignments with secondaries. All secondaries must have
    /// is_secondary set already.
    virtual void emit_mapped_singles(vector<vector<Alignment>>&& alns_batch);
    /// Emit a batch of pairs of Alignments.
    virtual void emit_pairs(vector<Alignment>&& aln1_batch, vector<Alignment>&& aln2_batch,
        vector<int64_t>&& tlen_limit_batch);
    /// Emit the mappings of a batch of pairs of Alignments. All secondaries
    /// must have is_secondary set already.
    virtual void emit_mapped_pairs(vector<vector<Alignment>>&& alns1_batch,
        vector<vector<Alignment>>&& alns2_batch, vector<int64_t>&& tlen_limit_batch);
    
private:

    /// If we are doing output to a file, this will hold the open file. Otherwise (for stdout) it will be empty.
    unique_ptr<ofstream> out_file;
    
    /// This holds a StreamMultiplexer on the output stream, for sharing it
    /// between threads.
    vg::io::StreamMultiplexer multiplexer;

    /// Emit single alignment as TSV.
    /// This is all we use; we don't do anything for pairing.
    void emit(Alignment&& aln_batch);
};

/**
 * Emit Alignments to a stream in SAM/BAM/CRAM format.
 * Thread safe.
 */
class HTSAlignmentEmitter : public AlignmentEmitter {
public:
    /// Create an HTSAlignmentEmitter writing to the given file (or "-") in the
    /// given HTS format ("SAM", "BAM", "CRAM"). path_length must map from
    /// contig name to length to include in the header. Sample names and read
    /// groups for the header will be guessed from the first reads. HTSlib
    /// positions will be read from the alignments' refpos, and the alignments
    /// must be surjected.
    HTSAlignmentEmitter(const string& filename, const string& format, const map<string, int64_t>& path_length, size_t max_threads);
    
    /// Tear down an HTSAlignmentEmitter and destroy HTSlib structures.
    ~HTSAlignmentEmitter();
    
    // Not copyable or movable
    HTSAlignmentEmitter(const HTSAlignmentEmitter& other) = delete;
    HTSAlignmentEmitter& operator=(const HTSAlignmentEmitter& other) = delete;
    HTSAlignmentEmitter(HTSAlignmentEmitter&& other) = delete;
    HTSAlignmentEmitter& operator=(HTSAlignmentEmitter&& other) = delete;
    
    
    /// Emit a batch of Alignments.
    virtual void emit_singles(vector<Alignment>&& aln_batch);
    /// Emit a batch of Alignments with secondaries. All secondaries must have
    /// is_secondary set already.
    virtual void emit_mapped_singles(vector<vector<Alignment>>&& alns_batch);
    /// Emit a batch of pairs of Alignments.
    virtual void emit_pairs(vector<Alignment>&& aln1_batch, vector<Alignment>&& aln2_batch,
        vector<int64_t>&& tlen_limit_batch);
    /// Emit the mappings of a batch of pairs of Alignments. All secondaries
    /// must have is_secondary set already.
    virtual void emit_mapped_pairs(vector<vector<Alignment>>&& alns1_batch,
        vector<vector<Alignment>>&& alns2_batch, vector<int64_t>&& tlen_limit_batch);
    
private:

    /// We hack about with htslib's BGZF EOF footers, so we need to know how long they are.
    static const size_t BGZF_FOOTER_LENGTH;

    /// If we are doing output to a file, this will hold the open file. Otherwise (for stdout) it will be empty.
    unique_ptr<ofstream> out_file;
    /// This holds a StreamMultiplexer on the output stream, for sharing it
    /// between threads.
    vg::io::StreamMultiplexer multiplexer;
    
    /// This holds our format name, for later error messages.
    string format;
    
    /// Store the path length map until the header can be made.
    map<string, int64_t> path_length;
    
    /// To back our samFile*s, we need the hFILE* objects wrapping our C++
    /// streams. We need to manually flush these after HTS headers are written,
    /// since bgzf_flush, which samtools calls, closes a BGZF block and sends
    /// the data to the hFILE* but does not actually flush the hFILE*.
    /// These will be pointers to the hFILE* for each thread's samFile*. We may
    /// only use them while the samFile* they belong to is still open; closing
    /// the samFile* will free the hFILE* but not null it out of this vector.
    vector<hFILE*> backing_files;
    
    /// We make one samFile* per thread, on each thread's output stream form
    /// the multiplexer. As soon as we create them, we show them the header, so
    /// they are initialized properly. If they have not yet been filled in
    /// (because the header is not ready yet), they are null.
    vector<samFile*> sam_files;
    
    /// We need a header
    atomic<bam_hdr_t*> atomic_header;
    /// We also need a header string.
    /// Not atomic, because by the time we read it we know the header is ready
    /// and nobody is writing to it.
    string sam_header;
    /// If the header isn't present when we want to write, we need a mutex to control creating it.
    mutex header_mutex;
    
    /// Remember if we are outputting BGZF-compressed data or not.
    /// If we are, we trim off spurious EOF markers and append our own.
    bool output_is_bgzf;
    
    /// Remember the HTSlib mode string we need to open our files.
    string hts_mode;
    
    virtual void convert_alignment(const Alignment& aln, vector<pair<int, char>>& cigar, bool& pos_rev, int64_t& pos, string& path_name) const;
    
    /// Convert an unpaired alignment to HTS format.
    /// Header must have been created already.
    void convert_unpaired(Alignment& aln, vector<bam1_t*>& dest);
    /// Convert a paired alignment to HTS format.
    /// Header must have been created already.
    void convert_paired(Alignment& aln1, Alignment& aln2, int64_t tlen_limit, vector<bam1_t*>& dest);
    
    /// Write and deallocate a bunch of BAM records. Takes care of locking the
    /// file. Header must have been written already.
    void save_records(bam_hdr_t* header, vector<bam1_t*>& records, size_t thread_number);
    
    /// Make sure that the HTS header has been written, and the samFile* in
    /// sam_files has been created for the given thread.
    ///
    /// If the header has not been written, blocks until it has been written.
    ///
    /// If we end up being the thread to write it, sniff header information
    /// from the given alignment.
    ///
    /// Returns the header pointer, so we don't have to do another atomic read
    /// later.
    bam_hdr_t* ensure_header(const Alignment& sniff, size_t thread_number);
    
    /// Given a header and a thread number, make sure the samFile* for that
    /// thread is initialized and ready to have alignments written to it. If
    /// true, actually writes the given header into the output file created by
    /// the multiplexer. If the samFile* was already initialized, flushes it
    /// out and makes a breakpoint.
    void initialize_sam_file(bam_hdr_t* header, size_t thread_number, bool keep_header = false);
    
};

class SplicedHTSAlignmentEmitter : public HTSAlignmentEmitter {
    
public:
    
    SplicedHTSAlignmentEmitter(const string& filename, const string& format,
                               const PathPositionHandleGraph& graph,
                               size_t max_threads);
    
    ~SplicedHTSAlignmentEmitter() = default;
    
    /// The minimum length of a deletion relative to the path that will be coded as a splice junction in the CIGAR
    size_t min_splice_length = 20;
    
private:

    /// Helper for constructor, makes path length map for parent class
    map<string, int64_t> make_path_length_index() const;
    
    /// Override for convert alignment that converts splices implicitly
    void convert_alignment(const Alignment& aln, vector<pair<int, char>>& cigar, bool& pos_rev, int64_t& pos, string& path_name) const;
    
    /// Convert a spliced alignment against a path to a cigar. The alignment must be
    /// colinear along a path and contain only mappings on the path, but it can have
    /// deletions relative to the path that follow edges in the graph.
    vector<pair<int, char>> spliced_cigar_against_path(const Alignment& aln, const string& path_name,
                                                       int64_t pos, bool rev) const;
    
    /// Graph that alignments were aligned against
    const PathPositionHandleGraph& graph;
    
};

/**
 * Emit Alignments to a stream in GAM or JSON format.
 * Thread safe.
 *
 * TODO: Split into Protobuf and JSON versions?
 */
class VGAlignmentEmitter : public AlignmentEmitter {
public:
    /// Create a VGAlignmentEmitter writing to the given file (or "-") in the given
    /// non-HTS format ("JSON", "GAM").
    VGAlignmentEmitter(const string& filename, const string& format, size_t max_threads);
    
    /// Finish and drstroy a VGAlignmentEmitter.
    ~VGAlignmentEmitter();
    
    /// Emit a batch of Alignments.
    virtual void emit_singles(vector<Alignment>&& aln_batch);
    /// Emit a batch of Alignments with secondaries. All secondaries must have
    /// is_secondary set already.
    virtual void emit_mapped_singles(vector<vector<Alignment>>&& alns_batch);
    /// Emit a batch of pairs of Alignments.
    virtual void emit_pairs(vector<Alignment>&& aln1_batch, vector<Alignment>&& aln2_batch,
        vector<int64_t>&& tlen_limit_batch);
    /// Emit the mappings of a batch of pairs of Alignments. All secondaries
    /// must have is_secondary set already.
    virtual void emit_mapped_pairs(vector<vector<Alignment>>&& alns1_batch,
        vector<vector<Alignment>>&& alns2_batch, vector<int64_t>&& tlen_limit_batch);
    
private:

    /// If we are doing output to a file, this will hold the open file. Otherwise (for stdout) it will be empty.
    unique_ptr<ofstream> out_file;
    
    /// This holds a StreamMultiplexer on the output stream, for sharing it
    /// between threads.
    vg::io::StreamMultiplexer multiplexer;
    
    /// We also keep ProtobufEmitters, one per thread, if we are doing protobuf output.
    vector<unique_ptr<vg::io::ProtobufEmitter<Alignment>>> proto;
};

}


#endif
