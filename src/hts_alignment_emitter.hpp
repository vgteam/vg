#ifndef VG_HTS_ALIGNMENT_EMITTER_HPP_INCLUDED
#define VG_HTS_ALIGNMENT_EMITTER_HPP_INCLUDED

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
#include "vg/io/alignment_emitter.hpp"

namespace vg {
using namespace std;

using namespace vg::io;

/**
 * Flag enum for controlling the behavior of alignment emiotters behind get_alignment_emitter().
 */
enum alignment_emitter_flags_t {
    /// Value for no flags set.
    ALIGNMENT_EMITTER_FLAG_NONE = 0,
    /// Skip surjection, and expect pre-surjected alignments.
    ALIGNMENT_EMITTER_FLAG_HTS_RAW = 1,
    /// Use splicing-aware conversion to HTSlib formats:
    /// alignments are spliced at known splice sites (i.e. edges in the graph), so
    /// form spliced CIGAR strings
    ALIGNMENT_EMITTER_FLAG_HTS_SPLICED = 2,
    /// When surjecting, discard low-complexity anchors and realign more freely
    /// against the target path.
    ALIGNMENT_EMITTER_FLAG_HTS_PRUNE_SUSPICIOUS_ANCHORS = 4
};

/// Get an AlignmentEmitter that can emit to the given file (or "-") in the
/// given format. When writing HTSlib formats (SAM, BAM, CRAM), paths should
/// contain the paths in the linear reference in sequence dictionary order (see
/// get_sequence_dictionary), and a PathPositionHandleGraph must be provided.
/// When writing GAF, a HandleGraph must be provided for obtaining node lengths
/// and sequences. Other formats do not need a graph.
///
/// flags is an ORed together set of flags from alignment_emitter_flags_t.
///
/// Automatically applies per-thread buffering, but needs to know how many OMP
/// threads will be in use.
unique_ptr<AlignmentEmitter> get_alignment_emitter(const string& filename, const string& format, 
                                                   const vector<tuple<path_handle_t, size_t, size_t>>& paths, size_t max_threads,
                                                   const HandleGraph* graph = nullptr, int flags = ALIGNMENT_EMITTER_FLAG_NONE);

/**
 * Produce a list of path handles in a fixed order, suitable for use with
 * get_alignment_emitter_with_surjection(), by parsing a file. The file may be
 * an HTSlib-style "sequence dictionary" (consisting of SAM @SQ header lines),
 * or a plain list of sequence names (which do not start with "@SQ"). If the
 * file is not openable or contains no entries, reports an error and quits. If
 * the filename is itself an empty string, all non-alt-allele paths from the
 * graph will be collected in arbitrary order.
 *
 * TODO: Be able to generate the autosomes human-sort, X, Y, MT order typical
 * of references.
 *
 * The tuple is <path, path length in graph, base path length>
 * For a subpath (ie chr1[1000-10000]) the base path length would be that of chr1
 * This information needs to come from the user in order to be correct, but 
 * if it's not specified, it'll be guessed from the graph
 */
vector<tuple<path_handle_t, size_t, size_t>> get_sequence_dictionary(const string& filename, const PathPositionHandleGraph& graph);                                                   

/**
 * Given a list of path handles and size info (from get_sequence_dictionary), return two things:
 *  1) names and lengths of all of base paths in order.
 *  2) a mapping of path names to length (reflects paths in the graph including subpaths) 
 *
 * If subpath_support is set to false, there won't be a distinction. 
 */
pair<vector<pair<string, int64_t>>, unordered_map<string, int64_t>> extract_path_metadata(
    const vector<tuple<path_handle_t, size_t, size_t>>& paths,  const PathPositionHandleGraph& graph,
    bool subpath_support = false);

/*
 * A class that can write SAM/BAM/CRAM files from parallel threads
 */
class HTSWriter {
public:
    /// Create an HTSWriter writing to the given file (or "-") in the
    /// given HTS format ("SAM", "BAM", "CRAM"). path_order_and_length must give each
    /// contig name and length to include in the header. Sample names and read
    /// groups for the header will be guessed from the first reads. HTSlib
    /// positions will be read from the alignments' refpos, and the alignments
    /// must be surjected.
    HTSWriter(const string& filename, const string& format, const vector<pair<string, int64_t>>& path_order_and_length,
              const unordered_map<string, int64_t>& subpath_to_length, size_t max_threads);
    
    /// Tear down an HTSWriter and destroy HTSlib structures.
    ~HTSWriter();
    
    // Not copyable or movable
    HTSWriter(const HTSWriter& other) = delete;
    HTSWriter& operator=(const HTSWriter& other) = delete;
    HTSWriter(HTSWriter&& other) = delete;
    HTSWriter& operator=(HTSWriter&& other) = delete;
    
protected:
    
    /// We hack about with htslib's BGZF EOF footers, so we need to know how long they are.
    static const size_t BGZF_FOOTER_LENGTH;
    
    /// If we are doing output to a file, this will hold the open file. Otherwise (for stdout) it will be empty.
    unique_ptr<ofstream> out_file;
    /// This holds a StreamMultiplexer on the output stream, for sharing it
    /// between threads.
    vg::io::StreamMultiplexer multiplexer;
    
    /// This holds our format name, for later error messages.
    string format;
    
    /// Store the path names and lengths in the order to put them in the header.
    vector<pair<string, int64_t>> path_order_and_length;
    /// With subpath support, the above list will store base path inoformation for the header
    /// The actual path lengths go here:
    unordered_map<string, int64_t> subpath_to_length;
    
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
    bam_hdr_t* ensure_header(const string& read_group, const string& sample_name, size_t thread_number);
    
    /// Given a header and a thread number, make sure the samFile* for that
    /// thread is initialized and ready to have alignments written to it. If
    /// true, actually writes the given header into the output file created by
    /// the multiplexer. If the samFile* was already initialized, flushes it
    /// out and makes a breakpoint.
    void initialize_sam_file(bam_hdr_t* header, size_t thread_number, bool keep_header = false);
};

/**
 * Emit Alignments to a stream in SAM/BAM/CRAM format.
 * Thread safe.
 */
class HTSAlignmentEmitter : public AlignmentEmitter, public HTSWriter {
public:
    /// Create an HTSAlignmentEmitter writing to the given file (or "-") in the
    /// given HTS format ("SAM", "BAM", "CRAM"). path_order_and_length must give
    /// contig names and lengths to include in the header, in order. Sample
    /// names and read groups for the header will be guessed from the first
    /// reads. HTSlib positions will be read from the alignments' refpos, and
    /// the alignments must be surjected.
    HTSAlignmentEmitter(const string& filename, const string& format,
                        const vector<pair<string, int64_t>>& path_order_and_length,
                        const unordered_map<string, int64_t>& subpath_to_length, size_t max_threads);
    
    /// Tear down an HTSAlignmentEmitter and destroy HTSlib structures.
    ~HTSAlignmentEmitter() = default;
    
    // Not copyable or movable
    HTSAlignmentEmitter(const HTSAlignmentEmitter& other) = delete;
    HTSAlignmentEmitter& operator=(const HTSAlignmentEmitter& other) = delete;
    HTSAlignmentEmitter(HTSAlignmentEmitter&& other) = delete;
    HTSAlignmentEmitter& operator=(HTSAlignmentEmitter&& other) = delete;
    
    
    /// Emit a batch of Alignments.
    void emit_singles(vector<Alignment>&& aln_batch);
    /// Emit a batch of Alignments with secondaries. All secondaries must have
    /// is_secondary set already.
    void emit_mapped_singles(vector<vector<Alignment>>&& alns_batch);
    /// Emit a batch of pairs of Alignments.
    void emit_pairs(vector<Alignment>&& aln1_batch, vector<Alignment>&& aln2_batch,
        vector<int64_t>&& tlen_limit_batch);
    /// Emit the mappings of a batch of pairs of Alignments. All secondaries
    /// must have is_secondary set already.
    void emit_mapped_pairs(vector<vector<Alignment>>&& alns1_batch,
        vector<vector<Alignment>>&& alns2_batch, vector<int64_t>&& tlen_limit_batch);
    
private:
    
    virtual void convert_alignment(const Alignment& aln, vector<pair<int, char>>& cigar, bool& pos_rev, int64_t& pos, string& path_name) const;
    
    /// Convert an unpaired alignment to HTS format.
    /// Header must have been created already.
    void convert_unpaired(Alignment& aln, bam_hdr_t* header, vector<bam1_t*>& dest);
    /// Convert a paired alignment to HTS format.
    /// Header must have been created already.
    void convert_paired(Alignment& aln1, Alignment& aln2, bam_hdr_t* header, int64_t tlen_limit,
                        vector<bam1_t*>& dest);

};

/*
 * An HTSAlgnmentEmitter that tries to detect splice edges in
 * the input data so that they can be encoded as N CIGAR operations
 */
class SplicedHTSAlignmentEmitter : public HTSAlignmentEmitter {
    
public:
    
    SplicedHTSAlignmentEmitter(const string& filename, const string& format,
                               const vector<pair<string, int64_t>>& path_order_and_length,
                               const unordered_map<string, int64_t>& subpath_to_length,
                               const PathPositionHandleGraph& graph,
                               size_t max_threads);
    
    ~SplicedHTSAlignmentEmitter() = default;
    
    /// The minimum length of a deletion relative to the path that will be coded as a splice junction in the CIGAR
    size_t min_splice_length = 20;
    
private:
    
    /// Override for convert alignment that converts splices implicitly
    void convert_alignment(const Alignment& aln, vector<pair<int, char>>& cigar, bool& pos_rev, int64_t& pos, string& path_name) const;
    
    /// Convert a spliced alignment against a path to a cigar. The alignment must be
    /// colinear along a path and contain only mappings on the path, but it can have
    /// deletions relative to the path that follow edges in the graph.
    vector<pair<int, char>> spliced_cigar_against_path(const Alignment& aln, const string& path_name, int64_t pos,
                                                       bool rev) const;
    
    /// Graph that alignments were aligned against
    const PathPositionHandleGraph& graph;
    
};


}


#endif
