#ifndef VG_ALIGNMENT_EMITTER_HPP_INCLUDED
#define VG_ALIGNMENT_EMITTER_HPP_INCLUDED

/**
 * \file alignment_emitter.hpp
 *
 * Defines a system for emitting alignments and groups of alignments in multiple formats.
 */

#include <mutex>
#include <vector>

#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

#include "vg.pb.h"
#include "stream/protobuf_emitter.hpp"

namespace vg {
using namespace std;

/**
 * Base class for a sink that takes alignments, possibly with pairing/secondary
 * relationships, and writes them out somewhere.
 *
 * All implementations must be thread safe.
 */
class AlignmentEmitter {
public:
    /// Emit a single Alignment
    virtual void emit_single(Alignment&& aln) = 0;
    /// Emit a single Alignment with secondaries. All secondaries must have is_secondary set already.
    virtual void emit_mapped_single(vector<Alignment>&& alns) = 0;
    /// Emit a pair of Alignments. The tlen_limit, if specified, is the maximum
    /// pairing distance to flag properly paired, if the output format cares
    /// about such things. TODO: Move to a properly paired annotation that runs
    /// with the Alignment.
    virtual void emit_pair(Alignment&& aln1, Alignment&& aln2, int64_t tlen_limit = 0) = 0;
    /// Emit the mappings of a pair of Alignments. All secondaries must have
    /// is_secondary set already. The tlen_limit, if specified, is the maximum
    /// pairing distance to flag properly paired, if the output format cares
    /// about such things. TODO: Move to a properly paired annotation that runs
    /// with the Alignment.
    virtual void emit_mapped_pair(vector<Alignment>&& alns1, vector<Alignment>&& alns2, int64_t tlen_limit = 0) = 0;
    
    /// Allow destruction through base class pointer.
    virtual ~AlignmentEmitter() = default;
};

/// Get an AlignmentEmitter that can emit to the given file (or "-") in the
/// given format. A table of contig lengths is required for HTSlib formats.
/// Automatically applies buffering.
unique_ptr<AlignmentEmitter> get_alignment_emitter(const string& filename, const string& format, const map<string, int64_t>& path_length);

/**
 * Throws per-OMP-thread buffers over the top of a backing AlignmentEmitter, which it owns.
 */
class OMPThreadBufferedAlignmentEmitter : public AlignmentEmitter {
public:
    /// Create an OMPThreadBufferedAlignmentEmitter that emits alignments to
    /// the given backing AlignmentEmitter. The backing emitter will become
    /// owned by this one.
    OMPThreadBufferedAlignmentEmitter(AlignmentEmitter* backing);
    
    /// Destroy and flush all the buffers
    ~OMPThreadBufferedAlignmentEmitter();
    
    /// Emit a single Alignment
    virtual void emit_single(Alignment&& aln);
    /// Emit a single Alignment with secondaries. All secondaries must have is_secondary set already.
    virtual void emit_mapped_single(vector<Alignment>&& alns);
    /// Emit a pair of Alignments.
    virtual void emit_pair(Alignment&& aln1, Alignment&& aln2, int64_t tlen_limit = 0);
    /// Emit the mappings of a pair of Alignments. All secondaries must have is_secondary set already.
    virtual void emit_mapped_pair(vector<Alignment>&& alns1, vector<Alignment>&& alns2, int64_t tlen_limit = 0);
    
private:
    /// Save all the buffered alignments from the given thread.
    void flush(size_t thread);

    /// Keep a reference to the backing emitter
    unique_ptr<AlignmentEmitter> backing;
    
    // We have one buffer for each type of emit operation, for each thread
    vector<vector<Alignment>> single_buffer;
    vector<vector<vector<Alignment>>> mapped_single_buffer;
    vector<vector<tuple<Alignment, Alignment, size_t>>> pair_buffer;
    vector<vector<tuple<vector<Alignment>, vector<Alignment>, size_t>>> mapped_pair_buffer;

    const static size_t BUFFER_LIMIT = 1000;
};

/**
 * Emit a TSV table describing alignments.
 */
class TSVAlignmentEmitter : public AlignmentEmitter {
public:
    
    /// Create a TSVAlignmentEmitter writing to the given file (or "-")
    TSVAlignmentEmitter(const string& filename);

    /// The default destructor should clean up the open file, if any.
    ~TSVAlignmentEmitter() = default;
    
    /// Emit a single Alignment
    virtual void emit_single(Alignment&& aln);
    /// Emit a single Alignment with secondaries. All secondaries must have is_secondary set already.
    virtual void emit_mapped_single(vector<Alignment>&& alns);
    /// Emit a pair of Alignments.
    virtual void emit_pair(Alignment&& aln1, Alignment&& aln2, int64_t tlen_limit = 0);
    /// Emit the mappings of a pair of Alignments. All secondaries must have is_secondary set already.
    virtual void emit_mapped_pair(vector<Alignment>&& alns1, vector<Alignment>&& alns2, int64_t tlen_limit = 0);
    
private:

    /// If we are doing output to a file, this will hold the open file. Otherwise (for stdout) it will be empty.
    unique_ptr<ofstream> out_file;

    /// Access to the output stream is protected by this mutex
    mutex sync;

    /// Emit a single alignment as TSV
    void emit_single_internal(Alignment&& aln, const lock_guard<mutex>& lock);
    
    /// Emit a pair of alignments as TSV, in order.
    void emit_pair_internal(Alignment&& aln1, Alignment&& aln2, const lock_guard<mutex>& lock);
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
    HTSAlignmentEmitter(const string& filename, const string& format, const map<string, int64_t>& path_length);
    
    /// Tear down an HTSAlignmentEmitter and destroy HTSlib structures.
    ~HTSAlignmentEmitter();
    
    // Not copyable or movable
    HTSAlignmentEmitter(const HTSAlignmentEmitter& other) = delete;
    HTSAlignmentEmitter& operator=(const HTSAlignmentEmitter& other) = delete;
    HTSAlignmentEmitter(HTSAlignmentEmitter&& other) = delete;
    HTSAlignmentEmitter& operator=(HTSAlignmentEmitter&& other) = delete;
    
    
    /// Emit a single Alignment
    virtual void emit_single(Alignment&& aln);
    /// Emit a single Alignment with secondaries. All secondaries must have is_secondary set already.
    virtual void emit_mapped_single(vector<Alignment>&& alns);
    /// Emit a pair of Alignments.
    virtual void emit_pair(Alignment&& aln1, Alignment&& aln2, int64_t tlen_limit = 0);
    /// Emit the mappings of a pair of Alignments. All secondaries must have is_secondary set already.
    virtual void emit_mapped_pair(vector<Alignment>&& alns1, vector<Alignment>&& alns2, int64_t tlen_limit = 0);
    
private:
    
    /// We need a mutex to synchronize on
    mutex sync;

    /// Remember what format we are using.
    string format;
    
    /// Sorte the path length map until the header can be made.
    map<string, int64_t> path_length;
    
    /// We need a samFile
    samFile* sam_file = nullptr;
    /// We need a header
    bam_hdr_t* hdr = nullptr;
    /// We also need a header string
    string sam_header;
    
    /// Emit a single alignment, with a lock already held.
    void emit_single_internal(Alignment&& aln, const lock_guard<mutex>& lock);
    /// Emit a pair of alignments, with a lock already held.
    void emit_pair_internal(Alignment&& aln1, Alignment&& aln2, int64_t tlen_limit, const lock_guard<mutex>& lock);
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
    VGAlignmentEmitter(const string& filename, const string& format);
    
    /// Finish and drstroy a VGAlignmentEmitter.
    ~VGAlignmentEmitter();
    
    /// Emit a single Alignment
    virtual void emit_single(Alignment&& aln);
    /// Emit a single Alignment with secondaries. All secondaries must have is_secondary set already.
    virtual void emit_mapped_single(vector<Alignment>&& alns);
    /// Emit a pair of Alignments.
    virtual void emit_pair(Alignment&& aln1, Alignment&& aln2, int64_t tlen_limit = 0);
    /// Emit the mappings of a pair of Alignments. All secondaries must have is_secondary set already.
    virtual void emit_mapped_pair(vector<Alignment>&& alns1, vector<Alignment>&& alns2, int64_t tlen_limit = 0);
    
private:

    /// If we are doing output to a file, this will hold the open file. Otherwise (for stdout) it will be empty.
    unique_ptr<ofstream> out_file;

    /// If we are doing Protobuf output we need a backing emitter. If we are doing JSON out, this will be empty.
    unique_ptr<stream::ProtobufEmitter<Alignment>> proto;
    
    /// If we are doing JSON, we need to take care of our own stream locking.
    mutex stream_mutex;
};

}


#endif
