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
    /// Emit a pair of Alignments.
    virtual void emit_pair(Alignment&& aln1, Alignment&& aln2) = 0;
    /// Emit the mappings of a pair of Alignments. All secondaries must have is_secondary set already.
    virtual void emit_mapped_pair(vector<Alignment>&& alns1, vector<Alignment>&& alns2) = 0;
    
    /// Allow destruction through base class pointer.
    virtual ~AlignmentEmitter();
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
    HTSAlignmentEmitter(const string& filename, const string& format, map<string, int64_t>& path_length);
    
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
    virtual void emit_pair(Alignment&& aln1, Alignment&& aln2);
    /// Emit the mappings of a pair of Alignments. All secondaries must have is_secondary set already.
    virtual void emit_mapped_pair(vector<Alignment>&& alns1, vector<Alignment>&& alns2);
    
    /// If the distance between paired reads is above this limit, they will not
    /// be flagged as properly paired.
    size_t tlen_limit = 0;
    
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
    void emit_pair_internal(Alignment&& aln1, Alignment&& aln2, const lock_guard<mutex>& lock);
};

/**
 * Emit Alignments to a stream in GAM or JSON format.
 * Thread safe.
 *
 * TODO: Split into Protobuf and JSON versions?
 */
class VGAlignmentEmitter : public AlignmentEmitter {
public:
    /// Create an AlignmentEmitter writing to the given file (or "-") in the given
    /// non-HTS format ("JSON", "GAM").
    VGAlignmentEmitter(const string& filename, const string& format);
    
    /// Emit a single Alignment
    virtual void emit_single(Alignment&& aln);
    /// Emit a single Alignment with secondaries. All secondaries must have is_secondary set already.
    virtual void emit_mapped_single(vector<Alignment>&& alns);
    /// Emit a pair of Alignments.
    virtual void emit_pair(Alignment&& aln1, Alignment&& aln2);
    /// Emit the mappings of a pair of Alignments. All secondaries must have is_secondary set already.
    virtual void emit_mapped_pair(vector<Alignment>&& alns1, vector<Alignment>&& alns2);
    
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
