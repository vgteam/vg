#ifndef VG_MULTIPATH_ALIGNMENT_EMITTER_HPP_INCLUDED
#define VG_MULTIPATH_ALIGNMENT_EMITTER_HPP_INCLUDED

/**
 * \file multipath_alignment_emitter.hpp
 *
 * Defines a system for emitting multipath alignments and groups of multipath alignments in multiple formats.
 */

#include <mutex>
#include <iostream>
#include <sstream>

#include <vg/vg.pb.h>
#include <vg/io/protobuf_emitter.hpp>
#include <vg/io/stream_multiplexer.hpp>
#include "multipath_alignment.hpp"

namespace vg {
using namespace std;

/*
 * Class that handles multithreaded output for multipath alignments
 */
class MultipathAlignmentEmitter {
public:
    
    /// Initialize with the intended output stream and the maximum number of threads that
    /// will be outputting. Optionally convert to single path alignments instead of multipath
    /// alignments.
    MultipathAlignmentEmitter(ostream& out, int num_threads, bool emit_single_path = false);
    ~MultipathAlignmentEmitter();
    
    /// Choose a read group to apply to all emitted alignments
    void set_read_group(const string& read_group);
    
    /// Choose a sample name to apply to all emitted alignments
    void set_sample_name(const string& sample_name);
    
    /// Emit paired read mappings as interleaved protobuf messages
    void emit_pairs(const string& name_1, const string& name_2,
                    vector<pair<multipath_alignment_t, multipath_alignment_t>>&& mp_aln_pairs);
    
    /// Emit read mappings as protobuf messages
    void emit_singles(const string& name, vector<multipath_alignment_t>&& mp_alns);
    
private:
    
    void convert_multipath_alignment(const multipath_alignment_t& mp_aln, Alignment& aln,
                                     const string* prev_name = nullptr,
                                     const string* next_name = nullptr) const;
    
    /// the stream that everything is emitted into
    ostream& out;
    
    /// shared multiplexer across streams
    vg::io::StreamMultiplexer multiplexer;
    
    /// an Alignment emitter for each thread
    vector<unique_ptr<vg::io::ProtobufEmitter<Alignment>>> aln_emitters;
    
    /// a MultipathAlignment emitter for each thread
    vector<unique_ptr<vg::io::ProtobufEmitter<MultipathAlignment>>> mp_aln_emitters;
    
    /// read group applied to alignments
    string read_group;
    
    /// sample name applied to alignments
    string sample_name;
    
};

}


#endif
