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
#include "alignment.hpp"

namespace vg {
using namespace std;

/*
 * Class that handles multithreaded output for multipath alignments
 */
class MultipathAlignmentEmitter {
public:
    
    /// Initialize with the intended output stream and the maximum number of threads that
    /// will be outputting. Optionally converts to single path alignments in either "gam" or "gaf"
    /// format.
    MultipathAlignmentEmitter(ostream& out, int num_threads, const HandleGraph& graph,
                              const string out_format = "gamp");
    ~MultipathAlignmentEmitter();
    
    /// Emit paired read mappings as interleaved protobuf messages
    void emit_pairs(vector<pair<multipath_alignment_t, multipath_alignment_t>>&& mp_aln_pairs);
    
    /// Emit read mappings as protobuf messages
    void emit_singles(vector<multipath_alignment_t>&& mp_alns);
    
private:
    
    void convert_multipath_alignment(const multipath_alignment_t& mp_aln, Alignment& aln,
                                     const multipath_alignment_t* prev_frag = nullptr,
                                     const multipath_alignment_t* next_frag = nullptr) const;
    
    /// the graph we're aligning against
    const HandleGraph& graph;
    
    /// the stream that everything is emitted into
    ostream& out;
    
    /// shared multiplexer across streams
    vg::io::StreamMultiplexer multiplexer;
    
    /// an Alignment emitter for each thread
    vector<unique_ptr<vg::io::ProtobufEmitter<Alignment>>> aln_emitters;
    
    /// a MultipathAlignment emitter for each thread
    vector<unique_ptr<vg::io::ProtobufEmitter<MultipathAlignment>>> mp_aln_emitters;
    
};

}


#endif
