/**
 * \file multipath_alignment_emitter.cpp
 *
 * Implements a system for emitting multipath alignments and groups of multipath alignments in multiple formats.
 */

#include "multipath_alignment_emitter.hpp"
#include "json2pb.h"

namespace vg {
using namespace std;

MultipathAlignmentEmitter::MultipathAlignmentEmitter(ostream& out, int num_threads,
                                                     bool emit_single_path) :
    out(out),
    multiplexer(out, num_threads)
{
    // init the emitters for the correct output type
    if (emit_single_path) {
        aln_emitters.reserve(num_threads);
        for (int i = 0; i < num_threads; ++i) {
            aln_emitters.emplace_back(new vg::io::ProtobufEmitter<Alignment>(multiplexer.get_thread_stream(i)));
        }
    }
    else {
        mp_aln_emitters.reserve(num_threads);
        for (int i = 0; i < num_threads; ++i) {
            mp_aln_emitters.emplace_back(new vg::io::ProtobufEmitter<MultipathAlignment>(multiplexer.get_thread_stream(i)));
        }
    }
}

MultipathAlignmentEmitter::~MultipathAlignmentEmitter() {
    for (auto& emitter : aln_emitters) {
        // Flush each ProtobufEmitter
        emitter->flush();
        // Make it go away before the stream
        emitter.reset();
    }
    for (auto& emitter : mp_aln_emitters) {
        // Flush each ProtobufEmitter
        emitter->flush();
        // Make it go away before the stream
        emitter.reset();
    }
}

void MultipathAlignmentEmitter::emit_pairs(vector<pair<multipath_alignment_t, multipath_alignment_t>>&& mp_aln_pairs) {
    
    int thread_number = omp_get_thread_num();
    
    if (!mp_aln_emitters.empty()) {
        vector<MultipathAlignment> mp_alns_out(2 * mp_aln_pairs.size());
        for (size_t i = 0; i < mp_aln_pairs.size(); ++i) {
            to_proto_multipath_alignment(mp_aln_pairs[i].first, mp_alns_out[2 * i]);
            to_proto_multipath_alignment(mp_aln_pairs[i].second, mp_alns_out[2 * i + 1]);
        }
        
        mp_aln_emitters[thread_number]->write_many(std::move(mp_alns_out));
        
        if (multiplexer.want_breakpoint(thread_number)) {
            // The multiplexer wants our data.
            // Flush and create a breakpoint.
            mp_aln_emitters[thread_number]->flush();
            multiplexer.register_breakpoint(thread_number);
        }
    }
    else {
        vector<Alignment> alns_out(2 * mp_aln_pairs.size());
        for (size_t i = 0; i < mp_aln_pairs.size(); ++i) {
            convert_multipath_alignment(mp_aln_pairs[i].first, alns_out[2 * i],
                                        nullptr,
                                        &mp_aln_pairs[i].second);
            convert_multipath_alignment(mp_aln_pairs[i].second, alns_out[2 * i + 1],
                                        &mp_aln_pairs[i].first,
                                        nullptr);
        }
        
        aln_emitters[thread_number]->write_many(std::move(alns_out));
        
        if (multiplexer.want_breakpoint(thread_number)) {
            // The multiplexer wants our data.
            // Flush and create a breakpoint.
            aln_emitters[thread_number]->flush();
            multiplexer.register_breakpoint(thread_number);
        }
    }
}

void MultipathAlignmentEmitter::emit_singles(vector<multipath_alignment_t>&& mp_alns) {
    
    int thread_number = omp_get_thread_num();
    
    if (!mp_aln_emitters.empty()) {
        vector<MultipathAlignment> mp_alns_out(mp_alns.size());
        for (size_t i = 0; i < mp_alns.size(); ++i) {
            to_proto_multipath_alignment(mp_alns[i], mp_alns_out[i]);
        }
        
        mp_aln_emitters[thread_number]->write_many(std::move(mp_alns_out));
        
        if (multiplexer.want_breakpoint(thread_number)) {
            // The multiplexer wants our data.
            // Flush and create a breakpoint.
            mp_aln_emitters[thread_number]->flush();
            multiplexer.register_breakpoint(thread_number);
        }
    }
    else {
        vector<Alignment> alns_out(mp_alns.size());
        for (size_t i = 0; i < mp_alns.size(); ++i) {
            convert_multipath_alignment(mp_alns[i], alns_out[i]);
        }
        
        aln_emitters[thread_number]->write_many(std::move(alns_out));
        
        if (multiplexer.want_breakpoint(thread_number)) {
            // The multiplexer wants our data.
            // Flush and create a breakpoint.
            aln_emitters[thread_number]->flush();
            multiplexer.register_breakpoint(thread_number);
        }
    }
}

void MultipathAlignmentEmitter::convert_multipath_alignment(const multipath_alignment_t& mp_aln,
                                                            Alignment& aln,
                                                            const multipath_alignment_t* prev_pair,
                                                            const multipath_alignment_t* next_pair) const {
    optimal_alignment(mp_aln, aln);
    if (prev_pair) {
        aln.mutable_fragment_prev()->set_name(prev_pair->name());
    }
    if (next_pair) {
        aln.mutable_fragment_next()->set_name(next_pair->name());
    }
    aln.set_identity(identity(aln.path()));
}
}
