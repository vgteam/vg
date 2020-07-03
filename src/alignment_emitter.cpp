/**
 * \file alignment_emitter.cpp
 *
 * Implements a system for emitting alignments and groups of alignments in multiple formats.
 */

#include "alignment_emitter.hpp"
#include "alignment_io.hpp"
#include "json2pb.h"
#include <vg/io/hfile_cppstream.hpp>
#include <vg/io/stream.hpp>
#include <omp.h>

#include <sstream>

//#define debug

namespace vg {
using namespace std;

// Implement all the single-read methods in terms of one-read batches
void AlignmentEmitter::emit_single(Alignment&& aln) {
    vector<Alignment> batch = { aln };
    emit_singles(std::move(batch));
}
void AlignmentEmitter::emit_mapped_single(vector<Alignment>&& alns) {
    vector<vector<Alignment>> batch = { alns };
    emit_mapped_singles(std::move(batch));
}
void AlignmentEmitter::emit_pair(Alignment&& aln1, Alignment&& aln2, int64_t tlen_limit) {
    vector<Alignment> batch1 = { aln1 };
    vector<Alignment> batch2 = { aln2 };
    vector<int64_t> tlen_limit_batch(1, tlen_limit);
    emit_pairs(std::move(batch1), std::move(batch2), std::move(tlen_limit_batch));
}
void AlignmentEmitter::emit_mapped_pair(vector<Alignment>&& alns1, vector<Alignment>&& alns2, int64_t tlen_limit) {
    vector<vector<Alignment>> batch1 = { alns1 };
    vector<vector<Alignment>> batch2 = { alns2 };
    vector<int64_t> tlen_limit_batch(1, tlen_limit);
    emit_mapped_pairs(std::move(batch1), std::move(batch2), std::move(tlen_limit_batch));
}

unique_ptr<AlignmentEmitter> get_non_hts_alignment_emitter(const string& filename, const string& format,
    const map<string, int64_t>& path_length, size_t max_threads, const HandleGraph* graph) {

    // Make the backing, non-buffered emitter
    AlignmentEmitter* backing = nullptr;
    if (format == "GAM" || format == "JSON") {
        // Make an emitter that supports VG formats
        backing = new VGAlignmentEmitter(filename, format, max_threads);
    } else if (format == "GAF") {
        backing = new GafAlignmentEmitter(filename, format, *graph, max_threads);
    } else if (format == "TSV") {
        backing = new TSVAlignmentEmitter(filename, max_threads);
    } else {
        cerr << "error [vg::get_non_hts_alignment_emitter]: Unimplemented output format " << format << endl;
        exit(1);
    }
    
    // Wrap it in a unique_ptr that will delete it
    return unique_ptr<AlignmentEmitter>(backing);
}

TSVAlignmentEmitter::TSVAlignmentEmitter(const string& filename, size_t max_threads) :
    out_file(filename == "-" ? nullptr : new ofstream(filename)),
    multiplexer(out_file.get() != nullptr ? *out_file : cout, max_threads) {
    
    if (out_file.get() != nullptr && !*out_file) {
        // Make sure we opened a file if we aren't writing to standard output
        cerr << "[vg::TSVAlignmentEmitter] failed to open " << filename << " for writing" << endl;
        exit(1);
    }
}

void TSVAlignmentEmitter::emit_singles(vector<Alignment>&& aln_batch) {
    for (auto&& aln : aln_batch) {
        emit(std::move(aln));
    }
    multiplexer.register_breakpoint(omp_get_thread_num());
}

void TSVAlignmentEmitter::emit_mapped_singles(vector<vector<Alignment>>&& alns_batch) {
    for (auto&& alns : alns_batch) {
        for (auto&& aln : alns) {
            emit(std::move(aln));
        }
    }
    multiplexer.register_breakpoint(omp_get_thread_num());
}

void TSVAlignmentEmitter::emit_pairs(vector<Alignment>&& aln1_batch,
                                     vector<Alignment>&& aln2_batch, 
                                     vector<int64_t>&& tlen_limit_batch) {
    // Ignore the tlen limit.
    assert(aln1_batch.size() == aln2_batch.size());
    for (size_t i = 0; i < aln1_batch.size(); i++) {
        // Emit each pair in order as read 1, then read 2
        emit(std::move(aln1_batch[i]));
        emit(std::move(aln2_batch[i]));
    }
    multiplexer.register_breakpoint(omp_get_thread_num());
}


void TSVAlignmentEmitter::emit_mapped_pairs(vector<vector<Alignment>>&& alns1_batch,
                                            vector<vector<Alignment>>&& alns2_batch,
                                            vector<int64_t>&& tlen_limit_batch) {
    assert(alns1_batch.size() == alns2_batch.size());
    for (size_t i = 0; i < alns1_batch.size(); i++) {
        // For each pair
        assert(alns1_batch[i].size() == alns2_batch[i].size());
        for (size_t j = 0; j < alns1_batch[i].size(); j++) {
            // Emit read 1 and read 2 pairs, together
            emit(std::move(alns1_batch[i][j]));
            emit(std::move(alns2_batch[i][j]));
        }
    }
    multiplexer.register_breakpoint(omp_get_thread_num());
}

void TSVAlignmentEmitter::emit(Alignment&& aln) {
    Position refpos;
    if (aln.refpos_size()) {
        refpos = aln.refpos(0);
    }

    // Get the stream to write to
    ostream& out = multiplexer.get_thread_stream(omp_get_thread_num());

    out << aln.name() << "\t"
        << refpos.name() << "\t"
        << refpos.offset() << "\t"
        << aln.mapping_quality() << "\t"
        << aln.score() << "\n";
}

VGAlignmentEmitter::VGAlignmentEmitter(const string& filename, const string& format, size_t max_threads):
    out_file(filename == "-" ? nullptr : new ofstream(filename)),
    multiplexer(out_file.get() != nullptr ? *out_file : cout, max_threads) {
    
    // We only support GAM and JSON formats
    assert(format == "GAM" || format == "JSON");
    
#ifdef debug
    cerr << "Creating VGAlignmentEmitter for " << format << " format to output file " << filename << " @ " << out_file.get() << endl;
    if (out_file.get() != nullptr) {
        cerr << "Output stream is at " << out_file->tellp() << endl;
    } else {
        cerr << "Output stream is at " << cout.tellp() << endl;
    }
#endif
    
    if (filename != "-") {
        // Check the file
        if (!*out_file) {
            // We couldn't get it open
            cerr << "[vg::VGAlignmentEmitter] failed to open " << filename << " for writing " << format << " output" << endl;
            exit(1);
        }
    }
    
    if (format == "GAM") {
        // We need per-thread emitters
        proto.reserve(max_threads);
        for (size_t i = 0; i < max_threads; i++) {
            // Make an emitter for each thread.
            proto.emplace_back(new vg::io::ProtobufEmitter<Alignment>(multiplexer.get_thread_stream(i)));
        }
    }
    
    // We later infer our format and output destination from out_file and proto being empty/set.
}

VGAlignmentEmitter::~VGAlignmentEmitter() {
#ifdef debug
    cerr << "Destroying VGAlignmentEmitter" << endl;
    if (out_file.get() != nullptr) {
        cerr << "Output stream is at " << out_file->tellp() << endl;
    } else {
        cerr << "Output stream is at " << cout.tellp() << endl;
    }
#endif

    if (!proto.empty()) {
        for (auto& emitter : proto) {
            // Flush each ProtobufEmitter
            emitter->flush(); 
            // Make it go away before the stream
            emitter.reset();
        }
    }
    
#ifdef debug
    cerr << "Destroyed VGAlignmentEmitter" << endl;
#endif
}

void VGAlignmentEmitter::emit_singles(vector<Alignment>&& aln_batch) {
    size_t thread_number = omp_get_thread_num();
    if (!proto.empty()) {
        // Save in protobuf
#ifdef debug
        #pragma omp critical (cerr)
        cerr << "VGAlignmentEmitter emitting " << aln_batch.size() << " reads to Protobuf in thread " << thread_number << endl;
#endif
        proto[thread_number]->write_many(std::move(aln_batch));
        if (multiplexer.want_breakpoint(thread_number)) {
            // The multiplexer wants our data.
            // Flush and create a breakpoint.
            proto[thread_number]->flush();
            multiplexer.register_breakpoint(thread_number);
        }
    } else {
        // Serialize to a string in our thread
        stringstream data;
        for (auto& aln : aln_batch) {
            multiplexer.get_thread_stream(thread_number) << pb2json(aln) << endl;
        }
        // No need to flush, we can always register a breakpoint.
        multiplexer.register_breakpoint(thread_number);
    }
}

void VGAlignmentEmitter::emit_mapped_singles(vector<vector<Alignment>>&& alns_batch) {
    size_t thread_number = omp_get_thread_num();
    if (!proto.empty()) {
        // Count up alignments
        size_t count = 0;
        for (auto& alns : alns_batch) {
            count += alns.size();
        }

#ifdef debug
        #pragma omp critical (cerr)
        cerr << "VGAlignmentEmitter emitting " << count << " alignments to Protobuf in thread " << thread_number << endl;
#endif
        
        if (count == 0) {
            // Nothing to do
            return;
        }
        
        // Collate one big vector to write together
        vector<Alignment> all;
        all.reserve(count);
        for (auto&& alns : alns_batch) {
            std::move(alns.begin(), alns.end(), std::back_inserter(all));
        }
        
        // Save in protobuf
        proto[thread_number]->write_many(std::move(all));
        if (multiplexer.want_breakpoint(thread_number)) {
            // The multiplexer wants our data.
            // Flush and create a breakpoint.
            proto[thread_number]->flush();
            multiplexer.register_breakpoint(thread_number);
#ifdef debug
            cerr << "Sent breakpoint from thread " << thread_number << endl;
#endif
        }
    } else {
        // Serialize to a string in our thread
        for (auto& alns : alns_batch) {
            for (auto& aln : alns) {
                multiplexer.get_thread_stream(thread_number) << pb2json(aln) << endl;
            }
        }
        // No need to flush, we can always register a breakpoint.
        multiplexer.register_breakpoint(thread_number);
#ifdef debug
        cerr << "Sent " << alns_batch.size() << " batches from thread " << thread_number << " followed by a breakpoint" << endl;
#endif
    }
}

void VGAlignmentEmitter::emit_pairs(vector<Alignment>&& aln1_batch,
                                    vector<Alignment>&& aln2_batch,
                                    vector<int64_t>&& tlen_limit_batch) {
    // Sizes need to match up
    assert(aln1_batch.size() == aln2_batch.size());
    assert(aln1_batch.size() == tlen_limit_batch.size());
    
    size_t thread_number = omp_get_thread_num();
    
    if (!proto.empty()) {
        // Save in protobuf
        
#ifdef debug
        #pragma omp critical (cerr)
        cerr << "VGAlignmentEmitter emitting paired reads to Protobuf in thread " << thread_number << endl;
#endif
        
        // Arrange into a vector in collated order
        vector<Alignment> all;
        all.reserve(aln1_batch.size() + aln2_batch.size());
        for (size_t i = 0; i < aln1_batch.size(); i++) {
            all.emplace_back(std::move(aln1_batch[i]));
            all.emplace_back(std::move(aln2_batch[i]));
        }
        
        // Save in protobuf
        proto[thread_number]->write_many(std::move(all));
        if (multiplexer.want_breakpoint(thread_number)) {
            // The multiplexer wants our data.
            // Flush and create a breakpoint.
            proto[thread_number]->flush();
            multiplexer.register_breakpoint(thread_number);
        }
    } else {
        // Serialize to a string in our thread in collated order
        stringstream data;
        for (size_t i = 0; i < aln1_batch.size(); i++) {
            multiplexer.get_thread_stream(thread_number) << pb2json(aln1_batch[i]) << endl
                << pb2json(aln2_batch[i]) << endl;
        }
        // No need to flush, we can always register a breakpoint.
        multiplexer.register_breakpoint(thread_number);
    }
}

void VGAlignmentEmitter::emit_mapped_pairs(vector<vector<Alignment>>&& alns1_batch,
                                           vector<vector<Alignment>>&& alns2_batch,
                                           vector<int64_t>&& tlen_limit_batch) {
    // Sizes need to match up
    assert(alns1_batch.size() == alns2_batch.size());
    assert(alns1_batch.size() == tlen_limit_batch.size());
    
    size_t thread_number = omp_get_thread_num();
    
    if (!proto.empty()) {
        // Save in protobuf
        
        // Count up all the alignments
        size_t count = 0;
        for (size_t i = 0; i < alns1_batch.size(); i++) {
            assert(alns1_batch[i].size() == alns2_batch[i].size());
            count += alns1_batch[i].size() * 2;
        }
        
#ifdef debug
        #pragma omp critical (cerr)
        cerr << "VGAlignmentEmitter emitting " << count << " mapped pairs to Protobuf in thread " << thread_number << endl;
#endif
        
        if (count == 0) {
            // Nothing to do
            return;
        }
        
        // Arrange into an interleaved vector
        vector<Alignment> all;
        all.reserve(count);
        for (size_t i = 0; i < alns1_batch.size(); i++) {
            for (size_t j = 0; j < alns1_batch[i].size(); j++) {
                all.emplace_back(std::move(alns1_batch[i][j]));
                all.emplace_back(std::move(alns2_batch[i][j]));
            }
        }
        
        // Save in protobuf
        proto[thread_number]->write_many(std::move(all));
        if (multiplexer.want_breakpoint(thread_number)) {
            // The multiplexer wants our data.
            // Flush and create a breakpoint.
            proto[thread_number]->flush();
            multiplexer.register_breakpoint(thread_number);
        }
    } else {
        // Serialize to an interleaved string in our thread
        stringstream data;
        for (size_t i = 0; i < alns1_batch.size(); i++) {
            assert(alns1_batch[i].size() == alns1_batch[i].size());
            for (size_t j = 0; j < alns1_batch[i].size(); j++) {
                multiplexer.get_thread_stream(thread_number) << pb2json(alns1_batch[i][j]) << endl
                    << pb2json(alns2_batch[i][j]) << endl;
            }
        }
        
        // No need to flush, we can always register a breakpoint.
        multiplexer.register_breakpoint(thread_number);
    }
}

GafAlignmentEmitter::GafAlignmentEmitter(const string& filename,
                                         const string& format,
                                         const HandleGraph& _graph,
                                         size_t max_threads):
    out_file(filename == "-" ? nullptr : new ofstream(filename)),
    multiplexer(out_file.get() != nullptr ? *out_file : cout, max_threads),
    graph(_graph) {
    
    // We only support GAF format
    assert(format == "GAF");
    
#ifdef debug
    cerr << "Creating GafAlignmentEmitter for " << format << " format to output file " << filename << " @ " << out_file.get() << endl;
    if (out_file.get() != nullptr) {
        cerr << "Output stream is at " << out_file->tellp() << endl;
    } else {
        cerr << "Output stream is at " << cout.tellp() << endl;
    }
#endif
    
    if (filename != "-") {
        // Check the file
        if (!*out_file) {
            // We couldn't get it open
            cerr << "[vg::GafAlignmentEmitter] failed to open " << filename << " for writing " << format << " output" << endl;
            exit(1);
        }
    }
    
    // We later infer our format and output destination from out_file and proto being empty/set.
}

GafAlignmentEmitter::~GafAlignmentEmitter() {
#ifdef debug
    cerr << "Destroying GafAlignmentEmitter" << endl;
    if (out_file.get() != nullptr) {
        cerr << "Output stream is at " << out_file->tellp() << endl;
    } else {
        cerr << "Output stream is at " << cout.tellp() << endl;
    }
#endif

#ifdef debug
    cerr << "Destroyed GafAlignmentEmitter" << endl;
#endif
}

void GafAlignmentEmitter::emit_singles(vector<Alignment>&& aln_batch) {
    size_t thread_number = omp_get_thread_num();
    // Serialize to a string in our thread
    for (auto& aln : aln_batch) {
        multiplexer.get_thread_stream(thread_number) << alignment_to_gaf(graph, aln) << endl;
    }
    // No need to flush, we can always register a breakpoint.
    multiplexer.register_breakpoint(thread_number);
}

void GafAlignmentEmitter::emit_mapped_singles(vector<vector<Alignment>>&& alns_batch) {
    size_t thread_number = omp_get_thread_num();
    // Serialize to a string in our thread
    for (auto& alns : alns_batch) {
        for (auto& aln : alns) {
            multiplexer.get_thread_stream(thread_number) << alignment_to_gaf(graph, aln) << endl;
        }
    }
    // No need to flush, we can always register a breakpoint.
    multiplexer.register_breakpoint(thread_number);
#ifdef debug
    cerr << "Sent " << alns_batch.size() << " batches from thread " << thread_number << " followed by a breakpoint" << endl;
#endif
}

void GafAlignmentEmitter::emit_pairs(vector<Alignment>&& aln1_batch,
                                    vector<Alignment>&& aln2_batch,
                                    vector<int64_t>&& tlen_limit_batch) {
    // Sizes need to match up
    assert(aln1_batch.size() == aln2_batch.size());
    assert(aln1_batch.size() == tlen_limit_batch.size());
    
    size_t thread_number = omp_get_thread_num();
    
    // Serialize to a string in our thread in collated order
    for (size_t i = 0; i < aln1_batch.size(); i++) {
        multiplexer.get_thread_stream(thread_number) << alignment_to_gaf(graph, aln1_batch[i]) << endl
                                                     << alignment_to_gaf(graph, aln2_batch[i]) << endl;
    }
    // No need to flush, we can always register a breakpoint.
    multiplexer.register_breakpoint(thread_number);
}

void GafAlignmentEmitter::emit_mapped_pairs(vector<vector<Alignment>>&& alns1_batch,
                                           vector<vector<Alignment>>&& alns2_batch,
                                           vector<int64_t>&& tlen_limit_batch) {
    // Sizes need to match up
    assert(alns1_batch.size() == alns2_batch.size());
    assert(alns1_batch.size() == tlen_limit_batch.size());
    
    size_t thread_number = omp_get_thread_num();
    // Serialize to an interleaved string in our thread
    for (size_t i = 0; i < alns1_batch.size(); i++) {
        assert(alns1_batch[i].size() == alns1_batch[i].size());
        for (size_t j = 0; j < alns1_batch[i].size(); j++) {
            multiplexer.get_thread_stream(thread_number) << alignment_to_gaf(graph, alns1_batch[i][j]) << endl
                                                         << alignment_to_gaf(graph, alns2_batch[i][j]) << endl;
        }
    }
    // No need to flush, we can always register a breakpoint.
    multiplexer.register_breakpoint(thread_number);
}

}
