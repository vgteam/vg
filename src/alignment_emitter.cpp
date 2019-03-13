/**
 * \file alignment_emitter.cpp
 *
 * Implements a system for emitting alignments and groups of alignments in multiple formats.
 */

#include "alignment_emitter.hpp"
#include "alignment.hpp"
#include "json2pb.h"

#include <sstream>

namespace vg {
using namespace std;

unique_ptr<AlignmentEmitter> get_alignment_emitter(const string& filename, const string& format, const map<string, int64_t>& path_length) {

    // Make the backing, non-buffered emitter
    AlignmentEmitter* backing = nullptr;
    if (format == "GAM" || format == "JSON") {
        // Make an emitter that supports VG formats
        backing = new VGAlignmentEmitter(filename, format);
    } else if (format == "SAM" || format == "BAM" || format == "CRAM") {
        // Make an emitter that supports HTSlib formats
        backing = new HTSAlignmentEmitter(filename, format, path_length);
    } else if (format == "TSV") {
        backing = new TSVAlignmentEmitter(filename);
    } else {
        cerr << "error [vg::get_alignment_emitter]: Unimplemented output format " << format << endl;
        exit(1);
    }

    // Make a buffering emitter to own the backing one, and return it in a unique_ptr.
    return make_unique<OMPThreadBufferedAlignmentEmitter>(backing);
}

OMPThreadBufferedAlignmentEmitter::OMPThreadBufferedAlignmentEmitter(AlignmentEmitter* backing) : backing(backing) {
    #pragma omp parallel
    {
        #pragma omp single
        {
            size_t threads = omp_get_num_threads();
            single_buffer.resize(threads);
            pair_buffer.resize(threads);
            mapped_single_buffer.resize(threads);
            mapped_pair_buffer.resize(threads);
        }
    }
    
}

OMPThreadBufferedAlignmentEmitter::~OMPThreadBufferedAlignmentEmitter() {
    // Flush our buffers
    for (size_t i = 0; i < single_buffer.size(); i++) {
        flush(i);
    }
    // Then we can clean up our backing emitter
    backing.reset();
}

void OMPThreadBufferedAlignmentEmitter::flush(size_t thread) {
    // Flush all the buffers
    for (auto& aln : single_buffer[thread]) {
        backing->emit_single(std::move(aln));
    }
    single_buffer[thread].clear();
    for (auto& record : pair_buffer[thread]) {
        backing->emit_pair(std::move(get<0>(record)), std::move(get<1>(record)), get<2>(record));
    }
    pair_buffer[thread].clear();
    for (auto& alns : mapped_single_buffer[thread]) {
        backing->emit_mapped_single(std::move(alns));
    }
    mapped_single_buffer[thread].clear();
    for (auto& record : mapped_pair_buffer[thread]) {
        backing->emit_mapped_pair(std::move(get<0>(record)), std::move(get<1>(record)), get<2>(record));
    }
    mapped_pair_buffer[thread].clear();
}

void OMPThreadBufferedAlignmentEmitter::emit_single(Alignment&& aln) {
    size_t i = omp_get_thread_num();
    single_buffer[i].emplace_back(std::move(aln));
    if (single_buffer[i].size() >= BUFFER_LIMIT) {
        flush(i);
    }
}

void OMPThreadBufferedAlignmentEmitter::emit_mapped_single(vector<Alignment>&& alns) {
    size_t i = omp_get_thread_num();
    mapped_single_buffer[i].emplace_back(std::move(alns));
    if (mapped_single_buffer[i].size() >= BUFFER_LIMIT) {
        flush(i);
    }
}

void OMPThreadBufferedAlignmentEmitter::emit_pair(Alignment&& aln1, Alignment&& aln2, int64_t tlen_limit) {
    size_t i = omp_get_thread_num();
    pair_buffer[i].emplace_back(std::move(aln1), std::move(aln2), tlen_limit);
    if (pair_buffer[i].size() >= BUFFER_LIMIT) {
        flush(i);
    }
}

void OMPThreadBufferedAlignmentEmitter::emit_mapped_pair(vector<Alignment>&& alns1, vector<Alignment>&& alns2, int64_t tlen_limit) {
    size_t i = omp_get_thread_num();
    mapped_pair_buffer[i].emplace_back(std::move(alns1), std::move(alns2), tlen_limit);
    if (mapped_pair_buffer[i].size() >= BUFFER_LIMIT) {
        flush(i);
    }
}

TSVAlignmentEmitter::TSVAlignmentEmitter(const string& filename) : out_file(filename == "-" ? nullptr : new ofstream(filename)) {
    if (out_file.get() != nullptr && !*out_file) {
        // Make sure we opened a file if we aren't writing to standard output
        cerr << "[vg::TSVAlignmentEmitter] failed to open " << filename << " for writing" << endl;
        exit(1);
    }
}

void TSVAlignmentEmitter::emit_single(Alignment&& aln) {
    lock_guard<mutex> lock(sync);
    emit_single_internal(std::move(aln), lock);
}

void TSVAlignmentEmitter::emit_mapped_single(vector<Alignment>&& alns) {
    lock_guard<mutex> lock(sync);
    for (auto& aln : alns) {
        // Emit each mapping of the alignment in order in a single run
        emit_single_internal(std::move(aln), lock);
    }
}

void TSVAlignmentEmitter::emit_pair(Alignment&& aln1, Alignment&& aln2, int64_t tlen_limit) {
    lock_guard<mutex> lock(sync);
    // Emit the two paired alignments paired with each other.
    emit_pair_internal(std::move(aln1), std::move(aln2), lock);
}


void TSVAlignmentEmitter::emit_mapped_pair(vector<Alignment>&& alns1, vector<Alignment>&& alns2, int64_t tlen_limit) {
    lock_guard<mutex> lock(sync);
    // Make sure we have the same number of mappings on each side.
    assert(alns1.size() == alns2.size());
    for (size_t i = 0; i < alns1.size(); i++) {
        // Emit each pair as paired alignments
        emit_pair_internal(std::move(alns1[i]), std::move(alns2[i]), lock);
    }
}

void TSVAlignmentEmitter::emit_single_internal(Alignment&& aln, const lock_guard<mutex>& lock) {
    Position refpos;
    if (aln.refpos_size()) {
        refpos = aln.refpos(0);
    }

    // Work out if we are writing to standard output or an open file
    ostream& out = (out_file.get() == nullptr) ? cout : *out_file;

    out << aln.name() << "\t"
        << refpos.name() << "\t"
        << refpos.offset() << "\t"
        << aln.mapping_quality() << "\t"
        << aln.score() << "\n";
}

void TSVAlignmentEmitter::emit_pair_internal(Alignment&& aln1, Alignment&& aln2, const lock_guard<mutex>& lock) {
    emit_single_internal(std::move(aln1), lock);
    emit_single_internal(std::move(aln2), lock);
}

HTSAlignmentEmitter::HTSAlignmentEmitter(const string& filename, const string& format, const map<string, int64_t>& path_length) : 
    format(format), path_length(path_length) {
    
    // Make sure we have an HTS format
    assert(format == "SAM" || format == "BAM" || format == "CRAM");
    
    // Compute the file mode to send to HTSlib depending on output format
    char out_mode[5];
    string out_format = "";
    strcpy(out_mode, "w");
    if (format == "BAM") {
        out_format = "b";
    } else if (format == "CRAM") {
        out_format = "c";
    } else {
        // Must be SAM
        out_format = "";
    }
    strcat(out_mode, out_format.c_str());
    int compress_level = 9; // TODO: support other compression levels
    if (compress_level >= 0) {
        char tmp[2];
        tmp[0] = compress_level + '0'; tmp[1] = '\0';
        strcat(out_mode, tmp);
    }
    
    // Open the file
    sam_file = sam_open(filename.c_str(), out_mode);
    
    if (sam_file == nullptr) {
        // We couldn't open the output file.
        cerr << "[vg::HTSAlignmentEmitter] failed to open " << filename << " for writing " << format << " output" << endl;
        exit(1);
    }

}

HTSAlignmentEmitter::~HTSAlignmentEmitter() {
    if (hdr != nullptr) {
        bam_hdr_destroy(hdr);
    }
    sam_close(sam_file);
}

void HTSAlignmentEmitter::emit_single_internal(Alignment&& aln, const lock_guard<mutex>& lock) {
    // We are already locked, so we can touch the header.
    
    if (hdr == nullptr) {
        // Create and write the header
        
        // Sniff out the read group and sample, and map from RG to sample
        map<string, string> rg_sample;
        if (!aln.sample_name().empty() && !aln.read_group().empty()) {
            // We have a sample and a read group
            rg_sample[aln.read_group()] = aln.sample_name();
        }
        
        hdr = hts_string_header(sam_header, path_length, rg_sample);
        
        // write the header
        if (sam_hdr_write(sam_file, hdr) != 0) {
            cerr << "[vg::HTSAlignmentEmitter] error: failed to write the SAM header" << endl;
            exit(1);
        }
    }
    
    // Look up the stuff we need from the Alignment to express it in BAM.
    // We assume the position is available in refpos(0)
    assert(aln.refpos_size() == 1);
    
    size_t path_len = 0;
    if (aln.refpos(0).name() != "") {
        path_len = path_length.at(aln.refpos(0).name()); 
    }
    // Extract the position so that it could be adjusted by cigar_against_path if we decided to supperss softclips. Which we don't.
    // TODO: Separate out softclip suppression.
    int64_t pos = aln.refpos(0).offset();
    vector<pair<int, char>> cigar = cigar_against_path(aln, aln.refpos(0).is_reverse(), pos, path_len, 0);
    // TODO: We're passing along a text header so we can make a SAM file so
    // we can make a BAM record by re-reading it, which we can then
    // possibly output as SAM again. Make this less complicated.
    bam1_t* b = alignment_to_bam(sam_header,
                                 aln,
                                 aln.refpos(0).name(),
                                 pos,
                                 aln.refpos(0).is_reverse(),
                                 cigar);
    int r = 0;
    r = sam_write1(sam_file, hdr, b);
    if (r == 0) {
        cerr << "[vg::HTSAlignmentEmitter] error: writing to output file failed" << endl;
        exit(1);
    }
    bam_destroy1(b);
}

void HTSAlignmentEmitter::emit_pair_internal(Alignment&& aln1, Alignment&& aln2, int64_t tlen_limit, const lock_guard<mutex>& lock) {
    // We are already locked, so we can touch the header.
    
    if (hdr == nullptr) {
        // Create and write the header
        
        // Sniff out the read group and sample, and map from RG to sample
        map<string, string> rg_sample;
        if (!aln1.sample_name().empty() && !aln1.read_group().empty()) {
            // We have a sample and a read group
            rg_sample[aln1.read_group()] = aln1.sample_name();
        }
        
        hdr = hts_string_header(sam_header, path_length, rg_sample);
        
        // write the header
        if (sam_hdr_write(sam_file, hdr) != 0) {
            cerr << "[vg::HTSAlignmentEmitter] error: failed to write the SAM header" << endl;
            exit(1);
        }
    }
    
    // Look up the stuff we need from the Alignment to express it in BAM.
    // We assume the position is available in refpos(0)
    assert(aln1.refpos_size() == 1);
    assert(aln2.refpos_size() == 1);
    
    size_t path_len1 = 0;
    if (aln1.refpos(0).name() != "") {
        path_len1 = path_length.at(aln1.refpos(0).name()); 
    }
    size_t path_len2 = 0;
    if (aln2.refpos(0).name() != "") {
        path_len2 = path_length.at(aln2.refpos(0).name()); 
    }
    
    // Extract the position so that it could be adjusted by cigar_against_path if we decided to supperss softclips. Which we don't.
    // TODO: Separate out softclip suppression.
    int64_t pos1 = aln1.refpos(0).offset();
    int64_t pos2 = aln2.refpos(0).offset();
    vector<pair<int, char>> cigar1 = cigar_against_path(aln1, aln1.refpos(0).is_reverse(), pos1, path_len1, 0);
    vector<pair<int, char>> cigar2 = cigar_against_path(aln2, aln2.refpos(0).is_reverse(), pos2, path_len2, 0);
    
    // Determine the TLEN for each read.
    auto tlens = compute_template_lengths(pos1, cigar1, pos2, cigar2);
    
    // TODO: We're passing along a text header so we can make a SAM file so
    // we can make a BAM record by re-reading it, which we can then
    // possibly output as SAM again. Make this less complicated.
    bam1_t* b1 = alignment_to_bam(sam_header,
                                  aln1,
                                  aln1.refpos(0).name(),
                                  pos1,
                                  aln1.refpos(0).is_reverse(),
                                  cigar1,
                                  aln2.refpos(0).name(),
                                  pos2,
                                  aln2.refpos(0).is_reverse(),
                                  tlens.first,
                                  tlen_limit);
    bam1_t* b2 = alignment_to_bam(sam_header,
                                  aln2,
                                  aln2.refpos(0).name(),
                                  pos2,
                                  aln2.refpos(0).is_reverse(),
                                  cigar2,
                                  aln1.refpos(0).name(),
                                  pos1,
                                  aln1.refpos(0).is_reverse(),
                                  tlens.second,
                                  tlen_limit);
    
    int r = 0;
    r = sam_write1(sam_file, hdr, b1);
    if (r == 0) {
        cerr << "[vg::HTSAlignmentEmitter] error: writing to output file failed" << endl;
        exit(1);
    }
    bam_destroy1(b1);
    r = 0;
    r = sam_write1(sam_file, hdr, b2);
    if (r == 0) {
        cerr << "[vg::HTSAlignmentEmitter] error: writing to output file failed" << endl;
        exit(1);
    }
    bam_destroy1(b2);
}

void HTSAlignmentEmitter::emit_single(Alignment&& aln) {
    lock_guard<mutex> lock(sync);
    emit_single_internal(std::move(aln), lock);
}

void HTSAlignmentEmitter::emit_mapped_single(vector<Alignment>&& alns) {
    lock_guard<mutex> lock(sync);
    for (auto& aln : alns) {
        // Emit each mapping of the alignment in order in a single run
        emit_single_internal(std::move(aln), lock);
    }
}

void HTSAlignmentEmitter::emit_pair(Alignment&& aln1, Alignment&& aln2, int64_t tlen_limit) {
    lock_guard<mutex> lock(sync);
    // Emit the two paired alignments paired with each other.
    emit_pair_internal(std::move(aln1), std::move(aln2), tlen_limit, lock);
}


void HTSAlignmentEmitter::emit_mapped_pair(vector<Alignment>&& alns1, vector<Alignment>&& alns2, int64_t tlen_limit) {
    lock_guard<mutex> lock(sync);
    // Make sure we have the same number of mappings on each side.
    assert(alns1.size() == alns2.size());
    for (size_t i = 0; i < alns1.size(); i++) {
        // Emit each pair as paired alignments
        emit_pair_internal(std::move(alns1[i]), std::move(alns2[i]), tlen_limit, lock);
    }
}

VGAlignmentEmitter::VGAlignmentEmitter(const string& filename, const string& format) {
    // We only support GAM and JSON formats
    assert(format == "GAM" || format == "JSON");
    
    if (filename != "-") {
        // Open the file
        out_file = make_unique<ofstream>(filename);
        if (!*out_file) {
            // We couldn't get it open
            cerr << "[vg::VGAlignmentEmitter] failed to open " << filename << " for writing " << format << " output" << endl;
            exit(1);
        }
    }
    
    if (format == "GAM") {
        // We ened an emitter
        
        // Pick where to send the output
        ostream& out_stream = (filename != "-") ? *out_file : cout;
        
        // Point a ProtobufEmitter there
        proto = make_unique<stream::ProtobufEmitter<Alignment>>(out_stream);
    }
    
    // We later infer our format and output destination from out_file and proto being empty/set.
}

VGAlignmentEmitter::~VGAlignmentEmitter() {
    // TODO: We've had trouble with alignments going missing. Hopefully this fiexs that.

    if (proto.get() != nullptr) {
        // Flush the ProtobufEmitter
        proto->emit_group(); 
        // Make it go away before the stream
        proto.reset();
    }
    
    if (out_file.get() != nullptr) {
        // Flush our file stream
        out_file->flush();
        // Destroy it
        out_file.reset();
    } else {
        // Flush standard output
        cout.flush();
    }
}

void VGAlignmentEmitter::emit_single(Alignment&& aln) {
    if (proto.get() != nullptr) {
        // Save in protobuf
        proto->write(std::move(aln));
    } else {
        // Find a stream to write JSON to
        ostream& out = (out_file.get() != nullptr) ? *out_file : cout;
        
        // Serialize to a string in our thread
        string data = pb2json(aln);
        
        // Lock and emit the string
        lock_guard<mutex> lock(stream_mutex);
        out << data << endl;
    }
}

void VGAlignmentEmitter::emit_mapped_single(vector<Alignment>&& alns) {
    if (proto.get() != nullptr) {
        // Save in protobuf
        proto->write_many(std::move(alns));
    } else {
        // Find a stream to write JSON to
        ostream& out = (out_file.get() != nullptr) ? *out_file : cout;
        
        // Serialize to a string in our thread
        stringstream data;
        for (auto& aln : alns) {
            data << pb2json(aln) << endl;
        }
        
        // Lock and emit the string
        lock_guard<mutex> lock(stream_mutex);
        out << data.str();
    }
}

void VGAlignmentEmitter::emit_pair(Alignment&& aln1, Alignment&& aln2, int64_t tlen_limit) {
    if (proto.get() != nullptr) {
        // Save in protobuf
        
        // Arrange into a vector
        vector<Alignment> alns {std::move(aln1), std::move(aln2)};
        
        // Save as a single unit
        proto->write_many(std::move(alns));
    } else {
        // Find a stream to write JSON to
        ostream& out = (out_file.get() != nullptr) ? *out_file : cout;
        
        // Serialize to a string in our thread
        stringstream data;
        data << pb2json(aln1) << endl;
        data << pb2json(aln2) << endl;
        
        // Lock and emit the string
        lock_guard<mutex> lock(stream_mutex);
        out << data.str();
    }
}

void VGAlignmentEmitter::emit_mapped_pair(vector<Alignment>&& alns1, vector<Alignment>&& alns2, int64_t tlen_limit) {
    // Sizes need to match up
    assert(alns1.size() == alns2.size());
    
    if (proto.get() != nullptr) {
        // Save in protobuf
        
        // Arrange into an interleaved vector
        vector<Alignment> alns;
        alns.reserve(alns1.size() + alns2.size());
        for (size_t i = 0; i < alns1.size(); i++) {
            alns.emplace_back(std::move(alns1[i]));
            alns.emplace_back(std::move(alns2[i]));
        }
        
        // Save the interleaved vector
        proto->write_many(std::move(alns));
    } else {
        // Find a stream to write JSON to
        ostream& out = (out_file.get() != nullptr) ? *out_file : cout;
        
        // Serialize to an interleaved string in our thread
        stringstream data;
        for (size_t i = 0; i < alns1.size(); i++) {
            data << pb2json(alns1[i]) << endl;
            data << pb2json(alns2[i]) << endl;
        }
        
        // Lock and emit the string
        lock_guard<mutex> lock(stream_mutex);
        out << data.str();
    }
}

}
