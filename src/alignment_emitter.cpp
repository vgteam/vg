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
    return make_unique<OMPSingleWriterBufferedAlignmentEmitter>(backing);
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
        
        // Now make sure we don't have to re-allocate while buffering.
        // With batches we can go over our flush threshold, so allocate twice as much space.
        size_t i = omp_get_thread_num();
        single_buffer[i].reserve(FLUSH_THRESHOLD * 2);
        get<0>(pair_buffer[i]).reserve(FLUSH_THRESHOLD * 2);
        get<1>(pair_buffer[i]).reserve(FLUSH_THRESHOLD * 2);
        get<2>(pair_buffer[i]).reserve(FLUSH_THRESHOLD * 2);
        mapped_single_buffer[i].reserve(FLUSH_THRESHOLD * 2);
        get<0>(mapped_pair_buffer[i]).reserve(FLUSH_THRESHOLD * 2);
        get<1>(mapped_pair_buffer[i]).reserve(FLUSH_THRESHOLD * 2);
        get<2>(mapped_pair_buffer[i]).reserve(FLUSH_THRESHOLD * 2);
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
    backing->emit_singles(std::move(single_buffer[thread]));
    single_buffer[thread].clear();
    single_buffer[thread].reserve(FLUSH_THRESHOLD * 2);
    
    backing->emit_pairs(std::move(get<0>(pair_buffer[thread])),
        std::move(get<1>(pair_buffer[thread])),
        std::move(get<2>(pair_buffer[thread]))); 
    get<0>(pair_buffer[thread]).clear();
    get<1>(pair_buffer[thread]).clear();
    get<2>(pair_buffer[thread]).clear();
    get<0>(pair_buffer[thread]).reserve(FLUSH_THRESHOLD * 2);
    get<1>(pair_buffer[thread]).reserve(FLUSH_THRESHOLD * 2);
    get<2>(pair_buffer[thread]).reserve(FLUSH_THRESHOLD * 2);
    
    backing->emit_mapped_singles(std::move(mapped_single_buffer[thread]));
    mapped_single_buffer[thread].clear();
    mapped_single_buffer[thread].reserve(FLUSH_THRESHOLD * 2);
     
    backing->emit_mapped_pairs(std::move(get<0>(mapped_pair_buffer[thread])),
        std::move(get<1>(mapped_pair_buffer[thread])),
        std::move(get<2>(mapped_pair_buffer[thread])));
    get<0>(mapped_pair_buffer[thread]).clear();
    get<1>(mapped_pair_buffer[thread]).clear();
    get<2>(mapped_pair_buffer[thread]).clear();
    get<0>(mapped_pair_buffer[thread]).reserve(FLUSH_THRESHOLD * 2);
    get<1>(mapped_pair_buffer[thread]).reserve(FLUSH_THRESHOLD * 2);
    get<2>(mapped_pair_buffer[thread]).reserve(FLUSH_THRESHOLD * 2);
}

void OMPThreadBufferedAlignmentEmitter::emit_singles(vector<Alignment>&& aln_batch) {
    size_t i = omp_get_thread_num();
    std::move(aln_batch.begin(), aln_batch.end(), std::back_inserter(single_buffer[i]));
    if (single_buffer[i].size() >= FLUSH_THRESHOLD) {
        flush(i);
    }
}

void OMPThreadBufferedAlignmentEmitter::emit_mapped_singles(vector<vector<Alignment>>&& alns_batch) {
    size_t i = omp_get_thread_num();
    std::move(alns_batch.begin(), alns_batch.end(), std::back_inserter(mapped_single_buffer[i]));
    if (mapped_single_buffer[i].size() >= FLUSH_THRESHOLD) {
        flush(i);
    }
}

void OMPThreadBufferedAlignmentEmitter::emit_pairs(vector<Alignment>&& aln1_batch,
                                                   vector<Alignment>&& aln2_batch,
                                                   vector<int64_t>&& tlen_limit_batch) {
    
    size_t i = omp_get_thread_num();
    std::move(aln1_batch.begin(), aln1_batch.end(), std::back_inserter(get<0>(pair_buffer[i])));
    std::move(aln2_batch.begin(), aln2_batch.end(), std::back_inserter(get<1>(pair_buffer[i])));
    std::move(tlen_limit_batch.begin(), tlen_limit_batch.end(), std::back_inserter(get<2>(pair_buffer[i])));
    if (get<0>(pair_buffer[i]).size() >= FLUSH_THRESHOLD) {
        flush(i);
    }
}

void OMPThreadBufferedAlignmentEmitter::emit_mapped_pairs(vector<vector<Alignment>>&& alns1_batch,
                                                          vector<vector<Alignment>>&& alns2_batch,
                                                          vector<int64_t>&& tlen_limit_batch) {
    size_t i = omp_get_thread_num();
    std::move(alns1_batch.begin(), alns1_batch.end(), std::back_inserter(get<0>(mapped_pair_buffer[i])));
    std::move(alns2_batch.begin(), alns2_batch.end(), std::back_inserter(get<1>(mapped_pair_buffer[i])));
    std::move(tlen_limit_batch.begin(), tlen_limit_batch.end(), std::back_inserter(get<2>(mapped_pair_buffer[i])));
    if (get<0>(mapped_pair_buffer[i]).size() >= FLUSH_THRESHOLD) {
        flush(i);
    }
}

OMPSingleWriterBufferedAlignmentEmitter::OMPSingleWriterBufferedAlignmentEmitter(AlignmentEmitter* backing) : backing(backing), stop(false) {
    #pragma omp parallel
    {
        #pragma omp single
        {
            size_t threads = omp_get_num_threads();
            single_buffer.resize(threads);
            pair_buffer.resize(threads);
            mapped_single_buffer.resize(threads);
            mapped_pair_buffer.resize(threads);
            buffer_mutexes.resize(threads);
        }
    }
    
    writer = thread(&OMPSingleWriterBufferedAlignmentEmitter::writer_thread_function, this);
    
}

OMPSingleWriterBufferedAlignmentEmitter::~OMPSingleWriterBufferedAlignmentEmitter() {
    // Stop our writer thread
    stop.store(true);
    // Wait on it
    writer.join();
    
    // Flush our buffers that are still not empty
    for (size_t i = 0; i < single_buffer.size(); i++) {
        flush(i);
    }
    // Then we can clean up our backing emitter
    backing.reset();
}


void OMPSingleWriterBufferedAlignmentEmitter::writer_thread_function() {
    
    
#ifdef debug
    // Track max buffer fill when we empty anything.
    // This has different semantics between data types but we ignore that.
    size_t buffer_high_water = 0;
#endif
    
    while(!stop.load()) {
        // Until we are told to stop
        
        for (size_t i = 0; i < buffer_mutexes.size(); i++) {
            // For each buffer, lock it
            if (buffer_mutexes[i].try_lock()) {

#ifdef debug
                buffer_high_water = max(buffer_high_water, single_buffer[i].size());
                buffer_high_water = max(buffer_high_water, get<0>(pair_buffer[i]).size());
                buffer_high_water = max(buffer_high_water, mapped_single_buffer[i].size());
                buffer_high_water = max(buffer_high_water, get<0>(mapped_pair_buffer[i]).size());
#endif
                
                
                // We got it. Flush.
                // TODO: Double-buffer so we don't stall this other thread.
                flush(i);
                
                buffer_mutexes[i].unlock();
            }
        }
        
    }
    
#ifdef debug
    cerr << "Buffer high water mark: " << buffer_high_water << " items" << endl;
#endif
}

void OMPSingleWriterBufferedAlignmentEmitter::flush(size_t thread) {
    // Flush all the buffers
    backing->emit_singles(std::move(single_buffer[thread]));
    single_buffer[thread].clear();
    
    backing->emit_pairs(std::move(get<0>(pair_buffer[thread])),
        std::move(get<1>(pair_buffer[thread])),
        std::move(get<2>(pair_buffer[thread]))); 
    get<0>(pair_buffer[thread]).clear();
    get<1>(pair_buffer[thread]).clear();
    get<2>(pair_buffer[thread]).clear();
    
    backing->emit_mapped_singles(std::move(mapped_single_buffer[thread]));
    mapped_single_buffer[thread].clear();
     
    backing->emit_mapped_pairs(std::move(get<0>(mapped_pair_buffer[thread])),
        std::move(get<1>(mapped_pair_buffer[thread])),
        std::move(get<2>(mapped_pair_buffer[thread])));
    get<0>(mapped_pair_buffer[thread]).clear();
    get<1>(mapped_pair_buffer[thread]).clear();
    get<2>(mapped_pair_buffer[thread]).clear();
}

void OMPSingleWriterBufferedAlignmentEmitter::emit_singles(vector<Alignment>&& aln_batch) {
    size_t i = omp_get_thread_num();
    buffer_mutexes[i].lock();
    while(single_buffer[i].size() >= MAX_BUFFER_SIZE) {
        // Avoid OOM by letting someone else run.
        // TODO: this is a wasteful spin where we could use a condition variable for buffer emptied.
        buffer_mutexes[i].unlock();
        std::this_thread::yield();
        buffer_mutexes[i].lock();
    }
    std::move(aln_batch.begin(), aln_batch.end(), std::back_inserter(single_buffer[i]));
    buffer_mutexes[i].unlock();
}

void OMPSingleWriterBufferedAlignmentEmitter::emit_mapped_singles(vector<vector<Alignment>>&& alns_batch) {
    size_t i = omp_get_thread_num();
    buffer_mutexes[i].lock();
    while(mapped_single_buffer[i].size() >= MAX_BUFFER_SIZE) {
        // Avoid OOM by letting someone else run.
        // TODO: this is a wasteful spin where we could use a condition variable for buffer emptied.
        buffer_mutexes[i].unlock();
        std::this_thread::yield();
        buffer_mutexes[i].lock();
    }
    std::move(alns_batch.begin(), alns_batch.end(), std::back_inserter(mapped_single_buffer[i]));
    buffer_mutexes[i].unlock();
}

void OMPSingleWriterBufferedAlignmentEmitter::emit_pairs(vector<Alignment>&& aln1_batch,
                                                   vector<Alignment>&& aln2_batch,
                                                   vector<int64_t>&& tlen_limit_batch) {
    
    size_t i = omp_get_thread_num();
    buffer_mutexes[i].lock();
    while(get<0>(pair_buffer[i]).size() >= MAX_BUFFER_SIZE) {
        // Avoid OOM by letting someone else run.
        // TODO: this is a wasteful spin where we could use a condition variable for buffer emptied.
        buffer_mutexes[i].unlock();
        std::this_thread::yield();
        buffer_mutexes[i].lock();
    }
    std::move(aln1_batch.begin(), aln1_batch.end(), std::back_inserter(get<0>(pair_buffer[i])));
    std::move(aln2_batch.begin(), aln2_batch.end(), std::back_inserter(get<1>(pair_buffer[i])));
    std::move(tlen_limit_batch.begin(), tlen_limit_batch.end(), std::back_inserter(get<2>(pair_buffer[i])));
    buffer_mutexes[i].unlock();
}

void OMPSingleWriterBufferedAlignmentEmitter::emit_mapped_pairs(vector<vector<Alignment>>&& alns1_batch,
                                                          vector<vector<Alignment>>&& alns2_batch,
                                                          vector<int64_t>&& tlen_limit_batch) {
    size_t i = omp_get_thread_num();
    buffer_mutexes[i].lock();
    while(get<0>(mapped_pair_buffer[i]).size() >= MAX_BUFFER_SIZE) {
        // Avoid OOM by letting someone else run.
        // TODO: this is a wasteful spin where we could use a condition variable for buffer emptied.
        buffer_mutexes[i].unlock();
        std::this_thread::yield();
        buffer_mutexes[i].lock();
    }
    std::move(alns1_batch.begin(), alns1_batch.end(), std::back_inserter(get<0>(mapped_pair_buffer[i])));
    std::move(alns2_batch.begin(), alns2_batch.end(), std::back_inserter(get<1>(mapped_pair_buffer[i])));
    std::move(tlen_limit_batch.begin(), tlen_limit_batch.end(), std::back_inserter(get<2>(mapped_pair_buffer[i])));
    buffer_mutexes[i].unlock();
}

TSVAlignmentEmitter::TSVAlignmentEmitter(const string& filename) : out_file(filename == "-" ? nullptr : new ofstream(filename)) {
    if (out_file.get() != nullptr && !*out_file) {
        // Make sure we opened a file if we aren't writing to standard output
        cerr << "[vg::TSVAlignmentEmitter] failed to open " << filename << " for writing" << endl;
        exit(1);
    }
}

void TSVAlignmentEmitter::emit_singles(vector<Alignment>&& aln_batch) {
    lock_guard<mutex> lock(sync);
    for (auto&& aln : aln_batch) {
        emit(std::move(aln), lock);
    }
}

void TSVAlignmentEmitter::emit_mapped_singles(vector<vector<Alignment>>&& alns_batch) {
    lock_guard<mutex> lock(sync);
    for (auto&& alns : alns_batch) {
        for (auto&& aln : alns) {
            emit(std::move(aln), lock);
        }
    }
}

void TSVAlignmentEmitter::emit_pairs(vector<Alignment>&& aln1_batch,
                                     vector<Alignment>&& aln2_batch, 
                                     vector<int64_t>&& tlen_limit_batch) {
    lock_guard<mutex> lock(sync);
    // Ignore the tlen limit.
    assert(aln1_batch.size() == aln2_batch.size());
    for (size_t i = 0; i < aln1_batch.size(); i++) {
        // Emit each pair in order as read 1, then read 2
        emit(std::move(aln1_batch[i]), lock);
        emit(std::move(aln2_batch[i]), lock);
    }
}


void TSVAlignmentEmitter::emit_mapped_pairs(vector<vector<Alignment>>&& alns1_batch,
                                            vector<vector<Alignment>>&& alns2_batch,
                                            vector<int64_t>&& tlen_limit_batch) {
    lock_guard<mutex> lock(sync);
    assert(alns1_batch.size() == alns2_batch.size());
    for (size_t i = 0; i < alns1_batch.size(); i++) {
        // For each pair
        assert(alns1_batch[i].size() == alns2_batch[i].size());
        for (size_t j = 0; j < alns1_batch[i].size(); j++) {
            // Emit read 1 and read 2 pairs, together
            emit(std::move(alns1_batch[i][j]), lock);
            emit(std::move(alns2_batch[i][j]), lock);
        }
    }
}

void TSVAlignmentEmitter::emit(Alignment&& aln, const lock_guard<mutex>& lock) {
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

HTSAlignmentEmitter::HTSAlignmentEmitter(const string& filename, const string& format, const map<string, int64_t>& path_length) : 
    format(format), path_length(path_length), sam_file(nullptr), file_mutex(), atomic_header(nullptr) {
    
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
    // Note that the destructor runs in only one thread, and only when
    // destruction is safe. No need to lock the header.
    if (atomic_header.load() != nullptr) {
        bam_hdr_destroy(atomic_header.load());
    }
    sam_close(sam_file);
}

bam_hdr_t* HTSAlignmentEmitter::ensure_header(const Alignment& sniff) {
    bam_hdr_t* header = atomic_header.load();
    if (header == nullptr) {
        // The header does not exist.
        
        // Lock the header mutex so we have exclusive control of the header
        lock_guard<mutex> header_lock(header_mutex);
        
        bam_hdr_t* header = atomic_header.load();
        if (header != nullptr) {
            // Someone else beat us to creating the header.
            // Header is ready.
            return header;
        }
        
        // Otherwise it is our turn to make the header.
        
        // Sniff out the read group and sample, and map from RG to sample
        map<string, string> rg_sample;
        if (!sniff.sample_name().empty() && !sniff.read_group().empty()) {
            // We have a sample and a read group
            rg_sample[sniff.read_group()] = sniff.sample_name();
        }
        
        // Make the header
        header = hts_string_header(sam_header, path_length, rg_sample);
        
        {
            // Lock the file
            lock_guard<mutex> file_lock(file_mutex);
        
            // Write the header to the file.
            if (sam_hdr_write(sam_file, header) != 0) {
                cerr << "[vg::HTSAlignmentEmitter] error: failed to write the SAM header" << endl;
                exit(1);
            }
        }
        
        // Save back to the atomic only after the header has been written and
        // it is safe for other threads to use it.
        atomic_header.store(header);
    }
    
    return header;
}

void HTSAlignmentEmitter::convert_unpaired(Alignment& aln, vector<bam1_t*>& dest) {
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
    dest.emplace_back(alignment_to_bam(sam_header,
                                       aln,
                                       aln.refpos(0).name(),
                                       pos,
                                       aln.refpos(0).is_reverse(),
                                       cigar));
}

void HTSAlignmentEmitter::convert_paired(Alignment& aln1, Alignment& aln2, int64_t tlen_limit, vector<bam1_t*>& dest) {
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
    dest.emplace_back(alignment_to_bam(sam_header,
                                       aln1,
                                       aln1.refpos(0).name(),
                                       pos1,
                                       aln1.refpos(0).is_reverse(),
                                       cigar1,
                                       aln2.refpos(0).name(),
                                       pos2,
                                       aln2.refpos(0).is_reverse(),
                                       tlens.first,
                                       tlen_limit));
    dest.emplace_back(alignment_to_bam(sam_header,
                                       aln2,
                                       aln2.refpos(0).name(),
                                       pos2,
                                       aln2.refpos(0).is_reverse(),
                                       cigar2,
                                       aln1.refpos(0).name(),
                                       pos1,
                                       aln1.refpos(0).is_reverse(),
                                       tlens.second,
                                       tlen_limit));
    
}

void HTSAlignmentEmitter::save_records(bam_hdr_t* header, vector<bam1_t*>& records) {
    {
        // Lock the file for output
        lock_guard<mutex> file_lock(file_mutex);
        for (auto& b : records) {
            // Emit each record
            if (sam_write1(sam_file, header, b) == 0) {
                cerr << "[vg::HTSAlignmentEmitter] error: writing to output file failed" << endl;
                exit(1);
            }
        }
    }
    
    for (auto& b : records) {
        // After unlocking the file, deallocate all the records
        bam_destroy1(b);
    }
}

void HTSAlignmentEmitter::emit_singles(vector<Alignment>&& aln_batch) {
    if (aln_batch.empty()) {
        // Nothing to do
        return;
    }
    
    // Make sure header exists
    bam_hdr_t* header = ensure_header(aln_batch.front());
    
    vector<bam1_t*> records;
    records.reserve(aln_batch.size());
    
    for (auto& aln : aln_batch) {
        // Convert each alignment to HTS format
        convert_unpaired(aln, records);
    }
    
    // Save to the file, locking once.
    save_records(header, records);
}

void HTSAlignmentEmitter::emit_mapped_singles(vector<vector<Alignment>>&& alns_batch) {
    // Count the total alignments to do
    size_t count = 0;
    // And find an alignment to base the header on
    Alignment* sniff = nullptr;
    for (auto& alns : alns_batch) {
        count += alns.size();
        if (!alns.empty() && sniff == nullptr) {
            sniff = &alns.front();
        }
    }
    
    if (count == 0) {
        // Nothing to do
        return;
    }
    
    // Make sure header exists
    assert(sniff != nullptr);
    bam_hdr_t* header = ensure_header(*sniff);
    
    vector<bam1_t*> records;
    records.reserve(count);
    
    for (auto& alns : alns_batch) {
        for (auto& aln : alns) {
            // Convert each alignment to HTS format
            convert_unpaired(aln, records);
        }
    }
    
    // Save to the file, locking once.
    save_records(header, records);

}


void HTSAlignmentEmitter::emit_pairs(vector<Alignment>&& aln1_batch,
                                     vector<Alignment>&& aln2_batch,
                                     vector<int64_t>&& tlen_limit_batch) {
    assert(aln1_batch.size() == aln2_batch.size());
    assert(aln1_batch.size() == tlen_limit_batch.size());
    
    
    if (aln1_batch.empty()) {
        // Nothing to do
        return;
    }
    
    // Make sure header exists
    bam_hdr_t* header = ensure_header(aln1_batch.front());
    
    vector<bam1_t*> records;
    records.reserve(aln1_batch.size() * 2);
    
    for (size_t i = 0; i < aln1_batch.size(); i++) {
        // Convert each alignment pair to HTS format
        convert_paired(aln1_batch[i], aln2_batch[i], tlen_limit_batch[i], records);
    }
    
    // Save to the file, locking once.
    save_records(header, records);
}
    
    
void HTSAlignmentEmitter::emit_mapped_pairs(vector<vector<Alignment>>&& alns1_batch, 
                                            vector<vector<Alignment>>&& alns2_batch,
                                            vector<int64_t>&& tlen_limit_batch) {
                                            
    assert(alns1_batch.size() == alns2_batch.size());
    assert(alns1_batch.size() == tlen_limit_batch.size());
    
    // Count the total alignments to do
    size_t count = 0;
    // And find an alignment to base the header on
    Alignment* sniff = nullptr;
    for (size_t i = 0; i < alns1_batch.size(); i++) {
        // Go through all pairs
        
        // Make sure each end of the pair has the same number of mappings
        assert(alns1_batch[i].size() == alns2_batch[i].size());
        
        count += alns1_batch[i].size() * 2;
        if (!alns1_batch[i].empty() && sniff == nullptr) {
            sniff = &alns1_batch[i].front();
        }
    }
    
    if (count == 0) {
        // Nothing to do
        return;
    }
    
    // Make sure header exists
    assert(sniff != nullptr);
    bam_hdr_t* header = ensure_header(*sniff);
    
    vector<bam1_t*> records;
    records.reserve(count);
    
    for (size_t i = 0; i < alns1_batch.size(); i++) {
        for (size_t j = 0; j < alns1_batch[i].size(); j++) {
            // Convert each alignment pair to HTS format
            convert_paired(alns1_batch[i][j], alns2_batch[i][j], tlen_limit_batch[i], records);
        }
    }
    
    // Save to the file, locking once.
    save_records(header, records);
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
        proto = make_unique<vg::io::ProtobufEmitter<Alignment>>(out_stream);
    }
    
    // We later infer our format and output destination from out_file and proto being empty/set.
}

VGAlignmentEmitter::~VGAlignmentEmitter() {
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

void VGAlignmentEmitter::emit_singles(vector<Alignment>&& aln_batch) {
    if (proto.get() != nullptr) {
        // Save in protobuf with a single lock
        proto->write_many(std::move(aln_batch));
    } else {
        // Find a stream to write JSON to
        ostream& out = (out_file.get() != nullptr) ? *out_file : cout;
        
        // Serialize to a string in our thread
        stringstream data;
        for (auto& aln : aln_batch) {
            data << pb2json(aln) << endl;
        }
        
        // Lock and emit the string
        lock_guard<mutex> lock(stream_mutex);
        out << data.str();
    }
}

void VGAlignmentEmitter::emit_mapped_singles(vector<vector<Alignment>>&& alns_batch) {
    if (proto.get() != nullptr) {
        // Count up alignments
        size_t count = 0;
        for (auto& alns : alns_batch) {
            count += alns.size();
        }
        
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
        
        // Save in protobuf with a single lock
        proto->write_many(std::move(all));
    } else {
        // Find a stream to write JSON to
        ostream& out = (out_file.get() != nullptr) ? *out_file : cout;
        
        // Serialize to a string in our thread
        stringstream data;
        for (auto& alns : alns_batch) {
            for (auto& aln : alns) {
                data << pb2json(aln) << endl;
            }
        }
        
        // Lock and emit the string
        lock_guard<mutex> lock(stream_mutex);
        out << data.str();
    }
}

void VGAlignmentEmitter::emit_pairs(vector<Alignment>&& aln1_batch,
                                    vector<Alignment>&& aln2_batch,
                                    vector<int64_t>&& tlen_limit_batch) {
    // Sizes need to match up
    assert(aln1_batch.size() == aln2_batch.size());
    assert(aln1_batch.size() == tlen_limit_batch.size());
    if (proto.get() != nullptr) {
        // Save in protobuf
        
        // Arrange into a vector in collated order
        vector<Alignment> all;
        all.reserve(aln1_batch.size() + aln2_batch.size());
        for (size_t i = 0; i < aln1_batch.size(); i++) {
            all.emplace_back(std::move(aln1_batch[i]));
            all.emplace_back(std::move(aln2_batch[i]));
        }
        
        // Save as a single unit with a single lock
        proto->write_many(std::move(all));
    } else {
        // Find a stream to write JSON to
        ostream& out = (out_file.get() != nullptr) ? *out_file : cout;
        
        // Serialize to a string in our thread in collated order
        stringstream data;
        for (size_t i = 0; i < aln1_batch.size(); i++) {
            data << pb2json(aln1_batch[i]) << endl;
            data << pb2json(aln2_batch[i]) << endl;
        }
        
        // Lock and emit the string
        lock_guard<mutex> lock(stream_mutex);
        out << data.str();
    }
}

void VGAlignmentEmitter::emit_mapped_pairs(vector<vector<Alignment>>&& alns1_batch,
                                           vector<vector<Alignment>>&& alns2_batch,
                                           vector<int64_t>&& tlen_limit_batch) {
    // Sizes need to match up
    assert(alns1_batch.size() == alns2_batch.size());
    assert(alns1_batch.size() == tlen_limit_batch.size());
    
    if (proto.get() != nullptr) {
        // Save in protobuf
        
        // Count up all the alignments
        size_t count = 0;
        for (size_t i = 0; i < alns1_batch.size(); i++) {
            assert(alns1_batch[i].size() == alns2_batch[i].size());
            count += alns1_batch[i].size() * 2;
        }
        
        if (count == 0) {
            // Nothing to do
            return;
        }
        
        // Arrange into an interleaved vector
        vector<Alignment> all;
        all.reserve(count);
        for (size_t i = 0; i < alns1_batch.size(); i++) {
            for (size_t j = 0; i < alns1_batch[i].size(); j++) {
                all.emplace_back(std::move(alns1_batch[i][j]));
                all.emplace_back(std::move(alns2_batch[i][j]));
            }
        }
        
        // Save the interleaved vector with one lock
        proto->write_many(std::move(all));
    } else {
        // Find a stream to write JSON to
        ostream& out = (out_file.get() != nullptr) ? *out_file : cout;
        
        // Serialize to an interleaved string in our thread
        stringstream data;
        for (size_t i = 0; i < alns1_batch.size(); i++) {
            assert(alns1_batch[i].size() == alns1_batch[i].size());
            for (size_t j = 0; i < alns1_batch[i].size(); j++) {
                data << pb2json(alns1_batch[i][j]) << endl;
                data << pb2json(alns2_batch[i][j]) << endl;
            }
        }
        
        // Lock and emit the string
        lock_guard<mutex> lock(stream_mutex);
        out << data.str();
    }
}

}
