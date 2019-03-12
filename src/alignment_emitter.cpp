/**
 * \file alignment_emitter.cpp
 *
 * Implements a system for emitting alignments and groups of alignments in multiple formats.
 */

#include "alignment_emitter.hpp"
#include "alignment.hpp"

namespace vg {
using namespace std;

HTSAlignmentEmitter::HTSAlignmentEmitter(const string& filename, const string& format, map<string, int64_t>& path_length) : 
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
        cerr << "[vg::AlignmentEmitter] failed to open " << filename << " for writing HTS output" << endl;
        exit(1);
    }

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
            cerr << "[vg::AlignmentEmitter] error: failed to write the SAM header" << endl;
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
        cerr << "[vg::AlignmentEmitter] error: writing to output file failed" << endl;
        exit(1);
    }
    bam_destroy1(b);
}

void HTSAlignmentEmitter::emit_pair_internal(Alignment&& aln1, Alignment&& aln2, const lock_guard<mutex>& lock) {
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
            cerr << "[vg::AlignmentEmitter] error: failed to write the SAM header" << endl;
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
                                  this->tlen_limit);
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
                                  this->tlen_limit);
    
    int r = 0;
    r = sam_write1(sam_file, hdr, b1);
    if (r == 0) {
        cerr << "[vg::AlignmentEmitter] error: writing to output file failed" << endl;
        exit(1);
    }
    bam_destroy1(b1);
    r = 0;
    r = sam_write1(sam_file, hdr, b2);
    if (r == 0) {
        cerr << "[vg::AlignmentEmitter] error: writing to output file failed" << endl;
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

void HTSAlignmentEmitter::emit_pair(Alignment&& aln1, Alignment&& aln2) {
    lock_guard<mutex> lock(sync);
    // Emit the two paired alignments paired with each other.
    emit_pair_internal(std::move(aln1), std::move(aln2), lock);
}


void HTSAlignmentEmitter::emit_mapped_pair(vector<Alignment>&& alns1, vector<Alignment>&& alns2) {
    lock_guard<mutex> lock(sync);
    // Make sure we have the same number of mappings on each side.
    assert(alns1.size() == alns2.size());
    for (size_t i = 0; i < alns1.size(); i++) {
        // Emit each pair as paired alignments
        emit_pair_internal(std::move(alns1[i]), std::move(alns2[i]), lock);
    }
}

}
