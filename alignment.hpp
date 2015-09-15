#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <iostream>
#include <functional>
#include <zlib.h>
#include "utility.hpp"
#include "path.hpp"
#include "vg.pb.h"

#include "htslib/hfile.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"

namespace vg {

const char* const BAM_DNA_LOOKUP = "=ACMGRSVTWYHKDBN";

int hts_for_each(string& filename, function<void(Alignment&)> lambda);
int hts_for_each_parallel(string& filename, function<void(Alignment&)> lambda);
int fastq_for_each(string& filename, function<void(Alignment&)> lambda);
bool get_next_alignment_from_fastq(gzFile fp, char* buffer, size_t len, Alignment& alignment);
bool get_next_interleaved_alignment_pair_from_fastq(gzFile fp, char* buffer, size_t len, Alignment& mate1, Alignment& mate2);
bool get_next_alignment_pair_from_fastqs(gzFile fp1, gzFile fp2, char* buffer, size_t len, Alignment& mate1, Alignment& mate2);

size_t fastq_unpaired_for_each(string& filename, function<void(Alignment&)> lambda);
size_t fastq_paired_interleaved_for_each(string& filename, function<void(Alignment&, Alignment&)> lambda);
size_t fastq_paired_two_files_for_each(string& file1, string& file2, function<void(Alignment&, Alignment&)> lambda);
// parallel versions of above
size_t fastq_unpaired_for_each_parallel(string& filename, function<void(Alignment&)> lambda);
size_t fastq_paired_interleaved_for_each_parallel(string& filename, function<void(Alignment&, Alignment&)> lambda);
size_t fastq_paired_two_files_for_each_parallel(string& file1, string& file2, function<void(Alignment&, Alignment&)> lambda);

bam_hdr_t* hts_file_header(string& filename, string& header);
bam_hdr_t* hts_string_header(string& header,
                             map<string, int64_t>& path_length,
                             map<string, string>& rg_sample);
void write_alignments(std::ostream& out, vector<Alignment>& buf);

Alignment bam_to_alignment(const bam1_t *b, map<string, string>& rg_sample);

bam1_t* alignment_to_bam(const string& sam_header,
                         const Alignment& alignment,
                         const string& refseq,
                         const int32_t refpos,
                         const string& cigar,
                         const string& mateseq,
                         const int32_t matepos,
                         const int32_t tlen);

string alignment_to_sam(const Alignment& alignment,
                        const string& refseq,
                        const int32_t refpos,
                        const string& cigar,
                        const string& mateseq,
                        const int32_t matepos,
                        const int32_t tlen);

string cigar_against_path(const Alignment& alignment);

int32_t sam_flag(const Alignment& alignment);
short quality_char_to_short(char c);
char quality_short_to_char(short i);
string string_quality_char_to_short(const string& quality);
string string_quality_short_to_char(const string& quality);
void alignment_quality_char_to_short(Alignment& alignment);
void alignment_quality_short_to_char(Alignment& alignment);
void parse_rg_sample_map(char* hts_header, map<string, string>& rg_sample);
int alignment_to_length(const Alignment& a);
int alignment_from_length(const Alignment& a);
Alignment merge_alignments(const vector<Alignment>& alns, const vector<size_t>& overlap);
// merge is destructive so we copy
Alignment merge_alignments(Alignment a1, Alignment a2, size_t overlap);
Alignment strip_from_start(const Alignment& aln, size_t drop);
Alignment strip_from_end(const Alignment& aln, size_t drop);
int softclip_start(Alignment& alignment);
int softclip_end(Alignment& alignment);
size_t to_length_after_pos(const Alignment& aln, const Position& pos);
size_t from_length_after_pos(const Alignment& aln, const Position& pos);
size_t to_length_before_pos(const Alignment& aln, const Position& pos);
size_t from_length_before_pos(const Alignment& aln, const Position& pos);

// Invert the orientation in the alignment of all the nodes whose IDs are
// listed. It needs a callback to ask the length of any given node.
void flip_nodes(Alignment& a, set<int64_t> ids, std::function<size_t(int64_t)> node_length);

}

#endif
