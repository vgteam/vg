#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <iostream>
#include <functional>
#include <zlib.h>
#include "vg.hpp"
#include "index.hpp"
#include "path.hpp"

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
int to_length(const Mapping& m);
int from_length(const Mapping& m);
void merge_alignments(Alignment& a1, const Alignment& a2);

}

#endif
