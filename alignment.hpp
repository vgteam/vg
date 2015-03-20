#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <iostream>
#include <functional>
#include "vg.hpp"
#include "index.hpp"

#include "htslib/hfile.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"

namespace vg {

const char* const BAM_DNA_LOOKUP = "=ACMGRSVTWYHKDBN";

int hts_for_each(string& filename, function<void(Alignment&)> lambda);
int hts_for_each_parallel(string& filename, function<void(Alignment&)> lambda);
Alignment bam_to_alignment(const bam1_t *b, map<string, string>& rg_sample);
//bam1_t* alignment_to_bam(Alignment& alignment, Index& index, map<string, string>& rg_sample);
bool project_alignment(Alignment& alignment, Index& index, string& path_name);
void write_alignments(std::ostream& out, vector<Alignment>& buf);
short quality_char_to_short(char c);
char quality_short_to_char(short i);
void alignment_quality_char_to_short(Alignment& alignment);
void alignment_quality_short_to_char(Alignment& alignment);
void parse_rg_sample_map(char* hts_header, map<string, string>& rg_sample);
int to_length(Mapping& m);
int from_length(Mapping& m);

}

#endif
