#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <iostream>
#include <functional>
#include "vg.hpp"

#include "htslib/hfile.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"

namespace vg {

const char* const BAM_DNA_LOOKUP = "=ACMGRSVTWYHKDBN";

int sam_for_each(string& filename, function<void(Alignment&)> lambda);
Alignment bam_to_alignment(bam1_t* b);
void write_alignments(std::ostream& out, vector<Alignment>& buf);
short quality_char_to_short_int(char c);
char quality_int_to_char(short i);

}

#endif
