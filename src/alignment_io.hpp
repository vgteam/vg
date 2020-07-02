#ifndef VG_ALIGNMENT_IO_HPP_INCLUDED
#define VG_ALIGNMENT_IO_HPP_INCLUDED

#include <iostream>
#include <functional>
#include <zlib.h>
#include "utility.hpp"
#include "path.hpp"
#include "position.hpp"
#include <vg/vg.pb.h>
#include "edit.hpp"
#include "htslib/hfile.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "handle.hpp"
#include "gafkluge.hpp"

namespace vg {

const uint64_t DEFAULT_PARALLEL_BATCHSIZE = 512;

// general
size_t unpaired_for_each_parallel(function<bool(Alignment&)> get_read_if_available,
                                  function<void(Alignment&)> lambda,
                                  uint64_t batch_size = DEFAULT_PARALLEL_BATCHSIZE);

size_t paired_for_each_parallel_after_wait(function<bool(Alignment&, Alignment&)> get_pair_if_available,
                                           function<void(Alignment&, Alignment&)> lambda,
                                           function<bool(void)> single_threaded_until_true,
                                           uint64_t batch_size = DEFAULT_PARALLEL_BATCHSIZE);
// single gaf
bool get_next_alignment_from_gaf(const HandleGraph& graph, htsFile* fp, kstring_t& s_buffer, gafkluge::GafRecord& g_buffer,
                                 Alignment& alignment);
bool get_next_interleaved_alignment_pair_from_gaf(const HandleGraph& graph, htsFile* fp, kstring_t& s_buffer,
                                                  gafkluge::GafRecord& g_buffer, Alignment& mate1, Alignment& mate2);
size_t gaf_unpaired_for_each(const HandleGraph& graph, const string& filename, function<void(Alignment&)> lambda);
size_t gaf_paired_interleaved_for_each(const HandleGraph& graph, const string& filename,
                                       function<void(Alignment&, Alignment&)> lambda);
// parallel gaf
size_t gaf_unpaired_for_each_parallel(const HandleGraph& graph, const string& filename,
                                      function<void(Alignment&)> lambda,
                                      uint64_t batch_size = DEFAULT_PARALLEL_BATCHSIZE);
size_t gaf_paired_interleaved_for_each_parallel(const HandleGraph& graph, const string& filename,
                                                function<void(Alignment&, Alignment&)> lambda,
                                                uint64_t batch_size = DEFAULT_PARALLEL_BATCHSIZE);
size_t gaf_paired_interleaved_for_each_parallel_after_wait(const HandleGraph& graph, const string& filename,
                                                           function<void(Alignment&, Alignment&)> lambda,
                                                           function<bool(void)> single_threaded_until_true,
                                                           uint64_t batch_size = DEFAULT_PARALLEL_BATCHSIZE);
// gaf conversion
gafkluge::GafRecord alignment_to_gaf(const HandleGraph& graph, const Alignment& aln, bool cs_cigar = true, bool base_quals = true);
void gaf_to_alignment(const HandleGraph& graph, const gafkluge::GafRecord& gaf, Alignment& aln);

// utility
short quality_char_to_short(char c);
char quality_short_to_char(short i);
string string_quality_char_to_short(const string& quality);
string string_quality_short_to_char(const string& quality);
void alignment_quality_char_to_short(Alignment& alignment);
void alignment_quality_short_to_char(Alignment& alignment);

}
#endif
