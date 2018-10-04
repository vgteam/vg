#ifndef VG_ALIGNMENT_HPP_INCLUDED
#define VG_ALIGNMENT_HPP_INCLUDED

#include <iostream>
#include <functional>
#include <zlib.h>
#include "utility.hpp"
#include "path.hpp"
#include "position.hpp"
#include "vg.pb.h"
#include "xg.hpp"
#include "edit.hpp"
#include "htslib/hfile.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"

namespace vg {

const char* const BAM_DNA_LOOKUP = "=ACMGRSVTWYHKDBN";

int hts_for_each(string& filename, function<void(Alignment&)> lambda);
int hts_for_each_parallel(string& filename, function<void(Alignment&)> lambda);
int hts_for_each(string& filename, function<void(Alignment&)> lambda, xg::XG* xgindex);
int hts_for_each_parallel(string& filename, function<void(Alignment&)> lambda, xg::XG* xgindex);
int fastq_for_each(string& filename, function<void(Alignment&)> lambda);
bool get_next_alignment_from_fastq(gzFile fp, char* buffer, size_t len, Alignment& alignment);
bool get_next_interleaved_alignment_pair_from_fastq(gzFile fp, char* buffer, size_t len, Alignment& mate1, Alignment& mate2);
bool get_next_alignment_pair_from_fastqs(gzFile fp1, gzFile fp2, char* buffer, size_t len, Alignment& mate1, Alignment& mate2);

size_t fastq_unpaired_for_each(const string& filename, function<void(Alignment&)> lambda);
size_t fastq_paired_interleaved_for_each(const string& filename, function<void(Alignment&, Alignment&)> lambda);
size_t fastq_paired_two_files_for_each(const string& file1, const string& file2, function<void(Alignment&, Alignment&)> lambda);
// parallel versions of above
size_t fastq_unpaired_for_each_parallel(const string& filename,
                                        function<void(Alignment&)> lambda);
    
size_t fastq_paired_interleaved_for_each_parallel(const string& filename,
                                                  function<void(Alignment&, Alignment&)> lambda);
    
size_t fastq_paired_interleaved_for_each_parallel_after_wait(const string& filename,
                                                             function<void(Alignment&, Alignment&)> lambda,
                                                             function<bool(void)> single_threaded_until_true);
    
size_t fastq_paired_two_files_for_each_parallel(const string& file1, const string& file2,
                                                function<void(Alignment&, Alignment&)> lambda);
    
size_t fastq_paired_two_files_for_each_parallel_after_wait(const string& file1, const string& file2,
                                                           function<void(Alignment&, Alignment&)> lambda,
                                                           function<bool(void)> single_threaded_until_true);

bam_hdr_t* hts_file_header(string& filename, string& header);
bam_hdr_t* hts_string_header(string& header,
                             map<string, int64_t>& path_length,
                             map<string, string>& rg_sample);
void write_alignment_to_file(const Alignment& aln, const string& filename);

void mapping_cigar(const Mapping& mapping, vector<pair<int, char> >& cigar);
string cigar_string(vector<pair<int, char> >& cigar);
string mapping_string(const string& source, const Mapping& mapping);

void cigar_mapping(const bam1_t *b, Mapping& mapping, xg::XG* xgindex);

Alignment bam_to_alignment(const bam1_t *b, map<string, string>& rg_sample, const bam_hdr_t *bh, xg::XG* xgindex);
Alignment bam_to_alignment(const bam1_t *b, map<string, string>& rg_sample);

/**
 * Convert a paired Alignment to a BAM record. If the alignment is unmapped,
 * refpos must be -1. Otherwise, refpos must be the position on the reference
 * sequence to which the alignment is aligned. Similarly, refseq must be the
 * sequence aligned to, or "" if unaligned. The mateseq and matepos fields must
 * be set similarly for the mate. Note that mateseq must not be "=".
 *
 * Remember to clean up with bam_destroy1(b);
 */
bam1_t* alignment_to_bam(const string& sam_header,
                         const Alignment& alignment,
                         const string& refseq,
                         const int32_t refpos,
                         const bool refrev,
                         const string& cigar,
                         const string& mateseq,
                         const int32_t matepos,
                         const int32_t tlen);
                         
/**
 * Convert an unpaired Alignment to a BAM record. If the alignment is unmapped,
 * refpos must be -1. Otherwise, refpos must be the position on the reference
 * sequence to which the alignment is aligned. Similarly, refseq must be the
 * sequence aligned to, or "" if unaligned.
 *
 * Remember to clean up with bam_destroy1(b);
 */
bam1_t* alignment_to_bam(const string& sam_header,
                        const Alignment& alignment,
                        const string& refseq,
                        const int32_t refpos,
                        const bool refrev,
                        const string& cigar);
                         
/**
 * Convert a paired Alignment to a SAM record. If the alignment is unmapped,
 * refpos must be -1. Otherwise, refpos must be the position on the reference
 * sequence to which the alignment is aligned. Similarly, refseq must be the
 * sequence aligned to, or "" if unaligned. The mateseq and matepos fields must
 * be set similarly for the mate. Note that mateseq must not be "=".
 */
string alignment_to_sam(const Alignment& alignment,
                        const string& refseq,
                        const int32_t refpos,
                        const bool refrev,
                        const string& cigar,
                        const string& mateseq,
                        const int32_t matepos,
                        const int32_t tlen);
                        
/**
 * Convert an unpaired Alignment to a SAM record. If the alignment is unmapped,
 * refpos must be -1. Otherwise, refpos must be the position on the reference
 * sequence to which the alignment is aligned. Similarly, refseq must be the
 * sequence aligned to, or "" if unaligned.
 */
string alignment_to_sam(const Alignment& alignment,
                        const string& refseq,
                        const int32_t refpos,
                        const bool refrev,
                        const string& cigar);
                        


string cigar_against_path(const Alignment& alignment, bool on_reverse_strand, int64_t& pos, size_t path_len, size_t softclip_suppress);
void mapping_against_path(Alignment& alignment, const bam1_t *b, xg::XG* xgindex, bool on_reverse_strand);

int32_t sam_flag(const Alignment& alignment, bool on_reverse_strand, bool paired);
short quality_char_to_short(char c);
char quality_short_to_char(short i);
string string_quality_char_to_short(const string& quality);
string string_quality_short_to_char(const string& quality);
void alignment_quality_char_to_short(Alignment& alignment);
void alignment_quality_short_to_char(Alignment& alignment);
void parse_rg_sample_map(char* hts_header, map<string, string>& rg_sample);
int alignment_to_length(const Alignment& a);
int alignment_from_length(const Alignment& a);
// Adds a2 onto the end of a1, returns reference to a1
Alignment& extend_alignment(Alignment& a1, const Alignment& a2, bool debug=false);
// Merge a set of alignments into one
Alignment merge_alignments(const vector<Alignment>& alns);
// Merge two alignments end-to-end (could be "concat")
Alignment merge_alignments(const Alignment& a1, const Alignment& a2, bool debug=false);
Alignment strip_from_start(const Alignment& aln, size_t drop);
Alignment strip_from_end(const Alignment& aln, size_t drop);
Alignment trim_alignment(const Alignment& aln, const Position& pos1, const Position& pos2);
vector<Alignment> alignment_ends(const Alignment& aln, size_t len1, size_t len2);
Alignment alignment_middle(const Alignment& aln, int len);
// generate a digest of the alignmnet
const string hash_alignment(const Alignment& aln);
// Flip the alignment's sequence and is_reverse flag, and flip and re-order its
// Mappings to match. A function to get node lengths is needed because the
// Mappings in the alignment will need to give their positions from the opposite
// ends of their nodes. Offsets will be updated to count unused bases from node
// start when considering the node in its new orientation.
Alignment reverse_complement_alignment(const Alignment& aln, const function<int64_t(id_t)>& node_length);
void reverse_complement_alignment_in_place(Alignment* aln, const function<int64_t(id_t)>& node_length);
vector<Alignment> reverse_complement_alignments(const vector<Alignment>& alns, const function<int64_t(int64_t)>& node_length);
int non_match_start(const Alignment& alignment);
int non_match_end(const Alignment& alignment);
int softclip_start(const Alignment& alignment);
int softclip_end(const Alignment& alignment);
int query_overlap(const Alignment& aln1, const Alignment& aln2);
int edit_count(const Alignment& alignment);
size_t to_length_after_pos(const Alignment& aln, const Position& pos);
size_t from_length_after_pos(const Alignment& aln, const Position& pos);
size_t to_length_before_pos(const Alignment& aln, const Position& pos);
size_t from_length_before_pos(const Alignment& aln, const Position& pos);
string signature(const Alignment& aln);
pair<string, string> signature(const Alignment& aln1, const Alignment& aln2);
string middle_signature(const Alignment& aln, int len);
pair<string, string> middle_signature(const Alignment& aln1, const Alignment& aln2, int len);

// project the alignment's path back into a different ID space
void translate_nodes(Alignment& a, const unordered_map<id_t, pair<id_t, bool> >& ids, const std::function<size_t(int64_t)>& node_length);

// Invert the orientation in the alignment of all the nodes whose IDs are
// listed. It needs a callback to ask the length of any given node.
void flip_nodes(Alignment& a, const set<int64_t>& ids, const std::function<size_t(int64_t)>& node_length);

/// Simplifies the Path in the Alignment. Note that this removes deletions at
/// the start and end of Mappings, so code that handles simplified Alignments
/// needs to handle offsets on internal Mappings.
Alignment simplify(const Alignment& a, bool trim_internal_deletions = true);

// quality information; a kind of poor man's pileup
map<id_t, int> alignment_quality_per_node(const Alignment& aln);

/// Parse regions from the given BED file into Alignments in a vector.
/// Reads the optional name, is_reverse, and score fields if present, and populates the relevant Alignment fields.
/// Skips and warns about malformed or illegal BED records.
void parse_bed_regions(istream& bedstream, xg::XG* xgindex, vector<Alignment>* out_alignments);
void parse_gff_regions(istream& gtfstream, xg::XG* xgindex, vector<Alignment>* out_alignments);

Position alignment_start(const Alignment& aln);
Position alignment_end(const Alignment& aln);Position alignment_start(const Alignment& aln);

/// return the path offsets as cached in the alignment
map<string ,vector<pair<size_t, bool> > > alignment_refpos_to_path_offsets(const Alignment& aln);
/// annotate the first alignment with its minimum distance to the second in their annotated paths
void alignment_set_distance_to_correct(Alignment& aln, const Alignment& base);
void alignment_set_distance_to_correct(Alignment& aln, const map<string ,vector<pair<size_t, bool> > >& base_offsets);

}

#endif
