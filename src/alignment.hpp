#ifndef VG_ALIGNMENT_HPP_INCLUDED
#define VG_ALIGNMENT_HPP_INCLUDED

#include <iostream>
#include <functional>
#include <zlib.h>
#include "utility.hpp"
#include "path.hpp"
#include "position.hpp"
#include <vg/vg.pb.h>
#include "vg/io/edit.hpp"
#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include "handle.hpp"
#include "vg/io/alignment_io.hpp"

namespace vg {

const char* const BAM_DNA_LOOKUP = "=ACMGRSVTWYHKDBN";

int hts_for_each(string& filename, function<void(Alignment&)> lambda);
int hts_for_each_parallel(string& filename, function<void(Alignment&)> lambda);
int hts_for_each(string& filename, function<void(Alignment&)> lambda,
                 const PathPositionHandleGraph* graph);
int hts_for_each_parallel(string& filename, function<void(Alignment&)> lambda,
                          const PathPositionHandleGraph* graph);
int fastq_for_each(string& filename, function<void(Alignment&)> lambda);

// fastq
bool get_next_alignment_from_fastq(gzFile fp, char* buffer, size_t len, Alignment& alignment);
bool get_next_interleaved_alignment_pair_from_fastq(gzFile fp, char* buffer, size_t len, Alignment& mate1, Alignment& mate2);
bool get_next_alignment_pair_from_fastqs(gzFile fp1, gzFile fp2, char* buffer, size_t len, Alignment& mate1, Alignment& mate2);

size_t fastq_unpaired_for_each(const string& filename, function<void(Alignment&)> lambda);
size_t fastq_paired_interleaved_for_each(const string& filename, function<void(Alignment&, Alignment&)> lambda);
size_t fastq_paired_two_files_for_each(const string& file1, const string& file2, function<void(Alignment&, Alignment&)> lambda);
// parallel versions of above
size_t fastq_unpaired_for_each_parallel(const string& filename,
                                        function<void(Alignment&)> lambda,
                                        uint64_t batch_size = vg::io::DEFAULT_PARALLEL_BATCHSIZE);
    
size_t fastq_paired_interleaved_for_each_parallel(const string& filename,
                                                  function<void(Alignment&, Alignment&)> lambda,
                                                  uint64_t batch_size = vg::io::DEFAULT_PARALLEL_BATCHSIZE);
    
size_t fastq_paired_interleaved_for_each_parallel_after_wait(const string& filename,
                                                             function<void(Alignment&, Alignment&)> lambda,
                                                             function<bool(void)> single_threaded_until_true,
                                                             uint64_t batch_size = vg::io::DEFAULT_PARALLEL_BATCHSIZE);
    
size_t fastq_paired_two_files_for_each_parallel(const string& file1, const string& file2,
                                                function<void(Alignment&, Alignment&)> lambda,
                                                uint64_t batch_size = vg::io::DEFAULT_PARALLEL_BATCHSIZE);
    
size_t fastq_paired_two_files_for_each_parallel_after_wait(const string& file1, const string& file2,
                                                           function<void(Alignment&, Alignment&)> lambda,
                                                           function<bool(void)> single_threaded_until_true,
                                                           uint64_t batch_size = vg::io::DEFAULT_PARALLEL_BATCHSIZE);

bam_hdr_t* hts_file_header(string& filename, string& header);
bam_hdr_t* hts_string_header(string& header,
                             const map<string, int64_t>& path_length,
                             const map<string, string>& rg_sample);
bam_hdr_t* hts_string_header(string& header,
                             const vector<pair<string, int64_t>>& path_order_and_length,
                             const map<string, string>& rg_sample);
void write_alignment_to_file(const Alignment& aln, const string& filename);

void mapping_cigar(const Mapping& mapping, vector<pair<int, char> >& cigar);
string cigar_string(const vector<pair<int, char> >& cigar);
string mapping_string(const string& source, const Mapping& mapping);

void cigar_mapping(const bam1_t *b, Mapping& mapping);

Alignment bam_to_alignment(const bam1_t *b,
                           const map<string, string>& rg_sample,
                           const map<int, path_handle_t>& tid_path_handle,
                           const bam_hdr_t *bh,
                           const PathPositionHandleGraph* graph);
Alignment bam_to_alignment(const bam1_t *b, const map<string, string>& rg_sample, const map<int, path_handle_t>& tid_path_handle);

/**
 * Add a CIGAR operation to a vector representing the parsed CIGAR string.
 *
 * Coalesces adjacent operations of the same type. Coalesces runs of inserts
 * and deletes into a signle delete followed by a single insert.
 */
inline void append_cigar_operation(const int length, const char operation, vector<pair<int, char>>& cigar) {
    if (cigar.empty()) {
        // Always append to an empty CIGAR
        cigar.emplace_back(length, operation);
    } else if (operation != cigar.back().second) {
        // We have changed operations
        if (operation == 'D' && cigar.back().second == 'I') {
            // This deletion needs to come before the adjacent insertion
            if (cigar.size() > 1 && cigar[cigar.size() - 2].second == 'D') {
                // Add to the deletion that laready exists before the insertion
                cigar[cigar.size() - 2].first += length;
            } else {
                // Create a new deletion
                cigar.emplace_back(length, operation);
                // Put it under the insertion
                std::swap(cigar[cigar.size() - 2], cigar.back());
            }
        } else {
            // This is an ordinary change of operations.
            cigar.emplace_back(length, operation);
        }
    } else {
        cigar.back().first += length;
    }
}

/**
 * Convert a paired Alignment to a BAM record. If the alignment is unmapped,
 * refpos must be -1. Otherwise, refpos must be the position on the reference
 * sequence to which the alignment is aligned. Similarly, refseq must be the
 * sequence aligned to, or "" if unaligned. The mateseq and matepos fields must
 * be set similarly for the mate. Note that mateseq must not be "=". If
 * tlen_max is given, it is a limit on the magnitude of tlen to consider the
 * read properly paired.
 *
 * Remember to clean up with bam_destroy1(b);
 */
bam1_t* alignment_to_bam(bam_hdr_t* bam_header,
                         const Alignment& alignment,
                         const string& refseq,
                         const int32_t refpos,
                         const bool refrev,
                         const vector<pair<int, char>>& cigar,
                         const string& mateseq,
                         const int32_t matepos,
                         bool materev,
                         const int32_t tlen,
                         const int32_t tlen_max = 0);
                         
/**
 * Convert an unpaired Alignment to a BAM record. If the alignment is unmapped,
 * refpos must be -1. Otherwise, refpos must be the position on the reference
 * sequence to which the alignment is aligned. Similarly, refseq must be the
 * sequence aligned to, or "" if unaligned.
 *
 * Remember to clean up with bam_destroy1(b);
 */
bam1_t* alignment_to_bam(bam_hdr_t* bam_header,
                         const Alignment& alignment,
                         const string& refseq,
                         const int32_t refpos,
                         const bool refrev,
                         const vector<pair<int, char>>& cigar);
                         
/**
 * Convert a paired Alignment to a SAM record. If the alignment is unmapped,
 * refpos must be -1. Otherwise, refpos must be the position on the reference
 * sequence to which the alignment is aligned. Similarly, refseq must be the
 * sequence aligned to, or "" if unaligned. The mateseq and matepos fields must
 * be set similarly for the mate. Note that mateseq must not be "=". If
 * tlen_max is given, it is a limit on the magnitude of tlen to consider the
 * read properly paired.
 */
string alignment_to_sam(const Alignment& alignment,
                        const string& refseq,
                        const int32_t refpos,
                        const bool refrev,
                        const vector<pair<int, char>>& cigar,
                        const string& mateseq,
                        const int32_t matepos,
                        bool materev,
                        const int32_t tlen,
                        const int32_t tlen_max = 0);
                        
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
                        const vector<pair<int, char>>& cigar);
                        

/// Returns the SAM bit-coded flag for alignment with
int32_t determine_flag(const Alignment& alignment,
                       const string& refseq,
                       const int32_t refpos,
                       const bool refrev,
                       const string& mateseq,
                       const int32_t matepos,
                       bool materev,
                       const int32_t tlen,
                       bool paired,
                       const int32_t tlen_max);

/// Create a CIGAR from the given Alignment. If softclip_suppress is nonzero,
/// suppress softclips up to that length. This will necessitate adjusting pos,
/// which is why it is passed by reference.
vector<pair<int, char>> cigar_against_path(const Alignment& alignment, bool on_reverse_strand, int64_t& pos, size_t path_len, size_t softclip_suppress);

/// Merge runs of successive I/D operations into a single I and D, remove 0-length
/// operations, and merge adjacent operations of the same type
void simplify_cigar(vector<pair<int, char>>& cigar);


/// Translate the CIGAR in the given BAM record into mappings in the given
/// Alignment against the given path in the given graph.
void mapping_against_path(Alignment& alignment, const bam1_t *b,
                          const path_handle_t& path, const PathPositionHandleGraph* graph,
                          bool on_reverse_strand);

/// Work out the TLEN values for two reads. The magnitude is the distance
/// between the outermost aligned bases, and the sign is positive for the
/// leftmost read and negative for the rightmost.
pair<int32_t, int32_t> compute_template_lengths(const int64_t& pos1, const vector<pair<int, char>>& cigar1,
    const int64_t& pos2, const vector<pair<int, char>>& cigar2);

int32_t sam_flag(const Alignment& alignment, bool on_reverse_strand, bool paired);
/// Populate a mapping from read group to sample name, given the text BAM header.
void parse_rg_sample_map(char* hts_header, map<string, string>& rg_sample);
/// Populate a mapping from target ID number to path handle in the given graph,
/// given a parsed BAM header. The graph may be null. Missing target paths in
/// the graph produce no warning or error and no map entry.
void parse_tid_path_handle_map(const bam_hdr_t* hts_header, const PathHandleGraph* graph, map<int, path_handle_t>& tid_path_handle);
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
int softclip_trim(Alignment& alignment);
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

/// Returns true if the alignment sequence contains any U's and false if the alignment sequence contains
/// and T's. In the case that both T's and U's are included, responds according to whichever comes first.
/// If the sequence contains neither U's nor T's, returns false.
bool uses_Us(const Alignment& alignment);

/// Replaces any U's in the sequence or the Path with T's
void convert_Us_to_Ts(Alignment& alignment);

/// Replaces any T's in the sequence or the Path with U's
void convert_Ts_to_Us(Alignment& alignment);

/// Simplifies the Path in the Alignment. Note that this removes deletions at
/// the start and end of Mappings, so code that handles simplified Alignments
/// needs to handle offsets on internal Mappings.
Alignment simplify(const Alignment& a, bool trim_internal_deletions = true);
    
/// Merge adjacent edits of the same type and convert all N matches to mismatches.
void normalize_alignment(Alignment& alignment);

// quality information; a kind of poor man's pileup
map<id_t, int> alignment_quality_per_node(const Alignment& aln);

/// Parse regions from the given BED file and call the given callback with each.
/// Does *not* write them to standard output.
/// Reads the optional name, is_reverse, and score fields if present, and populates the relevant Alignment fields.
/// Skips and warns about malformed or illegal BED records.
void parse_bed_regions(istream& bedstream, const PathPositionHandleGraph* graph, const std::function<void(Alignment&)>& callback);
/// Parse regions from the given GFF file and call the given callback with each.
/// Does *not* write them to standard output.
void parse_gff_regions(istream& gtfstream, const PathPositionHandleGraph* graph, const std::function<void(Alignment&)>& callback);
/// Parse regions from the given BED file into the given vector.
/// Does *not* write them to standard output.
void parse_bed_regions(istream& bedstream, const PathPositionHandleGraph* graph, vector<Alignment>* out_alignments);
/// Parse regions from the given GFF file into the given vector.
/// Does *not* write them to standard output.
void parse_gff_regions(istream& gtfstream, const PathPositionHandleGraph* graph, vector<Alignment>* out_alignments);

Position alignment_start(const Alignment& aln);
Position alignment_end(const Alignment& aln);Position alignment_start(const Alignment& aln);

/// return the path offsets as cached in the alignment
map<string ,vector<pair<size_t, bool> > > alignment_refpos_to_path_offsets(const Alignment& aln);
/// Annotate the first alignment with its minimum distance to the second in
/// their annotated paths. If translation is set, replace path names in aln
/// using that mapping, if they are found in it.
void alignment_set_distance_to_correct(Alignment& aln, const Alignment& base, const unordered_map<string, string>* translation = nullptr);
void alignment_set_distance_to_correct(Alignment& aln, const map<string, vector<pair<size_t, bool>>>& base_offsets, const unordered_map<string, string>* translation = nullptr);

/**
 * Represents a report on whether an alignment makes sense in the context of a graph.
 */
struct AlignmentValidity {
    /// The different kinds of possible problems with alignments
    enum Problem {
        OK,
        NODE_MISSING,
        NODE_TOO_SHORT,
        SEQ_DOES_NOT_MATCH
    };
    
    /// The kind of problem with the alignment.
    Problem problem = OK;
    /// The mapping in the alignment's path at which the problem was encountered.
    size_t bad_mapping_index = 0;
    /// An explanation for the problem.
    std::string message = "";
    
    /// We are truthy if the alignment has no problem, and falsey otherwise.
    inline operator bool() const {
        return problem == OK;
    }
};

/// Check to make sure edits on the alignment's path don't assume incorrect
/// node lengths or ids. Result can be used like a bool or inspected for
/// further details. Does not log anything itself about bad alignments.
AlignmentValidity alignment_is_valid(const Alignment& aln, const HandleGraph* hgraph, bool check_sequence = false);
    
/// Make an Alignment corresponding to a subregion of a stored path.
/// Positions are 0-based, and pos2 is excluded.
/// Respects path circularity, so pos2 < pos1 is not a problem.
/// If pos1 == pos2, returns an empty alignment.
Alignment target_alignment(const PathPositionHandleGraph* graph, const path_handle_t& path, size_t pos1, size_t pos2,
                           const string& feature, bool is_reverse);
/// Same as above, but uses the given Mapping, translated directly form a CIGAR string, as a source of edits.
/// The edits are inserted into the generated Alignment, cut as necessary to fit into the Alignment's Mappings.
Alignment target_alignment(const PathPositionHandleGraph* graph, const path_handle_t& path, size_t pos1, size_t pos2,
                           const string& feature, bool is_reverse, Mapping& cigar_mapping);

}

#endif
