#ifndef VG_READFILTER_HPP_INCLUDED
#define VG_READFILTER_HPP_INCLUDED

#include <vector>
#include <cstdlib>
#include <iostream>
#include <string>
#include <regex>
#include "vg.hpp"
#include "handle.hpp"
#include <vg/vg.pb.h>

/** \file
 * Provides a way to filter and transform reads, implementing the bulk of the
 * `vg filter` command.
 *
 */
namespace vg{

using namespace std;

class ReadFilter{
public:
    
    // Filtering parameters
    
    /// Actually take the complement of the filter
    bool complement_filter = false;
    /// Read name must have one of these prefixes, if any are present.
    /// TODO: This should be a trie but I don't have one handy.
    /// Must be sorted for vaguely efficient search.
    vector<string> name_prefixes;
    /// Read must not have a refpos set with a contig name containing a match to any of these
    vector<regex> excluded_refpos_contigs;
    /// If a read has one of the features in this set as annotations, the read
    /// is filtered out.
    unordered_set<string> excluded_features;
    double min_secondary = 0.;
    double min_primary = 0.;
    /// Should we rescore each alignment with default parameters and no e.g.
    /// haplotype info?
    bool rescore = false;
    bool frac_score = false;
    bool sub_score = false;
    int max_overhang = 99999;
    int min_end_matches = 0;
    bool verbose = false;
    double min_mapq = 0.;
    int repeat_size = 0;
    /// Should we drop split reads that follow edges not in the graph?
    bool drop_split = false;
    
    /// We can also pseudorandomly drop reads. What's the probability that we keep a read?
    double downsample_probability = 1.0;
    /// Samtools-compatible internal seed mask, for deciding which read pairs to keep.
    /// To be generated with rand() after srand() from the user-visible seed.
    uint32_t downsample_seed_mask = 0;
    
    /// How far in from the end should we look for ambiguous end alignment to
    /// clip off?
    int defray_length = 0;
    /// Limit defray recursion to visit this many nodes
    int defray_count = 99999;
    
    /// Number of threads from omp
    int threads = -1;
    /// GAM output buffer size
    int buffer_size = 512;
    /// Sometimes we only want a report, and not a filtered gam.  toggling off output
    /// speeds things up considerably.
    bool write_output = true;
    /// A HandleGraph is required for some filters (Note: ReadFilter doesn't own/free this)
    const HandleGraph* graph = nullptr;
    /// Interleaved input
    bool interleaved = false;
    /// When outputting paired reads, fail the pair only if both (all) reads
    /// fail (true) instead of if either (any) read fails (false)
    bool filter_on_all = false;
    
    // minimum base quality as PHRED score
    int min_base_quality = 0;
    // minimum fraction of bases in reads that must have quality at least <min_base_quality>
    double min_base_quality_fraction = 0.0;
                     
    // Keep some basic counts for when verbose mode is enabled
    struct Counts {
        enum FilterName { read = 0, wrong_name, wrong_refpos, excluded_feature, min_score, min_sec_score, max_overhang,
                          min_end_matches, min_mapq, split, repeat, defray, defray_all, random, min_base_qual, filtered,
                          last };
        vector<size_t> counts;
        Counts () : counts(FilterName::last, 0) {}
        Counts& operator+=(const Counts& other) {
            for (int i = 0; i < FilterName::last; ++i) {
                counts[i] += other.counts[i];
            }
            return *this;
        }
        /// If any read was filtered, count the other read as filtered
        Counts& set_paired_any() {
            for (int i = 0; i < FilterName::last; ++i) {
                counts[i] = counts[i] == 1 ? 2 : counts[i];
            }
            return *this;
        }
        /// If not all reads were filtered, count filtered ones as unfiltered.
        Counts& set_paired_all() {
            // We know that the pair as a whole was filtered out if counts[FilterName::filtered] == 2, and that it was kept otherwise.
            if (counts[FilterName::filtered] != 2) {
                // The read pair was not discarded, so clear out all the fail counts (including FilterName::filtered).
                for (int i = 0; i < FilterName::last; ++i) {
                    counts[i] = 0;
                }
            }
            // Otherwise the read pair was discarded, so leave all the
            // counts alone. There is at least one fail on each side to be
            // responsible for it (and if not verbose, only one fail on
            // each side).
            return *this;
        }
        void reset() {
            std::fill(counts.begin(), counts.end(), 0);
        }
        bool keep() {
            return counts[FilterName::filtered] == 0;
        }
    };
    
    bool append_regions = false;

    /**
     * Run all the filters on an alignment. The alignment may get modified in-place by the defray filter
     */
    Counts filter_alignment(Alignment& aln);
    
    /**
     * Filter the alignments available from the given stream, placing them on
     * standard output or in the appropriate file. Returns 0 on success, exit
     * code to use on error.
     *
     */
    int filter(istream* alignment_stream);
    
    /**
     * Look at either end of the given alignment, up to k bases in from the end.
     * See if that tail of the alignment is mapped such that another embedding
     * in the given graph can produce the same sequence as the sequence along
     * the embedding that the read actually has, and if so trim back the read.
     *
     * In the case of softclips, the aligned portion of the read is considered,
     * and if trimmign is required, the softclips are hard-clipped off.
     *
     * Returns true if the read had to be modified, and false otherwise.
     *
     * MUST NOT be called with a null index.
     */
    bool trim_ambiguous_ends(Alignment& alignment, int k);
    
private:

    /**
     * quick and dirty filter to see if removing reads that can slip around
     * and still map perfectly helps vg call.  returns true if at either
     * end of read sequence, at least k bases are repetitive, checking repeats
     * of up to size 2k
     */
    bool has_repeat(Alignment& aln, int k);
    
    /**
     * Trim only the end of the given alignment, leaving the start alone. Two
     * calls of this implement trim_ambiguous_ends above.
     */
    bool trim_ambiguous_end(Alignment& alignment, int k);
    
    /**
     * Return false if the read only follows edges in the graph, and true if
     * the read is split (or just incorrect) and takes edges not in the index.
     *
     * Throws an error if no graph is specified.
     */
    bool is_split(Alignment& alignment);
    
    /**
     * Based on the read name and paired-ness, compute the SAM-style QNAME and
     * use that and the configured sampling probability and seed in the
     * Samtools read sampling algorithm, to determine if the read should be
     * kept. Returns true if the read should stay, and false if it should be
     * removed. Always accepts or rejects paired reads together.
     */
    bool sample_read(const Alignment& read);
};
ostream& operator<<(ostream& os, const ReadFilter::Counts& counts);
}

#endif
