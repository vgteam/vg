#ifndef VG_READFILTER_HPP_INCLUDED
#define VG_READFILTER_HPP_INCLUDED

#include <vector>
#include <cstdlib>
#include <iostream>
#include <string>
#include <regex>
#include <fstream>
#include <sstream>

#include "vg.hpp"
#include "handle.hpp"
#include "IntervalTree.h"
#include "annotation.hpp"
#include "multipath_alignment_emitter.hpp"
#include <vg/io/alignment_emitter.hpp>
#include <vg/vg.pb.h>
#include <vg/io/stream.hpp>

#include <htslib/khash.h>

/** \file
 * Provides a way to filter and transform reads, implementing the bulk of the
 * `vg filter` command.
 *
 */
namespace vg{

using namespace std;

struct Counts;

template<typename Read>
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
    /// Read must contain at least one of these strings as a subsequence
    vector<string> subsequences;
    /// If a read has one of the features in this set as annotations, the read
    /// is filtered out.
    unordered_set<string> excluded_features;
    double min_secondary = numeric_limits<double>::lowest();
    double min_primary = numeric_limits<double>::lowest();
    /// Should we rescore each alignment with default parameters and no e.g.
    /// haplotype info?
    bool rescore = false;
    bool frac_score = false;
    bool sub_score = false;
    int max_overhang = numeric_limits<int>::max() / 2;
    int min_end_matches = numeric_limits<int>::min() / 2;
    bool verbose = false;
    double min_mapq = numeric_limits<double>::lowest();
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
    
    /// Filter to proper pairs
    bool only_proper_pairs = true;
    
    /// Filter to only mapped reads
    bool only_mapped = true;
    
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
    int min_base_quality = numeric_limits<int>::min() / 2;
    // minimum fraction of bases in reads that must have quality at least <min_base_quality>
    double min_base_quality_fraction = numeric_limits<double>::lowest();
      
    /**
     * Run all the filters on an alignment. The alignment may get modified in-place by the defray filter
     */
    Counts filter_alignment(Read& aln);
    
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
    bool trim_ambiguous_ends(Read& read, int k) const;
    
private:

    /**
     * quick and dirty filter to see if removing reads that can slip around
     * and still map perfectly helps vg call.  returns true if at either
     * end of read sequence, at least k bases are repetitive, checking repeats
     * of up to size 2k
     */
    bool has_repeat(const Read& read, int k) const;
    
    /**
     * Check if the alignment includes any aligned bases
     */
    bool is_mapped(const Read& read) const;
    
    /**
     * Trim only the end of the given alignment, leaving the start alone. Two
     * calls of this implement trim_ambiguous_ends above.
     */
    bool trim_ambiguous_end(Alignment& alignment, int k) const;
    
    /**
     * Return false if the read only follows edges in the graph, and true if
     * the read is split (or just incorrect) and takes edges not in the index.
     *
     * Throws an error if no graph is specified.
     */
    bool is_split(const Read& read) const;
    
    /**
     * Based on the read name and paired-ness, compute the SAM-style QNAME and
     * use that and the configured sampling probability and seed in the
     * Samtools read sampling algorithm, to determine if the read should be
     * kept. Returns true if the read should stay, and false if it should be
     * removed. Always accepts or rejects paired reads together.
     */
    bool sample_read(const Read& read) const;
    
    /**
     * Convert a multipath alignment to a single path
     */
    Alignment to_alignment(const MultipathAlignment& multipath_aln) const;
    
    /**
     * Get the score indicated by the params
     */
    double get_score(const Read& read) const;
    
    /**
     * Does the read name have one of the indicated prefixes?
     */
    bool matches_name(const Read& read) const;
    
    /**
     * Does the read match one of the excluded refpos contigs?
     */
    bool has_excluded_refpos(const Read& read) const;
    
    /**
     * Is the read annotated with any of the excluded features
     */
    bool has_excuded_feature(const Read& read) const;
    
    /**
     * Is the read a secondary alignments?
     */
    bool is_secondary(const Read& read) const;
    
    /**
     * How long are the read overhangs?
     */
    int get_overhang(const Read& read) const;
    
    /**
     * Internal helper for get_overhang
     */
    int alignment_overhang(const Alignment& aln) const;
    
    /**
     * What is the shortest run of end matches on either end?
     */
    int get_end_matches(const Read& read) const;
    
    /**
     * internal helper for get_end_matches
     */
    int alignment_end_matches(const Alignment& aln) const;
    
    /**
     * What is the read's mapping quality
     */
    int get_mapq(const Read& read) const;
    
    /**
     * What fraction of the base qualities are at least as large as the min base quality
     */
    double get_min_base_qual_fraction(const Read& read) const;
    
    /**
     * Is the read paired?
     */
    bool get_is_paired(const Read& read) const;
    
    /**
     * Is the read in a proper-mapped pair?
     */
    bool is_proper_pair(const Read& read) const;
    
    /**
     * Does the read contain at least one of the indicated sequences
     */
    bool contains_subsequence(const Read& read) const;
    
    /**
     * Write the read to stdout
     */
    void emit(Read& read);
    
    /**
     * Write a read pair to stdout
     */
    void emit(Read& read1, Read& read2);
    
    
    /// The twp specializations have different writing infrastructure
    unique_ptr<AlignmentEmitter> aln_emitter;
    unique_ptr<MultipathAlignmentEmitter> mp_aln_emitter;
    
    /// Helper function for filter
    void filter_internal(istream* in);
};

// Keep some basic counts for when verbose mode is enabled
struct Counts {
    // note: "last" must be kept as the final value in this enum
    enum FilterName { read = 0, wrong_name, wrong_refpos, excluded_feature, min_score, min_sec_score, max_overhang,
        min_end_matches, min_mapq, split, repeat, defray, defray_all, random, min_base_qual, subsequence, filtered,
        proper_pair, unmapped, last};
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
ostream& operator<<(ostream& os, const Counts& counts);


/**
 * Template implementations
 */

template <typename Read>
void ReadFilter<Read>::filter_internal(istream* in) {
    
    // keep counts of what's filtered to report (in verbose mode)
    vector<Counts> counts_vec(threads);
    
    function<void(Read&)> lambda = [&](Read& read) {
#ifdef debug
        cerr << "Encountered read named \"" << read.name() << "\" with " << read.sequence().size()
        << " bp sequence and " << read.quality().size() << " quality values" << endl;
#endif
        Counts read_counts = filter_alignment(read);
        counts_vec[omp_get_thread_num()] += read_counts;
        if ((read_counts.keep() != complement_filter) && write_output) {
            emit(read);
        }
    };
    
    function<void(Read&, Read&)> pair_lambda = [&](Read& read1, Read& read2) {
        Counts read_counts = filter_alignment(read1);
        read_counts += filter_alignment(read2);
        if (filter_on_all) {
            // Unless both reads were filtered out (total filtered count == 2), keep the read.
            read_counts.set_paired_all();
        } else {
            // Either read failing is sufficient to scuttle the pair.
            // So if we filter out one end for any reason, we filter out the other as well.
            read_counts.set_paired_any();
        }
        counts_vec[omp_get_thread_num()] += read_counts;
        if ((read_counts.keep() != complement_filter) && write_output) {
            emit(read1, read2);
        }
    };
    
    if (interleaved) {
        vg::io::for_each_interleaved_pair_parallel(*in, pair_lambda);
    } else {
        vg::io::for_each_parallel(*in, lambda);
    }
    
    if (verbose) {
        Counts& counts = counts_vec[0];
        for (int i = 1; i < counts_vec.size(); ++i) {
            counts += counts_vec[i];
        }
        cerr << counts;
    }
}

template<>
inline int ReadFilter<Alignment>::filter(istream* alignment_stream) {
    
    if(defray_length > 0 && graph == nullptr) {
        cerr << "HandleGraph (e.g. XG) required for end de-fraying" << endl;
        return 1;
    }
    
    if (write_output) {
        // Keep an AlignmentEmitter to multiplex output from multiple threads.
        aln_emitter = get_non_hts_alignment_emitter("-", "GAM", map<string, int64_t>(), get_thread_count());
    }
    
    filter_internal(alignment_stream);
    
    return 0;
}

template<>
inline int ReadFilter<MultipathAlignment>::filter(istream* alignment_stream) {
    
    if (defray_length > 0) {
        cerr << "Cannot defray multipath alignments" << endl;
        return 1;
    }
    if (!excluded_refpos_contigs.empty()) {
        cerr << "Cannot filter multipath alignments by ref pos" << endl;
        return 1;
    }
    
    if (write_output) {
        // Keep an AlignmentEmitter to multiplex output from multiple threads.
        mp_aln_emitter = unique_ptr<MultipathAlignmentEmitter>(new MultipathAlignmentEmitter("-", get_thread_count()));
    }
    
    filter_internal(alignment_stream);
    
    return 0;
}

template<>
inline void ReadFilter<Alignment>::emit(Alignment& aln) {
    aln_emitter->emit_single(std::move(aln));
}

template<>
inline void ReadFilter<Alignment>::emit(Alignment& aln1, Alignment& aln2) {
    aln_emitter->emit_pair(std::move(aln1), std::move(aln2));
}

template<>
inline void ReadFilter<MultipathAlignment>::emit(MultipathAlignment& mp_aln) {
    vector<multipath_alignment_t> emit_vec(1);
    from_proto_multipath_alignment(mp_aln, emit_vec.front());
    mp_aln_emitter->emit_singles(mp_aln.name(), std::move(emit_vec));
}

template<>
inline void ReadFilter<MultipathAlignment>::emit(MultipathAlignment& mp_aln1, MultipathAlignment& mp_aln2) {
    vector<pair<multipath_alignment_t, multipath_alignment_t>> emit_vec(1);
    from_proto_multipath_alignment(mp_aln1, emit_vec.front().first);
    from_proto_multipath_alignment(mp_aln2, emit_vec.front().second);
    mp_aln_emitter->emit_pairs(mp_aln1.name(), mp_aln2.name(), std::move(emit_vec));
}

template<typename Read>
Counts ReadFilter<Read>::filter_alignment(Read& read) {
    Counts counts;
    
    ++counts.counts[Counts::FilterName::read];
    bool keep = true;
    // filter (current) alignment
    if (!name_prefixes.empty()) {
        if (!matches_name(read)) {
            // There are prefixes and we don't match any, so drop the read.
            ++counts.counts[Counts::FilterName::wrong_name];
            keep = false;
        }
    }
    if ((keep || verbose) && !subsequences.empty()) {
        if (!contains_subsequence(read)) {
            // There are subsequences and we don't match any, so drop the read.
            ++counts.counts[Counts::FilterName::subsequence];
            keep = false;
        }
    }
    if ((keep || verbose) && only_proper_pairs) {
        if (!is_proper_pair(read)) {
            ++counts.counts[Counts::FilterName::proper_pair];
            keep = false;
        }
    }
    if ((keep || verbose) && !excluded_refpos_contigs.empty()) {
        if (has_excluded_refpos(read)) {
            ++counts.counts[Counts::FilterName::wrong_refpos];
            keep = false;
        }
    }
    if ((keep || verbose) && !excluded_features.empty()) {
        if (has_excuded_feature(read)) {
            ++counts.counts[Counts::FilterName::excluded_feature];
            keep = false;
        }
    }
    double score = get_score(read);
    bool secondary = is_secondary(read);
    if ((keep || verbose) && !secondary && score < min_primary) {
        ++counts.counts[Counts::FilterName::min_score];
        keep = false;
    }
    if ((keep || verbose) && secondary && score < min_secondary) {
        ++counts.counts[Counts::FilterName::min_sec_score];
        keep = false;
    }
    if ((keep || verbose) && max_overhang > 0) {
        if (get_overhang(read) > max_overhang) {
            ++counts.counts[Counts::FilterName::max_overhang];
            keep = false;
        }
    }
    if ((keep || verbose) && min_end_matches > 0) {
        if (get_end_matches(read) < min_end_matches) {
            ++counts.counts[Counts::FilterName::min_end_matches];
            keep = false;
        }
    }
    if ((keep || verbose) && min_mapq > 0) {
        if (get_mapq(read) < min_mapq) {
            ++counts.counts[Counts::FilterName::min_mapq];
            keep = false;
        }
    }
    if ((keep || verbose) && min_base_quality > 0 && min_base_quality_fraction > 0.0) {
        if (get_min_base_qual_fraction(read) < min_base_quality_fraction) {
            ++counts.counts[Counts::FilterName::min_base_qual];
            keep = false;
        }
    }
    if ((keep || verbose) && drop_split) {
        if (is_split(read)) {
            ++counts.counts[Counts::FilterName::split];
            keep = false;
        }
    }
    if ((keep || verbose) && repeat_size > 0) {
        if (has_repeat(read, repeat_size)) {
            ++counts.counts[Counts::FilterName::repeat];
            keep = false;
        }
    }
    if ((keep || verbose) && defray_length) {
        ++counts.counts[Counts::FilterName::defray];
        if (trim_ambiguous_ends(read, defray_length)) {
            // We keep these, because the alignments get modified.
            // Unless the *entire* read gets trimmed
            if (read.sequence().empty()) {
                keep = false;
                ++counts.counts[Counts::FilterName::defray_all];
            }
        }
    }
    if ((keep || verbose) && only_mapped) {
        if (!is_mapped(read)) {
            ++counts.counts[Counts::FilterName::unmapped];
            keep = false;
        }
    }
    if ((keep || verbose) && downsample_probability != 1.0) {
        if (!sample_read(read)) {
            ++counts.counts[Counts::FilterName::random];
            keep = false;
        }
    }
    
    if (!keep) {
        ++counts.counts[Counts::FilterName::filtered];
    }
    
    return counts;
}

template<typename Read>
Alignment ReadFilter<Read>::to_alignment(const MultipathAlignment& multipath_aln) const {
    multipath_alignment_t mp_aln;
    from_proto_multipath_alignment(multipath_aln, mp_aln);
    Alignment aln;
    optimal_alignment(mp_aln, aln);
    return aln;
}

template<>
inline double ReadFilter<Alignment>::get_score(const Alignment& aln) const {
    double score = (double)aln.score();
    double denom = aln.sequence().length();
    // toggle substitution score
    if (sub_score == true) {
        // hack in ident to replace old counting logic.
        score = aln.identity() * aln.sequence().length();
        assert(score <= denom);
    } else if (rescore == true) {
        // We need to recalculate the score with the base aligner always
        const static Aligner unadjusted;
        GSSWAligner* aligner = (GSSWAligner*)&unadjusted;
        // Also use the score
        score = aligner->score_contiguous_alignment(aln);
    }
    
    // toggle absolute or fractional score
    if (frac_score) {
        if (denom > 0.) {
            score /= denom;
        }
        else {
            assert(score == 0.);
        }
    }
    return score;
}

template<>
inline double ReadFilter<MultipathAlignment>::get_score(const MultipathAlignment& read) const {
    multipath_alignment_t mp_aln;
    from_proto_multipath_alignment(read, mp_aln);
    double score;
    if (!sub_score && !rescore) {
        score = optimal_alignment_score(mp_aln);
    }
    else {
        Alignment aln;
        optimal_alignment(mp_aln, aln);
        if (sub_score) {
            score = identity(aln.path());
        }
        else {
            const static Aligner unadjusted;
            GSSWAligner* aligner = (GSSWAligner*)&unadjusted;
            score = aligner->score_contiguous_alignment(aln);
        }
    }
    
    if (frac_score && read.sequence().size()) {
        score /= read.sequence().size();
    }
    return score;
}

template<typename Read>
bool ReadFilter<Read>::matches_name(const Read& aln) const {
    bool keep = true;
    // filter (current) alignment
    if (!name_prefixes.empty()) {
        // Make sure we match at least one name prefix
        
        bool found = false;
        
        // Do a binary search for the closest prefix and see if all of any prefix exists.
        // We assume the prefixes are sorted.
        size_t left_bound = 0;
        size_t left_match = 0;
        while (left_match < name_prefixes[left_bound].size() &&
               left_match < aln.name().size() &&
               name_prefixes[left_bound][left_match] == aln.name()[left_match]) {
            // Scan all the matches at the start
            left_match++;
        }
        
        size_t right_bound = name_prefixes.size() - 1;
        size_t right_match = 0;
        while (right_match < name_prefixes[right_bound].size() &&
               right_match < aln.name().size() &&
               name_prefixes[right_bound][right_match] == aln.name()[right_match]) {
            // Scan all the matches at the end
            right_match++;
        }
        
        if (left_match == name_prefixes[left_bound].size() || right_match == name_prefixes[right_bound].size()) {
            // We found a match already
            found = true;
        } else {
            while (left_bound + 1 < right_bound) {
                // Until we run out of unexamined prefixes, do binary search
                size_t center = (left_bound + right_bound) / 2;
                // No need to re-check any common prefix
                size_t center_match = min(left_match, right_match);
                
                while (center_match < name_prefixes[center].size() &&
                       center_match < aln.name().size() &&
                       name_prefixes[center][center_match] == aln.name()[center_match]) {
                    // Scan all the matches here
                    center_match++;
                }
                
                if (center_match == name_prefixes[center].size()) {
                    // We found a hit!
                    found = true;
                    break;
                }
                
                if (center_match == aln.name().size() ||
                    name_prefixes[center][center_match] > aln.name()[center_match]) {
                    // The match, if it exists, must be before us
                    right_bound = center;
                    right_match = center_match;
                }
                else {
                    // The match, if it exists, must be after us.
                    left_bound = center;
                    left_match = center_match;
                }
            }
        }
        
        if (!found) {
            // There are prefixes and we don't match any, so drop the read.
            keep = false;
        }
    }
    return keep;
}

template<>
inline bool ReadFilter<MultipathAlignment>::has_excluded_refpos(const MultipathAlignment& read) const {
    // TODO: multipath alignments don't record refpos
    return false;
}

template<>
inline bool ReadFilter<Alignment>::has_excluded_refpos(const Alignment& aln) const {
    bool found_match = false;
    if (!excluded_refpos_contigs.empty() && aln.refpos_size() != 0) {
        // We have refpos exclusion filters and a refpos is set.
        // We need to bang every refpos anme against every filter.
        
        for (auto& expression : excluded_refpos_contigs) {
            for (auto& refpos : aln.refpos()) {
                if (regex_search(refpos.name(), expression)) {
                    // We don't want this read because of this match
                    found_match = true;
                    break;
                }
            }
            if (found_match) {
                break;
            }
        }
    }
    return found_match;
}

template<typename Read>
bool ReadFilter<Read>::has_excuded_feature(const Read& read) const {
    bool found_match = false;
    vector<string> features(get_annotation<vector<string>>(read, "features"));
    
    for (auto& feature : features) {
        if (excluded_features.count(feature)) {
            // If the read has any banned features, fail it.
            found_match = true;
            break;
        }
    }
    return found_match;
}

template<>
inline bool ReadFilter<MultipathAlignment>::is_secondary(const MultipathAlignment& mp_aln) const {
    return get_annotation<bool>(mp_aln, "secondary");
}

template<>
inline bool ReadFilter<Alignment>::is_secondary(const Alignment& aln) const {
    return aln.is_secondary();
}

template<typename Read>
int ReadFilter<Read>::alignment_overhang(const Alignment& aln) const {
    int overhang = 0;
    if (aln.path().mapping_size() > 0) {
        const auto& left_mapping = aln.path().mapping(0);
        if (left_mapping.edit_size() > 0) {
            overhang = left_mapping.edit(0).to_length() - left_mapping.edit(0).from_length();
        }
        const auto& right_mapping = aln.path().mapping(aln.path().mapping_size() - 1);
        if (right_mapping.edit_size() > 0) {
            const auto& edit = right_mapping.edit(right_mapping.edit_size() - 1);
            overhang = max(overhang, edit.to_length() - edit.from_length());
        }
    }
    else {
        overhang = aln.sequence().size();
    }
    return overhang;
}

template<>
inline int ReadFilter<Alignment>::get_overhang(const Alignment& aln) const {
    return alignment_overhang(aln);
}

template<>
inline int ReadFilter<MultipathAlignment>::get_overhang(const MultipathAlignment& mp_aln) const {
    return alignment_overhang(to_alignment(mp_aln));
}

template<typename Read>
int ReadFilter<Read>::alignment_end_matches(const Alignment& aln) const {
    // compute end matches.
    int left_end_matches = 0;
    // from the left
    for (int i = 0; i < aln.path().mapping_size() && left_end_matches < min_end_matches; ++i) {
        for (int j = 0; j < aln.path().mapping(i).edit_size() && left_end_matches < min_end_matches; ++j) {
            const Edit& edit = aln.path().mapping(i).edit(j);
            if (edit.from_length() == edit.to_length() && edit.sequence().empty()) {
                left_end_matches += edit.to_length();
            } else {
                i = aln.path().mapping_size();
                break;
            }
        }
    }
    int right_end_matches = 0;
    // from the right
    for (int i = aln.path().mapping_size() - 1; i >= 0 && right_end_matches < min_end_matches; --i) {
        for (int j = aln.path().mapping(i).edit_size() - 1; j >= 0 && right_end_matches < min_end_matches; --j) {
            const Edit& edit = aln.path().mapping(i).edit(j);
            if (edit.from_length() == edit.to_length() && edit.sequence().empty()) {
                right_end_matches += edit.to_length();
            } else {
                i = -1;
                break;
            }
        }
    }
    return min(left_end_matches, right_end_matches);
}

template<>
inline int ReadFilter<MultipathAlignment>::get_end_matches(const MultipathAlignment& read) const {
    return alignment_end_matches(to_alignment(read));
}

template<>
inline int ReadFilter<Alignment>::get_end_matches(const Alignment& read) const {
    return alignment_end_matches(read);
}

template<typename Read>
int ReadFilter<Read>::get_mapq(const Read& read) const {
    return read.mapping_quality();
}

template<typename Read>
double ReadFilter<Read>::get_min_base_qual_fraction(const Read& read) const {
    int mq_count = 0;
    const string& base_qualities = read.quality();
    for (int i = 0; i < base_qualities.length(); ++i) {
        if (short(base_qualities[i]) >= min_base_quality) {
            ++mq_count;
        }
    }
    return (double)mq_count / (double)base_qualities.size() < min_base_quality_fraction;
}

template<>
inline bool ReadFilter<Alignment>::is_split(const Alignment& alignment) const {
    if(graph == nullptr) {
        // Can't tell if the read is split.
        throw runtime_error("HandleGraph (e.g. XG) required to check for split reads");
    }
    
    handle_t prev;
    for(size_t i = 0; i + 1 < alignment.path().mapping_size(); i++) {
        if (i == 0) {
            const auto& pos = alignment.path().mapping(i).position();
            prev = graph->get_handle(pos.node_id(), pos.is_reverse());
        }
        const auto& pos = alignment.path().mapping(i + 1).position();
        handle_t here = graph->get_handle(pos.node_id(), pos.is_reverse());
        
        // Can we find the same articulation of the edge as the alignment uses
        
        if(!graph->has_edge(prev, here)) {
            // We found a skip!
            if(verbose) {
                cerr << "Warning: read " << alignment.name() << " has an unknown edge "
                << graph->get_id(prev) << (graph->get_is_reverse(prev) ? "-" : "+") << " -> " << graph->get_id(here) << (graph->get_is_reverse(here) ? "-" : "+")
                << ". Removing!" << endl;
            }
            return true;
        }
        
        prev = here;
    }
    
    // No wandering jumps between nodes found
    return false;
}

template<>
inline bool ReadFilter<MultipathAlignment>::is_split(const MultipathAlignment& mp_aln) const {
    bool found_connection = false;
    for (const auto& subpath : mp_aln.subpath()) {
        if (subpath.connection_size() != 0) {
            found_connection = true;
            break;
        }
    }
    return found_connection;
}

// quick and dirty filter to see if removing reads that can slip around
// and still map perfectly helps vg call.  returns true if at either
// end of read sequence, at least k bases are repetitive, checking repeats
// of up to size 2k
template<typename Read>
bool ReadFilter<Read>::has_repeat(const Read& read, int k) const {
    if (k == 0) {
        return false;
    }
    const string& s = read.sequence();
    for (int i = 1; i <= 2 * k; ++i) {
        int covered = 0;
        bool ffound = true;
        bool bfound = true;
        for (int j = 1; (ffound || bfound) && (j + 1) * i < s.length(); ++j) {
            ffound = ffound && s.substr(0, i) == s.substr(j * i, i);
            bfound = bfound && s.substr(s.length() - i, i) == s.substr(s.length() - i - j * i, i);
            if (ffound || bfound) {
                covered += i;
            }
        }
        if (covered >= k) {
            return true;
        }
    }
    return false;
}

template<>
inline bool ReadFilter<Alignment>::trim_ambiguous_ends(Alignment& alignment, int k) const {
    assert(graph != nullptr);
    
    // Define a way to get node length, for flipping alignments
    function<int64_t(id_t)> get_node_length = [&](id_t node) {
        return graph->get_length(graph->get_handle(node));
    };
    
    // Because we need to flip the alignment, make sure it is new-style and
    // doesn't have any Mappings with no Edits.
    for(size_t i = 0; i < alignment.path().mapping_size(); i++) {
        if(alignment.path().mapping(i).edit_size() == 0) {
            // Complain!
            throw runtime_error("Found mapping with no edits in " + pb2json(alignment));
        }
    }
    
    // TODO: we're going to flip the alignment twice! This is a waste of time!
    // Define some kind of oriented view or something, or just two duplicated
    // trimming functions, so we can just trim once without flipping.
    
    // Trim the end
    bool end_changed = trim_ambiguous_end(alignment, k);
    // Flip and trim the start
    
    Alignment flipped = reverse_complement_alignment(alignment, get_node_length);
    
    if(trim_ambiguous_end(flipped, k)) {
        // The start needed trimming
        
        // Flip the trimmed flipped alignment back
        alignment = reverse_complement_alignment(flipped, get_node_length);
        // We definitely changed something
        return true;
    }
    
    // We maybe changed something
    return end_changed;
}

template<>
inline bool ReadFilter<MultipathAlignment>::trim_ambiguous_ends(MultipathAlignment& aln, int k) const {
    // TODO: apply this filter to mp alns?
    return false;
}

template<typename Read>
bool ReadFilter<Read>::trim_ambiguous_end(Alignment& alignment, int k) const {
    // What mapping in the alignment is the leftmost one starting in the last k
    // bases? (Except when that would be the first mapping, we use the second.)
    // Start out with it set to the past-the-end value.
    size_t trim_start_mapping = alignment.path().mapping_size();
    
    // How many real non-softclip bases have we seen reading in from the end of
    // the read?
    size_t real_base_count = 0;
    // How many softclip bases have we seen in from the end of the read?
    size_t softclip_base_count = 0;
    for(size_t i = alignment.path().mapping_size() - 1; i != -1 && i != 0; i--) {
        // Scan in from the end of the read.
        
        auto* mapping = alignment.mutable_path()->mutable_mapping(i);
        
        // We should always have edits in our mappings.
        assert(mapping->edit_size() > 0);
        
        for(int j = mapping->edit_size() - 1; j != -1; j--) {
            // Visit every edit in the mapping
            auto& edit = mapping->edit(j);
            
            
            if(real_base_count == 0 && edit.from_length() == 0) {
                // This is a trailing insert. Put it as a softclip
                softclip_base_count += edit.to_length();
            } else {
                // This is some other kind of thing. Record it as real bases.
                real_base_count += edit.to_length();
            }
        }
        
        if(real_base_count <= k) {
            // This mapping starts fewer than k non-softclipped alignment
            // bases from the end of the read.
            trim_start_mapping = i;
        } else {
            // This mapping starts more than k in from the end. So the
            // previous one, if we had one, must be the right one.
            break;
        }
    }
    
    if(trim_start_mapping == alignment.path().mapping_size()) {
        // No mapping was found that starts within the last k non-softclipped
        // bases. So there's nothing to do.
        return false;
    }
    
    if(real_base_count == 0) {
        // We have an anchoring mapping, but all the mappings we could trim are
        // softclips, so there's no point. TODO: will we ever get softclips
        // placed as the only thing on a node?
        return false;
    }
    
    // Which is the last assumed-non-ambiguous mapping from which we can anchor
    // our search?
    size_t root_mapping = trim_start_mapping - 1;
    
    // What's the sequence, including that root node, that we are looking for?
    // We need the sequence of the nodes, rather than the read's sequence,
    // because you can still be ambiguous even if you have a SNP on top of the
    // ambiguous thing.
    
    // We need to ignore all the offsets and from_lengths, except for the from
    // length on the last node to let us know if we end early. It's sort of
    // nonsense to have offsets and non-full from_lengths on internal mappings,
    // and everything is easiest if we use the full length sequence of the root
    // node.
    stringstream target_sequence_stream;
    for(size_t i = root_mapping; i < alignment.path().mapping_size(); i++) {
        // Collect the appropriately oriented from sequence from each mapping
        auto& mapping = alignment.path().mapping(i);
        handle_t handle = graph->get_handle(mapping.position().node_id(),
                                            mapping.position().is_reverse());
        string sequence = graph->get_sequence(handle);
        
        if(i == root_mapping) {
            // Use the full length of the node and ignore any offset
            target_sequence_stream << sequence;
        } else {
            // Use the offset plus the total from_length of all the
            // edits (in case we're the last node and ending early). We made
            // sure all non-root nodes had edits earlier.
            
            size_t from_length = mapping.position().offset();
            for(size_t j = 0; j < mapping.edit_size(); j++) {
                from_length += mapping.edit(j).from_length();
            }
            
            // Put in the sequence that the mapping visits
            target_sequence_stream << sequence.substr(0, from_length);
        }
    }
    string target_sequence = target_sequence_stream.str();
    
#ifdef debug
#pragma omp critical(cerr)
    cerr << "Need to look for " << target_sequence << " right of mapping " << root_mapping << endl;
#endif
    
    // We're not going to recurse hundreds of nodes deep, so we can use the real
    // stack and a real recursive function.
    
    // Do the DFS into the given node, after already having matched the given
    // number of bases of the target sequence. See if you can match any more
    // bases of the target sequence.
    
    // Return the total number of leaves in all subtrees that match the full
    // target sequence, and the depth in bases of the shallowest point at which
    // multiple subtrees with full lenght matches are unified.
    
    // We keep a maximum number of visited nodes here, just to prevent recursion
    // from going on forever in worst-case-type graphs
    size_t dfs_visit_count = 0;
    function<pair<size_t, size_t>(const handle_t&, size_t)> do_dfs =
    [&](const handle_t& handle, size_t matched) -> pair<size_t, size_t> {
        
        ++dfs_visit_count;
        
        // Grab the node sequence and match more of the target sequence.
        string node_sequence = graph->get_sequence(handle);
        
#ifdef debug
#pragma omp critical(cerr)
        cerr << "Node " << graph->get_id(handle) <<  " " << (graph->get_is_reverse(handle) ? "rev" : "fwd") << ": "
        << node_sequence << " at offset " << matched << " in " << target_sequence << endl;
#endif
        
        // Now count up the new matches between this node and the target sequence.
        size_t new_matches;
        for(
            // Start with no matches
            new_matches = 0;
            // Keep going as long as we're inside both strings and haven't had a mismatch
            new_matches < node_sequence.size() &&
            matched + new_matches < target_sequence.size() &&
            node_sequence[new_matches] == target_sequence[matched + new_matches];
            // Count up all the matches we find
            new_matches++
            );
        
        if(matched + new_matches == target_sequence.size()) {
            // We found a tail end of a complete match of the target sequence
            // on this node.
            
#ifdef debug
#pragma omp critical(cerr)
            cerr << "Node " << node_id << " is a matching leaf" << endl;
#endif
            
            // Return one match and unification at full length (i.e. nothing can
            // be discarded).
            return make_pair(1, target_sequence.size());
        }
        
        if(new_matches < node_sequence.size()) {
            // We didn't make it all the way through this node, nor did we
            // finish the target sequence; there's a mismatch between the node
            // and the target sequence.
            
#ifdef debug
#pragma omp critical(cerr)
            cerr << "Node " << node_id << " has a mismatch" << endl;
#endif
            
            // If we mismatch, return 0 matches and unification at full length.
            return make_pair(0, target_sequence.size());
        }
        
        // If we get through the whole node sequence without mismatching or
        // running out of target sequence, keep going.
        
#ifdef debug
#pragma omp critical(cerr)
        cerr << "Node " << graph->get_id(handle) << " has " << new_matches << " internal new matches" << endl;
#endif
        
        // We're going to call all the children and collect the results, and
        // then aggregate them. It might be slightly faster to aggregate while
        // calling, but that might be less clear.
        vector<pair<size_t, size_t>> child_results;
        
        graph->follow_edges(handle, false, [&](const handle_t& next) {
            if (dfs_visit_count < defray_count) {
                child_results.push_back(do_dfs(next, matched + node_sequence.size()));
                return true;
            }
            else {
#ifdef debug
#pragma omp critical(cerr)
                cerr << "Aborting read filter DFS at node " << graph->get_id(next) << " after " << dfs_visit_count << " visited" << endl;
#endif
                return false;
            }
        });
        
        // Sum up the total leaf matches, which will be our leaf match count.
        size_t total_leaf_matches = 0;
        // If we don't find multiple children with leaf matches, report
        // unification at the min unification depth of any subtree (and there
        // will only be one that isn't at full length).
        size_t children_with_leaf_matches = 0;
        size_t unification_depth = target_sequence.size();
        
        for(auto& result : child_results) {
            total_leaf_matches += result.first;
            if(result.first > 0) {
                children_with_leaf_matches++;
            }
            unification_depth = min(unification_depth, result.second);
        }
        if(children_with_leaf_matches > 1) {
            // If multiple children have nonzero leaf match counts, report
            // unification at the end of this node.
            unification_depth = matched + node_sequence.size();
        }
        
        return make_pair(total_leaf_matches, unification_depth);
    };
    
    // Search from the root mapping's node looking right in its orientation in
    // the mapping
    const auto& root_pos = alignment.path().mapping(root_mapping).position();
    auto result = do_dfs(graph->get_handle(root_pos.node_id(), root_pos.is_reverse()), 0);
    
#ifdef debug
#pragma omp critical(cerr)
    cerr << "Found " << result.first << " matching leaves with closest unification at " << result.second << endl;
#endif
    
    // We keep this much of the target sequence.
    size_t target_sequence_to_keep = result.second;
    
    if(target_sequence_to_keep == target_sequence.size()) {
        // Nothing to trim!
        return false;
    }
    
    // Figure out how many mappings we need to keep from the root in order to
    // get that much sequence; we know it is either full length or at a mapping
    // boundary. We handle the root special because it's always full length and
    // we have to cut after its end.
    size_t kept_sequence_accounted_for = graph->get_length(graph->get_handle(alignment.path().mapping(root_mapping).position().node_id()));
    size_t first_mapping_to_drop;
    for(first_mapping_to_drop = root_mapping + 1;
        first_mapping_to_drop < alignment.path().mapping_size();
        first_mapping_to_drop++) {
        // Consider starting dropping at each mapping after the root.
        if(kept_sequence_accounted_for == target_sequence_to_keep) {
            // OK this mapping really is the first one to drop.
            break;
        } else {
            // Keep going. Account for the sequence from this mapping.
            auto& mapping = alignment.path().mapping(first_mapping_to_drop);
            
            // We know it's not the root mapping, and it can't be the non-full-
            // length end mapping (because we would have kept the full length
            // target sequence and not had to cut). So assume full node is used.
            kept_sequence_accounted_for += graph->get_length(graph->get_handle(mapping.position().node_id()));
        }
    }
    
    // OK we know the first mapping to drop. We need to work out the to_size,
    // including all softclips, from there to the end, so we know how much to
    // trim off of the sequence and quality.
    size_t to_length_to_remove = 0;
    for(size_t i = first_mapping_to_drop; i < alignment.path().mapping_size(); i++) {
        // Go through all the mappings
        auto& mapping = alignment.path().mapping(i);
        for(size_t j = 0; j < mapping.edit_size(); j++) {
            // Add up the to_length of all the edits
            to_length_to_remove += mapping.edit(j).to_length();
        }
    }
    
#ifdef debug
    cerr << "Want to trim " << alignment.sequence().size() << " bp sequence and " << alignment.quality().size()
    << " quality values to remove " << to_length_to_remove << endl;
#endif
    
    // Make sure we have at least enough to trim.
    // Note that we allow the entire alignment to be trimmed away!
    assert(alignment.sequence().size() >= to_length_to_remove);
    assert(alignment.quality().empty() || alignment.quality().size() >= to_length_to_remove);
    // And that we made sence before trimming
    assert(alignment.quality().empty() || alignment.quality().size() == alignment.sequence().size());
    
    // Trim sequence
    alignment.set_sequence(alignment.sequence().substr(0, alignment.sequence().size() - to_length_to_remove));
    
    // Trim quality
    if(!alignment.quality().empty()) {
        alignment.set_quality(alignment.quality().substr(0, alignment.quality().size() - to_length_to_remove));
    }
    
    // Now we can discard the extra mappings
    size_t to_delete = alignment.path().mapping_size() - first_mapping_to_drop;
    alignment.mutable_path()->mutable_mapping()->DeleteSubrange(first_mapping_to_drop, to_delete);
    
    // Now the alignment is fixed!
    return true;
}

template<>
inline bool ReadFilter<Alignment>::get_is_paired(const Alignment& aln) const {
    return (!aln.fragment_prev().name().empty() || aln.fragment_prev().path().mapping_size() != 0 ||
            !aln.fragment_next().name().empty() || aln.fragment_next().path().mapping_size() != 0);
}

template<>
inline bool ReadFilter<MultipathAlignment>::get_is_paired(const MultipathAlignment& mp_aln) const {
    return !mp_aln.paired_read_name().empty();
}

template<typename Read>
bool ReadFilter<Read>::is_proper_pair(const Read& read) const {
    if (!has_annotation(read, "proper_pair")) {
        return false;
    }
    return get_annotation<bool>(read, "proper_pair");
}
    
template<typename Read>
bool ReadFilter<Read>::contains_subsequence(const Read& read) const {
    bool found = false;
    for (const string& seq : subsequences) {
        if (read.sequence().find(seq) != string::npos) {
            found = true;
            break;
        }
    }
    return found;
}

template<typename Read>
bool ReadFilter<Read>::sample_read(const Read& read) const {
    // Decide if the alignment is paired.
    // It is paired if fragment_next or fragment_prev point to something.
    bool is_paired = get_is_paired(read);
    
    // Compute the QNAME that samtools would use
    string qname;
    if (is_paired) {
        // Strip pair end identifiers like _1 or /2 that vg uses at the end of the name.
        qname = regex_replace(read.name(), regex("[/_][12]$"), "");
    } else {
        // Any _1 in the name is part of the actual read name.
        qname = read.name();
    }
    
    // Now treat it as samtools would.
    // See https://github.com/samtools/samtools/blob/60138c42cf04c5c473dc151f3b9ca7530286fb1b/sam_view.c#L101-L104
    
    // Hash that with __ac_X31_hash_string from htslib and XOR against the seed mask
    auto masked_hash = __ac_X31_hash_string(qname.c_str()) ^ downsample_seed_mask;
    
    // Hash that again with __ac_Wang_hash from htslib, apparently to mix the bits.
    uint32_t mixed_hash = __ac_Wang_hash(masked_hash);
    
    // Take the low 24 bits and compute a double from 0 to 1
    const int32_t LOW_24_BITS = 0xffffff;
    double sample = ((double)(mixed_hash & LOW_24_BITS)) / (LOW_24_BITS + 1);
    
    // If the result is >= the portion to downsample to, discard the read.
    // Otherwise, keep it.
    return (sample < downsample_probability);
}
    
}

#endif
