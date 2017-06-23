#ifndef VG_READFILTER_HPP
#define VG_READFILTER_HPP

#include <vector>
#include <cstdlib>
#include <iostream>
#include <string>
#include "vg.hpp"
#include "xg.hpp"
#include "vg.pb.h"

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
    double min_secondary = 0.;
    double min_primary = 0.;
    bool frac_score = false;
    bool sub_score = false;
    int max_overhang = 99999;
    int min_end_matches = 0;
    int context_size = 0;
    bool verbose = false;
    double min_mapq = 0.;
    int repeat_size = 0;
    // How far in from the end should we look for ambiguous end alignment to
    // clip off?
    int defray_length = 0;
    // Limit defray recursion to visit this many nodes
    int defray_count = 99999;
    // Should we drop split reads that follow edges not in the graph?
    bool drop_split = false;
    // default to 1 thread (as opposed to all)
    int threads = 1;

    // Keep some basic counts for when verbose mode is enabled
    struct Counts {
        vector<size_t> read;
        vector<size_t> filtered;
        vector<size_t> min_score;
        vector<size_t> max_overhang;
        vector<size_t> min_end_matches;
        vector<size_t> min_mapq;
        vector<size_t> split;
        vector<size_t> repeat;
        vector<size_t> defray;
        Counts() : read(2, 0), filtered(2, 0), min_score(2, 0),
                   max_overhang(2, 0), min_end_matches(2, 0), min_mapq(2, 0),
                   split(2, 0), repeat(2, 0), defray(2, 0) {}
        Counts& operator+=(const Counts& other) {
            for (int i = 0; i < 2; ++i) {
                read[i] += other.read[i];
                filtered[i] += other.filtered[i];
                min_score[i] += other.min_score[i];
                max_overhang[i] += other.max_overhang[i];
                min_end_matches[i] += other.min_end_matches[i];
                min_mapq[i] += other.min_mapq[i];
                split[i] += other.split[i];
                repeat[i] += other.repeat[i];
                defray[i] += other.defray[i];
            }
            return *this;
        }
    };
    
    // Extra filename things we need for chunking. TODO: refactor that somehow
    // to maybe be a different class?
    string regions_file;
    string outbase;
    bool append_regions = false;
    
    /**
     * Filter the alignments available from the given stream, placing them on
     * standard output or in the appropriate file. Returns 0 on success, exit
     * code to use on error.
     *
     * If an XG index is required, use the specified one. If one is required and
     * not provided, the function will complain and return nonzero.
     *
     * TODO: Refactor to be less CLI-aware and more modular-y.
     */
    int filter(istream* alignment_stream, xg::XG* xindex = nullptr);
    
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
    bool trim_ambiguous_ends(xg::XG* index, Alignment& alignment, int k);
    
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
    bool trim_ambiguous_end(xg::XG* index, Alignment& alignment, int k);
    
    /**
     * Return false if the read only follows edges in the xg index, and true if
     * the read is split (or just incorrect) and takes edges not in the index.
     *
     * Throws an error if no XG index is specified.
     */
    bool is_split(xg::XG* index, Alignment& alignment);
    
};
}

#endif
