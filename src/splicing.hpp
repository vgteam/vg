/**
 * \file splicing.hpp
 *
 * Defines splicing related objects and functions
 *
 */
#ifndef VG_SPLICE_REGION_HPP_INCLUDED
#define VG_SPLICE_REGION_HPP_INCLUDED

#include "dinucleotide_machine.hpp"
#include "incremental_subgraph.hpp"
#include "aligner.hpp"

namespace vg {

using namespace std;

/*
 * Object that represents the acceptable splice motifs and their
 * scores
 */
class SpliceMotifs {
public:
    // Defaults the three canonical human splicing motifs with the frequencies
    // reported in Burset, et al. (2000):
    // - GT-AG: 0.9924
    // - GC-AG: 0.0069
    // - AT-AC: 0.0005
    SpliceMotifs(const GSSWAligner& scorer);
    
    // Construct with triples of (5' dinucleotide, 3' dinucleotide, frequency)
    // and possibly non-default scoring
    SpliceMotifs(const vector<tuple<string, string, double>>& motifs,
                 const GSSWAligner& scorer);
    
    // the number of splicing motifs
    size_t size() const;
    // the dinucleotide motif on one side in the order that the nucleotides
    // are encountered when traversing into the intron
    const string& oriented_motif(size_t motif_num, bool left_side) const;
    // the score associated with a splicing motif
    int32_t score(size_t motif_num) const;
    // must be called if scoring parameters are changed
    void update_scoring(const GSSWAligner& scorer);
private:
    
    // internal function for the constructor
    void init(const vector<tuple<string, string, double>>& motifs,
              const GSSWAligner& scorer);
    
    // the splice table stored in its most useful form
    vector<tuple<string, string, int32_t>> data;
    // the original input data
    vector<tuple<string, string, double>> unaltered_data;
};


/*
 * Object that identifies possible splice sites in a small region of
 * the graph and answers queries about them.
 */
class SpliceRegion {
public:
    
    SpliceRegion(const pos_t& seed_pos, bool search_left, int64_t search_dist,
                 const HandleGraph& graph,
                 const DinucleotideMachine& dinuc_machine,
                 const SpliceMotifs& splice_motifs);
    SpliceRegion() = default;
    ~SpliceRegion() = default;
    
    // returns the locations in the of a splice motif and the distance of those location
    // from the search position. crashes if given a motif that was not provided
    // to the constructor.
    const vector<pair<pos_t, int64_t>>& candidate_splice_sites(size_t motif_num) const;
    
private:

    vector<vector<pair<pos_t, int64_t>>> motif_matches;
    
};

// return the position of the base when trimming a given length from either the start
// or the end of the alignment. softclips do not contributed to the total, and extra
// will be trimmed to avoid splitting an insertion. also returns the total amount of
// sequence trimmed (including softclip) and the score of the sub-alignment that would
// be removed.
tuple<pos_t, int64_t, int32_t> trimmed_end(const Alignment& aln, int64_t len, bool from_end,
                                           const HandleGraph& graph, const GSSWAligner& aligner);

}

#endif // VG_SPLICE_REGION_HPP_INCLUDED
