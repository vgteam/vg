/**
 * \file splice_region.hpp
 *
 * Defines SpliceRegion and some other splice related functions
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
 * Object that identifies possible splice sites in a small region of
 * the graph and answers queries about them.
 */
class SpliceRegion {
public:
    
    SpliceRegion(const pos_t& seed_pos, bool search_left, int64_t search_dist,
                 const HandleGraph& graph,
                 const DinucleotideMachine& dinuc_machine,
                 const vector<string>& splice_motifs);
    SpliceRegion() = default;
    ~SpliceRegion() = default;
    
    // returns the locations in the of a splice motif and the distance of those location
    // from the search position. crashes if given a motif that was not provided
    // to the constructor.
    const vector<pair<pos_t, int64_t>>& candidate_splice_sites(const string& motif) const;
    
    pos_t motif_pos_to_splice_pos();
    
private:

    unordered_map<string, vector<pair<pos_t, int64_t>>> motif_matches;
    
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
