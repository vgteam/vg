/**
 * \file splice_region.hpp
 *
 * Defines an object that finds candidate splice sites within the area around a given position
 *
 */
#ifndef VG_SPLICE_REGION_HPP_INCLUDED
#define VG_SPLICE_REGION_HPP_INCLUDED

#include "dinucleotide_machine.hpp"
#include "incremental_subgraph.hpp"

namespace vg {

using namespace std;

/*
 * Object that identifies possible splice sites in a small region of
 * the graph and answers queries about them.
 */
class SpliceRegion {
public:
    
    // note: search distance is relative to the *start* of the motif, whereas
    // the returned distances are relative to the *end* of the motif
    SpliceRegion(const pos_t& seed_pos, bool search_left, int64_t search_dist,
                 const HandleGraph& graph,
                 const DinucleotideMachine& dinuc_machine,
                 const vector<string>& splice_motifs);
    SpliceRegion() = default;
    ~SpliceRegion() = default;
    
    // returns the locations in the graph just past the location of a splice
    // motif (in the direction of the intron) and the distance of that location
    // from the search position. crashes if given a motif that was not provided
    // to the constructor.
    const vector<pair<pos_t, int64_t>>& candidate_splice_sites(const string& motif) const;
    
private:

    //
    unordered_map<string, vector<pair<pos_t, int64_t>>> motif_matches;
    
};

}

#endif // VG_SPLICE_REGION_HPP_INCLUDED
