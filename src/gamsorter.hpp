#ifndef VG_GAMSORTER_HPP_INCLUDED
#define VG_GAMSORTER_HPP_INCLUDED

#include "vg.pb.h"
#include "stream.hpp"
#include <string>
#include <queue>
#include <sstream>
#include <functional>
#include <algorithm>
#include <iostream>
#include <set>
#include <vector>
#include <unordered_map>
#include <tuple>

/**
 * \file gamsorter.hpp
 * GAM sorting tools.
 */
using namespace std;
namespace vg
{


/// Provides the ability to sort a GAM, either "dumbly" (in memory), or
/// "streaming" into temporary files. Paired alignments are not necessarily
/// going to end up next to each other, so if sorting by position make sure to
/// set the position cross-references first if you want to be able to find
/// them.
class GAMSorter
{

  
  public:

    //////////////////
    // Main entry points
    //////////////////
    
    /// Sort a stream of GAM-format data, using temporary files.
    void stream_sort(istream& gam_in, ostream& gam_out) const;
    
    /// Sort a stream of GAM-format data, loading it all into memory and doing
    /// a single giant sort operation.
    void dumb_sort(istream& gam_in, ostream& gam_out) const;
    
    //////////////////
    // Supporting API
    //////////////////

    /// Sort a vector of alignments, in place.
    void sort(vector<Alignment>& alns) const;

    /// Return true if out of Alignments a and b, alignment a must come before alignment b, and false otherwise.
    bool less_than(const Alignment& a, const Alignment& b) const;

    /// Determine the minumum Position visited by an Alignment. The minumum
    /// Position is the lowest node ID visited by the alignment, with the
    /// lowest offset visited on that node ID as the offset, and the
    /// orientation set to false if the forward strand is visited, and true if
    /// only the reverse strand is visited.
    Position get_min_position(const Alignment& aln) const;

    /// Determine the minimum position visited by a Path, as for an Alignment.
    Position get_min_position(const Path& path) const;

    /// Return True if the given Position values are equal, and false otherwise.
    bool equal_to(const Position& a, const Position& b) const;

    /// Return True if position A is less than position B in our sort, and false otherwise.
    /// Position order is defined first by node ID, then by strand (forward first), and then by offset within the strand.
    /// We can't sort by actual base on the forward strand, because we need to be able to sort without knowing the graph's node lengths.
    bool less_than(const Position& a, const Position& b) const;

    /// Return True if position A is greater than position B in our sort, and false otherwise.
    bool greater_than(const Position& a, const Position& b) const;

  private:
    int max_buf_size = 10; //20000;
};
}
#endif
