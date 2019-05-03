#ifndef VG_POSITION_HPP_INCLUDED
#define VG_POSITION_HPP_INCLUDED

#include <vg/vg.pb.h>
#include "types.hpp"
#include "utility.hpp"
#include "json2pb.h"
#include <gcsa/gcsa.h>

/** \file 
 * Functions for working with Positions and `pos_t`s.
 */

namespace vg {

using namespace std;
/// Reverse a Position and get a Position at the same **point between bases**, going the other direction.
/// To get a Position to the same *base*, subtract 1 from the resulting offset.
Position reverse(const Position& pos, size_t node_length);
/// Convert a Position to a (much smaller) pos_t.
pos_t make_pos_t(const Position& pos);
/// Create a pos_t from a gcsa node
pos_t make_pos_t(gcsa::node_type node);
/// Convert a pos_t to a Position.
Position make_position(const pos_t& pos);
/// Create a Position from a Node ID, an orientation flag, and an offset along that strand of the node.
Position make_position(id_t id, bool is_rev, off_t off);
/// Make a Position from a gcsa node
Position make_position(gcsa::node_type node);
/// Find the min distance in the path offsets where the path orientation is the same and different
pair<int64_t, int64_t> min_oriented_distances(const map<string, vector<pair<size_t, bool> > >& path_offsets1,
                                              const map<string, vector<pair<size_t, bool> > >& path_offsets2);

}

#endif
