#ifndef VG_POSITION_HPP_INCLUDED
#define VG_POSITION_HPP_INCLUDED

#include "vg.pb.h"
#include "types.hpp"
#include "utility.hpp"
#include "json2pb.h"
#include <gcsa/gcsa.h>
#include <iostream>

/** \file 
 * Functions for working with Positions and `pos_t`s.
 */

namespace vg {

using namespace std;

/// Return true if a pos_t is unset.
bool is_empty(const pos_t& pos);
/// Extract the id of the node a pos_t is on.
id_t id(const pos_t& pos);
/// Return true if a pos_t is on the reverse strand of its node.
bool is_rev(const pos_t& pos);
/// Get the offset from a pos_t.
off_t offset(const pos_t& pos);
/// Get a reference to the Node ID of a pos_t.
id_t& get_id(pos_t& pos);
/// Get a reference to the reverse flag of a pos_t.
bool& get_is_rev(pos_t& pos);
/// Get a reference to the offset field of a pos_t.
off_t& get_offset(pos_t& pos);
/// Reverse a pos_t and get a pos_t at the same base, going the other direction.
pos_t reverse(const pos_t& pos, size_t node_length);
/// Reverse a Position and get a Position at the same base, going the orther direction.
Position reverse(const Position& pos, size_t node_length);
/// Print a pos_t to a stream.
ostream& operator<<(ostream& out, const pos_t& pos);
/// Convert a Position to a (much smaller) pos_t.
pos_t make_pos_t(const Position& pos);
/// Create a pos_t from a Node ID, an orientation flag, and an offset.
pos_t make_pos_t(id_t id, bool is_rev, off_t off);
/// Create a pos_t from a gcsa node
pos_t make_pos_t(gcsa::node_type node);
/// Convert a pos_t to a Position.
Position make_position(const pos_t& pos);
/// Create a Position from a Node ID, an orientation flag, and an offset.
Position make_position(id_t id, bool is_rev, off_t off);
/// Make a Position from a gcsa node
Position make_position(gcsa::node_type node);
/// Find the min distance in the path offsets where the path orientation is the same and different
pair<int64_t, int64_t> min_oriented_distances(const map<string, vector<pair<size_t, bool> > >& path_offsets1,
                                              const map<string, vector<pair<size_t, bool> > >& path_offsets2);

}

#endif
