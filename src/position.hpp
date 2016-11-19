#ifndef VG_POS_H
#define VG_POS_H

#include "vg.pb.h"
#include "types.hpp"
#include "xg.hpp"
#include "lru_cache.h"
#include "utility.hpp"
#include "json2pb.h"
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
/// Convert a pos_t to a Position.
Position make_position(const pos_t& pos);
/// Create a Position from a Node ID, an orientation flag, and an offset.
Position make_position(id_t id, bool is_rev, off_t off);

// xg/position traversal helpers with caching
// used by the Sampler and by the Mapper
/// Get the length of a Node from an xg::XG index, with cacheing of deserialized nodes.
size_t xg_cached_node_length(id_t id, xg::XG* xgidx, LRUCache<id_t, Node>& node_cache);
/// Get the character at a position in an xg::XG index, with cacheing of deserialized nodes.
char xg_cached_pos_char(pos_t pos, xg::XG* xgidx, LRUCache<id_t, Node>& node_cache);
/// Get the characters at positions after the given position from an xg::XG index, with cacheing of deserialized nodes.
map<pos_t, char> xg_cached_next_pos_chars(pos_t pos, xg::XG* xgidx, LRUCache<id_t, Node>& node_cache);

}

#endif
