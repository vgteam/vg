#ifndef VG_POS_HPP_INCLUDED
#define VG_POS_HPP_INCLUDED

#include "vg.pb.h"
#include "types.hpp"
#include "xg.hpp"
#include "lru_cache.h"
#include "utility.hpp"
#include "json2pb.h"
#include "gcsa.h"
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

// xg/position traversal helpers with caching
// used by the Sampler and by the Mapper
string xg_cached_node_sequence(id_t id, xg::XG* xgidx, LRUCache<id_t, Node>& node_cache);
/// Get the length of a Node from an xg::XG index, with cacheing of deserialized nodes.
size_t xg_cached_node_length(id_t id, xg::XG* xgidx, LRUCache<id_t, Node>& node_cache);
/// Get the node start position in the sequence vector
int64_t xg_cached_node_start(id_t id, xg::XG* xgidx, LRUCache<id_t, int64_t>& node_start_cache);
/// Get the character at a position in an xg::XG index, with cacheing of deserialized nodes.
char xg_cached_pos_char(pos_t pos, xg::XG* xgidx, LRUCache<id_t, Node>& node_cache);
/// Get the characters at positions after the given position from an xg::XG index, with cacheing of deserialized nodes.
map<pos_t, char> xg_cached_next_pos_chars(pos_t pos, xg::XG* xgidx, LRUCache<id_t, Node>& node_cache, LRUCache<id_t, vector<Edge> >& edge_cache);
set<pos_t> xg_cached_next_pos(pos_t pos, bool whole_node, xg::XG* xgidx, LRUCache<id_t, Node>& node_cache, LRUCache<id_t, vector<Edge> >& edge_cache);
int xg_cached_distance(pos_t pos1, pos_t pos2, int maximum, xg::XG* xgidx, LRUCache<id_t, Node>& node_cache, LRUCache<id_t, vector<Edge> >& edge_cache);
set<pos_t> xg_cached_positions_bp_from(pos_t pos, int distance, bool rev, xg::XG* xgidx, LRUCache<id_t, Node>& node_cache, LRUCache<id_t, vector<Edge> >& edge_cache);
//void xg_cached_graph_context(VG& graph, const pos_t& pos, int length, xg::XG* xgidx, LRUCache<id_t, Node>& node_cache, LRUCache<id_t, vector<Edge> >& edge_cache);
Node xg_cached_node(id_t id, xg::XG* xgidx, LRUCache<id_t, Node>& node_cache);
vector<Edge> xg_cached_edges_of(id_t id, xg::XG* xgidx, LRUCache<id_t, vector<Edge> >& edge_cache);
vector<Edge> xg_cached_edges_on_start(id_t id, xg::XG* xgidx, LRUCache<id_t, vector<Edge> >& edge_cache);
vector<Edge> xg_cached_edges_on_end(id_t id, xg::XG* xgidx, LRUCache<id_t, vector<Edge> >& edge_cache);

}

//namespace std {
///// hash function for pos_t
//template<>
//struct hash<const vg::pos_t> {
//    inline size_t operator()(const vg::pos_t& pos) const {
//        return hash<pair<pair<vg::id_t, bool>, vg::off_t > >()(make_pair(make_pair(vg::id(pos), vg::is_rev(pos)), vg::offset(pos)));
//    }
//};
//}

#endif
