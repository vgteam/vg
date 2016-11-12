#ifndef VG_POS_H
#define VG_POS_H

#include "vg.pb.h"
#include "types.hpp"
#include "xg.hpp"
#include "lru_cache.h"
#include "utility.hpp"
#include "json2pb.h"
#include <iostream>

namespace vg {

using namespace std;

bool is_empty(const pos_t& pos);
id_t id(const pos_t& pos);
bool is_rev(const pos_t& pos);
off_t offset(const pos_t& pos);
id_t& get_id(pos_t& pos);
bool& get_is_rev(pos_t& pos);
off_t& get_offset(pos_t& pos);
pos_t reverse(const pos_t& pos, size_t node_length);
Position reverse(const Position& pos, size_t node_length);
ostream& operator<<(ostream& out, const pos_t& pos);
pos_t make_pos_t(const Position& pos);
pos_t make_pos_t(id_t id, bool is_rev, off_t off);
Position make_position(const pos_t& pos);
Position make_position(id_t id, bool is_rev, off_t off);

// xg/position traversal helpers with caching
// used by the Sampler and by the Mapper
string xg_cached_node_sequence(id_t id, xg::XG* xgidx, LRUCache<id_t, Node>& node_cache);
size_t xg_cached_node_length(id_t id, xg::XG* xgidx, LRUCache<id_t, Node>& node_cache);
char xg_cached_pos_char(pos_t pos, xg::XG* xgidx, LRUCache<id_t, Node>& node_cache);
map<pos_t, char> xg_cached_next_pos_chars(pos_t pos, xg::XG* xgidx, LRUCache<id_t, Node>& node_cache);
int xg_cached_distance(pos_t pos1, pos_t pos2, xg::XG* xgidx, LRUCache<id_t, Node>& node_cache, int maximum);
set<pos_t> xg_cached_positions_bp_from(pos_t pos, int distance, bool rev, xg::XG* xgidx, LRUCache<id_t, Node>& node_cache);

}

#endif
