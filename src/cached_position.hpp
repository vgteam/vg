#ifndef VG_CACHED_POS_HPP_INCLUDED
#define VG_CACHED_POS_HPP_INCLUDED

#include <vg/vg.pb.h>
#include "types.hpp"
#include "xg.hpp"
#include "lru_cache.h"
#include "utility.hpp"
#include "json2pb.h"
#include <gcsa/gcsa.h>
#include <iostream>

/** \file 
 * Functions for working with cached Positions and `pos_t`s.
 */

namespace vg {

using namespace std;

// xg/position traversal helpers with caching
// used by the Sampler and by the Mapper
string xg_cached_node_sequence(id_t id, XG* xgidx, LRUCache<id_t, Node>& node_cache);
/// Get the length of a Node from an XG index, with cacheing of deserialized nodes.
size_t xg_cached_node_length(id_t id, XG* xgidx, LRUCache<id_t, Node>& node_cache);
/// Get the node start position in the sequence vector
int64_t xg_cached_node_start(id_t id, XG* xgidx, LRUCache<id_t, int64_t>& node_start_cache);
/// Get the character at a position in an XG index, with cacheing of deserialized nodes.
char xg_cached_pos_char(pos_t pos, XG* xgidx, LRUCache<id_t, Node>& node_cache);
/// Get the characters at positions after the given position from an XG index, with cacheing of deserialized nodes.
map<pos_t, char> xg_cached_next_pos_chars(pos_t pos, XG* xgidx, LRUCache<id_t, Node>& node_cache, LRUCache<id_t, vector<Edge> >& edge_cache);
set<pos_t> xg_cached_next_pos(pos_t pos, bool whole_node, XG* xgidx, LRUCache<id_t, Node>& node_cache, LRUCache<id_t, vector<Edge> >& edge_cache);
int64_t xg_cached_distance(pos_t pos1, pos_t pos2, int64_t maximum, XG* xgidx, LRUCache<id_t, Node>& node_cache, LRUCache<id_t, vector<Edge> >& edge_cache);
set<pos_t> xg_cached_positions_bp_from(pos_t pos, int64_t distance, bool rev, XG* xgidx, LRUCache<id_t, Node>& node_cache, LRUCache<id_t, vector<Edge> >& edge_cache);
//void xg_cached_graph_context(VG& graph, const pos_t& pos, int length, XG* xgidx, LRUCache<id_t, Node>& node_cache, LRUCache<id_t, vector<Edge> >& edge_cache);
Node xg_cached_node(id_t id, XG* xgidx, LRUCache<id_t, Node>& node_cache);
vector<Edge> xg_cached_edges_of(id_t id, XG* xgidx, LRUCache<id_t, vector<Edge> >& edge_cache);
vector<Edge> xg_cached_edges_on_start(id_t id, XG* xgidx, LRUCache<id_t, vector<Edge> >& edge_cache);
vector<Edge> xg_cached_edges_on_end(id_t id, XG* xgidx, LRUCache<id_t, vector<Edge> >& edge_cache);

}

#endif
