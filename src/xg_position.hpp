#ifndef VG_XG_POS_HPP_INCLUDED
#define VG_XG_POS_HPP_INCLUDED

#include <vg/vg.pb.h>
#include "types.hpp"
#include "xg.hpp"
#include "lru_cache.h"
#include "utility.hpp"
#include "json2pb.h"
#include "path.hpp"
#include <gcsa/gcsa.h>
#include <iostream>
#include "algorithms/nearest_offsets_in_paths.hpp"

/** \file 
 * Functions for working with cached Positions and `pos_t`s.
 */

namespace vg {

using namespace std;

// xg/position traversal helpers with caching
// used by the Sampler and by the Mapper
string xg_node_sequence(id_t id, const PathPositionHandleGraph* xgidx);
/// Get the length of a Node from an PathPositionHandleGraph index, with cacheing of deserialized nodes.
size_t xg_node_length(id_t id, const PathPositionHandleGraph* xgidx);
/// Get the node start position in the sequence vector
int64_t xg_node_start(id_t id, const PathPositionHandleGraph* xgidx);
/// Get the character at a position in an PathPositionHandleGraph index, with cacheing of deserialized nodes.
char xg_pos_char(pos_t pos, const PathPositionHandleGraph* xgidx);
/// Get the characters at positions after the given position from an PathPositionHandleGraph index, with cacheing of deserialized nodes.
map<pos_t, char> xg_next_pos_chars(pos_t pos, const PathPositionHandleGraph* xgidx);
set<pos_t> xg_next_pos(pos_t pos, bool whole_node, const PathPositionHandleGraph* xgidx);
int64_t xg_distance(pos_t pos1, pos_t pos2, int64_t maximum, const PathPositionHandleGraph* xgidx);
set<pos_t> xg_positions_bp_from(pos_t pos, int64_t distance, bool rev, const PathPositionHandleGraph* xgidx);
//void xg_graph_context(VG& graph, const pos_t& pos, int length, PathPositionHandleGraph* xgidx);
Node xg_node(id_t id, const PathPositionHandleGraph* xgidx);
vector<Edge> xg_edges_on_start(id_t id, const PathPositionHandleGraph* xgidx);
vector<Edge> xg_edges_on_end(id_t id, const PathPositionHandleGraph* xgidx);

/// Get a map from path name to a list of positions on that path touched by or
/// near to the given Alignment. If nearby is set, search off the part of the
/// graph actually covered by the alignment, to try and find nearby positions
/// along nearby paths. Otherwise, we only search off the actually touched part
/// of the graph if we can't find any paths on it. If just_min is set, produce
/// only one position per path, which will be the lowest graph position
/// *touched* by the Alignment (could be a Mapping start *or* end). Otherwise,
/// produces one position per occurrence on the path of the position of each
/// Mapping in the Alignment, which will be the position at which the Mapping
/// *starts*.
///
/// During the search for nearby positions, walk up to search_limit bases, or
/// the length of the Alignment's sequence if search_limit is 0.
unordered_map<path_handle_t, vector<pair<size_t, bool> > > xg_alignment_path_offsets(const PathPositionHandleGraph* xgidx, const Alignment& aln,
    bool just_min, bool nearby, size_t search_limit = 0);

/// Annotate the given Alignment in place with the earliest touched positions,
/// as produced by xg_alignment_path_offsets, as refpos values. Always uses min
/// positions. Only resorts to nearby positions if no positions on a path are
/// touched. During the search for nearby positions, walk up to search_limit
/// bases, or the length of the Alignment's sequence if search_limit is 0.
void xg_annotate_with_initial_path_positions(const PathPositionHandleGraph* xgidx, Alignment& aln, size_t search_limit = 0);

void xg_neighborhood(const PathPositionHandleGraph& xgidx,int64_t id, size_t dist, Graph& g, bool use_steps);
void xg_add_paths_to_graph(const PathPositionHandleGraph& xgidx, map<int64_t, Node*>& nodes, Graph& g);
void xg_get_id_range(const PathPositionHandleGraph& xgidx, int64_t id1, int64_t id2, Graph& g);


}

#endif
