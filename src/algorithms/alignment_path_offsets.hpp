#pragma once

#include "../handle.hpp"
#include "../position.hpp"
#include "../path.hpp"
#include "../multipath_alignment.hpp"
#include <vg/vg.pb.h>
#include <functional>
#include "nearest_offsets_in_paths.hpp"

namespace vg {
namespace algorithms {

using namespace std;

/// Gives the path positions of the alignment. If just_min is set, gives the
/// minimum position on each path. Else, gives all Mapping start positions on
/// each path. If nearby is set, will search for a nearby path. Will recurse
/// with nearby set if it is not set on initial call and no positions are
/// found. Respects search_limit in bp in that case. If search_limit is 0, read
/// length is used.
///
/// If path_filter is set, and it returns false for a path, that path is not
/// used to annotate the read.
unordered_map<path_handle_t, vector<pair<size_t, bool> > >
alignment_path_offsets(const PathPositionHandleGraph& graph,
                       const Alignment& aln,
                       bool just_min,
                       bool nearby,
                       size_t search_limit = 0,
                       const std::function<bool(const path_handle_t&)>* path_filter = nullptr);

/// Find the position of a multipath alignment on paths. Returns the lowest offset
/// position on a path for each contiguous stretch of the multipath alignment, but
/// multiple positions on the same path may be returned if the multipath alignment
/// is disconnected or fans out toward the sources or sinks.
///
/// If path_filter is set, and it returns false for a path, that path is not
/// used to annotate the read.
unordered_map<path_handle_t, vector<pair<size_t, bool> > >
multipath_alignment_path_offsets(const PathPositionHandleGraph& graph,
                                 const multipath_alignment_t& mp_aln,
                                 const std::function<bool(const path_handle_t&)>* path_filter = nullptr);

/// Use the graph to annotate an Alignment with the first
/// position it touches on each reference path. Thread safe.
///
/// search_limit gives the maximum distance to search for a path if the
/// alignment does not actually touch any paths. If 0, the alignment's
/// sequence length is used.
///
/// If path_filter is set, and it returns false for a path, that path is not
/// used to annotate the read.
void annotate_with_initial_path_positions(const PathPositionHandleGraph& graph, Alignment& aln, size_t search_limit = 0, const std::function<bool(const path_handle_t&)>* path_filter = nullptr);

/// Use the graph to annotate an Alignment with the first
/// position it touches on each node it visits in each reference path. Thread
/// safe. If no nodes on any path are visited, searches for a nearby path
/// position to use.
///
/// search_limit gives the maximum distance to search for a path if the
/// alignment does not actually touch any paths. If 0, the alignment's
/// sequence length is used.
///
/// If path_filter is set, and it returns false for a path, that path is not
/// used to annotate the read.
void annotate_with_node_path_positions(const PathPositionHandleGraph& graph, Alignment& aln, size_t search_limit = 0, const std::function<bool(const path_handle_t&)>* path_filter = nullptr);

/// Use the graph to annotate an Alignment with positions on each reference
/// path. Thread safe.
///
/// If just_min is set, gives the minimum position on each path. Else, gives
/// all Mapping start positions on each path. If no positions on the path are
/// found, looks for nearby path positions in graph space. Respects
/// search_limit in bp in that case. If search_limit is 0, read length is used.
///
/// If path_filter is set, and it returns false for a path, that path is not
/// used to annotate the read.
void annotate_with_path_positions(const PathPositionHandleGraph& graph, Alignment& aln, bool just_min, size_t search_limit = 0, const std::function<bool(const path_handle_t&)>* path_filter = nullptr);

/// Use the graph annotate Alignments with the first position
/// they touch on each reference path. Thread safe.
///
/// search_limit gives the maximum distance to search for a path if the
/// alignment does not actually touch any paths. If 0, the alignment's
/// sequence length is used.
///
/// If path_filter is set, and it returns false for a path, that path is not
/// used to annotate the read.
void annotate_with_initial_path_positions(const PathPositionHandleGraph& graph, vector<Alignment>& aln, size_t search_limit = 0, const std::function<bool(const path_handle_t&)>* path_filter = nullptr);


}

}
