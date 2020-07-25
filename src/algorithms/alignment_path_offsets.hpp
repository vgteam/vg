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

/// gives just the path positions of the alignment
unordered_map<path_handle_t, vector<pair<size_t, bool> > >
alignment_path_offsets(const PathPositionHandleGraph& graph,
                       const Alignment& aln,
                       bool just_min,
                       bool nearby,
                       size_t search_limit = 0);

/// Find the position of a multipath alignment on paths. Returns the lowest offset
/// position on a path for each contiguous stretch of the multipath alignment, but
/// multiple positions on the same path may be returned if the multipath alignment
/// is disconnected or fans out toward the sources or sinks.
unordered_map<path_handle_t, vector<pair<size_t, bool> > >
multipath_alignment_path_offsets(const PathPositionHandleGraph& graph,
                                 const multipath_alignment_t& mp_aln);

/// Use the graph to annotate an Alignment with the first
/// position it touches on each reference path. Thread safe.
///
/// search_limit gives the maximum distance to search for a path if the
/// alignment does not actually touch any paths. If 0, the alignment's
/// sequence length is used.
void annotate_with_initial_path_positions(const PathPositionHandleGraph& graph, Alignment& aln, size_t search_limit = 0);

/// Use the graph annotate Alignments with the first position
/// they touch on each reference path. Thread safe.
///
/// search_limit gives the maximum distance to search for a path if the
/// alignment does not actually touch any paths. If 0, the alignment's
/// sequence length is used.
void annotate_with_initial_path_positions(const PathPositionHandleGraph& graph, vector<Alignment>& aln, size_t search_limit = 0);


}

}
