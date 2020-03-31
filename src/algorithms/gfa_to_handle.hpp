#ifndef VG_ALGORITHMS_GFA_TO_HANDLE_HPP_INCLUDED
#define VG_ALGORITHMS_GFA_TO_HANDLE_HPP_INCLUDED

/**
 * \file gfa_to_handle.hpp
 *
 * Defines algorithms for copying data from GFA files into handle graphs
 */

#include <iostream>
#include <gfakluge.hpp>
#include <cctype>

#include "../handle.hpp"

namespace vg {
namespace algorithms {
using namespace std;

/// Read a GFA file for a blunt-ended graph into a HandleGraph. Give "-" as a filename for stdin.
///
/// Optionally tries read the GFA from disk without creating an in-memory representation (defaults to
/// in-memory algorithm if input stream is not seekable).
///
/// Also optionally provides a hint about the node ID range to the handle graph implementation before
/// constructing it (defaults to no hint if input stream is not seekable).
void gfa_to_handle_graph(const string& filename,
                         MutableHandleGraph* graph,
                         bool try_from_disk = true,
                         bool try_id_increment_hint = false);

/// Same as above but also adds path elements from the GFA to the graph
void gfa_to_path_handle_graph(const string& filename,
                              MutablePathMutableHandleGraph* graph,
                              bool try_from_disk = true,
                              bool try_id_increment_hint = false);

}
}

#endif
