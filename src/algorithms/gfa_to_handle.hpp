// SPDX-FileCopyrightText: 2014 Erik Garrison
//
// SPDX-License-Identifier: MIT

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

/// This exception will be thrown if the GFA data is not acceptable.
struct GFAFormatError : std::runtime_error {
    // Keep the constructor from a message
    using std::runtime_error::runtime_error;
};

/// Read a GFA file for a blunt-ended graph into a HandleGraph. Give "-" as a filename for stdin.
///
/// Optionally tries read the GFA from disk without creating an in-memory representation (defaults to
/// in-memory algorithm if reading from stdin).
///
/// Also optionally provides a hint about the node ID range to the handle graph implementation before
/// constructing it (defaults to no hint if reading from stdin).
///
/// Throws GFAFormatError if the GFA file is not acceptable, and
/// std::ios_base::failure if an IO operation fails. Throws invalid_argument if
/// otherwise misused.
void gfa_to_handle_graph(const string& filename,
                         MutableHandleGraph* graph,
                         bool try_from_disk = true,
                         bool try_id_increment_hint = false,
                         const string& translation_filename = "");

/// Same as gfa_to_handle_graph but also adds path elements from the GFA to the graph
void gfa_to_path_handle_graph(const string& filename,
                              MutablePathMutableHandleGraph* graph,
                              bool try_from_disk = true,
                              bool try_id_increment_hint = false,
                              int64_t max_rgfa_rank = numeric_limits<int64_t>::max(),
                              const string& translation_filename = "");
                              
/// Same as above but operating on a stream. Assumed to be non-seekable; all conversion happens in memory.
/// Always streaming. Doesn't support ID increment hints.
void gfa_to_path_handle_graph_in_memory(istream& in,
                                        MutablePathMutableHandleGraph* graph,
                                        int64_t max_rgfa_rank = numeric_limits<int64_t>::max());


}
}

#endif
