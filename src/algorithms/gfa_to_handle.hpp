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

/// Operate on a stream line by line.  This can only work if the GFA is sorted.  If the GFA isn't
/// sorted, dump it to a temp file, and fall back on gfa_to_path_handle_graph()
void gfa_to_path_handle_graph_stream(istream& in,
                                     MutablePathMutableHandleGraph* graph,
                                     int64_t max_rgfa_rank = numeric_limits<int64_t>::max());


/// gfakluge can't parse line by line, which we need for streaming
/// ideally, it needs to be entirely replaced.  here's a bare minimum for parsing lines
/// in the meantime.  they return the fields as strings, don't support overlaps, and
/// optional tags get read as strings in the vectors. 
tuple<string, string, vector<string>> parse_gfa_s_line(const string& s_line);
tuple<string, bool, string, bool, vector<string>> parse_gfa_l_line(const string& l_line);
/// visit_step takes {path-name, rank (-1 if path empty), step id, step reversed}
/// and returns true if it wants to keep iterating (false means stop)
void parse_gfa_p_line(const string& p_line,
                      function<bool(const string&, int64_t, const string&, bool)> visit_step);


}
}

#endif
