#ifndef VG_ALGORITHMS_GFA_TO_HANDLE_HPP_INCLUDED
#define VG_ALGORITHMS_GFA_TO_HANDLE_HPP_INCLUDED

/**
 * \file gfa_to_handle.hpp
 *
 * Defines algorithms for copying data from GFA files into handle graphs
 */

#include <iostream>
#include <cctype>
#include <vector>

#include "../handle.hpp"

namespace vg {
namespace algorithms {
using namespace std;

/// This exception will be thrown if the GFA data is not acceptable.
struct GFAFormatError : std::runtime_error {
    // Keep the constructor from a message
    using std::runtime_error::runtime_error;
};

/**
 * Stores ID information for a graph imported from a GFA.
 * Either all IDs are numerically equivalent to their GFA string IDs, or they
 * are stored in the name_to_id map.
 */
struct GFAIDMapInfo : public NamedNodeBackTranslation {
    /// If true, GFA string IDs are just graph numerical IDs.
    bool numeric_mode = true;
    /// This holds the max node ID yet used.
    nid_t max_id = 0;
    /// This maps from GFA string ID to graph numerical ID.
    /// This is behind a unique_ptr so it can be safely pointed into.
    unique_ptr<unordered_map<string, nid_t>> name_to_id = std::make_unique<unordered_map<string, nid_t>>();
    
    /// This inverts the name to ID map, and is populated when
    /// invert_translation is called, so it can be accessed thread-safely.
    unique_ptr<unordered_map<nid_t, const std::string*>> id_to_name;
    
    /**
     * Prepare the backing data structures for get_back_graph_node_name(). Call after name_to_id is complete.
     */
    void invert_translation();
    
    /**
     * Back-translation of node ranges. Is a no-op for imported GFA graphs that
     * haven't been modified, since the GFA graph is itself the backing graph.
     */
    std::vector<oriented_node_range_t> translate_back(const oriented_node_range_t& range) const;
    
    /**
     * Get the GFA sequence name of a node, given its ID.
     * Assumes this will never be called until after name_to_id is fully populated.
     */
    std::string get_back_graph_node_name(const nid_t& back_node_id) const;

};

/// Read a GFA file for a blunt-ended graph into a HandleGraph. Give "-" as a filename for stdin.
///
/// Throws GFAFormatError if the GFA file is not acceptable, and
/// std::ios_base::failure if an IO operation fails. Throws invalid_argument if
/// otherwise misused.
void gfa_to_handle_graph(const string& filename,
                         MutableHandleGraph* graph,
                         GFAIDMapInfo* translation = nullptr);

/// Overload which serializes its translation to a file internally.
void gfa_to_handle_graph(const string& filename,
                         MutableHandleGraph* graph,
                         const string& translation_filename);

/// Same as gfa_to_handle_graph but also adds path elements from the GFA to the graph
void gfa_to_path_handle_graph(const string& filename,
                              MutablePathMutableHandleGraph* graph,
                              GFAIDMapInfo* translation = nullptr,
                              int64_t max_rgfa_rank = numeric_limits<int64_t>::max());

/// Overload which serializes its translation to a file internally.
void gfa_to_path_handle_graph(const string& filename,
                              MutablePathMutableHandleGraph* graph,
                              int64_t max_rgfa_rank,
                              const string& translation_filename);
                              
/// Load a GFA from a stream (assumed not to be seekable or reopenable) into a PathHandleGraph.
/// Does not give max ID hints, and so might be very slow when loading into an ODGI graph.
void gfa_to_path_handle_graph(istream& in,
                              MutablePathMutableHandleGraph* graph,
                              GFAIDMapInfo* translation = nullptr,
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
