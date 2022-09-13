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
/// Does not give max ID hints, and so might be very slow when loading into an ODGI graph.
void gfa_to_handle_graph(const string& filename,
                         MutableHandleGraph* graph,
                         GFAIDMapInfo* translation = nullptr);

/// Overload which serializes its translation to a file internally.
void gfa_to_handle_graph(const string& filename,
                         MutableHandleGraph* graph,
                         const string& translation_filename);

/// Load a GFA from a stream (assumed not to be seekable or reopenable) into a HandleGraph.
void gfa_to_handle_graph(istream& in,
                         MutableHandleGraph* graph,
                         GFAIDMapInfo* translation = nullptr);

/// Same as gfa_to_handle_graph but also adds path elements from the GFA to the graph.
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
void gfa_to_path_handle_graph(istream& in,
                              MutablePathMutableHandleGraph* graph,
                              GFAIDMapInfo* translation = nullptr,
                              int64_t max_rgfa_rank = numeric_limits<int64_t>::max());

/**
 * Lower-level tools for parsing GFA elements.
 *
 * Parsing functions return the fields as strings, and don't support overlaps.
 * Optional tags get read as strings in the vectors.
 *
 * Allows you to register "listeners" for different kinds of GFA file items, by
 * adding functions to the various *_listeners vectors. These listeners can
 * raise GFAFormatError or its subclasses if they do not like what the GFA is
 * saying. Some types of GFAFormatError can be caught internally and
 * processing of the file will continue with the next line, but *not* with the
 * next listener for that line, so the user is responsible for worring about
 * what happens if some but not all listeners for something end up getting
 * called because one failed.
 */
class GFAParser {
public:
    
    // We are going to split up existing line buffers.
    // So we need a cursor into one.
    using cursor_t = string::const_iterator;
    // And a range in one
    using chars_t = pair<cursor_t, cursor_t>;
    // And a way to get the string value for one
    inline static string extract(const chars_t& range) {
        return string(range.first, range.second);
    }
    // And a way to get the length of one
    inline static size_t length(const chars_t& range) {
        return range.second - range.first;
    }
    // And a way to tell if one is empty
    inline static bool empty(const chars_t& range) {
        return range.second == range.first;
    }
    // And a type for a collection of GFA tags.
    // This could become a range or list of ranges if we wanted to copy less.
    using tag_list_t = vector<string>;
    
    /**
     * Parse tags out from a possibly empty range to a vector of tag strings.
     */
    static tag_list_t parse_tags(const chars_t& tag_range);
    
    /**
     * Parse an H line to tags
     */
    static tuple<tag_list_t> parse_h(const string& h_line);
    
    /**
     * Parse an S line to name, sequence, and tags
     */
    static tuple<string, chars_t, tag_list_t> parse_s(const string& s_line);
    
    /**
     * Parse an L line to name, is_reverse, name, is_reverse, overlap, and tags
     */
    static tuple<string, bool, string, bool, chars_t, tag_list_t> parse_l(const string& l_line);
    
    /**
     * Parse a P line into name, visits, overlaps, and tags.
     */
    static tuple<string, chars_t, chars_t, tag_list_t> parse_p(const string& p_line);
    
    /**
     * Parse a W line into sample, haplotype, sequence, range (start and end), walk, and tags.
     * If some or all of the range is missing, uses NO_SUBRANGE and NO_END_POSITION form PathMetadata.
     * Doesn't include an end position if a start position isn't set.
     */
    static tuple<string, size_t, string, pair<int64_t, int64_t>, chars_t, tag_list_t> parse_w(const string& p_line);
    
    /**
     * Scan visits extracted from a P line.
     * Calls a callback with all the steps.
     * visit_step takes {rank (-1 if path empty), step node name, step reversed}
     * and returns true if it wants to keep iterating (false means stop).
     */
    static void scan_p_visits(const chars_t& visit_range,
                              function<bool(int64_t rank, const chars_t& node_name, bool is_reverse)> visit_step);
                              
    /**
     * Scan visits extracted from a W line.
     * Calls a callback with all the steps.
     * visit_step takes {rank (-1 if path empty), step node name, step reversed}
     * and returns true if it wants to keep iterating (false means stop).
     */
    static void scan_w_visits(const chars_t& visit_range,
                              function<bool(int64_t rank, const chars_t& node_name, bool is_reverse)> visit_step);
                              
    
    /**
     * Scan visits extracted from a P or W line, as specified in line_type.
     * Calls a callback with all the steps.
     * visit_step takes {rank (-1 if path empty), step node name, step reversed}
     * and returns true if it wants to keep iterating (false means stop).
     */
    static void scan_visits(const chars_t& visit_range, char line_type,
                            function<bool(int64_t rank, const chars_t& node_name, bool is_reverse)> visit_step);
   
    /**
     * Decode rGFA tags from the given list of tags from an S line.
     * Stores rGFA parameters at the given locations if set.
     * Returns true if a complete set of tags was found.
     */
    static bool decode_rgfa_tags(const tag_list_t& tags,
                                 string* out_name = nullptr,
                                 int64_t* out_offset = nullptr,
                                 int64_t* out_rank = nullptr);
    
    /**
     * Parse a GFA name into a numeric id.
     *
     * If all ids are numeric, they will be converted directly with stol.
     *
     * If all ids are non-numeric, they will get incrementing ids beginning
     * with 1, in order they are visited.
     *
     * If they are a mix of numeric and non-numeric, the numberic ones will be
     * converted with stol until the first non-numeric one is found, then it
     * will revert to using max-id.
     *
     * Since non-numeric ids are dependent on the order the nodes are scanned,
     * there is the unfortunate side effect that they will be different
     * sepending on whether the GFA is processed in lexicographic order or file
     * order.
     *
     * If the string ID has been seen before, returns 0.
     */
    static nid_t assign_new_sequence_id(const string& str, GFAIDMapInfo& id_map_info);
    
    /**
     * Find the existing sequence ID for the given node name, or 0 if it has not been seen yet.
     */
    static nid_t find_existing_sequence_id(const string& str, GFAIDMapInfo& id_map_info);
    
    // To actually parse GFA, we stick event listeners on here and then we go
    // through the GFA. It is the parser's job to make sure events aren't fired
    // before events they depend on (so a path is delayed until all the nodes
    // in it are parsed).
    
    // We can either use an internal ID map here
    unique_ptr<GFAIDMapInfo> internal_id_map;
    // Or have this pointed at an external one before we start parsing.
    GFAIDMapInfo* external_id_map = nullptr;
    
    /// Get the ID map we should be using for parsing.
    inline GFAIDMapInfo& id_map();
    
    /// These listeners are called for the header line(s), if any.
    vector<std::function<void(const tag_list_t& tags)>> header_listeners;
    /// These listeners will be called with information for all nodes.
    /// Listeners are protected from duplicate node IDs. 
    vector<std::function<void(nid_t id, const chars_t& sequence, const tag_list_t& tags)>> node_listeners;
    /// These listeners will be called with information for all edges, after
    /// the node listeners for the involved nodes.
    /// Listeners are not protected from duplicate edges.
    vector<std::function<void(nid_t from, bool from_is_reverse, nid_t to, bool to_is_reverse, const chars_t& overlap, const tag_list_t& tags)>> edge_listeners;
    /// These listeners will be called with information for all P line paths,
    /// after the listeners for all involved nodes, and for the first header if any.
    /// Listeners are not protected from duplicate path names.
    vector<std::function<void(const string& name, const chars_t& visits, const chars_t& overlaps, const tag_list_t& tags)>> path_listeners;
    /// These listeners will be called with information for all W line paths,
    /// after the listeners for all involved nodes, and for the first header if any.
    /// Listeners are not protected from duplicate path metadata.
    vector<std::function<void(const string& sample_name, int64_t haplotype, const string& contig_name, const pair<int64_t, int64_t>& subrange, const chars_t& visits, const tag_list_t& tags)>> walk_listeners;
    /// These listeners will be called with each visit of an rGFA path to a
    /// node, after the node listeners for the involved node, but in an
    /// unspecified order with respect to listeners for headers. They will be
    /// called in order along each path. The listener is responsible for
    /// detecting any gaps in the offset space and producing multiple subpaths
    /// if necessary.
    /// Listeners are protected from duplicate paths with the same name and
    /// different ranks, but not from overlaps of nodes in path offset space.
    vector<std::function<void(nid_t id, int64_t offset, size_t length, const string& path_name, int64_t path_rank)>> rgfa_listeners;
    
    /// Include paths from rGFA tags at this rank or lower. Set to -1 to ignore rGFA tags.
    int64_t max_rgfa_rank = -1;
    
    /// Set to true to treat duplicate paths as errors. Otherwise, they will be
    /// treated as warnings and the duplicated will be discarded. Some GFA
    /// files, like the first HPRC graph releases, include duplicate paths.
    bool stop_on_duplicate_paths = false;
    
    /**
     * Parse GFA from the given stream.
     */
    void parse(istream& in);
    
};

/// This exception will be thrown if the GFA data is not acceptable.
struct GFAFormatError : public std::runtime_error {
    /// We can make one from a message
    GFAFormatError(const string& message);
    /// We can also make one with a position and a possibly null parsing state
    GFAFormatError(const string& message, const GFAParser::cursor_t& position, const char* parsing_state = nullptr);
    
    // The error may or may not have a position in a buffer attached.
    bool has_position = false;
    GFAParser::cursor_t position;
    
    // The error also can be annotated file location information when it makes
    // it up the stack to where those things are known.
    // These are all 1-based
    size_t pass_number = 0;
    size_t line_number = 0;
    size_t column_number = 0;
    string file_name = "";
    
    // For making what() messages we need our own message buffer.
    mutable string message_buffer;
    
    /// Return a pointer to a string describing this exception.
    /// Not thread safe.
    virtual const char* what() const noexcept;
};

/// This exception will be thrown if the GFA data includes multiple copies of
/// what we take to be the same path. We need to be able to tolerate this
/// situation at least in some cases because it is true of the HPRC first
/// release graphs, which duplicate paths across P lines and rGFA tags. 
struct GFADuplicatePathError : public GFAFormatError {
    GFADuplicatePathError(const std::string& path_name);
    GFADuplicatePathError(const std::string& path_name, const GFAParser::cursor_t& position, const char* parsing_state = nullptr);
};

}
}

#endif
