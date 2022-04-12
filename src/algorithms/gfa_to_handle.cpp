#include "gfa_to_handle.hpp"
#include "../path.hpp"

#include <bdsg/odgi.hpp>

namespace vg {
namespace algorithms {

void GFAIDMapInfo::invert_translation() {
    if (!numeric_mode) {
        // Make the mapping
        id_to_name.reset(new unordered_map<nid_t, const std::string*>());
        // And then populate it
        for (auto mapping = name_to_id->begin(); mapping != name_to_id->end(); ++mapping) {
            id_to_name->emplace(mapping->second, &mapping->first); 
        }
    }
}

std::vector<oriented_node_range_t> GFAIDMapInfo::translate_back(const oriented_node_range_t& range) const {
    // Nodes haven't been split.
    return {range};
}

std::string GFAIDMapInfo::get_back_graph_node_name(const nid_t& back_node_id) const {
    if (numeric_mode) {
        // Just use string version of number
        return std::to_string(back_node_id);
    }
    // We must have filled in the relevant mapping otherwise.
    assert(id_to_name);
    // Go look up and dereference the name string.
    return *id_to_name->at(back_node_id);
}

static void write_gfa_translation(const GFAIDMapInfo& id_map_info, const string& translation_filename) {
    // don't write anything unless we have both an output file and at least one non-trivial mapping
    if (!translation_filename.empty() && !id_map_info.numeric_mode) {
        ofstream trans_file(translation_filename);
        if (!trans_file) {
            throw runtime_error("error:[gfa_to_handle_graph] Unable to open output translation file: " + translation_filename);
        }
        for (const auto& mapping : *id_map_info.name_to_id) {
            trans_file << "T\t" << mapping.first << "\t" << mapping.second << "\n";
        }
    }
}

/// Add listeners which let a GFA parser fill in a handle graph with nodes and edges.
static void add_graph_listeners(GFAParser& parser, MutableHandleGraph* graph) {
    parser.node_listeners.push_back([&parser, graph](nid_t id, const GFAParser::chars_t& sequence, const GFAParser::tag_list_t& tags) {
        graph->create_handle(GFAParser::extract(sequence), id);
    });
    parser.edge_listeners.push_back([&parser, graph](nid_t from, bool from_is_reverse, nid_t to, bool to_is_reverse, const GFAParser::chars_t& overlap, const GFAParser::tag_list_t& tags) {
        static const string not_blunt = ("error:[gfa_to_handle_graph] Can only load blunt-ended GFAs. "
            "Try \"bluntifying\" your graph with a tool like <https://github.com/vgteam/GetBlunted>, or "
            "transitively merge overlaps with a pipeline of <https://github.com/ekg/gimbricate> and "
            "<https://github.com/ekg/seqwish>.");
        if (GFAParser::length(overlap) > 0) {
            string overlap_text = GFAParser::extract(overlap);
            if (overlap_text != "0M" && overlap_text != "*") {
                // This isn't an allowed overlap value.
                throw GFAFormatError(not_blunt + " Found edge with a non-null alignment '" + overlap_text + "'.");
            }
        }
        
        graph->create_edge(graph->get_handle(from, from_is_reverse),
                           graph->get_handle(to, to_is_reverse));
        
    });
}

/// Add listeners which let a GFA parser fill in a path handle graph with paths.
static void add_path_listeners(GFAParser& parser, MutablePathMutableHandleGraph* graph) {
    parser.path_listeners.push_back([&parser, graph](const string& name, const GFAParser::chars_t& visits, const GFAParser::chars_t& overlaps, const GFAParser::tag_list_t& tags) {
        if (graph->has_path(name)) {
            // Prohibit duplicates
            throw GFAFormatError("Duplicate path name: " + name);
        }
        
        auto path_handle = graph->create_path_handle(name);
        
        // Overlaps are pre-checked in scan_p
        // TODO: Do it in a better place.
        
        GFAParser::scan_p_visits(visits, [&](int64_t step_rank,
                                             const GFAParser::chars_t& step_name,
                                             bool step_is_reverse) {
            if (step_rank >= 0) {
                // Not an empty path sentinel.
                // Find the node ID to visit.
                nid_t n = GFAParser::find_existing_sequence_id(GFAParser::extract(step_name), parser.id_map());
                // And add the step.
                graph->append_step(path_handle, graph->get_handle(n, step_is_reverse));
            }
            // Don't stop.
            return true;
        });
    });
    
    // For rGFA we need to have some state. Use a smart pointer to smuggle it
    // into the closure.
    // For each path name, we remember handle and expected next starting
    // position. If we aren't at the expected next starting position, there's a
    // gap and we need to make a new path.
    // TODO: This duplicates some work with the parser, which also caches rGFA path expected offsets.
    // TODO: Come up with a better listener interface that announces breaks and lets you keep the path handy?
    using rgfa_cache_t = unordered_map<string, pair<path_handle_t, int64_t>>;
    std::shared_ptr<rgfa_cache_t> rgfa_cache = std::make_shared<rgfa_cache_t>();
    
    parser.rgfa_listeners.push_back([&parser, graph, rgfa_cache](nid_t id, int64_t offset, size_t length, const string& path_name, int64_t path_rank) {
        auto found = rgfa_cache->find(path_name);
        if (found != rgfa_cache->end() && found->second.second != offset) {
            // This path already exists, but there's a gap. We need to drop it
            // from the cache and make a new one with the right subpath info.
            rgfa_cache->erase(found);
            found = rgfa_cache->end();
        }
        if (found == rgfa_cache->end()) {
            // Need to make a new path, possibly with subrange start info.
            
            std::pair<int64_t, int64_t> subrange;
            if (offset == 0) {
                // Don't send a subrange
                subrange = PathMetadata::NO_SUBRANGE;
            } else {
                // Start later than 0
                subrange = std::pair<int64_t, int64_t>(offset, PathMetadata::NO_END_POSITION);
            }
            
            // TODO: See if we can split up the path name into a sample/haplotype/etc. to give it a ref sense.
            path_handle_t path = graph->create_path(PathMetadata::SENSE_GENERIC,
                                                    PathMetadata::NO_SAMPLE_NAME,
                                                    path_name, 
                                                    PathMetadata::NO_HAPLOTYPE,
                                                    PathMetadata::NO_PHASE_BLOCK,
                                                    subrange);
            // Then cache it
            found = rgfa_cache->emplace_hint(found, path_name, std::make_pair(path, offset));
        }
        
        // Add the step to the path
        auto& path = found->second.first;
        // rGFA paths always visit sequences forward.
        handle_t step = graph->get_handle(id, false);
        graph->append_step(path, step); 
        
        // Increment the expected next offset
        auto& next_offset = found->second.second;
        next_offset += length;
    });
}

void gfa_to_handle_graph(const string& filename, MutableHandleGraph* graph,
                         GFAIDMapInfo* translation) {
                         
    get_input_file(filename, [&](istream& in) {
       gfa_to_handle_graph(in, graph, translation);
    });
}

void gfa_to_handle_graph(const string& filename, MutableHandleGraph* graph,
                         const string& translation_filename) {

    
    GFAIDMapInfo id_map_info;
    gfa_to_handle_graph(filename, graph, &id_map_info);
    write_gfa_translation(id_map_info, translation_filename);
}

void gfa_to_handle_graph(istream& in, MutableHandleGraph* graph,
                         GFAIDMapInfo* translation) {
                         
    GFAParser parser;
    if (translation) {
        // Use the given external translation so the caller can keep it around.
        parser.external_id_map = translation;
    }
    add_graph_listeners(parser, graph);
    
    parser.parse(in);
}


void gfa_to_path_handle_graph(const string& filename, MutablePathMutableHandleGraph* graph,
                              GFAIDMapInfo* translation, int64_t max_rgfa_rank) {
    
    get_input_file(filename, [&](istream& in) {
        gfa_to_path_handle_graph(in, graph, translation, max_rgfa_rank);
    });
}

void gfa_to_path_handle_graph(const string& filename, MutablePathMutableHandleGraph* graph,
                              int64_t max_rgfa_rank, const string& translation_filename) {

    GFAIDMapInfo id_map_info;
    gfa_to_path_handle_graph(filename, graph, &id_map_info, max_rgfa_rank);
    write_gfa_translation(id_map_info, translation_filename);

}

void gfa_to_path_handle_graph(istream& in,
                              MutablePathMutableHandleGraph* graph,
                              GFAIDMapInfo* translation,
                              int64_t max_rgfa_rank) {
    
    // TODO: Deduplicate this setup code with gfa_to_handle_graph more.
    GFAParser parser;
    if (translation) {
        // Use the given external translation so the caller can keep it around.
        parser.external_id_map = translation;
    }
    add_graph_listeners(parser, graph);
    
    // Set up for path input
    parser.max_rgfa_rank = max_rgfa_rank;
    add_path_listeners(parser, graph);
    
    parser.parse(in);
}

/// Read a range, stopping before any end character in the given null-terminated string.
/// Throws if the range would be empty or none of the characters are encountered.
static GFAParser::chars_t take_range_until(GFAParser::cursor_t& cursor, const GFAParser::cursor_t& end, const char* end_chars, const char* parsing_state = nullptr) {
    auto start = cursor;
    while (cursor != end) {
        for (const char* stop_char = end_chars; *stop_char; ++stop_char) {
            if (*cursor == *stop_char) {
                // We found a stop character
                if (cursor == start) {
                     throw GFAFormatError("Expected nonempty value", cursor, parsing_state);
                }
                return GFAParser::chars_t(start, cursor);
            }
        }
        ++cursor;
    }
    throw GFAFormatError("Expected terminator in " + std::string(end_chars), cursor, parsing_state);
}

/// Read a range, stopping at tab or end of line.
static GFAParser::chars_t take_optional_range(GFAParser::cursor_t& cursor, const GFAParser::cursor_t& end, const char* parsing_state = nullptr) {
    auto start = cursor;
    while (cursor != end && *cursor != '\t') {
        ++cursor;
    }
    return GFAParser::chars_t(start, cursor);
}

/// Read a range, stopping at tab or end of line.
/// Throw if it is empty.
static GFAParser::chars_t take_range(GFAParser::cursor_t& cursor, const GFAParser::cursor_t& end, const char* parsing_state = nullptr) {
    GFAParser::chars_t value = take_optional_range(cursor, end);
    if (GFAParser::empty(value)) {
        throw GFAFormatError("Expected nonempty value", cursor, parsing_state);
    }
    return value;
}

/// Read a string, stopping at tab or end of line.
static string take_optional_string(GFAParser::cursor_t& cursor, const GFAParser::cursor_t& end, const char* parsing_state = nullptr) {
    string value;
    while (cursor != end && *cursor != '\t') {
        value.push_back(*cursor);
        ++cursor;
    }
    return value;
}

/// Read a string, stopping at tab or end of line.
/// Throw if it is empty.
static string take_string(GFAParser::cursor_t& cursor, const GFAParser::cursor_t& end, const char* parsing_state = nullptr) {
    string value = take_optional_string(cursor, end);
    if (value.empty()) {
        throw GFAFormatError("Expected nonempty value", cursor, parsing_state);
    }
    return value;
}

/// Advance past a tab character. If it's not there, return false. If something
/// else is there, return GFAFormatError.
static bool take_optional_tab(GFAParser::cursor_t& cursor, const GFAParser::cursor_t& end, const char* parsing_state = nullptr) {
    if (cursor == end) {
        return false;
    }
    if (*cursor != '\t') {
        throw GFAFormatError("Expected tab", cursor, parsing_state); 
    }
    ++cursor;
    return true;
}

/// Take the given character. Throw an error if it isn't there.
static void take_character(GFAParser::cursor_t& cursor, const GFAParser::cursor_t& end, char value, const char* parsing_state = nullptr) {
    if (cursor == end || *cursor != value) {
        throw GFAFormatError("Expected " + value, cursor, parsing_state); 
    }
    ++cursor;
}

/// Take one character of two options. Return true if it is the first one,
/// false if it is the second, and throw an error if it is absent or something
/// else.
static bool take_flag_character(GFAParser::cursor_t& cursor, const GFAParser::cursor_t& end, char true_value, char false_value, const char* parsing_state = nullptr) {
    if (cursor != end) {
        if (*cursor == true_value) {
            ++cursor;
            return true;
        }
        if (*cursor == false_value) {
            ++cursor;
            return false;
        }
    }
    // Composing the error is tricky because of the bare characters.
    stringstream ss;
    ss << "Expected " << true_value << " or " << false_value;
    throw GFAFormatError(ss.str(), cursor, parsing_state);
}

/// Advance past a tab character that must exist.
static void take_tab(GFAParser::cursor_t& cursor, const GFAParser::cursor_t& end, const char* parsing_state = nullptr) {
    if (!take_optional_tab(cursor, end)) {
        throw GFAFormatError("Expected tab", cursor, parsing_state);
    }
}

GFAParser::tag_list_t GFAParser::parse_tags(const chars_t& tag_range) {
    tag_list_t tags;
    auto cursor = tag_range.first;
    auto& end = tag_range.second;
    
    while (cursor != end) {
        // Scan out a tag of non-tab characters
        string tag = take_string(cursor, end, "parsing tags");
        if (!tag.empty()) {
            // We found a tag. Save it.
            tags.emplace_back(std::move(tag));
        }
        take_optional_tab(cursor, end, "parsing tags");
    }
    
    return tags;
}

tuple<string, GFAParser::chars_t, GFAParser::tag_list_t> GFAParser::parse_s(const string& s_line) {
    auto cursor = s_line.begin();
    auto end = s_line.end();
    
    // Make sure we start with S
    take_character(cursor, end, 'S', "parsing S line start");
    take_tab(cursor, end, "parsing S line");
    
    // Parse out the name
    string name = take_string(cursor, end, "parsing sequence name");
    take_tab(cursor, end, "parsing end of sequence name");
    
    // Parse out the sequence
    chars_t sequence = take_range(cursor, end, "parsing sequence");
    take_optional_tab(cursor, end, "parsing end of sequence");
    
    // Now we're either at the end or at the tab before the tags. Parse the tags.
    auto tags = GFAParser::parse_tags(chars_t(cursor, end));
    
    return make_tuple(std::move(name), std::move(sequence), std::move(tags));
}

tuple<string, bool, string, bool, GFAParser::chars_t, GFAParser::tag_list_t> GFAParser::parse_l(const string& l_line) {
    auto cursor = l_line.begin();
    auto end = l_line.end();
    
    // Make sure we start with L
    take_character(cursor, end, 'L', "parsing L line start");
    take_tab(cursor, end, "parsing L line");
    
    // Parse out the first node name
    string n1 = take_string(cursor, end, "parsing first node name");
    take_tab(cursor, end, "parsing end of first node name");
    
    // Parse the first orientation
    bool n1_reverse = take_flag_character(cursor, end, '-', '+', "parsing first node orientation");
    take_tab(cursor, end, "parsing end of first node orientation");
    
    // Parse out the second node name
    string n2 = take_string(cursor, end, "parsing second node name");
    take_tab(cursor, end, "parsing end of sencod node name");
    
    // Parse the second orientation
    bool n2_reverse = take_flag_character(cursor, end, '-', '+', "parsing second node orientation");
    take_tab(cursor, end, "parsing end of second node orientation");
    
    // Parse out the overlaps
    chars_t overlaps = take_range(cursor, end, "parsing overlaps");
    take_optional_tab(cursor, end, "parsing end of overlaps");
    
    // Now we're either at the end or at the tab before the tags. Parse the tags.
    auto tags = GFAParser::parse_tags(chars_t(cursor, end));
    
    return make_tuple(std::move(n1), n1_reverse, std::move(n2), n2_reverse, std::move(overlaps), std::move(tags));
}

tuple<string, GFAParser::chars_t, GFAParser::chars_t, GFAParser::tag_list_t> GFAParser::parse_p(const string& p_line) {
    auto cursor = p_line.begin();
    auto end = p_line.end();
    
    // Make sure we start with P
    take_character(cursor, end, 'P', "parsing P line start");
    take_tab(cursor, end, "parsing P line");
    
    // Grab the path name
    string path_name = take_string(cursor, end, "parsing path name");
    take_tab(cursor, end, "parsing end of path name");
    
    // Parse out the visits
    chars_t visits = take_range(cursor, end, "parsing path visits");
    take_tab(cursor, end, "parsing end of path visits");
    
    // Parse out the overlaps
    chars_t overlaps = take_range(cursor, end, "parsing overlaps");
    take_optional_tab(cursor, end, "parsing end of overlaps");
    
    // Now we're either at the end or at the tab before the tags. Parse the tags.
    auto tags = GFAParser::parse_tags(chars_t(cursor, end));
    
    return make_tuple(std::move(path_name), std::move(visits), std::move(overlaps), std::move(tags));
}

void GFAParser::scan_p_visits(const chars_t& visit_range,
                              function<bool(int64_t rank, const chars_t& node_name, bool is_reverse)> visit_step) {
    
    auto cursor = visit_range.first;
    auto& end = visit_range.second;
    int64_t rank = 0;
    
    while (cursor != end) {
        // Until we run out of visit list range
        
        // Make a range for the visit node name
        chars_t name_range = take_range_until(cursor, end, "+-", "parsing name of visited node");
        bool is_reverse = take_flag_character(cursor, end, '-', '+', "parsing orientation of visited node");
        
        if (!visit_step(rank, name_range, is_reverse)) {
            // We should stop looping
            return;
        }
        
        if (cursor != visit_range.second) {
            // Go past the comma separator
            take_character(cursor, end, ',', "parsing visit separator");
        }
        // And advance the rank for the next visit
        ++rank;
    }
    
    if (rank == 0) {
        // Nothing was visited. Show an empty path.
        visit_step(-1, chars_t(visit_range.second, visit_range.second), false);
    }
}

void GFAParser::scan_p(const string& p_line,
                       function<bool(const string& path_name, int64_t rank, const chars_t& node_name, bool is_reverse)> visit_step) {
    
    tuple<string, chars_t, chars_t, tag_list_t> p_parse = GFAParser::parse_p(p_line);
    auto& path_name = get<0>(p_parse);
    auto& path_visits = get<1>(p_parse);
    auto& overlaps = get<2>(p_parse);
    
    for(auto it = overlaps.first; it != overlaps.second; ++it) {
        if (*it != '*' && *it != ',' && *it != 'M' && (*it < '0' || *it > '9')) {
            // This overlap isn't just * or a list of * or a list of matches with numbers.
            // We can't handle it
            throw GFAFormatError("Path " + path_name + " has nontrivial overlaps and can't be handled");
        }
    }
    
    scan_p_visits(path_visits, [&](int64_t rank, const chars_t& node_name, bool is_reverse) {
        return visit_step(path_name, rank, node_name, is_reverse);
    });
    
}

bool GFAParser::decode_rgfa_tags(const tag_list_t& tags,
                                 string* out_name,
                                 int64_t* out_offset,
                                 int64_t* out_rank) {
    bool has_sn = false;
    bool has_so = false;
    bool has_sr = false;
    for (auto& tag : tags) {
        // Try and parse each tag.
        // TODO: maybe check for SN:Z:, SO:i:, SR:i: as prefixes?
        size_t sep1 = tag.find(':');
        if (sep1 != string::npos) {
            string tag_name = tag.substr(0, sep1);
            if (tag_name == "SN" || tag_name == "SO" || tag_name == "SR") {
                size_t sep2 = tag.find(':', sep1 + 1);
                if (sep2 != string::npos) {
                    string tag_type = tag.substr(sep1 + 1, sep2 - sep1 - 1);
                    if (tag_name == "SN" && tag_type == "Z") {
                        // We found a string name
                        has_sn = true;
                        if (out_name) {
                            *out_name = tag.substr(sep2 + 1);
                        }
                    } else if (tag_name == "SO" && tag_type == "i") {
                        // We found an integer offset along the path
                        has_so = true;
                        if (out_offset) {
                            *out_offset = stoll(tag.substr(sep2 + 1));
                        }
                    } else if (tag_name == "SR" && tag_type == "i") {
                        // We found an integer rank for the path
                        has_sr = true;
                        if (out_rank) {
                            *out_rank = stoll(tag.substr(sep2 + 1));
                        }
                    }
                }
            }
        }
        if (has_sn && has_so && has_sr) {
            break;
        }
    }
    return has_sn && has_so && has_sr;
}

nid_t GFAParser::assign_new_sequence_id(const string& str, GFAIDMapInfo& id_map_info) {
    
    auto found = id_map_info.name_to_id->find(str);
    if (found != id_map_info.name_to_id->end()) {
        // already in map, so bail out
        return 0;
    }
    
    nid_t node_id = -1;
    if (id_map_info.numeric_mode) {
        if (any_of(str.begin(), str.end(), [](char c) { return !isdigit(c); })) {
            // non-numeric: use max id and add to map
            id_map_info.numeric_mode = false;
        } else {
            node_id = stoll(str);
            if (node_id <= 0) {
                // treat <= 0 as non-numeric
                id_map_info.numeric_mode = false;
            }
        }
    }

    // if numeric, the id was set to stoll above, otherwise we take it from current max
    if (!id_map_info.numeric_mode) {
        node_id = id_map_info.max_id + 1;
    }
    
    id_map_info.max_id = std::max(node_id, id_map_info.max_id);
    id_map_info.name_to_id->emplace_hint(found, str, node_id);

    return node_id;
}

nid_t GFAParser::find_existing_sequence_id(const string& str, GFAIDMapInfo& id_map_info) {
    auto found = id_map_info.name_to_id->find(str);
    if (found != id_map_info.name_to_id->end()) {
        // already in map, just return
        return found->second;
    }
    // Otherwise just fail
    return 0;
}

void GFAParser::parse(istream& in) {
    if (!in) {
        throw std::ios_base::failure("error:[GFAParser] Couldn't open input stream");
    }
    
    // Check if stream is seekable
    in.clear();
    // This tracks the position of the current line
    std::streampos in_pos = in.tellg();
    // And this is what we use when we want to say to read to the end of the file.
    // We will fill it in with a real EOF position if the stream is seekable.
    std::streampos eof_pos = -1;
    bool stream_is_seekable = false;
    if (in_pos >= 0 && in.good()) {
        // Input stream is seekable.
        stream_is_seekable = true;
        // Find EOF
        in.seekg(0, std::ios_base::end);
        eof_pos = in.tellg();
        in.seekg(in_pos);
        if (!in.good()) {
            throw std::runtime_error("Could not get end of GFA file");
        }
    }
    // Reset error flags
    in.clear();
    
#ifdef debug
    std::cerr << "Stream seekable? " << stream_is_seekable << std::endl;
#endif
    
    bool has_rgfa_tags = false;
    string line_buffer; // can be quite big
    
    // We should be able to parse in 2 passes. One to make all the nodes, and
    // one to make all the things that reference nodes we hadn't seen yet.
    // We track pass number 1-based.
    size_t pass_number = 1;
    // And within a pass we remember the line number. Also 1-based.
    size_t line_number = 1;
    
    // We buffer lines we can't actually handle right now.
    
    // If we can seek back to them, we keep this collection of ranges of
    // unprocessed lines. If the second field is the max value, it extends to EOF.
    // Third field is starting line number.
    vector<tuple<std::streampos, std::streampos, size_t>> unprocessed_ranges;
    // And we keep this flag for knowing when we need to close a range.
    bool last_line_handled = true;
    
    // If we can't seek to them, we put them into this temporary file.
    string buffer_name;
    ofstream buffer_out_stream;
    
    // We call this to save a line either in the buffer or the collection of unprocessed ranges,
    // until some time after the given node is observed.
    auto save_line_until_node = [&](const string& missing_node_name) {
        if (pass_number > 1) {
            // We should only hit this on the first pass. If we hit it later we are missing a node.
            throw GFAFormatError("GFA file references missing node " + missing_node_name);
        }
        if (stream_is_seekable) {
            // We should be able to get back here.
            if (unprocessed_ranges.empty() || get<1>(unprocessed_ranges.back()) != eof_pos) {
                // There's not currently a run of unprocessed lines that we are a part of. We need to start a new run.
                unprocessed_ranges.emplace_back(in_pos, eof_pos, line_number);
#ifdef debug
                std::cerr << "Started new unprocessed range at " << in_pos << std::endl;
#endif
                // Run will be closed when a line is handled.
            }
        } else {
            if (buffer_name.empty()) {
                // Make sure the buffer is available
                buffer_name = temp_file::create();
                buffer_out_stream.open(buffer_name);
                if (!buffer_out_stream) {
                    throw runtime_error("error:[GFAParser] Could not open fallback gfa temp file: " + buffer_name);
                }
                // Tell the user that we're having to use a buffer; the tests want us to.
                std::cerr << "warning:[GFAParser] Streaming GFA file references node " << missing_node_name << " before it is defined. "
                          << "Buffering lines in " << buffer_name << " until we are ready for them." << std::endl;
            }
            
            // Store the line into it so we can move on to the next line
            buffer_out_stream << line_buffer << "\n";
        }
    };
    
    // We want to warn about unrecognized line types, but each only once.
    set<char> warned_line_types;
    
    // And we need to buffer all the rGFA visits until we have seen all the nodes.
    // This is a visit at a path offset to a node. We also need the length so
    // we can know when it abuts later visits.
    using rgfa_visit_t = tuple<int64_t, nid_t, size_t>;
    // This is a heap of those, in order
    using visit_queue_t = std::priority_queue<rgfa_visit_t, vector<rgfa_visit_t>, std::greater<rgfa_visit_t>>;
    // This holds rGFA paths we have heard of, mapping from name to rank, start
    // position of next visit that is safe to announce, and list of buffered
    // visits in a min-heap
    unordered_map<string, tuple<int64_t, size_t, visit_queue_t>> rgfa_path_cache;
    
    // We call this to handle the current line if it is ready to be handled, or
    // buffer it if it can't. Return false if we are not ready for the line right
    // now and we saved it, and true otherwise.
    auto handle_line_if_ready = [&]() {
        if (!line_buffer.empty()) {
            switch(line_buffer[0]) {
            case 'H':
                // Header lines don't need anything done
                // TODO: Warn if we see version 0.1; it may have non-standard split P lines.
                break;
            case 'S':
                // Sequence lines can always be handled right now
                {
                    tuple<string, GFAParser::chars_t, tag_list_t> s_parse = GFAParser::parse_s(line_buffer);
                    auto& node_name = get<0>(s_parse);
                    auto& sequence_range = get<1>(s_parse);
                    auto& tags = get<2>(s_parse);
                    nid_t assigned_id = GFAParser::assign_new_sequence_id(node_name, this->id_map());
                    if (assigned_id == 0) {
                        // This name has been used already!
                        throw GFAFormatError("Duplicate sequence name: " + node_name);
                    }
                    for (auto& listener : this->node_listeners) {
                        // Tell all the listener functions
                        listener(assigned_id, sequence_range, tags);
                    }
                    if (this->max_rgfa_rank >= 0 && tags.size() >= 3) {
                        // We'll check for the 3 rGFA optional tags.
                        string rgfa_path_name;
                        int64_t rgfa_offset_on_path;
                        int64_t rgfa_path_rank;
                        if (decode_rgfa_tags(tags, &rgfa_path_name, &rgfa_offset_on_path, &rgfa_path_rank) &&
                            rgfa_path_rank <= this->max_rgfa_rank) {
                            
                            // We need to remember this rGFA path visit
                            auto found = rgfa_path_cache.find(rgfa_path_name);
                            if (found == rgfa_path_cache.end()) {
                                // This is a completely new path, so record its rank
                                found = rgfa_path_cache.emplace_hint(found, rgfa_path_name, std::make_tuple(rgfa_path_rank, 0, visit_queue_t()));
                            } else {
                                // This path existed already. Make sure we aren't showing a conflicting rank
                                if (rgfa_path_rank != get<0>(found->second)) {
                                    throw GFAFormatError("rGFA path " + rgfa_path_name + " has conflicting ranks " + std::to_string(rgfa_path_rank) + " and " + std::to_string(get<0>(found->second)));
                                }
                            }
                            auto& visit_queue = get<2>(found->second);
                            auto& next_offset = get<1>(found->second);
                            auto node_length = GFAParser::length(sequence_range);
                            if (next_offset == rgfa_offset_on_path) {
                                // It's safe to dispatch this visit right now since it's the next one expected along the path.
                                for (auto& listener : this->rgfa_listeners) {
                                    // Tell all the listener functions about this visit
                                    listener(assigned_id, rgfa_offset_on_path, node_length, rgfa_path_name, rgfa_path_rank);
                                }
                                // Advance the offset by the sequence length;
                                next_offset += node_length;
                                while (!visit_queue.empty() && next_offset == get<0>(visit_queue.top())) {
                                    // The lowest-offset queued visit can be handled now because it abuts what we just did.
                                    // Grab the visit.
                                    auto& visit = visit_queue.top();
                                    for (auto& listener : this->rgfa_listeners) {
                                        // Tell all the listener functions about this visit
                                        listener(get<1>(visit), get<0>(visit), get<2>(visit), rgfa_path_name, rgfa_path_rank);
                                    }
                                    // Advance the offset by the sequence length;
                                    next_offset += get<2>(visit);
                                    // And pop the visit off
                                    visit_queue.pop();
                                }
                            } else {
                                // Add this visit to the heap so we can handle it when we find the missing visits.
                                visit_queue.emplace(rgfa_offset_on_path, assigned_id, node_length);
                            }
                        }
                    }
                    return true;
                }
                break;
            case 'L':
                // Edges can be handled if their nodes exist already
                {
                    tuple<string, bool, string, bool, chars_t, tag_list_t> l_parse = GFAParser::parse_l(line_buffer);
                    
                    // We only get these IDs if they have been seen already as nodes
                    nid_t n1 = GFAParser::find_existing_sequence_id(get<0>(l_parse), this->id_map());
                    if (!n1) {
                        save_line_until_node(get<0>(l_parse));
                        return false;
                    }
                    nid_t n2 = GFAParser::find_existing_sequence_id(get<2>(l_parse), this->id_map());
                    if (!n2) {
                        save_line_until_node(get<2>(l_parse));
                        return false;
                    }
                    
                    for (auto& listener : this->edge_listeners) {
                        // Tell all the listener functions
                        listener(n1, get<1>(l_parse), n2, get<3>(l_parse), get<4>(l_parse), get<5>(l_parse));
                    }
                }
                break;
            case 'P':
                // Paths can be handled if all their nodes have been seen
                {
                    bool missing = false;
                    string missing_name;
                    
                    // TODO: we don't check for duplicate path lines here.
                    // Listeners might.
                    
                    // Pass 1: Make sure we have all the nodes in the graph
                    GFAParser::scan_p(line_buffer, [&](const string& path_name,
                                                       int64_t step_rank,
                                                       const GFAParser::chars_t& step_id,
                                                       bool step_is_reverse) {
                         if (step_rank >= 0) {
                            string step_string = GFAParser::extract(step_id);
                            nid_t n = GFAParser::find_existing_sequence_id(step_string, this->id_map());
                            if (!n) {
                                missing = true;
                                missing_name = std::move(step_string);
                                return false;
                            }
                         }
                         return true;
                    });
                    if (missing) {
                        save_line_until_node(missing_name);
                        return false;
                    }
                    
                    // Pass 2: Re-parse the pieces of the line.
                    // TODO: Deduplicate with scan
                    tuple<string, chars_t, chars_t, tag_list_t> p_parse = GFAParser::parse_p(line_buffer);
                    for (auto& listener : this->path_listeners) {
                        // Tell all the listener functions
                        listener(get<0>(p_parse), get<1>(p_parse), get<2>(p_parse), get<3>(p_parse));
                    }
                }
                break;
            default:
                if (!warned_line_types.count(line_buffer[0])) {
                    // Warn once about this weird line type.
                    warned_line_types.insert(line_buffer[0]);
                    cerr << "warning:[GFAParser] Ignoring unrecognized " << line_buffer[0] << " line type" << endl;
                }
            }
        }
        return true;
    };
    
    // We have a function to do a pass over a file of lines. It stops at the
    // given max offset, if set.
    auto process_lines_in_stream = [&](istream& in_stream, std::streampos max_offset) {
        if (stream_is_seekable) {
            // Keep our position in the input stream up to date.
            in_pos = in_stream.tellg();
        }
        while ((!stream_is_seekable || in_pos < max_offset) && getline(in_stream, line_buffer)) {
            // For each line in the input file, before the max offset
            if (!line_buffer.empty()) {
                // Handle all lines in the stream that we can handle now.
                if (handle_line_if_ready()) {
                    // If we handled the line, we need to mark the end of any unhandled range that might be running.
                    if (pass_number== 1 && stream_is_seekable && !unprocessed_ranges.empty() &&
                        get<1>(unprocessed_ranges.back()) == eof_pos) {
                        // the unprocessed range ends where this line started.
                        get<1>(unprocessed_ranges.back()) = in_pos;
#ifdef debug
                        std::cerr << "Ended unprocessed range at " << in_pos << std::endl;
#endif
                    }
                }
            }
            if (stream_is_seekable) {
                // Keep our position in the original input stream up to date.
                in_pos = in_stream.tellg();
            }
            line_number++;
        }
#ifdef debug
        std::cerr << "Stop processing run at " << in_pos << "/" << max_offset << std::endl;
#endif
    };
    
    try {
        
        pass_number = 1;
        line_number = 1;
        process_lines_in_stream(in, eof_pos);
        
        if (stream_is_seekable) {
            if (!unprocessed_ranges.empty()) {
                // Handle unprocessed ranges of the file by seeking back to them.
                
                // Make sure to clear out EOF.
                in.clear();
                
                pass_number = 2;
                for (auto& range : unprocessed_ranges) {
                    in.seekg(get<0>(range));
                    if (!in.good()) {
                        throw std::runtime_error("Unable to seek in GFA stream that should be seekable");
                    }
                    line_number = get<2>(range);
                    process_lines_in_stream(in, get<1>(range));
                }
            }
        } else {
            if (!buffer_name.empty()) {
                // We also have lines in the buffer to handle.
                buffer_out_stream.close();
                ifstream buffer_in_stream(buffer_name);
                pass_number = 2;
                // We forget the original line numbers and restart.
                line_number = 1;
                process_lines_in_stream(buffer_in_stream, eof_pos);
                
                // Clean up buffer before returning
                // TODO: Who takes care of this on GFA format error?
                buffer_in_stream.close();
                unlink(buffer_name.c_str());
            }
        }
        
        
        // Run through any rGFA paths that don't start at 0 or have gaps. 
        for (auto& kv : rgfa_path_cache) {
            auto& rgfa_path_name = kv.first;
            auto& rgfa_path_rank = get<0>(kv.second);
            auto& visit_queue = get<2>(kv.second);
            
            while (!visit_queue.empty()) {
                // Grab the visit.
                auto& visit = visit_queue.top();
                for (auto& listener : this->rgfa_listeners) {
                    // Tell all the listener functions about this visit
                    listener(get<1>(visit), get<0>(visit), get<2>(visit), rgfa_path_name, rgfa_path_rank);
                }
                // And pop the visit off
                visit_queue.pop();
            }
        }
    } catch (GFAFormatError& e) {
        // Inform the error of where exactly it is.
        // We have line buffer in scope and positions are relative to the line buffer.
        
        e.pass_number = pass_number;
        if (!stream_is_seekable && pass_number > 1) {
            // We're working on this temp file. Report it so line numbers make sense.
            e.file_name = buffer_name;
        }
        e.line_number = line_number;
        if (e.has_position) {
            // We can find the column we were at within the line.
            e.column_number = 1 + (e.position - line_buffer.begin());
        }
        
        // Re-throw the new and improved error
        throw e;
    }
}

GFAIDMapInfo& GFAParser::id_map() {
    if (external_id_map) {
        return *external_id_map;
    }
    if (!internal_id_map) {
        internal_id_map = make_unique<GFAIDMapInfo>();
    }
    return *internal_id_map;
};


GFAFormatError::GFAFormatError(const string& message) : std::runtime_error(message) {
    // Nothing to do!
}

GFAFormatError::GFAFormatError(const string& message, const GFAParser::cursor_t& position, const char* parsing_state) : std::runtime_error(message + (parsing_state ? (" while " + std::string(parsing_state)) : "")), has_position(true), position(position) {
    // Nothing to do!
}

const char* GFAFormatError::what() const noexcept {
    if (message_buffer.empty()) {
        // We need to generate the message
        stringstream ss;
        ss << "GFA format error: ";
        
        if (pass_number != 0) {
            // We do the pass first because we might need to report a buffer
            // line number instead of an original line number.
            ss << "On pass " << pass_number << ": ";
        }
        if (!file_name.empty()) {
            ss << "In file " << file_name << ": ";
        }
        if (line_number != 0) {
            ss << "On line " << line_number << ": ";
        }
        if (column_number != 0) {
            ss << "At column " << column_number << ": ";
        }
        
        // Add on the message from the base class
        ss << std::runtime_error::what();
        
        // Save the composed message
        message_buffer = ss.str();
    }
    return message_buffer.c_str();
}


}
}
