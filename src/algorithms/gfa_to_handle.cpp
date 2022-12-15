#include "gfa_to_handle.hpp"
#include "../path.hpp"

#include <gbwtgraph/utils.h>

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

    // For rGFA we need to have some state. Use a smart pointer to smuggle it
    // into the closure.
    // For each path name, we remember handle and expected next starting
    // position. If we aren't at the expected next starting position, there's a
    // gap and we need to make a new path.
    // TODO: This duplicates some work with the parser, which also caches rGFA path expected offsets.
    // TODO: Come up with a better listener interface that announces breaks and lets you keep the path handy?
    using rgfa_cache_t = unordered_map<string, pair<path_handle_t, int64_t>>;
    std::shared_ptr<rgfa_cache_t> rgfa_cache = std::make_shared<rgfa_cache_t>();
    
    // We also need some shared state for making reference sample (RS) tags on the header apply to P and W lines later
    std::shared_ptr<unordered_set<string>> reference_samples = std::make_shared<unordered_set<string>>();
    
    parser.header_listeners.push_back([&parser, reference_samples](const GFAParser::tag_list_t& tags) {
        for (const std::string& tag : tags) {
            if (tag.size() >= 5 &&
                std::equal(gbwtgraph::REFERENCE_SAMPLE_LIST_GFA_TAG.begin(), gbwtgraph::REFERENCE_SAMPLE_LIST_GFA_TAG.end(), tag.begin()) &&
                tag[2] == ':' &&
                tag[3] == 'Z' &&
                tag[4] == ':') {
             
                // This is a reference samples tag like GBWTGraph's GFA parser knows how to parse.
                // Parse the tag's value
                *reference_samples = gbwtgraph::parse_reference_samples_tag(tag.substr(5));
            }
        }
    });

    parser.path_listeners.push_back([&parser, graph, reference_samples](const string& name,
                                                                        const GFAParser::chars_t& visits,
                                                                        const GFAParser::chars_t& overlaps,
                                                                        const GFAParser::tag_list_t& tags) {
        // For P lines, we add the path.
        
        // Parse out the path name's metadata
        PathSense sense;
        string sample;
        string locus;
        size_t haplotype;
        size_t phase_block;
        subrange_t subrange;
        PathMetadata::parse_path_name(name,
                                      sense,
                                      sample,
                                      locus,
                                      haplotype,
                                      phase_block,
                                      subrange);
                                      
        if (sense == PathSense::HAPLOTYPE && reference_samples->count(sample)) {
            // This P line is about a sample that looks like a haplotype but
            // actually wants to be a reference.
            sense = PathSense::REFERENCE;
        } else if (sense == PathSense::REFERENCE && haplotype != PathMetadata::NO_HAPLOTYPE && !reference_samples->count(sample)) {
            // Mimic the GBWTGraph behavior of parsing full PanSN names
            // (sample, haplotype number, contig) as haplotypes by default,
            // even though we use PanSN names in vg to indicate reference
            // sense.
            // TODO: This is super ugly, can we just change the way the
            // metadata name format works, or use a dedicated PanSN parser here
            // instead?
            // TODO: Can we use GBWTGraph's regex priority system?
            sense = PathSense::HAPLOTYPE;
            if (phase_block == PathMetadata::NO_PHASE_BLOCK) {
                // Assign a phase block if none is specified, since haplotypes need one.
                phase_block = 0;
            }
        }
        
        // Compose what we think the path ought to be named.
        // TODO: When we get a has_path that takes fully specified metadata, use that instead.
        string implied_path_name = PathMetadata::create_path_name(sense,
                                                                  sample,
                                                                  locus,
                                                                  haplotype,
                                                                  phase_block,
                                                                  subrange);
        if (graph->has_path(implied_path_name)) {
            // This is a duplicate.
            throw GFADuplicatePathError(implied_path_name);
        }
        
        // Create the path.
        auto path_handle = graph->create_path(sense,
                                              sample,
                                              locus,
                                              haplotype,
                                              phase_block,
                                              subrange);
        
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
    
    parser.walk_listeners.push_back([&parser, graph, reference_samples](const string& sample_name,
                                                                        int64_t haplotype,
                                                                        const string& contig_name,
                                                                        const subrange_t& subrange,
                                                                        const GFAParser::chars_t& visits,
                                                                        const GFAParser::tag_list_t& tags) {
        // For W lines, we add the path with a bit more metadata.
        
        // By default this is interpreted as a haplotype
        PathSense sense;
        
        // We need to determine a phase block
        size_t phase_block;
        // And a haplotype.
        size_t assigned_haplotype = (size_t) haplotype;
        
        string assigned_sample_name;
        if (sample_name == "*") {
            // The sample name is elided from the walk.
            // This walk must be a generic path.
            sense = PathSense::GENERIC;
            // We don't send a sample name.
            assigned_sample_name = PathMetadata::NO_SAMPLE_NAME;
            if (assigned_haplotype != 0) {
                // We can't have multiple haplotypes for a generic path
                throw GFAFormatError("Generic path on omitted (*) sample has nonzero haplotype");
            }
            assigned_haplotype = PathMetadata::NO_HAPLOTYPE;
            phase_block = PathMetadata::NO_PHASE_BLOCK;
        } else {
            // This is probably a sample name we can use
            
            if (reference_samples->count(sample_name)) {
                // This sample is supposed to be reference.
                sense = PathSense::REFERENCE;
                phase_block = PathMetadata::NO_PHASE_BLOCK;
            } else {
                // We're a haplotype
                sense = PathSense::HAPLOTYPE;
                // GFA doesn't really encode phase blocks. Always use the 0th one.
                phase_block = 0;
            }
            
            // Keep the sample name
            assigned_sample_name = sample_name;
        }
        
        // Drop the subrange completely if it starts at 0.
        // TODO: Detect if there are going to be multiple walks describing
        // different subranges, and keep the subrange on the first one even if
        // it starts at 0, because then we know it's really a partial walk.
        subrange_t assigned_subrange = (subrange.first == 0) ? PathMetadata::NO_SUBRANGE : subrange;
        
        // Compose what we think the path ought to be named.
        // TODO: When we get a has_path that takes fully specified metadata, use that instead.
        string implied_path_name = PathMetadata::create_path_name(sense,
                                                                  assigned_sample_name,
                                                                  contig_name,
                                                                  assigned_haplotype,
                                                                  phase_block,
                                                                  assigned_subrange);
        if (graph->has_path(implied_path_name)) {
            // This is a duplicate.
            throw GFADuplicatePathError(implied_path_name);
        }
        
        // Create the path.
        auto path_handle = graph->create_path(sense,
                                              assigned_sample_name,
                                              contig_name,
                                              assigned_haplotype,
                                              phase_block,
                                              assigned_subrange);
        
        GFAParser::scan_w_visits(visits, [&](int64_t step_rank,
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
    
    
    
    parser.rgfa_listeners.push_back([&parser, graph, rgfa_cache](nid_t id,
                                                                 int64_t offset,
                                                                 size_t length,
                                                                 const string& path_name,
                                                                 int64_t path_rank) {
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
            path_handle_t path = graph->create_path(PathSense::GENERIC,
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

/// Read a range, stopping before any end character in the given null-terminated string,
/// or at the end of the input.
/// Throws if the range would be empty.
static GFAParser::chars_t take_range_until_optional(GFAParser::cursor_t& cursor, const GFAParser::cursor_t& end, const char* end_chars, const char* parsing_state = nullptr) {
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
    return GFAParser::chars_t(start, cursor);
}

/// Read a range, stopping before any end character in the given null-terminated string.
/// Throws if the range would be empty or none of the characters are encountered.
static GFAParser::chars_t take_range_until(GFAParser::cursor_t& cursor, const GFAParser::cursor_t& end, const char* end_chars, const char* parsing_state = nullptr) {
    GFAParser::chars_t range = take_range_until_optional(cursor, end, end_chars, parsing_state);
    if (cursor == end) {
        // We didn't find any of the terminators to stop before
        throw GFAFormatError("Expected terminator in " + std::string(end_chars), cursor, parsing_state);
    }
    return range;
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

/// Read a non-negative integer, stopping at tab or end of line.
/// Throw if it is empty. If it is '*', return the given default value.
static int64_t take_number(GFAParser::cursor_t& cursor, const GFAParser::cursor_t& end, int64_t default_value, const char* parsing_state = nullptr) {
    int64_t value = 0;
    if (cursor == end || !((*cursor >= '0' && *cursor <= '9') || *cursor == '*')) {
        // Number is empty and not properly elided
        throw GFAFormatError("Expected natural number", cursor, parsing_state);
    }
    if (*cursor == '*') {
        // Take the * and use the default value
        ++cursor;
        value = default_value;
    } else {
        while (cursor != end && *cursor >= '0' && *cursor <= '9') {
            // Read the base 10 number digit by digit
            value *= 10;
            value += (*cursor - '0');
            ++cursor;
        }
    }
    if (cursor != end && *cursor != '\t' && *cursor != '\n') {
        throw GFAFormatError("Unexpected data at end of number", cursor, parsing_state);
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

tuple<GFAParser::tag_list_t> GFAParser::parse_h(const string& h_line) {
    auto cursor = h_line.begin();
    auto end = h_line.end();
    
    // Make sure we start with H
    take_character(cursor, end, 'H', "parsing H line start");
    take_tab(cursor, end, "parsing H line");
    
    // Now we're either at the end or at the tab before the tags. Parse the tags.
    auto tags = GFAParser::parse_tags(chars_t(cursor, end));
    
    return make_tuple(std::move(tags));
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

tuple<string, size_t, string, pair<int64_t, int64_t>, GFAParser::chars_t, GFAParser::tag_list_t> GFAParser::parse_w(const string& w_line) {
    auto cursor = w_line.begin();
    auto end = w_line.end();
    
    // Make sure we start with W
    take_character(cursor, end, 'W', "parsing W line start");
    take_tab(cursor, end, "parsing W line");
    
    // Grab the sample name
    string sample_name = take_string(cursor, end, "parsing sample name");
    take_tab(cursor, end, "parsing end of sample name");
    
    // Grab the haplotype number
    int64_t haplotype_number = take_number(cursor, end, -1, "parsing haplotype number");
    if (haplotype_number == -1) {
        // This field is required
        throw GFAFormatError("Missing haplotype number in W line", cursor);
    }
    take_tab(cursor, end, "parsing end of haplotype number");
    
    // Grab the sequence/contig/locus name
    string sequence_name = take_string(cursor, end, "parsing sequence name");
    take_tab(cursor, end, "parsing end of sequence name");
    
    // Grab the start and end positions
    int64_t range_start = take_number(cursor, end, -1, "parsing subrange start");
    take_tab(cursor, end, "parsing end of subrange start");
    int64_t range_end = take_number(cursor, end, -1, "parsing subrange end");
    take_tab(cursor, end, "parsing end of subrange end");
    
    // Parse out the visits
    chars_t visits = take_range(cursor, end, "parsing walk visits");
    take_optional_tab(cursor, end, "parsing end of walk visits");
    
    // Now we're either at the end or at the tab before the tags. Parse the tags.
    auto tags = GFAParser::parse_tags(chars_t(cursor, end));
    
    // Process the path subrange a bit. Compose it into the sort of subrange
    // PathMetadata uses.
    pair<int64_t, int64_t> range = PathMetadata::NO_SUBRANGE;
    if (range_start != -1) {
        range.first = range_start;
        if (range_end != -1) {
            range.second = range_end;
        }
    }
    
    return make_tuple(std::move(sample_name), std::move(haplotype_number), std::move(sequence_name), std::move(range), std::move(visits), std::move(tags));
}

void GFAParser::scan_p_visits(const chars_t& visit_range,
                              function<bool(int64_t rank, const chars_t& node_name, bool is_reverse)> visit_step) {
    
    return GFAParser::scan_visits(visit_range, 'P', visit_step);
}

void GFAParser::scan_w_visits(const chars_t& visit_range,
                              function<bool(int64_t rank, const chars_t& node_name, bool is_reverse)> visit_step) {
    
    return GFAParser::scan_visits(visit_range, 'W', visit_step);
}

void GFAParser::scan_visits(const chars_t& visit_range, char line_type,
                            function<bool(int64_t rank, const chars_t& node_name, bool is_reverse)> visit_step) {
    
    auto cursor = visit_range.first;
    auto& end = visit_range.second;
    int64_t rank = 0;
    
    while (cursor != end) {
        // Until we run out of visit list range
        
        bool is_reverse;
        chars_t name_range;
        
        // Parse name and orientation as appropriate for line type
        if (line_type == 'P') {
            // Parse like a path line
            name_range = take_range_until(cursor, end, "+-", "parsing name of visited node");
            is_reverse = take_flag_character(cursor, end, '-', '+', "parsing orientation of visited node");
        } else if (line_type == 'W') {
            // Parse like a walk line
            is_reverse = take_flag_character(cursor, end, '<', '>', "parsing orientation of visited node");
            name_range = take_range_until_optional(cursor, end, "><\t\n", "parsing name of visited node");
        } else {
            throw std::runtime_error("Unimplemented line type for scanning visits");
        }
        
        
        if (!visit_step(rank, name_range, is_reverse)) {
            // We should stop looping
            return;
        }
        
        if (line_type == 'P') {
            // P lines might have comma separators
            if (cursor != visit_range.second) {
                // Go past the comma separator
                take_character(cursor, end, ',', "parsing visit separator");
            }
        }
        // And advance the rank for the next visit
        ++rank;
    }
    
    if (rank == 0) {
        // Nothing was visited. Show an empty path.
        visit_step(-1, chars_t(visit_range.second, visit_range.second), false);
    }
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
    size_t pass_number;
    // And within a pass we remember the line number. Also 1-based.
    size_t line_number;
    
    // We don't want to process any paths until we've seen the header, or we
    // know there isn't one, because it affects interpretation of path lines.
    bool awaiting_header;
    
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
    
    // For error handling, we need a way to tell an error where it is
    auto annotate_error = [&](GFAFormatError& e) {
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
    };
    
    // We call this to save a line either in the buffer or the collection of
    // unprocessed ranges, until all lines have been seen once.
    auto save_line_until_next_pass = [&]() {
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
            }
            
            // Store the line into it so we can move on to the next line
            buffer_out_stream << line_buffer << "\n";
        }
    };
    
    // We call this to save a line either in the buffer or the collection of unprocessed ranges,
    // until some time after the given node is observed.
    auto save_line_until_node = [&](const string& missing_node_name) {
        if (pass_number > 1) {
            // We should only hit this on the first pass. If we hit it later we are missing a node.
            throw GFAFormatError("GFA file references missing node " + missing_node_name);
        }
        if (!stream_is_seekable && buffer_name.empty()) {
            // Warn that we are missing this node because it is the first missing node. The tests want us to.
            #pragma omp critical (cerr)
            std::cerr << "warning:[GFAParser] Streaming GFA file references node " << missing_node_name << " before it is defined. "
                      << "GFA lines will be buffered in a temporary file." << std::endl;
        }
        // TODO: We could be more efficient if we could notice as soon as the node arrives and handle the line.
        // Sadly we can't yet, so just save it for the next pass.
        save_line_until_next_pass();
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
        if (line_buffer.empty()) {
            // No line to handle.
            return true;
        }
        try {
            switch(line_buffer[0]) {
            case 'H':
                // Header lines need tags examoned
                {
                    tuple<tag_list_t> h_parse = GFAParser::parse_h(line_buffer);
                    auto& tags = get<0>(h_parse);
                    for (auto& listener : this->header_listeners) {
                        // Tell all the listener functions
                        listener(tags);
                    }
                    awaiting_header = false;
                }
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
                // Paths can be handled if all their nodes have been seen, and we know enough about the header. 
                {
                    if (awaiting_header) {
                        save_line_until_next_pass();
                        return false;
                    }
                    
                    bool missing = false;
                    string missing_name;
                    
                    // TODO: we don't check for duplicate path lines here.
                    // Listeners might.
                    
                    // Parse out the path pieces: name, visits, overlaps, tags
                    tuple<string, chars_t, chars_t, tag_list_t> p_parse = GFAParser::parse_p(line_buffer);
                    auto& path_name = get<0>(p_parse);
                    auto& visits = get<1>(p_parse);
                    auto& overlaps = get<2>(p_parse);
                    auto& tags = get<3>(p_parse);
                    
                    for(auto it = overlaps.first; it != overlaps.second; ++it) {
                        if (*it != '*' && *it != ',' && *it != 'M' && (*it < '0' || *it > '9')) {
                            // This overlap isn't just * or a list of * or a list of matches with numbers.
                            // We can't handle it
                            throw GFAFormatError("Path " + path_name + " has nontrivial overlaps and can't be handled", it);
                        }
                    }
                    
                    // Make sure we have all the nodes in the graph
                    GFAParser::scan_p_visits(visits, [&](int64_t step_rank,
                                                         const GFAParser::chars_t& step_id,
                                                         bool step_is_reverse) {
                         if (step_rank == -1) {
                            // Nothing to do for empty paths
                            return true;
                        }
                         
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
                    
                    for (auto& listener : this->path_listeners) {
                        // Tell all the listener functions
                        listener(path_name, visits, overlaps, tags);
                    }
                }
                break;
            case 'W':
                // Walks can be handled if all their nodes have been seen, and we know enough about the hneader.
                {
                    if (awaiting_header) {
                        save_line_until_next_pass();
                        return false;
                    }
                
                    bool missing = false;
                    string missing_name;
                    
                    // Fins the pieces of the walk line
                    tuple<string, size_t, string, pair<int64_t, int64_t>, chars_t, tag_list_t> w_parse = GFAParser::parse_w(line_buffer);
                    auto& sample_name = get<0>(w_parse);
                    auto& haplotype = get<1>(w_parse);
                    auto& contig_name = get<2>(w_parse);
                    auto& subrange = get<3>(w_parse);
                    auto& visits = get<4>(w_parse);
                    auto& tags = get<5>(w_parse);
                    
                    GFAParser::scan_w_visits(visits, [&](int64_t step_rank,
                                                         const GFAParser::chars_t& step_id,
                                                         bool step_is_reverse) {
                        if (step_rank == -1) {
                            // Nothing to do for empty paths
                            return true;
                        }
                        
                        // For every node the walk visits, make sure we have seen it.
                        string step_string = GFAParser::extract(step_id);
                        nid_t n = GFAParser::find_existing_sequence_id(step_string, this->id_map());
                        if (!n) {
                            missing = true;
                            missing_name = std::move(step_string);
                            return false;
                        }
                        return true;
                    });
                    
                    if (missing) {
                        save_line_until_node(missing_name);
                        return false;
                    }
                    
                    for (auto& listener : this->walk_listeners) {
                        // Tell all the listener functions
                        listener(sample_name, haplotype, contig_name, subrange, visits, tags);
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
        } catch (GFADuplicatePathError& e) {
            // We couldn't do what this line said because we already have this path.
            if (stop_on_duplicate_paths) {
                // That's bad. Stop parsing.
                throw;
            } else {
                // We can tolerate this. Just move on to the next line.
                
                // Make sure the error is annotated with position info.
                // TODO: Invent some cool kind of stack-based context?
                annotate_error(e);
                
                // And report it as a warning.
                #pragma omp critical (cerr)
                std::cerr << "warning:[GFAParser] Skipping GFA " << line_buffer[0]
                    << " line: " << e.what() << std::endl;
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
        awaiting_header = true;
        process_lines_in_stream(in, eof_pos);
        
        if (stream_is_seekable) {
            if (!unprocessed_ranges.empty()) {
                // Handle unprocessed ranges of the file by seeking back to them.
                
                // Make sure to clear out EOF.
                in.clear();
                
                pass_number = 2;
                // There can't be any new headers.
                awaiting_header = false;
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
                // There can't be any new headers.
                awaiting_header = false;
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
        // Tell the error where it happened
        annotate_error(e);
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

GFADuplicatePathError::GFADuplicatePathError(const std::string& path_name) : GFAFormatError("Duplicate path " + path_name + " exists in graph") {
    // Nothing to do!
}

GFADuplicatePathError::GFADuplicatePathError(const std::string& path_name, const GFAParser::cursor_t& position, const char* parsing_state) : GFAFormatError("Duplicate path " + path_name + " exists in graph", position, parsing_state) {
    // Nothing to do!
}


}
}
