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

static void validate_gfa_edge(const gfak::edge_elem& e) {
    static const string not_blunt = ("error:[gfa_to_handle_graph] Can only load blunt-ended GFAs. "
        "Try \"bluntifying\" your graph with a tool like <https://github.com/vgteam/GetBlunted>, or "
        "transitively merge overlaps with a pipeline of <https://github.com/ekg/gimbricate> and "
        "<https://github.com/ekg/seqwish>.");
    if (e.source_begin != e.source_end || e.sink_begin != 0 || e.sink_end != 0) {
        throw GFAFormatError(not_blunt + " Found edge with an overlay: " + e.source_name + "[" + to_string(e.source_begin) + ":" + to_string(e.source_end) + "] -> " + e.sink_name + "[" + to_string(e.sink_begin) + ":" + to_string(e.sink_end) + "]");
    }
    if (!(e.alignment == "0M" || e.alignment == "*" || e.alignment.empty())) {
        throw GFAFormatError(not_blunt + " Found edge with a non-null alignment '" + e.alignment + "'.");
    }
    if (e.source_name.empty()) {
        throw GFAFormatError("error:[gfa_to_handle_graph] Found edge record with missing source name");
    }
    if (e.sink_name.empty()) {
        throw GFAFormatError("error:[gfa_to_handle_graph] Found edge record with missing sink name");
    }
}

static string process_raw_gfa_path_name(const string& path_name_raw)  {
    string processed = path_name_raw;
    processed.erase(remove_if(processed.begin(), processed.end(),
                              [](char c) { return isspace(c); }),
                    processed.end());
    return processed;
}

/// return whether a gfa node has all 3 rGFA tags
/// optionally parse them
static bool gfa_sequence_parse_rgfa_tags(const gfak::sequence_elem& s,
                                         string* out_name = nullptr,
                                         int64_t* out_offset = nullptr,
                                         int64_t* out_rank = nullptr) {
    bool has_sn = false;
    bool has_so = false;
    bool has_sr = false; 
    for (size_t i = 0 ; i < s.opt_fields.size() && (!has_sn || !has_so || !has_sr); ++i) {
        if (s.opt_fields[i].key == "SN" && s.opt_fields[i].type == "Z") {
            has_sn = true;
            if (out_name) {
                *out_name = s.opt_fields[i].val;
            }
        } else if (s.opt_fields[i].key == "SO" && s.opt_fields[i].type == "i") {
            has_so = true;
            if (out_offset) {
                *out_offset = stol(s.opt_fields[i].val);
            }
        } else if (s.opt_fields[i].key == "SR" && s.opt_fields[i].type == "i") {
            has_sr = true;
            if (out_rank) {
                *out_rank = stol(s.opt_fields[i].val);
            }
        }
    }
    return has_sn && has_so && has_sr;
}

static bool gfa_to_handle_graph_in_memory(istream& in, MutableHandleGraph* graph,
                                          gfak::GFAKluge& gg, GFAIDMapInfo& id_map_info) {
    if (!in) {
        throw std::ios_base::failure("error:[gfa_to_handle_graph] Couldn't open input stream");
    }
    gg.parse_gfa_file(in);

    // create nodes
    bool has_rgfa_tags = false;
    for (const auto& seq_record : gg.get_name_to_seq()) {
        graph->create_handle(seq_record.second.sequence, GFAParser::parse_sequence_id(seq_record.first, id_map_info));
        has_rgfa_tags = has_rgfa_tags || gfa_sequence_parse_rgfa_tags(seq_record.second);
    }
    
    // create edges
    for (const auto& links_record : gg.get_seq_to_edges()) {
        for (const auto& edge : links_record.second) {
            validate_gfa_edge(edge);
            nid_t a_id = GFAParser::parse_sequence_id(edge.source_name, id_map_info);
            if (!graph->has_node(a_id)) {
                throw GFAFormatError("error:[gfa_to_handle_graph] GFA edge starts at nonexistent GFA node \"" + edge.source_name + "\"");
            }
            nid_t b_id = GFAParser::parse_sequence_id(edge.sink_name, id_map_info);
            if (!graph->has_node(b_id)) {
                throw GFAFormatError("error:[gfa_to_handle_graph] GFA edge ends at nonexistent GFA node \"" + edge.sink_name + "\"");
            }
            
            // note: we're counting on implementations de-duplicating edges
            handle_t a = graph->get_handle(a_id, !edge.source_orientation_forward);
            handle_t b = graph->get_handle(b_id, !edge.sink_orientation_forward);
            graph->create_edge(a, b);
        }
    }
    return has_rgfa_tags;
}

static bool gfa_to_handle_graph_on_disk(const string& filename, MutableHandleGraph* graph,
                                        gfak::GFAKluge& gg, GFAIDMapInfo& id_map_info) {
    
    // adapted from
    // https://github.com/vgteam/odgi/blob/master/src/gfa_to_handle.cpp
    
    if (dynamic_cast<bdsg::ODGI*>(graph)) {
        // This kind of graph needs a hint about IDs to be efficient. 
        
        // find the minimum ID
        nid_t min_id = numeric_limits<nid_t>::max();
        gg.for_each_sequence_line_in_file(filename.c_str(), [&](gfak::sequence_elem s) {
                min_id = std::min(min_id, GFAParser::parse_sequence_id(s.name, id_map_info));
            });
        
        if (min_id != numeric_limits<nid_t>::max()) {
            // we found the min, set it as the increment
            graph->set_id_increment(min_id);
        }
    }
    
    // add in all nodes
    bool has_rgfa_tags = false;
    gg.for_each_sequence_line_in_file(filename.c_str(), [&](gfak::sequence_elem s) {        
            graph->create_handle(s.sequence, GFAParser::parse_sequence_id(s.name, id_map_info));
            has_rgfa_tags = has_rgfa_tags || gfa_sequence_parse_rgfa_tags(s);
    });
    
    // add in all edges
    gg.for_each_edge_line_in_file(filename.c_str(), [&](gfak::edge_elem e) {
        validate_gfa_edge(e);
        nid_t a_id = GFAParser::parse_sequence_id(e.source_name, id_map_info);
        if (!graph->has_node(a_id)) {
            throw GFAFormatError("error:[gfa_to_handle_graph] GFA edge starts at nonexistent GFA node \"" + e.source_name + "\"");
        }
        nid_t b_id = GFAParser::parse_sequence_id(e.sink_name, id_map_info);
        if (!graph->has_node(b_id)) {
            throw GFAFormatError("error:[gfa_to_handle_graph] GFA edge ends at nonexistent GFA node \"" + e.sink_name + "\"");
        }
        
        handle_t a = graph->get_handle(a_id, !e.source_orientation_forward);
        handle_t b = graph->get_handle(b_id, !e.sink_orientation_forward);
        graph->create_edge(a, b);
    });
    return has_rgfa_tags;
}

/// Parse nodes and edges and load them into the given GFAKluge.
/// If the input is a seekable file, filename will be filled in and unseekable will be nullptr.
/// If the input is not a seekable file, filename may be filled in, and unseekable will be set to a stream to read from.
/// Returns true if any "SN" rGFA tags are found in the graph nodes
static bool gfa_to_handle_graph_load_graph(const string& filename, istream* unseekable, MutableHandleGraph* graph,
                                           gfak::GFAKluge& gg, GFAIDMapInfo& id_map_info) {
    
    if (graph->get_node_count() > 0) {
        throw invalid_argument("error:[gfa_to_handle_graph] Must parse GFA into an empty graph");
    }
    bool has_rgfa_tags = false;
    if (!unseekable) {
        // Do the from-disk path
        has_rgfa_tags = gfa_to_handle_graph_on_disk(filename, graph, gg, id_map_info);
    } else {
        // Do the path for streams
        
        if (dynamic_cast<bdsg::ODGI*>(graph)) {
            // This kind of graph needs a hint about IDs to be efficient. 
            // But, the ID increment hint can't be done.
            cerr << "warning:[gfa_to_handle_graph] Skipping node ID increment hint because input stream for GFA does not support seeking. "
                 << "If performance suffers, consider using an alternate graph implementation or reading GFA from hard disk." << endl;
        }
        
        has_rgfa_tags = gfa_to_handle_graph_in_memory(*unseekable, graph, gg, id_map_info);
    }
    return has_rgfa_tags;
}

/// After the given GFAKluge has been populated with nodes and edges, load path information.
/// If the input is a seekable file, filename will be filled in and unseekable will be nullptr.
/// If the input is not a seekable file, filename may be filled in, and unseekable will be set to a stream to read from.
static void gfa_to_handle_graph_add_paths(const string& filename, istream* unseekable, MutablePathHandleGraph* graph,
                                          gfak::GFAKluge& gg, GFAIDMapInfo& id_map_info) {
                                   
                                   
    if (!unseekable) {
        // Input is from a seekable file on disk.
        
        // add in all paths
        gg.for_each_path_element_in_file(filename.c_str(), [&](const string& path_name_raw,
                                                               const string& node_id,
                                                               bool is_rev,
                                                               const string& cigar,
                                                               bool is_empty,
                                                               bool is_circular) {
            // remove white space in path name
            // TODO: why?
            string path_name = process_raw_gfa_path_name(path_name_raw);
            
            // get either the existing path handle or make a new one
            path_handle_t path;
            if (!graph->has_path(path_name)) {
                path = graph->create_path_handle(path_name, is_circular);
            } else {
                path = graph->get_path_handle(path_name);
            }
            
            // add the step
            nid_t target_node_id = GFAParser::parse_sequence_id(node_id, id_map_info);
            if (!graph->has_node(target_node_id)) {
                // We need to make sure the GFA isn't lying about the nodes
                // that exist or we will fail with weird errors in get_handle
                // or even later.
                throw GFAFormatError("error:[gfa_to_handle_graph] GFA path " + path_name_raw + " visits nonexistent GFA node \"" + node_id + "\"");
            }
            handle_t step = graph->get_handle(target_node_id, is_rev);
            graph->append_step(path, step);
        });
    } else {
        
        // gg will have parsed the GFA file in the non-path part of the algorithm
        // No reading to do.
        
        // create paths
        for (const auto& path_record : gg.get_name_to_path()) {
            
            // process this to match the disk backed implementation
            // TODO: why?
            string path_name = process_raw_gfa_path_name(path_record.first);
            path_handle_t path = graph->create_path_handle(path_name);
            
            for (size_t i = 0; i < path_record.second.segment_names.size(); ++i) {
                const string& node_id = path_record.second.segment_names.at(i);
                nid_t target_node_id = GFAParser::parse_sequence_id(node_id, id_map_info);
                if (!graph->has_node(target_node_id)) {
                    // The GFA wants to go somewhere that doesn't exist
                    throw GFAFormatError("error:[gfa_to_handle_graph] GFA path " + path_record.first + " visits nonexistent GFA node \"" + node_id + "\"");
                }
                handle_t step = graph->get_handle(target_node_id, !path_record.second.orientations.at(i));
                graph->append_step(path, step);
            }
        }
    }
    
    
}

/// add paths from the optional rgfa tags on sequence nodes (SN: name SO: offset SR: rank)
/// max_rank selects which ranks to consider. Usually, only rank-0 paths are full paths while rank > 0 are subpaths
static void gfa_to_handle_graph_add_rgfa_paths(const string filename, istream* unseekable, vector<gfak::sequence_elem>* rgfa_seq_elems,
                                               MutablePathHandleGraph* graph,
                                               gfak::GFAKluge& gg, GFAIDMapInfo& id_map_info,
                                               int64_t max_rank) {

    // build up paths in memory using a plain old stl structure
    // maps path-name to <rank, vector<node_id, offset>>
    unordered_map<string, pair<int64_t, vector<pair<nid_t, int64_t>>>> path_map;

    string rgfa_name;
    int64_t rgfa_offset;
    int64_t rgfa_rank;

    function<void(const gfak::sequence_elem&)> update_rgfa_path = [&](const gfak::sequence_elem& s) {
        if (gfa_sequence_parse_rgfa_tags(s, &rgfa_name, &rgfa_offset, &rgfa_rank) && rgfa_rank <= max_rank) {
            pair<int64_t, vector<pair<nid_t, int64_t>>>& val = path_map[rgfa_name];
            if (!val.second.empty()) {
                if (val.first != rgfa_rank) {
                    cerr << "warning:[gfa_to_handle_graph] Ignoring rGFA tags for sequence " << s.name
                         << " because they identify it as being on path " << rgfa_name << " with rank " << rgfa_rank
                         << " but a path with that name has already been found with a different rank (" << val.first << ")" << endl;
                    return;
                }
            } else {
                val.first = rgfa_rank;
            }
            nid_t seq_id = GFAParser::parse_sequence_id(s.name, id_map_info);
            // We can assume the nodes exist because we're looking at sequence lines already here.
            val.second.push_back(make_pair(seq_id, rgfa_offset));
        }
    };

    if (rgfa_seq_elems != nullptr) {
        // Input is a list
        for (const auto& rgfa_seq_elem : *rgfa_seq_elems) {
            update_rgfa_path(rgfa_seq_elem);
        }
    } else if (!unseekable) {
        // Input is from a seekable file on disk.
        gg.for_each_sequence_line_in_file(filename.c_str(), [&](gfak::sequence_elem s) {
                update_rgfa_path(s);
            });
    } else {
        // gg will have parsed the GFA file in the non-path part of the algorithm
        // No reading to do.
        for (const auto& seq_record : gg.get_name_to_seq()) {
            update_rgfa_path(seq_record.second);
        }
    }

    for (auto& path_offsets : path_map) {
        const string& name = path_offsets.first;
        int64_t rank = path_offsets.second.first;
        vector<pair<nid_t, int64_t>>& node_offsets = path_offsets.second.second;
        if (graph->has_path(name)) {
            cerr << "warning:[gfa_to_handle_graph] Ignoring rGFA tags for path " << name << " as a path with that name "
                 << "has already been imported from a P-line" << endl;
            continue;
        }
        // sort by offset
        sort(node_offsets.begin(), node_offsets.end(),
             [](const pair<nid_t, int64_t>& o1, const pair<nid_t, int64_t>& o2) { return o1.second < o2.second;});

        // make a path for each run of contiguous offsets
        int64_t prev_sequence_size;
        path_handle_t path_handle;
        for (int64_t i = 0; i < node_offsets.size(); ++i) {
            bool contiguous = i > 0 && node_offsets[i].second == node_offsets[i-1].second + prev_sequence_size;
            if (!contiguous) {
                // should probably detect and throw errors if overlap (as opposed to gap)
                string path_chunk_name = node_offsets[i].second == 0 ? name : Paths::make_subpath_name(name, node_offsets[i].second);
                path_handle = graph->create_path_handle(path_chunk_name);
            }
            handle_t step = graph->get_handle(node_offsets[i].first, false);
            graph->append_step(path_handle, step);
            prev_sequence_size = graph->get_length(step);
        }
    }    
}

static vector<gfak::sequence_elem> gfa_to_path_handle_graph_stream(istream& in, MutablePathMutableHandleGraph* graph,
                                                                   GFAIDMapInfo& id_map_info,
                                                                   int64_t max_rank) {
    if (!in) {
        throw std::ios_base::failure("error:[gfa_to_handle_graph] Couldn't open input stream");
    }

    bool has_rgfa_tags = false;
    string line_buffer; // can be quite big

    // store up rgfa nodes (without sequence), as there's no current way to avoid a second pass
    // to support them
    vector<gfak::sequence_elem> rgfa_seq_elems;

    function<void()> fall_back_to_disk = [&]() {
        string fb_name = temp_file::create();
        cerr << "warning:[gfa] Unable to stream GFA as it's not in canonical order.  Buffering to " << fb_name << endl;
        ofstream fb_file(fb_name);
        if (!fb_file) {
            throw runtime_error("error:[gfa] Could not open fallback gfa temp file: " + fb_name);
        }
        // put that last line back
        fb_file << line_buffer << "\n";
        // copy the rest of the file
        std::copy(istreambuf_iterator<char>(in),
                  istreambuf_iterator<char>(),
                  ostreambuf_iterator<char>(fb_file));
        fb_file.close();
        // read the file from disk
        gfak::GFAKluge gg;
        bool ret = gfa_to_handle_graph_on_disk(fb_name, graph, gg, id_map_info);
        gfa_to_handle_graph_add_paths(fb_name, nullptr, graph, gg, id_map_info);
        has_rgfa_tags = has_rgfa_tags || ret;
    };
    
    while (getline(in, line_buffer)) {
        if (!line_buffer.empty()) {
            // We mimic gfakluge behaviour by silently ignoring lines we don't parse
            if (line_buffer[0] == 'S') {
                tuple<string, string, tag_list_t> s_parse = GFAParser::parse_s(line_buffer);
                graph->create_handle(get<1>(s_parse), GFAParser::parse_sequence_id(get<0>(s_parse), id_map_info));
                if (get<2>(s_parse).size() >= 3) {
                    // We'll check for the 3 rGFA optional tags.  For now that means
                    // re-using some code based on gfakluge structures, unfortunately
                    // note: we only copy the name and tags, not the sequence
                    gfak::sequence_elem seq_elem;
                    seq_elem.name = get<0>(s_parse);
                    seq_elem.length = get<1>(s_parse).length();
                    for (const string& opt_tag : get<2>(s_parse)) {
                        vector<string> toks = split_delims(opt_tag, ":");
                        if (toks.size() == 3) {
                            gfak::opt_elem opt;
                            opt.key = toks[0];
                            opt.type = toks[1];
                            opt.val = toks[2];
                            seq_elem.opt_fields.push_back(opt);
                        }
                    }
                    int64_t rgfa_rank;
                    if (gfa_sequence_parse_rgfa_tags(seq_elem, nullptr, &rgfa_rank, nullptr) &&
                        rgfa_rank <= max_rank) {
                        rgfa_seq_elems.push_back(seq_elem);
                    }
                }                    
            } else if (line_buffer[0] == 'L') {
                tuple<string, bool, string, bool, range_t, tag_list_t> l_parse = GFAParser::parse_l(line_buffer);
                nid_t n1 = GFAParser::parse_sequence_id(get<0>(l_parse), id_map_info);
                nid_t n2 = GFAParser::parse_sequence_id(get<2>(l_parse), id_map_info);
                if (!graph->has_node(n1) || !graph->has_node(n2)) {
                    fall_back_to_disk();
                    break;
                }
                graph->create_edge(graph->get_handle(n1, get<1>(l_parse)),
                                   graph->get_handle(n2, get<3>(l_parse)));
            } else if (line_buffer[0] == 'P') {
                bool missing = false;
                // pass 1: make sure we have all the nodes in the graph
                GFAParser::scan_p(line_buffer, [&](const string& path_name,
                                                   int64_t step_rank,
                                                   const string& step_id,
                                                   bool step_is_reverse) {
                     if (step_rank >= 0) {
                         nid_t n = GFAParser::parse_sequence_id(step_id, id_map_info);
                         if (!graph->has_node(n)) {
                             missing = true;
                             return false;
                         }
                     }
                     return true;
                 });
                if (missing) {
                    fall_back_to_disk();
                    break;
                }
                path_handle_t path_handle;
                // pass 2: make the path
                GFAParser::scan_p(line_buffer, [&](const string& path_name,
                                                   int64_t step_rank,
                                                   const string& step_id,
                                                   bool step_is_reverse) {
                     if (step_rank <= 0) {
                         path_handle = graph->create_path_handle(path_name);
                     }
                     if (step_rank >= 0) {
                         nid_t n = GFAParser::parse_sequence_id(step_id, id_map_info);
                         graph->append_step(path_handle, graph->get_handle(n, step_is_reverse));
                     }
                     return true;
                 });
                
            }
        }
    }
    return rgfa_seq_elems;
}

/// Add listeners which let a GFA parser fill in a handle graph with nodes and edges.
static void add_graph_listeners(GFAParser& parser, MutableHandleGraph* graph) {
    parser.node_listeners.push_back([&parser, graph](nid_t id, const range_t& sequence, const GFAParser::tag_list_t& tags) {
        graph->create_handle(GFAParser::extract(sequence), id);
    });
    parser.edge_listeners.push_back([&parser, graph](nid_t from, bool from_is_reverse, nid_t to, bool to_is_reverse, const range_t& overlap, const tag_list_t& tags) {
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
static void add_path_listeners(GFAParser& parser, MutableHandleGraph* graph) {
    parser.path_listeners.push_back([&parser, graph](const string& name, const range_t& visits, const range_t& overlaps, const GFAParser::tag_list_t& tags) {
         path_handle = graph->create_path_handle(name);
         // TODO: Make sure overlaps has nothing important
         
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
static GFAParser::range_t take_range_until(GFAParser::cursor_t& cursor, const GFAParser::cursor_t& end, const char* end_chars, const char* parsing_state = nullptr) {
    auto start = cursor;
    while (cursor != end) {
        for (const char* stop_char = end_chars; *stop_char; ++stop_char) {
            if (*cursor == *stop_char) {
                // We found a stop character
                if (cursor == start) {
                     throw GFAFormatError("Expected nonempty value", cursor, parsing_state);
                }
                return range_t(start, cursor);
            }
        }
        ++cursor;
    }
    throw GFAFormatError("Expected terminator in " + std::string(end_chars), cursor, parsing_state);
}

/// Read a range, stopping at tab or end of line.
static GFAParser::range_t take_optional_range(GFAParser::cursor_t& cursor, const GFAParser::cursor_t& end, const char* parsing_state = nullptr) {
    auto start = cursor;
    while (cursor != end && *cursor != '\t') {
        ++cursor;
    }
    return range_t(start, cursor);
}

/// Read a range, stopping at tab or end of line.
/// Throw if it is empty.
static GFAParser::range_t take_range(GFAParser::cursor_t& cursor, const GFAParser::cursor_t& end, const char* parsing_state = nullptr) {
    GFAParser::range_t value = take_optional_range(cursor, end);
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
        if (*cursor == valse_value) {
            ++cursor;
            return false;
        }
    }
    throw GFAFormatError("Expected " + true_value + " or " + false_value, cursor, parsing_state);
}

/// Advance past a tab character that must exist.
static void take_tab(GFAParser::cursor_t& cursor, const GFAParser::cursor_t& end, const char* parsing_state = nullptr) {
    if (!take_optional_tab(cursor, end)) {
        throw GFAFormatError("Expected tab", cursor, parsing_state);
    }
}

GFAParser::tag_list_t GFAParser::parse_tags(const range_t& tag_range) {
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

tuple<string, range_t, tag_list_t> GFAParser::parse_s(const string& s_line) {
    auto cursor = s_line.begin();
    auto end = s_line.end();
    
    // Make sure we start with S
    take_character(cursor, end, 'S', "parsing S line start");
    take_tab(cursor, end "parsing S line");
    
    // Parse out the name
    string name = take_string(cursor, end, "parsing sequence name");
    take_tab(cursor, end "parsing end of sequence name");
    
    // Parse out the sequence
    range_t sequence = take_range(cursor, end, "parsing sequence");
    take_optional_tab(cursor, end, "parsing end of sequence");
    
    // Now we're either at the end or at the tab before the tags. Parse the tags.
    auto tags = GFAParser::parse_tags(range_t(cursor, end));
    
    return make_tuple(std::move(name), std::move(sequence), std::move(tags));
}

tuple<string, bool, string, bool, range_t, tag_list_t> GFAParser::parse_l(const string& l_line) {
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
    range_t overlaps = take_range(cursor, end, "parsing overlaps");
    take_optional_tab(cursor, end, "parsing end of overlaps");
    
    // Now we're either at the end or at the tab before the tags. Parse the tags.
    auto tags = GFAParser::parse_tags(range_t(cursor, end));
    
    return make_tuple(std::move(n1), n1_reverse, std::move(n2), n2_reverse, std::move(overlaps), std::move(tags));
}

tuple<string, range_t, range_t, tag_list_t> GFAParser::parse_p(const string& p_line) {
    auto cursor = p_line.begin();
    auto end = p_line.end();
    
    // Make sure we start with P
    take_character(cursor, end, 'P', "parsing P line start");
    take_tab(cursor, end, "parsing P line");
    
    // Grab the path name
    string path_name = take_string(cursor, end, "parsing path name");
    take_tab(cursor, end, "parsing end of path name");
    
    // Parse out the visits
    range_t visits = take_range(cursor, end, "parsing path visits");
    take_tab(cursor, end, "parsing end of path visits");
    
    // Parse out the overlaps
    range_t overlaps = take_range(cursor, end, "parsing overlaps");
    take_optional_tab(cursor, end, "parsing end of overlaps");
    
    // Now we're either at the end or at the tab before the tags. Parse the tags.
    auto tags = GFAParser::parse_tags(range_t(cursor, end));
    
    return make_tuple(std::move(path_name), std::move(visits), std::move(overlaps), std::move(tags));
}

void GFAParser::scan_p_visits(const range_t& visit_range,
                              function<bool(int64_t rank, const range_t& node_name, bool is_reverse)> visit_step) {
    
    auto cursor = visit_range.first;
    auto& end = visit_range.second;
    int64_t rank = 0;
    
    while (cursor != end) {
        // Until we run out of visit list range
        
        // Make a range for the visit node name
        range_t name_range = take_range_until(cursor, end, "+-", "parsing name of visited node");
        bool is_reverse = take_flag_character(cursor, end, '-', '+', "parsing orientation of visited node");
        
        if (!visit_step(rank, name_range, is_reverse)) {
            // We should stop looping
            return;
        }
        
        if (cursor != visit_range.end) {
            // Go past the comma separator
            take_character(cursor, end, ',', "parsing visit separator");
        }
        // And advance the rank for the next visit
        ++rank;
    }
    
    if (rank == 0) {
        // Nothing was visited. Show an empty path.
        visit_step(-1, range_t(visit_range.second, visit_range.second), false);
    }
}

void GFAParser::scan_p(const string& p_line,
                       function<bool(const string& path_name, int64_t rank, const range_t& node_name, bool is_reverse)> visit_step) {
    
    tuple<string, range_t, range_t, tag_list_t> p_parse = GFAParser::parse_p(p_line);
    auto& path_name = get<0>(p_parse);
    auto& path_visits = get<1>(p_parse);
    auto& overlaps = get<2>(p_parse);
    
    for(auto it = overlaps.first; it != overlaps.second; ++it) {
        if (*it != '*' && *it != ',') {
            // This overlap isn't just * or a list of *.
            // We can't handle it
            throw GFAFormatError("Path " + path_name + " has nontrivial overlaps and can't be handled");
        }
    }
    
    scan_p_visits(path_visits, [&](int64_t rank, const range_t& node_name, bool is_reverse) {
        return visit_step(path_name, rank, node_name, is_reverse);
    });
    
}

bool GFAParser::decode_rgfa_tags(const tag_list_t& tags,
                                 string* out_name = nullptr,
                                 int64_t* out_offset = nullptr,
                                 int64_t* out_rank = nullptr) {
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
                string tag_type = tag_name.substr(sep1 + 1, sep2);
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
        if (has_sn && has_so && has_sr) {
            break;
        }
    }
    return has_sn && has_so && has_sr;
}

nid_t GFAParser::parse_sequence_id(const string& str, GFAIDMapInfo& id_map_info) {
    
    auto found = id_map_info.name_to_id->find(str);
    if (found != id_map_info.name_to_id->end()) {
        // already in map, just return
        return found->second;
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
    id_map_info.name_to_id->emplace(str, node_id);

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

    bool has_rgfa_tags = false;
    string line_buffer; // can be quite big
    
    // We should be able to parse in 2 passes. One to make all the nodes, and
    // one to make all the things that reference nodes.
    int pass_number = 0;
    
    // We buffer lines we can't actually handle right now into this temporary file.
    string buffer_name;
    ofstream buffer_out_stream;
    auto save_line_until_node = [&](const string& missing_node_name) {
        if (pass_number > 0) {
            // We should only hit this on the first pass. If we hit it later we are missing a node.
            throw GFAFormatError("GFA file references missing node " + missing_node_name);
        }
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
    // buffer it if it can't.
    auto handle_line_if_ready = [&]() {
        if (!line_buffer.empty()) {
            switch(line_buffer[0]) {
            case 'H':
                // Header lines don't need anything done
                break;
            case 'S':
                // Sequence lines can always be handled right now
                {
                    tuple<string, GFAParser::range_t, tag_list_t> s_parse = GFAParser::parse_s(line_buffer);
                    auto& node_name = get<0>(s_parse);
                    auto& sequence_range = get<1>(s_parse);
                    auto& tags = get<2>(s_parse);
                    // TODO: enforce ID uniqueness here
                    nid_t assigned_id = GFAParser::parse_sequence_id(node_name, this->id_map());
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
                            // Buffer this visit.
                            auto& visit_queue = get<2>(found->second);
                            auto& next_offset = get<1>(found->second);
                            if (next_offset == rgfa_offset_on_path) {
                                // It's safe to dispatch this visit right now since it's the next one expected along the path.
                                for (auto& rgfa_listeners : this->path_listeners) {
                                    // Tell all the listener functions about this visit
                                    listener(assigned_id, rgfa_offset_on_path, rgfa_path_name, rgfa_path_rank);
                                }
                                // Advance the offset by the sequence length;
                                next_offset += GFAParser::length(sequence_range);
                                while (!visit_queue.empty() && next_offset == get<0>(visit_queue.top())) {
                                    // The lowest-offset queued visit can be handled now because it abuts what we just did.
                                    // Grab the visit.
                                    auto& visit = visit_queue.top();
                                    for (auto& rgfa_listeners : this->path_listeners) {
                                        // Tell all the listener functions about this visit
                                        listener(get<1>(visit), get<0>(visit), rgfa_path_name, rgfa_path_rank);
                                    }
                                    // Advance the offset by the sequence length;
                                    next_offset += get<2>(visit);
                                    // And pop the visit off
                                    visit_queue.pop();
                                }
                            } else {
                                // Add this visit to the heap so we can handle it when we find the missing visits.
                                visit_buffer.emplace(rgfa_offset_on_path, assigned_id, GFAParser::length(sequence_range));
                            }
                        }
                    }
                }
                break;
            case 'L':
                // Edges can be handled if their nodes exist already
                {
                    tuple<string, bool, string, bool, range_t, tag_list_t> l_parse = GFAParser::parse_l(line_buffer);
                    
                    // We only get these IDs if they have been seen already as nodes
                    nid_t n1 = GFAParser::find_existing_sequence_id(get<0>(l_parse), this->id_map());
                    if (!n1) {
                        save_line_until_node(get<0>(l_parse));
                        break;
                    }
                    nid_t n2 = GFAParser::find_existing_sequence_id(get<2>(l_parse), this->id_map());
                    if (!n2) {
                        save_line_until_node(get<2>(l_parse));
                        break;
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
                    // pass 1: make sure we have all the nodes in the graph
                    GFAParser::scan_p(line_buffer, [&](const string& path_name,
                                                       int64_t step_rank,
                                                       const string& step_id,
                                                       bool step_is_reverse) {
                         if (step_rank >= 0) {
                             nid_t n = GFAParser::find_existing_sequence_id(step_id, this->id_map());
                             if (!n) {
                                missing = true;
                                missing_name = step_id;
                                return false;
                             }
                         }
                         return true;
                    });
                    if (missing) {
                        save_line_until_node(missing_name);
                        break;
                    }
                    
                    // Re-parse the pieces of the line.
                    // TODO: Deduplicate with scan
                    tuple<string, range_t, range_t, tag_list_t> p_parse = GFAParser::parse_p(line_buffer);
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
    };
    
    // We have a function to do a pass over a file of lines.
    auto process_lines_in_stream = [&](istream& in_stream) {
        while (getline(in, line_buffer)) {
            // For each line in the input file
            if (!line_buffer.empty()) {
                // Handle all lines in the stream that we can handle now.
                handle_line_if_ready(line_buffer); 
            }
        }
    }
    
    process_lines_in_stream(in);
    
    if (!buffer_name.empty()) {
        // We also have lines in the buffer to handle.
        buffer_out_stream.close();
        pass_number++;
        ifstream buffer_in_stream(buffer_name);
        process_lines_in_stream(buffer_in_stream);
        
        // Clean up buffer before returning
        // TODO: Who takes care of this on GFA format error?
        buffer_in_stream.close();
        unlink(buffer_name);
    }
    
    // Run through any rGFA paths that don't start at 0 or have gaps. 
    for (auto& kv : rgfa_path_cache) {
        auto& rgfa_path_name = kv.first;
        auto& rgfa_path_rank = get<0>(kv.second);
        auto& visit_queue = get<2>(kv.second);
        
        while (!visit_queue.empty()) {
            // Grab the visit.
            auto& visit = visit_queue.top();
            for (auto& rgfa_listeners : this->path_listeners) {
                // Tell all the listener functions about this visit
                listener(get<1>(visit), get<0>(visit), rgfa_path_name, rgfa_path_rank);
            }
            // And pop the visit off
            visit_queue.pop();
        }
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

GFAFormatError::GFAFormatError(const GFAParser::cursor_t& position, const string& message, const char* parsing_state) : GFAFormatError(message + (parsing_state ? (" while " + std::string(parsing_state)) : "")), has_position(true), position(position) {
    // Nothing to do!
}


}
}
