#include "gfa_to_handle.hpp"

namespace vg {
namespace algorithms {

nid_t parse_gfa_sequence_id(const string& str) {
    if (any_of(str.begin(), str.end(), [](char c) { return !isdigit(c); })) {
        throw GFAFormatError("error:[gfa_to_handle_graph] Could not parse sequence ID '" + str + "'. GFA sequence IDs must be integers >= 1.");
    }
    nid_t node_id = stoll(str);
    if (node_id <= 0) {
        throw GFAFormatError("error:[gfa_to_handle_graph] Could not parse sequence ID '" + str + "'. GFA sequence IDs must be integers >= 1.");
    }
    return node_id;
}


void validate_gfa_edge(const gfak::edge_elem& e) {
    string not_blunt = ("error:[gfa_to_handle_graph] Can only load blunt-ended GFAs. "
        "Try \"bluntifying\" your graph with a tool like <https://github.com/hnikaein/stark>, or "
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

string process_raw_gfa_path_name(const string& path_name_raw)  {
    string processed = path_name_raw;
    processed.erase(remove_if(processed.begin(), processed.end(),
                              [](char c) { return isspace(c); }),
                    processed.end());
    return processed;
}

void gfa_to_handle_graph_in_memory(istream& in, MutableHandleGraph* graph,
                                   gfak::GFAKluge& gg) {
    if (!in) {
        throw std::ios_base::failure("error:[gfa_to_handle_graph] Couldn't open input stream");
    }
    gg.parse_gfa_file(in);
    
    // create nodes
    for (const auto& seq_record : gg.get_name_to_seq()) {
        graph->create_handle(seq_record.second.sequence, parse_gfa_sequence_id(seq_record.first));
    }
    
    // create edges
    for (const auto& links_record : gg.get_seq_to_edges()) {
        for (const auto& edge : links_record.second) {
            validate_gfa_edge(edge);
            // note: we're counting on implementations de-duplicating edges
            handle_t a = graph->get_handle(parse_gfa_sequence_id(edge.source_name), !edge.source_orientation_forward);
            handle_t b = graph->get_handle(parse_gfa_sequence_id(edge.sink_name), !edge.sink_orientation_forward);
            graph->create_edge(a, b);
        }
    }
}

void gfa_to_handle_graph_on_disk(const string& filename, MutableHandleGraph* graph,
                                 bool try_id_increment_hint, gfak::GFAKluge& gg) {
    
    // adapted from
    // https://github.com/vgteam/odgi/blob/master/src/gfa_to_handle.cpp
    
    if (try_id_increment_hint) {
        
        // find the minimum ID
        nid_t min_id = numeric_limits<nid_t>::max();
        gg.for_each_sequence_line_in_file(filename.c_str(), [&](gfak::sequence_elem s) {
            min_id = std::min(min_id, parse_gfa_sequence_id(s.name));
        });
        
        if (min_id != numeric_limits<nid_t>::max()) {
            // we found the min, set it as the increment
            graph->set_id_increment(min_id);
        }
    }
    
    // add in all nodes
    gg.for_each_sequence_line_in_file(filename.c_str(), [&](gfak::sequence_elem s) {
        
        graph->create_handle(s.sequence, parse_gfa_sequence_id(s.name));
        
    });
    
    // add in all edges
    gg.for_each_edge_line_in_file(filename.c_str(), [&](gfak::edge_elem e) {
        validate_gfa_edge(e);
        handle_t a = graph->get_handle(parse_gfa_sequence_id(e.source_name), !e.source_orientation_forward);
        handle_t b = graph->get_handle(parse_gfa_sequence_id(e.sink_name), !e.sink_orientation_forward);
        graph->create_edge(a, b);
    });
}


/// Parse nodes and edges and load them into the given GFAKluge.
/// If the input is a seekable file, filename will be filled in and unseekable will be nullptr.
/// If the input is not a seekable file, filename may be filled in, and unseekable will be set to a stream to read from.
void gfa_to_handle_graph_load_graph(const string& filename, istream* unseekable, MutableHandleGraph* graph,
                                    bool try_id_increment_hint, gfak::GFAKluge& gg) {
    
    if (graph->get_node_count() > 0) {
        throw invalid_argument("error:[gfa_to_handle_graph] Must parse GFA into an empty graph");
    }
    
    if (!unseekable) {
        // Do the from-disk path
        gfa_to_handle_graph_on_disk(filename, graph, try_id_increment_hint, gg);
    } else {
        // Do the path for streams
        
        if (try_id_increment_hint) {
            // The ID increment hint can't be done.
            cerr << "warning:[gfa_to_handle_graph] Skipping node ID increment hint because input stream for GFA does not support seeking. "
                 << "If performance suffers, consider using an alternate graph implementation or reading GFA from hard disk." << endl;
        }
        
        gfa_to_handle_graph_in_memory(*unseekable, graph, gg);
    }
}

/// After the given GFAKluge has been populated with nodes and edges, load path information.
/// If the input is a seekable file, filename will be filled in and unseekable will be nullptr.
/// If the input is not a seekable file, filename may be filled in, and unseekable will be set to a stream to read from.
void gfa_to_handle_graph_add_paths(const string& filename, istream* unseekable, MutablePathHandleGraph* graph,
                                   gfak::GFAKluge& gg) {
                                   
                                   
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
            handle_t step = graph->get_handle(parse_gfa_sequence_id(node_id), is_rev);
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
                handle_t step = graph->get_handle(parse_gfa_sequence_id(path_record.second.segment_names.at(i)),
                                                  !path_record.second.orientations.at(i));
                graph->append_step(path, step);
            }
        }
    }
    
    
}

void gfa_to_handle_graph(const string& filename, MutableHandleGraph* graph,
                         bool try_from_disk, bool try_id_increment_hint) {

    // What stream should we read from (isntead of opening the file), if any?
    istream* unseekable = nullptr;
    
    // If we open a file, it will live here.
    unique_ptr<ifstream> opened;
    
    if (filename == "-") {
        // Read from standard input
        unseekable = &cin;
    } else if (!try_from_disk) {
        // The file may be seekable actually, but we don't want to use the
        // seekable-file codepath for some reason.
        opened = make_unique<ifstream>(filename);
        if (!opened) {
            throw std::ios_base::failure("error:[gfa_to_handle_graph] Couldn't open file " + filename);
        }
        unseekable = opened.get();
    }
    
    gfak::GFAKluge gg;
    gfa_to_handle_graph_load_graph(filename, unseekable, graph, try_id_increment_hint, gg);
}


void gfa_to_path_handle_graph(const string& filename, MutablePathMutableHandleGraph* graph,
                              bool try_from_disk, bool try_id_increment_hint) {
    
    
    // What stream should we read from (isntead of opening the file), if any?
    istream* unseekable = nullptr;
    
    // If we open a file, it will live here.
    unique_ptr<ifstream> opened;
    
    if (filename == "-") {
        // Read from standard input
        unseekable = &cin;
    } else if (!try_from_disk) {
        // The file may be seekable actually, but we don't want to use the
        // seekable-file codepath for some reason.
        opened = make_unique<ifstream>(filename);
        if (!opened) {
            throw std::ios_base::failure("error:[gfa_to_handle_graph] Couldn't open file " + filename);
        }
        unseekable = opened.get();
    }
    
    gfak::GFAKluge gg;
    gfa_to_handle_graph_load_graph(filename, unseekable, graph, try_id_increment_hint, gg);
    
    // TODO: Deduplicate everything other than this line somehow.
    gfa_to_handle_graph_add_paths(filename, unseekable, graph, gg);
}

void gfa_to_path_handle_graph_in_memory(istream& in,
                                        MutablePathMutableHandleGraph* graph) {
    gfak::GFAKluge gg;
    gfa_to_handle_graph_load_graph("", &in, graph, false, gg);
    gfa_to_handle_graph_add_paths("", &in, graph, gg);
    
}

}
}
