#include "gfa_to_handle.hpp"

namespace vg {
namespace algorithms {

nid_t parse_gfa_sequence_id(const string& str) {
    if (any_of(str.begin(), str.end(), [](char c) { return !isdigit(c); })) {
        throw runtime_error("error:[gfa_to_handle_graph] Could not parse sequence ID '" + str + "'. GFA sequence IDs must be integers >= 1.");
    }
    nid_t node_id = stoll(str);
    if (node_id <= 0) {
        throw runtime_error("error:[gfa_to_handle_graph] Could not parse sequence ID '" + str + "'. GFA sequence IDs must be integers >= 1.");
    }
    return node_id;
}


void validate_gfa_edge(const gfak::edge_elem& e) {
    if (e.source_begin != e.source_end || e.sink_begin != 0 || e.sink_end != 0) {
        throw runtime_error("error:[gfa_to_handle_graph] Only can load blunt ended GFAs. Found edge with an overlay: " + e.source_name + "[" + to_string(e.source_begin) + ":" + to_string(e.source_end) + "] -> " + e.sink_name + "[" + to_string(e.sink_begin) + ":" + to_string(e.sink_end) + "]");
    }
    if (!(e.alignment == "0M" || e.alignment == "*" || e.alignment.empty())) {
        throw runtime_error("error:[gfa_to_handle_graph] Only can load blunt ended GFAs. Found edge with a non-null alignment '" + e.alignment + "'.");
    }
    if (e.source_name.empty()) {
        throw runtime_error("error:[gfa_to_handle_graph] Found edge record with missing source name");
    }
    if (e.sink_name.empty()) {
        throw runtime_error("error:[gfa_to_handle_graph] Found edge record with missing sink name");
    }
}

string process_raw_gfa_path_name(const string& path_name_raw)  {
    string processed = path_name_raw;
    processed.erase(remove_if(processed.begin(), processed.end(),
                              [](char c) { return isspace(c); }),
                    processed.end());
    return processed;
}

void gfa_to_handle_graph_in_memory(const string& filename, MutableHandleGraph* graph,
                                   gfak::GFAKluge& gg) {
    if (filename == "-") {
        gg.parse_gfa_file(cin);
    }
    else {
        ifstream in(filename);
        if (!in) {
            throw runtime_error("error:[gfa_to_handle_graph] Couldn't open file " + filename);
        }
        gg.parse_gfa_file(in);
    }
    
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

void gfa_to_handle_graph_internal(const string& filename, MutableHandleGraph* graph,
                                  bool try_from_disk, bool try_id_increment_hint,
                                  gfak::GFAKluge& gg) {
    
    if (graph->get_node_count() > 0) {
        throw runtime_error("error:[gfa_to_handle_graph] Must parse GFA into an empty graph");
    }
    
    if (try_id_increment_hint && filename == "-") {
        cerr << "warning:[gfa_to_handle_graph] Skipping node ID increment hint because input stream for GFA does not support seeking. If performance suffers, consider using an alternate graph implementation or reading GFA from hard disk." << endl;
    }
    
    if (try_from_disk && filename == "-") {
        cerr << "warning:[gfa_to_handle_graph] From-disk GFA reading was selected, but input stream is not seekable. Using in-memory conversion algorithm instead. If memory performance suffers, consider reading GFA from hard disk." << endl;
    }
    
    if (try_from_disk && filename != "-") {
        gfa_to_handle_graph_on_disk(filename, graph, try_id_increment_hint, gg);
    }
    else {
        gfa_to_handle_graph_in_memory(filename, graph, gg);
    }
}

void gfa_to_handle_graph(const string& filename, MutableHandleGraph* graph,
                         bool try_from_disk, bool try_id_increment_hint) {

    gfak::GFAKluge gg;
    gfa_to_handle_graph_internal(filename, graph, try_from_disk, try_id_increment_hint, gg);
}


void gfa_to_path_handle_graph(const string& filename, MutablePathMutableHandleGraph* graph,
                              bool try_from_disk, bool try_id_increment_hint) {
    
    
    if (graph->get_path_count() > 0) {
        throw runtime_error("error:[gfa_to_handle_graph] Must parse GFA into an empty graph");
    }
    
    gfak::GFAKluge gg;
    gfa_to_handle_graph_internal(filename, graph, try_from_disk, try_id_increment_hint, gg);
    
    if (try_from_disk && filename != "-") {
        
        
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
    }
    else {
        
        // gg will have parsed the GFA file in the non-path part of the algorithm
        
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

}
}
