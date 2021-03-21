#include "gfa.hpp"
#include <gfakluge.hpp>
#include "utility.hpp"
#include "path.hpp"
#include <sstream>

namespace vg {

using namespace std;
using namespace gfak;

static bool write_w_line(const PathHandleGraph* graph, ostream& out, const string& wline_sep,
                         path_handle_t path_handle);

void graph_to_gfa(const PathHandleGraph* graph, ostream& out, const set<string>& rgfa_paths,
                  bool rgfa_pline, const string& wline_sep) {
    GFAKluge gg;
    gg.set_version(1.0);
    for (auto h : gg.get_header()){
        out << h.second.to_string();
    }

    // TODO moving to GFAKluge
    // problem: protobuf longs don't easily go to strings....
    
    //Compute the rGFA tags of given paths (todo: support non-zero ranks)
    unordered_map<nid_t, pair<path_handle_t, size_t>> node_offsets;
    for (const string& path_name : rgfa_paths) {
        path_handle_t path_handle = graph->get_path_handle(path_name);
        size_t offset = 0;
        graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
                handle_t handle = graph->get_handle_of_step(step_handle);
                nid_t node_id = graph->get_id(handle);
                if (graph->get_is_reverse(handle)) {
                    stringstream ss;
                    ss << "error [gfa]: unable to write rGFA tags for path " << path_name << " because node "
                       << node_id << " is traversed on its reverse strand.  rGFA only supports the forward strand." << endl;
                    throw runtime_error(ss.str());
                }
                if (node_offsets.count(node_id)) {
                    cerr << "warning [gfa]: multiple selected rgfa paths found on node " << node_id << ": keeping tags for "
                         << graph->get_path_name(node_offsets[node_id].first) << " and ignoring those for " << path_name << endl;
                } else {
                    node_offsets[node_id] = make_pair(path_handle, offset);
                }
                offset += graph->get_length(handle);
            });
    }
  
    //Go through each node in the graph
    graph->for_each_handle([&](const handle_t& h) {
        sequence_elem s_elem;
        // Fill seq element for a node
        nid_t node_id = graph->get_id(h);
        s_elem.name = to_string(node_id);
        s_elem.sequence = graph->get_sequence(h);
        out << s_elem.to_string_1();
        //gg.add_sequence(s_elem);
        auto it = node_offsets.find(node_id);
        if (it != node_offsets.end()) {
            // add rGFA tags
            out << "\t" << "SN:Z:" << graph->get_path_name(it->second.first)
                << "\t" << "SO:i:" << it->second.second
                << "\t" << "SR:i:0"; // todo: support non-zero ranks?
        }
        out << "\n"; // Writing `std::endl` would flush the buffer.
        return true;
    });
    
    //Go through each path
    graph->for_each_path_handle([&](const path_handle_t& h) {
        path_elem p_elem;
        p_elem.name = graph->get_path_name(h);
        if (rgfa_pline || !rgfa_paths.count(p_elem.name)) {
            bool wrote_w_line = write_w_line(graph, out, wline_sep, h);
            if (!wrote_w_line) {
                graph->for_each_step_in_path(h, [&](const step_handle_t& ph) {
            
                    handle_t step_handle = graph->get_handle_of_step(ph);

                    p_elem.segment_names.push_back( std::to_string(graph->get_id(step_handle)) );
                    p_elem.orientations.push_back( !graph->get_is_reverse(step_handle) );
                    return true;
                });
                p_elem.overlaps.push_back("*");
                //gg.add_path(p_elem.name, p_elem);
                out << p_elem.to_string_1() << "\n";
            }
        }
    });

    graph->for_each_edge([&](const edge_t& h) {
        edge_elem ee;
        ee.type = 1;
        //TODO: I'm guessing this is what it wants? 
        ee.source_name = to_string(graph->get_id(h.first));
        ee.sink_name = to_string(graph->get_id(h.second));
        ee.source_orientation_forward = ! graph->get_is_reverse(h.first);
        ee.sink_orientation_forward =  ! graph->get_is_reverse(h.second);

        ee.alignment = "*";
        
        if (graph->get_is_reverse(h.first) && (graph->get_is_reverse(h.second) || graph->get_id(h.second) < graph->get_id(h.first))) {
            // Canonicalize edges to be + orientation first if possible, and
            // then low-ID to high-ID if possible, for testability. This edge
            // needs to flip.
            
            // Swap the nodes
            std::swap(ee.source_name, ee.sink_name);
            // Swap the orientations
            std::swap(ee.source_orientation_forward, ee.sink_orientation_forward);
            // Reverse the orientations
            ee.source_orientation_forward = !ee.source_orientation_forward;
            ee.sink_orientation_forward = !ee.sink_orientation_forward;
        }
        
        out << ee.to_string_1() << "\n"; // Writing `std::endl` would flush the buffer.
        return true;
        //gg.add_edge(ee.source_name, ee);
        //link_elem l;
        //l.source_name = to_string(e->from());
        //l.sink_name = to_string(e->to());
        //l.source_orientation_forward = ! e->from_start();
        //l.sink_orientation_forward =  ! e->to_end();
        //l.cigar = std::to_string(e->overlap()) + "M";
        //gg.add_link(l.source_name, l);
    }, false);
    //gg.output_to_stream(cout);
}

bool write_w_line(const PathHandleGraph* graph, ostream& out, const string& wline_sep,
                  path_handle_t path_handle) {
    if (wline_sep.empty()) {
        return false;
    }
    string path_name = graph->get_path_name(path_handle);
    vector<string> toks = split_delims(path_name, wline_sep);
    if (toks.size() < 3) {
        return false;
    }
    string& sample = toks[0];
    size_t hap_index;
    try {
        hap_index = stol(toks[1]);
    } catch(...) {
        return false;
    }
    auto subpath_parse = Paths::parse_subpath_name(toks[2]);
    string contig;
    size_t start_offset = 0;
    size_t end_offset = 0;
    if (get<0>(subpath_parse) == true) {
        contig = get<1>(subpath_parse);
        start_offset = get<2>(subpath_parse);
        end_offset = get<3>(subpath_parse);
    } else {
        contig = toks[2];
    }

    size_t path_length = 0 ;
    graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
            path_length += graph->get_length(graph->get_handle_of_step(step_handle));
        });

    if (end_offset != 0 && start_offset + path_length != end_offset) {
        cerr << "[gfa] warning: incorrect end offset (" << end_offset << ") extracted from from path name " << path_name
             << ", using " << (start_offset + path_length) << " instead" << endl;
    }

    out << "W\t" << sample << "\t" << hap_index << "\t" << contig << "\t" << start_offset << "\t" << (start_offset + path_length) << "\t";

    graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
            handle_t handle = graph->get_handle_of_step(step_handle);
            out << (graph->get_is_reverse(handle) ? "<" : ">") << graph->get_id(handle);
        });
    out << "\n";
    return true;
}

}
