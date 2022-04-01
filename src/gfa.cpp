#include "gfa.hpp"
#include <gfakluge.hpp>
#include "utility.hpp"
#include "path.hpp"
#include <sstream>

namespace vg {

using namespace std;
using namespace gfak;

/// Determine if a path should be written as a GFA W line or a GFA P line.
static bool should_write_as_w_line(const PathHandleGraph* graph, path_handle_t path_handle);
/// Write out a W line for a path. Uses a map to keep track of fake offset
/// ranges used to distinguish multiple phase blocks on a haplotype, since GFA
/// doesn't support them.
static void write_w_line(const PathHandleGraph* graph, ostream& out, path_handle_t path_handle, unordered_map<tuple<string, int64_t, string>, size_t>& last_phase_block_end);

void graph_to_gfa(const PathHandleGraph* graph, ostream& out, const set<string>& rgfa_paths,
                  bool rgfa_pline, bool use_w_lines) {
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
    
    // Sort the paths by name, making sure to treat subpath coordinates numerically
    vector<path_handle_t> path_handles;
    graph->for_each_path_handle([&](const path_handle_t& h) {
            path_handles.push_back(h);
        });
    std::sort(path_handles.begin(), path_handles.end(), [&](const path_handle_t& p1, const path_handle_t& p2) {
            string n1 = graph->get_path_name(p1);
            string n2 = graph->get_path_name(p2);
            auto s1 = Paths::parse_subpath_name(n1);
            auto s2 = Paths::parse_subpath_name(n2);
            if (!get<0>(s1) || !get<0>(s2)) {
                return n1 < n2;
            } else if (get<1>(s1) < get<1>(s2)) {
                return true;
            } else if (get<1>(s1) == get<1>(s2)) {
                if (get<2>(s1) < get<2>(s2)) {
                    return true;
                } else if (get<2>(s1) == get<2>(s2)) {
                    return get<3>(s1) < get<3>(s2);
                }
            }
            return false;
        });

    vector<path_handle_t> w_line_paths;

    // Paths as P-lines
    for (const path_handle_t& h : path_handles) {
        path_elem p_elem;
        p_elem.name = graph->get_path_name(h);
        if (rgfa_pline || !rgfa_paths.count(p_elem.name)) {
            if (use_w_lines && should_write_as_w_line(graph, h)) {
                w_line_paths.push_back(h);
            } else {
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
    }
    
    // Paths as W-lines
    {
        unordered_map<tuple<string, int64_t, string>, size_t> last_phase_block_end;
        for (const path_handle_t& h : w_line_paths) {
            write_w_line(graph, out, h, last_phase_block_end);
        }
    }

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

bool should_write_as_w_line(const PathHandleGraph* graph, path_handle_t path_handle) {
    auto sense = graph->get_sense(path_handle);
    // Haplotype and reference sense paths both have good W line descriptions.
    // TODO: how to tell them apart?
    return sense == PathMetadata::SENSE_HAPLOTYPE || sense == PathMetadata::SENSE_REFERENCE;
}

void write_w_line(const PathHandleGraph* graph, ostream& out, path_handle_t path_handle, unordered_map<tuple<string, int64_t, string>, size_t>& last_phase_block_end) {
    // Extract the path metadata
    string sample = graph->get_sample_name(path_handle);
    string contig = graph->get_locus_name(path_handle);
    int64_t hap_index = graph->get_haplotype(path_handle);
    int64_t phase_block = graph->get_phase_block(path_handle);
    auto subrange = graph->get_subrange(path_handle);
    size_t start_offset = 0;
    size_t end_offset = 0;
    if (subrange != PathMetadata::NO_SUBRANGE) {
        start_offset = subrange.first;
        if (subrange.second != PathMetadata::NO_END_POSITION) {
            end_offset = subrange.second;
        }
    }
    
    if (hap_index == PathMetadata::NO_HAPLOTYPE) {
        // No haplotype is actually assigned here.
        // We really shouldn't have paths with it assigned and not assigned, so assign it 0 and make the sample haploid.
        // TODO: check for collisions somehoe?
        hap_index = 0;
    }
     
    // Get the path length.
    // TODO: sniff if the graph has this cached somehow?
    size_t path_length = 0 ;
    graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
            path_length += graph->get_length(graph->get_handle_of_step(step_handle));
        });

    if (end_offset != 0 && start_offset + path_length != end_offset) {
        cerr << "[gfa] warning: incorrect end offset (" << end_offset << ") extracted from from path name " << path_name
             << ", using " << (start_offset + path_length) << " instead" << endl;
    }
    
    // See if we need to bump along the start offset to avoid collisions of phase blocks
    auto key = std::make_tuple<string, int64_t, string>(sample, hap_index, contig);
    auto& phase_block_end_cursor = last_phase_block_end[key];
    if (phase_block_end_cursor != 0) {
        if (start_offset != 0) {
            cerr << "[gfa] error: cannot write multiple phase blocks on a sample, haplotyope, and contig in GFA format"
                 << " when paths already have subranges. Fix path " << graph->get_path_name(path_handle) << endl;
            exit(1);
        }
        // Budge us to after the last thing and budge the cursor to after us, plus a gigabase.
        // TODO: How are we ever going to round-trip this back to phase blocks?
        // Encode some very high bits in the offsets even for the first path?
        start_offset += phase_block_end_cursor + 1000000000;
        phase_block_end_cursor += path_length;
    }

    out << "W\t" << sample << "\t" << hap_index << "\t" << contig << "\t" << start_offset << "\t" << (start_offset + path_length) << "\t";

    graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
            handle_t handle = graph->get_handle_of_step(step_handle);
            out << (graph->get_is_reverse(handle) ? "<" : ">") << graph->get_id(handle);
        });
    out << "\n";
}

}
