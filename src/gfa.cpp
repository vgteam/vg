#include "gfa.hpp"
#include <gfakluge.hpp>

namespace vg {

using namespace std;
using namespace gfak;

void graph_to_gfa(const PathHandleGraph* graph, ostream& out, const set<string>& rgfa_paths) {
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
        out << endl;
        return true;
    });
    
    //Go through each path
    graph->for_each_path_handle([&](const path_handle_t& h) {
        path_elem p_elem;
        p_elem.name = graph->get_path_name(h);
        if (!rgfa_paths.count(p_elem.name)) {
            graph->for_each_step_in_path(h, [&](const step_handle_t& ph) {

                    handle_t step_handle = graph->get_handle_of_step(ph);

                    p_elem.segment_names.push_back( std::to_string(graph->get_id(step_handle)) );
                    p_elem.orientations.push_back( !graph->get_is_reverse(step_handle) );
                    stringstream cigaro;
                    //cigaro << n->sequence().size() << (p.mapping(m_ind.position().is_reverse()) ? "M" : "M");
                    cigaro << graph->get_sequence(step_handle).size() << (graph->get_is_reverse(step_handle) ? "M" : "M");
                    p_elem.overlaps.push_back( cigaro.str() );
                    return true;
                });
            //gg.add_path(p_elem.name, p_elem);
            out << p_elem.to_string_1() << endl;
        }
        return true;
    });


    graph->for_each_edge([&](const edge_t& h) {
        edge_elem ee;
        ee.type = 1;
        //TODO: I'm guessing this is what it wants? 
        ee.source_name = to_string(graph->get_id(h.first));
        ee.sink_name = to_string(graph->get_id(h.second));
        ee.source_orientation_forward = ! graph->get_is_reverse(h.first);
        ee.sink_orientation_forward =  ! graph->get_is_reverse(h.second);

        ee.alignment = "0M";
        
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
        
        out << ee.to_string_1() << endl;;
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

}
