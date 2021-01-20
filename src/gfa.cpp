#include "gfa.hpp"
#include <gfakluge.hpp>

namespace vg {

using namespace std;
using namespace gfak;

void graph_to_gfa(const PathHandleGraph* graph, ostream& out) {
  GFAKluge gg;
  gg.set_version(1.0);
  for (auto h : gg.get_header()){
    out << h.second.to_string();
  }

    // TODO moving to GFAKluge
    // problem: protobuf longs don't easily go to strings....

    
    //Go through each node in the graph
    graph->for_each_handle([&](const handle_t& h) {
        sequence_elem s_elem;
        // Fill seq element for a node
        s_elem.name = to_string(graph->get_id(h));
        s_elem.sequence = graph->get_sequence(h);
        out << s_elem.to_string_1() << endl;
        //gg.add_sequence(s_elem);
        return true;
    });
    
    //Go through each path
    graph->for_each_path_handle([&](const path_handle_t& h) {
        path_elem p_elem;
        p_elem.name = graph->get_path_name(h);
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
