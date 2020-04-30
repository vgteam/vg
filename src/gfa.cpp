#include "gfa.hpp"
#include <gfakluge.hpp>

// Use sonLib pinch graphs
#include <sonLib/stPinchGraphs.h>

#include <structures/union_find.hpp>

namespace vg {

using namespace std;
using namespace gfak;

void graph_to_gfa(const VG* graph, ostream& out) {
  GFAKluge gg;
  gg.set_version(1.0);
  for (auto h : gg.get_header()){
    out << h.second.to_string();
  }

    // TODO moving to GFAKluge
    // problem: protobuf longs don't easily go to strings....
    
    graph->for_each_node([&](const Node* n) {
        sequence_elem s_elem;
        // Fill seq element for a node
        s_elem.name = to_string(n->id());
        s_elem.sequence = n->sequence();
        out << s_elem.to_string_1() << endl;
        //gg.add_sequence(s_elem);
    });
    
    auto& pathmap = graph->paths._paths;
    for (auto p : pathmap){
        path_elem p_elem;
        p_elem.name = p.first;
        for (auto m : p.second){
            p_elem.segment_names.push_back( std::to_string(m.node_id()) );
            p_elem.orientations.push_back( !m.is_reverse() );
            const Node* n = graph->get_node( m.node_id() );
            stringstream cigaro;
            //cigaro << n->sequence().size() << (p.mapping(m_ind.position().is_reverse()) ? "M" : "M");
            cigaro << n->sequence().size() << (m.is_reverse() ? "M" : "M");
            p_elem.overlaps.push_back( cigaro.str() );
        }
        out << p_elem.to_string_1() << endl;
        //gg.add_path(p_elem.name, p_elem);
    }

    graph->for_each_edge([&](const Edge* e) {
        edge_elem ee;
        ee.type = 1;
        ee.source_name = to_string(e->from());
        ee.sink_name = to_string(e->to());
        ee.source_orientation_forward = ! e->from_start();
        ee.sink_orientation_forward =  ! e->to_end();
        ee.alignment = std::to_string(e->overlap()) + "M";
        
        if (e->from_start() && (e->to_end() || e->to() < e->from())) {
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
        //gg.add_edge(ee.source_name, ee);
        //link_elem l;
        //l.source_name = to_string(e->from());
        //l.sink_name = to_string(e->to());
        //l.source_orientation_forward = ! e->from_start();
        //l.sink_orientation_forward =  ! e->to_end();
        //l.cigar = std::to_string(e->overlap()) + "M";
        //gg.add_link(l.source_name, l);
    });
    //gg.output_to_stream(cout);
}

}
