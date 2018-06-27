#include "gfa.hpp"
#include <gfakluge.hpp>

namespace vg {

using namespace std;
using namespace gfak;

void gfa_to_graph(istream& in, VG* graph, bool only_perfect_match) {
    // c++... split...
    // for line in stdin
    string line;
    auto too_many_fields = [&line]() {
        cerr << "[vg] error: too many fields in line " << endl << line << endl;
        exit(1);
    };

    bool reduce_overlaps = false;
    GFAKluge gg;
    gg.parse_gfa_file(in);

    map<string, sequence_elem, custom_key> name_to_seq = gg.get_name_to_seq();
    map<std::string, vector<edge_elem> > seq_to_edges = gg.get_seq_to_edges();
    map<string, sequence_elem>::iterator it;
    id_t curr_id = 1;
    map<string, id_t> id_names;
    std::function<id_t(const string&)> get_add_id = [&](const string& name) -> id_t {
        if (is_number(name)) {
            return std::stol(name);
        } else {
            auto id = id_names.find(name);
            if (id == id_names.end()) {
                id_names[name] = curr_id;
                return curr_id++;
            } else {
                return id->second;
            }
        }
    };
    for (it = name_to_seq.begin(); it != name_to_seq.end(); it++){
        auto source_id = get_add_id((it->second).name);
        //Make us some nodes
        Node n;
        n.set_sequence((it->second).sequence);
        n.set_id(source_id);
        n.set_name((it->second).name);
        graph->add_node(n);
        // Now some edges. Since they're placed in this map
        // by their from_node, it's no big deal to just iterate
        // over them.
        for (edge_elem l : seq_to_edges[(it->second).name]){
            auto sink_id = get_add_id(l.sink_name);
            Edge e;
            e.set_from(source_id);
            e.set_to(sink_id);
            e.set_from_start(!l.source_orientation_forward);
            e.set_to_end(!l.sink_orientation_forward);
            // get the cigar
            auto cigar_elems = vcflib::splitCigar(l.alignment);
            if (cigar_elems.size() == 1
                && cigar_elems.front().first > 0
                && cigar_elems.front().second == "M") {
                    reduce_overlaps = true;
                    e.set_overlap(cigar_elems.front().first);
            }
            graph->add_edge(e);
        }
        // for (path_elem p: seq_to_paths[(it->second).name]){
        //     paths.append_mapping(p.name, source_id, p.rank ,p.is_reverse);
        // }
        // remove overlapping sequences from the graph
    }
    map<string, path_elem> n_to_p = gg.get_name_to_path();
    for (auto name_path : n_to_p){
        for (int np = 0; np < name_path.second.segment_names.size(); np++){
            graph->paths.append_mapping(name_path.first, stol(name_path.second.segment_names[np]), np + 1, !name_path.second.orientations[np]);
        }
    }
    if (reduce_overlaps) {
        graph->bluntify();
    }
}

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
        out << p_elem.to_string() << endl;
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
