#include "to_gfa.hpp"

namespace vg {
namespace algorithms {

using namespace std;

void to_gfa(const PathHandleGraph& graph, ostream& out) {
    graph.for_each_handle([&out,&graph](const handle_t& h) {
            out << "S\t" << graph.get_id(h) << "\t" << graph.get_sequence(h) << std::endl;
            graph.follow_edges(h, true, [&out,&graph,&h](const handle_t& n) {
                    out << "L\t" << graph.get_id(h) << "\t"
                        << (graph.get_is_reverse(h)?"-":"+")
                        << "\t" << graph.get_id(n) << "\t"
                        << (graph.get_is_reverse(n)?"-":"+")
                        << "\t0M" << std::endl;
                });
        });
    graph.for_each_path_handle([&out,&graph](const path_handle_t& p) {
            //step_handle_t step = path_begin(p);
            out << "P\t" << graph.get_path_name(p) << "\t";
            graph.for_each_step_in_path(p, [&out,&graph](const step_handle_t& step) {
                    handle_t h = graph.get_handle_of_step(step);
                    out << graph.get_id(h) << (graph.get_is_reverse(h)?"-":"+");
                    if (graph.has_next_step(step)) out << ",";
                });
            out << "\t";
            graph.for_each_step_in_path(p, [&out,&graph](const step_handle_t& step) {
                    out << graph.get_length(graph.get_handle_of_step(step)) << "M";
                    if (graph.has_next_step(step)) out << ",";
                });
            out << std::endl;
        });
}
    
}
}
