#include "append_graph.hpp"

#include "topological_sort.hpp"
#include "copy_graph.hpp"

using namespace std;

namespace vg {
namespace algorithms {

void append_handle_graph(const HandleGraph* from, MutableHandleGraph* into) {

    // get the heads and tails before copying
    vector<handle_t> heads = algorithms::head_nodes(from);
    vector<handle_t> tails = algorithms::tail_nodes(into);

    // copy over the nodes and edges
    // beware: no checking duplicates
    copy_handle_graph(from, into);

    // connect all the tips
    for (handle_t head : heads) {
        handle_t into_head = into->get_handle(from->get_id(head), from->get_is_reverse(head));
        for (handle_t tail : tails) {
            into->create_edge(tail, into_head);
        }
    }
}
    
void append_path_handle_graph(const PathHandleGraph* from, MutablePathMutableHandleGraph* into,
                              bool only_connect_path_tips) {

    // get the heads and tails before copying
    vector<handle_t> heads;
    vector<handle_t> tails;
    if (!only_connect_path_tips) {
        heads = algorithms::head_nodes(from);
        tails = algorithms::tail_nodes(into);
    }

    // copy over the nodes and edges
    // beware: no checking duplicates
    copy_handle_graph(from, into);
    
    // connect all the tips
    for (handle_t head : heads) {
        handle_t into_head = into->get_handle(from->get_id(head), from->get_is_reverse(head));
        for (handle_t tail : tails) {
            into->create_edge(tail, into_head);
        }
    }

    // append the paths, adding linking edges
    from->for_each_path_handle([&](const path_handle_t& path_handle) {
            string path_name = from->get_path_name(path_handle);
            if (!into->has_path(path_name)) {
                // just copy it over if it's not there
                copy_path(from, path_handle, into);
            } else {
                // else we append it
                path_handle_t into_path_handle = into->get_path_handle(path_name);
                
                // add the missing edge
                handle_t tail = into->get_handle_of_step(into->path_back(into_path_handle));
                handle_t head = from->get_handle_of_step(from->path_begin(path_handle));
                if (!into->has_edge(tail, head)) {
                    into->create_edge(tail, head);
                }
                
                // append the steps
                from->for_each_step_in_path(path_handle, [&](const step_handle_t& step_handle) {
                        handle_t handle = from->get_handle_of_step(step_handle);
                        handle_t into_handle = into->get_handle(from->get_id(handle), from->get_is_reverse(handle));
                        into->append_step(into_path_handle, into_handle);
                    });
            }                    
        });
    
}
}
}
