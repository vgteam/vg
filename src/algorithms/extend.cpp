#include "extend.hpp"

namespace vg {
namespace algorithms {

using namespace std;

    void extend(const HandleGraph* source, MutableHandleGraph* into) {
        
        // add any nodes that are not yet present
        source->for_each_handle([&](const handle_t& handle) {
            if (!into->has_node(source->get_id(handle))) {
                into->create_handle(source->get_sequence(handle), source->get_id(handle));
            }
            // always keep going
            return true;
        });
        
        // add any edges that are not yet present
        source->for_each_edge([&](const edge_t& edge) {
            edge_t into_edge(into->get_handle(source->get_id(edge.first), source->get_is_reverse(edge.first)),
                             into->get_handle(source->get_id(edge.second), source->get_is_reverse(edge.second)));
            if (!into->has_edge(into_edge)) {
                into->create_edge(into_edge);
            }
            // always keep going
            return true;
        });
    }
}
}
