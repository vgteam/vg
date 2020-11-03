#include "next_pos_chars.hpp"


namespace vg {
namespace algorithms {

map<pos_t, char> next_pos_chars(const PathPositionHandleGraph& graph, pos_t pos) {
    map<pos_t, char> nexts;
    handle_t handle = graph.get_handle(id(pos), is_rev(pos));
    if (offset(pos) < graph.get_length(handle)-1) {
        ++get_offset(pos);
        char c = graph.get_base(handle, offset(pos));
        nexts[pos] = c;
    } else {
        graph.follow_edges(handle, false, [&](const handle_t& next) {
                char c = graph.get_base(next, 0);
                nexts[make_pos_t(graph.get_id(next), graph.get_is_reverse(next), 0)] = c;
            });
    }
    return nexts;
}

}
}
