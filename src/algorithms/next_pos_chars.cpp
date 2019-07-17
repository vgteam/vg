#include "next_pos_chars.hpp"


namespace vg {
namespace algorithms {

map<pos_t, char> next_pos_chars(const PathPositionHandleGraph& graph, pos_t pos) {
    map<pos_t, char> nexts;
    handle_t handle = xgidx->get_handle(id(pos), is_rev(pos));
    if (offset(pos) < xgidx->get_length(handle)-1) {
        ++get_offset(pos);
        char c = xgidx->get_base(handle, offset(pos));
        nexts[pos] = c;
    } else {
        xgidx->follow_edges(handle, false, [&](const handle_t& next) {
                char c = xgidx->get_base(next, 0);
                nexts[make_pos_t(xgidx->get_id(next), xgidx->get_is_reverse(next), 0)] = c;
            });
    }
    return nexts;
}

}
}
