#include "remove_high_degree.hpp"

namespace vg {
namespace algorithms {

using namespace std;

void remove_high_degree_nodes(DeletableHandleGraph& g, int max_degree) {
    vector<handle_t> to_remove;
    g.for_each_handle([&](const handle_t& h) {
            int edge_count = 0;
            g.follow_edges(h, false, [&](const handle_t& ignored) {
                    ++edge_count;
                });
            g.follow_edges(h, true, [&](const handle_t& ignored) {
                    ++edge_count;
                });
            if (edge_count > max_degree) {
                to_remove.push_back(h);
            }
        });
    // now destroy the high degree nodes
    for (auto& h : to_remove) {
        g.destroy_handle(h);
    }
}

}
}
