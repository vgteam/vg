#include "id_sort.hpp"

#include "apply_bulk_modifications.hpp"

#include <algorithm>
#include <vector>

namespace vg {
namespace algorithms {

using namespace std;

vector<handle_t> id_order(HandleGraph* g) {
    // We will fill and sort this
    vector<handle_t> to_return;
    g->for_each_handle([&](const handle_t& handle) {
        // Collect all the handles
        to_return.push_back(handle);
    });
    
    std::sort(to_return.begin(), to_return.end(), [&](const handle_t& a, const handle_t& b) {
        // Sort in ID order with the standard algorithm
        return g->get_id(a) < g->get_id(b);
    });
    
    return to_return;
}

void id_sort(MutableHandleGraph* g) {
    if (g->node_size() <= 1) {
        // A graph with <2 nodes has only one sort.
        return;
    }
    
    // Order handles by ID and then apply the order to the graph
    apply_ordering(g, id_order(g));
}
    
}
}
