#include "disjoint_components.hpp"

namespace vg {
namespace algorithms {

list<bdsg::HashGraph> disjoint_components(const HandleGraph& graph) {
    
    vector<unordered_set<nid_t>> weak_comps = handlealgs::weakly_connected_components(&graph);
    
    list<bdsg::HashGraph> comps;
    for (const auto& weak_comp : weak_comps) {
        comps.emplace_back();
        auto& comp = comps.back();
        for (auto node_id : weak_comp) {
            comp.create_handle(graph.get_sequence(graph.get_handle(node_id)), node_id);
        }
        comp.for_each_handle([&](const handle_t& handle) {
            handle_t original_handle = graph.get_handle(comp.get_id(handle));
            // TODO: this will create duplicate edges if we ever decide
            // to switch back to not deduplicating on the fly
            graph.follow_edges(handle, true, [&](const handle_t& prev) {
                comp.create_edge(comp.get_handle(graph.get_id(prev), graph.get_is_reverse(prev)),
                                 handle);
            });
            graph.follow_edges(handle, false, [&](const handle_t& next) {
                comp.create_edge(handle,
                                 comp.get_handle(graph.get_id(next), graph.get_is_reverse(next)));
            });
        });
    }
    
    return comps;
}

}
}
