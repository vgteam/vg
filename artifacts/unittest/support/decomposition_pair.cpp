/**
 * \file decomposition_pair.cpp
 * Implements DecompositionPair.
 */

#include "decomposition_pair.hpp"

namespace vg {
namespace unittest {

DecompositionPair::DecompositionPair(const HandleGraph& graph, const HandleGraphSnarlFinder& finder) :
    graph(graph),
    events(capture_events(finder, graph)),
    manager(ReplaySnarlFinder(&graph, events).find_snarls()),
    adapter(manager, graph) {

    ReplaySnarlFinder replay(&graph, events);
    fill_in_distance_index(&distance_index, &graph, &replay);
}

bool DecompositionPair::has_known_difference() const {
    bool found = false;
    manager.for_each_snarl_preorder([&](const Snarl* snarl) {
        found = found || snarl->start().node_id() == snarl->end().node_id() ||
            (snarl->type() != ULTRABUBBLE && manager.is_leaf(snarl) &&
             manager.shallow_contents(snarl, graph, false).first.empty());
    });
    manager.for_each_chain([&](const Chain* chain) {
        found = found || get_start_of(*chain) == get_end_of(*chain);
    });
    return found;
}

}
}
