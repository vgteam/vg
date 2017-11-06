#include "genome_state.hpp"

namespace vg {

using namespace std;

GenomeState::GenomeState(const SnarlManager& manager, const HandleGraph* graph,
    const unordered_set<pair<const Snarl*, const Snarl*>> telomeres) : telomeres(telomeres), manager(manager) {

    manager.for_each_snarl_preorder([&](const Snarl* snarl) {
        // For each snarl
        
        // Make an empty state for it
        state.emplace(snarl, SnarlState());
        
        // Make a net graph for it. TODO: we're not considering internal
        // connectivity, but what we really should do is consider internal
        // connectivity but only allowing for start to end traversals (but
        // including in unary snarls)
        net_graphs.emplace(snarl, manager.net_graph_of(snarl, graph, false));
    });
}


}
