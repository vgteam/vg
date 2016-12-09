#include "phase_duplicator.hpp"

namespace vg {

using namespace std;

PhaseDuplicator::PhaseDuplicator(const xg::XG& index) : index(index) {
    // Nothing to do!
}

pair<Graph, vector<Translation>> PhaseDuplicator::duplicate(set<id_t> subgraph, id_t& next_id) const {

    // Allocate space for the result
    pair<Graph, vector<Translation>> to_return;
    // Name the pair members
    Graph& duplicated = to_return.first;
    vector<Translation>& translations = to_return.second;
    
    // TODO: actually do the duplication
    
    // Ship out the result
    return to_return;
}

}
