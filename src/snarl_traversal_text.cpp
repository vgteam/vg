#include "snarl_traversal_text.hpp"

namespace vg {

handle_t visit_handle(const HandleGraph& graph, const Visit& visit) {
    return graph.get_handle(visit.node_id(), visit.backward());
}

string snarl_name(const Snarl& snarl, const NodeTranslation* translation, bool in_brackets) {
    return site_name(snarl.start().node_id(), snarl.start().backward(), snarl.end().node_id(),
                     snarl.end().backward(), translation, in_brackets);
}

string trav_string(const HandleGraph& graph, const SnarlTraversal& trav,
                   const NodeTranslation* translation) {
    string seq;
    for (int i = 0; i < trav.visit_size(); ++i) {
        const Visit& visit = trav.visit(i);
        if (visit.node_id() > 0) {
            seq += graph.get_sequence(graph.get_handle(visit.node_id(), visit.backward()));
        } else {
            seq += snarl_name(visit.snarl(), translation, true);
        }
    }
    return seq;
}

vector<SnarlTraversal> traversals_of(const vector<SnarlTraversal>& travs,
                                     const vector<int>& numbers) {
    vector<SnarlTraversal> picked;
    picked.reserve(numbers.size());
    for (int number : numbers) {
        picked.push_back(number >= 0 ? travs[number] : SnarlTraversal());
    }
    return picked;
}

}
