/** \file
 * Implements SnarlManagerDecomposition.
 */

#include "snarl_manager_decomposition.hpp"

#include <algorithm>
#include <stdexcept>

namespace vg {

using namespace std;

// A net handle packs the type in its low 3 bits, the start and end endpoints
// of the traversal in the next 2 bits each, and the payload in the rest.
static const uint64_t TYPE_BITS = 3;
static const uint64_t ENDPOINT_BITS = 2;
static const uint64_t PAYLOAD_SHIFT = TYPE_BITS + 2 * ENDPOINT_BITS;

SnarlManagerDecomposition::SnarlManagerDecomposition(const SnarlManager& manager, const HandleGraph& graph) :
    manager(manager), graph(graph), snarls_by_number(manager.num_snarls(), nullptr), hidden(manager.num_snarls(), false) {

    manager.for_each_snarl_unindexed([&](const Snarl* snarl) {
        size_t number = manager.snarl_number(snarl);
        snarls_by_number[number] = snarl;
        hidden[number] = manager.is_trivial(snarl, graph);
    });
}

net_handle_t SnarlManagerDecomposition::pack(NetType type, endpoint_t start, endpoint_t end, uint64_t payload) {
    uint64_t packed = (payload << PAYLOAD_SHIFT) | ((uint64_t) end << (TYPE_BITS + ENDPOINT_BITS)) |
        ((uint64_t) start << TYPE_BITS) | (uint64_t) type;
    return handlegraph::as_net_handle(packed);
}

SnarlManagerDecomposition::NetType SnarlManagerDecomposition::type_of(const net_handle_t& net) {
    return (NetType) (as_integer(net) & ((1 << TYPE_BITS) - 1));
}

uint64_t SnarlManagerDecomposition::payload_of(const net_handle_t& net) {
    return as_integer(net) >> PAYLOAD_SHIFT;
}

net_handle_t SnarlManagerDecomposition::with_traversal(const net_handle_t& net, endpoint_t start, endpoint_t end) {
    return pack(type_of(net), start, end, payload_of(net));
}

handle_t SnarlManagerDecomposition::to_handle(const Visit& visit) const {
    return graph.get_handle(visit.node_id(), visit.backward());
}

net_handle_t SnarlManagerDecomposition::node_net(const handle_t& handle) const {
    return get_net(handle, &graph);
}

net_handle_t SnarlManagerDecomposition::snarl_net(const Snarl* snarl, bool forward) const {
    return pack(NetType::SNARL, forward ? START : END, forward ? END : START, manager.snarl_number(snarl));
}

net_handle_t SnarlManagerDecomposition::chain_net(const Chain& chain, endpoint_t start, endpoint_t end) const {
    return pack(NetType::CHAIN, start, end, manager.snarl_number(chain.front().first));
}

const Chain& SnarlManagerDecomposition::chain_of_net(const net_handle_t& net) const {
    return *manager.chain_of(snarls_by_number.at(payload_of(net)));
}

bool SnarlManagerDecomposition::is_hidden(const Snarl* snarl) const {
    return hidden[manager.snarl_number(snarl)];
}

bool SnarlManagerDecomposition::is_cyclic(const Chain& chain) const {
    return chain_node(chain, 0) == chain_node(chain, chain.size());
}

size_t SnarlManagerDecomposition::node_count(const Chain& chain) const {
    return is_cyclic(chain) ? chain.size() : chain.size() + 1;
}

handle_t SnarlManagerDecomposition::chain_node(const Chain& chain, size_t index) const {
    if (index < chain.size()) {
        // The node we read into this snarl through, going along the chain.
        const pair<const Snarl*, bool>& entry = chain[index];
        return entry.second ? graph.flip(to_handle(entry.first->end())) : to_handle(entry.first->start());
    } else {
        // The node we read out of the last snarl through.
        const pair<const Snarl*, bool>& entry = chain.back();
        return entry.second ? graph.flip(to_handle(entry.first->start())) : to_handle(entry.first->end());
    }
}

vector<net_handle_t> SnarlManagerDecomposition::chain_children(const Chain& chain) const {
    vector<net_handle_t> children;
    for (size_t i = 0; i < chain.size(); i++) {
        children.push_back(node_net(chain_node(chain, i)));
        if (!is_hidden(chain[i].first)) {
            children.push_back(snarl_net(chain[i].first, !chain[i].second));
        }
    }
    if (!is_cyclic(chain)) {
        children.push_back(node_net(chain_node(chain, chain.size())));
    }
    return children;
}

bool SnarlManagerDecomposition::locate(const handle_t& handle, ChainPlace& place) const {
    for (bool flipped : {false, true}) {
        handle_t probe = flipped ? graph.flip(handle) : handle;
        const Snarl* into = manager.into_which_snarl(graph.get_id(probe), graph.get_is_reverse(probe));
        if (into == nullptr) {
            continue;
        }
        place.chain = manager.chain_of(into);
        size_t rank = manager.chain_rank_of(into);
        // The probe reads into the snarl either through the node before it
        // along the chain, or backward through the node after it.
        bool probe_forward = (probe == chain_node(*place.chain, rank));
        place.index = probe_forward ? rank : rank + 1;
        if (place.index == node_count(*place.chain)) {
            // In a cyclic chain the node after the last snarl is the first node.
            place.index = 0;
        }
        place.forward = (probe_forward != flipped);
        return true;
    }
    return false;
}

net_handle_t SnarlManagerDecomposition::parent_snarl_net(const Snarl* parent) const {
    return parent == nullptr ? get_root() : snarl_net(parent, true);
}

vector<nid_t> SnarlManagerDecomposition::trivial_chain_nodes(const Snarl* snarl) const {
    vector<nid_t> nodes;
    for (const id_t& id : manager.shallow_contents(snarl, graph, false).first) {
        if (manager.into_which_snarl(id, false) == nullptr && manager.into_which_snarl(id, true) == nullptr) {
            // This node bounds no child snarl, so it is a chain by itself.
            nodes.push_back(id);
        }
    }
    std::sort(nodes.begin(), nodes.end());
    return nodes;
}

void SnarlManagerDecomposition::index_trivial_chains() const {
    call_once(trivial_chains_indexed, [&]() {
        for (size_t number = 0; number < snarls_by_number.size(); number++) {
            if (hidden[number]) {
                continue;
            }
            for (const nid_t& child : trivial_chain_nodes(snarls_by_number[number])) {
                trivial_chain_parents[child] = snarls_by_number[number];
            }
        }
        graph.for_each_handle([&](const handle_t& handle) {
            nid_t child = graph.get_id(handle);
            if (!trivial_chain_parents.count(child) && manager.into_which_snarl(child, false) == nullptr &&
                manager.into_which_snarl(child, true) == nullptr) {
                // This node is in no snarl and bounds none.
                trivial_chain_parents[child] = nullptr;
                root_trivial_chain_nodes.push_back(child);
            }
        });
        std::sort(root_trivial_chain_nodes.begin(), root_trivial_chain_nodes.end());
    });
}

const Snarl* SnarlManagerDecomposition::trivial_chain_parent(nid_t id) const {
    index_trivial_chains();
    return trivial_chain_parents.at(id);
}

net_handle_t SnarlManagerDecomposition::get_snarl_net(const Snarl* snarl) const {
    return snarl_net(snarl, true);
}

const Snarl* SnarlManagerDecomposition::get_snarl(const net_handle_t& net) const {
    if (type_of(net) != NetType::SNARL && type_of(net) != NetType::SENTINEL) {
        throw runtime_error("error: net handle does not refer to a snarl");
    }
    return snarls_by_number.at(payload_of(net));
}

net_handle_t SnarlManagerDecomposition::get_root() const {
    return pack(NetType::ROOT, TIP, TIP, 0);
}

bool SnarlManagerDecomposition::is_root(const net_handle_t& net) const {
    return type_of(net) == NetType::ROOT;
}

bool SnarlManagerDecomposition::is_snarl(const net_handle_t& net) const {
    return type_of(net) == NetType::SNARL;
}

bool SnarlManagerDecomposition::is_chain(const net_handle_t& net) const {
    return type_of(net) == NetType::CHAIN || type_of(net) == NetType::TRIVIAL_CHAIN;
}

bool SnarlManagerDecomposition::is_node(const net_handle_t& net) const {
    return type_of(net) == NetType::NODE;
}

bool SnarlManagerDecomposition::is_sentinel(const net_handle_t& net) const {
    return type_of(net) == NetType::SENTINEL;
}

net_handle_t SnarlManagerDecomposition::get_net(const handle_t& handle, const HandleGraph* graph) const {
    bool reverse = graph->get_is_reverse(handle);
    return pack(NetType::NODE, reverse ? END : START, reverse ? START : END, graph->get_id(handle));
}

handle_t SnarlManagerDecomposition::get_handle(const net_handle_t& net, const HandleGraph* graph) const {
    switch (type_of(net)) {
    case NetType::NODE:
    case NetType::TRIVIAL_CHAIN:
        return graph->get_handle(payload_of(net), starts_at(net) == END);
    case NetType::SENTINEL:
        {
            // The first endpoint names the bound and the second which way it
            // faces: toward the other bound is in, and toward itself is out.
            const Snarl* snarl = snarls_by_number.at(payload_of(net));
            bool facing_in = starts_at(net) != ends_at(net);
            if (starts_at(net) == START) {
                handle_t in = graph->get_handle(snarl->start().node_id(), snarl->start().backward());
                return facing_in ? in : graph->flip(in);
            } else {
                handle_t out = graph->get_handle(snarl->end().node_id(), snarl->end().backward());
                return facing_in ? graph->flip(out) : out;
            }
        }
    default:
        throw runtime_error("error: trying to get a handle from a snarl, chain, or root");
    }
}

net_handle_t SnarlManagerDecomposition::get_parent(const net_handle_t& child) const {
    endpoint_t start = starts_at(child);
    endpoint_t end = ends_at(child);
    bool directed = (start == START && end == END) || (start == END && end == START);
    switch (type_of(child)) {
    case NetType::ROOT:
        throw runtime_error("error: trying to find the parent of the root");
    case NetType::SENTINEL:
        {
            // A sentinel facing in has the snarl going the same way; one
            // facing out has the snarl going out through that bound.
            endpoint_t snarl_start = start;
            if (start == end) {
                snarl_start = (start == START) ? END : START;
            }
            return pack(NetType::SNARL, snarl_start, snarl_start == START ? END : START, payload_of(child));
        }
    case NetType::SNARL:
        {
            const Snarl* snarl = snarls_by_number.at(payload_of(child));
            const Chain& chain = *manager.chain_of(snarl);
            if (!directed) {
                return chain_net(chain, START, END);
            }
            // The chain goes the way the snarl goes, after accounting for the
            // snarl's orientation in it.
            bool backward_in_chain = chain[manager.chain_rank_of(snarl)].second;
            bool chain_forward = (start == START) != backward_in_chain;
            return chain_net(chain, chain_forward ? START : END, chain_forward ? END : START);
        }
    case NetType::NODE:
        {
            handle_t handle = graph.get_handle(payload_of(child), start == END);
            ChainPlace place;
            if (!locate(handle, place)) {
                // The node is a trivial chain by itself, going the same way.
                return pack(NetType::TRIVIAL_CHAIN, start, end, payload_of(child));
            }
            return chain_net(*place.chain, place.forward ? START : END, place.forward ? END : START);
        }
    case NetType::CHAIN:
        return parent_snarl_net(manager.parent_of(chain_of_net(child).front().first));
    case NetType::TRIVIAL_CHAIN:
        return parent_snarl_net(trivial_chain_parent(payload_of(child)));
    }
    throw runtime_error("error: unknown net handle type");
}

net_handle_t SnarlManagerDecomposition::get_bound(const net_handle_t& snarl, bool get_end, bool face_in) const {
    switch (type_of(snarl)) {
    case NetType::SNARL:
        {
            // The first endpoint names the bound; the second is the other
            // bound if the sentinel faces in, and the same bound if it faces
            // out.
            endpoint_t bound = get_end ? END : START;
            endpoint_t other = get_end ? START : END;
            return pack(NetType::SENTINEL, bound, face_in ? other : bound, payload_of(snarl));
        }
    case NetType::CHAIN:
        {
            const Chain& chain = chain_of_net(snarl);
            handle_t node = get_end ? chain_node(chain, chain.size()) : chain_node(chain, 0);
            // The start faces in and the end faces out when read along the chain.
            return node_net(face_in == get_end ? graph.flip(node) : node);
        }
    case NetType::TRIVIAL_CHAIN:
        {
            handle_t node = graph.get_handle(payload_of(snarl), false);
            return node_net(face_in == get_end ? graph.flip(node) : node);
        }
    case NetType::ROOT:
        throw runtime_error("error: trying to get the bounds of the root");
    default:
        throw runtime_error("error: trying to get the bounds of a node or sentinel");
    }
}

net_handle_t SnarlManagerDecomposition::flip(const net_handle_t& net) const {
    if (type_of(net) == NetType::SENTINEL) {
        // Keep the bound and turn which way it faces.
        return with_traversal(net, starts_at(net), ends_at(net) == START ? END : START);
    }
    return with_traversal(net, ends_at(net), starts_at(net));
}

net_handle_t SnarlManagerDecomposition::canonical(const net_handle_t& net) const {
    switch (type_of(net)) {
    case NetType::ROOT:
        return get_root();
    case NetType::SENTINEL:
        // The bound facing in.
        return with_traversal(net, starts_at(net), starts_at(net) == START ? END : START);
    default:
        {
            // The first realizable traversal, which is start-to-end whenever
            // that is realizable.
            net_handle_t found = with_traversal(net, START, END);
            for_each_traversal_impl(net, [&](const net_handle_t& traversal) {
                found = traversal;
                return false;
            });
            return found;
        }
    }
}

handlegraph::SnarlDecomposition::endpoint_t SnarlManagerDecomposition::starts_at(const net_handle_t& traversal) const {
    return (endpoint_t) ((as_integer(traversal) >> TYPE_BITS) & ((1 << ENDPOINT_BITS) - 1));
}

handlegraph::SnarlDecomposition::endpoint_t SnarlManagerDecomposition::ends_at(const net_handle_t& traversal) const {
    return (endpoint_t) ((as_integer(traversal) >> (TYPE_BITS + ENDPOINT_BITS)) & ((1 << ENDPOINT_BITS) - 1));
}

net_handle_t SnarlManagerDecomposition::get_parent_traversal(const net_handle_t& traversal_start,
                                                             const net_handle_t& traversal_end) const {
    if (is_sentinel(traversal_start) && is_sentinel(traversal_end)) {
        return pack(NetType::SNARL, starts_at(traversal_start), starts_at(traversal_end), payload_of(traversal_start));
    }
    if (!is_node(traversal_start) || !is_node(traversal_end)) {
        throw runtime_error("error: the snarl manager does not record internal tips");
    }
    handle_t first = get_handle(traversal_start, &graph);
    handle_t last = get_handle(traversal_end, &graph);
    ChainPlace start_place;
    ChainPlace end_place;
    if (!locate(first, start_place)) {
        // A trivial chain is traversed the way its node is read.
        return pack(NetType::TRIVIAL_CHAIN, starts_at(traversal_start), ends_at(traversal_start), graph.get_id(first));
    }
    if (!locate(last, end_place) || end_place.chain != start_place.chain) {
        throw runtime_error("error: traversal bounds are not in the same chain");
    }
    // Reading forward from the first node starts at the chain's start, and
    // reading backward from the last node starts at its end. A cyclic
    // chain's first node is also its last.
    const Chain& chain = *start_place.chain;
    size_t last_index = is_cyclic(chain) ? 0 : node_count(chain) - 1;
    if (start_place.index != (start_place.forward ? 0 : last_index) ||
        end_place.index != (end_place.forward ? last_index : 0)) {
        throw runtime_error("error: traversal bounds are not the bounds of their chain");
    }
    return chain_net(chain, start_place.forward ? START : END, end_place.forward ? END : START);
}

bool SnarlManagerDecomposition::for_each_child_impl(const net_handle_t& traversal,
                                                    const function<bool(const net_handle_t&)>& iteratee) const {
    switch (type_of(traversal)) {
    case NetType::ROOT:
        {
            for (const Chain& chain : manager.chains_of(nullptr)) {
                if (!iteratee(chain_net(chain, START, END))) {
                    return false;
                }
            }
            index_trivial_chains();
            for (const nid_t& id : root_trivial_chain_nodes) {
                if (!iteratee(pack(NetType::TRIVIAL_CHAIN, START, END, id))) {
                    return false;
                }
            }
            return true;
        }
    case NetType::SNARL:
        {
            const Snarl* snarl = snarls_by_number.at(payload_of(traversal));
            for (const Chain& chain : manager.chains_of(snarl)) {
                if (!iteratee(chain_net(chain, START, END))) {
                    return false;
                }
            }
            for (const nid_t& id : trivial_chain_nodes(snarl)) {
                if (!iteratee(pack(NetType::TRIVIAL_CHAIN, START, END, id))) {
                    return false;
                }
            }
            return true;
        }
    case NetType::CHAIN:
        {
            const Chain& chain = chain_of_net(traversal);
            vector<net_handle_t> children = chain_children(chain);
            if (starts_at(traversal) == END && ends_at(traversal) == START) {
                // Read the chain backward, still starting at its start node
                // if it is cyclic.
                std::reverse(children.begin(), children.end());
                if (is_cyclic(chain)) {
                    std::rotate(children.begin(), children.end() - 1, children.end());
                }
                for (net_handle_t& child : children) {
                    child = flip(child);
                }
            } else if (!(starts_at(traversal) == START && ends_at(traversal) == END)) {
                // Not a directed traversal, so produce everything start-to-end.
                for (net_handle_t& child : children) {
                    child = with_traversal(child, START, END);
                }
            }
            for (const net_handle_t& child : children) {
                if (!iteratee(child)) {
                    return false;
                }
            }
            return true;
        }
    case NetType::TRIVIAL_CHAIN:
        {
            bool reverse = (starts_at(traversal) == END && ends_at(traversal) == START);
            return iteratee(pack(NetType::NODE, reverse ? END : START, reverse ? START : END, payload_of(traversal)));
        }
    default:
        // Nodes and sentinels have no children.
        return true;
    }
}

bool SnarlManagerDecomposition::for_each_traversal_impl(const net_handle_t& item,
                                                        const function<bool(const net_handle_t&)>& iteratee) const {
    vector<pair<endpoint_t, endpoint_t>> kinds;
    switch (type_of(item)) {
    case NetType::ROOT:
        kinds.emplace_back(TIP, TIP);
        break;
    case NetType::SENTINEL:
        kinds.emplace_back(starts_at(item), ends_at(item));
        break;
    case NetType::NODE:
    case NetType::TRIVIAL_CHAIN:
        kinds.emplace_back(START, END);
        kinds.emplace_back(END, START);
        break;
    case NetType::SNARL:
        {
            const Snarl* snarl = snarls_by_number.at(payload_of(item));
            if (snarl->start_end_reachable()) {
                kinds.emplace_back(START, END);
                kinds.emplace_back(END, START);
            }
            if (snarl->start_self_reachable()) {
                kinds.emplace_back(START, START);
            }
            if (snarl->end_self_reachable()) {
                kinds.emplace_back(END, END);
            }
        }
        break;
    case NetType::CHAIN:
        {
            // A chain can be crossed if all its snarls can be. It can be
            // entered and left through its start if some snarl can turn
            // around on the side facing the start, and all snarls before it
            // can be crossed; likewise for its end.
            const Chain& chain = chain_of_net(item);
            auto turns_toward = [&](const pair<const Snarl*, bool>& entry, bool toward_start) {
                // Can this snarl turn around on its side facing the chain's start or end?
                bool on_snarl_start = (toward_start != entry.second);
                return on_snarl_start ? entry.first->start_self_reachable() : entry.first->end_self_reachable();
            };
            bool through = true;
            bool start_start = false;
            for (const pair<const Snarl*, bool>& entry : chain) {
                if (through && turns_toward(entry, true)) {
                    start_start = true;
                }
                through = through && entry.first->start_end_reachable();
            }
            bool end_end = false;
            bool through_from_end = true;
            for (auto it = chain.rbegin(); it != chain.rend(); ++it) {
                if (through_from_end && turns_toward(*it, false)) {
                    end_end = true;
                }
                through_from_end = through_from_end && it->first->start_end_reachable();
            }
            if (through) {
                kinds.emplace_back(START, END);
                kinds.emplace_back(END, START);
            }
            if (start_start) {
                kinds.emplace_back(START, START);
            }
            if (end_end) {
                kinds.emplace_back(END, END);
            }
        }
        break;
    }
    for (const pair<endpoint_t, endpoint_t>& kind : kinds) {
        if (!iteratee(with_traversal(item, kind.first, kind.second))) {
            return false;
        }
    }
    return true;
}

bool SnarlManagerDecomposition::follow_net_edges_impl(const net_handle_t& here, const HandleGraph* graph, bool go_left,
                                                      const function<bool(const net_handle_t&)>& iteratee) const {
    // Going left from here is going right from its flip, and what we reach
    // is read the way we travel.
    return follow_right(go_left ? flip(here) : here, iteratee);
}

bool SnarlManagerDecomposition::follow_right(const net_handle_t& here,
                                             const function<bool(const net_handle_t&)>& iteratee) const {
    NetType type = type_of(here);
    if (type == NetType::ROOT) {
        return true;
    }

    if (type == NetType::NODE || type == NetType::SNARL) {
        // Step to the neighbouring child in the parent chain, if there is one.
        if (type == NetType::SNARL) {
            if (ends_at(here) == TIP) {
                return true;
            }
            const Snarl* snarl = snarls_by_number.at(payload_of(here));
            const Chain& chain = *manager.chain_of(snarl);
            size_t rank = manager.chain_rank_of(snarl);
            net_handle_t forward_snarl = snarl_net(snarl, !chain[rank].second);
            // Leaving through the side that faces the chain's end goes along
            // the chain, to the node after the snarl.
            bool along_chain = (ends_at(here) == ends_at(forward_snarl));
            handle_t node = chain_node(chain, along_chain ? rank + 1 : rank);
            return iteratee(node_net(along_chain ? node : graph.flip(node)));
        }
        ChainPlace place;
        if (!locate(get_handle(here, &graph), place)) {
            // A trivial chain has nothing beside its node.
            return true;
        }
        const Chain& chain = *place.chain;
        // Find the snarl beside the node in the direction we go.
        size_t rank;
        if (place.forward) {
            if (place.index == chain.size()) {
                return true;
            }
            rank = place.index;
        } else {
            if (place.index == 0 && !is_cyclic(chain)) {
                return true;
            }
            rank = (place.index == 0 ? chain.size() : place.index) - 1;
        }
        const pair<const Snarl*, bool>& entry = chain[rank];
        if (is_hidden(entry.first)) {
            // Step over a hidden snarl to the node on its other side.
            handle_t node = chain_node(chain, place.forward ? rank + 1 : rank);
            return iteratee(node_net(place.forward ? node : graph.flip(node)));
        }
        // Cross the snarl the way we travel.
        net_handle_t forward_snarl = snarl_net(entry.first, !entry.second);
        return iteratee(place.forward ? forward_snarl : flip(forward_snarl));
    }
    
    // Otherwise we are a child of a snarl or the root, or a bound of a snarl,
    // and we follow graph edges to other children of the same parent.
    const Snarl* parent;
    handle_t exit;
    switch (type) {
    case NetType::SENTINEL:
        if (starts_at(here) == ends_at(here)) {
            // Facing out of the snarl.
            return true;
        }
        parent = snarls_by_number.at(payload_of(here));
        exit = get_handle(here, &graph);
        break;
    case NetType::CHAIN:
        {
            const Chain& chain = chain_of_net(here);
            parent = manager.parent_of(chain.front().first);
            if (ends_at(here) == TIP) {
                return true;
            }
            exit = ends_at(here) == END ? chain_node(chain, chain.size()) : graph.flip(chain_node(chain, 0));
        }
        break;
    default:
        // A trivial chain
        parent = trivial_chain_parent(payload_of(here));
        exit = graph.get_handle(payload_of(here), ends_at(here) == START);
        break;
    }

    return graph.follow_edges(exit, false, [&](const handle_t& next) {
        if (parent != nullptr) {
            // We may reach the parent's bounds, reading out of the snarl.
            if (next == graph.flip(to_handle(parent->start()))) {
                return iteratee(pack(NetType::SENTINEL, START, START, manager.snarl_number(parent)));
            }
            if (next == to_handle(parent->end())) {
                return iteratee(pack(NetType::SENTINEL, END, END, manager.snarl_number(parent)));
            }
            if (graph.get_id(next) == parent->start().node_id() || graph.get_id(next) == parent->end().node_id()) {
                // This edge leaves the snarl.
                return true;
            }
        }
        ChainPlace place;
        if (!locate(next, place)) {
            // A trivial chain, entered the way the node is read.
            net_handle_t node = node_net(next);
            return iteratee(pack(NetType::TRIVIAL_CHAIN, starts_at(node), ends_at(node), graph.get_id(next)));
        }
        // Cross the chain the way we travel.
        if (place.index == 0 && place.forward) {
            return iteratee(chain_net(*place.chain, START, END));
        }
        if (!place.forward && (place.index == node_count(*place.chain) - 1 || is_cyclic(*place.chain))) {
            return iteratee(chain_net(*place.chain, END, START));
        }
        // An edge into the middle of a chain does not come from a sibling.
        return true;
    });
}

bool SnarlManagerDecomposition::for_each_tippy_child_impl(const net_handle_t& parent,
                                                          const function<bool(const net_handle_t&)>& iteratee) const {
    throw runtime_error("error: the snarl manager does not record internal tips");
}

}
