/** \file
 * Implements SnarlManagerDecomposition.
 */

#include "snarl_manager_decomposition.hpp"

#include <algorithm>
#include <cassert>
#include <iterator>
#include <stdexcept>

namespace vg {

using namespace std;

namespace {

using endpoint_t = handlegraph::SnarlDecomposition::endpoint_t;
using enum handlegraph::SnarlDecomposition::endpoint_t;

/**
 * What a net handle refers to.
 *
 * The payload is the SnarlManager::snarl_number() of a snarl or a sentinel's
 * snarl, the snarl number of a chain's first snarl, and the node ID of a node
 * or a trivial chain.
 *
 * A sentinel's start endpoint names its bound. Its end endpoint is the other
 * bound if the sentinel faces into the snarl, and the same bound if it faces
 * out.
 */
enum class NetType : uint64_t {
    ROOT = 0,
    SNARL,
    CHAIN,
    TRIVIAL_CHAIN,
    NODE,
    SENTINEL
};

// A net handle packs the type in its low bits, then the start and end
// endpoints of the traversal, then the payload.
constexpr uint64_t TYPE_BITS = 3;
constexpr uint64_t TYPE_MASK = (uint64_t{1} << TYPE_BITS) - 1;
constexpr uint64_t ENDPOINT_BITS = 2;
constexpr uint64_t ENDPOINT_MASK = (uint64_t{1} << ENDPOINT_BITS) - 1;
constexpr uint64_t START_SHIFT = TYPE_BITS;
constexpr uint64_t END_SHIFT = TYPE_BITS + ENDPOINT_BITS;
constexpr uint64_t PAYLOAD_SHIFT = TYPE_BITS + 2 * ENDPOINT_BITS;
static_assert(static_cast<uint64_t>(NetType::SENTINEL) <= TYPE_MASK);
static_assert(static_cast<uint64_t>(TIP) <= ENDPOINT_MASK);

net_handle_t pack(NetType type, endpoint_t start, endpoint_t end, uint64_t payload) {
    assert(payload >> (64 - PAYLOAD_SHIFT) == 0);
    return handlegraph::as_net_handle((payload << PAYLOAD_SHIFT) |
                                      (static_cast<uint64_t>(end) << END_SHIFT) |
                                      (static_cast<uint64_t>(start) << START_SHIFT) |
                                      static_cast<uint64_t>(type));
}

NetType type_of(const net_handle_t& net) {
    return static_cast<NetType>(handlegraph::as_integer(net) & TYPE_MASK);
}

endpoint_t start_of(const net_handle_t& net) {
    return static_cast<endpoint_t>((handlegraph::as_integer(net) >> START_SHIFT) & ENDPOINT_MASK);
}

endpoint_t end_of(const net_handle_t& net) {
    return static_cast<endpoint_t>((handlegraph::as_integer(net) >> END_SHIFT) & ENDPOINT_MASK);
}

uint64_t payload_of(const net_handle_t& net) {
    return handlegraph::as_integer(net) >> PAYLOAD_SHIFT;
}

/// Make a net handle to the same thing with a different traversal.
net_handle_t with_traversal(const net_handle_t& net, endpoint_t start, endpoint_t end) {
    return pack(type_of(net), start, end, payload_of(net));
}

/// Swap START and END. Must not be given TIP.
endpoint_t opposite(endpoint_t endpoint) {
    return endpoint == START ? END : START;
}

/// Make a net handle to the same thing, traversed start-to-end if forward
/// and end-to-start otherwise.
net_handle_t directed(const net_handle_t& net, bool forward) {
    endpoint_t start = forward ? START : END;
    return with_traversal(net, start, opposite(start));
}

bool is_forward(const net_handle_t& net) {
    return start_of(net) == START && end_of(net) == END;
}

bool is_backward(const net_handle_t& net) {
    return start_of(net) == END && end_of(net) == START;
}

bool faces_in(const net_handle_t& sentinel) {
    return start_of(sentinel) != end_of(sentinel);
}

/// Get a node's net handle, read the way the handle reads.
net_handle_t node_net(const HandleGraph& graph, const handle_t& handle) {
    return directed(pack(NetType::NODE, START, END, graph.get_id(handle)), !graph.get_is_reverse(handle));
}

handle_t to_handle(const HandleGraph& graph, const Visit& visit) {
    return graph.get_handle(visit.node_id(), visit.backward());
}

/// The kinds of traversal a snarl or chain allows between its bounds.
struct Reachability {
    bool start_to_end;
    bool start_to_start;
    bool end_to_end;
};

Reachability reachability(const Snarl& snarl) {
    return {snarl.start_end_reachable(), snarl.start_self_reachable(), snarl.end_self_reachable()};
}

/**
 * Return true if a traversal that enters a chain through one endpoint can
 * leave through the same endpoint: it can cross each snarl it meets until
 * one snarl lets it turn around on the side it entered that snarl from.
 * Takes the chain's (snarl, backward) pairs in the order the traversal meets
 * them.
 */
template<std::input_iterator Iterator>
bool can_turn_back(Iterator begin, Iterator end, bool from_chain_start) {
    for (Iterator it = begin; it != end; ++it) {
        const auto& [snarl, backward] = *it;
        // The traversal enters the snarl's start if it and the snarl both go
        // along the chain or both go against it.
        if (from_chain_start != backward ? snarl->start_self_reachable() : snarl->end_self_reachable()) {
            return true;
        }
        if (!snarl->start_end_reachable()) {
            return false;
        }
    }
    return false;
}

Reachability reachability(const Chain& chain) {
    return {
        std::ranges::all_of(chain, [](const pair<const Snarl*, bool>& entry) {
            return entry.first->start_end_reachable();
        }),
        can_turn_back(chain.begin(), chain.end(), true),
        can_turn_back(chain.rbegin(), chain.rend(), false)
    };
}

/// List the start and end of each traversal kind that a reachability allows.
vector<pair<endpoint_t, endpoint_t>> traversal_kinds(const Reachability& reachable) {
    vector<pair<endpoint_t, endpoint_t>> kinds;
    if (reachable.start_to_end) {
        kinds.emplace_back(START, END);
        kinds.emplace_back(END, START);
    }
    if (reachable.start_to_start) {
        kinds.emplace_back(START, START);
    }
    if (reachable.end_to_end) {
        kinds.emplace_back(END, END);
    }
    return kinds;
}

}

////////////////////////////////////////////////////////////////////////////
// ChainLayout
////////////////////////////////////////////////////////////////////////////

SnarlManagerDecomposition::ChainLayout::ChainLayout(const SnarlManager& manager, const HandleGraph& graph) :
    manager(manager), graph(graph) {
}

handle_t SnarlManagerDecomposition::ChainLayout::node(const Chain& chain, size_t index) const {
    if (index == chain.size()) {
        return to_handle(graph, get_end_of(chain));
    }
    // The node read into the snarl with this rank, going along the chain.
    const auto& [snarl, backward] = chain[index];
    return to_handle(graph, backward ? reverse(snarl->end()) : snarl->start());
}

bool SnarlManagerDecomposition::ChainLayout::is_cyclic(const Chain& chain) const {
    return get_start_of(chain) == get_end_of(chain);
}

size_t SnarlManagerDecomposition::ChainLayout::node_count(const Chain& chain) const {
    return is_cyclic(chain) ? chain.size() : chain.size() + 1;
}

size_t SnarlManagerDecomposition::ChainLayout::last_index(const Chain& chain) const {
    return is_cyclic(chain) ? 0 : chain.size();
}

optional<SnarlManagerDecomposition::ChainLayout::Place>
SnarlManagerDecomposition::ChainLayout::locate(const handle_t& handle) const {
    for (bool flipped : {false, true}) {
        handle_t probe = flipped ? graph.flip(handle) : handle;
        const Snarl* into = manager.into_which_snarl(graph.get_id(probe), graph.get_is_reverse(probe));
        if (into == nullptr) {
            continue;
        }
        const Chain* chain = manager.chain_of(into);
        size_t rank = manager.chain_rank_of(into);
        // The probe reads into the snarl either through the node before it
        // along the chain, or backward through the node after it.
        bool probe_along = (probe == node(*chain, rank));
        size_t index = probe_along ? rank : rank + 1;
        if (index == node_count(*chain)) {
            // In a cyclic chain the node after the last snarl is the first node.
            index = 0;
        }
        return Place {chain, index, probe_along != flipped};
    }
    return nullopt;
}

optional<size_t> SnarlManagerDecomposition::ChainLayout::snarl_ahead(const Place& place) const {
    const Chain& chain = *place.chain;
    if (place.along) {
        return place.index < chain.size() ? optional<size_t>(place.index) : nullopt;
    }
    if (place.index > 0) {
        return place.index - 1;
    }
    // Reading backward out of the first node wraps around a cyclic chain.
    return is_cyclic(chain) ? optional<size_t>(chain.size() - 1) : nullopt;
}

optional<endpoint_t> SnarlManagerDecomposition::ChainLayout::entry_endpoint(const Place& place) const {
    if (place.along && place.index == 0) {
        return START;
    }
    if (!place.along && place.index == last_index(*place.chain)) {
        return END;
    }
    return nullopt;
}

optional<endpoint_t> SnarlManagerDecomposition::ChainLayout::exit_endpoint(const Place& place) const {
    if (place.along && place.index == last_index(*place.chain)) {
        return END;
    }
    if (!place.along && place.index == 0) {
        return START;
    }
    return nullopt;
}

////////////////////////////////////////////////////////////////////////////
// TrivialChains
////////////////////////////////////////////////////////////////////////////

SnarlManagerDecomposition::TrivialChains::TrivialChains(const SnarlManager& manager, const HandleGraph& graph) :
    manager(manager), graph(graph) {
}

bool SnarlManagerDecomposition::TrivialChains::bounds_no_snarl(nid_t id) const {
    return manager.into_which_snarl(id, false) == nullptr && manager.into_which_snarl(id, true) == nullptr;
}

vector<nid_t> SnarlManagerDecomposition::TrivialChains::nodes_in(const Snarl* parent) const {
    if (parent == nullptr) {
        index_parents();
        return root_nodes;
    }
    vector<nid_t> nodes;
    for (const id_t& id : manager.shallow_contents(parent, graph, false).first) {
        if (bounds_no_snarl(id)) {
            nodes.push_back(id);
        }
    }
    std::ranges::sort(nodes);
    return nodes;
}

void SnarlManagerDecomposition::TrivialChains::index_parents() const {
    call_once(parents_indexed, [this]() {
        manager.for_each_snarl_unindexed([this](const Snarl* snarl) {
            for (const nid_t& id : nodes_in(snarl)) {
                parents[id] = snarl;
            }
        });
        graph.for_each_handle([this](const handle_t& handle) {
            nid_t id = graph.get_id(handle);
            if (!parents.contains(id) && bounds_no_snarl(id)) {
                // This node is in no snarl and bounds none.
                parents[id] = nullptr;
                root_nodes.push_back(id);
            }
        });
        std::ranges::sort(root_nodes);
    });
}

const Snarl* SnarlManagerDecomposition::TrivialChains::parent_of(nid_t id) const {
    index_parents();
    return parents.at(id);
}

////////////////////////////////////////////////////////////////////////////
// SnarlManagerDecomposition
////////////////////////////////////////////////////////////////////////////

SnarlManagerDecomposition::SnarlManagerDecomposition(const SnarlManager& manager, const HandleGraph& graph) :
    manager(manager), graph(graph), layout(manager, graph), trivial_chains(manager, graph),
    snarls_by_number(manager.num_snarls(), nullptr), hidden(manager.num_snarls(), false) {

    manager.for_each_snarl_unindexed([this](const Snarl* snarl) {
        if (snarl->start().node_id() == snarl->end().node_id() &&
            snarl->start().backward() != snarl->end().backward()) {
            throw invalid_argument("error: SnarlManagerDecomposition cannot present unary snarls");
        }
        size_t number = this->manager.snarl_number(snarl);
        snarls_by_number[number] = snarl;
        hidden[number] = this->manager.is_trivial(snarl, this->graph);
    });
}

net_handle_t SnarlManagerDecomposition::snarl_net(const Snarl* snarl, bool forward) const {
    return directed(pack(NetType::SNARL, START, END, manager.snarl_number(snarl)), forward);
}

net_handle_t SnarlManagerDecomposition::chain_net(const Chain& chain) const {
    return pack(NetType::CHAIN, START, END, manager.snarl_number(chain.front().first));
}

const Chain& SnarlManagerDecomposition::chain_of_net(const net_handle_t& net) const {
    return *manager.chain_of(snarls_by_number.at(payload_of(net)));
}

net_handle_t SnarlManagerDecomposition::parent_snarl_net(const Snarl* parent) const {
    return parent == nullptr ? get_root() : snarl_net(parent, true);
}

bool SnarlManagerDecomposition::is_hidden(const Snarl* snarl) const {
    return hidden[manager.snarl_number(snarl)];
}

net_handle_t SnarlManagerDecomposition::chain_node_net(const Chain& chain, size_t index, bool along) const {
    handle_t node = layout.node(chain, index);
    return node_net(graph, along ? node : graph.flip(node));
}

net_handle_t SnarlManagerDecomposition::chain_snarl_net(const Chain& chain, size_t rank, bool along) const {
    const auto& [snarl, backward] = chain[rank];
    return snarl_net(snarl, along != backward);
}

vector<net_handle_t> SnarlManagerDecomposition::snarl_children(const Snarl* parent) const {
    vector<net_handle_t> children;
    for (const Chain& chain : manager.chains_of(parent)) {
        children.push_back(chain_net(chain));
    }
    for (const nid_t& id : trivial_chains.nodes_in(parent)) {
        children.push_back(pack(NetType::TRIVIAL_CHAIN, START, END, id));
    }
    return children;
}

vector<net_handle_t> SnarlManagerDecomposition::chain_children(const net_handle_t& traversal) const {
    const Chain& chain = chain_of_net(traversal);
    vector<net_handle_t> children;
    for (size_t rank = 0; rank < chain.size(); rank++) {
        children.push_back(chain_node_net(chain, rank, true));
        if (!is_hidden(chain[rank].first)) {
            children.push_back(chain_snarl_net(chain, rank, true));
        }
    }
    if (!layout.is_cyclic(chain)) {
        children.push_back(chain_node_net(chain, chain.size(), true));
    }

    if (is_backward(traversal)) {
        // Read the chain backward. A cyclic chain's list must still begin
        // with node index 0, which is both its bounds.
        std::ranges::reverse(children);
        if (layout.is_cyclic(chain)) {
            std::ranges::rotate(children, children.end() - 1);
        }
        for (net_handle_t& child : children) {
            child = flip(child);
        }
    } else if (!is_forward(traversal)) {
        // For a traversal that is neither start-to-end nor end-to-start, the
        // interface requires children start-to-end.
        for (net_handle_t& child : children) {
            child = with_traversal(child, START, END);
        }
    }
    return children;
}

optional<net_handle_t> SnarlManagerDecomposition::next_in_chain(const net_handle_t& here) const {
    if (is_snarl(here)) {
        if (ends_at(here) == TIP) {
            return nullopt;
        }
        const Snarl* snarl = get_snarl(here);
        const Chain& chain = *manager.chain_of(snarl);
        size_t rank = manager.chain_rank_of(snarl);
        // Leaving through the side that faces the chain's end goes along the
        // chain, to the node after the snarl.
        bool along = (ends_at(here) == ends_at(chain_snarl_net(chain, rank, true)));
        return chain_node_net(chain, along ? rank + 1 : rank, along);
    }

    optional<ChainLayout::Place> place = layout.locate(get_handle(here, &graph));
    if (!place) {
        // A trivial chain has nothing beside its node.
        return nullopt;
    }
    optional<size_t> rank = layout.snarl_ahead(*place);
    if (!rank.has_value()) {
        return nullopt;
    }
    const Chain& chain = *place->chain;
    if (is_hidden(chain[*rank].first)) {
        // Step over a hidden snarl to the node on its other side.
        return chain_node_net(chain, place->along ? *rank + 1 : *rank, place->along);
    }
    return chain_snarl_net(chain, *rank, place->along);
}

optional<net_handle_t> SnarlManagerDecomposition::sibling_entered_by(const handle_t& handle,
                                                                     const Snarl* parent) const {
    if (parent != nullptr) {
        net_handle_t parent_net = snarl_net(parent, true);
        if (handle == graph.flip(to_handle(graph, parent->start()))) {
            return get_bound(parent_net, false, false);
        }
        if (handle == to_handle(graph, parent->end())) {
            return get_bound(parent_net, true, false);
        }
        // From inside the snarl, an edge reaches a bound node only by reading
        // out of the snarl.
        assert(graph.get_id(handle) != parent->start().node_id() && graph.get_id(handle) != parent->end().node_id());
    }
    optional<ChainLayout::Place> place = layout.locate(handle);
    if (!place) {
        // A trivial chain, entered the way the node is read.
        return directed(pack(NetType::TRIVIAL_CHAIN, START, END, graph.get_id(handle)), !graph.get_is_reverse(handle));
    }
    optional<endpoint_t> entry = layout.entry_endpoint(*place);
    if (!entry) {
        // Only the chain's own snarls reach into its middle.
        return nullopt;
    }
    return directed(chain_net(*place->chain), *entry == START);
}

vector<net_handle_t> SnarlManagerDecomposition::siblings_ahead(const net_handle_t& here) const {
    vector<net_handle_t> siblings;
    const Snarl* parent = nullptr;
    handle_t exit;
    if (is_sentinel(here)) {
        if (!faces_in(here)) {
            return siblings;
        }
        parent = get_snarl(here);
        exit = get_handle(here, &graph);
    } else {
        if (ends_at(here) == TIP) {
            return siblings;
        }
        if (type_of(here) == NetType::CHAIN && layout.is_cyclic(chain_of_net(here))) {
            // Every edge on the node at a cyclic chain's bounds leads into the
            // chain's own snarls.
            return siblings;
        }
        net_handle_t up = get_parent(here);
        parent = is_root(up) ? nullptr : get_snarl(up);
        exit = get_handle(get_bound(here, ends_at(here) == END, false), &graph);
    }
    graph.follow_edges(exit, false, [&siblings, parent, this](const handle_t& next) {
        if (optional<net_handle_t> sibling = sibling_entered_by(next, parent)) {
            siblings.push_back(*sibling);
        }
    });
    return siblings;
}

net_handle_t SnarlManagerDecomposition::get_snarl_net(const Snarl* snarl) const {
    return snarl_net(snarl, true);
}

const Snarl* SnarlManagerDecomposition::get_snarl(const net_handle_t& net) const {
    if (type_of(net) != NetType::SNARL && type_of(net) != NetType::SENTINEL) {
        throw invalid_argument("error: net handle does not refer to a snarl");
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

net_handle_t SnarlManagerDecomposition::get_net(const handle_t& handle, const HandleGraph* handle_graph) const {
    return node_net(*handle_graph, handle);
}

handle_t SnarlManagerDecomposition::get_handle(const net_handle_t& net, const HandleGraph* handle_graph) const {
    using enum NetType;
    switch (type_of(net)) {
    case NODE:
    case TRIVIAL_CHAIN:
        return handle_graph->get_handle(payload_of(net), starts_at(net) == END);
    case SENTINEL:
        {
            const Snarl* snarl = get_snarl(net);
            if (starts_at(net) == START) {
                // The start bound's Visit reads into the snarl.
                handle_t in = to_handle(*handle_graph, snarl->start());
                return faces_in(net) ? in : handle_graph->flip(in);
            }
            // The end bound's Visit reads out of the snarl.
            handle_t out = to_handle(*handle_graph, snarl->end());
            return faces_in(net) ? handle_graph->flip(out) : out;
        }
    default:
        throw invalid_argument("error: trying to get a handle from a snarl, chain, or root");
    }
}

net_handle_t SnarlManagerDecomposition::get_parent(const net_handle_t& child) const {
    using enum NetType;
    switch (type_of(child)) {
    case ROOT:
        throw invalid_argument("error: trying to find the parent of the root");
    case SENTINEL:
        // The interface wants a snarl parent start-to-end.
        return pack(SNARL, START, END, payload_of(child));
    case SNARL:
        {
            const Snarl* snarl = get_snarl(child);
            const Chain& chain = *manager.chain_of(snarl);
            if (!is_forward(child) && !is_backward(child)) {
                return chain_net(chain);
            }
            // The chain reads the way the snarl reads, unless the snarl is
            // backward in the chain.
            bool backward_in_chain = chain[manager.chain_rank_of(snarl)].second;
            return directed(chain_net(chain), is_forward(child) != backward_in_chain);
        }
    case NODE:
        {
            optional<ChainLayout::Place> place = layout.locate(get_handle(child, &graph));
            if (!place) {
                // The node is a trivial chain by itself, traversed the way the
                // node is.
                return pack(TRIVIAL_CHAIN, starts_at(child), ends_at(child), payload_of(child));
            }
            return directed(chain_net(*place->chain), place->along);
        }
    case CHAIN:
        return parent_snarl_net(manager.parent_of(chain_of_net(child).front().first));
    case TRIVIAL_CHAIN:
        return parent_snarl_net(trivial_chains.parent_of(payload_of(child)));
    }
    throw invalid_argument("error: unknown net handle type");
}

net_handle_t SnarlManagerDecomposition::get_bound(const net_handle_t& snarl, bool get_end, bool face_in) const {
    using enum NetType;
    switch (type_of(snarl)) {
    case SNARL:
        {
            endpoint_t bound = get_end ? END : START;
            return pack(SENTINEL, bound, face_in ? opposite(bound) : bound, payload_of(snarl));
        }
    case CHAIN:
        {
            const Chain& chain = chain_of_net(snarl);
            handle_t node = layout.node(chain, get_end ? chain.size() : 0);
            // Along the chain, the start node reads in and the end node reads
            // out.
            return node_net(graph, face_in == get_end ? graph.flip(node) : node);
        }
    case TRIVIAL_CHAIN:
        {
            handle_t node = graph.get_handle(payload_of(snarl), false);
            return node_net(graph, face_in == get_end ? graph.flip(node) : node);
        }
    case ROOT:
        throw invalid_argument("error: trying to get the bounds of the root");
    default:
        throw invalid_argument("error: trying to get the bounds of a node or sentinel");
    }
}

net_handle_t SnarlManagerDecomposition::flip(const net_handle_t& net) const {
    if (is_sentinel(net)) {
        // Keep the bound and turn which way it faces.
        return with_traversal(net, starts_at(net), opposite(ends_at(net)));
    }
    return with_traversal(net, ends_at(net), starts_at(net));
}

net_handle_t SnarlManagerDecomposition::canonical(const net_handle_t& net) const {
    using enum NetType;
    switch (type_of(net)) {
    case ROOT:
        return get_root();
    case SENTINEL:
        // The bound facing in.
        return with_traversal(net, starts_at(net), opposite(starts_at(net)));
    default:
        {
            // The first traversal that for_each_traversal_impl() lists, which
            // is start-to-end whenever that is realizable, or start-to-end if
            // none is.
            net_handle_t found = with_traversal(net, START, END);
            for_each_traversal_impl(net, [&found](const net_handle_t& traversal) {
                found = traversal;
                return false;
            });
            return found;
        }
    }
}

handlegraph::SnarlDecomposition::endpoint_t SnarlManagerDecomposition::starts_at(const net_handle_t& traversal) const {
    return start_of(traversal);
}

handlegraph::SnarlDecomposition::endpoint_t SnarlManagerDecomposition::ends_at(const net_handle_t& traversal) const {
    return end_of(traversal);
}

net_handle_t SnarlManagerDecomposition::get_parent_traversal(const net_handle_t& traversal_start,
                                                             const net_handle_t& traversal_end) const {
    if (is_sentinel(traversal_start) && is_sentinel(traversal_end)) {
        if (payload_of(traversal_end) != payload_of(traversal_start)) {
            throw invalid_argument("error: traversal bounds are not sentinels of the same snarl");
        }
        return pack(NetType::SNARL, starts_at(traversal_start), starts_at(traversal_end), payload_of(traversal_start));
    }
    if (!is_node(traversal_start) || !is_node(traversal_end)) {
        throw invalid_argument("error: traversal bounds must be two sentinels or two nodes");
    }
    optional<ChainLayout::Place> start_place = layout.locate(get_handle(traversal_start, &graph));
    if (!start_place) {
        if (payload_of(traversal_end) != payload_of(traversal_start)) {
            throw invalid_argument("error: traversal bounds are not in the same trivial chain");
        }
        // A trivial chain is entered the way its start bound reads and left
        // the way its end bound reads.
        return pack(NetType::TRIVIAL_CHAIN, starts_at(traversal_start), ends_at(traversal_end),
                    payload_of(traversal_start));
    }
    optional<ChainLayout::Place> end_place = layout.locate(get_handle(traversal_end, &graph));
    if (!end_place || end_place->chain != start_place->chain) {
        throw invalid_argument("error: traversal bounds are not in the same chain");
    }
    optional<endpoint_t> entry = layout.entry_endpoint(*start_place);
    optional<endpoint_t> exit = layout.exit_endpoint(*end_place);
    if (!entry || !exit) {
        throw invalid_argument("error: traversal bounds are not the bounds of their chain");
    }
    return with_traversal(chain_net(*start_place->chain), *entry, *exit);
}

bool SnarlManagerDecomposition::for_each_child_impl(const net_handle_t& traversal,
                                                    const function<bool(const net_handle_t&)>& iteratee) const {
    using enum NetType;
    switch (type_of(traversal)) {
    case ROOT:
        return std::ranges::all_of(snarl_children(nullptr), iteratee);
    case SNARL:
        return std::ranges::all_of(snarl_children(get_snarl(traversal)), iteratee);
    case CHAIN:
        return std::ranges::all_of(chain_children(traversal), iteratee);
    case TRIVIAL_CHAIN:
        return iteratee(directed(pack(NODE, START, END, payload_of(traversal)), !is_backward(traversal)));
    default:
        // Nodes and sentinels have no children.
        return true;
    }
}

bool SnarlManagerDecomposition::for_each_traversal_impl(const net_handle_t& item,
                                                        const function<bool(const net_handle_t&)>& iteratee) const {
    vector<pair<endpoint_t, endpoint_t>> kinds;
    using enum NetType;
    switch (type_of(item)) {
    case ROOT:
        kinds.emplace_back(TIP, TIP);
        break;
    case SENTINEL:
        kinds.emplace_back(starts_at(item), ends_at(item));
        break;
    case NODE:
    case TRIVIAL_CHAIN:
        kinds = traversal_kinds({true, false, false});
        break;
    case SNARL:
        kinds = traversal_kinds(reachability(*get_snarl(item)));
        break;
    case CHAIN:
        kinds = traversal_kinds(reachability(chain_of_net(item)));
        break;
    }
    return std::ranges::all_of(kinds, [&item, &iteratee](const pair<endpoint_t, endpoint_t>& kind) {
        return iteratee(with_traversal(item, kind.first, kind.second));
    });
}

bool SnarlManagerDecomposition::follow_net_edges_impl(const net_handle_t& here, const HandleGraph*, bool go_left,
                                                      const function<bool(const net_handle_t&)>& iteratee) const {
    // Going left from here is going right from its flip.
    net_handle_t ahead = go_left ? flip(here) : here;
    using enum NetType;
    switch (type_of(ahead)) {
    case ROOT:
        return true;
    case NODE:
    case SNARL:
        {
            optional<net_handle_t> next = next_in_chain(ahead);
            return !next.has_value() || iteratee(*next);
        }
    default:
        return std::ranges::all_of(siblings_ahead(ahead), iteratee);
    }
}

bool SnarlManagerDecomposition::for_each_tippy_child_impl(const net_handle_t&,
                                                          const function<bool(const net_handle_t&)>&) const {
    throw domain_error("error: the snarl manager does not record internal tips");
}

}
