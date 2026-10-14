/** \file
 * Implements SnarlManagerDecomposition.
 */

#include "snarl_manager_decomposition.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <deque>
#include <iterator>
#include <stdexcept>
#include <unordered_set>

namespace vg {

using namespace std;

namespace {

using endpoint_t = handlegraph::SnarlDecomposition::endpoint_t;
using enum handlegraph::SnarlDecomposition::endpoint_t;

////////////////////////////////////////////////////////////////////////////
// Net handle layout
////////////////////////////////////////////////////////////////////////////

/*
 * A net handle is a 64-bit integer holding three fields, from the low bits
 * up:
 *
 * - The element type, a NetType, in TYPE_BITS bits.
 * - The traversal: the endpoint where it starts and the endpoint where it
 *   ends, each an endpoint_t in ENDPOINT_BITS bits.
 * - The payload, which says which element of that type the handle is for:
 *   - for a snarl, the snarl's SnarlManager::snarl_number();
 *   - for a sentinel, the snarl number of the snarl it bounds;
 *   - for a chain, the snarl number of the chain's first snarl;
 *   - for a node or a trivial chain, the node ID;
 *   - for the root, 0.
 *
 * A sentinel's traversal says which bound it is and which way it faces: it
 * starts at the endpoint of its bound, and ends at the other endpoint if it
 * faces into the snarl, or at the same endpoint if it faces out.
 */

/// @brief The type of net graph element a net handle is for.
enum class NetType : uint64_t {
    ROOT = 0,
    SNARL,
    CHAIN,
    /// @brief A chain holding a single free node.
    TRIVIAL_CHAIN,
    NODE,
    /// @brief A snarl's bound, facing into or out of the snarl.
    SENTINEL
};

constexpr uint64_t TYPE_BITS = 3;
constexpr uint64_t TYPE_MASK = (uint64_t{1} << TYPE_BITS) - 1;
constexpr uint64_t ENDPOINT_BITS = 2;
constexpr uint64_t ENDPOINT_MASK = (uint64_t{1} << ENDPOINT_BITS) - 1;
constexpr uint64_t START_SHIFT = TYPE_BITS;
constexpr uint64_t END_SHIFT = TYPE_BITS + ENDPOINT_BITS;
constexpr uint64_t PAYLOAD_SHIFT = TYPE_BITS + 2 * ENDPOINT_BITS;
static_assert(static_cast<uint64_t>(NetType::SENTINEL) <= TYPE_MASK);
static_assert(static_cast<uint64_t>(TIP) <= ENDPOINT_MASK);

/// @brief Pack the fields of a net handle.
/// @pre The payload fits in the bits above the traversal.
net_handle_t pack(NetType type, endpoint_t start, endpoint_t end, uint64_t payload) {
    assert(payload >> (64 - PAYLOAD_SHIFT) == 0);
    return handlegraph::as_net_handle((payload << PAYLOAD_SHIFT) |
                                      (static_cast<uint64_t>(end) << END_SHIFT) |
                                      (static_cast<uint64_t>(start) << START_SHIFT) |
                                      static_cast<uint64_t>(type));
}

/// @brief Unpack the element type of a net handle.
NetType type_of(const net_handle_t& net) {
    return static_cast<NetType>(handlegraph::as_integer(net) & TYPE_MASK);
}

/// @brief Unpack the endpoint where a net handle's traversal starts.
endpoint_t start_of(const net_handle_t& net) {
    return static_cast<endpoint_t>((handlegraph::as_integer(net) >> START_SHIFT) & ENDPOINT_MASK);
}

/// @brief Unpack the endpoint where a net handle's traversal ends.
endpoint_t end_of(const net_handle_t& net) {
    return static_cast<endpoint_t>((handlegraph::as_integer(net) >> END_SHIFT) & ENDPOINT_MASK);
}

/// @brief Unpack the payload of a net handle.
uint64_t payload_of(const net_handle_t& net) {
    return handlegraph::as_integer(net) >> PAYLOAD_SHIFT;
}

/// @brief Make a net handle for the same element with a different traversal.
net_handle_t with_traversal(const net_handle_t& net, endpoint_t start, endpoint_t end) {
    return pack(type_of(net), start, end, payload_of(net));
}

////////////////////////////////////////////////////////////////////////////
// Traversals and handles
////////////////////////////////////////////////////////////////////////////

/// @brief Swap START and END.
/// @pre The endpoint is not TIP.
endpoint_t opposite(endpoint_t endpoint) {
    assert(endpoint != TIP);
    return endpoint == START ? END : START;
}

/// @brief Make a net handle for the same element, traversed start-to-end if
/// forward and end-to-start otherwise.
net_handle_t directed(const net_handle_t& net, bool forward) {
    endpoint_t start = forward ? START : END;
    return with_traversal(net, start, opposite(start));
}

/// @brief Return true if a net handle's traversal is start-to-end.
bool is_forward(const net_handle_t& net) {
    return start_of(net) == START && end_of(net) == END;
}

/// @brief Return true if a net handle's traversal is end-to-start.
bool is_backward(const net_handle_t& net) {
    return start_of(net) == END && end_of(net) == START;
}

/// @brief Return true if a sentinel faces into its snarl.
bool faces_in(const net_handle_t& sentinel) {
    return start_of(sentinel) != end_of(sentinel);
}

/// @brief Get a node's net handle, read the way the graph handle reads.
net_handle_t node_net(const HandleGraph& graph, const handle_t& handle) {
    return directed(pack(NetType::NODE, START, END, graph.get_id(handle)), !graph.get_is_reverse(handle));
}

/// @brief Convert a Visit to a node into a graph handle.
handle_t to_handle(const HandleGraph& graph, const Visit& visit) {
    return graph.get_handle(visit.node_id(), visit.backward());
}

////////////////////////////////////////////////////////////////////////////
// Chains
////////////////////////////////////////////////////////////////////////////

/// @brief Return true if a chain's last snarl leads back into its first chain
/// node, read the same way.
bool is_cyclic(const Chain& chain) {
    return get_start_of(chain) == get_end_of(chain);
}

/// @brief Count a chain's distinct chain nodes: one more than its snarls, or
/// as many as its snarls if it is cyclic.
size_t chain_node_count(const Chain& chain) {
    return is_cyclic(chain) ? chain.size() : chain.size() + 1;
}

/// @brief Get the position of a chain's end bound, which for a cyclic chain is
/// its start bound at position 0.
size_t last_position(const Chain& chain) {
    return is_cyclic(chain) ? 0 : chain.size();
}

/**
 * @brief Get the chain node at a position, read along the chain.
 *
 * @param position A position up to chain.size(). For a cyclic chain,
 * chain.size() gives the node at position 0 again.
 */
handle_t chain_node(const HandleGraph& graph, const Chain& chain, size_t position) {
    if (position == chain.size()) {
        return to_handle(graph, get_end_of(chain));
    }
    const auto& [snarl, backward] = chain[position];
    return to_handle(graph, backward ? reverse(snarl->end()) : snarl->start());
}

/**
 * @brief Get the rank of the snarl that a chain node reads into.
 *
 * @param is_reverse True if the node is read against the chain.
 * @return The rank, or nullopt if the node is read out of the chain at one of
 * its bounds.
 */
optional<size_t> snarl_ahead(const Chain& chain, size_t position, bool is_reverse) {
    if (!is_reverse) {
        return position < chain.size() ? optional<size_t>(position) : nullopt;
    }
    if (position > 0) {
        return position - 1;
    }
    // Reading in reverse out of position 0 wraps around a cyclic chain.
    return is_cyclic(chain) ? optional<size_t>(chain.size() - 1) : nullopt;
}

/**
 * @brief Get the chain's endpoint that a chain node, read one way, enters
 * the chain through.
 *
 * @param is_reverse True if the node is read against the chain.
 * @return START for the start bound read along the chain, END for the end
 * bound read in reverse, and nullopt otherwise.
 */
optional<endpoint_t> entry_endpoint(const Chain& chain, size_t position, bool is_reverse) {
    if (!is_reverse && position == 0) {
        return START;
    }
    if (is_reverse && position == last_position(chain)) {
        return END;
    }
    return nullopt;
}

/**
 * @brief Get the chain's endpoint that a chain node, read one way, leaves
 * the chain through.
 *
 * @param is_reverse True if the node is read against the chain.
 * @return END for the end bound read along the chain, START for the start
 * bound read in reverse, and nullopt otherwise.
 */
optional<endpoint_t> exit_endpoint(const Chain& chain, size_t position, bool is_reverse) {
    if (!is_reverse && position == last_position(chain)) {
        return END;
    }
    if (is_reverse && position == 0) {
        return START;
    }
    return nullopt;
}

/// @brief The ways a walk can pass through a chain, judged from its snarls'
/// flags in the manager.
struct ChainPassage {
    /// @brief True if a walk can cross the chain from one bound to the other.
    bool crosses;
    /// @brief True if a walk entering at the start bound can leave there.
    bool turns_at_start;
    /// @brief True if a walk entering at the end bound can leave there.
    bool turns_at_end;
};

/**
 * @brief Return true if a walk that enters a chain at one bound can leave
 * through the same bound.
 *
 * @details The walk crosses each snarl it meets until one lets it turn
 * around on the side it entered that snarl from.
 *
 * @param begin The chain's (snarl, backward) pairs in the order the walk
 * meets them.
 * @param from_chain_start True if the walk enters at the chain's start bound.
 */
template<std::input_iterator Iterator>
bool can_turn_back(Iterator begin, Iterator end, bool from_chain_start) {
    for (Iterator it = begin; it != end; ++it) {
        const auto& [snarl, backward] = *it;
        // The walk enters the snarl's start if it and the snarl both go along
        // the chain or both go against it.
        if (from_chain_start != backward ? snarl->start_self_reachable() : snarl->end_self_reachable()) {
            return true;
        }
        if (!snarl->start_end_reachable()) {
            return false;
        }
    }
    return false;
}

/// @brief Work out the ways a walk can pass through a chain.
ChainPassage passage_through(const Chain& chain) {
    return {
        std::ranges::all_of(chain, [](const pair<const Snarl*, bool>& entry) {
            return entry.first->start_end_reachable();
        }),
        can_turn_back(chain.begin(), chain.end(), true),
        can_turn_back(chain.rbegin(), chain.rend(), false)
    };
}

/// @brief Return true if a node has a dead end: a side with no edges on it.
bool has_dead_end(const HandleGraph& graph, nid_t id) {
    handle_t handle = graph.get_handle(id, false);
    return graph.get_degree(handle, false) == 0 || graph.get_degree(handle, true) == 0;
}

}

////////////////////////////////////////////////////////////////////////////
// NodeDirectory
////////////////////////////////////////////////////////////////////////////

SnarlManagerDecomposition::NodeDirectory::NodeDirectory(const SnarlManager& manager, const HandleGraph& graph) {
    manager.for_each_chain([this, &manager, &graph](const Chain* chain) {
        const Snarl* parent = manager.parent_of(chain->front().first);
        for (size_t position = 0; position < chain_node_count(*chain); position++) {
            handle_t node = chain_node(graph, *chain, position);
            entries[graph.get_id(node)] = Entry {parent, ChainPosition {chain, position, graph.get_is_reverse(node)}};
        }
    });
    manager.for_each_snarl_unindexed([this, &manager, &graph](const Snarl* snarl) {
        for (const id_t& id : manager.shallow_contents(snarl, graph, false).first) {
            if (!entries.contains(id)) {
                entries[id] = Entry {snarl, nullopt};
                free_nodes_by_snarl[snarl].push_back(id);
            }
        }
    });
    graph.for_each_handle([this, &graph](const handle_t& handle) {
        nid_t id = graph.get_id(handle);
        if (!entries.contains(id)) {
            entries[id] = Entry {nullptr, nullopt};
            free_nodes_by_snarl[nullptr].push_back(id);
        }
    });
    for (auto& [snarl, nodes] : free_nodes_by_snarl) {
        std::ranges::sort(nodes);
    }
}

optional<SnarlManagerDecomposition::ChainPosition>
SnarlManagerDecomposition::NodeDirectory::chain_position(nid_t id) const {
    return entries.at(id).chain_position;
}

const Snarl* SnarlManagerDecomposition::NodeDirectory::containing_snarl(nid_t id) const {
    return entries.at(id).containing_snarl;
}

const vector<nid_t>& SnarlManagerDecomposition::NodeDirectory::free_nodes(const Snarl* snarl) const {
    static const vector<nid_t> none;
    auto found = free_nodes_by_snarl.find(snarl);
    return found == free_nodes_by_snarl.end() ? none : found->second;
}

////////////////////////////////////////////////////////////////////////////
// TraversalSet
////////////////////////////////////////////////////////////////////////////

void SnarlManagerDecomposition::TraversalSet::add(endpoint_t start, endpoint_t end) {
    bits |= uint16_t{1} << (3 * start + end);
}

void SnarlManagerDecomposition::TraversalSet::add_both_ways(endpoint_t one, endpoint_t other) {
    add(one, other);
    add(other, one);
}

bool SnarlManagerDecomposition::TraversalSet::contains(endpoint_t start, endpoint_t end) const {
    return (bits >> (3 * start + end)) & 1;
}

SnarlManagerDecomposition::TraversalSet
SnarlManagerDecomposition::TraversalSet::for_chain(const vector<pair<TraversalSet, bool>>& shown_snarls) {
    TraversalSet traversals;
    if (std::ranges::all_of(shown_snarls, [](const pair<TraversalSet, bool>& entry) {
            return entry.first.contains(START, END);
        })) {
        traversals.add_both_ways(START, END);
    }
    for (const auto& [snarl_traversals, backward] : shown_snarls) {
        // The side of the snarl that faces the chain's start.
        endpoint_t toward_start = backward ? END : START;
        endpoint_t toward_end = opposite(toward_start);
        // TODO: A distance index lets a chain turn around through any of its
        // snarls, without checking that the snarl can be reached from the
        // chain's bound, so we do too.
        if (snarl_traversals.contains(toward_start, toward_start)) {
            traversals.add(START, START);
        }
        if (snarl_traversals.contains(toward_end, toward_end)) {
            traversals.add(END, END);
        }
    }
    if (reaches_tip(shown_snarls, true)) {
        traversals.add_both_ways(START, TIP);
    }
    if (reaches_tip(shown_snarls, false)) {
        traversals.add_both_ways(END, TIP);
    }
    if (links_tips(shown_snarls)) {
        traversals.add(TIP, TIP);
    }
    return traversals;
}

bool SnarlManagerDecomposition::TraversalSet::reaches_tip(const vector<pair<TraversalSet, bool>>& shown_snarls,
                                                          bool from_start) {
    for (size_t i = 0; i < shown_snarls.size(); i++) {
        const auto& [snarl_traversals, backward] = shown_snarls[from_start ? i : shown_snarls.size() - 1 - i];
        // The walk enters the snarl's start if it and the snarl both go along
        // the chain or both go against it.
        if (endpoint_t entered = (backward == from_start) ? END : START; snarl_traversals.contains(entered, TIP)) {
            return true;
        }
        if (!snarl_traversals.contains(START, END)) {
            return false;
        }
    }
    return false;
}

bool SnarlManagerDecomposition::TraversalSet::links_tips(const vector<pair<TraversalSet, bool>>& shown_snarls) {
    // True when a walk from a tip in an earlier snarl can reach the snarl
    // being looked at, going along the chain.
    // TODO: Walks that turn around between two tips are not found.
    bool leaving_tip = false;
    for (const auto& [snarl_traversals, backward] : shown_snarls) {
        endpoint_t toward_start = backward ? END : START;
        endpoint_t toward_end = opposite(toward_start);
        if (snarl_traversals.contains(TIP, TIP) || (leaving_tip && snarl_traversals.contains(toward_start, TIP))) {
            return true;
        }
        leaving_tip = snarl_traversals.contains(TIP, toward_end) ||
            (leaving_tip && snarl_traversals.contains(START, END));
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////
// SnarlSurveyor
////////////////////////////////////////////////////////////////////////////

SnarlManagerDecomposition::SnarlSurveyor::SnarlSurveyor(const HandleGraph& graph, const NodeDirectory& directory) :
    graph(graph), directory(directory) {
}

vector<handle_t> SnarlManagerDecomposition::SnarlSurveyor::exits_of_child(const handle_t& entering) const {
    optional<ChainPosition> place = directory.chain_position(graph.get_id(entering));
    if (!place) {
        return {entering};
    }
    const Chain& chain = *place->chain;
    optional<endpoint_t> entry = entry_endpoint(chain, place->position,
                                                graph.get_is_reverse(entering) != place->is_reverse);
    assert(entry.has_value());
    ChainPassage passage = passage_through(chain);
    handle_t start_out = graph.flip(chain_node(graph, chain, 0));
    handle_t end_out = chain_node(graph, chain, chain.size());
    vector<handle_t> exits;
    if (passage.crosses) {
        exits.push_back(*entry == START ? end_out : start_out);
    }
    if (*entry == START ? passage.turns_at_start : passage.turns_at_end) {
        exits.push_back(*entry == START ? start_out : end_out);
    }
    return exits;
}

bool SnarlManagerDecomposition::SnarlSurveyor::reaches_tip(const Snarl* snarl, const vector<handle_t>& starts) const {
    unordered_set<handle_t> seen(starts.begin(), starts.end());
    deque<handle_t> to_leave(starts.begin(), starts.end());
    while (!to_leave.empty()) {
        handle_t here = to_leave.front();
        to_leave.pop_front();
        bool stopped = !graph.follow_edges(here, false, [this, snarl, &seen, &to_leave](const handle_t& next) {
            nid_t id = graph.get_id(next);
            if (id == snarl->start().node_id() || id == snarl->end().node_id()) {
                // The walk leaves the snarl here.
                return true;
            }
            if (!directory.chain_position(id) && has_dead_end(graph, id)) {
                // Coming in through the side with an edge, a tip node's dead
                // end is ahead.
                return false;
            }
            for (const handle_t& exit : exits_of_child(next)) {
                if (seen.insert(exit).second) {
                    to_leave.push_back(exit);
                }
            }
            return true;
        });
        if (stopped) {
            return true;
        }
    }
    return false;
}

SnarlManagerDecomposition::TraversalSet SnarlManagerDecomposition::SnarlSurveyor::survey(const Snarl* snarl) const {
    TraversalSet traversals;
    if (snarl->start_end_reachable()) {
        traversals.add_both_ways(START, END);
    }
    if (snarl->start_self_reachable()) {
        traversals.add(START, START);
    }
    if (snarl->end_self_reachable()) {
        traversals.add(END, END);
    }
    if (reaches_tip(snarl, {to_handle(graph, snarl->start())})) {
        traversals.add_both_ways(START, TIP);
    }
    if (reaches_tip(snarl, {graph.flip(to_handle(graph, snarl->end()))})) {
        traversals.add_both_ways(END, TIP);
    }
    // A walk from a tip starts by reading a tip node away from its dead end.
    vector<handle_t> from_tips;
    for (const nid_t& id : directory.free_nodes(snarl)) {
        for (bool is_reverse : {false, true}) {
            handle_t handle = graph.get_handle(id, is_reverse);
            if (graph.get_degree(handle, true) == 0 && graph.get_degree(handle, false) != 0) {
                from_tips.push_back(handle);
            }
        }
    }
    if (!from_tips.empty() && reaches_tip(snarl, from_tips)) {
        traversals.add(TIP, TIP);
    }
    return traversals;
}

////////////////////////////////////////////////////////////////////////////
// SnarlManagerDecomposition
////////////////////////////////////////////////////////////////////////////

SnarlManagerDecomposition::SnarlManagerDecomposition(const SnarlManager& manager, const HandleGraph& graph) :
    manager(manager), graph(graph), directory(manager, graph), snarls(manager.num_snarls()) {

    manager.for_each_snarl_unindexed([this](const Snarl* snarl) {
        if (snarl->start().node_id() == snarl->end().node_id() &&
            snarl->start().backward() != snarl->end().backward()) {
            throw invalid_argument("error: SnarlManagerDecomposition cannot present unary snarls");
        }
        SnarlRecord& record = snarls[this->manager.snarl_number(snarl)];
        record.snarl = snarl;
        // TODO: A distance index hides every snarl with no children, even one
        // that is not an ultrabubble, so we do too.
        record.hidden = this->manager.chains_of(snarl).empty() && directory.free_nodes(snarl).empty();
    });
    SnarlSurveyor surveyor(graph, directory);
    for (SnarlRecord& record : snarls) {
        if (!record.hidden) {
            record.traversals = surveyor.survey(record.snarl);
        }
    }
    manager.for_each_chain([this](const Chain* chain) {
        // TODO: A distance index treats every hidden snarl as crossable, with
        // no turnarounds, so we look only at shown snarls.
        vector<pair<TraversalSet, bool>> shown_snarls;
        for (const auto& [snarl, backward] : *chain) {
            if (const SnarlRecord& record = snarls[this->manager.snarl_number(snarl)]; !record.hidden) {
                shown_snarls.emplace_back(record.traversals, backward);
            }
        }
        chains[chain] = TraversalSet::for_chain(shown_snarls);
    });
}

const SnarlManagerDecomposition::SnarlRecord& SnarlManagerDecomposition::record_of(const net_handle_t& net) const {
    return snarls.at(payload_of(net));
}

net_handle_t SnarlManagerDecomposition::snarl_net(const Snarl* snarl, bool forward) const {
    return directed(pack(NetType::SNARL, START, END, manager.snarl_number(snarl)), forward);
}

net_handle_t SnarlManagerDecomposition::chain_net(const Chain& chain) const {
    return pack(NetType::CHAIN, START, END, manager.snarl_number(chain.front().first));
}

const Chain& SnarlManagerDecomposition::chain_of_net(const net_handle_t& net) const {
    return *manager.chain_of(record_of(net).snarl);
}

net_handle_t SnarlManagerDecomposition::chain_node_net(const Chain& chain, size_t position, bool is_reverse) const {
    if (position == chain_node_count(chain)) {
        // TODO: A distance index reaches the start of a cyclic chain, going
        // around the chain, as its end bound facing out, so we do too.
        return get_bound(chain_net(chain), true, false);
    }
    handle_t node = chain_node(graph, chain, position);
    return node_net(graph, is_reverse ? graph.flip(node) : node);
}

net_handle_t SnarlManagerDecomposition::chain_snarl_net(const Chain& chain, size_t rank, bool is_reverse) const {
    const auto& [snarl, backward] = chain[rank];
    return snarl_net(snarl, is_reverse == backward);
}

vector<net_handle_t> SnarlManagerDecomposition::snarl_children(const Snarl* parent) const {
    vector<net_handle_t> children;
    for (const Chain& chain : manager.chains_of(parent)) {
        children.push_back(chain_net(chain));
    }
    for (const nid_t& id : directory.free_nodes(parent)) {
        children.push_back(pack(NetType::TRIVIAL_CHAIN, START, END, id));
    }
    return children;
}

vector<net_handle_t> SnarlManagerDecomposition::chain_children(const net_handle_t& traversal) const {
    const Chain& chain = chain_of_net(traversal);
    vector<net_handle_t> children;
    for (size_t rank = 0; rank < chain.size(); rank++) {
        children.push_back(chain_node_net(chain, rank, false));
        if (!snarls[manager.snarl_number(chain[rank].first)].hidden) {
            children.push_back(chain_snarl_net(chain, rank, false));
        }
    }
    if (!is_cyclic(chain)) {
        children.push_back(chain_node_net(chain, chain.size(), false));
    }

    if (is_backward(traversal)) {
        // Read the chain in reverse. A cyclic chain's list must still begin
        // with its position 0 node, which is both its bounds.
        std::ranges::reverse(children);
        if (is_cyclic(chain)) {
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
            // The walk stops at the tip.
            return nullopt;
        }
        const Snarl* snarl = get_snarl(here);
        const Chain& chain = *manager.chain_of(snarl);
        size_t rank = manager.chain_rank_of(snarl);
        // Leaving through the snarl's side that faces the chain's start goes
        // against the chain, to the chain node before the snarl.
        bool is_reverse = (ends_at(here) == ends_at(chain_snarl_net(chain, rank, true)));
        return chain_node_net(chain, is_reverse ? rank : rank + 1, is_reverse);
    }

    handle_t handle = get_handle(here, &graph);
    optional<ChainPosition> place = directory.chain_position(graph.get_id(handle));
    if (!place) {
        // A trivial chain's node has nothing beside it.
        return nullopt;
    }
    const Chain& chain = *place->chain;
    bool is_reverse = graph.get_is_reverse(handle) != place->is_reverse;
    optional<size_t> rank = snarl_ahead(chain, place->position, is_reverse);
    if (!rank.has_value()) {
        return nullopt;
    }
    if (snarls[manager.snarl_number(chain[*rank].first)].hidden) {
        // Step over a hidden snarl to the chain node on its other side.
        return chain_node_net(chain, is_reverse ? *rank : *rank + 1, is_reverse);
    }
    return chain_snarl_net(chain, *rank, is_reverse);
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
    nid_t id = graph.get_id(handle);
    optional<ChainPosition> place = directory.chain_position(id);
    if (!place) {
        // A trivial chain, entered the way the node is read.
        return directed(pack(NetType::TRIVIAL_CHAIN, START, END, id), !graph.get_is_reverse(handle));
    }
    bool is_reverse = graph.get_is_reverse(handle) != place->is_reverse;
    optional<endpoint_t> entry = entry_endpoint(*place->chain, place->position, is_reverse);
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
        parent = record_of(here).snarl;
        exit = get_handle(here, &graph);
    } else {
        if (ends_at(here) == TIP) {
            return siblings;
        }
        if (type_of(here) == NetType::CHAIN && is_cyclic(chain_of_net(here))) {
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
    if (snarls.at(manager.snarl_number(snarl)).hidden) {
        throw invalid_argument("error: snarl is hidden from the decomposition");
    }
    return snarl_net(snarl, true);
}

const Snarl* SnarlManagerDecomposition::get_snarl(const net_handle_t& net) const {
    if (type_of(net) != NetType::SNARL) {
        throw invalid_argument("error: net handle is not for a snarl");
    }
    return record_of(net).snarl;
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
        return handle_graph->get_handle(payload_of(net), ends_at(net) == START);
    case SENTINEL:
        {
            const Snarl* snarl = record_of(net).snarl;
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
            optional<ChainPosition> place = directory.chain_position(payload_of(child));
            if (!place) {
                // The node is a trivial chain by itself, traversed the way the
                // node is.
                return pack(TRIVIAL_CHAIN, starts_at(child), ends_at(child), payload_of(child));
            }
            if (!is_forward(child) && !is_backward(child)) {
                return chain_net(*place->chain);
            }
            bool is_reverse = graph.get_is_reverse(get_handle(child, &graph)) != place->is_reverse;
            return directed(chain_net(*place->chain), !is_reverse);
        }
    case CHAIN:
        {
            const Snarl* parent = manager.parent_of(chain_of_net(child).front().first);
            return parent == nullptr ? get_root() : snarl_net(parent, true);
        }
    case TRIVIAL_CHAIN:
        {
            const Snarl* parent = directory.containing_snarl(payload_of(child));
            return parent == nullptr ? get_root() : snarl_net(parent, true);
        }
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
            if (get_end && is_cyclic(chain)) {
                // TODO: A distance index gives a cyclic chain's end bound as
                // its start bound's node read end-to-end, or start-to-start
                // if the chain reads that node in reverse, so we do too.
                endpoint_t out = ends_at(node_net(graph, chain_node(graph, chain, 0)));
                endpoint_t end = face_in ? opposite(out) : out;
                return pack(NODE, end, end, graph.get_id(chain_node(graph, chain, 0)));
            }
            handle_t node = chain_node(graph, chain, get_end ? chain.size() : 0);
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
        return with_traversal(net, starts_at(net), opposite(starts_at(net)));
    default:
        {
            // Use the first realizable traversal in the order a distance
            // index prefers.
            static const array<pair<endpoint_t, endpoint_t>, 6> preferred {{
                {START, END}, {START, TIP}, {END, TIP}, {START, START}, {END, END}, {TIP, TIP}
            }};
            TraversalSet realizable = traversals_of(net);
            for (const auto& [start, end] : preferred) {
                if (realizable.contains(start, end)) {
                    return with_traversal(net, start, end);
                }
            }
            return with_traversal(net, START, END);
        }
    }
}

handlegraph::SnarlDecomposition::endpoint_t SnarlManagerDecomposition::starts_at(const net_handle_t& traversal) const {
    return start_of(traversal);
}

handlegraph::SnarlDecomposition::endpoint_t SnarlManagerDecomposition::ends_at(const net_handle_t& traversal) const {
    return end_of(traversal);
}

pair<net_handle_t, endpoint_t> SnarlManagerDecomposition::traversal_bound(const net_handle_t& net, bool is_start) const {
    using enum NetType;
    switch (type_of(net)) {
    case SENTINEL:
        return {pack(SNARL, START, END, payload_of(net)), starts_at(net)};
    case SNARL:
    case CHAIN:
    case TRIVIAL_CHAIN:
        {
            net_handle_t parent = get_parent(with_traversal(net, START, END));
            if (is_root(parent)) {
                throw invalid_argument("error: the root has no traversals that start or end at a child");
            }
            return {with_traversal(parent, START, END), TIP};
        }
    case NODE:
        {
            optional<ChainPosition> place = directory.chain_position(payload_of(net));
            if (!place) {
                // A free node's trivial chain is entered where the node's
                // traversal starts and left where it ends.
                return {pack(TRIVIAL_CHAIN, START, END, payload_of(net)), is_start ? starts_at(net) : ends_at(net)};
            }
            const Chain& chain = *place->chain;
            bool is_reverse = graph.get_is_reverse(get_handle(net, &graph)) != place->is_reverse;
            optional<endpoint_t> endpoint = is_start ? entry_endpoint(chain, place->position, is_reverse)
                                                     : exit_endpoint(chain, place->position, is_reverse);
            if (!endpoint) {
                throw invalid_argument("error: node traversal does not read into or out of its chain at a bound");
            }
            return {chain_net(chain), *endpoint};
        }
    default:
        throw invalid_argument("error: the root cannot start or end a traversal");
    }
}

net_handle_t SnarlManagerDecomposition::get_parent_traversal(const net_handle_t& traversal_start,
                                                             const net_handle_t& traversal_end) const {
    auto [start_parent, start] = traversal_bound(traversal_start, true);
    auto [end_parent, end] = traversal_bound(traversal_end, false);
    if (start_parent != end_parent) {
        throw invalid_argument("error: traversal bounds are not in the same snarl or chain");
    }
    if (!traversals_of(start_parent).contains(start, end)) {
        throw invalid_argument("error: traversal between the given bounds is not realizable");
    }
    return with_traversal(start_parent, start, end);
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

SnarlManagerDecomposition::TraversalSet SnarlManagerDecomposition::traversals_of(const net_handle_t& item) const {
    TraversalSet traversals;
    using enum NetType;
    switch (type_of(item)) {
    case ROOT:
        traversals.add(TIP, TIP);
        break;
    case SENTINEL:
        traversals.add(starts_at(item), ends_at(item));
        break;
    case NODE:
    case TRIVIAL_CHAIN:
        traversals.add_both_ways(START, END);
        break;
    case SNARL:
        traversals = record_of(item).traversals;
        break;
    case CHAIN:
        traversals = chains.at(&chain_of_net(item));
        break;
    }
    return traversals;
}

bool SnarlManagerDecomposition::for_each_traversal_impl(const net_handle_t& item,
                                                        const function<bool(const net_handle_t&)>& iteratee) const {
    return traversals_of(item).for_each([&item, &iteratee](endpoint_t start, endpoint_t end) {
        return iteratee(with_traversal(item, start, end));
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
            if (!next.has_value()) {
                return true;
            }
            if (go_left && is_node(here) && is_node(*next)) {
                // TODO: A distance index reads a chain node it reaches by going
                // left from a neighbouring chain node the way the traversal
                // reads, so we do too.
                return iteratee(flip(*next));
            }
            return iteratee(*next);
        }
    default:
        return std::ranges::all_of(siblings_ahead(ahead), iteratee);
    }
}

}
