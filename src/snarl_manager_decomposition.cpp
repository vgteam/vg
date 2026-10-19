/** \file
 * The net handle encoding of SnarlManagerDecomposition, and the walks that
 * survey which traversals of each snarl and chain are realizable.
 */

#include "snarl_manager_decomposition.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <deque>
#include <set>
#include <stdexcept>
#include <unordered_set>

namespace vg {

using namespace std;

namespace {

using endpoint_t = handlegraph::SnarlDecomposition::endpoint_t;
using enum handlegraph::SnarlDecomposition::endpoint_t;

////////////////////////////////////////////////////////////////////////////
// Net handle encoding
////////////////////////////////////////////////////////////////////////////

/*
 * A net handle is a 64-bit integer holding three fields, from the low bits
 * up:
 *
 * - The element type, a NetType, in TYPE_BITS bits.
 * - The traversal, as the endpoint where it starts and the endpoint where it
 *   ends, each an endpoint_t in ENDPOINT_BITS bits.
 * - The payload, which says which element of that type the handle is for:
 *   - for a snarl, its SnarlManager::snarl_number();
 *   - for a sentinel, the snarl number of the snarl it bounds;
 *   - for a chain, the snarl number of the first snarl in the chain;
 *   - for a node or a trivial chain, the node ID;
 *   - for the root, 0.
 *
 * The traversal of a sentinel says which bound it is and which way it faces.
 * It starts at the endpoint of its bound. It ends at the other endpoint if
 * it faces into the snarl, and at the same endpoint if it faces out.
 */

/// The type of net graph element a net handle is for.
enum class NetType : uint64_t {
    ROOT = 0,
    SNARL,
    CHAIN,
    /// The trivial chain of a free node.
    TRIVIAL_CHAIN,
    NODE,
    /// A bound of a snarl, facing into or out of the snarl.
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

/**
 * Pack the fields of a net handle.
 *
 * @pre The payload fits in the bits above the traversal.
 */
net_handle_t pack(NetType type, endpoint_t start, endpoint_t end, uint64_t payload) {
    assert(payload >> (64 - PAYLOAD_SHIFT) == 0);
    return handlegraph::as_net_handle((payload << PAYLOAD_SHIFT) |
                                      (static_cast<uint64_t>(end) << END_SHIFT) |
                                      (static_cast<uint64_t>(start) << START_SHIFT) |
                                      static_cast<uint64_t>(type));
}

/// Unpack the element type of a net handle.
NetType type_of(const net_handle_t& net) {
    return static_cast<NetType>(handlegraph::as_integer(net) & TYPE_MASK);
}

/// Unpack the endpoint where the traversal of a net handle starts.
endpoint_t start_of(const net_handle_t& net) {
    return static_cast<endpoint_t>((handlegraph::as_integer(net) >> START_SHIFT) & ENDPOINT_MASK);
}

/// Unpack the endpoint where the traversal of a net handle ends.
endpoint_t end_of(const net_handle_t& net) {
    return static_cast<endpoint_t>((handlegraph::as_integer(net) >> END_SHIFT) & ENDPOINT_MASK);
}

/// Unpack the payload of a net handle.
uint64_t payload_of(const net_handle_t& net) {
    return handlegraph::as_integer(net) >> PAYLOAD_SHIFT;
}

/// Make a net handle for the same element with a different traversal.
net_handle_t with_traversal(const net_handle_t& net, endpoint_t start, endpoint_t end) {
    return pack(type_of(net), start, end, payload_of(net));
}

////////////////////////////////////////////////////////////////////////////
// Traversals and handles
////////////////////////////////////////////////////////////////////////////

/**
 * Swap START and END.
 *
 * @pre The endpoint is START or END.
 */
endpoint_t opposite(endpoint_t endpoint) {
    assert(endpoint != TIP);
    return endpoint == START ? END : START;
}

/// Make a net handle for the same element, traversed start-to-end if forward
/// and end-to-start otherwise.
net_handle_t directed(const net_handle_t& net, bool forward) {
    endpoint_t start = forward ? START : END;
    return with_traversal(net, start, opposite(start));
}

/// Return true if the traversal of a net handle is start-to-end.
bool is_forward(const net_handle_t& net) {
    return start_of(net) == START && end_of(net) == END;
}

/// Return true if the traversal of a net handle is end-to-start.
bool is_backward(const net_handle_t& net) {
    return start_of(net) == END && end_of(net) == START;
}

/// Return true if a sentinel faces into its snarl.
bool faces_in(const net_handle_t& sentinel) {
    return start_of(sentinel) != end_of(sentinel);
}

/// Get the net handle for the node of a graph handle, start-to-end if the
/// handle is forward and end-to-start if it is reverse.
net_handle_t node_net(const HandleGraph& graph, const handle_t& handle) {
    return directed(pack(NetType::NODE, START, END, graph.get_id(handle)), !graph.get_is_reverse(handle));
}

/// Convert a Visit to a node into a graph handle.
handle_t to_handle(const HandleGraph& graph, const Visit& visit) {
    return graph.get_handle(visit.node_id(), visit.backward());
}

/// Return true if the left side of a handle is a dead end.
bool left_is_dead_end(const HandleGraph& graph, const handle_t& handle) {
    return graph.get_degree(handle, true) == 0;
}

////////////////////////////////////////////////////////////////////////////
// Chains
////////////////////////////////////////////////////////////////////////////

/// Return true if a chain is a looping chain.
bool is_looping(const Chain& chain) {
    return get_start_of(chain) == get_end_of(chain);
}

/// Count the distinct chain nodes of a chain. A looping chain has as many as
/// it has snarls, and any other chain has one more.
size_t chain_node_count(const Chain& chain) {
    return is_looping(chain) ? chain.size() : chain.size() + 1;
}

/// Get the position of the end bound of a chain.
size_t last_position(const Chain& chain) {
    return is_looping(chain) ? 0 : chain.size();
}

/**
 * Get the chain node at a position, oriented forward along the chain.
 *
 * @param position A position up to chain.size(). For a looping chain,
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
 * Get the rank of the snarl that a walk enters on leaving a chain node.
 *
 * @param position The position of the chain node.
 * @param is_reverse True if the walk goes in reverse along the chain.
 * @return The rank of the snarl in the chain, or nullopt if the walk leaves
 * the chain at one of its bounds.
 */
optional<size_t> snarl_ahead(const Chain& chain, size_t position, bool is_reverse) {
    if (!is_reverse) {
        return position < chain.size() ? optional<size_t>(position) : nullopt;
    }
    if (position > 0) {
        return position - 1;
    }
    // Going in reverse from position 0, a walk on a looping chain enters the
    // last snarl.
    return is_looping(chain) ? optional<size_t>(chain.size() - 1) : nullopt;
}

/**
 * Get the endpoint of a chain where a walk that visits a chain node enters
 * the chain.
 *
 * @param position The position of the chain node.
 * @param is_reverse True if the walk goes in reverse along the chain.
 * @return START for a walk going forward along the chain from its start
 * bound, END for a walk going in reverse along the chain from its end bound,
 * and nullopt for a walk that enters the chain node from inside the chain.
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
 * Get the endpoint of a chain where a walk that visits a chain node leaves
 * the chain.
 *
 * @param position The position of the chain node.
 * @param is_reverse True if the walk goes in reverse along the chain.
 * @return END for a walk going forward along the chain from its end bound,
 * START for a walk going in reverse along the chain from its start bound,
 * and nullopt for a walk that leaves the chain node into the chain.
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

/**
 * List the snarls of a manager, indexed by SnarlManager::snarl_number().
 *
 * @throws std::invalid_argument If a snarl is unary.
 */
vector<const Snarl*> number_snarls(const SnarlManager& manager) {
    vector<const Snarl*> snarls(manager.num_snarls());
    manager.for_each_snarl_unindexed([&](const Snarl* snarl) {
        if (snarl->start().node_id() == snarl->end().node_id() &&
            snarl->start().backward() != snarl->end().backward()) {
            throw invalid_argument("error: SnarlManagerDecomposition cannot present unary snarls");
        }
        snarls[manager.snarl_number(snarl)] = snarl;
    });
    return snarls;
}

}

////////////////////////////////////////////////////////////////////////////
// NodeDirectory
////////////////////////////////////////////////////////////////////////////

SnarlManagerDecomposition::NodeDirectory::NodeDirectory(const SnarlManager& manager, const HandleGraph& graph) {
    // Place the chain nodes first, so that every node not yet placed after
    // that is a free node.
    manager.for_each_chain([this, &manager, &graph](const Chain* chain) {
        const Snarl* parent = manager.parent_of(chain->front().first);
        for (size_t position = 0; position < chain_node_count(*chain); position++) {
            handle_t node = chain_node(graph, *chain, position);
            entries[graph.get_id(node)] = Entry {parent, ChainPosition {chain, position, graph.get_is_reverse(node)}};
        }
    });
    // The shallow contents of a snarl are the chain nodes of its child
    // chains, which already have places, and its free nodes.
    manager.for_each_snarl_unindexed([this, &manager, &graph](const Snarl* snarl) {
        for (const id_t& id : manager.shallow_contents(snarl, graph, false).first) {
            if (!entries.contains(id)) {
                entries[id] = Entry {snarl, nullopt};
                free_nodes_by_snarl[snarl].push_back(id);
            }
        }
    });
    // Each node in no snarl is a free node of the root.
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
    static const vector<nid_t> no_nodes;
    auto found = free_nodes_by_snarl.find(snarl);
    return found == free_nodes_by_snarl.end() ? no_nodes : found->second;
}

////////////////////////////////////////////////////////////////////////////
// TraversalSet
////////////////////////////////////////////////////////////////////////////

void SnarlManagerDecomposition::TraversalSet::add(endpoint_t start, endpoint_t end) {
    bits |= uint16_t{1} << (3 * start + end);
}

bool SnarlManagerDecomposition::TraversalSet::contains(endpoint_t start, endpoint_t end) const {
    return (bits >> (3 * start + end)) & 1;
}

////////////////////////////////////////////////////////////////////////////
// TraversalSurvey
////////////////////////////////////////////////////////////////////////////

SnarlManagerDecomposition::TraversalSurvey::TraversalSurvey(const SnarlManager& manager, const HandleGraph& graph,
                                                            const NodeDirectory& directory) {
    Surveyor surveyor(manager, graph, directory, *this);

    // Survey children before parents. In reverse preorder, every snarl comes
    // after the snarls inside it, and the child chains of a snarl hold only
    // snarls inside it. The top-level chains go last.
    vector<const Snarl*> preorder;
    manager.for_each_snarl_preorder([&preorder](const Snarl* snarl) {
        preorder.push_back(snarl);
    });
    for (auto it = preorder.rbegin(); it != preorder.rend(); ++it) {
        for (const Chain& chain : manager.chains_of(*it)) {
            chains[&chain] = surveyor.survey_chain(chain);
        }
        snarls[*it] = surveyor.survey_snarl(*it);
    }
    for (const Chain& chain : manager.chains_of(nullptr)) {
        chains[&chain] = surveyor.survey_chain(chain);
    }
}

const SnarlManagerDecomposition::TraversalSet&
SnarlManagerDecomposition::TraversalSurvey::of_snarl(const Snarl* snarl) const {
    return snarls.at(snarl);
}

const SnarlManagerDecomposition::TraversalSet&
SnarlManagerDecomposition::TraversalSurvey::of_chain(const Chain* chain) const {
    return chains.at(chain);
}

SnarlManagerDecomposition::TraversalSurvey::Surveyor::Surveyor(const SnarlManager& manager, const HandleGraph& graph,
                                                           const NodeDirectory& directory,
                                                           const TraversalSurvey& survey) :
    manager(manager), graph(graph), directory(directory), survey(survey) {
}

bool SnarlManagerDecomposition::TraversalSurvey::Surveyor::enter_child(const handle_t& entering, vector<handle_t>& exits) const {
    optional<ChainPosition> place = directory.chain_position(graph.get_id(entering));
    if (!place) {
        // A walk visits a free node and leaves through its other side.
        exits.push_back(entering);
        return false;
    }
    const Chain& chain = *place->chain;
    optional<endpoint_t> entry = entry_endpoint(chain, place->position,
                                                graph.get_is_reverse(entering) != place->is_reverse);
    if (!entry) {
        throw logic_error("error: an edge enters a chain between its bounds, so the snarl decomposition is broken");
    }
    const TraversalSet& traversals = survey.of_chain(&chain);
    handle_t start_out = graph.flip(chain_node(graph, chain, 0));
    handle_t end_out = chain_node(graph, chain, chain.size());
    if (traversals.contains(*entry, opposite(*entry))) {
        exits.push_back(*entry == START ? end_out : start_out);
    }
    if (traversals.contains(*entry, *entry)) {
        exits.push_back(*entry == START ? start_out : end_out);
    }
    return traversals.contains(*entry, TIP);
}

SnarlManagerDecomposition::TraversalSurvey::Surveyor::WalkEnds
SnarlManagerDecomposition::TraversalSurvey::Surveyor::walk_snarl(const Snarl* snarl, const vector<handle_t>& starts) const {
    const Chain& chain = *manager.chain_of(snarl);
    size_t rank = manager.chain_rank_of(snarl);
    handle_t start_out = graph.flip(chain_node(graph, chain, rank));
    handle_t end_out = chain_node(graph, chain, rank + 1);

    // Search breadth-first over the handles that a walk leaves through their
    // right sides.
    WalkEnds ends;
    unordered_set<handle_t> seen(starts.begin(), starts.end());
    deque<handle_t> to_leave(starts.begin(), starts.end());
    auto go = [&seen, &to_leave](const handle_t& exit) {
        if (seen.insert(exit).second) {
            to_leave.push_back(exit);
        }
    };
    while (!to_leave.empty()) {
        handle_t here = to_leave.front();
        to_leave.pop_front();
        bool has_edges = false;
        graph.follow_edges(here, false, [&](const handle_t& next) {
            has_edges = true;
            if (next == end_out) {
                ends.end = true;
            } else if (next == start_out) {
                ends.start = true;
            } else {
                vector<handle_t> exits;
                ends.tip |= enter_child(next, exits);
                for (const handle_t& exit : exits) {
                    go(exit);
                }
            }
        });
        // The right side of here is a dead end. Unless here is on a bound of
        // the snarl, that dead end is a tip of the snarl.
        nid_t id = graph.get_id(here);
        if (!has_edges && id != snarl->start().node_id() && id != snarl->end().node_id()) {
            ends.tip = true;
        }
    }
    return ends;
}

SnarlManagerDecomposition::TraversalSurvey::Surveyor::WalkEnds
SnarlManagerDecomposition::TraversalSurvey::Surveyor::walk_chain(const Chain& chain,
                                                       const vector<pair<size_t, bool>>& starts) const {
    // Search breadth-first over the chain nodes that a walk leaves, and the
    // way it goes along the chain when it leaves each.
    WalkEnds ends;
    set<pair<size_t, bool>> seen(starts.begin(), starts.end());
    deque<pair<size_t, bool>> to_leave(starts.begin(), starts.end());
    auto go = [&seen, &to_leave](size_t position, bool is_reverse) {
        if (seen.emplace(position, is_reverse).second) {
            to_leave.emplace_back(position, is_reverse);
        }
    };
    while (!to_leave.empty()) {
        auto [position, is_reverse] = to_leave.front();
        to_leave.pop_front();
        if (!is_reverse) {
            if (position == chain.size()) {
                ends.end = true;
                continue;
            }
            // Enter the snarl ahead at its start.
            const TraversalSet& traversals = survey.of_snarl(chain[position].first);
            if (traversals.contains(START, END)) {
                go(position + 1, false);
            }
            if (traversals.contains(START, START)) {
                go(position, true);
            }
            ends.tip |= traversals.contains(START, TIP);
        } else {
            if (position == 0) {
                ends.start = true;
                continue;
            }
            // Enter the snarl behind at its end.
            const TraversalSet& traversals = survey.of_snarl(chain[position - 1].first);
            if (traversals.contains(END, START)) {
                go(position - 1, true);
            }
            if (traversals.contains(END, END)) {
                go(position, false);
            }
            ends.tip |= traversals.contains(END, TIP);
        }
    }
    return ends;
}

void SnarlManagerDecomposition::TraversalSurvey::Surveyor::add_walks_from_bounds(TraversalSet& traversals,
                                                                                 const WalkEnds& from_start,
                                                                                 const WalkEnds& from_end) {
    // A walk taken backward is a walk, so a traversal from a tip is
    // realizable when the same traversal taken backward is.
    for (const auto& [from, ends] : {pair(START, from_start), pair(END, from_end)}) {
        if (ends.start) {
            traversals.add(from, START);
        }
        if (ends.end) {
            traversals.add(from, END);
        }
        if (ends.tip) {
            traversals.add(from, TIP);
            traversals.add(TIP, from);
        }
    }
}

SnarlManagerDecomposition::TraversalSet SnarlManagerDecomposition::TraversalSurvey::Surveyor::survey_snarl(const Snarl* snarl) const {
    const Chain& chain = *manager.chain_of(snarl);
    size_t rank = manager.chain_rank_of(snarl);
    handle_t start_in = chain_node(graph, chain, rank);
    handle_t end_in = graph.flip(chain_node(graph, chain, rank + 1));

    TraversalSet traversals;
    WalkEnds from_start = walk_snarl(snarl, {start_in});
    WalkEnds from_end = walk_snarl(snarl, {end_in});
    add_walks_from_bounds(traversals, from_start, from_end);

    // A walk from a tip starts at a dead end in the snarl. The dead end is a
    // side of a free node, the outer side of a bound of a child chain, or a
    // tip inside a child chain.
    vector<handle_t> from_tips;
    bool tip_to_tip = false;
    for (const nid_t& id : directory.free_nodes(snarl)) {
        for (bool is_reverse : {false, true}) {
            if (handle_t handle = graph.get_handle(id, is_reverse); left_is_dead_end(graph, handle)) {
                from_tips.push_back(handle);
            }
        }
    }
    for (const Chain& child : manager.chains_of(snarl)) {
        const TraversalSet& child_traversals = survey.of_chain(&child);
        handle_t child_start_in = chain_node(graph, child, 0);
        handle_t child_end_in = graph.flip(chain_node(graph, child, child.size()));
        tip_to_tip |= child_traversals.contains(TIP, TIP);
        if (child_traversals.contains(TIP, START)) {
            from_tips.push_back(graph.flip(child_start_in));
        }
        if (child_traversals.contains(TIP, END)) {
            from_tips.push_back(graph.flip(child_end_in));
        }
        // A walk from a dead end on the outer side of a bound enters the
        // child chain there.
        for (const handle_t& entering : {child_start_in, child_end_in}) {
            if (left_is_dead_end(graph, entering)) {
                tip_to_tip |= enter_child(entering, from_tips);
            }
        }
    }
    if (tip_to_tip || (!from_tips.empty() && walk_snarl(snarl, from_tips).tip)) {
        traversals.add(TIP, TIP);
    }
    return traversals;
}

SnarlManagerDecomposition::TraversalSet SnarlManagerDecomposition::TraversalSurvey::Surveyor::survey_chain(const Chain& chain) const {
    TraversalSet traversals;
    WalkEnds from_start = walk_chain(chain, {{0, false}});
    WalkEnds from_end = walk_chain(chain, {{chain.size(), true}});
    add_walks_from_bounds(traversals, from_start, from_end);

    // A walk from a tip starts in a snarl of the chain. It leaves that snarl
    // through its start or its end, or ends at a tip in the same snarl.
    vector<pair<size_t, bool>> from_tips;
    bool tip_to_tip = false;
    for (size_t rank = 0; rank < chain.size(); rank++) {
        const TraversalSet& snarl_traversals = survey.of_snarl(chain[rank].first);
        tip_to_tip |= snarl_traversals.contains(TIP, TIP);
        if (snarl_traversals.contains(TIP, START)) {
            from_tips.emplace_back(rank, true);
        }
        if (snarl_traversals.contains(TIP, END)) {
            from_tips.emplace_back(rank + 1, false);
        }
    }
    if (tip_to_tip || (!from_tips.empty() && walk_chain(chain, from_tips).tip)) {
        traversals.add(TIP, TIP);
    }
    return traversals;
}

////////////////////////////////////////////////////////////////////////////
// SnarlManagerDecomposition
////////////////////////////////////////////////////////////////////////////

SnarlManagerDecomposition::SnarlManagerDecomposition(const SnarlManager& manager, const HandleGraph& graph) :
    manager(manager), graph(graph), snarls_by_number(number_snarls(manager)), directory(manager, graph),
    survey(manager, graph, directory) {
}

handle_t SnarlManagerDecomposition::snarl_bound(const Snarl* snarl, bool get_end) const {
    return chain_node(graph, *manager.chain_of(snarl), manager.chain_rank_of(snarl) + (get_end ? 1 : 0));
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

net_handle_t SnarlManagerDecomposition::chain_node_net(const Chain& chain, size_t position, bool is_reverse) const {
    handle_t node = chain_node(graph, chain, position);
    return node_net(graph, is_reverse ? graph.flip(node) : node);
}

bool SnarlManagerDecomposition::for_each_snarl_child(const Snarl* parent,
                                                     const function<bool(const net_handle_t&)>& iteratee) const {
    for (const Chain& chain : manager.chains_of(parent)) {
        if (!iteratee(chain_net(chain))) {
            return false;
        }
    }
    for (const nid_t& id : directory.free_nodes(parent)) {
        if (!iteratee(pack(NetType::TRIVIAL_CHAIN, START, END, id))) {
            return false;
        }
    }
    return true;
}

bool SnarlManagerDecomposition::for_each_chain_child(const net_handle_t& traversal,
                                                     const function<bool(const net_handle_t&)>& iteratee) const {
    const Chain& chain = chain_of_net(traversal);
    // Number the children forward along the chain. The chain node at position
    // i is child 2i, and the snarl after it is child 2i + 1. A looping chain
    // has no child after its last snarl, since that would be its first chain
    // node again.
    size_t child_count = 2 * chain.size() + (is_looping(chain) ? 0 : 1);
    auto child = [this, &chain](size_t index, bool is_reverse) {
        return index % 2 == 0 ? chain_node_net(chain, index / 2, is_reverse)
                              : snarl_net(chain[index / 2].first, !is_reverse);
    };
    if (is_backward(traversal)) {
        // An end-to-start traversal of a looping chain starts at its first
        // chain node, which is also its end bound.
        size_t first = is_looping(chain) ? 0 : child_count - 1;
        if (!iteratee(child(first, true))) {
            return false;
        }
        for (size_t index = child_count - 1; index > 0; index--) {
            if (index != first && !iteratee(child(index, true))) {
                return false;
            }
        }
        return first == 0 || iteratee(child(0, true));
    }
    for (size_t index = 0; index < child_count; index++) {
        net_handle_t next = child(index, false);
        if (!is_forward(traversal)) {
            next = with_traversal(next, START, END);
        }
        if (!iteratee(next)) {
            return false;
        }
    }
    return true;
}

optional<net_handle_t> SnarlManagerDecomposition::next_in_chain(const net_handle_t& here) const {
    if (is_snarl(here)) {
        if (ends_at(here) == TIP) {
            return nullopt;
        }
        const Snarl* snarl = get_snarl(here);
        size_t rank = manager.chain_rank_of(snarl);
        // Leaving a snarl through its end, a walk next visits the chain node
        // after it, oriented forward along the chain. Leaving through its
        // start, the walk next visits the chain node before it, oriented in
        // reverse along the chain.
        bool is_reverse = ends_at(here) == START;
        return chain_node_net(*manager.chain_of(snarl), is_reverse ? rank : rank + 1, is_reverse);
    }

    handle_t handle = get_handle(here, &graph);
    optional<ChainPosition> place = directory.chain_position(graph.get_id(handle));
    if (!place) {
        // A free node is the only child of its trivial chain.
        return nullopt;
    }
    bool is_reverse = graph.get_is_reverse(handle) != place->is_reverse;
    optional<size_t> rank = snarl_ahead(*place->chain, place->position, is_reverse);
    if (!rank.has_value()) {
        return nullopt;
    }
    return snarl_net((*place->chain)[*rank].first, !is_reverse);
}

net_handle_t SnarlManagerDecomposition::sibling_entered_by(const handle_t& handle,
                                                                     const Snarl* parent) const {
    // A handle on a bound of the parent, oriented out of the parent, leaves
    // the parent.
    if (parent != nullptr) {
        if (handle == snarl_bound(parent, true)) {
            return get_bound(snarl_net(parent, true), true, false);
        }
        if (handle == graph.flip(snarl_bound(parent, false))) {
            return get_bound(snarl_net(parent, true), false, false);
        }
    }
    nid_t id = graph.get_id(handle);
    optional<ChainPosition> place = directory.chain_position(id);
    if (!place) {
        // The handle enters the trivial chain of a free node, going the way
        // the handle is oriented.
        return directed(pack(NetType::TRIVIAL_CHAIN, START, END, id), !graph.get_is_reverse(handle));
    }
    bool is_reverse = graph.get_is_reverse(handle) != place->is_reverse;
    optional<endpoint_t> entry = entry_endpoint(*place->chain, place->position, is_reverse);
    if (!entry) {
        throw logic_error("error: an edge enters a chain between its bounds, so the snarl decomposition is broken");
    }
    return directed(chain_net(*place->chain), *entry == START);
}

bool SnarlManagerDecomposition::for_each_sibling_ahead(const net_handle_t& here,
                                                       const function<bool(const net_handle_t&)>& iteratee) const {
    const Snarl* parent = nullptr;
    handle_t exit;
    if (is_sentinel(here)) {
        if (!faces_in(here)) {
            // A sentinel facing out leads out of the snarl.
            return true;
        }
        parent = snarls_by_number.at(payload_of(here));
        exit = get_handle(here, &graph);
    } else {
        if (ends_at(here) == TIP) {
            return true;
        }
        if (type_of(here) == NetType::CHAIN && is_looping(chain_of_net(here))) {
            // Every edge of the boundary node of a looping chain leads into a
            // snarl of the chain.
            return true;
        }
        net_handle_t up = get_parent(here);
        parent = is_root(up) ? nullptr : get_snarl(up);
        exit = get_handle(get_bound(here, ends_at(here) == END, false), &graph);
    }
    return graph.follow_edges(exit, false, [&](const handle_t& next) {
        return iteratee(sibling_entered_by(next, parent));
    });
}

bool SnarlManagerDecomposition::for_each_traversal_entered(const net_handle_t& entered,
                                                           const function<bool(const net_handle_t&)>& iteratee) const {
    if (!is_snarl(entered) && type_of(entered) != NetType::CHAIN) {
        // A node, a trivial chain or a sentinel has one traversal from each
        // place a walk can enter it.
        return iteratee(entered);
    }
    return traversals_of(entered).for_each([this, &entered, &iteratee](endpoint_t start, endpoint_t end) {
        return start != starts_at(entered) || iteratee(with_traversal(entered, start, end));
    });
}

bool SnarlManagerDecomposition::follow_right(const net_handle_t& here,
                                             const function<bool(const net_handle_t&)>& iteratee) const {
    auto each_traversal = [this, &iteratee](const net_handle_t& entered) {
        return for_each_traversal_entered(entered, iteratee);
    };
    using enum NetType;
    switch (type_of(here)) {
    case ROOT:
        return true;
    case NODE:
    case SNARL:
        {
            optional<net_handle_t> next = next_in_chain(here);
            return !next || each_traversal(*next);
        }
    default:
        return for_each_sibling_ahead(here, each_traversal);
    }
}

net_handle_t SnarlManagerDecomposition::get_snarl_net(const Snarl* snarl) const {
    return snarl_net(snarl, true);
}

const Snarl* SnarlManagerDecomposition::get_snarl(const net_handle_t& net) const {
    if (type_of(net) != NetType::SNARL) {
        throw invalid_argument("error: net handle is not for a snarl");
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
        return handle_graph->get_handle(payload_of(net), ends_at(net) == START);
    case SENTINEL:
        {
            const Snarl* snarl = snarls_by_number.at(payload_of(net));
            // We have the start bound oriented into the snarl, or the end
            // bound oriented out of it. Flip it if the sentinel faces the
            // other way.
            handle_t bound = snarl_bound(snarl, starts_at(net) == END);
            if (bool oriented_in = starts_at(net) == START; faces_in(net) != oriented_in) {
                bound = graph.flip(bound);
            }
            return handle_graph->get_handle(graph.get_id(bound), graph.get_is_reverse(bound));
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
        return pack(SNARL, START, END, payload_of(child));
    case SNARL:
        {
            // A snarl is oriented forward along its chain, so the parent
            // traversal of the chain goes the same way as the snarl traversal.
            const Chain& chain = *manager.chain_of(get_snarl(child));
            return directed(chain_net(chain), !is_backward(child));
        }
    case NODE:
        {
            optional<ChainPosition> place = directory.chain_position(payload_of(child));
            if (!place) {
                // A free node is the only child of its trivial chain, which
                // goes the way the node traversal goes.
                return pack(TRIVIAL_CHAIN, starts_at(child), ends_at(child), payload_of(child));
            }
            if (!is_forward(child) && !is_backward(child)) {
                // The interface gives the parent start-to-end when the child
                // has no direction along it.
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
            // Oriented forward along the chain, the start boundary node faces
            // into the chain and the end boundary node faces out of it.
            handle_t node = chain_node(graph, chain, get_end ? chain.size() : 0);
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
            static const array<pair<endpoint_t, endpoint_t>, 6> preferred {{
                {START, END}, {START, TIP}, {END, TIP}, {START, START}, {END, END}, {TIP, TIP}
            }};
            TraversalSet realizable = traversals_of(net);
            for (const auto& [start, end] : preferred) {
                if (realizable.contains(start, end)) {
                    return with_traversal(net, start, end);
                }
            }
            // TODO: No traversal of this net graph element is realizable, but
            // the interface wants the canonical traversal to be realizable.
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

pair<net_handle_t, endpoint_t> SnarlManagerDecomposition::parent_endpoint(const net_handle_t& net, bool is_start) const {
    using enum NetType;
    switch (type_of(net)) {
    case SENTINEL:
        if (faces_in(net) != is_start) {
            throw invalid_argument(is_start ? "error: a traversal cannot start at a sentinel facing out"
                                            : "error: a traversal cannot end at a sentinel facing in");
        }
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
                // A traversal of a free node enters its trivial chain where it
                // starts, and leaves it where it ends.
                return {pack(TRIVIAL_CHAIN, START, END, payload_of(net)), is_start ? starts_at(net) : ends_at(net)};
            }
            const Chain& chain = *place->chain;
            bool is_reverse = graph.get_is_reverse(get_handle(net, &graph)) != place->is_reverse;
            optional<endpoint_t> endpoint = is_start ? entry_endpoint(chain, place->position, is_reverse)
                                                     : exit_endpoint(chain, place->position, is_reverse);
            if (!endpoint) {
                throw invalid_argument("error: node traversal does not enter or leave its chain at a bound");
            }
            return {chain_net(chain), *endpoint};
        }
    default:
        throw invalid_argument("error: the root cannot start or end a traversal");
    }
}

net_handle_t SnarlManagerDecomposition::get_parent_traversal(const net_handle_t& traversal_start,
                                                             const net_handle_t& traversal_end) const {
    auto [start_parent, start] = parent_endpoint(traversal_start, true);
    auto [end_parent, end] = parent_endpoint(traversal_end, false);
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
        return for_each_snarl_child(nullptr, iteratee);
    case SNARL:
        return for_each_snarl_child(get_snarl(traversal), iteratee);
    case CHAIN:
        return for_each_chain_child(traversal, iteratee);
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
        traversals.add(START, END);
        traversals.add(END, START);
        break;
    case SNARL:
        traversals = survey.of_snarl(get_snarl(item));
        break;
    case CHAIN:
        traversals = survey.of_chain(&chain_of_net(item));
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
    if (!go_left) {
        return follow_right(here, iteratee);
    }
    // In a bidirected graph, a walk taken backward is a walk. So the
    // traversals from which going right gives here are the flips of the
    // traversals that going right from the flip of here gives.
    return follow_right(flip(here), [this, &iteratee](const net_handle_t& next) {
        return iteratee(flip(next));
    });
}

bool SnarlManagerDecomposition::for_each_tippy_child_impl(const net_handle_t& parent,
                                                          const function<bool(const net_handle_t&)>& iteratee) const {
    // Give each realizable traversal of a child that starts at an endpoint.
    auto each_from = [&](const net_handle_t& child, endpoint_t start) {
        return traversals_of(child).for_each([&](endpoint_t from, endpoint_t to) {
            return from != start || iteratee(with_traversal(child, from, to));
        });
    };
    using enum NetType;
    switch (type_of(parent)) {
    case ROOT:
    case SNARL:
    {
        // A tip of the parent is a tip in a child chain, the outer side of a
        // boundary node of a child chain with no edges on it, or a dead end on
        // a free node.
        const Snarl* snarl = is_root(parent) ? nullptr : get_snarl(parent);
        for (const Chain& chain : manager.chains_of(snarl)) {
            net_handle_t child = chain_net(chain);
            if (!each_from(child, TIP)) {
                return false;
            }
            if (left_is_dead_end(graph, chain_node(graph, chain, 0)) && !each_from(child, START)) {
                return false;
            }
            if (left_is_dead_end(graph, graph.flip(chain_node(graph, chain, chain.size()))) &&
                !each_from(child, END)) {
                return false;
            }
        }
        for (const nid_t& id : directory.free_nodes(snarl)) {
            for (bool is_reverse : {false, true}) {
                if (left_is_dead_end(graph, graph.get_handle(id, is_reverse)) &&
                    !iteratee(directed(pack(TRIVIAL_CHAIN, START, END, id), !is_reverse))) {
                    return false;
                }
            }
        }
        return true;
    }
    case CHAIN:
        for (const auto& [snarl, backward] : chain_of_net(parent)) {
            if (!each_from(snarl_net(snarl, true), TIP)) {
                return false;
            }
        }
        return true;
    default:
        return true;
    }
}

}
