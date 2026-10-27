/** \file
 * Implements EditSnarlDecompositionOverlay.
 */

#include "edit_snarl_decomposition_overlay.hpp"

#include <stdexcept>
#include <string>

#include <handlegraph/util.hpp>

namespace vg {

using namespace std;

namespace {

/// The lowest bit of an overlay net handle, set for an edit.
constexpr uint64_t EDIT_BIT = 1;
/// The bit above it in an edit's net handle, set when the edit is read end to start.
constexpr uint64_t REVERSE_BIT = 2;
/// How far an edit's number is shifted up in its net handle.
constexpr uint64_t EDIT_NUMBER_SHIFT = 2;

}

EditSnarlDecompositionOverlay::EditSnarlDecompositionOverlay(const SnarlDecomposition& inner,
                                                             const HandleGraph& graph) :
    inner(inner), graph(graph) {
    // Nothing to do
}

net_handle_t EditSnarlDecompositionOverlay::from_inner(const net_handle_t& inner_net) {
    const uint64_t value = handlegraph::as_integer(inner_net);
    if (value >> 63 != 0) {
        throw std::logic_error("EditSnarlDecompositionOverlay: an inner net handle has its top bit "
                               "set, so it cannot be shifted up");
    }
    return handlegraph::as_net_handle(value << 1);
}

net_handle_t EditSnarlDecompositionOverlay::to_inner(const net_handle_t& net) {
    return handlegraph::as_net_handle(handlegraph::as_integer(net) >> 1);
}

net_handle_t EditSnarlDecompositionOverlay::edit_net(size_t number, bool reverse) {
    return handlegraph::as_net_handle(((uint64_t)number << EDIT_NUMBER_SHIFT)
                                      | (reverse ? REVERSE_BIT : 0) | EDIT_BIT);
}

size_t EditSnarlDecompositionOverlay::edit_number(const net_handle_t& net) {
    return handlegraph::as_integer(net) >> EDIT_NUMBER_SHIFT;
}

bool EditSnarlDecompositionOverlay::edit_is_reverse(const net_handle_t& net) {
    return (handlegraph::as_integer(net) & REVERSE_BIT) != 0;
}

EditSnarlDecompositionOverlay::Edit
EditSnarlDecompositionOverlay::edit_of(const net_handle_t& net) const {
    std::lock_guard<std::mutex> lock(edits_mutex);
    return edits.at(edit_number(net));
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Edits
////////////////////////////////////////////////////////////////////////////////////////////////////

net_handle_t EditSnarlDecompositionOverlay::lowest_holder(const handle_t& a,
                                                          const handle_t& b) const {
    // The canonical inner net handles of the elements holding a node, from the node up to the root.
    auto ancestors = [&](const handle_t& node) {
        vector<net_handle_t> out;
        net_handle_t here = inner.canonical(inner.get_net(graph.forward(node), &graph));
        out.push_back(here);
        while (!inner.is_root(here)) {
            here = inner.canonical(inner.get_parent(here));
            out.push_back(here);
        }
        return out;
    };
    const vector<net_handle_t> from_a = ancestors(a);
    const vector<net_handle_t> from_b = ancestors(b);
    // The two lists end at the root. Walk back from there while they agree; the last element they
    // share is the lowest one holding both nodes.
    size_t shared = 0;
    while (shared < from_a.size() && shared < from_b.size()
           && from_a[from_a.size() - 1 - shared] == from_b[from_b.size() - 1 - shared]) {
        ++shared;
    }
    return from_inner(from_a[from_a.size() - shared]);
}

net_handle_t EditSnarlDecompositionOverlay::add_edit(const handle_t& start, size_t start_offset,
                                                     const handle_t& end, size_t end_offset) {
    if (start_offset > graph.get_length(start) || end_offset > graph.get_length(end)) {
        throw std::invalid_argument("EditSnarlDecompositionOverlay: an edit's offset is past the "
                                    "end of its node");
    }
    Edit edit{start, start_offset, end, end_offset, lowest_holder(start, end)};
    std::lock_guard<std::mutex> lock(edits_mutex);
    size_t number;
    if (!free_numbers.empty()) {
        number = free_numbers.back();
        free_numbers.pop_back();
        edits[number] = edit;
    } else {
        number = edits.size();
        edits.push_back(edit);
    }
    children_of[edit.parent].push_back(number);
    return edit_net(number, false);
}

void EditSnarlDecompositionOverlay::remove_edits(const net_handle_t& element) {
    // Find every element at or below `element` that has edits, by walking up from each.
    std::lock_guard<std::mutex> lock(edits_mutex);
    const net_handle_t top = canonical(element);
    vector<net_handle_t> emptied;
    for (const auto& [parent, numbers] : children_of) {
        net_handle_t here = parent;
        while (true) {
            if (here == top) {
                emptied.push_back(parent);
                break;
            }
            if (is_root(here)) {
                break;
            }
            here = canonical(from_inner(inner.get_parent(to_inner(here))));
        }
    }
    for (const net_handle_t& parent : emptied) {
        for (size_t number : children_of[parent]) {
            free_numbers.push_back(number);
        }
        children_of.erase(parent);
    }
}

bool EditSnarlDecompositionOverlay::is_edit(const net_handle_t& net) const {
    return (handlegraph::as_integer(net) & EDIT_BIT) != 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Bounds of sites, for snarls and edits
////////////////////////////////////////////////////////////////////////////////////////////////////

handle_t EditSnarlDecompositionOverlay::get_site_start(const net_handle_t& site) const {
    if (is_edit(site)) {
        const Edit edit = edit_of(site);
        return edit_is_reverse(site) ? graph.flip(edit.end) : edit.start;
    }
    // get_bound() ignores the traversal, so read the bound the traversal starts from.
    const net_handle_t inner_site = to_inner(site);
    if (inner.starts_at(inner_site) == END) {
        return graph.flip(inner.get_handle(inner.get_bound(inner_site, true, false), &graph));
    }
    return inner.get_handle(inner.get_bound(inner_site, false, true), &graph);
}

handle_t EditSnarlDecompositionOverlay::get_site_end(const net_handle_t& site) const {
    if (is_edit(site)) {
        const Edit edit = edit_of(site);
        return edit_is_reverse(site) ? graph.flip(edit.start) : edit.end;
    }
    const net_handle_t inner_site = to_inner(site);
    if (inner.starts_at(inner_site) == END) {
        return graph.flip(inner.get_handle(inner.get_bound(inner_site, false, true), &graph));
    }
    return inner.get_handle(inner.get_bound(inner_site, true, false), &graph);
}

size_t EditSnarlDecompositionOverlay::get_site_start_offset(const net_handle_t& site) const {
    if (is_edit(site)) {
        const Edit edit = edit_of(site);
        // Read the other way, the end boundary becomes the start, and the bases before it along
        // the flipped node are the ones after it along the node.
        return edit_is_reverse(site) ? graph.get_length(edit.end) - edit.end_offset
                                     : edit.start_offset;
    }
    return graph.get_length(get_site_start(site));
}

size_t EditSnarlDecompositionOverlay::get_site_end_offset(const net_handle_t& site) const {
    if (is_edit(site)) {
        const Edit edit = edit_of(site);
        return edit_is_reverse(site) ? graph.get_length(edit.start) - edit.start_offset
                                     : edit.end_offset;
    }
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// SnarlDecomposition, passed through to the inner decomposition except for edits
////////////////////////////////////////////////////////////////////////////////////////////////////

net_handle_t EditSnarlDecompositionOverlay::get_root() const {
    return from_inner(inner.get_root());
}

bool EditSnarlDecompositionOverlay::is_root(const net_handle_t& net) const {
    return !is_edit(net) && inner.is_root(to_inner(net));
}

bool EditSnarlDecompositionOverlay::is_snarl(const net_handle_t& net) const {
    return !is_edit(net) && inner.is_snarl(to_inner(net));
}

bool EditSnarlDecompositionOverlay::is_chain(const net_handle_t& net) const {
    return !is_edit(net) && inner.is_chain(to_inner(net));
}

bool EditSnarlDecompositionOverlay::is_node(const net_handle_t& net) const {
    return !is_edit(net) && inner.is_node(to_inner(net));
}

bool EditSnarlDecompositionOverlay::is_sentinel(const net_handle_t& net) const {
    return !is_edit(net) && inner.is_sentinel(to_inner(net));
}

net_handle_t EditSnarlDecompositionOverlay::get_net(const handle_t& handle,
                                                    const HandleGraph* graph) const {
    return from_inner(inner.get_net(handle, graph));
}

handle_t EditSnarlDecompositionOverlay::get_handle(const net_handle_t& net,
                                                   const HandleGraph* graph) const {
    if (is_edit(net)) {
        throw std::logic_error("EditSnarlDecompositionOverlay: get_handle on an edit");
    }
    return inner.get_handle(to_inner(net), graph);
}

net_handle_t EditSnarlDecompositionOverlay::get_parent(const net_handle_t& child) const {
    if (is_edit(child)) {
        return edit_of(child).parent;
    }
    return from_inner(inner.get_parent(to_inner(child)));
}

net_handle_t EditSnarlDecompositionOverlay::get_bound(const net_handle_t& snarl, bool get_end,
                                                      bool face_in) const {
    if (is_edit(snarl)) {
        throw std::logic_error("EditSnarlDecompositionOverlay: get_bound on an edit");
    }
    return from_inner(inner.get_bound(to_inner(snarl), get_end, face_in));
}

net_handle_t EditSnarlDecompositionOverlay::flip(const net_handle_t& net) const {
    if (is_edit(net)) {
        return edit_net(edit_number(net), !edit_is_reverse(net));
    }
    return from_inner(inner.flip(to_inner(net)));
}

net_handle_t EditSnarlDecompositionOverlay::canonical(const net_handle_t& net) const {
    if (is_edit(net)) {
        return edit_net(edit_number(net), false);
    }
    return from_inner(inner.canonical(to_inner(net)));
}

SnarlDecomposition::endpoint_t
EditSnarlDecompositionOverlay::starts_at(const net_handle_t& traversal) const {
    if (is_edit(traversal)) {
        return edit_is_reverse(traversal) ? END : START;
    }
    return inner.starts_at(to_inner(traversal));
}

SnarlDecomposition::endpoint_t
EditSnarlDecompositionOverlay::ends_at(const net_handle_t& traversal) const {
    if (is_edit(traversal)) {
        return edit_is_reverse(traversal) ? START : END;
    }
    return inner.ends_at(to_inner(traversal));
}

net_handle_t EditSnarlDecompositionOverlay::get_parent_traversal(
    const net_handle_t& traversal_start, const net_handle_t& traversal_end) const {
    if (is_edit(traversal_start) || is_edit(traversal_end)) {
        throw std::logic_error("EditSnarlDecompositionOverlay: get_parent_traversal on an edit");
    }
    return from_inner(inner.get_parent_traversal(to_inner(traversal_start), to_inner(traversal_end)));
}

bool EditSnarlDecompositionOverlay::for_each_child_impl(
    const net_handle_t& traversal, const std::function<bool(const net_handle_t&)>& iteratee) const {
    if (is_edit(traversal)) {
        // An edit has no children.
        return true;
    }
    // A node has no children in the inner decomposition, which may refuse to be asked.
    if (!is_node(traversal)
        && !inner.for_each_child(to_inner(traversal), [&](const net_handle_t& child) {
               return iteratee(from_inner(child));
           })) {
        return false;
    }
    // Copy the numbers out, so that the iteratee can add edits without the lock held.
    vector<size_t> numbers;
    {
        std::lock_guard<std::mutex> lock(edits_mutex);
        auto found = children_of.find(canonical(traversal));
        if (found != children_of.end()) {
            numbers = found->second;
        }
    }
    for (size_t number : numbers) {
        if (!iteratee(edit_net(number, false))) {
            return false;
        }
    }
    return true;
}

bool EditSnarlDecompositionOverlay::for_each_traversal_impl(
    const net_handle_t& item, const std::function<bool(const net_handle_t&)>& iteratee) const {
    if (is_edit(item)) {
        return iteratee(edit_net(edit_number(item), false))
               && iteratee(edit_net(edit_number(item), true));
    }
    return inner.for_each_traversal(to_inner(item), [&](const net_handle_t& traversal) {
        return iteratee(from_inner(traversal));
    });
}

bool EditSnarlDecompositionOverlay::follow_net_edges_impl(
    const net_handle_t& here, const HandleGraph* graph, bool go_left,
    const std::function<bool(const net_handle_t&)>& iteratee) const {
    if (is_edit(here)) {
        // An edit is not a step of any chain, so no net edges lead from it.
        return true;
    }
    return inner.follow_net_edges(to_inner(here), graph, go_left, [&](const net_handle_t& next) {
        return iteratee(from_inner(next));
    });
}

}
