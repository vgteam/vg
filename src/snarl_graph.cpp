#include <atomic>
#include <charconv>
#include <chrono>
#include <cstdio>
#include <limits>

#include <omp.h>

#include "snarl_graph.hpp"
#include "symbolic_allele.hpp"
#include "read_likelihood_caller.hpp"
#include "algorithms/expand_context.hpp"
#include "annotation.hpp"
#include "gref.hpp"
#include "traversal_clusters.hpp"
#include "utility.hpp"

//#define debug

namespace vg {


SnarlGraph::SnarlGraph(const HandleGraph* backing_graph, SnarlManager& snarl_manager, vector<const Snarl*> snarls) :
    backing_graph(backing_graph),
    snarl_manager(snarl_manager) {
    for (const Snarl* snarl : snarls) {
        if (!snarl_manager.is_trivial(snarl, *backing_graph)) {
            this->snarls[backing_graph->get_handle(snarl->start().node_id(), snarl->start().backward())] =
                make_pair(backing_graph->get_handle(snarl->end().node_id(), snarl->end().backward()), true);
            this->snarls[backing_graph->get_handle(snarl->end().node_id(), !snarl->end().backward())] =
                make_pair(backing_graph->get_handle(snarl->start().node_id(), !snarl->start().backward()), false);
        }
    }
}

pair<bool, handle_t> SnarlGraph::node_to_snarl(handle_t handle) const {
    auto i = snarls.find(handle);
    if (i != snarls.end()) {
        return make_pair(true, i->second.first);
    } else {
        return make_pair(false, handle);
    }
}

tuple<bool, handle_t, edge_t> SnarlGraph::edge_to_snarl_edge(edge_t edge) const {
    auto i = snarls.find(edge.first);
    edge_t out_edge;
    handle_t out_node;
    bool out_found = false;
    if (i != snarls.end()) {
        // edge is from snarl start to after snarl end
        out_edge.first = i->second.first;
        out_edge.second = edge.second;
        out_node = edge.first;
        out_found = true;
    } else {
        // reverse of above
        i = snarls.find(backing_graph->flip(edge.second));
        if (i != snarls.end()) {
            out_edge.first = edge.first;
            out_edge.second = backing_graph->flip(i->second.first);
            out_node = edge.second;
            out_found = true;
        }
    }
    // note that we only have those two cases since our internal map contains
    // both orientations of the snarl.

    return make_tuple(out_found, out_node, out_edge);
}

void SnarlGraph::embed_snarl(Visit& visit) {
    handle_t handle = backing_graph->get_handle(visit.node_id(), visit.backward());
    auto it = snarls.find(handle);
    if (it != snarls.end()) {
        // edit the Visit in place to replace id, with the full snarl
        Snarl* snarl = visit.mutable_snarl();
        snarl->mutable_start()->set_node_id(visit.node_id());
        snarl->mutable_start()->set_backward(visit.backward());
        handle_t other = it->second.first;
        snarl->mutable_end()->set_node_id(backing_graph->get_id(other));
        snarl->mutable_end()->set_backward(backing_graph->get_is_reverse(other));
        if (it->second.second == false) {
            // put the snarl in an orientation consisten with other indexes
            swap(*snarl->mutable_start(), *snarl->mutable_end());
            snarl->mutable_start()->set_backward(!snarl->start().backward());
            snarl->mutable_end()->set_backward(!snarl->end().backward());
        }
        visit.set_node_id(0);
    }    
}

void SnarlGraph::embed_snarls(SnarlTraversal& traversal) {
    for (size_t i = 0; i < traversal.visit_size(); ++i) {
        Visit& visit = *traversal.mutable_visit(i);
        if (visit.node_id() > 0) { 
            embed_snarl(visit);
        }
    }
}

void SnarlGraph::embed_ref_path_snarls(SnarlTraversal& traversal) {
    vector<Visit> out_trav;
    size_t snarl_count = 0;
    bool in_snarl = false;
    handle_t snarl_end;
    for (size_t i = 0; i < traversal.visit_size(); ++i) {
        Visit& visit = *traversal.mutable_visit(i);
        handle_t handle = backing_graph->get_handle(visit.node_id(), visit.backward());
        if (in_snarl) {
            // nothing to do if we're in a snarl except check for the end and come out
            if (handle == snarl_end) {
                in_snarl = false;
            } 
        } else {
            // if we're not in a snarl, check for a new one
            auto it = snarls.find(handle);
            if (it != snarls.end()) {
                embed_snarl(visit);
                snarl_end = it->second.first;
                in_snarl = true;
                ++snarl_count;
            }
            out_trav.push_back(visit);
        }
    }

    // switch in the updated traversal
    if (snarl_count > 0) {
        traversal.clear_visit();
        for (Visit& visit : out_trav) {
            *traversal.add_visit() = visit;
        }
    }
}

bool SnarlGraph::follow_edges_impl(const handle_t& handle, bool go_left, const std::function<bool(const handle_t&)>& iteratee) const {
    if (!go_left) {        
        auto i = snarls.find(handle);
        if (i == snarls.end()) {
            return backing_graph->follow_edges(handle, go_left, iteratee);
        } else {
            return backing_graph->follow_edges(i->second.first, go_left, iteratee);
        }
    } else {
        return this->follow_edges_impl(backing_graph->flip(handle), !go_left, iteratee);
    }
}

// a lot of these don't strictly make sense.  ex, we would want has_node to
// hide stuff inside snarls.  but... we don't want to pay the cost of maintining
// structures for functions that aren't used..
bool SnarlGraph::has_node(nid_t node_id) const {
    return backing_graph->has_node(node_id);
}
handle_t SnarlGraph::get_handle(const nid_t& node_id, bool is_reverse) const {
    return backing_graph->get_handle(node_id, is_reverse);
}
nid_t SnarlGraph::get_id(const handle_t& handle) const {
    return backing_graph->get_id(handle);
}
bool SnarlGraph::get_is_reverse(const handle_t& handle) const {
    return backing_graph->get_is_reverse(handle);
}
handle_t SnarlGraph::flip(const handle_t& handle) const {
    return backing_graph->flip(handle);
}
size_t SnarlGraph::get_length(const handle_t& handle) const {
    return backing_graph->get_length(handle);
}
std::string SnarlGraph::get_sequence(const handle_t& handle) const {
    return backing_graph->get_sequence(handle);
}
size_t SnarlGraph::get_node_count() const {
    return backing_graph->get_node_count();
}
nid_t SnarlGraph::min_node_id() const {
    return backing_graph->min_node_id();
}
nid_t SnarlGraph::max_node_id() const {
    return backing_graph->max_node_id();
}
bool SnarlGraph::for_each_handle_impl(const std::function<bool(const handle_t&)>& iteratee, bool parallel) const {
    return backing_graph->for_each_handle(iteratee, parallel);
}

}

