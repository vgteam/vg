/**
 * \file incremental_subgraph.cpp: contains the implementation of IncrementalSubgraph
 */


#include "incremental_subgraph.hpp"

//#define debug_incremental_subgraph

namespace vg {

using namespace std;

IncrementalSubgraph::IncrementalSubgraph(const HandleGraph& graph,
                                         const pos_t& start_pos,
                                         bool extract_left,
                                         int64_t max_distance) :
    graph(&graph), extract_left(extract_left), max_distance(max_distance)
{
    handle_t start = graph.get_handle(id(start_pos), is_rev(start_pos));
    int64_t dist = extract_left ? offset(start_pos) - graph.get_length(start) : -offset(start_pos);
    extracted.emplace_back(start, vector<size_t>(), vector<size_t>(), dist, dist);
    extracted_index[start] = 0;
    
    // initialize the frontier
    int64_t dist_thru = dist + graph.get_length(start);
    if (dist_thru < max_distance) {
        graph.follow_edges(start, extract_left, [&](const handle_t& next) {
            // we start by covering one edge (the one we're traversing right now)
            size_t back_degree = graph.get_degree(next, !extract_left) - 1;
            auto entry = frontier.emplace(back_degree, dist_thru, next);
            frontier_index[next] = entry.first;
        });
    }
        
}

bool IncrementalSubgraph::is_extendable() const {
    return !frontier.empty();
}

handle_t IncrementalSubgraph::extend() {

    // get frontier group with fewest uncovered edges (breaking ties
    auto front = frontier.begin();
    auto nearest = *front;
    
#ifdef debug_incremental_subgraph
    cerr << "adding frontier node " << graph->get_id(get<2>(nearest)) << " " << graph->get_is_reverse(get<2>(nearest)) << endl;
#endif
    
    // remove the node from corresponding index
    frontier_index.erase(get<2>(nearest));
    frontier.erase(nearest);
    
    // the distance to the end of this node
    int64_t dist_thru = get<1>(nearest) + graph->get_length(get<2>(nearest));
    
    graph->follow_edges(get<2>(nearest), extract_left, [&](const handle_t& next) {
        auto it = frontier_index.find(next);
        if (it != frontier_index.end()) {
            // we've already added this node into the frontier, but now we're covering
            // one more of its edges is covered, move it into the next bucket
            auto record = *it->second;
            frontier.erase(it->second);
            auto entry = frontier.emplace(get<0>(record) - 1,
                                          min(dist_thru, get<1>(record)),
                                          get<2>(record));
            it->second = entry.first;
#ifdef debug_incremental_subgraph
            cerr << "decrementing incoming count of frontier node " << graph->get_id(next) << " " << graph->get_is_reverse(next) << " to " << get<0>(record) - 1 << endl;
#endif
        }
        else if (dist_thru < max_distance) {
            // we haven't added this node to the frontier, but it's still within the
            // distance limits
#ifdef debug_incremental_subgraph
            cerr << "adding " << graph->get_id(next) << " " << graph->get_is_reverse(next) << " to the frontier" << endl;
#endif
            
            // we start by covering one edge (the one we're traversing right now)
            size_t back_degree = graph->get_degree(next, !extract_left) - 1;
            auto entry = frontier.emplace(back_degree, dist_thru, next);
            frontier_index[next] = entry.first;
        }
    });
    
    // make an entry in the extracted
    extracted_index[get<2>(nearest)] = extracted.size();
    extracted.emplace_back();
    auto& extracted_record = extracted.back();
    get<0>(extracted_record) = get<2>(nearest);
    get<3>(extracted_record) = get<1>(nearest);
    get<4>(extracted_record) = 0;
    // add edges to predecessors that have already been taken out of the frontier
    graph->follow_edges(get<2>(nearest), !extract_left, [&](const handle_t& prev) {
        // TODO: is it sufficient to only add edges to the last copy in the presence
        // of cycles?
        auto it = extracted_index.find(prev);
        if (it != extracted_index.end()) {
            get<1>(extracted_record).push_back(it->second);
            get<2>(extracted[it->second]).push_back(extracted.size() - 1);
            get<4>(extracted_record) = max<int64_t>(get<4>(extracted_record),
                                                    get<4>(extracted[it->second])
                                                    + graph->get_length(it->first));
#ifdef debug_incremental_subgraph
            cerr << "add edge " << it->second << " -> " << extracted.size() - 1 << endl;
            cerr << "update max dist to " << get<4>(extracted_record) << endl;
#endif
        }
    });
    return get_handle(extracted.size(), false);
}

size_t IncrementalSubgraph::order_of(const handle_t& handle) const {
    return handlegraph::number_bool_packing::unpack_number(handle);
}

handle_t IncrementalSubgraph::handle_at_order(size_t i) const {
    return handlegraph::number_bool_packing::pack(i, false);
}

int64_t IncrementalSubgraph::min_distance_from_start(const handle_t& handle) const {
    return get<3>(extracted[order_of(handle)]);
}

int64_t IncrementalSubgraph::max_distance_from_start(const handle_t& handle) const {
    return get<4>(extracted[order_of(handle)]);
}

bool IncrementalSubgraph::extracting_left() const {
    return extract_left;
}

bool IncrementalSubgraph::has_node(id_t node_id) const {
    return node_id > 0 && node_id <= extracted.size();
}

handle_t IncrementalSubgraph::get_handle(const id_t& node_id, bool is_reverse) const {
    return handlegraph::number_bool_packing::pack(node_id - 1, is_reverse);
}

id_t IncrementalSubgraph::get_id(const handle_t& handle) const {
    return order_of(handle) + 1;
}

bool IncrementalSubgraph::get_is_reverse(const handle_t& handle) const {
    return handlegraph::number_bool_packing::unpack_bit(handle);
}

handle_t IncrementalSubgraph::flip(const handle_t& handle) const {
    return handlegraph::number_bool_packing::toggle_bit(handle);
}

size_t IncrementalSubgraph::get_length(const handle_t& handle) const {
    return graph->get_length(get<0>(extracted[order_of(handle)]));
}

string IncrementalSubgraph::get_sequence(const handle_t& handle) const {
    return graph->get_sequence(get_underlying_handle(handle));
}

bool IncrementalSubgraph::follow_edges_impl(const handle_t& handle, bool go_left,
                                            const function<bool(const handle_t&)>& iteratee) const {
    bool left_edges = (go_left != get_is_reverse(handle)) != extract_left;
    auto& edges = left_edges ? get<1>(extracted[order_of(handle)]) : get<2>(extracted[order_of(handle)]);
    bool keep_going = true;
    for (size_t i = 0; i < edges.size() && keep_going; ++i) {
        keep_going = iteratee(handlegraph::number_bool_packing::pack(edges[i],
                                                                     get_is_reverse(handle)));
    }
    return keep_going;
}

bool IncrementalSubgraph::for_each_handle_impl(const function<bool(const handle_t&)>& iteratee,
                                               bool parallel) const {
    bool keep_going = true;
    for (size_t i = 0; i < extracted.size() && keep_going; ++i) {
        keep_going = iteratee(handlegraph::number_bool_packing::pack(i, false));
    }
    // not doing parallel, never expect to use it
    return keep_going;
}

size_t IncrementalSubgraph::get_node_count() const {
    return extracted.size();
}

id_t IncrementalSubgraph::min_node_id() const {
    return 1;
}

id_t IncrementalSubgraph::max_node_id() const {
    return extracted.size();
}

size_t IncrementalSubgraph::get_degree(const handle_t& handle, bool go_left) const {
    bool left_edges = (go_left != get_is_reverse(handle)) != extract_left;
    return (left_edges ? get<1>(extracted[order_of(handle)]) : get<2>(extracted[order_of(handle)])).size();
}

size_t IncrementalSubgraph::get_edge_count() const {
    size_t count = 0;
    for (const auto& record : extracted) {
        count += get<1>(record).size();
    }
    return count;
}

char IncrementalSubgraph::get_base(const handle_t& handle, size_t index) const {
    return graph->get_base(get_underlying_handle(handle), index);
}

string IncrementalSubgraph::get_subsequence(const handle_t& handle, size_t index, size_t size) const {
    return graph->get_subsequence(get_underlying_handle(handle), index, size);
}

handle_t IncrementalSubgraph::get_underlying_handle(const handle_t& handle) const {
    auto underlying = get<0>(extracted[order_of(handle)]);
    return get_is_reverse(handle) ? graph->flip(underlying) : underlying;
}
}

