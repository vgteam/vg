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
#ifdef debug_incremental_subgraph
    cerr << "initializing incremental graph from " << start_pos << " in " << (extract_left ? "left" : "right") << " direction up to distance " << max_distance << endl;
#endif
    
    handle_t start = graph.get_handle(id(start_pos), is_rev(start_pos));
    int64_t dist = extract_left ? offset(start_pos) - graph.get_length(start) : -offset(start_pos);
    
    extracted.emplace_back(start, vector<size_t>(), vector<size_t>(), dist, dist);
    
    // initialize the frontier
    int64_t dist_thru = dist + graph.get_length(start);
    if (dist_thru < max_distance) {
        graph.follow_edges(start, extract_left, [&](const handle_t& next) {
            
            // add all the back edges except the one we're using
            auto unseen_back_edges = new unordered_set<handle_t>();
            graph.follow_edges(next, !extract_left, [&](const handle_t& prev) {
                if (prev != start) {
                    // mark all edges unseen except the one we're traversing
                    unseen_back_edges->emplace(prev);
                }
            });
            auto seen_back_edges = new vector<size_t>(1, 0);
            // add the frontier and its random access index
            auto entry = frontier.emplace(dist_thru, next, unseen_back_edges, seen_back_edges);
            graph.follow_edges(next, !extract_left, [&](const handle_t& prev) {
                if (prev != start) {
                    // mark all edges unseen except the one we're traversing
                    frontier_index[next][prev].emplace(entry.first);
                }
            });
            
#ifdef debug_incremental_subgraph
            cerr << "initialized frontier with node " << graph.get_id(next) << " " << graph.get_is_reverse(next) << " at " << &(*entry.first) << ", which has unseen backward edges to:" << endl;
            for (auto h : *unseen_back_edges) {
                cerr << "\t" << graph.get_id(h) << " " << graph.get_is_reverse(h) << endl;
            }
            cerr << "allocated unseen edges at " << unseen_back_edges << endl;
#endif
        });
    }
        
}

IncrementalSubgraph::~IncrementalSubgraph() {
    for (auto& record : frontier) {
        delete get<2>(record);
        delete get<3>(record);
    }
}

bool IncrementalSubgraph::is_extendable() const {
    return !frontier.empty();
}

handle_t IncrementalSubgraph::extend() {
    // get frontier group with fewest uncovered edges (breaking ties arbitrarily)
    auto it = frontier.begin();
    auto nearest = *it;
    
//    cerr << "frontier state:" << endl;
//    for (const auto& val : frontier) {
//        cerr << "\tnode " << graph->get_id(get<1>(val)) << " " << graph->get_is_reverse(get<1>(val)) << ", dist " << get<0>(val) << ", num unseen " << get<2>(val)->size() << ", num seen " << get<3>(val)->size() << endl;
//    }
    
    if ((get_node_count() + 1) % 1000 == 0 ) {
        cerr << "extracting to size " << get_node_count() + 1 << ", frontier size " << frontier.size() << endl;
        
    }
    
#ifdef debug_incremental_subgraph
    cerr << "####" << endl;
    cerr << "extracting a copy of " << graph->get_id(get<1>(nearest)) << " " << graph->get_is_reverse(get<1>(nearest)) << " at " << &(*it) << ", which has remaining unseen backward edges at " << get<2>(nearest) << " to:" << endl;
    for (auto h : *get<2>(nearest)) {
        cerr << "\t" << graph->get_id(h) << " " << graph->get_is_reverse(h) << endl;
    }
    cerr << "and seen backward edges to:" << endl;
    for (auto i : *get<3>(nearest)) {
        auto h = get<0>(extracted[i]);
        cerr << "\t" << i << " (" << graph->get_id(h) << " " << graph->get_is_reverse(h) << ")" << endl;
    }
#endif
    
    // remove the node from the frontier
    frontier.erase(it);
    
    // the index for copies of this node in the frontier
    auto& target_index = frontier_index.at(get<1>(nearest));
    for (auto prev : *get<2>(nearest)) {
        // the index for all unseen predecessors of this frontier node
        auto& source_index = target_index.at(prev);
        // the frontier node we're removing should always be the first one
        source_index.erase(source_index.begin());
    }
    
    // make a node in the extracted subgraph
    extracted.emplace_back();
    auto& extracted_record = extracted.back();
    get<0>(extracted_record) = get<1>(nearest);
    
    // add the edges
    get<1>(extracted_record) = move(*get<3>(nearest));
    for (size_t prev : get<1>(extracted_record)) {
        get<2>(extracted[prev]).push_back(extracted.size() - 1);
    }
    
    // retrieve minimum distance
    get<3>(extracted_record) = get<0>(nearest);
    // figure out the maximum distance
    get<4>(extracted_record) = 0;
    for (size_t prev : get<1>(extracted_record)) {
        auto& prev_record = extracted[prev];
        get<4>(extracted_record) = max<int64_t>(get<4>(extracted_record),
                                                get<4>(prev_record) + graph->get_length(get<0>(prev_record)));
    }
    
#ifdef debug_incremental_subgraph
    cerr << "min distance: " << get<3>(extracted_record) << ", max distance: " << get<4>(extracted_record) << endl;
#endif
    
    // the distance to the end of this node
    int64_t dist_thru = get<0>(nearest) + graph->get_length(get<1>(nearest));
    
    graph->follow_edges(get<1>(nearest), extract_left, [&](const handle_t& next) {
//        cerr << "handling next node " << graph->get_id(next) << " " << graph->get_is_reverse(next) << endl;
        // see if we can mark this edge on one of the copies of the next node
        // that are currently in the frontier
        bool marked_edge = false;
        auto index_iter = frontier_index.find(next);
        if (index_iter != frontier_index.end()) {
//            cerr << "has existed in frontier" << endl;
            // this node exists in the frontier
            auto source_iter = index_iter->second.find(get<1>(nearest));
            if (source_iter != index_iter->second.end() && !source_iter->second.empty()) {
//                cerr << "a frontier copy needs this edge marked" << endl;
                // there are copies of this node in the frontier that haven't had
                // this edge marked yet, get the highest priority copy
                auto frontier_iter = *source_iter->second.begin();
                auto frontier_entry = *frontier_iter;
                
                
//                cerr << "frontier state:" << endl;
//                for (const auto& val : frontier) {
//                    cerr << "node " << graph->get_id(get<1>(val)) << " " << graph->get_is_reverse(get<1>(val)) << ", dist " << get<0>(val) << ", edges " << get<2>(val) << " " << get<3>(val) << ", location " << &val << endl;
//                    for (auto n : *get<2>(val)) {
//                        cerr << "\t" << graph->get_id(n) << " " << graph->get_is_reverse(n) << endl;
//                    }
//                }
//                cerr << "frontier index state:" << endl;
//                for (const auto& outer_rec : frontier_index) {
//                    cerr << "target node " << graph->get_id(outer_rec.first) << " " << graph->get_is_reverse(outer_rec.first) << endl;
//                    for (const auto& inner_rec : outer_rec.second) {
//                        cerr << "\tsource node " << graph->get_id(inner_rec.first) << " " << graph->get_is_reverse(inner_rec.first) << endl;
//                        for (auto fiter : inner_rec.second) {
//                            cerr << "\t\t" << &(*fiter) << endl;
//                        }
//                    }
//                }
                // remove this frontier entry from everywhere that the index is holding it
//                cerr << "removing this frontier entry from the index, prevs at " << get<2>(frontier_entry) << endl;
                for (auto prev : *get<2>(frontier_entry)) {
                    
                    index_iter->second.at(prev).erase(frontier_iter);
                }
                frontier.erase(frontier_iter);
//                cerr << "updating the frontier entry" << endl;
                // mark the unseen edge and move it to the seen edges
                get<2>(frontier_entry)->erase(get<1>(nearest));
                get<3>(frontier_entry)->push_back(extracted.size() - 1);
                // possibly update the minimum distance
                get<0>(frontier_entry) = min(get<0>(frontier_entry), dist_thru);
                // put it back in the frontier
                auto new_entry_iter = frontier.emplace(frontier_entry);
                // reinsert the new iterator everwhere it needs to go in the
                // random access index
                for (auto prev : *get<2>(frontier_entry)) {
                    index_iter->second.at(prev).emplace(new_entry_iter.first);
                }
                
#ifdef debug_incremental_subgraph
                cerr << "updated frontier node " << graph->get_id(next) << " " << graph->get_is_reverse(next) << ", which has unseen backward edges to:" << endl;
                for (auto h : *get<2>(frontier_entry)) {
                    cerr << "\t" << graph->get_id(h) << " " << graph->get_is_reverse(h) << endl;
                }
#endif
                
                marked_edge = true;
            }
            
//
//            for (auto copy_iter = iter->second.begin(), end = iter->second.end(); copy_iter != end; ++copy_iter) {
//                auto frontier_iter = *copy_iter;
//                // TODO: another data structure might make it possible to not iterate
//                // over all copies of this node in the frontier
//                auto incoming = get<2>(*frontier_iter);
//                auto incoming_iter = incoming->find(get<1>(nearest));
//                if (incoming_iter != incoming->end()) {
//                    // this frontier node has an unmarked edge from the node
//                    // we are extracting
//                    // remove from the frontier
//                    auto frontier_entry = *frontier_iter;
//                    frontier.erase(frontier_iter);
//                    // mark the unseen edge and move it to the seen edges
//                    incoming->erase(incoming_iter);
//                    get<3>(frontier_entry)->push_back(extracted.size() - 1);
//                    // possibly update the minimum distance
//                    get<0>(frontier_entry) = min(get<0>(frontier_entry), dist_thru);
//                    // put it back in the frontier
//                    auto new_entry_iter = frontier.emplace(frontier_entry);
//                    // and update our random access index
//                    iter->second.erase(copy_iter);
//                    iter->second.emplace(new_entry_iter.first);
//
//                    marked_edge = true;
//                    break;
//                }
//            }
        }
        
//        cerr << "frontier state:" << endl;
//        for (const auto& val : frontier) {
//            cerr << "node " << graph->get_id(get<1>(val)) << " " << graph->get_is_reverse(get<1>(val)) << ", dist " << get<0>(val) << ", edges " << get<2>(val) << " " << get<3>(val) << ", location " << &val << endl;
//            for (auto n : *get<2>(val)) {
//                cerr << "\t" << graph->get_id(n) << " " << graph->get_is_reverse(n) << endl;
//            }
//        }
//        cerr << "frontier index state:" << endl;
//        for (const auto& outer_rec : frontier_index) {
//            cerr << "target node " << graph->get_id(outer_rec.first) << " " << graph->get_is_reverse(outer_rec.first) << endl;
//            for (const auto& inner_rec : outer_rec.second) {
//                cerr << "\tsource node " << graph->get_id(inner_rec.first) << " " << graph->get_is_reverse(inner_rec.first) << endl;
//                for (auto fiter : inner_rec.second) {
//                    cerr << "\t\t" << &(*fiter) << endl;
//                }
//            }
//        }
        
        if (!marked_edge && dist_thru < max_distance) {
            // we need to add a new copy of this node to the frontier
//            cerr << "adding a new frontier copy" << endl;
            
            auto unseen_back_edges = new unordered_set<handle_t>();
//            cerr << "identifying unseen back edges" << endl;
            graph->follow_edges(next, !extract_left, [&](const handle_t& prev) {
                if (prev != get<1>(nearest)) {
                    // mark all edges unseen except the one we're traversing
                    unseen_back_edges->emplace(prev);
                }
            });
            auto seen_back_edges = new vector<size_t>(1, extracted.size() - 1);
            // add the frontier and its random access index
//            cerr << "enter it into the frontier" << endl;
            auto entry = frontier.emplace(dist_thru, next, unseen_back_edges, seen_back_edges);
//            cerr << "add the frontier node into the frontier index" << endl;
            auto& successor_index = frontier_index[next];
            graph->follow_edges(next, !extract_left, [&](const handle_t& prev) {
                if (prev != get<1>(nearest)) {
                    // mark all edges unseen except the one we're traversing
//                    cerr << "doing it for prev node " << graph->get_id(prev) << " " << graph->get_is_reverse(prev) << endl;
//                    cerr << "exists for next? " << frontier_index.count(next) << endl;
//                    cerr << "exists for prev? " << frontier_index[next].count(prev) << endl;
                    successor_index[prev].emplace(entry.first);
                }
            });
#ifdef debug_incremental_subgraph
            cerr << "added new frontier copy of node " << graph->get_id(next) << " " << graph->get_is_reverse(next) << " at " << &(*entry.first) << ", which has unseen backward edges to:" << endl;
            for (auto h : *unseen_back_edges) {
                cerr << "\t" << graph->get_id(h) << " " << graph->get_is_reverse(h) << endl;
            }
            cerr << "allocated unseen edges at " << unseen_back_edges << endl;
#endif
        }
    });
    
    // clean up the heap objects
//    cerr << "deleting unseen edges at " << get<2>(nearest) << endl;
    delete get<2>(nearest);
    delete get<3>(nearest);
    
//    cerr << "frontier state:" << endl;
//    for (const auto& val : frontier) {
//        cerr << "node " << graph->get_id(get<1>(val)) << " " << graph->get_is_reverse(get<1>(val)) << ", dist " << get<0>(val) << ", edges " << get<2>(val) << " " << get<3>(val) << ", location " << &val << endl;
//        for (auto n : *get<2>(val)) {
//            cerr << "\t" << graph->get_id(n) << " " << graph->get_is_reverse(n) << endl;
//        }
//    }
//    cerr << "frontier index state:" << endl;
//    for (const auto& outer_rec : frontier_index) {
//        cerr << "target node " << graph->get_id(outer_rec.first) << " " << graph->get_is_reverse(outer_rec.first) << endl;
//        for (const auto& inner_rec : outer_rec.second) {
//            cerr << "\tsource node " << graph->get_id(inner_rec.first) << " " << graph->get_is_reverse(inner_rec.first) << endl;
//            for (auto fiter : inner_rec.second) {
//                cerr << "\t\t" << &(*fiter) << endl;
//            }
//        }
//    }
    
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

