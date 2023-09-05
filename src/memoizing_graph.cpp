/**
 * \file memoizing_graph.cpp: contains the implementation of MemoizingGraph
 */


#include "memoizing_graph.hpp"


namespace vg {

using namespace std;

    MemoizingGraph::MemoizingGraph(const PathPositionHandleGraph* graph) : graph(graph) {
        
    }
    
    bool MemoizingGraph::has_node(id_t node_id) const {
        bool found_node = false;
        if (graph->has_node(node_id)) {
            if (graph->get_length(graph->get_handle(node_id)) > 0) {
                found_node = true;
            }
        }
        return found_node;
    }
    
    handle_t MemoizingGraph::get_handle(const id_t& node_id, bool is_reverse) const {
        // we have to do some ugly stuff to keep libhandlegraph's const requirements while still
        // updating memos
        auto& memo = const_cast<MemoizingGraph*>(this)->get_handle_memo;
        
        handle_t to_return;
        auto it = memo.find(node_id);
        if (it != memo.end()) {
            to_return = is_reverse ? graph->flip(it->second) : it->second;
        }
        else if (memo.size() < max_handle_memo_size) {
            handle_t handle = graph->get_handle(node_id);
            memo[node_id] = handle;
            to_return = is_reverse ? graph->flip(handle) : handle;
        }
        else {
            to_return = graph->get_handle(node_id, is_reverse);
        }
        return to_return;
    }
    
    id_t MemoizingGraph::get_id(const handle_t& handle) const {
        return graph->get_id(handle);
    }
    
    bool MemoizingGraph::get_is_reverse(const handle_t& handle) const {
        return graph->get_is_reverse(handle);
    }
    
    handle_t MemoizingGraph::flip(const handle_t& handle) const {
        return graph->flip(handle);
    }
    
    size_t MemoizingGraph::get_length(const handle_t& handle) const {
        return graph->get_length(handle);
    }
    
    string MemoizingGraph::get_sequence(const handle_t& handle) const {
        return graph->get_sequence(handle);
    }
    
    bool MemoizingGraph::follow_edges_impl(const handle_t& handle, bool go_left,
                                           const function<bool(const handle_t&)>& iteratee) const {
        
        return graph->follow_edges(handle, go_left, [&](const handle_t& next) {
            return iteratee(next);
        });
    }
    
    bool MemoizingGraph::for_each_handle_impl(const function<bool(const handle_t&)>& iteratee, bool parallel) const {
        return graph->for_each_handle([&](const handle_t& handle) {
            return iteratee(handle);
        }, parallel);
    }
    
    size_t MemoizingGraph::get_node_count() const {
        return graph->get_node_count();
    }
    
    id_t MemoizingGraph::min_node_id() const {
        return graph->min_node_id();
    }
    
    id_t MemoizingGraph::max_node_id() const {
        return graph->max_node_id();
    }
    
    size_t MemoizingGraph::get_path_count() const {
        return graph->get_path_count();
    }
    
    bool MemoizingGraph::has_path(const std::string& path_name) const {
        return graph->has_path(path_name);
    }
    
    path_handle_t MemoizingGraph::get_path_handle(const std::string& path_name) const {
        return graph->get_path_handle(path_name);
    }
    
    std::string MemoizingGraph::get_path_name(const path_handle_t& path_handle) const {
        return graph->get_path_name(path_handle);
    }
    
    bool MemoizingGraph::get_is_circular(const path_handle_t& path_handle) const {
        return graph->get_is_circular(path_handle);
    }
    
    size_t MemoizingGraph::get_step_count(const path_handle_t& path_handle) const {
        return graph->get_step_count(path_handle);
    }
    
    handle_t MemoizingGraph::get_handle_of_step(const step_handle_t& step_handle) const {
        return graph->get_handle_of_step(step_handle);
    }
    
    path_handle_t MemoizingGraph::get_path_handle_of_step(const step_handle_t& step_handle) const {
        return graph->get_path_handle_of_step(step_handle);
    }
    
    step_handle_t MemoizingGraph::path_begin(const path_handle_t& path_handle) const {
        return graph->path_begin(path_handle);
    }
    
    step_handle_t MemoizingGraph::path_end(const path_handle_t& path_handle) const {
        return graph->path_end(path_handle);
    }
    
    step_handle_t MemoizingGraph::path_back(const path_handle_t& path_handle) const {
        return graph->path_back(path_handle);
    }
    
    step_handle_t MemoizingGraph::path_front_end(const path_handle_t& path_handle) const {
        return graph->path_front_end(path_handle);
    }
    
    bool MemoizingGraph::has_next_step(const step_handle_t& step_handle) const {
        return graph->has_next_step(step_handle);
    }
    
    bool MemoizingGraph::has_previous_step(const step_handle_t& step_handle) const {
        return graph->has_previous_step(step_handle);
    }
    
    step_handle_t MemoizingGraph::get_next_step(const step_handle_t& step_handle) const {
        return graph->get_next_step(step_handle);
    }
    
    step_handle_t MemoizingGraph::get_previous_step(const step_handle_t& step_handle) const {
        return graph->get_previous_step(step_handle);
    }

    bool MemoizingGraph::for_each_path_handle_impl(const std::function<bool(const path_handle_t&)>& iteratee) const {
        return graph->for_each_path_handle(iteratee);
    }
    
    bool MemoizingGraph::for_each_step_on_handle_impl(const handle_t& handle,
                                                      const std::function<bool(const step_handle_t&)>& iteratee) const {
        return graph->for_each_step_on_handle(handle, iteratee);
    }
    
    std::vector<step_handle_t> MemoizingGraph::steps_of_handle(const handle_t& handle,
                                                               bool match_orientation) const {
        
        // we have to do some ugly stuff to keep libhandlegraph's const requirements while still
        // updating memos
        auto& memo = const_cast<MemoizingGraph*>(this)->steps_of_handle_memo;

        vector<step_handle_t> to_return;
        auto it = memo.find(forward(handle));
        if (it != memo.end()) {
            if (match_orientation) {
                for (const step_handle_t& step : it->second) {
                    if (graph->get_is_reverse(graph->get_handle_of_step(step)) == graph->get_is_reverse(handle)) {
                        to_return.push_back(step);
                    }
                }
            }
            else {
                to_return = it->second;
            }
        }
        else if (memo.size() < max_steps_of_handle_memo_size) {
            memo[forward(handle)] = graph->steps_of_handle(handle);
            if (match_orientation) {
                for (const step_handle_t& step : memo[forward(handle)]) {
                    if (graph->get_is_reverse(graph->get_handle_of_step(step)) == graph->get_is_reverse(handle)) {
                        to_return.push_back(step);
                    }
                }
            }
            else {
                to_return = memo[forward(handle)];
            }
        }
        else {
            to_return = graph->steps_of_handle(handle, match_orientation);
        }
        return to_return;
    }
    
    bool MemoizingGraph::is_empty(const path_handle_t& path_handle) const {
        return graph->is_empty(path_handle);
    }
    
    size_t MemoizingGraph::get_path_length(const path_handle_t& path_handle) const {
        return graph->get_path_length(path_handle);
    }
    
    size_t MemoizingGraph::get_position_of_step(const step_handle_t& step) const {
        return graph->get_position_of_step(step);
    }
    
    step_handle_t MemoizingGraph::get_step_at_position(const path_handle_t& path,
                                                       const size_t& position) const {
        return graph->get_step_at_position(path, position);
    }
    
    bool MemoizingGraph::for_each_step_position_on_handle(const handle_t& handle,
                                                          const std::function<bool(const step_handle_t&, const bool&, const size_t&)>& iteratee) const {
        return graph->for_each_step_position_on_handle(handle, iteratee);
    }
}

