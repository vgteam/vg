#include <atomic>
#include "subgraph_overlay.hpp"

#include <handlegraph/util.hpp>

namespace vg {

using namespace std;
using namespace handlegraph;

SubgraphOverlay::SubgraphOverlay(const HandleGraph* backing, const unordered_set<nid_t>* node_subset) :
    backing_graph(backing),
    node_subset(node_subset) {
    if (!node_subset->empty()) {
        auto minmax_nodes = std::minmax_element(node_subset->begin(), node_subset->end());
        min_node = *minmax_nodes.first;
        max_node = *minmax_nodes.second;
    } 
}

SubgraphOverlay::~SubgraphOverlay() {
    
}

bool SubgraphOverlay::has_node(nid_t node_id) const {
    return node_subset->count(node_id);
}
   
handle_t SubgraphOverlay::get_handle(const nid_t& node_id, bool is_reverse) const {
    if (has_node(node_id)) {
        return backing_graph->get_handle(node_id, is_reverse);
    } else {
        throw runtime_error("Node " + std::to_string(node_id) + " not in subgraph overlay");
    }
}
    
nid_t SubgraphOverlay::get_id(const handle_t& handle) const {
    return backing_graph->get_id(handle);
}
    
bool SubgraphOverlay::get_is_reverse(const handle_t& handle) const {
    return backing_graph->get_is_reverse(handle);
}

handle_t SubgraphOverlay::flip(const handle_t& handle) const {
    return backing_graph->flip(handle);
}
    
size_t SubgraphOverlay::get_length(const handle_t& handle) const {
    return backing_graph->get_length(handle);
}

std::string SubgraphOverlay::get_sequence(const handle_t& handle) const {
    return backing_graph->get_sequence(handle);
}
    
size_t SubgraphOverlay::get_node_count() const {
    return node_subset->size();
}

nid_t SubgraphOverlay::min_node_id() const {
    return min_node;
}
    
nid_t SubgraphOverlay::max_node_id() const {
    return max_node;
}

bool SubgraphOverlay::follow_edges_impl(const handle_t& handle, bool go_left, const std::function<bool(const handle_t&)>& iteratee) const {
    std::function<bool(const handle_t&)> subgraph_iteratee = [&](const handle_t& handle) {
        if (has_node(backing_graph->get_id(handle))) {
            if (iteratee(handle) == false) {
                return false;
            }
        }
        return true;
    };
    if (has_node(backing_graph->get_id(handle))) {
        return backing_graph->follow_edges(handle, go_left, subgraph_iteratee);
    }
    return true;
}
    
bool SubgraphOverlay::for_each_handle_impl(const std::function<bool(const handle_t&)>& iteratee, bool parallel) const {

    if (!parallel) {
        bool keep_going = true;
        for (auto node_it = node_subset->begin(); keep_going && node_it != node_subset->end(); ++node_it) {
            keep_going = iteratee(get_handle(*node_it, false));
        }
        return keep_going;
    } else {
        // copy them into something easy to iterate with omp
        vector<nid_t> node_vec(node_subset->begin(), node_subset->end());
        std::atomic<bool> keep_going(true);
#pragma omp parallel for
        for (size_t i = 0; i < node_vec.size(); ++i) {
            keep_going = keep_going && iteratee(backing_graph->get_handle(node_vec[i]));
        }
        return keep_going;
    }
}

PathSubgraphOverlay::PathSubgraphOverlay(const PathHandleGraph* backing, const unordered_set<nid_t>* node_subset) :
    SubgraphOverlay(backing, node_subset),
    backing_path_graph(backing) {

    backing->for_each_path_handle([&](const path_handle_t& path_handle) {
            bool fully_contained = true;
            backing->for_each_step_in_path(path_handle, [&](const step_handle_t& step_handle) -> bool {
                    if (!has_node(backing->get_id(backing->get_handle_of_step(step_handle)))) {
                        fully_contained = false;
                        return false;
                    }
                    return true;
                });
            if (fully_contained) {
                path_subset.insert(path_handle);
            }
        });
}

PathSubgraphOverlay::~PathSubgraphOverlay() {
}

size_t PathSubgraphOverlay::get_path_count() const {
    return path_subset.size();
}
    
bool PathSubgraphOverlay::has_path(const std::string& path_name) const {
    return backing_path_graph->has_path(path_name) &&
        path_subset.count(backing_path_graph->get_path_handle(path_name));
}
    
path_handle_t PathSubgraphOverlay::get_path_handle(const std::string& path_name) const {
    if (!has_path(path_name)) {
        throw runtime_error("Path " + path_name + " not in subgraph overlay");
    } else {
        return backing_path_graph->get_path_handle(path_name);
    }
}

std::string PathSubgraphOverlay::get_path_name(const path_handle_t& path_handle) const {
    return backing_path_graph->get_path_name(path_handle);
}
    
bool PathSubgraphOverlay::get_is_circular(const path_handle_t& path_handle) const {
    return backing_path_graph->get_is_circular(path_handle);
}
    
size_t PathSubgraphOverlay::get_step_count(const path_handle_t& path_handle) const {
    return backing_path_graph->get_step_count(path_handle);
}
    
handle_t PathSubgraphOverlay::get_handle_of_step(const step_handle_t& step_handle) const {
    return backing_path_graph->get_handle_of_step(step_handle);
}
    
path_handle_t PathSubgraphOverlay::get_path_handle_of_step(const step_handle_t& step_handle) const {
    return backing_path_graph->get_path_handle_of_step(step_handle);
}

step_handle_t PathSubgraphOverlay::path_begin(const path_handle_t& path_handle) const {
    return backing_path_graph->path_begin(path_handle);
}
    
step_handle_t PathSubgraphOverlay::path_end(const path_handle_t& path_handle) const {
    return backing_path_graph->path_end(path_handle);
}
    
step_handle_t PathSubgraphOverlay::path_back(const path_handle_t& path_handle) const {
    return backing_path_graph->path_back(path_handle);
}
    
step_handle_t PathSubgraphOverlay::path_front_end(const path_handle_t& path_handle) const {
    return backing_path_graph->path_front_end(path_handle);
}

bool PathSubgraphOverlay::has_next_step(const step_handle_t& step_handle) const {
    return backing_path_graph->has_next_step(step_handle);
}

bool PathSubgraphOverlay::has_previous_step(const step_handle_t& step_handle) const {
    return backing_path_graph->has_previous_step(step_handle);
}
    
step_handle_t PathSubgraphOverlay::get_next_step(const step_handle_t& step_handle) const {
    return backing_path_graph->get_next_step(step_handle);
}
    
step_handle_t PathSubgraphOverlay::get_previous_step(const step_handle_t& step_handle) const {
    return backing_path_graph->get_previous_step(step_handle);
}

bool PathSubgraphOverlay::for_each_path_handle_impl(const std::function<bool(const path_handle_t&)>& iteratee) const {
    bool keep_going = true;
    for (auto path_it = path_subset.begin(); keep_going && path_it != path_subset.end(); ++path_it) {
        keep_going = iteratee(*path_it);
    }

    return keep_going;
}

bool PathSubgraphOverlay::for_each_step_on_handle_impl(const handle_t& handle,
                                                       const std::function<bool(const step_handle_t&)>& iteratee) const {
    return backing_path_graph->for_each_step_on_handle(handle, iteratee);
}

}
