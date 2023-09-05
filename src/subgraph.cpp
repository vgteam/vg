/**
 * \file subgraph.cpp: contains the implementation of SubHandleGraph
 */


#include "subgraph.hpp"

#include <atomic>

namespace vg {

using namespace std;

    SubHandleGraph::SubHandleGraph(const HandleGraph* super) : super(super) {
        // nothing to do
    }
    
    void SubHandleGraph::add_handle(const handle_t& handle) {
        
        id_t node_id = super->get_id(handle);
        
        min_id = min(node_id, min_id);
        max_id = max(node_id, max_id);
        
        contents.insert(node_id);
    }
    
    bool SubHandleGraph::has_node(id_t node_id) const {
        return contents.count(node_id);
    }
    
    handle_t SubHandleGraph::get_handle(const id_t& node_id, bool is_reverse) const {
        if (!contents.count(node_id)) {
            cerr << "error:[SubHandleGraph] subgraph does not contain node with ID " << node_id << endl;
            exit(1);
        }
        return super->get_handle(node_id, is_reverse);
    }
    
    id_t SubHandleGraph::get_id(const handle_t& handle) const {
        return super->get_id(handle);
    }
    
    bool SubHandleGraph::get_is_reverse(const handle_t& handle) const {
        return super->get_is_reverse(handle);
    }
    
    handle_t SubHandleGraph::flip(const handle_t& handle) const {
        return super->flip(handle);
    }
    
    size_t SubHandleGraph::get_length(const handle_t& handle) const {
        return super->get_length(handle);
    }
    
    string SubHandleGraph::get_sequence(const handle_t& handle) const {
        return super->get_sequence(handle);
    }
    
    bool SubHandleGraph::follow_edges_impl(const handle_t& handle, bool go_left, const function<bool(const handle_t&)>& iteratee) const {
        // only let it travel along edges whose endpoints are in the subgraph
        bool keep_going = true;
        super->follow_edges(handle, go_left, [&](const handle_t& handle) {
            if (contents.count(super->get_id(handle))) {
                keep_going = iteratee(handle);
            }
            return keep_going;
        });
        return keep_going;
    }
    
    bool SubHandleGraph::for_each_handle_impl(const function<bool(const handle_t&)>& iteratee, bool parallel) const {
        if (parallel) {
            atomic<bool> keep_going(true);
            // do parallelism taskwise inside the iteration
#pragma omp parallel
            {
#pragma omp single
                {
                    for(auto iter = contents.begin(); keep_going && iter != contents.end(); iter++) {
#pragma omp task
                        {
                            if (!iteratee(super->get_handle(*iter))) {
                                keep_going = false;
                            }
                        }
                    }
                }
            }
            return keep_going;
        }
        else {
            // non-parallel
            for (id_t node_id : contents) {
                if (!iteratee(super->get_handle(node_id))) {
                    return false;
                }
            }
            return true;
        }
    }
    
    size_t SubHandleGraph::get_node_count() const {
        return contents.size();
    }
    
    id_t SubHandleGraph::min_node_id() const {
        return min_id;
    }
    
    id_t SubHandleGraph::max_node_id() const {
        return max_id;
    }

    handle_t SubHandleGraph::get_underlying_handle(const handle_t& handle) const {
        return handle;
    }
}

