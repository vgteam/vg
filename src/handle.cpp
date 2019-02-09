#include "handle.hpp"
#include "snarls.hpp"

/** \file handle.cpp
 * Implement handle graph utility methods.
 */

namespace vg {

using namespace std;

handle_t HandleGraph::get_handle(const Visit& visit) const {
    return get_handle(visit.node_id(), visit.backward());
}

size_t HandleGraph::get_degree(const handle_t& handle, bool go_left) const {
    size_t count = 0;
    follow_edges(handle, go_left, [&](const handle_t& ignored) {
        // Just manually count every edge we get by looking at the handle in that orientation
        count++;
    });
    return count;
}

bool HandleGraph::has_edge(const handle_t& left, const handle_t& right) const {
    bool not_seen = true;
    follow_edges(left, false, [&](const handle_t& next) {
        not_seen = (next != right);
        return not_seen;
    });
    return !not_seen;
}

Visit HandleGraph::to_visit(const handle_t& handle) const {
    return vg::to_visit(this->get_id(handle), this->get_is_reverse(handle));
}

handle_t HandleGraph::forward(const handle_t& handle) const {
    return this->get_is_reverse(handle) ? this->flip(handle) : handle;
}

pair<handle_t, handle_t> HandleGraph::edge_handle(const handle_t& left, const handle_t& right) const {
    // The degeneracy is between any pair and a pair of the same nodes but reversed in order and orientation.
    // We compare those two pairs and construct the smaller one.
    auto flipped_right = this->flip(right);
    
    if (as_integer(left) > as_integer(flipped_right)) {
        // The other orientation would be smaller.
        return make_pair(flipped_right, this->flip(left));
    } else if(as_integer(left) == as_integer(flipped_right)) {
        // Our left and the flipped pair's left would be equal.
        auto flipped_left = this->flip(left);
        if (as_integer(right) > as_integer(flipped_left)) {
            // And our right is too big, so flip.
            return make_pair(flipped_right, flipped_left);
        } else {
            // No difference or we're smaller.
            return make_pair(left, right);
        }
    } else {
        // We're smaller
        return make_pair(left, right);
    }
}

handle_t HandleGraph::traverse_edge_handle(const edge_t& edge, const handle_t& left) const {
    if (left == edge.first) {
        // The cannonical orientation is the one we want
        return edge.second;
    } else if (left == this->flip(edge.second)) {
        // We really want the other orientation
        return this->flip(edge.first);
    } else {
        // This isn't either handle that the edge actually connects. Something has gone wrong.
        throw runtime_error("Cannot view edge " +
            to_string(this->get_id(edge.first)) + " " + to_string(this->get_is_reverse(edge.first)) + " -> " +
            to_string(this->get_id(edge.second)) + " " + to_string(this->get_is_reverse(edge.second)) +
            " from non-participant " + to_string(this->get_id(left)) + " " + to_string(this->get_is_reverse(left)));
    }
}
    
void HandleGraph::for_each_edge(const function<bool(const edge_t&)>& iteratee, bool parallel) const {
    for_each_handle([&](const handle_t& handle){
        bool keep_going = true;
        // filter to edges where this node is lower ID or any rightward self-loops
        follow_edges(handle, false, [&](const handle_t& next) {
            if (get_id(handle) <= get_id(next)) {
                keep_going = iteratee(edge_handle(handle, next));
            }
            return keep_going;
        });
        if (keep_going) {
            // filter to edges where this node is lower ID or leftward reversing
            // self-loop
            follow_edges(handle, true, [&](const handle_t& prev) {
                if (get_id(handle) < get_id(prev) ||
                    (get_id(handle) == get_id(prev) && !get_is_reverse(prev))) {
                    keep_going = iteratee(edge_handle(prev, handle));
                }
                return keep_going;
            });
        }
    }, parallel);
}
    
bool PathHandleGraph::is_empty(const path_handle_t& path_handle) const {
    // By default, we can answer emptiness queries with the length query.
    // But some implementations may have an expensive length query and a cheaper emptiness one
    return get_occurrence_count(path_handle) == 0;
}

void PathHandleGraph::for_each_occurrence_in_path(const path_handle_t& path, const function<void(const occurrence_handle_t&)>& iteratee) const {
    
#ifdef debug
    cerr << "Go over occurrences in path " << get_path_name(path) << " with " <<  get_occurrence_count(path) << " occurrences in it" << endl;
#endif

    if (is_empty(path)) {
        // Nothing to do!
#ifdef debug
        cerr << "Path is empty!" << endl;
#endif
        return;
    }
    
    // Otherwise the path is nonempty so it is safe to try and grab a first occurrence
    auto here = get_first_occurrence(path);
    // Run for the first occurrence
    iteratee(here);
    while (has_next_occurrence(here)) {
        // Run for all subsequent occurrences on the path
        here = get_next_occurrence(here);
        iteratee(here);
    }
}

}


