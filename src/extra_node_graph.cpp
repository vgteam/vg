#include "extra_node_graph.hpp"
#include "utility.hpp"

#include <handlegraph/util.hpp>

//#define debug

namespace vg {

using namespace std;
using namespace handlegraph;

ExtraNodeGraph::ExtraNodeGraph(
    const HandleGraph* backing, 
    const string& sequence,
    const vector<handle_t>& edges_in,
    const vector<handle_t>& edges_out) :
    backing(backing), 
    in_from(edges_in.begin(), edges_in.end()),
    out_to(edges_out.begin(), edges_out.end()),
    added_id(backing->max_node_id() + 1),
    sequence(sequence) {
   
    // Nothing to do!
}

handle_t ExtraNodeGraph::get_created_handle() const {
    return added_fwd;
}

bool ExtraNodeGraph::has_node(id_t node_id) const {
    return node_id == added_id || backing->has_node(node_id);
}

handle_t ExtraNodeGraph::get_handle(const id_t& node_id, bool is_reverse) const {
    if (node_id == added_id) {
        // They asked for the added node
        return is_reverse ? added_rev : added_fwd;
    } else {
        // Otherwise they asked for something in the backing graph
        handle_t backing_handle = backing->get_handle(node_id, is_reverse);
        
        // Budge up to make room for the added node in each orientation
        return from_backing(backing_handle);
    }
}

id_t ExtraNodeGraph::get_id(const handle_t& handle) const {
    if (handle == added_fwd || handle == added_rev) {
        return added_id;
    } else {
        return backing->get_id(to_backing(handle));
    }
}

bool ExtraNodeGraph::get_is_reverse(const handle_t& handle) const {
    if (handle == added_fwd) {
        return false;
    } else if (handle == added_rev) {
        return true;
    } else {
        return backing->get_is_reverse(to_backing(handle));
    }
}

handle_t ExtraNodeGraph::flip(const handle_t& handle) const {
    if (is_ours(handle)) {
        // In our block of handles, orientation is the low bit
        return as_handle(as_integer(handle) ^ 1);
    } else {
        // Make the backing graph flip it
        return from_backing(backing->flip(to_backing(handle)));
    }
}

size_t ExtraNodeGraph::get_length(const handle_t& handle) const {
    if (is_ours(handle)) {
        // Our node is the same length in any orientation
        return sequence.length();
    } else {
        return backing->get_length(to_backing(handle));
    }
}

string ExtraNodeGraph::get_sequence(const handle_t& handle) const {
    if (handle == added_fwd) {
        return sequence;
    } else if (handle == added_rev) {
        return reverse_complement(sequence);
    } else {
        assert(!is_ours(handle));
        return backing->get_sequence(to_backing(handle));
    }
}

bool ExtraNodeGraph::follow_edges_impl(const handle_t& handle, bool go_left, const function<bool(const handle_t&)>& iteratee) const {
    if (is_ours(handle)) {
        // This is some orientation of our added node
        
        if (handle == added_fwd) {
            // We are in the forward orientation
            if (go_left) {
                // We are going left, so get everywhere we come from
                for (auto& prev : in_from) {
                    if (!iteratee(from_backing(prev))) {
                        return false;
                    }
                }
            } else {
                // We are going right, so get everything we go to
                for (auto& next : out_to) {
                    if (!iteratee(from_backing(next))) {
                        return false;
                    }
                }
            }
        } else if (handle == added_rev) {
            // We are going in the reverse orientation
            if (go_left) {
                // We actually want where we go to, backward
                for (auto& next : out_to) {
                    if (!iteratee(flip(from_backing(next)))) {
                        return false;
                    }
                }
            } else {
                // We actually want where we come from, backward
                for (auto& prev : in_from) {
                    if (!iteratee(flip(from_backing(prev)))) {
                        return false;
                    }
                }
            }
        }
        
        return true;
    } else {
        // The handle refers to a node in the backing graph
        auto backing_handle = to_backing(handle);
        
        // If we get through those, do the actual edges in the backing graph
        bool keep_going = backing->follow_edges(backing_handle, go_left, [&](const handle_t& found) -> bool {
            return iteratee(from_backing(found));
        });
        
        // Also handle edges to/from our added node.
        
        if (keep_going &&
            ((go_left && out_to.count(backing_handle)) ||
            (!go_left && in_from.count(backing_handle)))) {
            // Visit it forward
            keep_going &= iteratee(added_fwd);
        }
        
        if (keep_going &&
            ((go_left && out_to.count(backing->flip(backing_handle))) ||
            (!go_left && in_from.count(backing->flip(backing_handle))))) {
            // Visit it reverse
            keep_going &= iteratee(added_rev);
        }
        
        return keep_going;
        
    }
}

bool ExtraNodeGraph::for_each_handle_impl(const function<bool(const handle_t&)>& iteratee, bool parallel) const {
    
    // First do the node we added
    if (!iteratee(added_fwd)) {
        return false;
    }
    
#ifdef debug
    cerr << "Try backing graph " << (parallel ? "in parallel" : "") << endl;
#endif
    return backing->for_each_handle([&](const handle_t& backing_handle) -> bool {
        // Now do each backing node, possibly in parallel.
#ifdef debug
        cerr << "Invoke iteratee on " << backing->get_id(backing_handle) << endl;
#endif
        return iteratee(from_backing(backing_handle));
    }, parallel);
}

size_t ExtraNodeGraph::get_node_count() const {
    return backing->get_node_count() + 1;
}

id_t ExtraNodeGraph::min_node_id() const {
    return min(backing->min_node_id(), added_id);
}
    
id_t ExtraNodeGraph::max_node_id() const {
    return max(backing->max_node_id(), added_id);
}

size_t ExtraNodeGraph::get_degree(const handle_t& handle, bool go_left) const {
    if (is_ours(handle)) {
        if ((handle == added_fwd && !go_left) || (handle == added_rev && go_left)) {
            // Edges out of the added node
            return out_to.size();
        } else if ((handle == added_fwd && go_left) || (handle == added_rev && !go_left)) {
            // Edges into the added node
            return in_from.size();
        }
    } else {
        // We need to find the backing graph degree and possibly adjust it if this is a head or tail
        handle_t backing_handle = to_backing(handle);
        
        size_t degree = backing->get_degree(backing_handle, go_left);
        
        if ((go_left && out_to.count(backing_handle)) || (!go_left && in_from.count(backing_handle))) {
            // Forward version of this handle connects to the added node on this side.
            degree++;
        }
        
        if ((go_left && in_from.count(backing->flip(backing_handle))) || (!go_left && out_to.count(backing->flip(backing_handle)))) {
            // Reverse version of this handle connects to added node on this side.
            degree++;
        }
        
        return degree;
    }
    
    // We must return from one of the other branches
    throw runtime_error("Did not hit a return statement that should have been hit");
}
    

}
