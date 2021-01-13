#include "source_sink_overlay.hpp"

#include <handlegraph/util.hpp>

//#define debug

namespace vg {

using namespace std;
using namespace handlegraph;

SourceSinkOverlay::SourceSinkOverlay(const HandleGraph* backing, size_t length, id_t source_id, id_t sink_id,
    bool break_disconnected) : node_length(length), backing(backing), source_id(source_id), sink_id(sink_id) {
   
    // Both IDs or neither must be specified.
    assert((this->source_id == 0) == (this->sink_id == 0));
   
    if (this->source_id == 0 || this->sink_id == 0) {
        // We need to autodetect our source and sink IDs
        id_t backing_max_id = backing->get_node_count() > 0 ? backing->max_node_id() : 0;
        
        this->source_id = backing_max_id + 1;
        this->sink_id = this->source_id + 1;
    }
    
#ifdef debug
    cerr << "Make overlay for kmer size " << length << " with source " << this->source_id << " and sink " << this->sink_id << endl;
#endif
    
    // We have to divide the graph into connected components and get ahold of the tips.
    vector<pair<unordered_set<id_t>, vector<handle_t>>> components = handlealgs::weakly_connected_components_with_tips(backing);
    
    for (auto& component : components) {
        // Unpack each component
        auto& component_ids = component.first;
        auto& component_tips = component.second;
        
#ifdef debug
        cerr << "Weakly connected component of " << component_ids.size() << " has " << component_tips.size() << " tips:" << endl;
        for (auto& tip : component_tips) {
            cerr << "\t" << backing->get_id(tip) << " orientation " << backing->get_is_reverse(tip) << endl;
        }
#endif
        
        // All the components need to be nonempty
        assert(!component_ids.empty());
        
        for (auto& handle : component_tips) {
            // We need to cache the heads and tails as sets of handles, so we know to
            // make edges to all of them when reading out of our synthetic source and
            // sink nodes.
            
            if (backing->get_is_reverse(handle)) {
                // It's a tail. Insert it forward as a tail.
                backing_tails.insert(backing->flip(handle));
            } else {
                // It's a head
                backing_heads.insert(handle);
            }
            
        }
        
        if (component_tips.empty() && break_disconnected) {
            // If we're supposed to break open cycles, we also mix in an arbitrary node
            // from each tipless component as a head, and each handle that reads into
            // it as a tail.
            
            // Choose a fake head arbitrarily
            handle_t fake_head = backing->get_handle(*component_ids.begin(), false);
            backing_heads.insert(fake_head);
            
            // Find the fake tails that are to the left of it
            backing->follow_edges(fake_head, true, [&](const handle_t& fake_tail) {
                backing_tails.insert(fake_tail);
            });
        }
    }
    
    
}

handle_t SourceSinkOverlay::get_source_handle() const {
    return source_fwd;
}

handle_t SourceSinkOverlay::get_sink_handle() const {
    return sink_fwd;
}

bool SourceSinkOverlay::has_node(id_t node_id) const {
    return backing->has_node(node_id);
}

handle_t SourceSinkOverlay::get_handle(const id_t& node_id, bool is_reverse) const {
    if (node_id == source_id) {
        // They asked for the source node
        return is_reverse ? source_rev : source_fwd;
    } else if (node_id == sink_id) {
        // They asked for the sink node
        return is_reverse ? sink_rev : sink_fwd;
    } else {
        // Otherwise they asked for something in the backing graph
        handle_t backing_handle = backing->get_handle(node_id, is_reverse);
        
        // Budge up to make room for the source and sink in each orientation
        return as_handle(as_integer(backing_handle) + 4);
    }
}

id_t SourceSinkOverlay::get_id(const handle_t& handle) const {
    if (handle == source_fwd || handle == source_rev) {
        return source_id;
    } else if (handle == sink_fwd || handle == sink_rev) {
        return sink_id;
    } else {
        return backing->get_id(to_backing(handle));
    }
}

bool SourceSinkOverlay::get_is_reverse(const handle_t& handle) const {
    if (handle == source_fwd || handle == sink_fwd) {
        return false;
    } else if (handle == source_rev || handle == sink_rev) {
        return true;
    } else {
        return backing->get_is_reverse(to_backing(handle));
    }
}

handle_t SourceSinkOverlay::flip(const handle_t& handle) const {
    if (is_ours(handle)) {
        // In our block of two handles, orientation is the low bit
        return as_handle(as_integer(handle) ^ 1);
    } else {
        // Make the backing graph flip it
        return from_backing(backing->flip(to_backing(handle)));
    }
}

size_t SourceSinkOverlay::get_length(const handle_t& handle) const {
    if (is_ours(handle)) {
        // Both our fake nodes are the same length
        return node_length;
    } else {
        return backing->get_length(to_backing(handle));
    }
}

string SourceSinkOverlay::get_sequence(const handle_t& handle) const {
    if (handle == source_fwd || handle == sink_rev) {
        // Reading into the graph is all '#'
        return string(node_length, '#');
    } else if (handle == source_rev || handle == sink_fwd) {
        // Reading out of the graph is all '$'
        return string(node_length, '$');
    } else {
        assert(!is_ours(handle));
        return backing->get_sequence(to_backing(handle));
    }
}

bool SourceSinkOverlay::follow_edges_impl(const handle_t& handle, bool go_left, const function<bool(const handle_t&)>& iteratee) const {
    if (is_ours(handle)) {
        // We only care about the right of the source and the left of the sink
        if ((handle == source_fwd && !go_left) || (handle == source_rev && go_left)) {
            // We want the right of the source (every head node in the backing graph)
            // Make sure to put it in the appropriate orientation.
            
            for (const handle_t& backing_head : backing_heads) {
                // Feed each backing graph head to the iteratee, in the
                // appropriate orientation depending on which way we want to
                // go.
                if (!iteratee(from_backing(go_left ? backing->flip(backing_head) : backing_head))) {
                    // If they say to stop, stop
                    return false;
                }
            }
            
        } else if ((handle == sink_fwd && go_left) || (handle == sink_rev && !go_left)) {
            // We want the left of the sink (every tail node in the backing graph)
            // Make sure to put it in the appropriate orientation.
            
            for (const handle_t& backing_tail : backing_tails) {
                // Feed each backing graph tail to the iteratee, in the
                // appropriate orientation depending on which way we want to
                // go.
                if (!iteratee(from_backing(go_left ? backing_tail : backing->flip(backing_tail)))) {
                    // If they say to stop, stop
                    return false;
                }
            }
        }
        return true;
    } else {
        // The handle refers to a node in the backing graph
        auto backing_handle = to_backing(handle);
        
        if ((backing_heads.count(backing_handle) && go_left) || (backing_heads.count(backing->flip(backing_handle)) && !go_left)) {
            // We want to read left off a head (possibly in reverse) into the synthetic source
            if (!iteratee(go_left ? source_fwd : source_rev)) {
                // If they say stop, stop
                return false;
            }
        }
        
        if ((backing_tails.count(backing_handle) && !go_left) || (backing_tails.count(backing->flip(backing_handle)) && go_left)) {
            // We want to read right off a tail (possibly in reverse) into the synthetic sink
            if (!iteratee(go_left ? sink_rev : sink_fwd)) {
                // If they say stop, stop
                return false;
            }
        }
    
        // If we get through those, do the actual edges in the backing graph
        return backing->follow_edges(backing_handle, go_left, [&](const handle_t& found) -> bool {
            return iteratee(from_backing(found));
        });
    }
}

bool SourceSinkOverlay::for_each_handle_impl(const function<bool(const handle_t&)>& iteratee, bool parallel) const {
    
    // First do the sourece and sink we added
    if (!iteratee(source_fwd)) {
        return false;
    }
    if (!iteratee(sink_fwd)) {
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

size_t SourceSinkOverlay::get_node_count() const {
    return backing->get_node_count() + 2;
}

id_t SourceSinkOverlay::min_node_id() const {
    return min(backing->min_node_id(), min(source_id, sink_id));
}
    
id_t SourceSinkOverlay::max_node_id() const {
    return max(backing->max_node_id(), max(source_id, sink_id));
}

size_t SourceSinkOverlay::get_degree(const handle_t& handle, bool go_left) const {
    if (is_ours(handle)) {
        if ((handle == source_fwd && !go_left) || (handle == source_rev && go_left)) {
            // We are reading into every graph head
            return backing_heads.size();
        } else if ((handle == sink_fwd && go_left) || (handle == sink_rev && !go_left)) {
            // We are reading into every graph tail
            return backing_tails.size();
        }
        // Otherwise we're reading off the outside ends of the source/sink nodes
        return 0;
    } else {
        // We need to find the backing graph degree and possibly adjust it if this is a head or tail
        handle_t backing_handle = to_backing(handle);
        
        size_t degree = backing->get_degree(backing_handle, go_left);
        
        if (backing_heads.count(backing->forward(backing_handle))) {
            // We are a head. Are we going off the left end when forward, or the right end when reverse?
            if (go_left != backing->get_is_reverse(backing_handle)) {
                // If so we count the synthetic edge.
                degree++;
            }
        }
        if (backing_tails.count(backing->forward(backing_handle))) {
            // We are a tial. Are we going off the left end when reverse, or the right end when forward?
            if (go_left != !backing->get_is_reverse(backing_handle)) {
                // If so we count the synthetic edge.
                degree++;
            }
        }
        
        return degree;
    }
}

handle_t SourceSinkOverlay::get_underlying_handle(const handle_t& handle) const {
    if (is_ours(handle)) {
        throw std::runtime_error("error:[SourceSinkOverlay] cannot request underlying handle of source or sink node");
    }
    return to_backing(handle);
}

}
