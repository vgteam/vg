/**
 * \file path_subgraph.cpp: contains the implementation of PathSubgraph
 */


#include "path_subgraph.hpp"
#include "path.hpp"
#include <handlegraph/util.hpp>

namespace vg {

using namespace std;

    PathSubgraph::PathSubgraph(const HandleGraph* base, const Path& path) : super(base), defining_path(path) {
        // Check our input
        assert(defining_path.mapping_size() > 0);
    }
    
    bool PathSubgraph::has_node(id_t node_id) const {
        return (node_id > 0 && node_id <= defining_path.mapping_size());
    }
    
    handle_t PathSubgraph::get_handle(const id_t& node_id, bool is_reverse) const {
        return handlegraph::number_bool_packing::pack(node_id, is_reverse);
    }
    
    id_t PathSubgraph::get_id(const handle_t& handle) const {
        return handlegraph::number_bool_packing::unpack_number(handle);
    }
    
    bool PathSubgraph::get_is_reverse(const handle_t& handle) const {
        return handlegraph::number_bool_packing::unpack_bit(handle);
    }
    
    handle_t PathSubgraph::flip(const handle_t& handle) const {
        return handlegraph::number_bool_packing::toggle_bit(handle);
    }
    
    size_t PathSubgraph::get_length(const handle_t& handle) const {
        size_t index = (size_t)get_id(handle) - 1;
        // No need to go back to the backing graph; the path knows lengths.
        return mapping_from_length(defining_path.mapping(index));
    }
    
    string PathSubgraph::get_sequence(const handle_t& handle) const {
        size_t index = (size_t)get_id(handle) - 1;
        bool backward = get_is_reverse(handle);
        auto& pos = defining_path.mapping(index).position();
        // Get the backing handle to the node we are visiting, in the orientation we are visiting it.
        handle_t backing_handle = super->get_handle(pos.node_id(), pos.is_reverse() != backward);

        // Grab its full sequence in the correct orientation.
        string backing_sequence = get_sequence(backing_handle);

        if (index == 0 && pos.offset() != 0) {
            // We need to trim off the start of the backing node in its orientation along the path.

            size_t desired_length = get_length(handle);

            if (backward) {
                // Its orientation along the path is opposite our visit orientation.
                // Trim the end of the string.
                return backing_sequence.substr(0, desired_length); 
            } else {
                // Its orientation along the path is our visit orientation.
                // Trim the start of the string.
                return backing_sequence.substr(backing_sequence.size() - desired_length);
            }
        } else if (index == defining_path.mapping_size()) {
            // We may need to trim off the end of the backing node in its orientation along the path.

            size_t desired_length = get_length(handle);

            if (backward) {
                // Its orientation along the path is opposite our visit orientation.
                // Trim the start of the string.
                return backing_sequence.substr(backing_sequence.size() - desired_length);
            } else {
                // Its orientation along the path is our visit orientation.
                // Trim the end of the string.
                return backing_sequence.substr(0, desired_length); 
            }

        } else {
            // Just send the sequence as is
            return backing_sequence;
        }
    }
    
    bool PathSubgraph::follow_edges_impl(const handle_t& handle, bool go_left, const function<bool(const handle_t&)>& iteratee) const {
        // There's only ever 0 or 1 edges
        size_t index = (size_t)get_id(handle) - 1;
        bool backward = get_is_reverse(handle);

        if (index == 0 && (go_left != backward)) {
            // Hit left edge
            return true;
        }

        if (index == defining_path.mapping_size() - 1 && (go_left == backward)) {
            // Hit right edge
            return true;
        }

        // Otherwise we can go somewhere
        if (go_left == backward) {
            // Going forward in path
            index++;
        } else {
            index--;
        }

        // Go there and return the bool flag
        return iteratee(get_handle(index + 1, backward));
    }
    
    bool PathSubgraph::for_each_handle_impl(const function<bool(const handle_t&)>& iteratee, bool parallel) const {
        // TODO: implement parallel mode.
        // Paths should be short so we shouldn't really need it.
        for (size_t i = 0; i < defining_path.mapping_size(); i++) {
            // Try visiting each path visit
            if (!iteratee(get_handle(i + 1, false))) {
                // Stop early
                return false;
            }
        }
        return true;
    }
    
    size_t PathSubgraph::node_size() const {
        return defining_path.mapping_size();
    }
    
    id_t PathSubgraph::min_node_id() const {
        return 1;
    }
    
    id_t PathSubgraph::max_node_id() const {
        return defining_path.mapping_size();
    }

}

