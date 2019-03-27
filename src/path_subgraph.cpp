/**
 * \file path_subgraph.cpp: contains the implementation of PathSubgraph
 */


#include "path_subgraph.hpp"
#include "path.hpp"
#include <handlegraph/util.hpp>
#include <iostream>

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
        assert(node_id >= 1 && node_id <= defining_path.mapping_size());
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
        assert(index >= 0 && index < defining_path.mapping_size());
        bool backward = get_is_reverse(handle);
        auto& pos = defining_path.mapping(index).position();
        // Get the backing handle to the node we are visiting, in the orientation we are visiting it.
        handle_t backing_handle = super->get_handle(pos.node_id(), pos.is_reverse() != backward);

        // Grab its full sequence in the correct orientation.
        string backing_sequence = super->get_sequence(backing_handle);

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
        assert(index >= 0 && index < defining_path.mapping_size());
        bool backward = get_is_reverse(handle);

        if (index == 0 && ((go_left && !backward) || (!go_left && backward))) {
            // Hit left edge
            return true;
        }

        if (index == defining_path.mapping_size() - 1 && ((go_left && backward) || (!go_left && !backward))) {
            // Hit right edge
            return true;
        }

        // Otherwise we can go somewhere
        if ((go_left && backward) || (!go_left && !backward)) {
            // Going forward in path
            index++;
        } else {
            index--;
        }

        assert(index >= 0 && index < defining_path.mapping_size());

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

    Path PathSubgraph::translate_down(const Path& path_against_subgraph) const {
        Path translated;

        for (auto& subgraph_mapping : path_against_subgraph.mapping()) {
            // Translate each mapping
            Mapping* translated_mapping = translated.add_mapping();

            // Look up the defining Mapping we are visiting
            auto& defining_mapping = defining_path.mapping(subgraph_mapping.position().node_id() - 1);

            // TODO: simplify out repeated code here once we're sure each case is really correct
            if (defining_mapping.position().is_reverse() == false && subgraph_mapping.position().is_reverse() == false) {
                // We're in the forward orientation all the way through.
                // If there's an offset in the defining mapping, we need to add that to the offset in the subgraph mapping.
                // If the defining mapping has a short length, we don't care because we know the subgraph mapping won't be longer.
                translated_mapping->mutable_position()->set_node_id(defining_mapping.position().node_id());
                translated_mapping->mutable_position()->set_offset(defining_mapping.position().offset() + subgraph_mapping.position().offset());
                // The result will be forward
            } else if (defining_mapping.position().is_reverse() == false && subgraph_mapping.position().is_reverse() == true) {
                // We're in the forward orientation agaisnt the backing graph but the reverse orientation against the path.
                // Any shortness in the path mapping from length needs to be turned into an offset and added to the backing path offset.
                // Any offset in it will be ignored.
                size_t shortness = mapping_from_length(defining_mapping) - mapping_from_length(subgraph_mapping) - subgraph_mapping.position().offset();

                translated_mapping->mutable_position()->set_node_id(defining_mapping.position().node_id());
                translated_mapping->mutable_position()->set_offset(defining_mapping.position().offset() + shortness);
                // We come out backward
                translated_mapping->mutable_position()->set_is_reverse(true);
            } else if (defining_mapping.position().is_reverse() == true && subgraph_mapping.position().is_reverse() == false) {
                // We're in the reverse orientation against the backing graph, and the mapping to the path agrees with that.
                // We need to sum the offsets and ignore shortness
                translated_mapping->mutable_position()->set_node_id(defining_mapping.position().node_id());
                translated_mapping->mutable_position()->set_offset(defining_mapping.position().offset() + subgraph_mapping.position().offset());
                // And we will stay reverse
                translated_mapping->mutable_position()->set_is_reverse(true);
            } else {
                // We're in the reverse orientation in the backing graph, but then flip back against that.
                // We need to add shortness to offset
                size_t shortness = mapping_from_length(defining_mapping) - mapping_from_length(subgraph_mapping) - subgraph_mapping.position().offset();
                
                translated_mapping->mutable_position()->set_node_id(defining_mapping.position().node_id());
                translated_mapping->mutable_position()->set_offset(defining_mapping.position().offset() + shortness);
                // We come out in the forward orientation
                translated_mapping->mutable_position()->set_is_reverse(true);
            }

            // The edits always stay the same
            for (auto& edit : subgraph_mapping.edit()) {
                *translated_mapping->add_edit() = edit;
            }
        }

        return translated;
    }

}

