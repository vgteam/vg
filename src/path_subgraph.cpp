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
        bool result = (node_id > 0 && node_id <= defining_path.mapping_size());
#ifdef debug
        cerr << "Have node " << node_id << ": " << result << endl;
#endif
        return result;
    }
    
    vector<handle_t> PathSubgraph::get_topological_order() const {
        vector<handle_t> order;
        order.reserve(defining_path.mapping_size());
        for (id_t i = 1; i <= defining_path.mapping_size(); i++) {
            // Make one handle per node in the path
            order.push_back(get_handle(i, false));
        }
        
#ifdef debug
        cerr << "Path: " << pb2json(defining_path) << endl;
        cerr << "Order:";
        for (auto& h : order) {
            cerr << " " << get_id(h) << (get_is_reverse(h) ? "-" : "+");
        }
        cerr << endl;
#endif
        
        return order;
    }
    
    handle_t PathSubgraph::get_handle(const id_t& node_id, bool is_reverse) const {
        assert(node_id >= 1 && node_id <= defining_path.mapping_size());
        handle_t handle = handlegraph::number_bool_packing::pack(node_id, is_reverse);
        assert(get_id(handle) == node_id);
        assert(get_is_reverse(handle) == is_reverse);
        return handle;
    }
    
    id_t PathSubgraph::get_id(const handle_t& handle) const {
        return handlegraph::number_bool_packing::unpack_number(handle);
    }
    
    bool PathSubgraph::get_is_reverse(const handle_t& handle) const {
        return handlegraph::number_bool_packing::unpack_bit(handle);
    }
    
    handle_t PathSubgraph::flip(const handle_t& handle) const {
        handle_t flipped = handlegraph::number_bool_packing::toggle_bit(handle);
        assert(get_is_reverse(flipped) != get_is_reverse(handle));
        assert(get_id(flipped) == get_id(handle));
        return flipped;
    }
    
    size_t PathSubgraph::get_length(const handle_t& handle) const {
        size_t index = (size_t)get_id(handle) - 1;
        // No need to go back to the backing graph; the path knows lengths.
        return mapping_from_length(defining_path.mapping(index));
    }
    
    string PathSubgraph::get_sequence(const handle_t& handle) const {
        // Find the backing node in its local forward orientation
        size_t index = (size_t)get_id(handle) - 1;
        assert(index >= 0 && index < defining_path.mapping_size());
        auto& pos = defining_path.mapping(index).position();
        handle_t backing_handle = super->get_handle(pos.node_id(), false);
        
        // Get its sequence in its local forward orientation
        string backing_sequence = super->get_sequence(backing_handle);
        
        // Work out what range of that sequence we want.
        size_t wanted_length = get_length(handle);
        size_t backing_first = 0;
        
#ifdef debug
        cerr << "Start selecting " << wanted_length << " bp starting at " << backing_first << endl;
#endif
        
        // For every offset, even 0
            
        // Work out whether we should do it from the
        // start or end of the backing sequence in its local forward orientation.
        // If the path visits the node forward, we cut from the start.
        // Otherwise, we cut from the end.
        bool cut_from_start = !pos.is_reverse();
        
        // Reposition the window
        // accordingly.
        size_t budge;
        if (cut_from_start) {
            // Account for the space at the start of the node consumed by the offset
            budge = pos.offset();
        } else {
            // Leave only the space at the end of the node consumed by the offset.
            // Budge by all the unwanted bases not consumed by the offset.
            budge = backing_sequence.size() - wanted_length - pos.offset();
        }
        
#ifdef debug
        cerr << "Budge by " << budge << endl;
#endif
        
        backing_first += budge;
        
#ifdef debug
        cerr << "End selecting " << wanted_length << " bp starting at " << backing_first << endl;
#endif
        
        // Pull out and reverse complement if necessary
        string wanted_sequence = backing_sequence.substr(backing_first, wanted_length);
        if (get_is_reverse(handle) != pos.is_reverse()) {
            // If we reverse the backing sequence only once, flip it.
            wanted_sequence = reverse_complement(wanted_sequence);
#ifdef debug
            cerr << "Flip it" << endl;
#endif
        }
        
#ifdef debug
        cerr << "Mapping " << pb2json(defining_path.mapping(index)) << " on sequence "
            << backing_sequence << " visited " << (get_is_reverse(handle) ? "rev" : "fwd")
            << " produces " << wanted_sequence << endl;
#endif
        
        // Return it
        return wanted_sequence;
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
    
    size_t PathSubgraph::get_node_count() const {
        return defining_path.mapping_size();
    }
    
    id_t PathSubgraph::min_node_id() const {
        return 1;
    }
    
    id_t PathSubgraph::max_node_id() const {
        return defining_path.mapping_size();
    }
    
    handle_t PathSubgraph::get_underlying_handle(const handle_t& handle) const {
        // Look up the defining Mapping we are visiting
        auto& defining_mapping = defining_path.mapping(get_id(handle) - 1);
        
        // Get the handle corresponding to this mapping in our path.
        return super->get_handle(defining_mapping.position().node_id(), defining_mapping.position().is_reverse());
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

