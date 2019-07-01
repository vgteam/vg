/**
 * \file tree_subgraph.cpp
 * Contains the implementation of the TreeSubgraph.
 */


#include "tree_subgraph.hpp"
#include <handlegraph/util.hpp>
#include <iostream>

namespace vg {

using namespace std;


TreeSubgraph::TreeSubgraph(const HandleGraph* super, vector<pair<int64_t, handle_t>>&& tree, size_t root_trim) : super(super),
    tree(tree), root_trim(root_trim), children(tree.size()) {
    
    for (size_t i = 1; i < tree.size(); i++) {
        // Populate children by scanning the tree vector.
        // Don't scan the root because it can't be anyone's child.
        
        // Tell the parent of this tree node that it has us for a child.
        children.at(tree[i].first).push_back(i);
    }
}

vector<handle_t> TreeSubgraph::get_topological_order() const {
    vector<handle_t> to_return;
    
    for (size_t i = 0; i < tree.size(); i++) {
        // The tree is already in topological order, so we just use that order.
        to_return.push_back(handlegraph::number_bool_packing::pack(i, false));
    }
    
    return to_return;
}

handle_t TreeSubgraph::get_root() const {
    if (tree.empty()) {
        throw runtime_error("Tree is empty and has no root");
    }
    
    // Return the handle for the 0th element in the tree.
    return handlegraph::number_bool_packing::pack(0, false); 
}

bool TreeSubgraph::has_node(id_t node_id) const {
    return node_id > 0 && node_id <= tree.size();
}

handle_t TreeSubgraph::get_handle(const id_t& node_id, bool is_reverse) const {
    return handlegraph::number_bool_packing::pack(node_id - 1, is_reverse);
}

id_t TreeSubgraph::get_id(const handle_t& handle) const {
    return handlegraph::number_bool_packing::unpack_number(handle) + 1;
}

bool TreeSubgraph::get_is_reverse(const handle_t& handle) const {
    return handlegraph::number_bool_packing::unpack_bit(handle);
}

handle_t TreeSubgraph::flip(const handle_t& handle) const {
    return handlegraph::number_bool_packing::toggle_bit(handle);
}

size_t TreeSubgraph::get_length(const handle_t& handle) const {
    // Get the length in the backing graph
    size_t length = super->get_length(get_underlying_handle(handle));
    
    if (get_id(handle) == 1 && root_trim != 0) {
        // Trim the root as necessary
        length -= root_trim;
    }
    
    return length;
}

string TreeSubgraph::get_sequence(const handle_t& handle) const {
    // TODO: use get_subsequence to efficiently trim the root

    // Get the full backing sequence in the correct orientation to return
    string sequence = super->get_sequence(get_underlying_handle(handle));

    if (get_id(handle) == 1 && root_trim != 0) {
        // Trim the root
        
        if (get_is_reverse(handle)) {
            // We need to cut off the end
            sequence = sequence.substr(0, sequence.size() - root_trim);
        } else {
            // We need to cut off the start
            sequence = sequence.substr(root_trim);
        }
    }

    return sequence;
}

bool TreeSubgraph::follow_edges_impl(const handle_t& handle, bool go_left, const function<bool(const handle_t&)>& iteratee) const {
    // Work out our index in the tree vector
    size_t index = get_id(handle) - 1;
    // Work out if we want the parent or the child
    bool find_parent = (go_left != get_is_reverse(handle));
    // Work out if we need to flip the result when we get there
    bool flip_result = get_is_reverse(handle);
    
    if (find_parent) {
        if (index == 0) {
            // No parent of the root
            return true;
        }
        
        // Otherwise, just visit the one parent.
        size_t result_index = tree.at(index).first;
        return iteratee(get_handle(result_index + 1, flip_result));
    } else {
        bool keep_going = true;
        for (const size_t& result_index : children.at(index)) {
            // Go through all the children
            
            keep_going &= iteratee(get_handle(result_index + 1, flip_result));
            
            if (!keep_going) {
                break;
            }
        }
        return keep_going;
    }
   
}

bool TreeSubgraph::for_each_handle_impl(const function<bool(const handle_t&)>& iteratee, bool parallel) const {
    for (size_t i = 0; i < tree.size(); i++) {
        // For each index, turn it into the corresponding ID and get the handle and run on it.
        if (!iteratee(get_handle(i + 1, false))) {
            // Asked to stop early.
            return false;
        }
    }
    return true;
}

size_t TreeSubgraph::get_node_count() const {
    return tree.size();
}

id_t TreeSubgraph::min_node_id() const {
    return 1;
}

id_t TreeSubgraph::max_node_id() const {
    return tree.size();
}

handle_t TreeSubgraph::get_underlying_handle(const handle_t& handle) const {
    size_t index = get_id(handle) - 1;
    bool flip = get_is_reverse(handle);
    
    // Find the backing graph handle
    handle_t to_return = tree.at(index).second;
    
    if (flip) {
        // Flip it if necessary
        to_return = super->flip(to_return);
    }
    
    return to_return;
}

Path TreeSubgraph::translate_down(const Path& path_against_subgraph) const {
    // Copy the whole path
    Path translated = path_against_subgraph;
    
    for (size_t i = 0; i < translated.mapping_size(); i++) {
        // Get the handle in ourselves
        handle_t visited = get_handle(translated.mapping(i).position().node_id(), translated.mapping(i).position().is_reverse());
        
        // Translate it down
        handle_t underlying = get_underlying_handle(visited);
        
        // Put its ID and orientation in.
        translated.mutable_mapping(i)->mutable_position()->set_node_id(super->get_id(underlying));
        translated.mutable_mapping(i)->mutable_position()->set_is_reverse(super->get_is_reverse(underlying));
        
        if (get_id(visited) == 1 && !get_is_reverse(visited) && root_trim != 0) {
            // We're on the forward strand of the root, and the root needs trimming, so adjust the offset.
            // If we are at 0 on the trimmed node, we are at root_trim on the untrimmed node.
            translated.mutable_mapping(i)->mutable_position()->set_offset(translated.mapping(i).position().offset() + root_trim);
        }
    }
    
    return translated;
}

}

