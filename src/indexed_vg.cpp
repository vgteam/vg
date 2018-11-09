/**
 * \file indexed_vg.cpp
 * Implementation for the IndexedVG class, which provides a HandleGraph interface to an on-disk VG file.
 */

#include "indexed_vg.hpp"
#include "utility.hpp"

namespace vg {

using namespace std;

IndexedVG::IndexedVG(string graph_filename) : vg_filename(graph_filename), index() {
    
    // Decide where the index ought to be stored
    string index_filename = vg_filename + ".vgi";
    
    ifstream index_in_stream(index_filename);
    if (index_in_stream.good()) {
        // We found the idnex, load it
        index.load(index_in_stream);
    } else {
        // We need to build the index
        
        // Get the file to write the index to
        ofstream index_out_stream(vg_filename);
        if (!index_out_stream.good()) {
            // We couldn't load the index and we can't save it
            throw runtime_error("Could not open index file " + vg_filename + " for reading or writing");
        }
        
        // TODO: Show progress as we do this?
        with_cursor([&](cursor_t& cursor) {
            // Get a cursor to the start of the file
            assert(cursor.seek_group(0));
            
            // Compute the index
            index.index(cursor);
            
            // Save the index
            index.save(index_out_stream);
        });
    }
    
}

// TODO: We ought to use some kind of handle packing that relates to file offsets for graph chunks contasining nodes.
// For now we just use the EasyHandlePacking and hit the index every time.

handle_t IndexedVG::get_handle(const id_t& node_id, bool is_reverse) const {
    return EasyHandlePacking::pack(node_id, is_reverse);
}

id_t IndexedVG::get_id(const handle_t& handle) const {
    return EasyHandlePacking::unpack_number(handle);
}

bool IndexedVG::get_is_reverse(const handle_t& handle) const {
    return EasyHandlePacking::unpack_bit(handle);
}

handle_t IndexedVG::flip(const handle_t& handle) const {
    return EasyHandlePacking::toggle_bit(handle);
}

size_t IndexedVG::get_length(const handle_t& handle) const {
    // We don't have a more efficient way to get the length than loading the sequence
    return get_sequence(handle).size();
}

string IndexedVG::get_sequence(const handle_t& handle) const {
    
    // Get the ID of the node we are looking for
    id_t id = get_id(handle);
    
    // We will pull the sequence out into this string.
    string found_sequence;
    
    with_cursor([&](cursor_t& cursor) {
        // Get ahold of the vg file cursor
        
        index.find(cursor, id, [&](const Graph& graph) {
            // For each relevant Graph (which may just have some edges to the node we are looking for)
            if (graph.node_size() > 0) {
                // We only care if it actually has nodes
                
                if (graph.node(0).id() > id) {
                    // This graph is out of range (too late)
                    return;
                }
                
                if (graph.node(graph.node_size() - 1).id() < id) {
                    // This graph is out of range (too early)
                    return;
                }
                
                // Otherwise binary search for the actual node
                size_t min_index = 0;
                size_t past_max_index = graph.node_size();
                
                while(true) {
                    size_t middle_index = (past_max_index + min_index) / 2;
                    
                    if (middle_index >= graph.node_size()) {
                        // Node was not found
                        throw runtime_error("Node " + to_string(id) + " not found in expected graph chunk! Is graph sorted correctly?");
                    }
                    
                    // Grab the actual node
                    auto& middle_node = graph.node(middle_index);
                    
                    if (middle_node.id() == id) {
                        // We found it
                        found_sequence = middle_node.sequence();
                        break;
                    } else if (middle_node.id() > id) {
                        // Only look left of here
                        past_max_index = middle_index;
                    } else {
                        // We must have found a node with too small an ID
                        // Only look right of here
                        min_index = middle_index + 1;
                    }
                    
                }
                
            }
        });
    });
    
    if (get_is_reverse(handle)) {
        // Reverse complement the sequence if necessary
        found_sequence = reverse_complement(found_sequence);
    }
    
    return found_sequence;
}

bool IndexedVG::follow_edges(const handle_t& handle, bool go_left, const function<bool(const handle_t&)>& iteratee) const {
    // TODO: implement stopping early in the backing index!
    
    if (go_left) {
        // Go right from our reverse version, and return flipped results
        return follow_edges(flip(handle), false, [&](const handle_t& other) -> bool {
            return iteratee(flip(other));
        });
    }
    
    // Now we only have to handle the going right case
    
    // If this is false, don't call the iteratee any more.
    // TODO: Pass the early stop on to the index scan.
    bool keep_going = true;
    
    with_cursor([&](cursor_t& cursor) {
        // Get ahold of the vg file cursor
        
        // Get the ID of the node we are looking for
        id_t id = get_id(handle);
        bool is_reverse = get_is_reverse(handle);
        
        index.find(cursor, id, [&](const Graph& graph) {
            if (!keep_going) {
                return;
            }
            for (auto& edge : graph.edge()) {
                if (edge.from() == id) {
                    if (
                        (!edge.from_start() && !is_reverse) || // Normal down an edge case
                        (edge.from_start() && is_reverse) // We also end up going down the edge the same way from the edge's point of view
                    ) {
                        keep_going &= iteratee(get_handle(edge.to(), edge.to_end()));
                    }
                }
                
                if (!keep_going) {
                    return;
                }
                
                if (edge.to() == id) {
                    if (
                        (edge.to_end() && !is_reverse) || // We read up this edge
                        (!edge.to_end() && is_reverse) // We also read up this edge
                    ) {
                        keep_going &= iteratee(get_handle(edge.from(), !edge.from_start()));
                    }
                }
                
                if (!keep_going) {
                    return;
                }
            }
        });
    });
    
    return keep_going;
}

void IndexedVG::for_each_handle(const function<bool(const handle_t&)>& iteratee, bool parallel) const {
    // We have to scan the whole graph for this
    // TODO: This ought to populate a cache for when we actually want to get at the nodes we find
    // TODO: Implement parallel mode
    
    // No need to query the index; we can start at the beginning
    
    // This holds if we should keep going or not
    bool keep_going = true;
    
    with_cursor([&](cursor_t& cursor) {
        // Get a cursor to the beginning
        assert(cursor.seek_group(0));
        
        while (keep_going && cursor.has_next()) {
            // Loop over each graph chunk
            Graph graph = move(cursor.take());
            
            for (auto& node : graph.node()) {
                // Loop over all the nodes in each chunk 
                if (!iteratee(get_handle(node.id(), false))) {
                    // After the iteratee looks at it, if it didn't like it, stop
                    keep_going = false;
                    break;
                }
            }
        }
    });
}

size_t IndexedVG::node_size() const {
    // TODO: Add total distinct node count to the index or cache it or something.
    // Right now we just scan.
    
    size_t count = 0;
    for_each_handle([&](const handle_t& ignored) {
        count++;
    });
    return count;
}

id_t IndexedVG::min_node_id() const {
    // We can just seek to the start and get the first node in the first chunk with nodes.
    
    // This holds the first real node ID we find
    id_t min_node_id = numeric_limits<id_t>::max();
    
    with_cursor([&](cursor_t& cursor) {
        // Get a cursor to the beginning
        assert(cursor.seek_group(0));
    
        while (min_node_id == numeric_limits<id_t>::max() && cursor.has_next()) {
            // Loop over each graph chunk until we find a real node
            Graph graph = move(cursor.take());
            
            for (auto& node : graph.node()) {
                // Remember the first node we come across
                min_node_id = node.id();
                break;
            }
        }
    });
    
    return min_node_id;
}

id_t IndexedVG::max_node_id() const {

    // This is the max node ID observed
    id_t max_observed = 0;
    
    with_cursor([&](cursor_t& cursor) {
        // Grab a cursor
        
        // Scan locally-forward but globally-backward through the vg file
        index.scan_backward([&](int64_t start_vo, int64_t past_end_vo) -> bool {
            // For each start point in backward order, seek to it
            assert(cursor.seek_group(start_vo));
            
            while(cursor.has_next() && cursor.tell_group() < past_end_vo) {
                // For each graph forward from there
                Graph graph = move(cursor.take());
                
                if (graph.node_size() > 0) {
                    // We have nodes. The last one will have the highest ID.
                    max_observed = graph.node(graph.node_size() - 1).id();
                }
                
                // We have to scan until the end of each group we visit
            }
            
            if (max_observed != 0) {
                // We found something. Stop going back towards the beginning.
                return false;
            } else {
                // We haven't seen any nodes yet. Keep searching
                return true;
            }
        });
    });
    
    return max_observed;
}

auto IndexedVG::with_cursor(function<void(cursor_t&)> callback) const -> void {
    // TODO: Use a cursor pool
    
    // Open the file
    ifstream vg_stream(vg_filename);
    assert(vg_stream.good());
    
    // Set up the cursor
    cursor_t cursor(vg_stream);
    
    // Let the callback use it
    callback(cursor);
}

    
}

