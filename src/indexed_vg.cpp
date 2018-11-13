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
        ofstream index_out_stream(index_filename);
        if (!index_out_stream.good()) {
            // We couldn't load the index and we can't save it
            throw runtime_error("Could not open index file " + index_filename + " for reading or writing");
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

void IndexedVG::print_report() const {
    cerr << cursor_streams.size() << " cursors outstanding, " << cursor_pool.size() << " cursors free" << endl;
    cerr << graph_cache.size() << " cache entries" << endl;
    cerr << cache_hits << " cache hits, " << cache_misses << " cache misses" << endl;
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
        
    find(id, [&](const CacheEntry& entry) -> bool {
        // For each relevant entry (which may just have some edges to the node we are looking for)
        auto found = entry.id_to_node_index.find(id);
        if (found != entry.id_to_node_index.end()) {
            // We found the node!
            // Copy out its sequence
            found_sequence = entry.merged_group.node(found->second).sequence();
            // Stop
            return false;
        }
        
        // Otherwise we don't have the node we want
        return true;
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

void IndexedVG::with_cursor(function<void(cursor_t&)> callback) const {
    // We'll fill this in with a cursor from the pool if we can get one, and a new one otherwise
    unique_ptr<cursor_t> obtained_cursor;
    {
        // Get ahold of the pool
        lock_guard<mutex> lock(cursor_pool_mutex);
        if (!cursor_pool.empty()) {
            // Grab a cursor from the pool
            obtained_cursor = move(cursor_pool.front());
            cursor_pool.pop_front();
        } else {
            // Open a new file stream
            cursor_streams.emplace_back(vg_filename);
            assert(cursor_streams.back().good());
            
            // Make a cursor around it
            obtained_cursor = unique_ptr<cursor_t>(new cursor_t(cursor_streams.back()));
        }
    }
    
    // Let the callback use it
    callback(*obtained_cursor.get());
    
    {
        // Get ahold of the pool
        lock_guard<mutex> lock(cursor_pool_mutex);
        
        // Put the cursor back in the pool
        cursor_pool.emplace_back(move(obtained_cursor));
    }
    
    // TODO: Does this moving unique_ptrs make sense or should we copy around indexes or real pointers
}

void IndexedVG::find(id_t id, const function<bool(const CacheEntry&)>& iteratee) const {
    
    /*auto group_found = node_group_cache.find(id);
    if (group_found != node_group_cache.end()) {
        // There's just the one entry for this node.
        // TODO: edges
        iteratee(graph_cache.at(group_found->second));
        return;
    }*/

    index.find(id, [&](int64_t run_start_vo, int64_t run_past_end_vo) -> bool {
        // Loop over the index and get all the VO run ranges
        
        // We will set this to false if the iteratee says to stop
        bool keep_going = true;
        
        // Scan through each run
        // Start at the run's start
        int64_t scan_vo = run_start_vo;
        
        while(keep_going && scan_vo < run_past_end_vo)  {
            // Get the cache entry for this run
            decltype(graph_cache)::const_iterator found;
            {
                lock_guard<mutex> lock(cache_mutex);
                found = graph_cache.find(scan_vo);
            }
            
            if (found != graph_cache.end()) {
                // If we have a cached group for this VO, use it
                // TODO: support the cached group being removed from the cache
                
                cache_hits++;
                
                // Iterate over the subgraph for the node from the cached run
                keep_going &= iteratee(found->second);
                
                if (!keep_going) {
                    break;
                }
                
                // Advance to the next group
                scan_vo = found->second.next_group;
            } else {
                // Otherwise, we need to actually access the file
                
                cache_misses++;
                
                int64_t cache_line_vo = scan_vo;
                
                // We're going to have a cache entry
                CacheEntry* cache_line = nullptr;
                 
                with_cursor([&](cursor_t& cursor) {
                    // Does nothing if we are already in the right place.
                    assert(cursor.seek_group(scan_vo));
                    
                    // Read the group into the cache
                    cache_line = new CacheEntry(cursor);
                    
                });
                
                // Iterate over the subgraph for the node from the cached run
                keep_going &= iteratee(*cache_line);
                
                // Remember where to look for the next group
                scan_vo = cache_line->next_group;
                
                /*for (auto& node : cache_line->merged_group.node()) {
                    node_group_cache[node.id()] = cache_line_vo;
                }*/
                
                // Save back to the cache
                {
                    lock_guard<mutex> lock(cache_mutex);
                    graph_cache.emplace(cache_line_vo, move(*cache_line));
                }
                
                delete cache_line;
                
                if (!keep_going) {
                    break;
                }
            }
            
            // When scan_vo hits or passes the past-end for the range, we will be done with it
        }
        
        // Keep looking if the iteratee wants to
        return keep_going;
    });
}

IndexedVG::CacheEntry::CacheEntry(cursor_t& cursor) {
    
    int64_t group_vo = cursor.tell_group();
    
    while (cursor.has_next() && cursor.tell_group() == group_vo) {
        // Merge the whole group together
        merged_group.MergeFrom(cursor.take());
    }

    if (!cursor.has_next()) {
        // We hit EOF
        next_group = numeric_limits<int64_t>::max();
    } else {
        // We found another group.
        // TODO: Can we avoid deserializing its first chunk?
        next_group = cursor.tell_group();
    }
    
    // Compute the indexes
    for (size_t i = 0; i < merged_group.node_size(); i++) {
        // Record the index of every node
        id_to_node_index[merged_group.node(i).id()] = i;
    }
    
    for (size_t i = 0; i < merged_group.edge_size(); i++) {
        // And of every edge by end node IDs
        auto& edge = merged_group.edge(i);
        id_to_edge_indices[edge.from()].push_back(i);
        id_to_edge_indices[edge.to()].push_back(i);
    }
}

Graph IndexedVG::CacheEntry::query(const id_t& id) const {
    Graph to_return;
    
    auto node_found = id_to_node_index.find(id);
    if (node_found != id_to_node_index.end()) {
        // We have the node in question, so send it
        *to_return.add_node() = merged_group.node(node_found->second);
    }
    
    auto edges_found = id_to_edge_indices.find(id);
    if (edges_found != id_to_edge_indices.end()) {
        // We have edges on it, so send them
        for (auto& edge_index : edges_found->second) {
            *to_return.add_edge() = merged_group.edge(edge_index);
        }
    }
    
    // TODO: Path visits
    
    return to_return;
}

    
}

