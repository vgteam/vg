/**
 * \file indexed_vg.cpp
 * Implementation for the IndexedVG class, which provides a HandleGraph interface to an on-disk VG file.
 */

#include "indexed_vg.hpp"
#include "utility.hpp"
#include "vg/io/json2pb.h"

#include <handlegraph/util.hpp>

#include <atomic>

namespace vg {

using namespace std;

IndexedVG::IndexedVG(string graph_filename) : vg_filename(graph_filename), index(),
    cursor_streams(), cursor_pool(), cursor_pool_mutex(), group_cache(100), cache_mutex() {
    
    // Decide where the index ought to be stored
    string index_filename = vg_filename + ".vgi";
    
    ifstream index_in_stream(index_filename);
    if (index_in_stream.good()) {
        // We found the index, load it
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
    cerr << group_cache.size() << " cache entries" << endl;
    // TODO: Cache hit/miss counts from the LRUcache do not appear to be
    // correct (hits seem to be counted as misses). So we don't report them
    // here.
}

bool IndexedVG::has_node(id_t node_id) const {
    bool id_in_graph = false;
    find(node_id, [&](const CacheEntry& entry) -> bool {
            // For each relevant entry (which may just have some edges to the node we are looking for)
            auto found = entry.id_to_node_index.find(node_id);
            if (found != entry.id_to_node_index.end()) {
                // We found the node!
                id_in_graph = true;
                // Stop
                return false;
            }
        
            // Otherwise we don't have the node we want
            return true;
        });

    return id_in_graph;
}

// TODO: We ought to use some kind of handle packing that relates to file offsets for graph chunks contasining nodes.
// For now we just use the handlegraph::number_bool_packing and hit the index every time.

handle_t IndexedVG::get_handle(const id_t& node_id, bool is_reverse) const {
    return handlegraph::number_bool_packing::pack(node_id, is_reverse);
}

id_t IndexedVG::get_id(const handle_t& handle) const {
    return handlegraph::number_bool_packing::unpack_number(handle);
}

bool IndexedVG::get_is_reverse(const handle_t& handle) const {
    return handlegraph::number_bool_packing::unpack_bit(handle);
}

handle_t IndexedVG::flip(const handle_t& handle) const {
    return handlegraph::number_bool_packing::toggle_bit(handle);
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

bool IndexedVG::follow_edges_impl(const handle_t& handle, bool go_left, const function<bool(const handle_t&)>& iteratee) const {
    // TODO: implement stopping early in the backing index!
    
#ifdef debug
    cerr << "Following edges " << (go_left ? "left" : "right") << " from " << get_id(handle) << " orientation " << get_is_reverse(handle) << endl;
#endif
    
    if (go_left) {
        // Go right from our reverse version, and return flipped results
        return follow_edges(flip(handle), false, [&](const handle_t& other) -> bool {
            return iteratee(flip(other));
        });
    }
    
    // Now we only have to handle the going right case
    
    // If this is false, don't call the iteratee any more.
    bool keep_going = true;

    // Get the ID of the node we are looking for
    id_t id = get_id(handle);
    bool is_reverse = get_is_reverse(handle);
    
    find(id, [&](const CacheEntry& entry) -> bool {
        // For each CacheEntry that describes a graph group that may have edges touching the ID we are looking for
        
#ifdef debug
        cerr << "Relevant cache entry: " << &entry << endl;
#endif
        
        // Find the list of edge indices that touch this node ID
        auto found = entry.id_to_edge_indices.find(id);
        
        if (found == entry.id_to_edge_indices.end()) {
            // No relevant edges in this cache entry. Get the next potentially relevant one.
            return true;
        }
        
#ifdef debug
        cerr << "Entry has " << found->second.size() << " edge indices that touch node " << id << endl;
#endif
        
        for (auto edge_index : found->second) {
            // Look up each relevant edge in the graph
            auto& edge = entry.merged_group.edge(edge_index);
            
#ifdef debug
            cerr << "Consider edge #" << edge_index << " in cache entry graph" << endl;
#endif
            
            if (edge.from() == id) {
                // This edge touches us on its from end
                if (
                    (!edge.from_start() && !is_reverse) || // Normal down an edge case
                    (edge.from_start() && is_reverse) // We also end up going down the edge the same way from the edge's point of view
                ) {
#ifdef debug
                    cerr << "Follow edge " << pb2json(edge) << " from from end" << endl;
#endif
                    keep_going &= iteratee(get_handle(edge.to(), edge.to_end()));
                }
            }
            
            if (!keep_going) {
                return false;
            }
            
            if (edge.to() == id) {
                // This edge touches us on its to end
                if (
                    (edge.to_end() && !is_reverse) || // We read up this edge
                    (!edge.to_end() && is_reverse) // We also read up this edge
                ) {
#ifdef debug
                    cerr << "Follow edge " << pb2json(edge) << " from to end" << endl;
#endif
                    keep_going &= iteratee(get_handle(edge.from(), !edge.from_start()));
                }
            }
            
            if (!keep_going) {
                return false;
            }
        }
        
        return keep_going;
    });
    
    if (!keep_going) {
#ifdef debug
        cerr << "Stopped early!" << endl;
#endif
    }
        
    return keep_going;
}

bool IndexedVG::for_each_handle_impl(const function<bool(const handle_t&)>& iteratee, bool parallel) const {
    // We have to scan the whole graph for this
    
    int64_t group_vo = 0;
    atomic<bool> keep_going(true);
    
    while(keep_going) {
        // Look up the cache entry here
        bool still_in_file = with_cache_entry(group_vo, [&](const CacheEntry& entry) {
            
            if (parallel) {
                // Handle each node in the cache entry in a task
                #pragma omp parallel for
                for (size_t i = 0; i < entry.merged_group.node_size(); i++) {
                    // Show a handle for every node to the iteratee
                    // stopping is best effort in multithreaded mode; we don't try very hard.
                    if (keep_going) {
                        if (!iteratee(get_handle(entry.merged_group.node(i).id(), false))) {
                            keep_going = false;
                        }
                    }
                }
            } else {
                // Do it single threaded
                for (auto& node : entry.merged_group.node()) {
                    // Show a handle for every node to the iteratee
                    if (!iteratee(get_handle(node.id(), false))) {
                        keep_going = false;
                    }
                    
                    if (!keep_going) {
                        // If it is done, stop.
                        break;
                    }
                }
            }
            
            // Move on to the next cache entry for the next group
            group_vo = entry.next_group;
        });
        
        if (!still_in_file) {
            // We hit EOF
            break;
        }
    }
    
    return keep_going;
}

size_t IndexedVG::get_node_count() const {
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
    
    // This is the virtual offset of the serialized graph group we are considering
    int64_t group_vo = 0;
    
    while (min_node_id == numeric_limits<id_t>::max()) {
        // Get graph groups in order
        bool still_in_file = with_cache_entry(group_vo, [&](const CacheEntry& entry) {
            if (entry.merged_group.node_size() > 0) {
                // This graph has a first node
                min_node_id = entry.merged_group.node(0).id();
            }
            
            // Move on to the next cache entry for the next group
            group_vo = entry.next_group;
        });
        
        if (!still_in_file) {
            // We hit EOF
            break;
        }
    }
    
    return min_node_id;
}

id_t IndexedVG::max_node_id() const {

    // This is the max node ID observed
    id_t max_observed = 0;
   
    // Scan locally-forward but globally-backward through the vg file
    index.scan_backward([&](int64_t start_vo, int64_t past_end_vo) -> bool {
        
        int64_t group_vo = start_vo;
        
        while (group_vo < past_end_vo) {
            bool still_in_file = with_cache_entry(group_vo, [&](const CacheEntry& entry) {
                if (entry.merged_group.node_size() > 0) {
                    // We have nodes. The last one will have the highest ID.
                    max_observed = entry.merged_group.node(entry.merged_group.node_size() - 1).id();
                }
                // Move on to the next cache entry for the next group
                group_vo = entry.next_group;
            });
            
            if (!still_in_file) {
                // We hit EOF
                break;
            }
        }
        
        if (max_observed != 0) {
            // We found something. Stop going back towards the beginning.
            return false;
        } else {
            // We haven't seen any nodes yet. Keep searching
            return true;
        }
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
    // We will set this to false if the iteratee says to stop
    bool keep_going = true;
    
    index.find(id, [&](int64_t run_start_vo, int64_t run_past_end_vo) -> bool {
        // Loop over the index and get all the VO run ranges relevant to the given ID.
        
        // Scan through each run
        // Start at the run's start
        int64_t scan_vo = run_start_vo;
        
        while(keep_going && scan_vo < run_past_end_vo)  {
            
            // We are working on the group that starts at scan_vo
            
            // Go get it, unless scan_vo is at EOF
            bool still_in_file = with_cache_entry(scan_vo, [&](const CacheEntry& entry) {
                // Now we have a cache entry for the group we were looking for.
                // Show it to the iteratee.
                keep_going &= iteratee(entry);
                
                // Advance to the next group
                scan_vo = entry.next_group;
            });
            
            // We should never hit EOF when operating on ranges from the index.
            // If we do, the index is invalid.
            assert(still_in_file);
            
            // When scan_vo hits or passes the past-end for the range, we will be done with it
        }
        
        // Keep looking if the iteratee wants to, and get a new range.
        return keep_going;
    });
}

bool IndexedVG::with_cache_entry(int64_t group_vo, const function<void(const CacheEntry&)>& callback) const {

    if (group_vo == numeric_limits<int64_t>::max()) {
        // We got the EOF sentinel. We can't seek there.
        return false;
    }

    // This will point to the cache entry for the group when we find or make it.
    shared_ptr<CacheEntry> cache_entry;
    
    {
        lock_guard<mutex> lock(cache_mutex);
        // See if it is cached. Gets a pair of the item (if found) and a flag for whether it was found
        auto cache_pair = group_cache.retrieve(group_vo);
        if (cache_pair.second) {
            // We found it
            cache_entry = move(cache_pair.first);
        }
    }

    if (!cache_entry) {
        // If it wasn't found, load it up. We could synchronize to do
        // this with the cache lock held, to stop all threads banging
        // on the disk until one of them caches it. But we probably
        // want to allow simultaneous reads from disk overall.
        
        with_cursor([&](cursor_t& cursor) {
            // Try to get to the VO we are supposed to go to
            auto pre_seek_group = cursor.tell_group();
            if (!cursor.seek_group(group_vo)) {
                cerr << "error[vg::IndexedVG]: Could not seek from group pos " << pre_seek_group
                    << " to group pos " << group_vo << endl;
                cerr << "Current position: group " << cursor.tell_group()
                    << " has_current: " << cursor.has_current() << endl;
                assert(false);
            }
            
            if (cursor.has_current()) {
                // We seeked to a real thing and not EOF
            
                // Read the group into a cache entry
                cache_entry = shared_ptr<CacheEntry>(new CacheEntry(cursor));
            }
        });
        
        if (cache_entry) {
            // We actually found a valid group.
            
            lock_guard<mutex> lock(cache_mutex);
            // Save a copy of the shared pointer into the cache
            group_cache.put(group_vo, cache_entry);
        }
        
    }
    
    if (cache_entry) {
        // We aren't at EOF or anything, so call the callback
        callback(*cache_entry);
        return true;
    }
    
    // We didn't find it in the file.
    return false;
    
}

IndexedVG::CacheEntry::CacheEntry(cursor_t& cursor) {
    
    // We want to cache the group we are pointed at
    int64_t group_vo = cursor.tell_group();
    
    while (cursor.has_current() && cursor.tell_group() == group_vo) {
        // Merge the whole group together
        merged_group.MergeFrom(cursor.take());
    }

    if (!cursor.has_current()) {
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
        if (edge.to() != edge.from()) {
            // If it's not a self loop we need to point to it from both ends
            id_to_edge_indices[edge.to()].push_back(i);
        }
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

