#ifndef VG_INDEXED_VG_HPP_INCLUDED
#define VG_INDEXED_VG_HPP_INCLUDED

/**
 * \file indexed_vg.hpp
 * Contains an implementation of a HandleGraph backed by a sorted, indexed .vg file
 */
 
#include <string>
#include <mutex>

#include "stream_index.hpp"
#include "handle.hpp"

namespace vg {

using namespace std;

/**
 * Use a .vg file on disk with a .vgi index to provide random access to the
 * graph data without loading the entire graph into memory. Sort of a
 * compromise between an XG and a VG, except unlike either we don't need the
 * whole graph in memory.
 *
 * We require that all nodes in the graph appear in ID order within their
 * chunks, and that all chunks appear in ID order. So all nodes are in ID order
 * in the file.
 *
 * Cannot be copied since internally it contains a ProtobufIterator wrapping an
 * open file. Can only be moved.
 *
 * All operations are thread-safe to call. Internally we can't be seeking a
 * cursor off to another location in the middle of looping over a run of
 * matchung chunks, but we handle that ourselves.
 */
class IndexedVG : public HandleGraph {

public:

    /// Open a .vg file. If the .vg has a .vg.vgi index, it wil be loaded. If
    /// not, an index will be generated and saved.
    IndexedVG(string graph_filename);
    
    // We are moveable
    IndexedVG(IndexedVG&& other) = default;
    IndexedVG& operator=(IndexedVG&& other) = default;
    
    void print_report() const;

private:
    // We are not copyable
    // TODO: Make us copyable or add a cursor pool to justify non-copyability
    IndexedVG(const IndexedVG& other) = delete;
    IndexedVG& operator=(const IndexedVG& other) = delete;
    
public:

    ///////////////
    // Handle Graph Interface
    ///////////////
    
    /// Look up the handle for the node with the given ID in the given orientation
    virtual handle_t get_handle(const id_t& node_id, bool is_reverse = false) const;
    
    /// Get the ID from a handle
    virtual id_t get_id(const handle_t& handle) const;
    
    /// Get the orientation of a handle
    virtual bool get_is_reverse(const handle_t& handle) const;
    
    /// Invert the orientation of a handle (potentially without getting its ID)
    virtual handle_t flip(const handle_t& handle) const;
    
    /// Get the length of a node
    virtual size_t get_length(const handle_t& handle) const;
    
    /// Get the sequence of a node, presented in the handle's local forward
    /// orientation.
    virtual string get_sequence(const handle_t& handle) const;
    
    /// Loop over all the handles to next/previous (right/left) nodes. Passes
    /// them to a callback which returns false to stop iterating and true to
    /// continue. Returns true if we finished and false if we stopped early.
    virtual bool follow_edges(const handle_t& handle, bool go_left, const function<bool(const handle_t&)>& iteratee) const;
    
    /// Loop over all the nodes in the graph in their local forward
    /// orientations, in their internal stored order. Stop if the iteratee
    /// returns false. Can be told to run in parallel, in which case stopping
    /// after a false return value is on a best-effort basis and iteration
    /// order is not defined.
    virtual void for_each_handle(const function<bool(const handle_t&)>& iteratee, bool parallel = false) const;
    
    /// Return the number of nodes in the graph
    virtual size_t node_size() const;
    
    /// Return the smallest ID in the graph, or some smaller number if the
    /// smallest ID is unavailable. Return value is unspecified if the graph is empty.
    virtual id_t min_node_id() const;
    
    /// Return the largest ID in the graph, or some larger number if the
    /// largest ID is unavailable. Return value is unspecified if the graph is empty.
    virtual id_t max_node_id() const;
    
    using HandleGraph::follow_edges;
    using HandleGraph::for_each_handle;
    using HandleGraph::get_handle;
    
    
    
protected:
    /// We store the graph filename, so we can have cursors to it created on demand.
    /// This is necessary to have e.g. random accesses to bits of the graph while looping over the graph as a whole.
    /// The downside is we lose BGZF block cacheing between different streams of access.
    string vg_filename;
    
    /// Index data about the vg file
    StreamIndex<Graph> index;
    
    /// Defien the type we use for cursors into the backing file.
    using cursor_t = StreamIndex<Graph>::cursor_t;
    
    /// Get temporary ownership of a cursor to the backing vg file.
    void with_cursor(function<void(cursor_t&)> callback) const;
    
    /// Input streams referenced by cursors live in this list that grows forever
    mutable list<ifstream> cursor_streams;
    /// Cursors live in this free pool
    mutable list<unique_ptr<cursor_t>> cursor_pool;
    /// Access is protected by this mutex
    mutable mutex cursor_pool_mutex;
    
    /// Represents an entry in the cache for a parsed run of graphs.
    /// Has its own indexes and the virtual offset of the next group.
    struct CacheEntry {
        /// Make a cache entry for a group by reading a cursor at that group's start
        CacheEntry(cursor_t& to_read);
    
        // Can be moved
        CacheEntry(CacheEntry&& other) = default;
        CacheEntry& operator=(CacheEntry&& other) = default;
        
        /// Pull out the subgraph for the given node
        Graph query(const id_t& id) const;
    
        /// All the graphs get merged into this one
        Graph merged_group;
        
        /// This maps from node ID to index in the merged graph.
        unordered_map<id_t, size_t> id_to_node_index;
        
        /// This maps from node ID to indexes of touched edges in the merged graph.
        unordered_map<id_t, vector<size_t>> id_to_edge_indices;
        
        // TODO: Path visits
        
        /// This is the virtual offset of the next group
        int64_t next_group;
    };
    
    /// Wrapper around the index's find, with cacheing. Supports stopping
    /// early, but doesn't do internal filtering of chunks/runs where the node
    /// being queried is in a hole. Runs the iteratee on CacheEntry objects for
    /// the runs that might have info on the requested node, in order.
    void find(id_t id, const function<bool(const CacheEntry&)>& iteratee) const;
    
    /// Cache that stores deserialized and indexed Graph object lists by group
    /// virtual offset. Also stores the virtual offset for the next group, or
    /// the max value for EOF, so we can do groups not in here.
    mutable unordered_map<int64_t, CacheEntry> graph_cache;
    /// This maps from node ID to the group it lives in in the graph cache
    mutable unordered_map<id_t, int64_t> node_group_cache;
    /// The cache is protected with this mutex
    mutable mutex cache_mutex;
    
    mutable size_t cache_hits;
    mutable size_t cache_misses;
};

}


#endif
