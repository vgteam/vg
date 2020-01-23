#ifndef VG_INDEXED_VG_HPP_INCLUDED
#define VG_INDEXED_VG_HPP_INCLUDED

/**
 * \file indexed_vg.hpp
 * Contains an implementation of a HandleGraph backed by a sorted, indexed .vg file
 */

#include <lru_cache.h>

#include <string>
#include <list>
#include <mutex>

#include "stream_index.hpp"
#include "handle.hpp"


namespace vg {

using namespace std;

/** Use a .vg file on disk with a .vgi index to provide random access to the
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
 *
 * Internally, we keep a pool of cursors into the backing graph file, and each
 * time we need to actually access the backing graph file we grab a cursor or
 * make one if we don't have a free one.
 *
 * Internally we also keep a least-recently-used cache of indexed
 * merged-together graph groups. The cache is keyed by group start VO. The
 * cache holds shared pointers to cache entries, so that one thread can be
 * evicting something from the cache while another is still working with it.
 */
class IndexedVG : public HandleGraph {

public:

    /// Open a .vg file. If the .vg has a .vg.vgi index, it wil be loaded. If
    /// not, an index will be generated and saved.
    IndexedVG(string graph_filename);
    
    // TODO: This gets implicitly deleted and generates warning because of the
    // StreamIndex member variable
    // We are moveable
    //IndexedVG(IndexedVG&& other) = default;
    
    // TODO: This gets implicitly deleted and generates warning because StreamIndex
    // member variable is not movable
    //IndexedVG& operator=(IndexedVG&& other) = default;
    
    void print_report() const;

private:
    // We are not copyable because we keep a pool of open files
    IndexedVG(const IndexedVG& other) = delete;
    
    // TODO: This gets implicitly deleted and generates warning because of the
    // StreamIndex member variable
    //IndexedVG& operator=(const IndexedVG& other) = delete;
    
public:

    ///////////////
    // Handle Graph Interface
    ///////////////

    /// Check if a node exists by ID
    virtual bool has_node(id_t node_id) const;

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
    virtual bool follow_edges_impl(const handle_t& handle, bool go_left, const function<bool(const handle_t&)>& iteratee) const;
    
    /// Loop over all the nodes in the graph in their local forward
    /// orientations, in their internal stored order. Stop if the iteratee
    /// returns false. Can be told to run in parallel, in which case stopping
    /// after a false return value is on a best-effort basis and iteration
    /// order is not defined.
    virtual bool for_each_handle_impl(const function<bool(const handle_t&)>& iteratee, bool parallel = false) const;
    
    /// Return the number of nodes in the graph
    virtual size_t get_node_count() const;
    
    /// Return the smallest ID in the graph, or some smaller number if the
    /// smallest ID is unavailable. Return value is unspecified if the graph is empty.
    virtual id_t min_node_id() const;
    
    /// Return the largest ID in the graph, or some larger number if the
    /// largest ID is unavailable. Return value is unspecified if the graph is empty.
    virtual id_t max_node_id() const;
    
protected:
    /// We store the graph filename, so we can have cursors to it created on demand.
    /// This is necessary to have e.g. random accesses to bits of the graph while looping over the graph as a whole.
    /// The downside is we lose BGZF block cacheing between different streams of access.
    string vg_filename;
    
    /// Index data about the vg file
    StreamIndex<Graph> index;
    
    /// Define the type we use for cursors into the backing file.
    using cursor_t = StreamIndex<Graph>::cursor_t;
    
    /// Get temporary ownership of a cursor to the backing vg file.
    void with_cursor(function<void(cursor_t&)> callback) const;
    
    /// Input streams referenced by cursors live in this list that grows forever
    mutable list<ifstream> cursor_streams;
    /// Cursors live in this free pool
    mutable list<unique_ptr<cursor_t>> cursor_pool;
    /// Access is protected by this mutex
    mutable mutex cursor_pool_mutex;
    
    /// Represents an entry in the cache for a parsed group of graphs.
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
        
        /// This is the virtual offset of the next group in the file.
        /// If this was the last group in the file, this is numeric_limits<int64_t>::max().
        int64_t next_group;
    };
    
    /// Wrapper around the index's find, with cacheing. Supports stopping
    /// early, but doesn't do internal filtering of chunks/runs where the node
    /// being queried is in a hole. Runs the iteratee on CacheEntry objects for
    /// the runs that might have info on the requested node, in order.
    /// Internally holds shared_ptr copies to the cache entries it is handing
    /// out references to. Users must do all everything they need the
    /// CacheEntry for within the callback as the reference may not be valid
    /// afterwards.
    void find(id_t id, const function<bool(const CacheEntry&)>& iteratee) const;
    
    /// Load or use the cached version of the CacheEntry for the given group
    /// start VO. If the EOF sentinel numeric_limits<int64_t>::max() is passed,
    /// the callback is not called and false is returned. (This is to enable
    /// easy looping to scan over CacheEntries.) Passing any other past-the-end
    /// VO is prohibited, and may produce an error. Handles locking the cache
    /// for updates and keeping the CacheEntry reference live while the
    /// callback is running.
    bool with_cache_entry(int64_t group_vo, const function<void(const CacheEntry&)>& callback) const;
    
    /// This is the cache that holds CacheEntries for groups we have already parsed and indexed.
    /// We can only access the cache from one thread at a time, but the shared pointers let us
    /// be working with the actual data in ther threads.
    mutable LRUCache<int64_t, shared_ptr<CacheEntry>> group_cache;
    /// The cache is protected with this mutex
    mutable mutex cache_mutex;
};

}


#endif
