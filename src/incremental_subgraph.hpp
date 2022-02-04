/** \file
 * incremental_subgraph.hpp: defines a subgraph that is extracted from the parent
 * graph on an as-needed basis
 */
#ifndef VG_INCREMENTAL_SUBGRAPH_HPP_INCLUDED
#define VG_INCREMENTAL_SUBGRAPH_HPP_INCLUDED

#include <deque>

#include "handle.hpp"

namespace vg {

using namespace std;

/**
 * A subgraph that is extracted, made single-stranded, and DAG-fied online on an as-needed
 * basis from the parent graph. It is restricted to subgraphs that extend from a single position
 * in the graph in one direction.
 */
class IncrementalSubgraph : public ExpandingOverlayGraph {
public:
    
    /// Initialize
    IncrementalSubgraph(const HandleGraph& graph,
                        const pos_t& starting_position,
                        bool extract_left,
                        int64_t max_distance = numeric_limits<int64_t>::max());
    
    /// Default constructor -- not actually functional
    IncrementalSubgraph() = default;
    
    /// Default destructor
    ~IncrementalSubgraph();
    
    //////////////////////////
    /// Specialized interface
    //////////////////////////
    
    // TODO: prune method that removes outgoing edges into the frontier
    // from an extracted node
    
    /// True if there are additional nodes
    bool is_extendable() const;
    
    /// Extract an additional node
    handle_t extend();
    
    /// The order of a node in a topological order of the extracted graph
    size_t order_of(const handle_t& handle) const;
    
    /// The node at a given position in the topological order
    handle_t handle_at_order(size_t i) const;
    
    /// The minimum distance from the start position
    int64_t min_distance_from_start(const handle_t& handle) const;
    
    /// The maximum distance from the start position
    int64_t max_distance_from_start(const handle_t& handle) const;
    
    /// Are we doing extraction left or right?
    bool extracting_left() const;
    
    //////////////////////////
    /// HandleGraph interface
    //////////////////////////
    
    /// Method to check if a node exists by ID
    bool has_node(id_t node_id) const;
    
    /// Look up the handle for the node with the given ID in the given orientation
    handle_t get_handle(const id_t& node_id, bool is_reverse = false) const;
    
    /// Get the ID from a handle
    id_t get_id(const handle_t& handle) const;
    
    /// Get the orientation of a handle
    bool get_is_reverse(const handle_t& handle) const;
    
    /// Invert the orientation of a handle (potentially without getting its ID)
    handle_t flip(const handle_t& handle) const;
    
    /// Get the length of a node
    size_t get_length(const handle_t& handle) const;
    
    /// Get the sequence of a node, presented in the handle's local forward
    /// orientation.
    string get_sequence(const handle_t& handle) const;
    
    /// Loop over all the handles to next/previous (right/left) nodes. Passes
    /// them to a callback which returns false to stop iterating and true to
    /// continue. Returns true if we finished and false if we stopped early.
    bool follow_edges_impl(const handle_t& handle, bool go_left, const function<bool(const handle_t&)>& iteratee) const;
    
    /// Loop over all the nodes in the graph in their local forward
    /// orientations, in their internal stored order. Stop if the iteratee
    /// returns false. Can be told to run in parallel, in which case stopping
    /// after a false return value is on a best-effort basis and iteration
    /// order is not defined.
    bool for_each_handle_impl(const function<bool(const handle_t&)>& iteratee, bool parallel = false) const;
    
    /// Return the number of nodes in the graph
    size_t get_node_count() const;
    
    /// Return the smallest ID in the graph, or some smaller number if the
    /// smallest ID is unavailable. Return value is unspecified if the graph is empty.
    id_t min_node_id() const;
    
    /// Return the largest ID in the graph, or some larger number if the
    /// largest ID is unavailable. Return value is unspecified if the graph is empty.
    id_t max_node_id() const;
    
    ///////////////////////////////////
    /// Optional HandleGraph interface
    ///////////////////////////////////
    
    /// Get the number of edges on the right (go_left = false) or left (go_left
    /// = true) side of the given handle. The default implementation is O(n) in
    /// the number of edges returned, but graph implementations that track this
    /// information more efficiently can override this method.
    size_t get_degree(const handle_t& handle, bool go_left) const;
    
    /// Return the total number of edges in the graph. If not overridden,
    /// counts them all in linear time.
    size_t get_edge_count() const;
    
    /// Returns one base of a handle's sequence, in the orientation of the
    /// handle.
    char get_base(const handle_t& handle, size_t index) const;
    
    /// Returns a substring of a handle's sequence, in the orientation of the
    /// handle. If the indicated substring would extend beyond the end of the
    /// handle's sequence, the return value is truncated to the sequence's end.
    /// By default O(n) in the size of the handle's sequence, but can be overriden.
    string get_subsequence(const handle_t& handle, size_t index, size_t size) const;
    
    ///////////////////////////////////
    /// ExpandingOverlayGraph interface
    ///////////////////////////////////
    
    /**
     * Returns the handle in the underlying graph that corresponds to a handle in the
     * overlay
     */
    handle_t get_underlying_handle(const handle_t& handle) const;
    
private:
    
    // comparator for the frontier, order first by number of unseen incoming edges
    // and then by distance and break ties arbitrarily based on handle values
    // TODO: rearrange the tuples to make this less awkward
    struct FCmp {
        inline bool operator()(const tuple<int64_t, handle_t, unordered_set<handle_t>*, vector<size_t>*>& a,
                               const tuple<int64_t, handle_t, unordered_set<handle_t>*, vector<size_t>*>& b) const {
            return (get<2>(a)->size() < get<2>(b)->size() ||
                    (get<2>(a)->size() == get<2>(b)->size() && a < b));
        }
    };
    
    pair<size_t, size_t> underlying_interval(size_t i) const;
    
    /// direction we're extracting from the start pos
    bool extract_left;
    
    /// farthest distance we will travel from the start pos
    int64_t max_distance;
    
    
    /// records of (underlying handle, left edges, right edges, min distance, max distance)
    vector<tuple<handle_t, vector<size_t>, vector<size_t>, int64_t, int64_t>> extracted;
    /// the number of times a handle has been extracted
    unordered_map<handle_t, size_t> num_extracted;
    
    /// records of (unseen edges going into, seen edges going into, distance, node).
    /// serves as an updateable priority queue for nodes that are adjacent to the extracted nodes
    /// the container classes are created on the heap so that we can do a remove-modify-replace
    /// update to frontier entries without deep-copying the containers
    set<tuple<int64_t, handle_t, unordered_set<handle_t>*, vector<size_t>*>, FCmp> frontier;
    /// provides random access into the frontier by handle
    unordered_map<handle_t, deque<decltype(frontier)::iterator>> frontier_index;
    
    /// The underlying graph
    const HandleGraph* graph = nullptr;
};

}

#endif // VG_INCREMENTAL_SUBGRAPH_HPP_INCLUDED
