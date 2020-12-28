/** \file
 * chopped_graph.hpp: defines a handle graph overlay that splits nodes
 * into multiple nodes to satisfy a maximum length
 */

#ifndef VG_CHOPPED_GRAPH_HPP_INCLUDED
#define VG_CHOPPED_GRAPH_HPP_INCLUDED

#include "handle.hpp"
#include "utility.hpp"

namespace vg {

using namespace std;

/**
 * A HandleGraph implementation that overlays some other handle graph and
 * splits its nodes to satisfy a maximum length
 */
template<uint64_t MaxNodeLength = 1024>
class ChoppedGraph : virtual public HandleGraph {
public:
    
    /// Initialize as the reverse version of another graph, optionally also
    /// complementing
    ChoppedGraph(const HandleGraph& graph);
    
    /// Default constructor -- not actually functional
    ChoppedGraph() = default;
    
    /// Default destructor
    virtual ~ChoppedGraph() = default;
    
    //////////////////////////
    /// HandleGraph interface
    //////////////////////////
    
    // Method to check if a node exists by ID
    bool has_node(nid_t node_id) const;
    
    /// Returns true if there is an edge that allows traversal from the left
    /// handle to the right handle. By default O(n) in the number of edges
    /// on left, but can be overridden with more efficient implementations.
    bool has_edge(const handle_t& left, const handle_t& right) const;
    
    /// Look up the handle for the node with the given ID in the given orientation
    handle_t get_handle(const nid_t& node_id, bool is_reverse = false) const;
    
    /// Get the ID from a handle
    nid_t get_id(const handle_t& handle) const;
    
    /// Get the orientation of a handle
    bool get_is_reverse(const handle_t& handle) const;
    
    /// Invert the orientation of a handle (potentially without getting its ID)
    handle_t flip(const handle_t& handle) const;
    
    /// Get the length of a node
    size_t get_length(const handle_t& handle) const;
    
    /// Get the sequence of a node, presented in the handle's local forward
    /// orientation.
    string get_sequence(const handle_t& handle) const;
    
    /// Returns a substring of a handle's sequence, in the orientation of the
    /// handle. If the indicated substring would extend beyond the end of the
    /// handle's sequence, the return value is truncated to the sequence's end.
    string get_subsequence(const handle_t& handle, size_t index, size_t size) const;
    
    /// Returns one base of a handle's sequence, in the orientation of the
    /// handle.
    char get_base(const handle_t& handle, size_t index) const;
    
    /// Get the number of edges on the right (go_left = false) or left (go_left
    /// = true) side of the given handle.
    size_t get_degree(const handle_t& handle, bool go_left) const;
    
    /// Return the number of nodes in the graph
    size_t get_node_count() const;
    
    /// Return the total length of all nodes in the graph, in bp. If not
    /// overridden, loops over all nodes in linear time.
    size_t get_total_length() const;
    
    /// Return the total number of edges in the graph. If not overridden,
    /// counts them all in linear time.
    size_t get_edge_count() const;
    
    /// Return the smallest ID in the graph, or some smaller number if the
    /// smallest ID is unavailable. Return value is unspecified if the graph is empty.
    nid_t min_node_id() const;
    
    /// Return the largest ID in the graph, or some larger number if the
    /// largest ID is unavailable. Return value is unspecified if the graph is empty.
    nid_t max_node_id() const;
    
protected:
    
    /// Loop over all the handles to next/previous (right/left) nodes. Passes
    /// them to a callback which returns false to stop iterating and true to
    /// continue. Returns true if we finished and false if we stopped early.
    bool follow_edges_impl(const handle_t& handle, bool go_left,
                           const function<bool(const handle_t&)>& iteratee) const;
    
    /// Loop over all the nodes in the graph in their local forward
    /// orientations, in their internal stored order. Stop if the iteratee
    /// returns false. Can be told to run in parallel, in which case stopping
    /// after a false return value is on a best-effort basis and iteration
    /// order is not defined.
    bool for_each_handle_impl(const function<bool(const handle_t&)>& iteratee,
                              bool parallel = false) const;
    
    /// The index among an underlying handle's chop segments this one is
    size_t chop_index(const handle_t& handle) const;
    
    /// The number of segments a node has been chopped into
    size_t chopped_size(const handle_t& underlying) const;
    
    /// The n-th chopped segment of a node from the underlying graph
    handle_t chopped_handle(const handle_t& underlying_handle, size_t chop_idx) const;
    
    /**
     * Returns the handle in the underlying graph that corresponds to a handle in the
     * overlay
     */
    handle_t get_underlying_handle(const handle_t& handle) const;
    
    /// The underlying graph we're making splitting
    const HandleGraph* graph = nullptr;
    
    /// TODO: this might be more efficient with a map in the other
    /// direction instead (although it would take more memory)
    /// this would a raw index encoding for handles that makes chopped index queries
    /// need a galloping search instead of get handle queries
    vector<nid_t> chopped_id_begin;
    
    /// The number of bits used to represent the index among chopped segments
    uint8_t chopped_index_bit_count = 0;
    /// A mask to pull the chop index out of a packed handle
    uint64_t chopped_index_mask = 0;
    /// The average ratio between a node ID and its underlying ID
    double find_multiplier = 0.0;
};

/**
 * The same as ChoppedGraph except with the PathHandleGraph interface added
 */
template<uint64_t MaxNodeLength = 1024>
class ChoppedPathGraph : public ChoppedGraph<MaxNodeLength>, virtual public PathHandleGraph {
public:
        
    ChoppedPathGraph(const PathHandleGraph& graph);
    ChoppedPathGraph() = default;
    virtual ~ChoppedPathGraph() = default;
    
    /// Returns the number of paths stored in the graph
    size_t get_path_count() const;
    
    /// Determine if a path name exists and is legal to get a path handle for.
    bool has_path(const string& path_name) const;
    
    /// Look up the path handle for the given path name.
    /// The path with that name must exist.
    path_handle_t get_path_handle(const string& path_name) const;
    
    /// Look up the name of a path from a handle to it
    string get_path_name(const path_handle_t& path_handle) const;
    
    /// Look up whether a path is circular
    bool get_is_circular(const path_handle_t& path_handle) const;
    
    /// Returns the number of node steps in the path
    size_t get_step_count(const path_handle_t& path_handle) const;
    
    /// Returns the number of node steps on a handle
    size_t get_step_count(const handle_t& handle) const;
    
    /// Returns true if the given path is empty, and false otherwise
    bool is_empty(const path_handle_t& path_handle) const;
    
    /// Get a node handle (node ID and orientation) from a handle to an step on a path
    handle_t get_handle_of_step(const step_handle_t& step_handle) const;
    
    /// Returns a handle to the path that an step is on
    path_handle_t get_path_handle_of_step(const step_handle_t& step_handle) const;
    
    /// Get a handle to the first step, which will be an arbitrary step in a circular path
    /// that we consider "first" based on our construction of the path. If the path is empty,
    /// then the implementation must return the same value as path_end().
    step_handle_t path_begin(const path_handle_t& path_handle) const;
    
    /// Get a handle to a fictitious position past the end of a path. This position is
    /// returned by get_next_step for the final step in a path in a non-circular path.
    /// Note: get_next_step will *NEVER* return this value for a circular path.
    step_handle_t path_end(const path_handle_t& path_handle) const;
    
    /// Get a handle to the last step, which will be an arbitrary step in a circular path that
    /// we consider "last" based on our construction of the path. If the path is empty
    /// then the implementation must return the same value as path_front_end().
    step_handle_t path_back(const path_handle_t& path_handle) const;
    
    /// Get a handle to a fictitious position before the beginning of a path. This position is
    /// return by get_previous_step for the first step in a path in a non-circular path.
    /// Note: get_previous_step will *NEVER* return this value for a circular path.
    step_handle_t path_front_end(const path_handle_t& path_handle) const;
    
    /// Returns true if the step is not the last step in a non-circular path.
    bool has_next_step(const step_handle_t& step_handle) const;
    
    /// Returns true if the step is not the first step in a non-circular path.
    bool has_previous_step(const step_handle_t& step_handle) const;
    
    /// Returns a handle to the next step on the path. If the given step is the final step
    /// of a non-circular path, this method has undefined behavior. In a circular path,
    /// the "last" step will loop around to the "first" step.
    step_handle_t get_next_step(const step_handle_t& step_handle) const;
    
    /// Returns a handle to the previous step on the path. If the given step is the first
    /// step of a non-circular path, this method has undefined behavior. In a circular path,
    /// it will loop around from the "first" step (i.e. the one returned by path_begin) to
    /// the "last" step.
    step_handle_t get_previous_step(const step_handle_t& step_handle) const;
    
protected:
    
    /// The step from the underlying graph that this one is derived from
    step_handle_t get_underlying_step_handle(const step_handle_t& step_handle) const;
    
    /// The chopped step handle that corresponds to the n-th chopped segment of the node
    step_handle_t chopped_step_handle(const step_handle_t& underlying_step_handle, size_t chop_idx) const;
    
    /// Which chopped segment of the node does this correspond to
    size_t chop_index(const step_handle_t& step_handle) const;
    
    /// Execute a function on each path in the graph. If it returns false, stop
    /// iteration. Returns true if we finished and false if we stopped early.
    bool for_each_path_handle_impl(const std::function<bool(const path_handle_t&)>& iteratee) const;
    
    /// Execute a function on each step of a handle in any path. If it
    /// returns false, stop iteration. Returns true if we finished and false if
    /// we stopped early.
    bool for_each_step_on_handle_impl(const handle_t& handle,
                                      const std::function<bool(const step_handle_t&)>& iteratee) const;
    
    /// For simplicity, we store the same pointer again with the more general type
    const PathHandleGraph* path_graph;
    
    /// This is nontrivial to compute on the fly, so we'll memoize these values
    spp::sparse_hash_map<path_handle_t, size_t> chopped_steps_count;
    
    
    /// Make some inherited members accesible without scope resolutions
    using ChoppedGraph<MaxNodeLength>::get_underlying_handle;
    using ChoppedGraph<MaxNodeLength>::chop_index;
    using ChoppedGraph<MaxNodeLength>::chopped_size;
    using ChoppedGraph<MaxNodeLength>::chopped_index_bit_count;
    using ChoppedGraph<MaxNodeLength>::chopped_index_mask;
    using ChoppedGraph<MaxNodeLength>::chopped_handle;
};

/**
 * The same as ChoppedPathGraph except with the PathPositionHandleGraph interface added
 */
template<uint64_t MaxNodeLength = 1024>
class ChoppedPathPositionGraph : public ChoppedPathGraph<MaxNodeLength>, virtual public PathPositionHandleGraph {
public:
    
    ChoppedPathPositionGraph(const PathPositionHandleGraph& graph);
    ChoppedPathPositionGraph() = default;
    virtual ~ChoppedPathPositionGraph() = default;
    
    /// Returns the length of a path measured in bases of sequence.
    size_t get_path_length(const path_handle_t& path_handle) const;
    
    /// Returns the position along the path of the beginning of this step measured in
    /// bases of sequence. In a circular path, positions start at the step returned by
    /// path_begin().
    size_t get_position_of_step(const step_handle_t& step) const;
    
    /// Returns the step at this position, measured in bases of sequence starting at
    /// the step returned by path_begin(). If the position is past the end of the
    /// path, returns path_end().
    step_handle_t get_step_at_position(const path_handle_t& path, const size_t& position) const;
    
protected:
    
    const PathPositionHandleGraph* path_position_graph;
    
    using ChoppedPathGraph<MaxNodeLength>::chopped_step_handle;
    using ChoppedPathGraph<MaxNodeLength>::chop_index;
    using ChoppedPathGraph<MaxNodeLength>::get_underlying_step_handle;
};








/**
 *
 * Template implementations
 *
 */

template<uint64_t MaxNodeLength>
ChoppedGraph<MaxNodeLength>::ChoppedGraph(const HandleGraph& graph) : graph(&graph) {
    
    // note: assumes dense IDs starting at 1 for efficiency
    // make a vector where every ID can also index before and after by 1
    chopped_id_begin.resize(graph.max_node_id() + 2, 0);
    size_t max_chopped_count = 0;
    graph.for_each_handle([&](const handle_t& handle) {
        // assign it an interval of IDs for the chopped nodes
        size_t chopped_count = max<size_t>(1, (graph.get_length(handle) + MaxNodeLength - 1) / MaxNodeLength);
        // temporarily assign only the count in the vector
        chopped_id_begin[graph.get_id(handle) + 1] = chopped_count;
        max_chopped_count = max(chopped_count, max_chopped_count);
    });
    
    // convert it into a prefix sum
    chopped_id_begin[0] = 1; // make IDs start at 1
    for (size_t i = 1; i < chopped_id_begin.size(); ++i) {
        chopped_id_begin[i] += chopped_id_begin[i - 1];
    }
    
    // inefficient way to find most significant bit, but whatever
    for (uint8_t i = 1; i <= 64; ++i) {
        if (max_chopped_count & (1 << (i - 1))) {
            chopped_index_bit_count = i;
        }
    }
    chopped_index_mask = numeric_limits<uint64_t>::max() >> (64 - chopped_index_bit_count);
    
    // memoize this constant that we frequently use for find
    // note: have to use double because otherwise risk of overflow
    // when > 2^32 nodes
    find_multiplier = double(chopped_id_begin.size()) / double(max_node_id());
}

template<uint64_t MaxNodeLength>
size_t ChoppedGraph<MaxNodeLength>::chop_index(const handle_t& handle) const {
    return handlegraph::as_integer(handle) & chopped_index_mask;
}

template<uint64_t MaxNodeLength>
size_t ChoppedGraph<MaxNodeLength>::chopped_size(const handle_t& underlying) const {
    // this is also implicitly available from the node length...
    nid_t node_id = graph->get_id(underlying);
    return chopped_id_begin[node_id + 1] - chopped_id_begin[node_id];
}

template<uint64_t MaxNodeLength>
handle_t ChoppedGraph<MaxNodeLength>::chopped_handle(const handle_t& underlying_handle, size_t chop_idx) const {
    return handlegraph::as_handle((handlegraph::as_integer(underlying_handle) << chopped_index_bit_count) | chop_idx);
}

template<uint64_t MaxNodeLength>
bool ChoppedGraph<MaxNodeLength>::has_node(nid_t node_id) const {
    // IDs form a contiguous range, so being less than the max is a sufficient condition
    return node_id <= max_node_id();
}

template<uint64_t MaxNodeLength>
bool ChoppedGraph<MaxNodeLength>::has_edge(const handle_t& left, const handle_t& right) const {
    
    handle_t under_left = get_underlying_handle(left);
    handle_t under_right = get_underlying_handle(right);
    
    if (((graph->get_is_reverse(under_left) && chop_index(left) == 0) ||
         (!graph->get_is_reverse(under_left) && chop_index(left) + 1 == chopped_size(under_left))) &&
        ((graph->get_is_reverse(under_right) && chop_index(right) + 1 == chopped_size(under_right)) ||
         (!graph->get_is_reverse(under_right) && chop_index(right) == 0))) {
        
        // at the boundary of unchopped nodes
        return graph->has_edge(under_left, under_right);
    }
    else {
        // not at the boundary of unchopped nodes, but possibly at a boundary internal to an unchopped node
        return (under_left == under_right &&
                ((graph->get_is_reverse(under_left) && chop_index(right) + 1 == chop_index(left)) ||
                 (!graph->get_is_reverse(under_left) && chop_index(right) == chop_index(left) + 1)));
    }
}

template<uint64_t MaxNodeLength>
handle_t ChoppedGraph<MaxNodeLength>::get_handle(const nid_t& node_id, bool is_reverse) const {
    // start searching around where we'd expect to find the node ID
    // if nodes were evenly distributed through the vector
    // TODO: do i actually need the min here?
    nid_t search_min = min<nid_t>(node_id * find_multiplier,
                                   chopped_id_begin.size() - 1);
    nid_t search_max = search_min + 1;
    
    // expand search windows out with galloping search in either direction
    // as necessary
    while (search_min != 0 &&
           chopped_id_begin[search_min] > node_id) {
        nid_t stride = 2 * search_max - search_min;
        search_max = search_min;
        search_min = search_min > stride ? search_min - stride : 0;
    }
    while (search_max != chopped_id_begin.size() &&
           chopped_id_begin[search_max] <= node_id) {
        nid_t stride = 2 * search_max - search_min;
        search_min = search_max;
        search_max = min<nid_t>(search_max + stride, chopped_id_begin.size());
    }
    
    // binary search to find the highest index that is at most
    // this node ID
    while (search_min + 1 < search_max) {
        nid_t middle = (search_min + search_max) / 2;
        if (chopped_id_begin[middle] > node_id) {
            search_max = middle;
        }
        else {
            search_min = middle;
        }
    }
    
    size_t chop_idx = node_id - chopped_id_begin[search_min];
    return chopped_handle(graph->get_handle(search_min, is_reverse), chop_idx);
}

template<uint64_t MaxNodeLength>
nid_t ChoppedGraph<MaxNodeLength>::get_id(const handle_t& handle) const {
    return chopped_id_begin[graph->get_id(get_underlying_handle(handle))] + chop_index(handle);
}

template<uint64_t MaxNodeLength>
bool ChoppedGraph<MaxNodeLength>::get_is_reverse(const handle_t& handle) const {
    return graph->get_is_reverse(get_underlying_handle(handle));
}

template<uint64_t MaxNodeLength>
handle_t ChoppedGraph<MaxNodeLength>::flip(const handle_t& handle) const {
    return chopped_handle(graph->flip(get_underlying_handle(handle)), chop_index(handle));
}

template<uint64_t MaxNodeLength>
size_t ChoppedGraph<MaxNodeLength>::get_length(const handle_t& handle) const {
    handle_t underlying = get_underlying_handle(handle);
    if (chop_index(handle) + 1 == chopped_size(underlying)) {
        size_t length = graph->get_length(underlying);
        return length ? ((length - 1) % MaxNodeLength) + 1 : 0;
    }
    else {
        return MaxNodeLength;
    }
}

template<uint64_t MaxNodeLength>
string ChoppedGraph<MaxNodeLength>::get_sequence(const handle_t& handle) const {
    handle_t underlying = get_underlying_handle(handle);
    string seq = graph->get_subsequence(graph->forward(underlying),
                                        chop_index(handle) * MaxNodeLength, MaxNodeLength);
    if (graph->get_is_reverse(underlying)) {
        reverse_complement_in_place(seq);
    }
    return seq;
}

template<uint64_t MaxNodeLength>
string ChoppedGraph<MaxNodeLength>::get_subsequence(const handle_t& handle,
                                                    size_t index, size_t size) const {
    handle_t underlying = get_underlying_handle(handle);
    size_t under_idx;
    if (graph->get_is_reverse(underlying)) {
        under_idx = max<int64_t>(0, graph->get_length(underlying) - (chop_index(handle) + 1) * MaxNodeLength) + index;
    }
    else {
        under_idx = chop_index(handle) * MaxNodeLength + index;
    }
    return graph->get_subsequence(underlying, under_idx, min<size_t>(MaxNodeLength - index, size));
}

template<uint64_t MaxNodeLength>
char ChoppedGraph<MaxNodeLength>::get_base(const handle_t& handle, size_t index) const {
    handle_t underlying = get_underlying_handle(handle);
    // TODO: repetetive with above
    size_t under_idx;
    if (graph->get_is_reverse(underlying)) {
        under_idx = max<int64_t>(0, graph->get_length(underlying) - (chop_index(handle) + 1) * MaxNodeLength) + index;
    }
    else {
        under_idx = chop_index(handle) * MaxNodeLength + index;
    }
    return graph->get_base(underlying, under_idx);
}

template<uint64_t MaxNodeLength>
size_t ChoppedGraph<MaxNodeLength>::get_degree(const handle_t& handle, bool go_left) const {
    handle_t underlying = get_underlying_handle(handle);
    bool leftward = (graph->get_is_reverse(underlying) != go_left);
    size_t chop_idx = chop_index(handle);
    if ((leftward && chop_idx == 0) || (!leftward && chop_idx + 1 == chopped_size(underlying))) {
        // across edges in the underlying graph
        return graph->get_degree(underlying, go_left);
    }
    else {
        // across a chop edge
        return 1;
    }
}

template<uint64_t MaxNodeLength>
bool ChoppedGraph<MaxNodeLength>::follow_edges_impl(const handle_t& handle, bool go_left,
                                     const function<bool(const handle_t&)>& iteratee) const {
    
    size_t chop_idx =  chop_index(handle);
    handle_t underlying = get_underlying_handle(handle);
    bool traverse_leftward = go_left != graph->get_is_reverse(underlying);
    if (traverse_leftward && chop_idx != 0) {
        // move left internally to the chopped node
        return iteratee(chopped_handle(underlying, chop_idx - 1));
    }
    else if (!traverse_leftward && chop_idx + 1 < chopped_size(underlying)) {
        // move right internally to the chopped node
        return iteratee(chopped_handle(underlying, chop_idx + 1));
    }
    else {
        // follow edges in the underlying graph
        return graph->follow_edges(underlying, go_left,
                                   [&](const handle_t& next) {
            return iteratee(chopped_handle(next,
                                           go_left != graph->get_is_reverse(next) ?
                                           chopped_size(next) - 1 : 0));
        });
    }
}

template<uint64_t MaxNodeLength>
bool ChoppedGraph<MaxNodeLength>::for_each_handle_impl(const function<bool(const handle_t&)>& iteratee,
                                        bool parallel) const {
    
    return graph->for_each_handle([&](const handle_t& handle) {
        bool keep_going = true;
        for (size_t i = 0, n = chopped_size(handle); i < n && keep_going; ++i) {
            keep_going = iteratee(chopped_handle(handle, i));
        }
        return keep_going;
    }, parallel);
}

template<uint64_t MaxNodeLength>
size_t ChoppedGraph<MaxNodeLength>::get_node_count() const {
    // nodes start at 1 and are allocated contiguously, so these are the same
    return max_node_id();
}

template<uint64_t MaxNodeLength>
size_t ChoppedGraph<MaxNodeLength>::get_total_length() const {
    // contains the same sequence
    return graph->get_total_length();
}

template<uint64_t MaxNodeLength>
size_t ChoppedGraph<MaxNodeLength>::get_edge_count() const {
    // one edge is added for each new node that's added
    return graph->get_edge_count() + (get_node_count() - graph->get_node_count());
}

template<uint64_t MaxNodeLength>
nid_t ChoppedGraph<MaxNodeLength>::min_node_id() const {
    return 1;
}

template<uint64_t MaxNodeLength>
nid_t ChoppedGraph<MaxNodeLength>::max_node_id() const {
    return chopped_id_begin.back() - 1;
}

template<uint64_t MaxNodeLength>
handle_t ChoppedGraph<MaxNodeLength>::get_underlying_handle(const handle_t& handle) const {
    return handlegraph::as_handle(handlegraph::as_integer(handle) >> chopped_index_bit_count);
}

template<uint64_t MaxNodeLength>
ChoppedPathGraph<MaxNodeLength>::ChoppedPathGraph(const PathHandleGraph& graph)
    : ChoppedGraph<MaxNodeLength>(graph), path_graph(&graph)
{
    // compute the step count for the chopped paths
    chopped_steps_count.reserve(graph.get_path_count());
    graph.for_each_path_handle([&](const path_handle_t& path_handle){
        size_t count = 0;
        for (handle_t handle : graph.scan_path(path_handle)) {
            count += chopped_size(handle);
        }
        chopped_steps_count[path_handle] = count;
    });
}


template<uint64_t MaxNodeLength>
step_handle_t ChoppedPathGraph<MaxNodeLength>::get_underlying_step_handle(const step_handle_t& step_handle) const {
    step_handle_t underlying;
    handlegraph::as_integers(underlying)[0] = handlegraph::as_integers(step_handle)[0] >> chopped_index_bit_count;
    handlegraph::as_integers(underlying)[1] = handlegraph::as_integers(step_handle)[1];
    return underlying;
}

template<uint64_t MaxNodeLength>
step_handle_t ChoppedPathGraph<MaxNodeLength>::chopped_step_handle(const step_handle_t& underlying_step_handle,
                                                                   size_t chop_idx) const {
    step_handle_t step;
    handlegraph::as_integers(step)[0] = handlegraph::as_integers(underlying_step_handle)[0] << chopped_index_bit_count | chop_idx;
    handlegraph::as_integers(step)[1] = handlegraph::as_integers(underlying_step_handle)[1];
    return step;
}

template<uint64_t MaxNodeLength>
size_t ChoppedPathGraph<MaxNodeLength>::chop_index(const step_handle_t& step_handle) const {
    return handlegraph::as_integers(step_handle)[0] & chopped_index_mask;
}

template<uint64_t MaxNodeLength>
size_t ChoppedPathGraph<MaxNodeLength>::get_path_count() const {
    return path_graph->get_path_count();
}

template<uint64_t MaxNodeLength>
bool ChoppedPathGraph<MaxNodeLength>::has_path(const string& path_name) const {
    return path_graph->has_path(path_name);
}

template<uint64_t MaxNodeLength>
path_handle_t ChoppedPathGraph<MaxNodeLength>::get_path_handle(const string& path_name) const {
    return path_graph->get_path_handle(path_name);
}

template<uint64_t MaxNodeLength>
string ChoppedPathGraph<MaxNodeLength>::get_path_name(const path_handle_t& path_handle) const {
    return path_graph->get_path_name(path_handle);
}

template<uint64_t MaxNodeLength>
bool ChoppedPathGraph<MaxNodeLength>::get_is_circular(const path_handle_t& path_handle) const {
    return path_graph->get_is_circular(path_handle);
}

template<uint64_t MaxNodeLength>
size_t ChoppedPathGraph<MaxNodeLength>::get_step_count(const path_handle_t& path_handle) const {
    return chopped_steps_count.at(path_handle);
}

template<uint64_t MaxNodeLength>
size_t ChoppedPathGraph<MaxNodeLength>::get_step_count(const handle_t& handle) const {
    return path_graph->get_step_count(ChoppedGraph<MaxNodeLength>::get_underlying_handle(handle));
}

template<uint64_t MaxNodeLength>
bool ChoppedPathGraph<MaxNodeLength>::is_empty(const path_handle_t& path_handle) const {
    return path_graph->is_empty(path_handle);
}

template<uint64_t MaxNodeLength>
handle_t ChoppedPathGraph<MaxNodeLength>::get_handle_of_step(const step_handle_t& step_handle) const {
    step_handle_t underlying = get_underlying_step_handle(step_handle);
    return this->chopped_handle(path_graph->get_handle_of_step(underlying), chop_index(step_handle));
}

template<uint64_t MaxNodeLength>
path_handle_t ChoppedPathGraph<MaxNodeLength>::get_path_handle_of_step(const step_handle_t& step_handle) const {
    return path_graph->get_path_handle_of_step(get_underlying_step_handle(step_handle));
}

template<uint64_t MaxNodeLength>
step_handle_t ChoppedPathGraph<MaxNodeLength>::path_begin(const path_handle_t& path_handle) const {
    step_handle_t underlying_begin = path_graph->path_begin(path_handle);
    if (underlying_begin == path_graph->path_end(path_handle)) {
        return underlying_begin;
    }
    else {
        handle_t handle = path_graph->get_handle_of_step(underlying_begin);
        if (path_graph->get_is_reverse(handle)) {
            return chopped_step_handle(underlying_begin, chopped_size(handle) - 1);
        }
        else {
            return chopped_step_handle(underlying_begin, 0);
        }
    }
}

template<uint64_t MaxNodeLength>
step_handle_t ChoppedPathGraph<MaxNodeLength>::path_end(const path_handle_t& path_handle) const {
    return path_graph->path_end(path_handle);
}

template<uint64_t MaxNodeLength>
step_handle_t ChoppedPathGraph<MaxNodeLength>::path_back(const path_handle_t& path_handle) const {
    step_handle_t underlying_back = path_graph->path_back(path_handle);
    if (underlying_back == path_graph->path_front_end(path_handle)) {
        return underlying_back;
    }
    else {
        handle_t handle = path_graph->get_handle_of_step(underlying_back);
        if (path_graph->get_is_reverse(handle)) {
            return chopped_step_handle(underlying_back, 0);
        }
        else {
            return chopped_step_handle(underlying_back, chopped_size(handle) - 1);
        }
    }
}

template<uint64_t MaxNodeLength>
step_handle_t ChoppedPathGraph<MaxNodeLength>::path_front_end(const path_handle_t& path_handle) const {
    return path_graph->path_front_end(path_handle);
}

template<uint64_t MaxNodeLength>
bool ChoppedPathGraph<MaxNodeLength>::has_next_step(const step_handle_t& step_handle) const {
    size_t chop_idx = chop_index(step_handle);
    step_handle_t underlying = get_underlying_step_handle(step_handle);
    bool leftward = path_graph->get_is_reverse(path_graph->get_handle_of_step(underlying));
    if ((leftward && chop_idx > 0) ||
        (!leftward && chop_idx + 1 < chopped_size(path_graph->get_handle_of_step(underlying)))) {
        return true;
    }
    else {
        return path_graph->has_next_step(underlying);
    }
}

template<uint64_t MaxNodeLength>
bool ChoppedPathGraph<MaxNodeLength>::has_previous_step(const step_handle_t& step_handle) const {
    size_t chop_idx = chop_index(step_handle);
    step_handle_t underlying = get_underlying_step_handle(step_handle);
    bool leftward = path_graph->get_is_reverse(path_graph->get_handle_of_step(underlying));
    if ((leftward && chop_idx + 1 < chopped_size(path_graph->get_handle_of_step(underlying))) ||
        (!leftward && chop_idx > 0)) {
        return true;
    }
    else {
        return path_graph->has_previous_step(underlying);
    }
}

template<uint64_t MaxNodeLength>
step_handle_t ChoppedPathGraph<MaxNodeLength>::get_next_step(const step_handle_t& step_handle) const {
    size_t chop_idx = chop_index(step_handle);
    step_handle_t underlying = get_underlying_step_handle(step_handle);
    bool leftward = path_graph->get_is_reverse(path_graph->get_handle_of_step(underlying));
    if ((leftward && chop_idx == 0) ||
        (!leftward && chop_idx + 1 == chopped_size(path_graph->get_handle_of_step(underlying)))) {
        step_handle_t next_underlying = path_graph->get_next_step(underlying);
        if (next_underlying == path_graph->path_end(path_graph->get_path_handle_of_step(underlying))) {
            return next_underlying;
        }
        else {
            handle_t next_handle = path_graph->get_handle_of_step(next_underlying);
            return chopped_step_handle(next_underlying,
                                       path_graph->get_is_reverse(next_handle) ?
                                       chopped_size(next_handle) - 1 : 0);
        }
        
    }
    else {
        return chopped_step_handle(underlying, leftward ? chop_idx - 1 : chop_idx + 1);
    }
}

template<uint64_t MaxNodeLength>
step_handle_t ChoppedPathGraph<MaxNodeLength>::get_previous_step(const step_handle_t& step_handle) const {
    
    size_t chop_idx = chop_index(step_handle);
    step_handle_t underlying = get_underlying_step_handle(step_handle);
    bool leftward = path_graph->get_is_reverse(path_graph->get_handle_of_step(underlying));
    
    if ((leftward && chop_idx + 1 == chopped_size(path_graph->get_handle_of_step(underlying))) ||
        (!leftward && chop_idx == 0)) {
        
        step_handle_t prev_underlying = path_graph->get_previous_step(underlying);
        if (prev_underlying == path_graph->path_front_end(path_graph->get_path_handle_of_step(underlying))) {
            return prev_underlying;
        }
        else {
            handle_t prev_handle = path_graph->get_handle_of_step(prev_underlying);
            return chopped_step_handle(prev_underlying,
                                       path_graph->get_is_reverse(prev_handle) ?
                                       0 : chopped_size(prev_handle) - 1);
        }
        
    }
    else {
        return chopped_step_handle(underlying, leftward ? chop_idx + 1 : chop_idx - 1);
    }
}

template<uint64_t MaxNodeLength>
bool ChoppedPathGraph<MaxNodeLength>::for_each_path_handle_impl(const std::function<bool(const path_handle_t&)>& iteratee) const {
    return path_graph->for_each_path_handle(iteratee);
}

template<uint64_t MaxNodeLength>
bool ChoppedPathGraph<MaxNodeLength>::for_each_step_on_handle_impl(const handle_t& handle,
                                  const std::function<bool(const step_handle_t&)>& iteratee) const {
    
    size_t chop_idx = chop_index(handle);
    return path_graph->for_each_step_on_handle(get_underlying_handle(handle), [&](const step_handle_t& step) {
        return iteratee(chopped_step_handle(step, chop_idx));
    });
}

template<uint64_t MaxNodeLength>
ChoppedPathPositionGraph<MaxNodeLength>::ChoppedPathPositionGraph(const PathPositionHandleGraph& graph)
    : ChoppedPathGraph<MaxNodeLength>(graph), path_position_graph(&graph)
{
    
}

template<uint64_t MaxNodeLength>
size_t ChoppedPathPositionGraph<MaxNodeLength>::get_path_length(const path_handle_t& path_handle) const {
    // chopping doesn't change path length
    return path_position_graph->get_path_length(path_handle);
}

template<uint64_t MaxNodeLength>
size_t ChoppedPathPositionGraph<MaxNodeLength>::get_position_of_step(const step_handle_t& step) const {
    step_handle_t underlying = get_underlying_step_handle(step);
    handle_t underlying_handle = path_position_graph->get_handle_of_step(underlying);
    if (!path_position_graph->get_is_reverse(underlying_handle)) {
        return path_position_graph->get_position_of_step(underlying) + chop_index(step) * MaxNodeLength;
    }
    else {
        return (path_position_graph->get_position_of_step(underlying)
                + max<int64_t>(0, path_position_graph->get_length(underlying_handle) - (chop_index(step) + 1) * MaxNodeLength));
    }
}

template<uint64_t MaxNodeLength>
step_handle_t ChoppedPathPositionGraph<MaxNodeLength>::get_step_at_position(const path_handle_t& path,
                                                                            const size_t& position) const {
    step_handle_t underlying = path_position_graph->get_step_at_position(path, position);
    size_t remaining = position - path_position_graph->get_position_of_step(underlying);
    handle_t underlying_handle = path_position_graph->get_handle_of_step(underlying);
    if (!path_position_graph->get_is_reverse(underlying_handle)) {
        return chopped_step_handle(underlying, remaining / MaxNodeLength);
    }
    else {
        return chopped_step_handle(underlying, (path_position_graph->get_length(underlying_handle) - remaining - 1) / MaxNodeLength);
    }
}


}

#endif
