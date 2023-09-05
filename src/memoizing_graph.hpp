#ifndef VG_MEMOIZING_GRAPH_HPP_INCLUDED
#define VG_MEMOIZING_GRAPH_HPP_INCLUDED

/** \file
 * memoizing_graph.hpp: defines a handle graph implementation memoizes the results
 * of certain handle operations
 */

#include "handle.hpp"

#include <unordered_map>

namespace vg {

using namespace std;

    /**
     * A PathPositionHandleGraph implementation that memoizes the results of get_handle
     * and steps_of_handle.
     */
    class MemoizingGraph : public PathPositionHandleGraph {
    public:
        
        /// Initialize with a pointer to graph we want to memoize operations for
        MemoizingGraph(const PathPositionHandleGraph* graph);
        
        /// Default constructor -- not actually functional
        MemoizingGraph() = default;
        
        /// Default destructor
        ~MemoizingGraph() = default;
        
        //////////////////////////
        /// HandleGraph interface
        //////////////////////////
        
        /// Method to check if a node exists by ID
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
        
        /// Return the number of nodes in the graph.
        virtual size_t get_node_count() const;
        
        /// Return the smallest ID in the graph, or some smaller number if the
        /// smallest ID is unavailable. Return value is unspecified if the graph is empty.
        virtual id_t min_node_id() const;
        
        /// Return the largest ID in the graph, or some larger number if the
        /// largest ID is unavailable. Return value is unspecified if the graph is empty.
        virtual id_t max_node_id() const;
        
        ////////////////////////////////////////////
        // Path handle graph interface
        ////////////////////////////////////////////
        
        /// Returns the number of paths stored in the graph
        virtual size_t get_path_count() const;
        
        /// Determine if a path name exists and is legal to get a path handle for.
        virtual bool has_path(const std::string& path_name) const;
        
        /// Look up the path handle for the given path name.
        /// The path with that name must exist.
        virtual path_handle_t get_path_handle(const std::string& path_name) const;
        
        /// Look up the name of a path from a handle to it
        virtual std::string get_path_name(const path_handle_t& path_handle) const;
        
        /// Look up whether a path is circular
        virtual bool get_is_circular(const path_handle_t& path_handle) const;
        
        /// Returns the number of node steps in the path
        virtual size_t get_step_count(const path_handle_t& path_handle) const;
        
        /// Get a node handle (node ID and orientation) from a handle to an step on a path
        virtual handle_t get_handle_of_step(const step_handle_t& step_handle) const;
        
        /// Returns a handle to the path that an step is on
        virtual path_handle_t get_path_handle_of_step(const step_handle_t& step_handle) const;
        
        /// Get a handle to the first step, which will be an arbitrary step in a circular path
        /// that we consider "first" based on our construction of the path. If the path is empty,
        /// then the implementation must return the same value as path_end().
        virtual step_handle_t path_begin(const path_handle_t& path_handle) const;
        
        /// Get a handle to a fictitious position past the end of a path. This position is
        /// returned by get_next_step for the final step in a path in a non-circular path.
        /// Note: get_next_step will *NEVER* return this value for a circular path.
        virtual step_handle_t path_end(const path_handle_t& path_handle) const;
        
        /// Get a handle to the last step, which will be an arbitrary step in a circular path that
        /// we consider "last" based on our construction of the path. If the path is empty
        /// then the implementation must return the same value as path_front_end().
        virtual step_handle_t path_back(const path_handle_t& path_handle) const;
        
        /// Get a handle to a fictitious position before the beginning of a path. This position is
        /// return by get_previous_step for the first step in a path in a non-circular path.
        /// Note: get_previous_step will *NEVER* return this value for a circular path.
        virtual step_handle_t path_front_end(const path_handle_t& path_handle) const;
        
        /// Returns true if the step is not the last step in a non-circular path.
        virtual bool has_next_step(const step_handle_t& step_handle) const;
        
        /// Returns true if the step is not the first step in a non-circular path.
        virtual bool has_previous_step(const step_handle_t& step_handle) const;
        
        /// Returns a handle to the next step on the path. If the given step is the final step
        /// of a non-circular path, this method has undefined behavior. In a circular path,
        /// the "last" step will loop around to the "first" step.
        virtual step_handle_t get_next_step(const step_handle_t& step_handle) const;
        
        /// Returns a handle to the previous step on the path. If the given step is the first
        /// step of a non-circular path, this method has undefined behavior. In a circular path,
        /// it will loop around from the "first" step (i.e. the one returned by path_begin) to
        /// the "last" step.
        virtual step_handle_t get_previous_step(const step_handle_t& step_handle) const;
        
    protected:
        
        /// Execute a function on each path in the graph. If it returns false, stop
        /// iteration. Returns true if we finished and false if we stopped early.
        virtual bool for_each_path_handle_impl(const std::function<bool(const path_handle_t&)>& iteratee) const;
        
        /// Execute a function on each step of a handle in any path. If it
        /// returns false, stop iteration. Returns true if we finished and false if
        /// we stopped early.
        virtual bool for_each_step_on_handle_impl(const handle_t& handle,
                                                  const std::function<bool(const step_handle_t&)>& iteratee) const;
        
    public:
        
        /// Returns a vector of all steps of a node on paths. Optionally restricts to
        /// steps that match the handle in orientation.
        virtual std::vector<step_handle_t> steps_of_handle(const handle_t& handle,
                                                           bool match_orientation = false) const;
        
        /// Returns true if the given path is empty, and false otherwise
        virtual bool is_empty(const path_handle_t& path_handle) const;
        
        
        ////////////////////////////////////////////////////////////////////////////
        // Path position handle graph interface
        ////////////////////////////////////////////////////////////////////////////
        
        /// Returns the length of a path measured in bases of sequence.
        virtual size_t get_path_length(const path_handle_t& path_handle) const;
        
        /// Returns the position along the path of the beginning of this step measured in
        /// bases of sequence. In a circular path, positions start at the step returned by
        /// path_begin().
        virtual size_t get_position_of_step(const step_handle_t& step) const;
        
        /// Returns the step at this position, measured in bases of sequence starting at
        /// the step returned by path_begin(). If the position is past the end of the
        /// path, returns path_end().
        virtual step_handle_t get_step_at_position(const path_handle_t& path,
                                                   const size_t& position) const;
        
    protected:
        
        /// Execute an itteratee on each step and its path relative position and orientation
        /// on a handle in any path. Iteration will stop early if the iteratee returns false.
        /// This method returns false if iteration was stopped early, else true.
        virtual bool for_each_step_position_on_handle(const handle_t& handle,
                                                      const std::function<bool(const step_handle_t&, const bool&, const size_t&)>& iteratee) const;
        
    public:
        
        /// The largest number of calls to get_handle we will memoize
        size_t max_handle_memo_size = 500;
        
        /// The largest number of calls to steps_of_handle we will memoize
        size_t max_steps_of_handle_memo_size = 500;
        
    private:
        /// The graph we're memoizing operations for
        const PathPositionHandleGraph* graph = nullptr;
        
        /// Memo for get_handle
        unordered_map<id_t, handle_t> get_handle_memo;
        
        /// Memo for steps_of_handle
        unordered_map<handle_t, vector<step_handle_t>> steps_of_handle_memo;
    };
}

#endif
