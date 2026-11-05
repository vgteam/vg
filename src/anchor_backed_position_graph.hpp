#ifndef VG_ANCHOR_BACKED_POSITION_GRAPH_HPP_INCLUDED
#define VG_ANCHOR_BACKED_POSITION_GRAPH_HPP_INCLUDED

/** \file
 * A path position graph that knows positions on one path from a few anchors,
 * so a Surjector can surject onto a haplotype path without a position index
 * over it.
 */

#include "handle.hpp"

#include <handlegraph/path_position_handle_graph.hpp>

#include <unordered_map>
#include <utility>
#include <vector>

namespace vg {

using namespace std;

/// A run of steps on the target path whose first and last steps have known
/// positions.
struct PrecomputedAnchor {
    /// The first step of the run.
    step_handle_t step_begin;
    /// The last step of the run, which may be step_begin.
    step_handle_t step_end;
    /// The 0-based offset on the target path of the first base of step_begin.
    size_t path_offset_step_begin = 0;
    /// The 0-based offset on the target path of the first base of step_end.
    size_t path_offset_step_end = 0;
};

/// Counts of the position queries on this thread and of the steps walked to
/// answer them. Reset it before a surjection and read it afterwards.
struct AnchorGraphStats {
    size_t pos_of_step_calls = 0;  ///< calls to get_position_of_step
    size_t pos_of_step_walk  = 0;  ///< steps walked by get_position_of_step
    size_t step_at_pos_calls = 0;  ///< calls to get_step_at_position
    size_t step_at_pos_walk  = 0;  ///< steps walked by get_step_at_position
};
extern thread_local AnchorGraphStats g_anchor_graph_stats;

/**
 * A PathPositionHandleGraph over a PathHandleGraph that knows positions on one
 * target path only
 *
 * The positions come from anchors: steps whose positions the caller already
 * knows. A position query walks the path from the nearest step with a known
 * position and remembers every step it passes. Everything else goes to the
 * base graph.
 *
 * The graph is built for surjecting one read onto the target path. Only the
 * steps of the target path inside the anchors' runs are reported as steps on
 * their nodes, so a Surjector ignores the many other haplotype paths on each
 * node.
 *
 * @pre The anchors' runs cover every step of the target path that a surjected
 * alignment visits.
 *
 * @note Position queries change the cache of known positions, so one object
 * must not be queried from two threads at once.
 */
class AnchorBackedPositionGraph : public PathPositionHandleGraph {
public:

    /**
     * Make a graph that knows positions on target_path from the anchors.
     *
     * @param base The graph to forward to. It is not owned.
     * @param target_path_length The length of target_path in bases, which
     * get_path_length returns and which must be correct.
     * @throws std::invalid_argument if base is null.
     */
    AnchorBackedPositionGraph(const PathHandleGraph* base,
                              const std::vector<PrecomputedAnchor>& anchors,
                              path_handle_t target_path,
                              size_t target_path_length);

    ~AnchorBackedPositionGraph() = default;

    AnchorBackedPositionGraph() = delete;
    AnchorBackedPositionGraph(const AnchorBackedPositionGraph&) = delete;
    AnchorBackedPositionGraph& operator=(const AnchorBackedPositionGraph&) = delete;

    //////////////////////////
    // HandleGraph interface, forwarded to the base graph
    //////////////////////////

    virtual bool has_node(nid_t node_id) const override;
    virtual handle_t get_handle(const nid_t& node_id, bool is_reverse = false) const override;
    virtual nid_t get_id(const handle_t& handle) const override;
    virtual bool get_is_reverse(const handle_t& handle) const override;
    virtual handle_t flip(const handle_t& handle) const override;
    virtual size_t get_length(const handle_t& handle) const override;
    virtual std::string get_sequence(const handle_t& handle) const override;
    virtual size_t get_node_count() const override;
    virtual nid_t min_node_id() const override;
    virtual nid_t max_node_id() const override;

    //////////////////////////
    // PathHandleGraph interface, forwarded to the base graph
    //////////////////////////

    virtual size_t get_path_count() const override;
    virtual bool has_path(const std::string& path_name) const override;
    virtual path_handle_t get_path_handle(const std::string& path_name) const override;
    virtual std::string get_path_name(const path_handle_t& path_handle) const override;
    virtual bool get_is_circular(const path_handle_t& path_handle) const override;
    virtual size_t get_step_count(const path_handle_t& path_handle) const override;
    virtual handle_t get_handle_of_step(const step_handle_t& step_handle) const override;
    virtual path_handle_t get_path_handle_of_step(const step_handle_t& step_handle) const override;
    virtual step_handle_t path_begin(const path_handle_t& path_handle) const override;
    virtual step_handle_t path_end(const path_handle_t& path_handle) const override;
    virtual step_handle_t path_back(const path_handle_t& path_handle) const override;
    virtual step_handle_t path_front_end(const path_handle_t& path_handle) const override;
    virtual bool has_next_step(const step_handle_t& step_handle) const override;
    virtual bool has_previous_step(const step_handle_t& step_handle) const override;
    virtual step_handle_t get_next_step(const step_handle_t& step_handle) const override;
    virtual step_handle_t get_previous_step(const step_handle_t& step_handle) const override;

    //////////////////////////
    // PathPositionHandleGraph interface, for the target path only
    //////////////////////////

    /// Return the length of the target path given to the constructor.
    /// @throws std::logic_error for any other path.
    virtual size_t get_path_length(const path_handle_t& path_handle) const override;

    /**
     * Return the offset of the first base of a step on the target path.
     *
     * We walk back to the nearest step with a known position, so a query costs
     * the distance to it, not to the path start.
     *
     * @throws std::logic_error if the step is on another path.
     */
    virtual size_t get_position_of_step(const step_handle_t& step) const override;

    /**
     * Return the step of the target path that covers a position.
     *
     * We walk to it from the closer of the known positions on either side.
     *
     * @returns path_end for a position at or past the end of the path.
     * @throws std::logic_error for any other path.
     */
    virtual step_handle_t get_step_at_position(const path_handle_t& path,
                                               const size_t& position) const override;

    //////////////////////////
    // Inspection
    //////////////////////////

    /// Get the number of steps with known positions.
    size_t cached_position_count() const { return step_pos_.size(); }

    /// Get the path that this graph knows positions on.
    path_handle_t target_path() const { return target_path_; }

protected:

    // Iteration, forwarded to the base graph except for
    // for_each_step_on_handle_impl, which reports target path steps only.
    virtual bool follow_edges_impl(const handle_t& handle, bool go_left,
                                   const std::function<bool(const handle_t&)>& iteratee) const override;
    virtual bool for_each_handle_impl(const std::function<bool(const handle_t&)>& iteratee,
                                      bool parallel = false) const override;
    virtual bool for_each_path_handle_impl(const std::function<bool(const path_handle_t&)>& iteratee) const override;
    virtual bool for_each_step_on_handle_impl(const handle_t& handle,
                                              const std::function<bool(const step_handle_t&)>& iteratee) const override;

private:

    /// The graph that everything is forwarded to. Not owned.
    const PathHandleGraph* base_ = nullptr;

    /// The path that this graph knows positions on.
    path_handle_t target_path_;

    /// The length of the target path in bases.
    size_t target_path_length_ = 0;

    /// The offset of the first base of each step with a known position. It
    /// holds the anchors' first and last steps, path_begin at 0, path_end at
    /// the path length, and every step a position query has walked past.
    mutable std::unordered_map<step_handle_t, size_t> step_pos_;

    /// The known positions from construction, as (position, step) pairs sorted
    /// by position, for get_step_at_position to binary search. Steps learned
    /// by later queries are only in step_pos_.
    std::vector<std::pair<size_t, step_handle_t>> sorted_known_;

    /// The steps of the target path inside the anchors' runs, by node ID.
    std::unordered_map<nid_t, std::vector<step_handle_t>> node_target_steps_;
};

} // namespace vg

#endif // VG_ANCHOR_BACKED_POSITION_GRAPH_HPP_INCLUDED
