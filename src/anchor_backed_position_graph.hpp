#ifndef VG_ANCHOR_BACKED_POSITION_GRAPH_HPP_INCLUDED
#define VG_ANCHOR_BACKED_POSITION_GRAPH_HPP_INCLUDED

/** \file
 * anchor_backed_position_graph.hpp
 *
 * A PathPositionHandleGraph implementation that answers position queries
 * from pre-computed anchors instead of a built positional index.
 *
 * Intended use: feed this to a Surjector instead of a bdsg::ReferencePathOverlay
 * built over the target haplotype. Surject's calls to `get_position_of_step`,
 * `get_step_at_position`, and `get_path_length` are answered from anchor
 * data; everything else delegates to the underlying gbwtgraph::GBWTGraph
 * (which only implements PathHandleGraph).
 *
 * Cost model:
 *   - Anchor-boundary step → O(1) hash-map lookup.
 *   - Interior step (between anchor boundaries) → O(k) walk via get_next_step
 *     from the nearest known step, accumulating node lengths. Memoized.
 *   - Position-to-step → O(log n) binary search + O(k) walk from the
 *     nearest known position.
 *
 * Restrictions:
 *   - Answers position queries ONLY for the single `target_path` passed to
 *     the constructor. Queries for any other path are programming errors;
 *     they assert in debug builds and throw in release builds.
 *   - Caller is responsible for the cached path length being correct.
 *     Typically the coordinate translator that built the anchors computes it
 *     by walking the target's GBWT path summing node lengths.
 */

#include "handle.hpp"

#include <handlegraph/path_position_handle_graph.hpp>

#include <unordered_map>
#include <utility>
#include <vector>

namespace vg {

using namespace std;

/// Minimal anchor payload used by AnchorBackedPositionGraph. Carries the
/// per-chunk step handles and their base offsets on the target path; this
/// is exactly what the position graph needs.
///
/// When the wider `surject_with_anchors` API is wired up, this struct can
/// grow additional fields (read range, graph_path, GBWT search states) —
/// AnchorBackedPositionGraph only reads `step_begin/end` and
/// `path_offset_step_begin/end`.
struct PrecomputedAnchor {
    /// First step on the target path that this anchor pins.
    step_handle_t step_begin;
    /// Last step (inclusive) on the target path that this anchor pins.
    step_handle_t step_end;
    /// Base position of `step_begin` on the target path (0-based, forward).
    size_t path_offset_step_begin = 0;
    /// Base position of `step_end` on the target path (0-based, forward).
    size_t path_offset_step_end = 0;
};

/// Per-thread diagnostics for the position-query cost. Reset before a surjection
/// and read afterwards to see how much walking the position queries did — the
/// difference between a healthy sparse-anchor run and an O(path^2) blow-up.
struct AnchorGraphStats {
    size_t pos_of_step_calls = 0;  ///< calls to get_position_of_step
    size_t pos_of_step_walk  = 0;  ///< total steps walked backward on cache miss
    size_t step_at_pos_calls = 0;  ///< calls to get_step_at_position
    size_t step_at_pos_walk  = 0;  ///< total steps walked forward on cache miss
};
extern thread_local AnchorGraphStats g_anchor_graph_stats;

/**
 * A PathPositionHandleGraph that delegates everything except position
 * queries to an underlying PathHandleGraph, and answers position queries
 * from a sparse set of pre-computed anchor step positions plus a precomputed
 * path length.
 */
class AnchorBackedPositionGraph : public PathPositionHandleGraph {
public:

    /// Construct over `base` (typically `&gbz.graph`), pinning a single
    /// target path. `anchors` carries the known step→position bindings
    /// (typically two per chunk: step_begin and step_end). `target_path_length`
    /// is the total base length of the target path, precomputed by the
    /// caller (e.g. by summing node lengths along the GBWT path).
    AnchorBackedPositionGraph(const PathHandleGraph* base,
                              const std::vector<PrecomputedAnchor>& anchors,
                              path_handle_t target_path,
                              size_t target_path_length);

    /// Default destructor.
    ~AnchorBackedPositionGraph() = default;

    // No default construction or copy — the cache and base pointer are
    // tightly coupled to the constructor arguments.
    AnchorBackedPositionGraph() = delete;
    AnchorBackedPositionGraph(const AnchorBackedPositionGraph&) = delete;
    AnchorBackedPositionGraph& operator=(const AnchorBackedPositionGraph&) = delete;

    //////////////////////////
    // HandleGraph interface — all delegate to base_.
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
    // PathHandleGraph interface — all delegate to base_.
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
    // PathPositionHandleGraph interface — the interesting overrides.
    //////////////////////////

    /// Returns the precomputed length for `target_path`. Asserts on any
    /// other path (we cannot answer length for unknown paths).
    virtual size_t get_path_length(const path_handle_t& path_handle) const override;

    /// Returns the position of `step` on the target path. Cache hit if step
    /// is at an anchor boundary; otherwise walks forward from the nearest
    /// known step accumulating node lengths, memoizes the result, and
    /// returns it. Asserts if step is not on the target path.
    virtual size_t get_position_of_step(const step_handle_t& step) const override;

    /// Returns the step containing `position` on the target path, or
    /// `path_end(target_path)` if `position >= path_length`. Asserts on any
    /// other path.
    virtual step_handle_t get_step_at_position(const path_handle_t& path,
                                               const size_t& position) const override;

    //////////////////////////
    // Inspection / testing helpers.
    //////////////////////////

    /// Number of step→position entries known (including memoized + sentinels).
    size_t cached_position_count() const { return step_pos_.size(); }

    /// The target path this adapter answers positions for.
    path_handle_t target_path() const { return target_path_; }

protected:

    /// Iteratee-impl methods required by the base class. Delegate to base_
    /// via its public template wrappers (which call back into base_'s impls).
    virtual bool follow_edges_impl(const handle_t& handle, bool go_left,
                                   const std::function<bool(const handle_t&)>& iteratee) const override;
    virtual bool for_each_handle_impl(const std::function<bool(const handle_t&)>& iteratee,
                                      bool parallel = false) const override;
    virtual bool for_each_path_handle_impl(const std::function<bool(const path_handle_t&)>& iteratee) const override;
    virtual bool for_each_step_on_handle_impl(const handle_t& handle,
                                              const std::function<bool(const step_handle_t&)>& iteratee) const override;

private:

    /// Underlying graph for handle / step / metadata queries. Not owned.
    const PathHandleGraph* base_ = nullptr;

    /// The single path we know positions for.
    path_handle_t target_path_;

    /// Total base length of `target_path_` — caller-precomputed.
    size_t target_path_length_ = 0;

    /// Known step → base-position table. Populated from anchors (two entries
    /// per anchor: step_begin / step_end) plus path_begin → 0 and
    /// path_end → target_path_length_ sentinels. Lazily extended via
    /// memoization when get_position_of_step walks to a previously-unknown
    /// step. Mutable because get_position_of_step is const.
    mutable std::unordered_map<step_handle_t, size_t> step_pos_;

    /// Initial known positions sorted ascending: (position, step). Built
    /// once at construction (sentinels + anchor entries). Used by
    /// get_step_at_position for binary search; NOT updated by memoization
    /// (the unordered_map handles cache hits regardless).
    std::vector<std::pair<size_t, step_handle_t>> sorted_known_;

    /// node id -> the target path's steps that visit it within the read's
    /// overlap region (the union of the anchors' [step_begin, step_end] spans),
    /// precomputed once at construction. for_each_step_on_handle returns these
    /// instead of enumerating EVERY haplotype's step on the node via the base
    /// graph — that enumeration is what made surjection cost
    /// O(#haplotypes-per-node x read-nodes) (~12s on dense graphs) even though
    /// the Surjector discards all non-target steps.
    std::unordered_map<nid_t, std::vector<step_handle_t>> node_target_steps_;

    /// Walk forward from `start` on the target path, accumulating node
    /// lengths starting from `start_pos`, until we reach `target` or pass
    /// `max_distance` bases. Returns the position when target is found, or
    /// numeric_limits<size_t>::max() if the walk exceeded max_distance.
    /// Memoizes every step it visits along the way.
    size_t walk_forward_to(step_handle_t start, size_t start_pos,
                           const step_handle_t& target,
                           size_t max_distance) const;
};

} // namespace vg

#endif // VG_ANCHOR_BACKED_POSITION_GRAPH_HPP_INCLUDED
