#include "anchor_backed_position_graph.hpp"

#include <algorithm>
#include <cassert>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace vg {

thread_local AnchorGraphStats g_anchor_graph_stats;

AnchorBackedPositionGraph::AnchorBackedPositionGraph(
    const PathHandleGraph* base,
    const std::vector<PrecomputedAnchor>& anchors,
    path_handle_t target_path,
    size_t target_path_length)
    : base_(base),
      target_path_(target_path),
      target_path_length_(target_path_length)
{
    if (!base_) {
        throw std::invalid_argument("AnchorBackedPositionGraph: base graph pointer is null");
    }

    // Seed the cache with anchor step→position pairs. We treat each anchor's
    // begin and end steps as known boundaries. Duplicate entries (e.g.
    // adjacent anchors sharing a step) collapse cleanly via map semantics.
    step_pos_.reserve(anchors.size() * 2 + 2);
    for (const auto& a : anchors) {
        step_pos_[a.step_begin] = a.path_offset_step_begin;
        step_pos_[a.step_end]   = a.path_offset_step_end;
    }

    // Path-begin sentinel. The first step on the target path is always at
    // base position 0; cache it so backward queries / position-zero walks
    // have an anchor.
    step_handle_t path_first = base_->path_begin(target_path_);
    step_pos_[path_first] = 0;

    // Path-end sentinel. path_end is past-the-end and not a real step; we
    // cache its "position" as the path length so walks toward the end can
    // detect termination without a separate iteration check.
    step_handle_t path_past_end = base_->path_end(target_path_);
    step_pos_[path_past_end] = target_path_length_;

    // Build sorted_known_ for binary search. Include all map entries.
    sorted_known_.reserve(step_pos_.size());
    for (const auto& kv : step_pos_) {
        sorted_known_.emplace_back(kv.second, kv.first);
    }
    std::sort(sorted_known_.begin(), sorted_known_.end(),
              [](const std::pair<size_t, step_handle_t>& a,
                 const std::pair<size_t, step_handle_t>& b) {
                  return a.first < b.first;
              });

    // Precompute node -> target-path steps over the read's overlap region.
    // We walk the target path from each anchor's step_begin to its step_end
    // (that is exactly the stretch of the target path the read covers) and
    // record every step. for_each_step_on_handle then returns just these,
    // instead of asking the base graph to enumerate every haplotype on the
    // node. The walk is bounded by the read's span, not the path length.
    for (const auto& a : anchors) {
        step_handle_t step = a.step_begin;
        // Safety cap so a bad anchor can never spin: the region can't be longer
        // than the whole path.
        for (size_t guard = 0; guard <= target_path_length_ + 1; ++guard) {
            nid_t node = base_->get_id(base_->get_handle_of_step(step));
            node_target_steps_[node].push_back(step);
            if (step == a.step_end) {
                break;
            }
            if (!base_->has_next_step(step)) {
                break;
            }
            step = base_->get_next_step(step);
        }
    }
}

// ============================================================
// PathPositionHandleGraph overrides
// ============================================================

size_t AnchorBackedPositionGraph::get_path_length(const path_handle_t& path_handle) const {
    if (!(path_handle == target_path_)) {
        std::ostringstream msg;
        msg << "AnchorBackedPositionGraph::get_path_length called for non-target path '"
            << base_->get_path_name(path_handle) << "' (target is '"
            << base_->get_path_name(target_path_) << "')";
        throw std::logic_error(msg.str());
    }
    return target_path_length_;
}

size_t AnchorBackedPositionGraph::get_position_of_step(const step_handle_t& step) const {
    ++g_anchor_graph_stats.pos_of_step_calls;

    // Cache hit?
    auto cached = step_pos_.find(step);
    if (cached != step_pos_.end()) {
        return cached->second;
    }

    // Verify the queried step is on our target path. If not, we cannot
    // answer (the underlying base_ is a plain PathHandleGraph; it doesn't
    // know positions either).
    if (!(base_->get_path_handle_of_step(step) == target_path_)) {
        std::ostringstream msg;
        msg << "AnchorBackedPositionGraph::get_position_of_step: step is not on the target path '"
            << base_->get_path_name(target_path_) << "'";
        throw std::logic_error(msg.str());
    }

    // Walk BACKWARD from `step` to the nearest step whose position we already
    // know — an anchor boundary, the path-begin sentinel (position 0), or a
    // previously-memoized step — accumulating node lengths as we go. Then
    // position(step) = known_pos + accumulated.
    //
    // This is O(distance to the nearest known step), NOT O(path). The previous
    // implementation re-walked from path_begin on every cache miss and ignored
    // the already-cached prefix, so a Surjector issuing many position queries
    // with sparse anchors cost O(path^2). Walking backward and backfilling every
    // step we pass keeps each query bounded by the gap to the nearest known
    // position, and that gap only shrinks as the cache fills.
    //
    // `accumulated` for a cursor is the summed length of steps in [cursor, step)
    // (0 when cursor == step). trail records each visited (cursor, accumulated)
    // so we can backfill all of them once we find a known anchor.
    std::vector<std::pair<step_handle_t, size_t>> trail;
    step_handle_t cursor = step;
    size_t accumulated = 0;
    size_t known_pos = 0;

    while (true) {
        auto it = step_pos_.find(cursor);
        if (it != step_pos_.end()) {
            known_pos = it->second;
            break;
        }
        trail.emplace_back(cursor, accumulated);
        if (!base_->has_previous_step(cursor)) {
            // First step on the path; its position is 0. (path_begin is normally
            // already cached at 0, so this is just a safety net.)
            known_pos = 0;
            break;
        }
        step_handle_t prev = base_->get_previous_step(cursor);
        accumulated += base_->get_length(base_->get_handle_of_step(prev));
        cursor = prev;
        ++g_anchor_graph_stats.pos_of_step_walk;
    }

    const size_t acc_end = accumulated;
    // Backfill every step we visited: position(s) = known_pos + (acc_end - acc_s).
    for (const auto& e : trail) {
        step_pos_.emplace(e.first, known_pos + (acc_end - e.second));
    }
    return known_pos + acc_end;
}

step_handle_t AnchorBackedPositionGraph::get_step_at_position(
    const path_handle_t& path, const size_t& position) const
{
    ++g_anchor_graph_stats.step_at_pos_calls;
    if (!(path == target_path_)) {
        std::ostringstream msg;
        msg << "AnchorBackedPositionGraph::get_step_at_position called for non-target path '"
            << base_->get_path_name(path) << "' (target is '"
            << base_->get_path_name(target_path_) << "')";
        throw std::logic_error(msg.str());
    }

    // Past the end → path_end sentinel, as the interface requires.
    if (position >= target_path_length_) {
        return base_->path_end(target_path_);
    }

    // Binary search sorted_known_ for the largest known position ≤ `position`.
    // `upper_bound(position)` gives the first entry with pos > position;
    // step back one to get the floor.
    auto it = std::upper_bound(
        sorted_known_.begin(), sorted_known_.end(),
        position,
        [](size_t value, const std::pair<size_t, step_handle_t>& entry) {
            return value < entry.first;
        });
    // sorted_known_ always contains (0, path_begin), so `it != begin()`
    // is guaranteed for any non-negative position.
    assert(it != sorted_known_.begin());
    --it;

    step_handle_t cursor = it->second;
    size_t cursor_pos = it->first;

    // Walk forward until the next step would exceed `position`. The current
    // step is the answer when (cursor_pos ≤ position < cursor_pos + node_len).
    while (true) {
        // Sentinels (path_begin and path_end) have associated positions but
        // aren't real steps for length lookup. path_end never matches because
        // its position == target_path_length_ > position (we returned above).
        // path_begin's handle IS a real step on the path, so it walks normally.
        handle_t h = base_->get_handle_of_step(cursor);
        size_t cursor_end = cursor_pos + base_->get_length(h);
        if (position < cursor_end) {
            // Memoize the walk-target so the inverse query is also free.
            step_pos_.emplace(cursor, cursor_pos);
            return cursor;
        }
        cursor_pos = cursor_end;
        if (!base_->has_next_step(cursor)) {
            // Walked off the end — should not happen because position < length.
            return base_->path_end(target_path_);
        }
        // Memoize every step we pass.
        step_pos_.emplace(cursor, cursor_pos - base_->get_length(h));
        cursor = base_->get_next_step(cursor);
        ++g_anchor_graph_stats.step_at_pos_walk;
    }
}

// ============================================================
// HandleGraph delegates
// ============================================================

bool AnchorBackedPositionGraph::has_node(nid_t node_id) const {
    return base_->has_node(node_id);
}

handle_t AnchorBackedPositionGraph::get_handle(const nid_t& node_id, bool is_reverse) const {
    return base_->get_handle(node_id, is_reverse);
}

nid_t AnchorBackedPositionGraph::get_id(const handle_t& handle) const {
    return base_->get_id(handle);
}

bool AnchorBackedPositionGraph::get_is_reverse(const handle_t& handle) const {
    return base_->get_is_reverse(handle);
}

handle_t AnchorBackedPositionGraph::flip(const handle_t& handle) const {
    return base_->flip(handle);
}

size_t AnchorBackedPositionGraph::get_length(const handle_t& handle) const {
    return base_->get_length(handle);
}

std::string AnchorBackedPositionGraph::get_sequence(const handle_t& handle) const {
    return base_->get_sequence(handle);
}

size_t AnchorBackedPositionGraph::get_node_count() const {
    return base_->get_node_count();
}

nid_t AnchorBackedPositionGraph::min_node_id() const {
    return base_->min_node_id();
}

nid_t AnchorBackedPositionGraph::max_node_id() const {
    return base_->max_node_id();
}

bool AnchorBackedPositionGraph::follow_edges_impl(const handle_t& handle, bool go_left,
                                                  const std::function<bool(const handle_t&)>& iteratee) const {
    return base_->follow_edges(handle, go_left, iteratee);
}

bool AnchorBackedPositionGraph::for_each_handle_impl(const std::function<bool(const handle_t&)>& iteratee,
                                                    bool parallel) const {
    return base_->for_each_handle(iteratee, parallel);
}

// ============================================================
// PathHandleGraph delegates
// ============================================================

size_t AnchorBackedPositionGraph::get_path_count() const {
    return base_->get_path_count();
}

bool AnchorBackedPositionGraph::has_path(const std::string& path_name) const {
    return base_->has_path(path_name);
}

path_handle_t AnchorBackedPositionGraph::get_path_handle(const std::string& path_name) const {
    return base_->get_path_handle(path_name);
}

std::string AnchorBackedPositionGraph::get_path_name(const path_handle_t& path_handle) const {
    return base_->get_path_name(path_handle);
}

bool AnchorBackedPositionGraph::get_is_circular(const path_handle_t& path_handle) const {
    return base_->get_is_circular(path_handle);
}

size_t AnchorBackedPositionGraph::get_step_count(const path_handle_t& path_handle) const {
    return base_->get_step_count(path_handle);
}

handle_t AnchorBackedPositionGraph::get_handle_of_step(const step_handle_t& step_handle) const {
    return base_->get_handle_of_step(step_handle);
}

path_handle_t AnchorBackedPositionGraph::get_path_handle_of_step(const step_handle_t& step_handle) const {
    return base_->get_path_handle_of_step(step_handle);
}

step_handle_t AnchorBackedPositionGraph::path_begin(const path_handle_t& path_handle) const {
    return base_->path_begin(path_handle);
}

step_handle_t AnchorBackedPositionGraph::path_end(const path_handle_t& path_handle) const {
    return base_->path_end(path_handle);
}

step_handle_t AnchorBackedPositionGraph::path_back(const path_handle_t& path_handle) const {
    return base_->path_back(path_handle);
}

step_handle_t AnchorBackedPositionGraph::path_front_end(const path_handle_t& path_handle) const {
    return base_->path_front_end(path_handle);
}

bool AnchorBackedPositionGraph::has_next_step(const step_handle_t& step_handle) const {
    return base_->has_next_step(step_handle);
}

bool AnchorBackedPositionGraph::has_previous_step(const step_handle_t& step_handle) const {
    return base_->has_previous_step(step_handle);
}

step_handle_t AnchorBackedPositionGraph::get_next_step(const step_handle_t& step_handle) const {
    return base_->get_next_step(step_handle);
}

step_handle_t AnchorBackedPositionGraph::get_previous_step(const step_handle_t& step_handle) const {
    return base_->get_previous_step(step_handle);
}

bool AnchorBackedPositionGraph::for_each_path_handle_impl(
    const std::function<bool(const path_handle_t&)>& iteratee) const {
    return base_->for_each_path_handle(iteratee);
}

bool AnchorBackedPositionGraph::for_each_step_on_handle_impl(
    const handle_t& handle,
    const std::function<bool(const step_handle_t&)>& iteratee) const {
    // Return ONLY the target path's steps on this node (precomputed over the
    // read's overlap region), not every haplotype's. The Surjector filters to
    // the target path anyway, so this is behavior-preserving for surjection but
    // avoids the all-haplotype GBWT enumeration that dominated the runtime.
    // A node with no entry is either off the target path or outside the read's
    // overlap region — in both cases the target path contributes no step the
    // Surjector would keep, so returning nothing matches the filtered result.
    auto it = node_target_steps_.find(base_->get_id(handle));
    if (it == node_target_steps_.end()) {
        return true;
    }
    for (const step_handle_t& step : it->second) {
        if (!iteratee(step)) {
            return false;
        }
    }
    return true;
}

// ============================================================
// Private helpers
// ============================================================

size_t AnchorBackedPositionGraph::walk_forward_to(
    step_handle_t start, size_t start_pos,
    const step_handle_t& target,
    size_t max_distance) const
{
    step_handle_t cursor = start;
    size_t cursor_pos = start_pos;
    size_t walked = 0;

    while (walked <= max_distance) {
        if (cursor == target) {
            step_pos_[target] = cursor_pos;
            return cursor_pos;
        }
        step_pos_.emplace(cursor, cursor_pos);

        handle_t h = base_->get_handle_of_step(cursor);
        size_t step_len = base_->get_length(h);
        cursor_pos += step_len;
        walked += step_len;

        if (!base_->has_next_step(cursor)) {
            break;
        }
        cursor = base_->get_next_step(cursor);
    }
    return std::numeric_limits<size_t>::max();
}

} // namespace vg
