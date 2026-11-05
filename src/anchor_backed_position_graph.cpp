#include "anchor_backed_position_graph.hpp"

#include <algorithm>
#include <cassert>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_set>
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

    // Learn the positions of the anchors' first and last steps. Anchors that
    // share a step must agree on its position.
    step_pos_.reserve(anchors.size() * 2 + 2);
    for (const auto& a : anchors) {
        step_pos_[a.step_begin] = a.path_offset_step_begin;
        step_pos_[a.step_end]   = a.path_offset_step_end;
    }

    // The path start is a known position, so a walk back always stops.
    step_handle_t path_first = base_->path_begin(target_path_);
    step_pos_[path_first] = 0;

    // path_end, one past the last step, sits at the path length.
    step_handle_t path_past_end = base_->path_end(target_path_);
    step_pos_[path_past_end] = target_path_length_;

    sorted_known_.reserve(step_pos_.size());
    for (const auto& kv : step_pos_) {
        sorted_known_.emplace_back(kv.second, kv.first);
    }
    std::sort(sorted_known_.begin(), sorted_known_.end(),
              [](const std::pair<size_t, step_handle_t>& a,
                 const std::pair<size_t, step_handle_t>& b) {
                  return a.first < b.first;
              });

    // Index the steps in each anchor's run by node, for
    // for_each_step_on_handle_impl. Runs can overlap, as when adjacent anchors
    // share an end step, so list each step once.
    std::unordered_set<step_handle_t> indexed;
    for (const auto& a : anchors) {
        step_handle_t step = a.step_begin;
        // A run has at most as many steps as the path has bases. The cap stops
        // the walk if step_end is not ahead of step_begin.
        for (size_t guard = 0; guard <= target_path_length_ + 1; ++guard) {
            if (indexed.insert(step).second) {
                nid_t node = base_->get_id(base_->get_handle_of_step(step));
                node_target_steps_[node].push_back(step);
            }
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

    auto cached = step_pos_.find(step);
    if (cached != step_pos_.end()) {
        return cached->second;
    }

    if (!(base_->get_path_handle_of_step(step) == target_path_)) {
        std::ostringstream msg;
        msg << "AnchorBackedPositionGraph::get_position_of_step: step is not on the target path '"
            << base_->get_path_name(target_path_) << "'";
        throw std::logic_error(msg.str());
    }

    // Walk back to the nearest step with a known position, adding up the
    // lengths of the steps between. Remembering every step passed makes later
    // queries nearby cheap.
    //
    // accumulated is the length of the steps from cursor up to, but not
    // including, step. trail holds (cursor, accumulated) for each step passed.
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
            // This is the first step of the path, at position 0.
            known_pos = 0;
            break;
        }
        step_handle_t prev = base_->get_previous_step(cursor);
        accumulated += base_->get_length(base_->get_handle_of_step(prev));
        cursor = prev;
        ++g_anchor_graph_stats.pos_of_step_walk;
    }

    const size_t acc_end = accumulated;
    // Remember the position of every step passed.
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

    if (position >= target_path_length_) {
        return base_->path_end(target_path_);
    }

    // Find the last known position at or before position.
    auto it = std::upper_bound(
        sorted_known_.begin(), sorted_known_.end(),
        position,
        [](size_t value, const std::pair<size_t, step_handle_t>& entry) {
            return value < entry.first;
        });
    // path_begin is known at 0, so there is always one.
    assert(it != sorted_known_.begin());
    --it;

    // Walk back from the next known position when it is closer. Before the
    // first anchor the last known position is the path start, so walking
    // forward would cross the whole path up to the read.
    auto next = it + 1;
    if (next != sorted_known_.end() && next->first < target_path_length_ &&
        next->first - position < position - it->first) {
        step_handle_t cursor = next->second;
        size_t cursor_pos = next->first;
        while (cursor_pos > position && base_->has_previous_step(cursor)) {
            cursor = base_->get_previous_step(cursor);
            cursor_pos -= base_->get_length(base_->get_handle_of_step(cursor));
            step_pos_.emplace(cursor, cursor_pos);
            ++g_anchor_graph_stats.step_at_pos_walk;
        }
        return cursor;
    }

    step_handle_t cursor = it->second;
    size_t cursor_pos = it->first;

    // Walk forward to the step that covers position. We never reach path_end,
    // since position is less than the path length.
    while (true) {
        handle_t h = base_->get_handle_of_step(cursor);
        size_t cursor_end = cursor_pos + base_->get_length(h);
        if (position < cursor_end) {
            step_pos_.emplace(cursor, cursor_pos);
            return cursor;
        }
        cursor_pos = cursor_end;
        if (!base_->has_next_step(cursor)) {
            // Only a wrong path length gets here.
            return base_->path_end(target_path_);
        }
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
    // Report only the target path's steps inside the anchors' runs. A
    // Surjector onto the target path keeps only target path steps, and the
    // runs cover every step the alignment visits, so it sees the same steps.
    // Listing every haplotype's steps on a node costs far more.
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

} // namespace vg
