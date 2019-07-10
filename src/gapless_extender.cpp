#include "gapless_extender.hpp"

#include <algorithm>
#include <queue>
#include <set>
#include <stack>

namespace vg {

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr size_t GaplessExtender::MAX_MISMATCHES;

//------------------------------------------------------------------------------

Position GaplessExtension::starting_position(const GBWTGraph& graph) const {
    Position position;
    if (this->empty()) {
        return position;
    }

    position.set_node_id(graph.get_id(this->path.front()));
    position.set_is_reverse(graph.get_is_reverse(this->path.front()));
    position.set_offset(this->offset);

    return position;
}

Position GaplessExtension::tail_position(const GBWTGraph& graph) const {
    Position position;
    if (this->empty()) {
        return position;
    }

    position.set_node_id(graph.get_id(this->path.back()));
    position.set_is_reverse(graph.get_is_reverse(this->path.back()));
    position.set_offset(this->tail_offset(graph));

    return position;
}

size_t GaplessExtension::tail_offset(const GBWTGraph& graph) const {
    size_t result = this->offset + this->length();
    for (size_t i = 0; i + 1 < this->path.size(); i++) {
        result -= graph.get_length(this->path[i]);
    }
    return result;
}

Path GaplessExtension::to_path(const GBWTGraph& graph, const std::string& sequence) const {

    Path result;

    auto mismatch = this->mismatch_positions.begin(); // The next mismatch.
    size_t read_offset = this->read_interval.first;   // Current offset in the read.
    size_t node_offset = this->offset;                // Current offset in the current node.
    for (size_t i = 0; i < this->path.size(); i++) {
        size_t limit = std::min(read_offset + graph.get_length(this->path[i]) - node_offset, this->read_interval.second);
        Mapping& mapping = *(result.add_mapping());
        mapping.mutable_position()->set_node_id(graph.get_id(this->path[i]));
        mapping.mutable_position()->set_offset(node_offset);
        mapping.mutable_position()->set_is_reverse(graph.get_is_reverse(this->path[i]));
        while (mismatch != this->mismatch_positions.end() && *mismatch < limit) {
            if (read_offset < *mismatch) {
                Edit& exact_match = *(mapping.add_edit());
                exact_match.set_from_length(*mismatch - read_offset);
                exact_match.set_to_length(*mismatch - read_offset);
            }
            Edit& edit = *(mapping.add_edit());
            edit.set_from_length(1);
            edit.set_to_length(1);
            edit.set_sequence(std::string(1, sequence[*mismatch]));
            read_offset = *mismatch + 1;
            ++mismatch;
        }
        if (read_offset < limit) {
            Edit& exact_match = *(mapping.add_edit());
            exact_match.set_from_length(limit - read_offset);
            exact_match.set_to_length(limit - read_offset);
            read_offset = limit;
        }
        mapping.set_rank(i + 1);
        node_offset = 0;
    }

    return result;
}

//------------------------------------------------------------------------------

GaplessExtender::GaplessExtender() :
    graph(nullptr), aligner(nullptr)
{
}

GaplessExtender::GaplessExtender(const GBWTGraph& graph, const Aligner& aligner) :
    graph(&graph), aligner(&aligner)
{
}

//------------------------------------------------------------------------------

template<class Element>
void in_place_subvector(std::vector<Element>& vec, size_t head, size_t tail) {
    if (head >= tail || tail > vec.size()) {
        vec.clear();
        return;
    }
    if (head > 0) {
        for (size_t i = head; i < tail; i++) {
            vec[i - head] = std::move(vec[i]);
        }
    }
    vec.resize(tail - head);
}

// Match the initial node, assuming that read_offset or node_offset is 0.
// Updates score, internal_score, and old_score but does not consider full-length bonuses.
void match_initial(GaplessExtension& match, const std::string& seq, std::pair<const char*, size_t> target, const Aligner* aligner) {
    size_t node_offset = match.offset;
    while (match.read_interval.second < seq.length() && node_offset < target.second) {
        if (seq[match.read_interval.second] != target.first[node_offset]) {
            match.score -= aligner->mismatch;
            match.internal_score++;
            match.old_score++;
        } else {
            match.score += aligner->match;
        }
        match.read_interval.second++;
        node_offset++;
    }
}

// Match forward but stop before the mismatch count reaches the limit.
// Updates score and internal_score but does not consider full-length bonuses.
// Returns the tail offset (the number of characters matched).
size_t match_forward(GaplessExtension& match, const std::string& seq, std::pair<const char*, size_t> target, uint32_t mismatch_limit, const Aligner* aligner) {
    size_t node_offset = 0;
    while (match.read_interval.second < seq.length() && node_offset < target.second) {
        if (seq[match.read_interval.second] != target.first[node_offset]) {
            if (match.internal_score + 1 >= mismatch_limit) {
                return node_offset;
            }
            match.score -= aligner->mismatch;
            match.internal_score++;
        } else {
            match.score += aligner->match;
        }
        match.read_interval.second++;
        node_offset++;
    }
    return node_offset;
}

// Match forward but stop before the mismatch count reaches the limit.
// Starts from the offset in the match and updates it.
// Updates score and internal_score but does not consider full-length bonuses.
void match_backward(GaplessExtension& match, const std::string& seq, std::pair<const char*, size_t> target, uint32_t mismatch_limit, const Aligner* aligner) {
    while (match.read_interval.first > 0 && match.offset > 0) {
        if (seq[match.read_interval.first - 1] != target.first[match.offset - 1]) {
            if (match.internal_score + 1 >= mismatch_limit) {
                return;
            }
            match.score -= aligner->mismatch;
            match.internal_score++;
        } else {
            match.score += aligner->match;
        }
        match.read_interval.first--;
        match.offset--;
    }
}

// Sort the extensions from left to right. Remove duplicates and empty extensions.
void remove_duplicates(std::vector<GaplessExtension>& result) {
    auto sort_order = [](const GaplessExtension& a, const GaplessExtension& b) -> bool {
        if (a.read_interval != b.read_interval) {
            return (a.read_interval < b.read_interval);
        }
        if (a.state.backward.node != b.state.backward.node) {
            return (a.state.backward.node < b.state.backward.node);
        }
        if (a.state.forward.node != b.state.forward.node) {
            return (a.state.forward.node < b.state.forward.node);
        }
        return (a.state.backward.range < b.state.backward.range);
    };
    std::sort(result.begin(), result.end(), sort_order);
    size_t tail = 0;
    for (size_t i = 0; i < result.size(); i++) {
        if (result[i].empty()) {
            continue;
        }
        if (tail == 0 || result[i].read_interval != result[tail - 1].read_interval || result[i].state != result[tail - 1].state) {
            if (i > tail) {
                result[tail] = std::move(result[i]);
            }
            tail++;
        }
    }
    result.resize(tail);
}

// Realign the extensions to find the mismatching positions.
void find_mismatches(const std::string& seq, const GBWTGraph& graph, std::vector<GaplessExtension>& result) {
    for (GaplessExtension& extension : result) {
        if (extension.internal_score == 0) {
            continue;
        }
        extension.mismatch_positions.reserve(extension.internal_score);
        size_t node_offset = extension.offset, read_offset = extension.read_interval.first;
        for (const handle_t& handle : extension.path) {
            std::pair<const char*, size_t> target = graph.get_sequence_view(handle);
            while (node_offset < target.second && read_offset < extension.read_interval.second) {
                if (target.first[node_offset] != seq[read_offset]) {
                    extension.mismatch_positions.push_back(read_offset);
                }
                node_offset++;
                read_offset++;
            }
            node_offset = 0;
        }
    }
}

//------------------------------------------------------------------------------

std::vector<GaplessExtension> GaplessExtender::extend(cluster_type& cluster, const std::string& sequence, size_t max_mismatches, bool trim_extensions) const {

    std::vector<GaplessExtension> result;
    if (this->graph == nullptr || this->aligner == nullptr || sequence.empty()) {
        return result;
    }

    // Find either the best extension for each seed or the best full-length alignment
    // for the entire cluster. If we have found a full-length alignment with
    // at most max_mismatches mismatches, we are no longer interested in extensions with
    // at least that many mismatches.
    bool full_length_found = false;
    uint32_t full_length_mismatches = std::numeric_limits<uint32_t>::max();
    for (seed_type seed : cluster) {
        GaplessExtension best_match {
            { }, static_cast<size_t>(0), gbwt::BidirectionalState(),
            { static_cast<size_t>(0), static_cast<size_t>(0) }, { },
            std::numeric_limits<int32_t>::min(), false, false,
            false, false, std::numeric_limits<uint32_t>::max(), std::numeric_limits<uint32_t>::max()
        };
        bool best_match_is_full_length = false;

        // Match the initial node and add it to the queue, unless we already have
        // at least as good full-length alignment,
        std::priority_queue<GaplessExtension> extensions;
        {
            size_t read_offset = get_read_offset(seed);
            size_t node_offset = get_node_offset(seed);
            GaplessExtension match {
                { seed.first }, node_offset, this->graph->get_bd_state(seed.first),
                { read_offset, read_offset }, { },
                static_cast<int32_t>(0), false, false,
                false, false, static_cast<uint32_t>(0), static_cast<uint32_t>(0)
            };
            match_initial(match, sequence, this->graph->get_sequence_view(seed.first), this->aligner);
            if (match.internal_score >= full_length_mismatches) {
                continue;
            }
            if (match.read_interval.first == 0) {
                match.left_full = true;
                match.left_maximal = true;
                match.score += this->aligner->full_length_bonus;
            }
            if (match.read_interval.second >= sequence.length()) {
                match.right_full = true;
                match.right_maximal = true;
                match.score += this->aligner->full_length_bonus;
            }
            extensions.push(match);
        }

        // Extend the most promising extensions first, using alignment scores for priority.
        // First make the extension right-maximal and then left-maximal.
        while (!extensions.empty()) {
            GaplessExtension curr = extensions.top();
            extensions.pop();
            if (curr.internal_score >= full_length_mismatches) {
                continue;
            }
            // Always allow at least max_mismatches / 2 mismatches in the current flank.
            uint32_t mismatch_limit = std::max(
                static_cast<uint32_t>(max_mismatches + 1),
                static_cast<uint32_t>(max_mismatches / 2 + curr.old_score + 1));
            mismatch_limit = std::min(mismatch_limit, full_length_mismatches);
            bool found_extension = false;

            // Case 1: Extend to the right.
            if (!curr.right_maximal) {
                this->graph->follow_paths(curr.state, false, [&](const gbwt::BidirectionalState& next_state) -> bool {
                    if (next_state.empty()) {
                        return true;
                    }
                    handle_t handle = GBWTGraph::node_to_handle(next_state.forward.node);
                    GaplessExtension next {
                        { }, curr.offset, next_state,
                        curr.read_interval, { },
                        curr.score, curr.left_full, curr.right_full,
                        curr.left_maximal, curr.right_maximal, curr.internal_score, curr.old_score
                    };
                    size_t node_offset = match_forward(next, sequence, this->graph->get_sequence_view(handle), mismatch_limit, this->aligner);
                    if (node_offset == 0) { // Did not match anything.
                        return true;
                    }
                    next.path.reserve(curr.path.size() + 1);
                    next.path.insert(next.path.end(), curr.path.begin(), curr.path.end());
                    next.path.push_back(handle);
                    // Did the extension become right-maximal?
                    if (next.read_interval.second >= sequence.length()) {
                        next.right_full = true;
                        next.right_maximal = true;
                        next.score += this->aligner->full_length_bonus;
                        next.old_score = next.internal_score;
                    } else if (node_offset < this->graph->get_length(handle)) {
                        next.right_maximal = true;
                        next.old_score = next.internal_score;
                    }
                    extensions.push(next);
                    found_extension = true;
                    return true;
                });
                if (!found_extension) {
                    curr.right_maximal = true;
                    curr.old_score = curr.internal_score;
                    extensions.push(curr);
                }
            }

            // Case 2: Extend to the left.
            else if (!curr.left_maximal) {
                this->graph->follow_paths(curr.state, true, [&](const gbwt::BidirectionalState& next_state) -> bool {
                    if (next_state.empty()) {
                        return true;
                    }
                    handle_t handle = GBWTGraph::node_to_handle(gbwt::Node::reverse(next_state.backward.node));
                    size_t node_length = this->graph->get_length(handle);
                    GaplessExtension next {
                        { }, node_length, next_state,
                        curr.read_interval, { },
                        curr.score, curr.left_full, curr.right_full,
                        curr.left_maximal, curr.right_maximal, curr.internal_score, curr.old_score
                    };
                    match_backward(next, sequence, this->graph->get_sequence_view(handle), mismatch_limit, this->aligner);
                    if (next.offset >= node_length) { // Did not match anything.
                        return true;
                    }
                    next.path.reserve(curr.path.size() + 1);
                    next.path.push_back(handle);
                    next.path.insert(next.path.end(), curr.path.begin(), curr.path.end());
                    // Did the extension become left-maximal?
                    if (next.read_interval.first == 0) {
                        next.left_full = true;
                        next.left_maximal = true;
                        next.score += this->aligner->full_length_bonus;
                        next.old_score = next.internal_score;
                    } else if (next.offset > 0) {
                        next.left_maximal = true;
                        next.old_score = next.internal_score;
                    }
                    extensions.push(next);
                    found_extension = true;
                    return true;
                });
                if (!found_extension) {
                    curr.left_maximal = true;
                    curr.old_score = curr.internal_score;
                    extensions.push(curr);
                }
            }

            // Case 3: Maximal extension with a better score than the best extension so far.
            else if (best_match < curr) {
                best_match = std::move(curr);
                if (best_match.full() && best_match.internal_score <= max_mismatches) {
                    full_length_found = true;
                    full_length_mismatches = best_match.internal_score;
                    best_match_is_full_length = true;
                }
            }
        }

        // Handle the best match.
        if (best_match.score > 0) {
            if (best_match_is_full_length) {
                result.clear();
                result.push_back(best_match);
                if (best_match.internal_score == 0) {
                    break;
                }
            } else if (!full_length_found && !best_match.empty()) {
                result.push_back(best_match);
            }
        }
    }

    // Remove duplicates, find mismatches, and trim mismatches to maximize score.
    // If we have a full-length alignment with sufficiently few mismatches, we do
    // not trim it.
    remove_duplicates(result);
    find_mismatches(sequence, *(this->graph), result);
    if (trim_extensions) {
        this->trim(result, max_mismatches);
    }

    return result;
}

//------------------------------------------------------------------------------

// Trim mismatches from the extension to maximize the score. Returns true if the
// extension was trimmed.
bool trim_mismatches(GaplessExtension& extension, const GBWTGraph& graph, const Aligner& aligner) {

    bool trimmed = false;

    // Trim mismatches from the tail.
    size_t best_trim = 0;
    std::pair<size_t, size_t> new_interval = extension.read_interval;
    int32_t best_score = extension.score;
    int32_t full_length = (extension.right_full ? aligner.full_length_bonus : 0);
    for (size_t trim = 1; trim <= extension.mismatches(); trim++) {
        size_t new_end = extension.mismatch_positions[extension.mismatches() - trim];
        size_t trim_length = extension.read_interval.second - new_end;
        int32_t score = extension.score - (trim_length - trim) * aligner.match + trim * aligner.mismatch - full_length;
        if (score > best_score) {
            best_trim = trim;
            new_interval.second = new_end;
            best_score = score;
        }
    }
    if (best_trim > 0) {
        extension.mismatch_positions.resize(extension.mismatches() - best_trim);
        extension.score = best_score;
        trimmed = true;
    }

    // Trim mismatches from the head.
    best_trim = 0;
    full_length = (extension.left_full ? aligner.full_length_bonus : 0);
    for (size_t trim = 1; trim <= extension.mismatches(); trim++) {
        size_t new_begin = extension.mismatch_positions[trim - 1] + 1;
        size_t trim_length = new_begin - extension.read_interval.first;
        int32_t score = extension.score - (trim_length - trim) * aligner.match + trim * aligner.mismatch - full_length;
        if (score > best_score) {
            best_trim = trim;
            new_interval.first = new_begin;
            best_score = score;
        }
    }
    if (best_trim > 0) {
        in_place_subvector(extension.mismatch_positions, best_trim, extension.mismatches());
        extension.score = best_score;
        trimmed = true;
    }

    // We have already updated the mismatch positions. Update all other fields
    // if necessary. If the extension became empty, simply clear the path.
    if (new_interval.first >= new_interval.second) {
        extension.read_interval = new_interval;
        extension.left_full = extension.right_full = false;
        extension.path.clear();
    } else if (new_interval != extension.read_interval) {
        if (new_interval.first > extension.read_interval.first) {
            extension.left_full = false;
        }
        if (new_interval.second < extension.read_interval.second) {
            extension.right_full = false;
        }
        size_t node_offset = extension.offset, read_offset = extension.read_interval.first;
        extension.read_interval = new_interval;
        size_t path_head = 0;
        while (path_head < extension.path.size()) {
            size_t node_length = graph.get_length(extension.path[path_head]);
            read_offset += node_length - node_offset;
            node_offset = 0;
            if (read_offset > new_interval.first) {
                extension.offset = node_length - (read_offset - new_interval.first);
                break;
            }
            path_head++;
        }
        size_t path_tail = path_head + 1;
        while (read_offset < new_interval.second) {
            read_offset += graph.get_length(extension.path[path_tail]);
            path_tail++;
        }
        // Adjust the path and determine the new search state.
        if (path_head > 0 || path_tail < extension.path.size()) {
            in_place_subvector(extension.path, path_head, path_tail);
            extension.state = graph.bd_find(extension.path);
        }
    }

    return trimmed;
}

void GaplessExtender::trim(std::vector<GaplessExtension>& extensions, size_t max_mismatches) const {
    bool trimmed = false;
    for (GaplessExtension& extension : extensions) {
        if (!extension.full() || extension.mismatches() > max_mismatches) {
            trimmed |= trim_mismatches(extension, *(this->graph), *(this->aligner));
        }
    }
    if (trimmed) {
        remove_duplicates(extensions);
    }
}

//------------------------------------------------------------------------------

} // namespace vg
