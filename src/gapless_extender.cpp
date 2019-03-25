#include "gapless_extender.hpp"

#include <queue>
#include <set>
#include <stack>

namespace vg {

//------------------------------------------------------------------------------

GaplessExtender::GaplessExtender() :
    graph(nullptr)
{
}

GaplessExtender::GaplessExtender(const GBWTGraph& graph) :
    graph(&graph)
{
}

//------------------------------------------------------------------------------

struct GaplessMatch {
    size_t score;
    size_t start, limit; // In the sequence.
    size_t offset; // In the initial node of the path.
    gbwt::BidirectionalState state;

    std::vector<handle_t> path;
    std::vector<size_t> mismatches; // Mismatching sequence positions in an arbitrary order.

    size_t length() const { return this->limit - this->start; }
    bool empty() const { return (this->length() == 0); }

    // Low-priority elements have many mismatches over a short range.
    bool operator<(const GaplessMatch& another) const {
        return ((this->score > another.score) ||
                (this->score == another.score && this->length() < another.length()));
    }
};

std::ostream& operator<<(std::ostream& out, const GaplessMatch& match) {
    out << "(" << match.score << ", (" << match.start << ", " << match.limit << "), ";
    gbwt::operator<<(out, match.state) << ")";
    return out;
}

// Match forward, starting from the given target offset.
void match_forward(const std::string& seq, std::pair<const char*, size_t> target, size_t target_offset,
                   GaplessMatch& match, size_t error_bound) {
    while (match.limit < seq.length() && target_offset < target.second) {
        if (seq[match.limit] != target.first[target_offset]) {
            match.score++;
            if (match.score < error_bound) {
                match.mismatches.push_back(match.limit);
            } else {
                return;
            }
        }
        match.limit++;
        target_offset++;
    }
}

// Match backward, starting from the end of the target and updating match offset.
void match_backward(const std::string& seq, std::pair<const char*, size_t> target,
                    GaplessMatch& match, size_t error_bound) {
    match.offset = target.second;
    while (match.start > 0 && match.offset > 0) {
        match.start--;
        match.offset--;
        if (seq[match.start] != target.first[match.offset]) {
            match.score++;
            if (match.score < error_bound) {
                match.mismatches.push_back(match.start);
            } else {
                return;
            }
        }
    }
}

// Convert the gapless match to a path.
Path gapless_match_to_path(GaplessMatch& match, const GBWTGraph& graph, const std::string& sequence) {

    Path result;
    if (match.empty()) {
        return result;
    }

    std::sort(match.mismatches.begin(), match.mismatches.end());
    size_t sequence_offset = 0; // Start of the unmapped part in the sequence.
    size_t mismatch_offset = 0; // In match.mismatches.
    size_t node_offset = match.offset; // Start of the alignment in the current node.
    for (size_t i = 0; i < match.path.size(); i++) {
        size_t limit = std::min(sequence_offset + graph.get_length(match.path[i]) - node_offset, match.limit);
        Mapping& mapping = *(result.add_mapping());
        mapping.mutable_position()->set_node_id(graph.get_id(match.path[i]));
        mapping.mutable_position()->set_offset(node_offset);
        mapping.mutable_position()->set_is_reverse(graph.get_is_reverse(match.path[i]));
        while (mismatch_offset < match.mismatches.size() && match.mismatches[mismatch_offset] < limit) {
            if (sequence_offset < match.mismatches[mismatch_offset]) {
                Edit& exact_match = *(mapping.add_edit());
                exact_match.set_from_length(match.mismatches[mismatch_offset] - sequence_offset);
                exact_match.set_to_length(match.mismatches[mismatch_offset] - sequence_offset);
            }
            Edit& mismatch = *(mapping.add_edit());
            mismatch.set_from_length(1);
            mismatch.set_to_length(1);
            mismatch.set_sequence(std::string(1, sequence[match.mismatches[mismatch_offset]]));
            sequence_offset = match.mismatches[mismatch_offset] + 1;
            mismatch_offset++;
        }
        if (sequence_offset < limit) {
            Edit& exact_match = *(mapping.add_edit());
            exact_match.set_from_length(limit - sequence_offset);
            exact_match.set_to_length(limit - sequence_offset);
            sequence_offset = limit;
        }
        mapping.set_rank(i + 1);
        node_offset = 0;
    }

    return result;
}

std::pair<Path, size_t> GaplessExtender::extend_seeds(std::vector<std::pair<size_t, pos_t>>& cluster, const std::string& sequence, size_t max_mismatches) const {

    GaplessMatch best_match {
        max_mismatches + 1,
        static_cast<size_t>(0), static_cast<size_t>(0),
        static_cast<size_t>(0),
        { },
        { }
    };
    if (this->graph == nullptr) {
        return std::make_pair(gapless_match_to_path(best_match, *(this->graph), sequence), best_match.score);
    }

    // Process the seeds in sorted order.
    std::pair<size_t, pos_t> prev(0, make_pos_t(0, false, 0));
    std::sort(cluster.begin(), cluster.end());
    for (size_t i = 0; i < cluster.size(); i++) {

        // Start matching as early in the initial node as possible.
        std::pair<size_t, pos_t> hit = cluster[i];
        size_t adjustment = std::min(static_cast<size_t>(offset(hit.second)), hit.first);
        get_offset(hit.second) -= adjustment;
        hit.first -= adjustment;
        if (hit == prev) {
            continue; // This seed was redundant.
        }
        prev = hit;

        // Match the initial node.
        std::priority_queue<GaplessMatch> forward, backward;
        {
            handle_t handle = GBWTGraph::node_to_handle(pos_to_gbwt(hit.second));
            GaplessMatch match {
                static_cast<size_t>(0),
                hit.first, hit.first,
                static_cast<size_t>(offset(hit.second)),
                this->graph->get_bd_state(handle),
                { },
                { }            
            };
            match_forward(sequence, this->graph->get_sequence_view(handle), match.offset, match, best_match.score);
            if (match.score >= best_match.score) { 
                continue;
            } else {
                match.path.push_back(handle);
            }
            if (match.limit >= sequence.length()) {
                if (match.start == 0) {
                    best_match = match;
                } else {
                    backward.push(match);
                }
            } else {
                forward.push(match);
            }
            if (best_match.score == 0) {
                return std::make_pair(gapless_match_to_path(best_match, *(this->graph), sequence), best_match.score);
            }
        }

        // Match forward over all paths.
        while (!forward.empty()) {
            GaplessMatch curr = forward.top();
            forward.pop();
            this->graph->follow_paths(curr.state, false, [&](const gbwt::BidirectionalState& next_state) -> bool {
                if (next_state.empty()) {
                    return true;
                }
                handle_t handle = GBWTGraph::node_to_handle(next_state.forward.node);
                GaplessMatch next {
                    curr.score,
                    curr.start, curr.limit,
                    curr.offset,
                    next_state,
                    { },
                    { }
                };
                match_forward(sequence, this->graph->get_sequence_view(handle), 0, next, best_match.score);
                if (next.score >= best_match.score) {
                    return true;
                } else {
                    next.path.reserve(curr.path.size() + 1);
                    next.path.insert(next.path.end(), curr.path.begin(), curr.path.end());
                    next.path.push_back(handle);
                    next.mismatches.insert(next.mismatches.end(), curr.mismatches.begin(), curr.mismatches.end());
                }
                if (next.limit >= sequence.length()) {
                    if (next.start == 0) {
                        best_match = next;
                    } else {
                        backward.push(next);
                    }
                } else {
                    forward.push(next);
                }
                return true;
            });
            if (best_match.score == 0) {
                return std::make_pair(gapless_match_to_path(best_match, *(this->graph), sequence), best_match.score);
            }
        }

        // Match backward over all paths.
        while (!backward.empty()) {
            GaplessMatch curr = backward.top();
            backward.pop();
            this->graph->follow_paths(curr.state, true, [&](const gbwt::BidirectionalState& next_state) -> bool {
                if (next_state.empty()) {
                    return true;
                }
                handle_t handle = GBWTGraph::node_to_handle(gbwt::Node::reverse(next_state.backward.node));
                GaplessMatch next {
                    curr.score,
                    curr.start, curr.limit,
                    curr.offset, // This will be replaced in match_backward().
                    next_state,
                    { },
                    { }
                };
                match_backward(sequence, this->graph->get_sequence_view(handle), next, best_match.score);
                if (next.score >= best_match.score) {
                    return true;
                } else {
                    next.path.reserve(curr.path.size() + 1);
                    next.path.push_back(handle);
                    next.path.insert(next.path.end(), curr.path.begin(), curr.path.end());
                    next.mismatches.insert(next.mismatches.end(), curr.mismatches.begin(), curr.mismatches.end());
                }
                if (next.start == 0) {
                    best_match = next;
                } else {
                    backward.push(next);
                }
                return true;
            });
            if (best_match.score == 0) {
                return std::make_pair(gapless_match_to_path(best_match, *(this->graph), sequence), best_match.score);
            }
        }
    }

    return std::make_pair(gapless_match_to_path(best_match, *(this->graph), sequence), best_match.score);
}

//------------------------------------------------------------------------------

struct MaximalGBWTMatch {
    size_t start, limit; // In the sequence.
    size_t offset; // In the initial node of the path.
    gbwt::BidirectionalState state;

    std::vector<handle_t> path;

    size_t length() const { return this->limit - this->start; }
    bool empty() const { return (this->length() == 0); }

    // Two matches are considered equivalent, if the substrings and search states are identical.
    bool operator<(const MaximalGBWTMatch& another) const {
        if (this->start != another.start) {
            return (this->start < another.start);
        }
        if (this->limit != another.limit) {
            return (this->limit < another.limit);
        }
        if (this->state.forward.node != another.state.forward.node) {
            return (this->state.forward.node < another.state.forward.node);
        }
        if (this->state.backward.node != another.state.backward.node) {
            return (this->state.backward.node < another.state.backward.node);
        }
        if (this->state.forward.range.first != another.state.forward.range.first) {
            return (this->state.forward.range.first < another.state.forward.range.first);
        }
        if (this->state.forward.range.second != another.state.forward.range.second) {
            return (this->state.forward.range.second < another.state.forward.range.second);
        }
        return false;
    }
};

// Match forward, starting from the given target offset.
// Returns the number of successfully matched characters.
size_t match_forward(const std::string& seq, std::pair<const char*, size_t> target, size_t target_offset,
                   MaximalGBWTMatch& match) {
    size_t initial_offset = target_offset;
    while (match.limit < seq.length() && target_offset < target.second) {
        if (seq[match.limit] != target.first[target_offset]) {
            break;
        }
        match.limit++;
        target_offset++;
    }
    return target_offset - initial_offset;
}

// Match backward, starting from the end of the target and updating match offset.
// Returns the number of successfully matched characters.
size_t match_backward(const std::string& seq, std::pair<const char*, size_t> target, MaximalGBWTMatch& match) {
    match.offset = target.second;
    while (match.start > 0 && match.offset > 0) {
        match.start--;
        match.offset--;
        if (seq[match.start] != target.first[match.offset]) {
            return target.second - match.offset - 1;
        }
    }
    return target.second - match.offset;
}

// Convert the exact match to a path.
Path maximal_match_to_path(const MaximalGBWTMatch& match, const GBWTGraph& graph, const std::string& sequence) {

    Path result;
    if (match.empty()) {
        return result;
    }

    size_t sequence_offset = match.start; // Start of the unmapped part in the sequence.
    size_t node_offset = match.offset; // Start of the alignment in the current node.
    for (size_t i = 0; i < match.path.size(); i++) {
        size_t limit = std::min(sequence_offset + graph.get_length(match.path[i]) - node_offset, match.limit);

        Mapping& mapping = *(result.add_mapping());
        mapping.mutable_position()->set_node_id(graph.get_id(match.path[i]));
        mapping.mutable_position()->set_offset(node_offset);
        mapping.mutable_position()->set_is_reverse(graph.get_is_reverse(match.path[i]));

        Edit& exact_match = *(mapping.add_edit());
        exact_match.set_from_length(limit - sequence_offset);
        exact_match.set_to_length(limit - sequence_offset);
        sequence_offset = limit;

        mapping.set_rank(i + 1);
        node_offset = 0;
    }

    return result;
}

std::vector<std::pair<Path, size_t>> GaplessExtender::seeds_to_mems(std::vector<std::pair<size_t, pos_t>>& cluster, const std::string& sequence) const {

    // Process the seeds in sorted order.
    handle_t prev_handle = GBWTGraph::node_to_handle(gbwt::ENDMARKER);
    size_t prev_start = 0, prev_limit = 0;
    std::sort(cluster.begin(), cluster.end());
    std::set<MaximalGBWTMatch> mems;
    for (size_t i = 0; i < cluster.size(); i++) {

        // Skip redundant seeds.
        std::pair<size_t, pos_t> hit = cluster[i];
        handle_t handle = GBWTGraph::node_to_handle(pos_to_gbwt(hit.second));
        size_t adjustment = std::min(static_cast<size_t>(offset(hit.second)), hit.first);
        if (handle == prev_handle && offset(hit.second) - adjustment == prev_start && offset(hit.second) < prev_limit) {
            continue;
        }
        prev_handle = handle;
        prev_start = offset(hit.second) - adjustment;
        // prev_limit is updated later

        // Match the initial node.
        std::stack<MaximalGBWTMatch> forward, backward;
        {
            // Maximal match within the initial node.
            MaximalGBWTMatch match {
                hit.first, hit.first,
                static_cast<size_t>(offset(hit.second)),
                this->graph->get_bd_state(handle),
                { handle }            
            };
            std::pair<const char*, size_t> node_view = this->graph->get_sequence_view(handle);
            match_forward(sequence, node_view, match.offset, match);
            prev_limit = match.offset + match.length();
            match_backward(sequence, node_view, match);
            if (match.limit >= sequence.length()) {
                if (match.start == 0) {
                    mems.insert(match);
                } else if (match.offset == 0) {
                    backward.push(match);
                }
            } else if (prev_limit >= node_view.second) {
                forward.push(match);
            }
        }

        // Match forward over all paths.
        while (!forward.empty()) {
            MaximalGBWTMatch curr = forward.top();
            forward.pop();
            size_t extend_total = 0; // Number of paths in successful extensions.
            this->graph->follow_paths(curr.state, false, [&](const gbwt::BidirectionalState& next_state) -> bool {
                if (next_state.empty()) {
                    return true;
                }
                handle_t handle = GBWTGraph::node_to_handle(next_state.forward.node);
                MaximalGBWTMatch next {
                    curr.start, curr.limit,
                    curr.offset,
                    next_state,
                    { }
                };
                std::pair<const char*, size_t> node_view = this->graph->get_sequence_view(handle);
                size_t matching_chars = match_forward(sequence, node_view, 0, next);
                if (matching_chars == 0) {
                    return true;
                }
                extend_total += next_state.size();
                next.path.reserve(curr.path.size() + 1);
                next.path.insert(next.path.end(), curr.path.begin(), curr.path.end());
                next.path.push_back(handle);
                if (next.limit >= sequence.length()) {
                    if (next.start == 0) {
                        mems.insert(next);
                    } else if (next.offset == 0) {
                        backward.push(next);
                    }
                } else if (matching_chars >= node_view.second) {
                    forward.push(next);
                }
                return true;
            });
            // We could not extend all paths, so some MEMs end at curr.
            if (extend_total < curr.state.size() && curr.offset == 0) {
                backward.push(curr);
            }
        }

        // Match backward over all paths.
        while (!backward.empty()) {
            MaximalGBWTMatch curr = backward.top();
            backward.pop();
            size_t extend_total = 0; // Number of paths in successful extensions.
            this->graph->follow_paths(curr.state, true, [&](const gbwt::BidirectionalState& next_state) -> bool {
                if (next_state.empty()) {
                    return true;
                }
                handle_t handle = GBWTGraph::node_to_handle(gbwt::Node::reverse(next_state.backward.node));
                MaximalGBWTMatch next {
                    curr.start, curr.limit,
                    curr.offset, // This will be replaced in match_backward().
                    next_state,
                    { }
                };
                std::pair<const char*, size_t> node_view = this->graph->get_sequence_view(handle);
                size_t matching_chars = match_backward(sequence, node_view, next);
                if (matching_chars == 0) {
                    return true;
                }
                extend_total += next_state.size();
                next.path.reserve(curr.path.size() + 1);
                next.path.push_back(handle);
                next.path.insert(next.path.end(), curr.path.begin(), curr.path.end());
                if (next.start == 0) {
                    mems.insert(next);
                } else if (next.offset == 0) {
                    backward.push(next);
                }
                return true;
            });
            // We could not extend all paths, so some MEMs start at curr.
            if (extend_total < curr.state.size()) {
                mems.insert(curr);
            }
        }
    }

    // FIXME this would be a good place to filter MEMs contained in other MEMs
    std::vector<std::pair<Path, size_t>> result;
    for(const MaximalGBWTMatch& match : mems) {
        result.emplace_back(maximal_match_to_path(match, *(this->graph), sequence), match.start);
    }
    return result;
}

//------------------------------------------------------------------------------

} // namespace vg
