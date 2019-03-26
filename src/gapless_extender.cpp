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
    size_t seq_start, seq_limit; // Sequence range.
    size_t node_start, node_limit; // In the initial/final nodes of the path.
    gbwt::BidirectionalState state;

    std::vector<handle_t> path;

    size_t length() const { return this->seq_limit - this->seq_start; }
    bool empty() const { return (this->length() == 0); }

    // Compare sequence ranges, start/end nodes, and node offsets.
    bool operator<(const MaximalGBWTMatch& another) const {
        if (this->seq_start != another.seq_start) {
            return (this->seq_start < another.seq_start);
        }
        if (this->seq_limit != another.seq_limit) {
            return (this->seq_limit < another.seq_limit);
        }
        if (this->state.backward.node != another.state.backward.node) {
            return (this->state.backward.node < another.state.backward.node);
        }
        if (this->node_start != another.node_start) {
            return (this->node_start < another.node_start);
        }
        if (this->state.forward.node != another.state.forward.node) {
            return (this->state.forward.node < another.state.forward.node);
        }
        if (this->node_limit != another.node_limit) {
            return (this->node_limit < another.node_limit);
        }
        return false;
    }
};

// Match forward as long as the characters match,
void match_forward(const std::string& seq, std::pair<const char*, size_t> target, MaximalGBWTMatch& match) {
    while (match.seq_limit < seq.length() && match.node_limit < target.second) {
        if (seq[match.seq_limit] != target.first[match.node_limit]) {
            break;
        }
        match.seq_limit++;
        match.node_limit++;
    }
}

// Match backward as long as the characters match.
void match_backward(const std::string& seq, std::pair<const char*, size_t> target, MaximalGBWTMatch& match) {
    while (match.seq_start > 0 && match.node_start > 0) {
        if (seq[match.seq_start - 1] != target.first[match.node_start - 1]) {
            break;
        }
        match.seq_start--;
        match.node_start--;
    }
}

// Convert the exact match to a path.
Path maximal_match_to_path(const MaximalGBWTMatch& match, const GBWTGraph& graph) {

    Path result;
    if (match.empty()) {
        return result;
    }

    for (size_t i = 0; i < match.path.size(); i++) {
        size_t start = (i == 0 ? match.node_start : 0);
        size_t limit = (i + 1 == match.path.size() ? match.node_limit : graph.get_length(match.path[i]));

        Mapping& mapping = *(result.add_mapping());
        mapping.mutable_position()->set_node_id(graph.get_id(match.path[i]));
        mapping.mutable_position()->set_offset(start);
        mapping.mutable_position()->set_is_reverse(graph.get_is_reverse(match.path[i]));

        Edit& exact_match = *(mapping.add_edit());
        exact_match.set_from_length(limit - start);
        exact_match.set_to_length(limit - start);

        mapping.set_rank(i + 1);
    }

    return result;
}

std::vector<std::pair<Path, size_t>> GaplessExtender::maximal_extensions(std::vector<std::pair<size_t, pos_t>>& cluster, const std::string& sequence) const {

    // Process the seeds in sorted order.
    std::pair<size_t, pos_t> prev(0, make_pos_t(0, false, 0));
    size_t prev_limit = 0; // Limit in the initial node.
    std::sort(cluster.begin(), cluster.end());
    std::set<MaximalGBWTMatch> matches;
    for (size_t i = 0; i < cluster.size(); i++) {

        // Skip redundant seeds.
        std::pair<size_t, pos_t> normalized = cluster[i];
        size_t adjustment = std::min(static_cast<size_t>(offset(normalized.second)), normalized.first);
        normalized.first -= adjustment;
        get_offset(normalized.second) -= adjustment;
        if (normalized == prev && offset(cluster[i].second) < prev_limit) {
            continue;
        }
        prev = normalized;
        // prev_limit is updated later when we match the first node.

        // Match the initial node.
        handle_t handle = GBWTGraph::node_to_handle(pos_to_gbwt(cluster[i].second));
        MaximalGBWTMatch match {
            cluster[i].first, cluster[i].first,
            static_cast<size_t>(offset(cluster[i].second)), static_cast<size_t>(offset(cluster[i].second)),
            this->graph->get_bd_state(handle),
            { handle }            
        };
        std::pair<const char*, size_t> node_view = this->graph->get_sequence_view(handle);
        match_forward(sequence, node_view, match);
        prev_limit = match.node_limit;
        match_backward(sequence, node_view, match);

        // Match forward.
        while (match.node_limit >= node_view.second && match.seq_limit < sequence.length()) {
            bool extension = false, ambiguous = false;
            gbwt::BidirectionalState successor;
            this->graph->follow_paths(match.state, false, [&](const gbwt::BidirectionalState& next_state) -> bool {
                if (ambiguous) {
                    return false;
                }
                if (next_state.empty()) {
                    return true;
                }
                handle_t next_handle = GBWTGraph::node_to_handle(next_state.forward.node);
                if (this->graph->starts_with(next_handle, sequence[match.seq_limit])) {
                    if (extension) {
                        ambiguous = true;
                        return false;
                    } else {
                        extension = true;
                        successor = next_state;
                        handle = next_handle;
                        return true;
                    }
                }
                return true;
            });
            if (extension && !ambiguous) {
                node_view = this->graph->get_sequence_view(handle);
                match.seq_limit++;
                match.node_limit = 1;
                match.state = successor;
                match.path.push_back(handle);
                match_forward(sequence, node_view, match);
            } else {
                break;
            }
        }

        // Match backward.
        while (match.node_start == 0 && match.seq_start > 0) {
            bool extension = false, ambiguous = false;
            gbwt::BidirectionalState successor;
            this->graph->follow_paths(match.state, true, [&](const gbwt::BidirectionalState& next_state) -> bool {
                if (ambiguous) {
                    return false;
                }
                if (next_state.empty()) {
                    return true;
                }
                handle_t next_handle = GBWTGraph::node_to_handle(gbwt::Node::reverse(next_state.backward.node));
                if (this->graph->ends_with(next_handle, sequence[match.seq_start - 1])) {
                    if (extension) {
                        ambiguous = true;
                        return false;
                    } else {
                        extension = true;
                        successor = next_state;
                        handle = next_handle;
                        return true;
                    }
                }
                return true;
            });
            if (extension && !ambiguous) {
                node_view = this->graph->get_sequence_view(handle);
                match.seq_start--;
                match.node_start = node_view.second - 1;
                match.state = successor;
                match.path.insert(match.path.begin(), handle);
                match_backward(sequence, node_view, match);
            } else {
                break;
            }
        }

        if (!match.empty()) {
            matches.insert(match);
        }
    }

    // Convert the matches to Path objects.
    std::vector<std::pair<Path, size_t>> result;
    for(const MaximalGBWTMatch& match : matches) {
        result.emplace_back(maximal_match_to_path(match, *(this->graph)), match.seq_start);
    }
    return result;
}

//------------------------------------------------------------------------------

} // namespace vg
