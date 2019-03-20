#include "gapless_extender.hpp"

#include <queue>
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

// Convert the gapless match to a path and a mismatch count.
std::pair<Path, size_t> gapless_match_to_path(GaplessMatch& match, const GBWTGraph& graph, const std::string& sequence) {

    std::pair<Path, size_t> result(Path(), match.score);
    if (match.empty()) {
        return result;
    }

    std::sort(match.mismatches.begin(), match.mismatches.end());
    size_t sequence_offset = 0; // Start of the unmapped part in the sequence.
    size_t mismatch_offset = 0; // In match.mismatches.
    size_t node_offset = match.offset; // Start of the alignment in the current node.
    for (size_t i = 0; i < match.path.size(); i++) {
        size_t limit = sequence_offset + graph.get_length(match.path[i]) - node_offset;
        Mapping& mapping = *(result.first.add_mapping());
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

std::pair<Path, size_t> GaplessExtender::extend_seeds(std::vector<std::pair<size_t, pos_t>>& cluster, const std::string& sequence, size_t max_mismatches) {

    GaplessMatch best_match {
        max_mismatches + 1,
        static_cast<size_t>(0), static_cast<size_t>(0),
        static_cast<size_t>(0),
        { },
        { }
    };
    if (this->graph == nullptr) {
        return gapless_match_to_path(best_match, *(this->graph), sequence);
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
                return gapless_match_to_path(best_match, *(this->graph), sequence);
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
                return gapless_match_to_path(best_match, *(this->graph), sequence);
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
                return gapless_match_to_path(best_match, *(this->graph), sequence);
            }
        }
    }

    return gapless_match_to_path(best_match, *(this->graph), sequence);
}

//------------------------------------------------------------------------------

} // namespace vg
