#include "gapless_extender.hpp"

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
    size_t result = this->offset + this->core_length();
    for (size_t i = 0; i + 1 < this->path.size(); i++) {
        result -= graph.get_length(this->path[i]);
    }
    return result;
}

Path GaplessExtension::to_path(const GBWTGraph& graph, const std::string& sequence) const {

    Path result;

    // Skip mismatches before the core interval.
    auto mismatch = this->mismatch_positions.begin();
    size_t sequence_offset = this->core_interval.first; // Start of the unmapped part in the sequence.
    while (mismatch != this->mismatch_positions.end() && *mismatch < sequence_offset) {
        ++mismatch;
    }

    size_t node_offset = this->offset; // Start of the alignment in the current node.
    for (size_t i = 0; i < this->path.size(); i++) {
        size_t limit = std::min(sequence_offset + graph.get_length(this->path[i]) - node_offset, this->core_interval.second);
        Mapping& mapping = *(result.add_mapping());
        mapping.mutable_position()->set_node_id(graph.get_id(this->path[i]));
        mapping.mutable_position()->set_offset(node_offset);
        mapping.mutable_position()->set_is_reverse(graph.get_is_reverse(this->path[i]));
        while (mismatch != this->mismatch_positions.end() && *mismatch < limit) {
            if (sequence_offset < *mismatch) {
                Edit& exact_match = *(mapping.add_edit());
                exact_match.set_from_length(*mismatch - sequence_offset);
                exact_match.set_to_length(*mismatch - sequence_offset);
            }
            Edit& edit = *(mapping.add_edit());
            edit.set_from_length(1);
            edit.set_to_length(1);
            edit.set_sequence(std::string(1, sequence[*mismatch]));
            sequence_offset = *mismatch + 1;
            ++mismatch;
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

std::vector<GaplessExtension> GaplessExtender::extend(cluster_type& cluster, const std::string& sequence, size_t max_mismatches) const {

    if (this->graph == nullptr) {
        return std::vector<GaplessExtension>();
    }

    // Try to find a full-length alignment.
    GaplessExtension full_length = this->extend_seeds(cluster, sequence, max_mismatches, false);
    if (full_length.full()) {
        return { full_length };
    }

    // Find maximal unambiguous extensions and extend their flanks.
    std::vector<GaplessExtension> result = this->maximal_extensions(cluster, sequence, true);
    this->extend_flanks(result, sequence, max_mismatches / 2);

    return result;
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

// Convert the GaplessMatch to GaplessExtension.
GaplessExtension match_to_extension(GaplessMatch& match, const std::string& sequence) {

    GaplessExtension result {
        { },
        match.offset,
        match.state,
        { match.start, match.limit },
        (match.start == 0 && match.limit == sequence.length()),
        { match.start, match.limit },
        { }
    };

    result.path.swap(match.path);
    std::sort(match.mismatches.begin(), match.mismatches.end());
    result.mismatch_positions.swap(match.mismatches);

    return result;
}

GaplessExtension GaplessExtender::extend_seeds(cluster_type& cluster, const std::string& sequence, size_t max_mismatches, bool cluster_is_sorted) const {

    GaplessMatch best_match {
        max_mismatches + 1,
        static_cast<size_t>(0), static_cast<size_t>(0),
        static_cast<size_t>(0),
        { },
        { }
    };
    if (this->graph == nullptr) {
        return match_to_extension(best_match, sequence);
    }

    // Process the seeds in sorted order.
    seed_type prev(0, make_pos_t(0, false, 0));
    if (!cluster_is_sorted) {
        std::sort(cluster.begin(), cluster.end());
    }
    for (size_t i = 0; i < cluster.size(); i++) {

        // Start matching as early in the initial node as possible.
        seed_type hit = cluster[i];
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
                return match_to_extension(best_match, sequence);
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
                return match_to_extension(best_match, sequence);
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
                return match_to_extension(best_match, sequence);
            }
        }
    }

    return match_to_extension(best_match, sequence);
}

//------------------------------------------------------------------------------

struct UnambiguousMatch {
    size_t seq_start, seq_limit; // Sequence range.
    size_t node_start, node_limit; // In the initial/final nodes of the path.
    gbwt::BidirectionalState state;

    std::vector<handle_t> path;

    size_t length() const { return this->seq_limit - this->seq_start; }
    bool empty() const { return (this->length() == 0); }

    // Compare sequence ranges, start/end nodes, and node offsets.
    bool operator<(const UnambiguousMatch& another) const {
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
void match_forward(const std::string& seq, std::pair<const char*, size_t> target, UnambiguousMatch& match) {
    while (match.seq_limit < seq.length() && match.node_limit < target.second) {
        if (seq[match.seq_limit] != target.first[match.node_limit]) {
            break;
        }
        match.seq_limit++;
        match.node_limit++;
    }
}

// Match backward as long as the characters match.
void match_backward(const std::string& seq, std::pair<const char*, size_t> target, UnambiguousMatch& match) {
    while (match.seq_start > 0 && match.node_start > 0) {
        if (seq[match.seq_start - 1] != target.first[match.node_start - 1]) {
            break;
        }
        match.seq_start--;
        match.node_start--;
    }
}

// Convert UnambiguousMatch to GaplessExtension.
GaplessExtension unambiguous_match_to_extension(const UnambiguousMatch& match, const std::string& sequence) {

    GaplessExtension result {
        match.path, // We have to copy the path, because we get a const reference from an std::set iterator.
        match.node_start,
        match.state,
        { match.seq_start, match.seq_limit },
        (match.seq_start == 0 && match.seq_limit == sequence.length()),
        { match.seq_start, match.seq_limit },
        { }
    };

    return result;
}

std::vector<GaplessExtension> GaplessExtender::maximal_extensions(cluster_type& cluster, const std::string& sequence, bool cluster_is_sorted) const {

    // Process the seeds in sorted order.
    seed_type prev(0, make_pos_t(0, false, 0));
    size_t prev_limit = 0; // Limit in the initial node.
    if (!cluster_is_sorted) {
        std::sort(cluster.begin(), cluster.end());
    }
    std::set<UnambiguousMatch> matches;
    for (size_t i = 0; i < cluster.size(); i++) {

        // Skip redundant seeds.
        seed_type normalized = cluster[i];
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
        UnambiguousMatch match {
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

    // Convert the matches to GaplessExtension objects.
    std::vector<GaplessExtension> result;
    result.reserve(matches.size());
    for(const UnambiguousMatch& match : matches) {
        result.emplace_back(unambiguous_match_to_extension(match, sequence));
    }
    return result;
}

//------------------------------------------------------------------------------

struct FlankState {
    gbwt::BidirectionalState state;
    size_t seq_start, seq_limit; // In the sequence.
    size_t match_start, match_limit; // Positions bounding the first/last matching character in the sequence.
    size_t node_start, node_limit; // In the initial/final nodes of the path.

    std::vector<size_t> head_mismatches, tail_mismatches;

    size_t length() const { return this->seq_limit - this->seq_start; }
    size_t mismatches() const { return this->head_mismatches.size() + this->tail_mismatches.size(); }

    bool left_maximal(size_t max_mismatches) const {
        return (this->head_mismatches.size() > max_mismatches || this->seq_start == 0);
    }

    bool right_maximal(size_t max_mismatches, size_t read_length) const {
        return (this->tail_mismatches.size() > max_mismatches || this->seq_limit >= read_length);
    }

    bool at_start() const {
        return (this->node_start == 0);
    }

    bool at_end(size_t node_length) const {
        return (this->node_limit >= node_length);
    }

    void trim_head() {
        this->seq_start = this->match_start;
        size_t tail = this->head_mismatches.size();
        while (tail > 0 && this->head_mismatches[tail - 1] < this->seq_start) {
            tail--;
        }
        this->head_mismatches.resize(tail);
    }

    void trim_tail() {
        this->seq_limit = this->match_limit;
        size_t tail = this->tail_mismatches.size();
        while (tail > 0 && this->tail_mismatches[tail - 1] >= this->seq_limit) {
            tail--;
        }
        this->tail_mismatches.resize(tail);
    }

    // Low-priority elements cover a shorter range or have more mismatches over the same range.
    bool operator<(const FlankState& another) const {
        return ((this->length() < another.length()) ||
                (this->length() == another.length() && this->mismatches() > another.mismatches()));
    }
};

// Match forward.
void match_forward(const std::string& seq, std::pair<const char*, size_t> target, FlankState& match, size_t error_bound) {
    while (match.seq_limit < seq.length() && match.node_limit < target.second) {
        if (seq[match.seq_limit] == target.first[match.node_limit]) {
            match.seq_limit++;
            match.match_limit = match.seq_limit;
            match.node_limit++;
        } else {
            match.tail_mismatches.push_back(match.seq_limit);
            match.seq_limit++;
            match.node_limit++;
            if (match.tail_mismatches.size() > error_bound) {
                return;
            }
        }
    }
}

// Match backward.
void match_backward(const std::string& seq, std::pair<const char*, size_t> target, FlankState& match, size_t error_bound) {
    while (match.seq_start > 0 && match.node_start > 0) {
        match.seq_start--;
        match.node_start--;
        if (seq[match.seq_start] == target.first[match.node_start]) {
            match.match_start = match.seq_start;
        } else {
            match.head_mismatches.push_back(match.seq_start);
            if (match.head_mismatches.size() > error_bound) {
                return;
            }
        }
    }
}

void GaplessExtender::extend_flanks(std::vector<GaplessExtension>& extensions, const std::string& sequence, size_t max_mismatches) const {
    if (this->graph == nullptr) {
        return;
    }

    for (GaplessExtension& extension : extensions) {
        if (!extension.exact() || extension.empty() || extension.full()) {
            continue;
        }

        FlankState best_match {
            extension.state,
            extension.core_interval.first, extension.core_interval.second,
            extension.core_interval.first, extension.core_interval.second,
            extension.offset, extension.tail_offset(*(this->graph)),
            { }, { }
        };

        // Match the initial/final nodes of the path.
        std::stack<FlankState> forward, backward;
        {
            handle_t backward_handle = GBWTGraph::node_to_handle(gbwt::Node::reverse(best_match.state.backward.node));
            std::pair<const char*, size_t> backward_view = this->graph->get_sequence_view(backward_handle);
            match_backward(sequence, backward_view, best_match, max_mismatches);

            handle_t forward_handle = GBWTGraph::node_to_handle(best_match.state.forward.node);
            std::pair<const char*, size_t> forward_view = this->graph->get_sequence_view(forward_handle);
            match_forward(sequence, forward_view, best_match, max_mismatches);

            if (best_match.right_maximal(max_mismatches, sequence.length())) {
                best_match.trim_tail();
                if (best_match.left_maximal(max_mismatches)) {
                    best_match.trim_head();
                } else if (best_match.at_start()) {
                    backward.push(best_match);
                }
            } else if (best_match.at_end(forward_view.second)) {
                forward.push(best_match);
            }
        }

        // Forward.
        while (!forward.empty()) {
            FlankState curr = forward.top();
            forward.pop();
            bool extension = false;
            this->graph->follow_paths(curr.state, false, [&](const gbwt::BidirectionalState& next_state) -> bool {
                if (next_state.empty()) {
                    return true;
                }
                extension = true;
                handle_t handle = GBWTGraph::node_to_handle(next_state.forward.node);
                std::pair<const char*, size_t> seq_view = this->graph->get_sequence_view(handle);
                FlankState next = curr;
                next.state = next_state;
                next.node_limit = 0;
                match_forward(sequence, seq_view, next, max_mismatches);
                if (next.right_maximal(max_mismatches, sequence.length())) {
                    next.trim_tail();
                    if (next.left_maximal(max_mismatches)) {
                        next.trim_head();
                        if (best_match < next) {
                            best_match = next;
                        }
                    } else if (next.at_start()) {
                        backward.push(next);
                    }
                } else if (next.at_end(seq_view.second)){
                    forward.push(next);
                }
                return true;
            });
            if (!extension) {
                curr.trim_tail();
                if (curr.left_maximal(max_mismatches)) {
                    curr.trim_head();
                    if (best_match < curr) {
                        best_match = curr;
                    }
                } else if (curr.at_start()) {
                    backward.push(curr);
                }
            }
        }

        // Backward.
        while (!backward.empty()) {
            FlankState curr = backward.top();
            backward.pop();
            bool extension = false;
            this->graph->follow_paths(curr.state, true, [&](const gbwt::BidirectionalState& next_state) -> bool {
                if (next_state.empty()) {
                    return true;
                }
                extension = true;
                handle_t handle = GBWTGraph::node_to_handle(gbwt::Node::reverse(next_state.backward.node));
                std::pair<const char*, size_t> seq_view = this->graph->get_sequence_view(handle);
                FlankState next = curr;
                next.state = next_state;
                next.node_start = seq_view.second;
                match_backward(sequence, seq_view, next, max_mismatches);
                if (next.left_maximal(max_mismatches)) {
                    next.trim_head();
                    if (best_match < next) {
                        best_match = next;
                    }
                } else if (next.at_start()) {
                    backward.push(next);
                }
                return true;
            });
            if (!extension) {
                curr.trim_head();
                if (best_match < curr) {
                    best_match = curr;
                }
            }
        }

        // Use best_match as the flanked extension.
        extension.flanked_interval.first = best_match.seq_start;
        extension.flanked_interval.second = best_match.seq_limit;
        extension.mismatch_positions.reserve(best_match.head_mismatches.size() + best_match.tail_mismatches.size());
        extension.mismatch_positions.insert(extension.mismatch_positions.end(),
                                            best_match.head_mismatches.rbegin(),
                                            best_match.head_mismatches.rend());
        extension.mismatch_positions.insert(extension.mismatch_positions.end(),
                                            best_match.tail_mismatches.begin(),
                                            best_match.tail_mismatches.end());
    }
}

//------------------------------------------------------------------------------

} // namespace vg
