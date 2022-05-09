#include "gbwt_extender.hpp"

#include <algorithm>
#include <array>
#include <cstring>
#include <queue>
#include <set>
#include <stack>

namespace vg {

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr size_t GaplessExtender::MAX_MISMATCHES;
constexpr double GaplessExtender::OVERLAP_THRESHOLD;

//------------------------------------------------------------------------------

bool GaplessExtension::contains(const HandleGraph& graph, seed_type seed) const {
    handle_t expected_handle = GaplessExtender::get_handle(seed);
    size_t expected_node_offset = GaplessExtender::get_node_offset(seed);
    size_t expected_read_offset = GaplessExtender::get_read_offset(seed);

    size_t read_offset = this->read_interval.first;
    size_t node_offset = this->offset;
    for (handle_t handle : this->path) {
        size_t len = graph.get_length(handle) - node_offset;
        read_offset += len;
        node_offset += len;
        if (handle == expected_handle && read_offset - expected_read_offset == node_offset - expected_node_offset) {
            return true;
        }
        node_offset = 0;
    }

    return false;
}

Position GaplessExtension::starting_position(const HandleGraph& graph) const {
    Position position;
    if (this->empty()) {
        return position;
    }

    position.set_node_id(graph.get_id(this->path.front()));
    position.set_is_reverse(graph.get_is_reverse(this->path.front()));
    position.set_offset(this->offset);

    return position;
}

Position GaplessExtension::tail_position(const HandleGraph& graph) const {
    Position position;
    if (this->empty()) {
        return position;
    }

    position.set_node_id(graph.get_id(this->path.back()));
    position.set_is_reverse(graph.get_is_reverse(this->path.back()));
    position.set_offset(this->tail_offset(graph));

    return position;
}

size_t GaplessExtension::tail_offset(const HandleGraph& graph) const {
    size_t result = this->offset + this->length();
    for (size_t i = 0; i + 1 < this->path.size(); i++) {
        result -= graph.get_length(this->path[i]);
    }
    return result;
}

size_t GaplessExtension::overlap(const HandleGraph& graph, const GaplessExtension& another) const {
    size_t result = 0;
    size_t this_pos = this->read_interval.first, another_pos = another.read_interval.first;
    auto this_iter = this->path.begin(), another_iter = another.path.begin();
    size_t this_offset = this->offset, another_offset = another.offset;
    while (this_pos < this->read_interval.second && another_pos < another.read_interval.second) {
        if (this_pos == another_pos && *this_iter == *another_iter && this_offset == another_offset) {
            size_t len = std::min({ graph.get_length(*this_iter) - this_offset,
                                    this->read_interval.second - this_pos,
                                    another.read_interval.second - another_pos });
            result += len;
            this_pos += len;
            another_pos += len;
            ++this_iter;
            ++another_iter;
            this_offset = 0;
            another_offset = 0;
        } else if (this_pos <= another_pos) {
            this_pos += graph.get_length(*this_iter) - this_offset;
            ++this_iter;
            this_offset = 0;
        } else {
            another_pos += graph.get_length(*another_iter) - another_offset;
            ++another_iter;
            another_offset = 0;
        }
    }
    return result;
}

Path GaplessExtension::to_path(const HandleGraph& graph, const std::string& sequence) const {

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

ReadMasker::ReadMasker(const std::string& valid_chars) : mask(256, 'X') {
    for (char c : valid_chars) {
        this->mask[static_cast<size_t>(c)] = c;
    }
}

void ReadMasker::operator()(std::string& sequence) const {
    for (char& c : sequence) {
        c = this->mask[static_cast<size_t>(c)];
    }
}

//------------------------------------------------------------------------------

GaplessExtender::GaplessExtender() :
    graph(nullptr), aligner(nullptr), mask("ACGT")
{
}

GaplessExtender::GaplessExtender(const gbwtgraph::GBWTGraph& graph, const Aligner& aligner) :
    graph(&graph), aligner(&aligner), mask("ACGT")
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

// Compute the score based on read_interval, internal_score, left_full, and right_full.
void set_score(GaplessExtension& extension, const Aligner* aligner) {
    // Assume that everything matches.
    extension.score = static_cast<int32_t>((extension.read_interval.second - extension.read_interval.first) * aligner->match);
    // Handle the mismatches.
    extension.score -= static_cast<int32_t>(extension.internal_score * (aligner->match + aligner->mismatch));
    // Handle full-length bonuses.
    extension.score += static_cast<int32_t>(extension.left_full * aligner->full_length_bonus);
    extension.score += static_cast<int32_t>(extension.right_full * aligner->full_length_bonus); 
}

// Match the initial node, assuming that read_offset or node_offset is 0.
// Updates internal_score and old_score; use set_score() to compute score.
void match_initial(GaplessExtension& match, const std::string& seq, gbwtgraph::view_type target) {
    size_t node_offset = match.offset;
    size_t left = std::min(seq.length() - match.read_interval.second, target.second - node_offset);
    while (left > 0) {
        size_t len = std::min(left, sizeof(std::uint64_t));
        std::uint64_t a = 0, b = 0;
        std::memcpy(&a, seq.data() + match.read_interval.second, len);
        std::memcpy(&b, target.first + node_offset, len);
        if (a == b) {
            match.read_interval.second += len;
            node_offset += len;
        } else {
            for (size_t i = 0; i < len; i++) {
                if (seq[match.read_interval.second] != target.first[node_offset]) {
                    match.internal_score++;
                }
                match.read_interval.second++;
                node_offset++;
            }
        }
        left -= len;
    }
    match.old_score = match.internal_score;
}

// Match forward but stop before the mismatch count reaches the limit.
// Updates internal_score; use set_score() to recompute score.
// Returns the tail offset (the number of characters matched).
size_t match_forward(GaplessExtension& match, const std::string& seq, gbwtgraph::view_type target, uint32_t mismatch_limit) {
    size_t node_offset = 0;
    size_t left = std::min(seq.length() - match.read_interval.second, target.second - node_offset);
    while (left > 0) {
        size_t len = std::min(left, sizeof(std::uint64_t));
        std::uint64_t a = 0, b = 0;
        std::memcpy(&a, seq.data() + match.read_interval.second, len);
        std::memcpy(&b, target.first + node_offset, len);
        if (a == b) {
            match.read_interval.second += len;
            node_offset += len;
        } else {
            for (size_t i = 0; i < len; i++) {
                if (seq[match.read_interval.second] != target.first[node_offset]) {
                    if (match.internal_score + 1 >= mismatch_limit) {
                        return node_offset;
                    }
                    match.internal_score++;
                }
                match.read_interval.second++;
                node_offset++;
            }
        }
        left -= len;
    }
    return node_offset;
}

// Match forward but stop before the mismatch count reaches the limit.
// Starts from the offset in the match and updates it.
// Updates internal_score; use set_score() to recompute score.
void match_backward(GaplessExtension& match, const std::string& seq, gbwtgraph::view_type target, uint32_t mismatch_limit) {
    size_t left = std::min(match.read_interval.first, match.offset);
    while (left > 0) {
        size_t len = std::min(left, sizeof(std::uint64_t));
        std::uint64_t a = 0, b = 0;
        std::memcpy(&a, seq.data() + match.read_interval.first - len, len);
        std::memcpy(&b, target.first + match.offset - len, len);
        if (a == b) {
            match.read_interval.first -= len;
            match.offset -= len;
        } else {
            for (size_t i = 0; i < len; i++) {
                if (seq[match.read_interval.first - 1] != target.first[match.offset - 1]) {
                    if (match.internal_score + 1 >= mismatch_limit) {
                        return;
                    }
                    match.internal_score++;
                }
                match.read_interval.first--;
                match.offset--;
            }
        }
        left -= len;
    }
}

// Sort full-length extensions by internal_score, remove ones that are not
// full-length alignments, remove duplicates, and return the best extensions
// that have sufficiently low overlap.
void handle_full_length(const HandleGraph& graph, std::vector<GaplessExtension>& result, double overlap_threshold) {
    std::sort(result.begin(), result.end(), [](const GaplessExtension& a, const GaplessExtension& b) -> bool {
        if (a.full() && b.full()) {
            return (a.internal_score < b.internal_score);
        }
        return a.full();
    });
    size_t tail = 0;
    for (size_t i = 0; i < result.size(); i++) {
        if (!(result[i].full())) {
            break; // No remaining full-length extensions.
        }
        bool overlap = false;
        for (size_t prev = 0; prev < tail; prev++) {
            if (result[i].overlap(graph, result[prev]) > overlap_threshold * result[prev].length()) {
                overlap = true;
                break;
            }
        }
        if (overlap) {
            continue;
        }
        if (i > tail) {
            result[tail] = std::move(result[i]);
        }
        tail++;
    }
    result.resize(tail);
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
        if (a.state.backward.range != b.state.backward.range) {
           return (a.state.backward.range < b.state.backward.range);
        }
        if (a.state.forward.range != b.state.forward.range) {
           return (a.state.forward.range < b.state.forward.range);
        }
        return (a.offset < b.offset);
    };
    std::sort(result.begin(), result.end(), sort_order);
    size_t tail = 0;
    for (size_t i = 0; i < result.size(); i++) {
        if (result[i].empty()) {
            continue;
        }
        if (tail == 0 || result[i] != result[tail - 1]) {
            if (i > tail) {
                result[tail] = std::move(result[i]);
            }
            tail++;
        }
    }
    result.resize(tail);
}

// Realign the extensions to find the mismatching positions.
void find_mismatches(const std::string& seq, const gbwtgraph::CachedGBWTGraph& graph, std::vector<GaplessExtension>& result) {
    for (GaplessExtension& extension : result) {
        if (extension.internal_score == 0) {
            continue;
        }
        extension.mismatch_positions.reserve(extension.internal_score);
        size_t node_offset = extension.offset, read_offset = extension.read_interval.first;
        for (const handle_t& handle : extension.path) {
            gbwtgraph::view_type target = graph.get_sequence_view(handle);
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

size_t interval_length(std::pair<size_t, size_t> interval) {
    return interval.second - interval.first;
}

std::vector<handle_t> get_path(const std::vector<handle_t>& first, handle_t second) {
    std::vector<handle_t> result;
    result.reserve(first.size() + 1);
    result.insert(result.end(), first.begin(), first.end());
    result.push_back(second);
    return result;
}

std::vector<handle_t> get_path(handle_t first, const std::vector<handle_t>& second) {
    std::vector<handle_t> result;
    result.reserve(second.size() + 1);
    result.push_back(first);
    result.insert(result.end(), second.begin(), second.end());
    return result;
}

std::vector<handle_t> get_path(const std::vector<handle_t>& first, gbwt::node_type second) {
    return get_path(first, gbwtgraph::GBWTGraph::node_to_handle(second));
}

std::vector<handle_t> get_path(gbwt::node_type reverse_first, const std::vector<handle_t>& second) {
    return get_path(gbwtgraph::GBWTGraph::node_to_handle(gbwt::Node::reverse(reverse_first)), second);
}

//------------------------------------------------------------------------------

// Trim mismatches from the extension to maximize the score. Returns true if the
// extension was trimmed.
bool trim_mismatches(GaplessExtension& extension, const gbwtgraph::CachedGBWTGraph& graph, const Aligner& aligner) {

    if (extension.exact()) {
        return false;
    }

    // Start with the initial run of matches.
    auto mismatch = extension.mismatch_positions.begin();
    std::pair<size_t, size_t> current_interval(extension.read_interval.first, *mismatch);
    int32_t current_score = interval_length(current_interval) * aligner.match;
    if (extension.left_full) {
        current_score += aligner.full_length_bonus;
    }

    // Process the alignment and keep track of the best interval we have seen so far.
    std::pair<size_t, size_t> best_interval = current_interval;
    int32_t best_score = current_score;
    while (mismatch != extension.mismatch_positions.end()) {
        // See if we should start a new interval after the mismatch.
        if (current_score >= aligner.mismatch) {
            current_interval.second++;
            current_score -= aligner.mismatch;
        } else {
            current_interval.first = current_interval.second = *mismatch + 1;
            current_score = 0;
        }
        ++mismatch;

        // Process the following run of matches.
        if (mismatch == extension.mismatch_positions.end()) {
            size_t length = extension.read_interval.second - current_interval.second;
            current_interval.second = extension.read_interval.second;
            current_score += length * aligner.match;
            if (extension.right_full) {
                current_score += aligner.full_length_bonus;
            }
        } else {
            size_t length = *mismatch - current_interval.second;
            current_interval.second = *mismatch;
            current_score += length * aligner.match;
        }

        // Update the best interval.
        if (current_score > best_score || (current_score > 0 && current_score == best_score && interval_length(current_interval) > interval_length(best_interval))) {
            best_interval = current_interval;
            best_score = current_score;
        }
    }

    // Special cases: no trimming or complete trimming.
    if (best_interval == extension.read_interval) {
        return false;
    }
    if (interval_length(best_interval) == 0) {
        extension.path.clear();
        extension.read_interval = best_interval;
        extension.mismatch_positions.clear();
        extension.score = 0;
        extension.left_full = extension.right_full = false;
        return true;
    }

    // Update alignment statistics.
    bool path_changed = false;
    if (best_interval.first > extension.read_interval.first) {
        extension.left_full = false;
    }
    if (best_interval.second < extension.read_interval.second) {
        extension.right_full = false;
    }
    size_t node_offset = extension.offset, read_offset = extension.read_interval.first;
    extension.read_interval = best_interval;
    extension.score = best_score;

    // Trim the path.
    size_t head = 0;
    while (head < extension.path.size()) {
        size_t node_length = graph.get_length(extension.path[head]);
        read_offset += node_length - node_offset;
        node_offset = 0;
        if (read_offset > extension.read_interval.first) {
            extension.offset = node_length - (read_offset - extension.read_interval.first);
            break;
        }
        head++;
    }
    size_t tail = head + 1;
    while (read_offset < extension.read_interval.second) {
        read_offset += graph.get_length(extension.path[tail]);
        tail++;
    }
    if (head > 0 || tail < extension.path.size()) {
        in_place_subvector(extension.path, head, tail);
        extension.state = graph.bd_find(extension.path);
    }

    // Trim the mismatches.
    head = 0;
    while (head < extension.mismatch_positions.size() && extension.mismatch_positions[head] < extension.read_interval.first) {
        head++;
    }
    tail = head;
    while (tail < extension.mismatch_positions.size() && extension.mismatch_positions[tail] < extension.read_interval.second) {
        tail++;
    }
    in_place_subvector(extension.mismatch_positions, head, tail);

    return true;
}

//------------------------------------------------------------------------------

std::vector<GaplessExtension> GaplessExtender::extend(cluster_type& cluster, std::string sequence, const gbwtgraph::CachedGBWTGraph* cache, size_t max_mismatches, double overlap_threshold) const {

    std::vector<GaplessExtension> result;
    if (this->graph == nullptr || this->aligner == nullptr || cluster.empty() || sequence.empty()) {
        return result;
    }
    result.reserve(cluster.size());
    this->mask(sequence);

    // Allocate a cache if we were not provided with one.
    bool free_cache = (cache == nullptr);
    if (free_cache) {
        cache = new gbwtgraph::CachedGBWTGraph(*(this->graph));
    }

    // Find the best extension starting from each seed.
    size_t best_alignment = std::numeric_limits<size_t>::max();
    for (seed_type seed : cluster) {

        // Check if the seed is contained in an exact full-length alignment.
        if (best_alignment < result.size() && result[best_alignment].internal_score == 0) {
            if (result[best_alignment].contains(*cache, seed)) {
                continue;
            }
        }

        GaplessExtension best_match {
            { }, static_cast<size_t>(0), gbwt::BidirectionalState(),
            { static_cast<size_t>(0), static_cast<size_t>(0) }, { },
            std::numeric_limits<int32_t>::min(), false, false,
            false, false, std::numeric_limits<uint32_t>::max(), std::numeric_limits<uint32_t>::max()
        };

        // Match the initial node and add it to the queue.
        std::priority_queue<GaplessExtension> extensions;
        {
            size_t read_offset = get_read_offset(seed);
            size_t node_offset = get_node_offset(seed);
            GaplessExtension match {
                { seed.first }, node_offset, cache->get_bd_state(seed.first),
                { read_offset, read_offset }, { },
                static_cast<int32_t>(0), false, false,
                false, false, static_cast<uint32_t>(0), static_cast<uint32_t>(0)
            };
            match_initial(match, sequence, cache->get_sequence_view(seed.first));
            if (match.read_interval.first == 0) {
                match.left_full = true;
                match.left_maximal = true;
            }
            if (match.read_interval.second >= sequence.length()) {
                match.right_full = true;
                match.right_maximal = true;
            }
            set_score(match, this->aligner);
            extensions.push(std::move(match));
        }

        // Extend the most promising extensions first, using alignment scores for priority.
        // First make the extension right-maximal and then left-maximal.
        while (!extensions.empty()) {
            GaplessExtension curr = std::move(extensions.top());
            extensions.pop();

            // Case 1: Extend to the right.
            if (!curr.right_maximal) {
                size_t num_extensions = 0;
                // Always allow at least max_mismatches / 2 mismatches in the current flank.
                uint32_t mismatch_limit = std::max(
                    static_cast<uint32_t>(max_mismatches + 1),
                    static_cast<uint32_t>(max_mismatches / 2 + curr.old_score + 1));
                cache->follow_paths(curr.state, false, [&](const gbwt::BidirectionalState& next_state) -> bool {
                    handle_t handle = gbwtgraph::GBWTGraph::node_to_handle(next_state.forward.node);
                    GaplessExtension next {
                        { }, curr.offset, next_state,
                        curr.read_interval, { },
                        curr.score, curr.left_full, curr.right_full,
                        curr.left_maximal, curr.right_maximal, curr.internal_score, curr.old_score
                    };
                    size_t node_offset = match_forward(next, sequence, cache->get_sequence_view(handle), mismatch_limit);
                    if (node_offset == 0) { // Did not match anything.
                        return true;
                    }
                    next.path = get_path(curr.path, handle);
                    // Did the extension become right-maximal?
                    if (next.read_interval.second >= sequence.length()) {
                        next.right_full = true;
                        next.right_maximal = true;
                        next.old_score = next.internal_score;
                    } else if (node_offset < cache->get_length(handle)) {
                        next.right_maximal = true;
                        next.old_score = next.internal_score;
                    }
                    set_score(next, this->aligner);
                    num_extensions += next.state.size();
                    extensions.push(std::move(next));
                    return true;
                });
                // We could not extend all threads in 'curr' to the right. The unextended ones
                // may have different left extensions, so we must consider 'curr' right-maximal.
                if (num_extensions < curr.state.size()) {
                    curr.right_maximal = true;
                    curr.old_score = curr.internal_score;
                    extensions.push(std::move(curr));
                }
                continue;
            }

            // Case 2: Extend to the left.
            if (!curr.left_maximal) {
                bool found_extension = false;
                // Always allow at least max_mismatches / 2 mismatches in the current flank.
                uint32_t mismatch_limit = std::max(
                    static_cast<uint32_t>(max_mismatches + 1),
                    static_cast<uint32_t>(max_mismatches / 2 + curr.old_score + 1));
                cache->follow_paths(curr.state, true, [&](const gbwt::BidirectionalState& next_state) -> bool {
                    handle_t handle = gbwtgraph::GBWTGraph::node_to_handle(gbwt::Node::reverse(next_state.backward.node));
                    size_t node_length = cache->get_length(handle);
                    GaplessExtension next {
                        { }, node_length, next_state,
                        curr.read_interval, { },
                        curr.score, curr.left_full, curr.right_full,
                        curr.left_maximal, curr.right_maximal, curr.internal_score, curr.old_score
                    };
                    match_backward(next, sequence, cache->get_sequence_view(handle), mismatch_limit);
                    if (next.offset >= node_length) { // Did not match anything.
                        return true;
                    }
                    next.path = get_path(handle, curr.path);
                    // Did the extension become left-maximal?
                    if (next.read_interval.first == 0) {
                        next.left_full = true;
                        next.left_maximal = true;
                        // No need to set old_score.
                    } else if (next.offset > 0) {
                        next.left_maximal = true;
                        // No need to set old_score.
                    }
                    set_score(next, this->aligner);
                    extensions.push(std::move(next));
                    found_extension = true;
                    return true;
                });
                if (!found_extension) {
                    curr.left_maximal = true;
                    // No need to set old_score.
                } else {
                    continue;
                }
            }

            // Case 3: Maximal extension with a better score than the best extension so far.
            if (best_match < curr) {
                best_match = std::move(curr);
            }
        }

        // Add the best match to the result and update the best_alignment offset.
        if (!best_match.empty()) {
            if (best_match.full() && (best_alignment >= result.size() || best_match.internal_score < result[best_alignment].internal_score)) {
                best_alignment = result.size();
            }
            result.emplace_back(std::move(best_match));
        }
    }

    // If we have a good enough full-length alignment, return the best sufficiently
    // distinct full-length alignments.
    if (best_alignment < result.size() && result[best_alignment].internal_score <= max_mismatches) {
        handle_full_length(*cache, result, overlap_threshold);
        find_mismatches(sequence, *cache, result);
    }

    // Otherwise remove duplicates, find mismatches, and trim the extensions to maximize
    // score.
    else {
        remove_duplicates(result);
        find_mismatches(sequence, *cache, result);
        bool trimmed = false;
        for (GaplessExtension& extension : result) {
            trimmed |= trim_mismatches(extension, *cache, *(this->aligner));
        }
        if (trimmed) {
            remove_duplicates(result);
        }
    }

    // Free the cache if we allocated it.
    if (free_cache) {
        delete cache;
        cache = nullptr;
    }

    return result;
}

//------------------------------------------------------------------------------

bool GaplessExtender::full_length_extensions(const std::vector<GaplessExtension>& result, size_t max_mismatches) {
    return (result.size() > 0 && result.front().full() && result.front().mismatches() <= max_mismatches);
}

//------------------------------------------------------------------------------

struct state_hash {
    size_t operator()(const gbwt::BidirectionalState& state) const {
        size_t result = wang_hash_64(state.forward.node);
        result ^= wang_hash_64(state.forward.range.first) + 0x9e3779b9 + (result << 6) + (result >> 2);
        result ^= wang_hash_64(state.forward.range.second) + 0x9e3779b9 + (result << 6) + (result >> 2);
        result ^= wang_hash_64(state.backward.node) + 0x9e3779b9 + (result << 6) + (result >> 2);
        result ^= wang_hash_64(state.backward.range.first) + 0x9e3779b9 + (result << 6) + (result >> 2);
        result ^= wang_hash_64(state.backward.range.second) + 0x9e3779b9 + (result << 6) + (result >> 2);
        return result;
    }
};

//------------------------------------------------------------------------------

WFAExtender::WFAExtender() :
    graph(nullptr), mask("ACGT"),
    mismatch(1), gap_open(1), gap_extend(1)
{
}

// Scoring parameters are from:
//   Eizenga, Paten:
//   Improving the time and space complexity of the WFA algorithm and generalizing its scoring.
//   bioRxiv, 2022.
WFAExtender::WFAExtender(const gbwtgraph::GBWTGraph& graph, const Aligner& aligner) :
    graph(&graph), mask("ACGT"),
    mismatch(2 * (aligner.match + aligner.mismatch)),
    gap_open(2 * aligner.gap_open),
    gap_extend(2 * aligner.gap_extension + aligner.match)
{
}

//------------------------------------------------------------------------------

// A position in an alignment between a sequence and a graph.
struct MatchPos {
    uint32_t seq_offset;
    uint32_t node_offset;
    std::stack<uint32_t> path; // Sequence of tree offsets from a leaf to the relevant node.

    // Creates an empty position.
    MatchPos() : seq_offset(0), node_offset(0) {}

    // Creates a position with the given offsets and path.
    MatchPos(uint32_t seq_offset, uint32_t node_offset, const std::stack<uint32_t>& path) : seq_offset(seq_offset), node_offset(node_offset), path(path) {}

    bool empty() const { return this->path.empty(); }
    bool at_last_node() const { return (this->path.size() == 1); }
    uint32_t node() const { return this->path.top(); }
    void pop() { this->path.pop(); }

    // Returns the non-empty position further on the diagonal (on the sequence).
    static const MatchPos& max(const MatchPos& a, const MatchPos& b) {
        if (a.empty()) { return b; }
        if (b.empty()) { return a; }
        return (a.seq_offset > b.seq_offset ? a : b);
    }
};

// A point in an WFA score matrix (for a specific node).
struct WFAPoint {
    int32_t  score;
    int32_t  diagonal; // seq_offset - target offset
    uint32_t seq_offset;
    uint32_t node_offset;

    constexpr static std::pair<uint32_t, uint32_t> NO_OFFSET {
        std::numeric_limits<uint32_t>::max(), std::numeric_limits<uint32_t>::max()
    };

    // Converts the point to an alignment position with the given path.
    MatchPos pos(const std::stack<uint32_t>& path) const {
        return MatchPos(this->seq_offset, this->node_offset, path);
    }

    // For ordering the points in WFANode.
    bool operator<(const WFAPoint& another) const {
        return (this->score < another.score || (this->score == another.score && this->diagonal < another.diagonal));
    }
};

//------------------------------------------------------------------------------

struct WFANode {
    gbwt::SearchState state;

    // Offsets in the vector of nodes.
    uint32_t parent;
    std::vector<uint32_t> children;

    // All haplotypes end here.
    bool dead_end;

    constexpr static size_t MATCHES = 0;
    constexpr static size_t INSERTIONS = 1; // characters in the sequence but not in the graph
    constexpr static size_t DELETIONS = 2;  // characters in the graph but not in the sequence

    // Points on the wavefronts are sorted by score, diagonal.
    std::array<std::vector<WFAPoint>, 3> wavefronts;

    WFANode(const gbwt::SearchState& state, uint32_t parent) :
        state(state),
        parent(parent), children(),
        dead_end(false),
        wavefronts() {
    }

    bool is_leaf() const { return (this->children.empty() || this->dead_end); }
    bool expanded() const { return (!this->children.empty() || this->dead_end); }

    bool same_node(pos_t pos) const {
        return (gbwt::Node::id(this->state.node) == id(pos) && gbwt::Node::is_reverse(this->state.node) == is_rev(pos));
    }

    size_t length(const gbwtgraph::GBWTGraph& graph) const {
        return graph.get_length(gbwtgraph::GBWTGraph::node_to_handle(this->state.node));
    }

    // Returns the position for the given score and diagonal with the given path, or an empty position if it does not exist.
    MatchPos find_pos(size_t type, int32_t score, int32_t diagonal, std::stack<uint32_t>& path) const {
        WFAPoint key { score, diagonal, 0, 0 };
        const std::vector<WFAPoint>& points = this->wavefronts[type];
        auto iter = std::lower_bound(points.begin(), points.end(), key);
        if (iter == points.end() || key < *iter) {
            return MatchPos();
        }
        return iter->pos(path);
    }

    // Update the WFA matrix with the given alignment position.
    void update(size_t type, int32_t score, int32_t diagonal, const MatchPos& pos) {
        this->update(type, score, diagonal, pos.seq_offset, pos.node_offset);
    }

    // Update the WFA matrix with the given offsets.
    void update(size_t type, int32_t score, int32_t diagonal, uint32_t seq_offset, uint32_t node_offset) {
        WFAPoint key { score, diagonal, seq_offset, node_offset };
        std::vector<WFAPoint>& points = this->wavefronts[type];
        auto iter = std::lower_bound(points.begin(), points.end(), key);
        if (iter == points.end() || key < *iter) {
            points.insert(iter, key);
        } else {
            *iter = key;
        }
    }

    // Returns a position at the first non-match after the given position.
    void match_forward(const std::string& sequence, const gbwtgraph::GBWTGraph& graph, MatchPos& pos) const {
        handle_t handle = gbwtgraph::GBWTGraph::node_to_handle(this->state.node);
        gbwtgraph::view_type node_seq = graph.get_sequence_view(handle);
        while (pos.seq_offset < sequence.length() && pos.node_offset < node_seq.second && sequence[pos.seq_offset] == node_seq.first[pos.node_offset]) {
            pos.seq_offset++; pos.node_offset++;
        }
    }

    // Returns a position at the start of the run of matches before the given position.
    void match_backward(const std::string& sequence, const gbwtgraph::GBWTGraph& graph, MatchPos& pos) const {
        handle_t handle = gbwtgraph::GBWTGraph::node_to_handle(this->state.node);
        gbwtgraph::view_type node_seq = graph.get_sequence_view(handle);
        while (pos.seq_offset > 0 && pos.node_offset > 0 && sequence[pos.seq_offset - 1] == node_seq.first[pos.node_offset - 1]) {
            pos.seq_offset--; pos.node_offset--;
        }
    }
};

//------------------------------------------------------------------------------

class WFATree {
public:
    const gbwtgraph::GBWTGraph& graph;
    const std::string& sequence;

    std::vector<WFANode> nodes;

    // Best alignment found so far. If we reached the destination in the graph,
    // the score includes the possible insertion at the end but the point itself
    // does not.
    WFAPoint candidate_point;
    uint32_t candidate_node;

    // Scores copied from the extender.
    int32_t mismatch, gap_open, gap_extend;

    // The closed range of diagonals reached with each score.
    std::vector<std::pair<int32_t, int32_t>> diagonals;

    // The overall closed range of diagonals reached.
    std::pair<int32_t, int32_t> max_diagonals;

    WFATree(const gbwtgraph::GBWTGraph& graph, const std::string& sequence, const gbwt::SearchState& root, int32_t mismatch, int32_t gap_open, int32_t gap_extend) :
        graph(graph), sequence(sequence),
        nodes(),
        candidate_point({ std::numeric_limits<int32_t>::max(), 0, 0, 0 }), candidate_node(0),
        mismatch(mismatch), gap_open(gap_open), gap_extend(gap_extend),
        diagonals({ { 0, 0 } }), max_diagonals(0, 0)
    {
        this->nodes.emplace_back(root, 0);
    }

    uint32_t size() const { return this->nodes.size(); }
    static bool is_root(uint32_t node) { return (node == 0); }
    uint32_t parent(uint32_t node) const { return this->nodes[node].parent; }

    // wf_extend() in the paper.
    // If we reach the end of a node, we continue to the start of the next node even
    // if we do not use any characters in it.
    void extend(int32_t score, pos_t to) {
        for (int32_t diagonal = this->max_diagonals.first; diagonal <= this->max_diagonals.second; diagonal++) {
            std::vector<uint32_t> leaves = this->get_leaves();
            this->extend_over(score, diagonal, to, leaves);
        }
    }

    // wf_next() in the paper.
    // If we reach the end of a node, we continue to the start of the next node even
    // if we do not use any characters in it.
    void next(int32_t score, pos_t to) {
        std::pair<int32_t, int32_t> diagonal_range = this->get_diagonals(score);
        for (int32_t diagonal = diagonal_range.first; diagonal <= diagonal_range.second; diagonal++) {
            std::vector<uint32_t> leaves = this->get_leaves();
            // Note that we may do the same update from multiple leaves.
            for (uint32_t leaf : leaves) {
                MatchPos ins_open = this->find_pos(WFANode::MATCHES, leaf, score - this->gap_open - this->gap_extend, diagonal - 1, true, false);
                MatchPos ins_extend = this->find_pos(WFANode::INSERTIONS, leaf, score - this->gap_extend, diagonal - 1, true, false);
                MatchPos ins = MatchPos::max(ins_open, ins_extend);
                if (!ins.empty()) {
                    ins.seq_offset++;
                    this->nodes[ins.node()].update(WFANode::INSERTIONS, score, diagonal, ins);
                }

                MatchPos del_open = this->find_pos(WFANode::MATCHES, leaf, score - this->gap_open - this->gap_extend, diagonal + 1, false, true);
                MatchPos del_extend = this->find_pos(WFANode::DELETIONS, leaf, score - this->gap_extend, diagonal + 1, false, true);
                MatchPos del = MatchPos::max(del_open, del_extend);
                if (!del.empty()) {
                    del.node_offset++;
                    this->nodes[del.node()].update(WFANode::DELETIONS, score, diagonal, del);
                    this->propagate(del, score, diagonal, WFANode::DELETIONS);
                }

                MatchPos subst = this->find_pos(WFANode::MATCHES, leaf, score - this->mismatch, diagonal, true, true);
                if (!subst.empty()) {
                    subst.seq_offset++;
                    subst.node_offset++;
                }
                subst = MatchPos::max(subst, ins);
                subst = MatchPos::max(subst, del);
                if (!subst.empty()) {
                    // If we are just past the end position, we get a candidate alignment by
                    // assuming that the rest of the sequence is an insertion.
                    // If `subst` is an insertion, we already got a better alignment from the
                    // match preceding it.
                    if (this->nodes[subst.node()].same_node(to) && subst.node_offset == offset(to) + 1) {
                        uint32_t gap_length = this->sequence.length() - subst.seq_offset;
                        int32_t new_score = score;
                        if (gap_length > 0) {
                            new_score += this->gap_open + gap_length * this->gap_extend;
                        }
                        if (new_score < this->candidate_point.score) {
                            this->candidate_point = { new_score, diagonal, subst.seq_offset, subst.node_offset };
                            this->candidate_node = subst.node();
                        }
                    }
                    this->nodes[subst.node()].update(WFANode::MATCHES, score, diagonal, subst);
                    this->propagate(subst, score, diagonal, WFANode::MATCHES);
                }
            }
        }
    }

private:

    // wf_extend() on a specific diagonal for the set of (local) haplotypes corresponding to
    // the given list of leaves in the tree of GBWT search states.
    void extend_over(int32_t score, int32_t diagonal, pos_t to, const std::vector<uint32_t>& leaves) {
        for (uint32_t leaf : leaves) {
            MatchPos pos = this->find_pos(WFANode::MATCHES, leaf, score, diagonal, false, false);
            if (pos.empty()) {
                continue; // This should never happen.
            }
            while (true) {
                bool may_reach_to = this->nodes[pos.node()].same_node(to) & (pos.node_offset <= offset(to));
                this->nodes[pos.node()].match_forward(this->sequence, this->graph, pos);
                // We got a match that covers the end position but may extend past it.
                // This gives us a candidate where the rest of the sequence is an insertion.
                if (may_reach_to && pos.node_offset > offset(to)) {
                    uint32_t gap_length = (this->sequence.length() - pos.seq_offset) + (pos.node_offset - offset(to) - 1);
                    int32_t new_score = score;
                    if (gap_length > 0) {
                        score += this->gap_open + gap_length * this->gap_extend;
                    }
                    if (new_score < this->candidate_point.score) {
                        this->candidate_point = { new_score, diagonal, pos.seq_offset, pos.node_offset };
                        this->candidate_node = pos.node();
                    }
                }
                this->nodes[pos.node()].update(WFANode::MATCHES, score, diagonal, pos);
                if (pos.node_offset < this->nodes[pos.node()].length(this->graph)) {
                    break;
                }
                if (pos.at_last_node()) {
                    if (this->propagate(pos, score, diagonal, WFANode::MATCHES)) {
                        this->extend_over(score, diagonal, to, this->nodes[leaf].children);
                    }
                    break;
                }
                pos.pop();
                pos.node_offset = 0;
            }
        }
    }

    std::vector<uint32_t> get_leaves() const {
        std::vector<uint32_t> leaves;
        for (uint32_t node = 0; node < this->size(); node++) {
            if (this->nodes[node].is_leaf()) {
                leaves.push_back(node);
            }
        }
        return leaves;
    }

    std::pair<int32_t, int32_t> update_range(std::pair<int32_t, int32_t> range, int32_t score) const {
        if (score > 0) {
            range.first = std::min(range.first, this->diagonals[score].first);
            range.second = std::max(range.second, this->diagonals[score].second);
        }
        return range;
    }

    // Determine the diagonal range for the given score and update diagonals, max_diagonals.
    std::pair<int32_t, int32_t> get_diagonals(int32_t score) {
        // Determine the diagonal range for the given score.
        std::pair<int32_t, int32_t> range(0, 0);
        range = this->update_range(range, score - this->mismatch); // Mismatch.
        range = this->update_range(range, score - this->gap_open - this->gap_extend); // New gap.
        range = this->update_range(range, score - this->gap_extend); // Extend an existing gap.
        range.first--; range.second++;

        this->max_diagonals.first = std::min(this->max_diagonals.first, range.first);
        this->max_diagonals.second = std::max(this->max_diagonals.second, range.second);
        this->diagonals.resize(score + 1, std::pair<int32_t, int32_t>(0, 0));
        this->diagonals.back() = range;

        return range;
    }

    // If we have reached the end of the path corresponding to the position, expand
    // the children of the final node if they have not been expanded yet and then
    // propagate the position to the WFA matrixes of the given type in the children.
    // Returns true if we actually propagated the position to the children.
    bool propagate(const MatchPos& pos, int32_t score, int32_t diagonal, size_t type) {
        if (!pos.at_last_node() || pos.node_offset < this->nodes[pos.node()].length(this->graph)) {
            return false;
        }
        if (!this->nodes[pos.node()].expanded()) {
            bool found = false;
            this->graph.follow_paths(this->nodes[pos.node()].state, [&](const gbwt::SearchState& child) -> bool {
                this->nodes[pos.node()].children.push_back(this->size());
                this->nodes.emplace_back(child, pos.node());
                found = true;
                return true;
            });
            if (!found) {
                this->nodes[pos.node()].dead_end = true;
            }
        }
        for (uint32_t child : this->nodes[pos.node()].children) {
            this->nodes[child].update(type, score, diagonal, pos.seq_offset, 0);
        }
        return !this->nodes[pos.node()].children.empty();
    }

    // Returns the furthest position in given WFA matrix for (score, diagonal) at the
    // specified node or its ancestors, or an empty position if it does not exist.
    // Returns an empty position if an extendable position is requested but the position
    // cannot be extended.
    MatchPos find_pos(size_t type, uint32_t node, int32_t score, int32_t diagonal, bool extendable_seq, bool extendable_graph) const {
        if (score < 0) {
            return MatchPos();
        }
        std::stack<uint32_t> path;
        while(true) {
            path.push(node);
            MatchPos pos = this->nodes[node].find_pos(type, score, diagonal, path);
            if (!pos.empty()) {
                if (extendable_seq && pos.seq_offset >= this->sequence.length()) {
                    return MatchPos();
                }
                if (extendable_graph && this->at_dead_end(pos)) {
                    return MatchPos();
                }
                return pos;
            }
            if (is_root(node)) {
                return MatchPos();
            }
            node = this->parent(node);
        }
    }

    // Assumes that the position is non-empty.
    bool at_dead_end(const MatchPos& pos) const {
        return (this->nodes[pos.node()].dead_end && pos.node_offset >= this->nodes[pos.node()].length(this->graph));
    }
};

//------------------------------------------------------------------------------

void WFAExtender::connect(std::string sequence, pos_t from, pos_t to) const {
    if (this->graph == nullptr || sequence.empty()) {
        return;
    }
    gbwt::SearchState root_state = this->graph->get_state(this->graph->get_handle(id(from), is_rev(from)));
    if (root_state.empty()) {
        return;
    }
    this->mask(sequence);

    WFATree tree(*(this->graph), sequence, root_state, this->mismatch, this->gap_open, this->gap_extend);

    int32_t score = 0;
    // FIXME stop condition?
    while (true) {
        tree.extend(score, to);
        if (tree.candidate_point.score <= score) {
            break;
        }
        score++; // FIXME skip impossible scores
        tree.next(score, to);
    }

    // FIXME Finally we have to backtrace the actual alignment.
    // FIXME note that we may have edits in a node where we do not match anything
}

void WFAExtender::suffix(const std::string& sequence, pos_t from) const {
    this->connect(sequence, from, pos_t(0, false, 0));
}

void WFAExtender::prefix(const std::string& sequence, pos_t to) const {
    if (this->graph == nullptr) {
        return;
    }

    // Flip the position, extend forward, and reverse the return value.
    to = reverse_base_pos(to, this->graph->get_length(this->graph->get_handle(id(to), is_rev(to))));
    this->connect(reverse_complement(sequence), to, pos_t(0, false, 0));

    // FIXME reverse the return value
}

//------------------------------------------------------------------------------

} // namespace vg
