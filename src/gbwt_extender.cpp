#include "gbwt_extender.hpp"

#include <algorithm>
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

struct WFAPoint {
    int32_t  score;
    int32_t  diagonal;
    uint32_t seq_offset;
    uint32_t node_offset;

    bool operator<(const WFAPoint& another) const {
        return (this->score < another.score || (this->score == another.score && this->diagonal < another.diagonal));
    }
};

struct WFANode {
    gbwt::SearchState state;

    // Offsets in the vector of nodes.
    size_t parent;
    std::vector<size_t> children;

    // All haplotypes end here.
    bool dead_end;

    // The end position has been reached.
    // FIXME or maybe we should just take the gap to the end and report it as the best alignment so far.
    bool at_end;

    // Wavefronts sorted by (score, diagonal).
    std::vector<WFAPoint> matches;
    std::vector<WFAPoint> insertions;
    std::vector<WFAPoint> deletions;
};

void
WFAExtender::extend(std::string sequence, pos_t from, pos_t to) const {
    if (this->graph == nullptr || sequence.empty()) {
        return;
    }
    this->mask(sequence);

    std::vector<WFANode> tree;
    // For M, I, and D with each score: maintain the min and max diagonal
    // We may want to maintain a vector of possible scores.

    // while true:
    //   wf_extend()
    //   if found: break
    //   advance s
    //   wf_next()

    // Both wf_extend() and wf_next() iterate over all leaves of the tree.
    // If the (diagonal, score) does not exist in this node, try parent.
    // A node where all haplotypes end also counts as a leaf.

    // If we reach the end position, we fork. One option is to take a gap
    // until the end of the read. Another is continue on the haplotypes
    // in case there is a loop and we reach the end again.

    // FIXME stop condition?

    // Finally we have to backtrace the actual alignment.
}

//------------------------------------------------------------------------------

} // namespace vg
