#include "gbwt_extender.hpp"
#include "hash_map.hpp"

#include <algorithm>
#include <array>
#include <cstring>
#include <queue>
#include <set>

#include <structures/immutable_list.hpp>

namespace vg {

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr size_t GaplessExtender::MAX_MISMATCHES;
constexpr double GaplessExtender::OVERLAP_THRESHOLD;

//------------------------------------------------------------------------------

bool GaplessExtension::for_each_read_interval(const HandleGraph& graph, const std::function<bool(size_t, size_t, const seed_type&)>& iteratee) const {
    // Track correspondign read and node offsets on the current node
    size_t read_offset = this->read_interval.first;
    size_t node_offset = this->offset;
    for (const handle_t& handle : this->path) {
        // For each node
        
        // How many bases of the node do we use? Either remaining node or remaining read if shorter.
        size_t len = std::min(graph.get_length(handle) - node_offset, this->read_interval.second - read_offset);
        if (!iteratee(read_offset, len, seed_type(handle, read_offset - node_offset))) {
            return false;
        }
        read_offset += len;
        node_offset = 0;
    }
    return true;
}

bool GaplessExtension::contains(const HandleGraph& graph, const seed_type& seed) const {
    // Scan all the seeds we represent to see if that one is one of them.
    bool found = false;
    for_each_read_interval(graph, [&](size_t read_offset, size_t len, const seed_type& our_seed) {
        if (our_seed == seed) {
            found = true;
            return false;
        }
        return true;
    });

    return found;
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

std::vector<GaplessExtension> GaplessExtender::extend(cluster_type& cluster, std::string sequence, const gbwtgraph::CachedGBWTGraph* cache, size_t max_mismatches, double overlap_threshold, bool trim) const {

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
        std::priority_queue<std::pair<GaplessExtension, size_t>> extensions;
        // Each extension gets a unique deduplicating number to break score
        // ties and make the queue order stable across different backing STL
        // implementations.
        size_t extension_number = 0;

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
            extensions.emplace(std::move(match), extension_number++);
        }

        // Extend the most promising extensions first, using alignment scores for priority.
        // First make the extension right-maximal and then left-maximal.
        while (!extensions.empty()) {
            GaplessExtension curr = std::move(extensions.top().first);
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
                    extensions.emplace(std::move(next), extension_number++);
                    return true;
                });
                // We could not extend all threads in 'curr' to the right. The unextended ones
                // may have different left extensions, so we must consider 'curr' right-maximal.
                if (num_extensions < curr.state.size()) {
                    curr.right_maximal = true;
                    curr.old_score = curr.internal_score;
                    extensions.emplace(std::move(curr), extension_number++);
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
                    extensions.emplace(std::move(next), extension_number++);
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
        if (trim) {
            // It's OK if out extensions don't include all matches between the
            // read and each node that are in phase with our seeds. Trim back
            // to maximize score.
            bool trimmed = false;
            for (GaplessExtension& extension : result) {
                trimmed |= trim_mismatches(extension, *cache, *(this->aligner));
            }
            if (trimmed) {
                remove_duplicates(result);
    
            }
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

WFAAlignment WFAAlignment::from_extension(const GaplessExtension& extension) {

    // Start by aggregate-initializing.
    WFAAlignment to_return {
        extension.path,
        {},
        (uint32_t)extension.offset,
        (uint32_t)extension.read_interval.first,
        (uint32_t)extension.length(),
        extension.score,
        true
    };

    // We need to make edits for all the mismatches.
    // This tracks the base after the last edit in edits, in the sequence space.
    size_t edits_made_up_to = to_return.seq_offset;
    for (auto& mismatch_at : extension.mismatch_positions) {
        // For each mismatch position
        if (!to_return.edits.empty() && edits_made_up_to == mismatch_at && to_return.edits.back().first == mismatch) {
            // If we can glom it onto an existing mismatch, do that.
            ++to_return.edits.back().second;
        } else {
            // Otherwise, we need some new edits
            if (edits_made_up_to < mismatch_at) {
                // Add a match for the intervening non-mismatch sequence
                to_return.edits.emplace_back(match, mismatch_at - edits_made_up_to);
            }
            // Add a new 1 base mismatch
            to_return.edits.emplace_back(mismatch, 1);
        }
        // Advance the cursor to through this mismatch
        edits_made_up_to = mismatch_at;
    }
    if (edits_made_up_to < to_return.seq_offset + to_return.length) {
        // Add any trailing match
        to_return.edits.emplace_back(match, (to_return.seq_offset + to_return.length) - edits_made_up_to);
    }

    return to_return;
}

WFAAlignment WFAAlignment::make_unlocalized_insertion(size_t sequence_offset, size_t length, int score) {
    // We can do it all by aggregate-initializing
    return {{}, {{insertion, length}}, 0, (uint32_t)sequence_offset, (uint32_t)length, score, true};
}

WFAAlignment WFAAlignment::make_empty() {
    // We can do it all by aggregate-initializing
    return {{}, {}, 0, 0, 0, 0, true};
}

bool WFAAlignment::unlocalized_insertion() const {
    return (
        this->ok &&
        this->path.empty() &&
        this->edits.size() == 1 &&
        this->edits.front().first == insertion
    );
}

int64_t WFAAlignment::final_offset(const gbwtgraph::GBWTGraph& graph) const {
    int64_t final_offset = this->node_offset;
    for (auto edit : this->edits) {
        if (edit.first != WFAAlignment::insertion) {
            final_offset += edit.second;
        }
    }
    for (size_t i = 0; i + 1 < this->path.size(); i++) {
        final_offset -= graph.get_length(this->path[i]);
    }
    return final_offset;
}

void WFAAlignment::flip(const gbwtgraph::GBWTGraph& graph, const std::string& sequence) {
    this->seq_offset = sequence.length() - this->seq_offset - this->length;

    if (this->path.empty()) {
        return;
    }
    this->node_offset = graph.get_length(this->path.back()) - this->final_offset(graph);

    // Reverse the path and the edits.
    std::reverse(this->path.begin(), this->path.end());
    for (size_t i = 0; i < this->path.size(); i++) {
        this->path[i] = graph.flip(this->path[i]);
    }
    std::reverse(this->edits.begin(), this->edits.end());
}

void WFAAlignment::append(Edit edit, uint32_t length) {
    if (length == 0) {
        return;
    }
    if (this->edits.empty() || this->edits.back().first != edit) {
        this->edits.push_back(std::make_pair(edit, length));
    } else {
        this->edits.back().second += length;
    }
}

//#define debug_join

void WFAAlignment::join(const WFAAlignment& second) {
#ifdef debug_join
    std::cerr << "Joining alignment of sequence " << seq_offset << " - " << (seq_offset + length)
        << " with alignment of " << second.seq_offset << " - " << (second.seq_offset + second.length) << std::endl;
    std::cerr << "Left alignment: ";
    print(std::cerr);
    std::cerr << std::endl;
    std::cerr << "Right alignment: ";
    second.print(std::cerr);
    std::cerr << std::endl;
#endif

    if (!ok) {
        throw std::runtime_error("Cannot join onto an alignment that is not OK");
    }

    if (!second.ok) {
        throw std::runtime_error("Cannot join an alignment that is not OK onto another alignment");
    }

    if (second.empty()) {
        // We are joining an empty alignment onto us. Do nothing.
        return;
    }

    if (empty()) {
        // We are ourselves empty. Just be replaced.
        *this = second;
        return;
    }

    // Otherwise there is actual splicing to do.

    // Do error checking
    if (seq_offset + length != second.seq_offset) {
        throw std::runtime_error("Cannot join alignments because past-end position " +
                                 std::to_string(seq_offset + length) +
                                 " is not at start position " +
                                 std::to_string(second.seq_offset));
    }
    if (path.empty() && ! unlocalized_insertion()) {
        throw std::runtime_error("Cannot join alignments because first alignment has no path");
    }
    if (second.path.empty() && ! second.unlocalized_insertion()) {
        throw std::runtime_error("Cannot join alignments because second alignment has no path");
    }
    if (edits.empty()) {
        throw std::runtime_error("Cannot join alignments because first alignment has no edits");
    }
    if (second.edits.empty()) {
        throw std::runtime_error("Cannot join alignments because second alignment has no edits");
    }
    
    if (!second.unlocalized_insertion()) {
        // The second alignment has a path
        if (unlocalized_insertion()) {
            // We don't, so include the first handle and copy the offset on it.
            node_offset = second.node_offset;
            path.push_back(second.path.front());
        } else if (second.node_offset == 0) {
            // Include the first handle from the second alignment because it can't be shared
            path.push_back(second.path.front());
        } else {
            // It must be shared with this alignment
            if (second.path.front() != path.back()) {
                throw std::runtime_error("Cannot join alignments because second alignment starts in the middle of a handle that first alignment doesn't end on");
            }
        }

        // Copy all the other path handles.
        std::copy(second.path.begin() + 1, second.path.end(), std::back_inserter(path));
    }

    for (auto& edit : second.edits) {
        // Copy over all the edits
        append(edit.first, edit.second);
    }

    // Offsets don't need to change.

    // Add the length
    length += second.length;

    // Add the score
    score += second.score;

    // And that's all the fields we have!
}

//#define debug_path

Path WFAAlignment::to_path(const HandleGraph& graph, const std::string& sequence) const {

    if (!*this) {
        throw std::runtime_error("WFAAlignment is not OK and cannot become a path");
    }

    if (this->seq_offset + this->length > sequence.size()) {
        throw std::runtime_error("WFAAlignment extends past end of sequence");
    }

    Path result;
    
    if (this->unlocalized_insertion()) {
        // Special handling for pure insertions with no position: make a Path
        // that has no position.
        vg::Mapping* m = result.add_mapping();
        vg::Edit* e = m->add_edit();
        e->set_to_length(this->edits.front().second);
        e->set_sequence(sequence.substr(this->seq_offset, this->edits.front().second));
        return result;
    }

    if (this->path.empty()) {
        // Not a pure insert, but has an empty path. We assume it's an empty WFAAlignment.
        // Nothing to do!
        return result;
    }

    // Walk through the sequence
    size_t sequence_cursor = this->seq_offset;

    // And each node along the path
    auto path_it = this->path.begin();

    // And each base along the current node
    size_t node_cursor = this->node_offset;
    size_t first_node_length = graph.get_length(*path_it);
    if (this->node_offset >= first_node_length) {
        throw std::runtime_error("WFAAlignment has offset to or past end of first node");
    }
    // When the base along the node hits this, we leave the node.
    // We track this ourselves to avoid repeated length queries.
    size_t node_end = first_node_length;

    // Walk through the edits
    auto edit_it = this->edits.begin();
    if (edit_it == this->edits.end()) {
        throw std::runtime_error("WFAAlignment has no edits");
    }
    // And track how much of the current edit has been resolved aleady.
    size_t current_edit_used = 0;

    // As we walk along, we build a mapping. Set up the first mapping
    Mapping* mapping_in_progress = result.add_mapping();
    // And set its position, with the offset
    mapping_in_progress->mutable_position()->set_node_id(graph.get_id(*path_it));
    mapping_in_progress->mutable_position()->set_is_reverse(graph.get_is_reverse(*path_it));
    mapping_in_progress->mutable_position()->set_offset(node_cursor);

    while (edit_it != this->edits.end()) {
        if (current_edit_used == edit_it->second) {
            // There's no edit left, but there ought to be as a loop invariant,
            // if no empty edits exist.
            throw std::runtime_error("WFAAlignment has empty edit");
        }

        // What kind of edit is it?
        auto& edit_type = edit_it->first;

        // And how much is left?
        size_t length_to_resolve = edit_it->second - current_edit_used;

        if (edit_type == match || edit_type == mismatch || edit_type == deletion) {
            // These edits consume some graph.
            // Make sure there is a graph node.
            if (path_it == this->path.end()) {
                throw std::runtime_error("WFAAlignment tried to go past end of path");
            }
            if (node_cursor == node_end) {
                // Make sure we aren't starting right at the end of the node.
                // We make sure this doesn't happen when we advance nodes, as
                // long as all nodes are nonempty.
                // TODO: Can we hit this from an anchored tail alignment somehow?
                throw std::runtime_error("WFAAlignment tried to go past end of node (" + std::to_string(node_end) + " bp)");
            }
            // Limit to length of graph node.
            length_to_resolve = std::min(length_to_resolve, node_end - node_cursor);
        }

        assert(length_to_resolve > 0);

#ifdef debug_path
        std::cerr << "Use " << length_to_resolve << " bp of " << edit_it->second << edit_it->first << " against node " << (path_it != this->path.end() ? graph.get_id(*path_it) : (nid_t)0) << " to go from edit " << (edit_it - edits.begin()) << " offset " << current_edit_used << " = path step " << (path_it - path.begin()) << " offset " << node_cursor;
#endif

        // Create a vg Edit to translate to
        vg::Edit* created = mapping_in_progress->add_edit();

        if (edit_type == match || edit_type == mismatch || edit_type == deletion) {
            // These edits consume some graph
            created->set_from_length(length_to_resolve);
            node_cursor += length_to_resolve;
        }
        if (edit_type == mismatch || edit_type == insertion) {
            // These edits carry sequence
            if (sequence_cursor + length_to_resolve > this->seq_offset + this->length) {
                throw std::runtime_error("WFAAlignment uses more sequence than provided");
            }
            created->set_sequence(sequence.substr(sequence_cursor, length_to_resolve));
        }
        if (edit_type == match || edit_type == mismatch || edit_type == insertion) {
            // These edits consume some sequence
            created->set_to_length(length_to_resolve);
            sequence_cursor += length_to_resolve;
        }

        // Now we've resolved at least part of this edit.
        current_edit_used += length_to_resolve;

        if (current_edit_used == edit_it->second) {
            // Finished the edit.
            // Reset to the start of the next edit.
            ++edit_it;
            current_edit_used = 0;
        }
        if (edit_type == match || edit_type == mismatch || edit_type == deletion) {
            // These edits consume some graph. So we may need to advance in the graph now.
            if (node_cursor == node_end) {
                // Finished the node.

#ifdef debug_path
                std::cerr << " (leave node at path step " << (path_it - path.begin()) << " offset " << node_cursor << ")";
#endif

                // We already checked above, and the path cursor isn't at the end.
                assert(path_it != this->path.end());

                // Reset to the start of the next node.
                node_cursor = 0;
                // And advance along the path if possible.
                ++path_it;
                if (path_it != this->path.end()) {
                    // We've reached a new node, so work out where the end is
                    node_end = graph.get_length(*path_it);

                    if (node_cursor == node_end) {
                        throw std::runtime_error("WFAAlignment has empty node " + std::to_string(graph.get_id(*path_it)));
                    }

                    // Also start a new Mapping
                    mapping_in_progress = result.add_mapping();
                    // And set its position
                    mapping_in_progress->mutable_position()->set_node_id(graph.get_id(*path_it));
                    mapping_in_progress->mutable_position()->set_is_reverse(graph.get_is_reverse(*path_it));
                    // The offset will always be 0 since we entered from somewhere.
                } else {
                    // No next node, so we should be at the end of what we use in the graph.
                    // If we try to use more graph, we will throw an error.
                    node_end = 0;
                }
            }
        }

#ifdef debug_path
        std::cerr << " to edit " << (edit_it - edits.begin()) << " offset " << current_edit_used << " = path step " << (path_it - path.begin()) << " offset " << node_cursor << std::endl;
#endif
    }

    return result;
}

std::ostream& WFAAlignment::print(const HandleGraph& graph, std::ostream& out) const {
    out << "{";
    if (!ok) {
        out << " NOT OK!";
    }
    out << " path = [";
    for (handle_t handle : this->path) {
        out << " (" << graph.get_id(handle) << ", " << graph.get_is_reverse(handle) << ")";
    }
    out << " ], edits = [ ";
    for (auto edit : this->edits) {
        out << edit.second << edit.first;
    }
    out << " ], node offset = " << this->node_offset;
    out << ", sequence range = [" << this->seq_offset << ", " << (this->seq_offset + this->length) << ")";
    out << ", score = " << this->score << " }";

    return out;
}

std::ostream& WFAAlignment::print(std::ostream& out) const {
    out << "{";
    if (!ok) {
        out << " NOT OK!";
    }
    out << " path = [";
    for (handle_t handle : this->path) {
        out << " (" << as_integer(handle) << ")";
    }
    out << " ], edits = [ ";
    // Print up to a manageable number of edits. Sometimes we can end up trying
    // to print apparently infinite edits and make many GB of logs.
    for (size_t i = 0; i < std::min((size_t) 100, this->edits.size()); i++) {
        auto edit = this->edits.at(i);
        out << edit.second << edit.first;
    }
    if (this->edits.size() > 100) {
        out << "...";
    }
    out << " ], node offset = " << this->node_offset;
    out << ", sequence range = [" << this->seq_offset << ", " << (this->seq_offset + this->length) << ")";
    out << ", score = " << this->score << " }";

    return out;
}

void WFAAlignment::check_lengths(const HandleGraph& graph) const {
    // Compute read and graph lengths from the edits
    size_t edit_graph_length = 0;
    size_t edit_read_length = 0;
    for (auto& e : edits) {
        if (e.first == match || e.first == mismatch || e.first == insertion) {
            // These edits use read sequence
            edit_read_length += e.second;
        }
        if (e.first == match || e.first == mismatch || e.first == deletion) {
            // These edits use graph sequence
            edit_graph_length += e.second;
        }
    }

    // Compute graph length from the path
    size_t path_graph_length = 0;
    for (auto& h : path) {
        path_graph_length += graph.get_length(h);
    }
    path_graph_length -= node_offset;

    if (edit_graph_length > path_graph_length) {
        // We want to use more graph than we got.
        print(graph, std::cerr);
        std::cerr << std::endl;
        throw std::runtime_error("WFAAlignment has path graph length " + std::to_string(path_graph_length) + " but edit graph length " + std::to_string(edit_graph_length));
    }
    if (edit_read_length != length) {
        // We want to use a different amount of read than we should.
        print(graph, std::cerr);
        std::cerr << std::endl;
        throw std::runtime_error("WFAAlignment has length " + std::to_string(length) + " but edit read length " + std::to_string(edit_read_length));
    }
}

std::ostream& operator<<(std::ostream& out, const WFAAlignment::Edit& edit) {
    return out << std::to_string(edit);
}

}

namespace std {

std::string to_string(const vg::WFAAlignment::Edit& edit) {
    switch (edit) {
    case vg::WFAAlignment::match:
        return "M";
        break;
    case vg::WFAAlignment::mismatch:
        return "X";
        break;
    case vg::WFAAlignment::insertion:
        return "I";
        break;
    case vg::WFAAlignment::deletion:
        return "D";
        break;
    default:
        throw std::runtime_error("Unknown edit operation");
    }
}

}

namespace vg {

//------------------------------------------------------------------------------

const WFAExtender::ErrorModel WFAExtender::default_error_model {
    WFAExtender::ErrorModel::default_mismatches(),
    WFAExtender::ErrorModel::default_gaps(),
    WFAExtender::ErrorModel::default_gap_length(),
    WFAExtender::ErrorModel::default_distance()
};

WFAExtender::WFAExtender() :
    graph(nullptr), mask("ACGT"), aligner(nullptr)
{
}

WFAExtender::WFAExtender(const gbwtgraph::GBWTGraph& graph, const Aligner& aligner, const ErrorModel& error_model) :
    graph(&graph), mask("ACGT"), aligner(&aligner), error_model(&error_model)
{
    // Check that the scoring parameters are reasonable.
    assert(this->aligner->match >= 0);
    assert(this->aligner->mismatch > 0);
    assert(this->aligner->gap_open >= this->aligner->gap_extension);
    assert(this->aligner->gap_extension > 0);

    // Check that the error model makes sense.
    for (auto& event : {
        this->error_model->mismatches,
        this->error_model->gaps,
        this->error_model->gap_length
    }) {
        assert(event.per_base >= 0);
        assert(event.min >= 0);
        assert(event.max >= event.min);
    }
}

//------------------------------------------------------------------------------

// A position in an alignment between a sequence and a graph.
struct MatchPos {
    uint32_t seq_offset;
    uint32_t node_offset;

    /// We need a stack-like type that can be copied or referred to quickly,
    /// since we consider MatchPos objects for each WFANode along a path when
    /// doing find_pos(), and they contain the path.
    /// This path type assumes that we only ever push and *then* pop.
    struct PathList {
        const static size_t NUM_INLINE = 4;
        size_t item_count = 0;
        uint32_t inline_items[NUM_INLINE];
        // TODO: There are probably some obscure situations where we could have a bug
        // because this is a shared pointer. Document them.
        std::shared_ptr<std::vector<uint32_t>> additional_items;

        void push(uint32_t value) {
            ++this->item_count;
            if (this->item_count == NUM_INLINE + 1) {
                this->additional_items.reset(new std::vector<uint32_t>());
            }
            if (this->item_count > NUM_INLINE) {
                this->additional_items->resize(this->item_count - NUM_INLINE);
            }
            this->top() = value;
        }

        bool empty() const {
            return this->item_count == 0;
        }

        size_t size() const {
            return this->item_count;
        }

        const uint32_t& top() const {
            if (this->item_count > NUM_INLINE) {
                return this->additional_items->at(this->item_count - NUM_INLINE - 1);
            } else {
                return this->inline_items[this->item_count - 1];
            }
        }

        uint32_t& top() {
            if (this->item_count > NUM_INLINE) {
                return this->additional_items->at(this->item_count - NUM_INLINE - 1);
            } else {
                return this->inline_items[this->item_count - 1];
            }
        }

        void pop() {
            --this->item_count;
            if (this->item_count == NUM_INLINE) {
                this->additional_items.reset();
            }
        }
    };

    PathList path; // Sequence of tree offsets from a leaf to the relevant node.

    // Creates an empty position.
    MatchPos() : seq_offset(0), node_offset(0) {}

    // Creates a position with the given offsets and path.
    MatchPos(uint32_t seq_offset, uint32_t node_offset, const PathList& path) : seq_offset(seq_offset), node_offset(node_offset), path(path) {}

    bool empty() const { return this->path.empty(); }
    bool at_last_node() const { return this->path.size() == 1; }
    uint32_t node() const { return this->path.top(); }
    void pop() { this->path.pop(); }

    // Returns the distance from the start of the alignment as seq_offset + target_offset,
    // assuming that we are on the given diagonal.
    std::int32_t distance(std::int32_t diagonal) {
        return 2 * std::int32_t(this->seq_offset) - diagonal;
    }

    // Positions are ordered by sequence offsets. Empty positions are smaller than
    // non-empty ones.
    bool operator<(const MatchPos& another) {
        if (this->empty()) {
            return !another.empty();
        }
        if (another.empty()) {
            return false;
        }
        return (this->seq_offset < another.seq_offset);
    }
};

// A point in an WFA score matrix (for a specific node). Must be aggregate-initializable.
struct WFAPoint {
    int32_t  score;
    int32_t  diagonal; // seq_offset - target offset
    uint32_t seq_offset;
    uint32_t node_offset;

    // Returns the offset in the target sequence.
    int32_t target_offset() const {
        return static_cast<int32_t>(this->seq_offset) - this->diagonal;
    }

    // Returns the four-parameter alignment score.
    int32_t alignment_score(const Aligner& aligner) const {
        return (static_cast<int32_t>(aligner.match) * (static_cast<int32_t>(this->seq_offset) + this->target_offset()) - this->score) / 2;
    }

    // Returns the four-parameter alignment score with an implicit final insertion.
    int32_t alignment_score(const Aligner& aligner, uint32_t final_insertion) const {
        return (static_cast<int32_t>(aligner.match) * (static_cast<int32_t>(this->seq_offset + final_insertion) + this->target_offset()) - this->score) / 2;
    }

    // Converts the point to an alignment position with the given path.
    MatchPos pos(const MatchPos::PathList& path) const {
        return MatchPos(this->seq_offset, this->node_offset, path);
    }

    // For ordering the points in WFANode.
    bool operator<(const WFAPoint& another) const {
        return (this->score < another.score || (this->score == another.score && this->diagonal < another.diagonal));
    }

    // We have a way to split up into a key and a value, and to be out back
    // together again.
    using key_type = std::pair<int32_t, int32_t>;
    using value_type = std::pair<uint32_t, uint32_t>;
    using map_entry_type = std::pair<key_type, value_type>;

    key_type key() const {
        return key_type(this->score, this->diagonal);
    }

    value_type value() const {
        return value_type(this->seq_offset, this->node_offset);
    }

    static WFAPoint from_key_value(const key_type& key, const value_type& value) {
        return {key.first, key.second, value.first, value.second};
    }

    static WFAPoint from_map_entry(const map_entry_type& entry) {
        return from_key_value(entry.first, entry.second);
    }
};

//------------------------------------------------------------------------------

/*
    WFANode corresponds to a path in a graph and a set of haplotypes that follow
    the path. We store the path itself, the GBWT search state at the end of the
    path, and the concatenation of node sequences. We additionally store links
    to the parent and child nodes in WFATree, using offsets in the vector of
    nodes.

    The path continues until the haplotypes branch, it exceeds 1024 bp, or it
    covers the target position.
*/
struct WFANode {
    // The path in the graph.
    gbwt::vector_type path;

    // GBWT search state at the end of the path.
    gbwt::SearchState state;

    // Concatenation of node sequences.
    std::string node_sequence;

    // Offsets in the vector of nodes.
    std::uint32_t parent;
    std::vector<std::uint32_t> children;

    // Offset for the target position in the sequence.
    std::uint32_t target_offset;

    // All haplotypes end here.
    bool dead_end;

    constexpr static size_t TARGET_LENGTH = 1024;

    constexpr static size_t MATCHES = 0;
    constexpr static size_t INSERTIONS = 1; // characters in the sequence but not in the graph
    constexpr static size_t DELETIONS = 2;  // characters in the graph but not in the sequence

    // Points on the wavefronts are indexed by score, diagonal.
    std::array<hash_map<WFAPoint::key_type, WFAPoint::value_type>, 3> wavefronts;

    WFANode(const gbwtgraph::CachedGBWTGraph& graph, const gbwt::SearchState& state, pos_t target, std::uint32_t parent) :
        parent(parent), target_offset(std::numeric_limits<std::uint32_t>::max()), dead_end(false)
    {
        if(this->append_node(graph, state, target)) {
            return;
        }

        while (this->node_sequence.length() < TARGET_LENGTH) {
            size_t successors = 0;
            gbwt::SearchState next_state;
            graph.follow_paths(this->state, [&](const gbwt::SearchState& next) {
                successors++;
                next_state = next;
                return true;
            });
            if (successors == 0) {
                this->dead_end = true;
                break;
            } else if (successors > 1) {
                break;
            }
            if (this->append_node(graph, next_state, target)) {
                break;
            }
        }
    }

    size_t length() const { return this->node_sequence.length(); }
    bool is_leaf() const { return (this->children.empty() || this->dead_end); }
    bool expanded() const { return (!this->children.empty() || this->dead_end); }

    // WFANode::find_pos
    // Returns the position for the given score and diagonal with the given path,
    // or an empty position if it does not exist.
    MatchPos find_pos(size_t type, int32_t score, int32_t diagonal, const MatchPos::PathList& path) const {
        WFAPoint::key_type key { score, diagonal };
        auto& points = this->wavefronts[type];
        // Find the item with the same score and diagonal
        auto iter = points.find(key);
        if (iter == points.end()) {
            // Nothing is stored for this score and diagonal
            return MatchPos();
        }
        // Otherwise, something is stored for this score and diagonal. Reconstitute it.
        WFAPoint point = WFAPoint::from_key_value(key, iter->second);
        return point.pos(path);
    }

    // Update the WFA matrix with the given alignment position.
    void update(size_t type, int32_t score, int32_t diagonal, const MatchPos& pos) {
        this->update(type, score, diagonal, pos.seq_offset, pos.node_offset);
    }

    // Update the WFA matrix with the given offsets.
    void update(size_t type, int32_t score, int32_t diagonal, uint32_t seq_offset, uint32_t node_offset) {
        // TODO: Make the WFAPoint and then split it up with its methods instead?
        WFAPoint::key_type key { score, diagonal };
        WFAPoint::value_type value { seq_offset, node_offset };
        auto& points = this->wavefronts[type];
        auto iter = points.find(key);
        if (iter == points.end()) {
            // This is a new score and diagonal
            points.emplace_hint(iter, key, value);
        } else {
            // This score and diagonal already exists, so overwrite the value
            iter->second = value;
        }
    }

    // Advances the position to the first non-match at or after the current position.
    void match_forward(const std::string& sequence, MatchPos& pos) const {
        while (
            pos.seq_offset < sequence.length() &&
            pos.node_offset < this->node_sequence.length() &&
            sequence[pos.seq_offset] == this->node_sequence[pos.node_offset]
        ) {
            pos.seq_offset++;
            pos.node_offset++;
        }
    }

private:
    // Proceed to the given search state and append the path. Returns true if we reached the target.
    bool append_node(const gbwtgraph::CachedGBWTGraph& graph, gbwt::SearchState next, pos_t target) {
        this->state = next;
        this->path.push_back(this->state.node);
        gbwtgraph::view_type view = graph.get_sequence_view(gbwtgraph::GBWTGraph::node_to_handle(this->state.node));
        this->node_sequence.append(view.first, view.second);
        if (gbwt::Node::encode(id(target), is_rev(target)) == this->state.node) {
            this->target_offset = this->node_sequence.length() - (view.second - offset(target));
            return true;
        }
        return false;
    }
};

//------------------------------------------------------------------------------

/*
    WFATree represents a trie of haplotypes starting from a given position in
    the graph. The tree is expanded lazily as needed, and WFA alignment is done
    over all haplotypes in the tree. Each node is a WFANode that represents a
    non-branching set of haplotypes over a path in the graph. 
*/
class WFATree {
public:
    const gbwtgraph::CachedGBWTGraph graph;
    const std::string& sequence;

    // Start and end positions in the graph (exclusive).
    pos_t from, to;

    // Node identifiers are offsets in this vector. Node 0 is the root.
    std::vector<WFANode> nodes;

    // Best alignment found so far. If we reached the target position in the
    // graph, the score includes the implicit insertion at the end but the
    // position itself does not.
    WFAPoint candidate_point;
    std::uint32_t candidate_node;

    // WFA score (penalty) parameters derived from the actual scoring parameters.
    std::int32_t mismatch, gap_open, gap_extend;

    // Stop if no alignment has been found with this score or less.
    std::int32_t score_bound;

    // Furthest distance (seq_offset + target_offset) reached so far.
    std::int32_t max_distance;

    // Minimum distance we are still considering.
    std::int32_t min_distance;

    struct ScoreProperties {
        std::int32_t min_diagonal;
        std::int32_t max_diagonal;
        bool reachable_with_gap;
    };

    // A set of possible scores and diagonals reached with them. Because we
    // derive the three WFA scoring parameters from four Aligner parameters
    // (that include a match bonus), many small scores are impossible and can
    // be skipped.
    std::map<std::int32_t, ScoreProperties> possible_scores;

    WFATree(
        const gbwtgraph::GBWTGraph& graph, const std::string& sequence,
        pos_t from, pos_t to,
        const Aligner& aligner, const WFAExtender::ErrorModel& error_model
    ) :
        graph(graph), sequence(sequence), from(from), to(to),
        nodes(),
        candidate_point({ std::numeric_limits<std::int32_t>::max(), 0, 0, 0 }), candidate_node(0),
        mismatch(2 * (aligner.match + aligner.mismatch)),
        gap_open(2 * (aligner.gap_open - aligner.gap_extension)),
        gap_extend(2 * aligner.gap_extension + aligner.match),
        score_bound(0), max_distance(0), min_distance(0),
        possible_scores()
    {
        // Create the root node based on the starting position. Because the start
        // is outside the alignment, we may already have exhausted the node.
        handle_t handle = this->graph.get_handle(id(this->from), is_rev(this->from));
        gbwt::SearchState state = this->graph.get_state(handle);
        WFANode root(this->graph, state, this->to, 0);
        root.update(WFANode::MATCHES, 0, 0, 0, offset(this->from) + 1);
        this->nodes.push_back(root);

        // Determine score bound based on the error model and sequence length.
        std::int32_t max_mismatches = error_model.mismatches.evaluate(sequence.length());
        std::int32_t max_gaps = error_model.gaps.evaluate(sequence.length());
        std::int32_t max_gap_length =  error_model.gap_length.evaluate(sequence.length());
        this->score_bound = max_mismatches * this->mismatch + max_gaps * this->gap_open + max_gap_length * this->gap_extend;

        this->possible_scores[0] = { 0, 0, false };
    }

    std::uint32_t size() const { return this->nodes.size(); }
    static bool is_root(std::uint32_t node) { return (node == 0); }
    uint32_t parent(std::uint32_t node) const { return this->nodes[node].parent; }

    // Assumes length > 0.
    std::int32_t gap_extend_penalty(std::uint32_t length) const {
        return std::int32_t(length) * this->gap_extend;
    }

    // Assumes length > 0.
    std::int32_t gap_penalty(std::uint32_t length) const {
        return this->gap_open + this->gap_extend_penalty(length);
    }

    // wf_extend() in the paper.
    // If we reach the end of a node, we continue to the start of the next node even
    // if we do not use any characters in it.
    void extend(std::int32_t score) {
        auto iter = this->possible_scores.find(score);
        if (iter == this->possible_scores.end()) {
            // TODO: This should not happen. Should there be a warning/error?
            return;
        }
        for (std::int32_t diagonal = iter->second.min_diagonal; diagonal <= iter->second.max_diagonal; diagonal++) {
            std::vector<std::uint32_t> leaves = this->get_leaves();
            this->extend_over(score, diagonal, leaves);
        }
    }

    // Returns the next possible score after the given score, assuming that the
    // given score has already been determined possible. Also updates the set
    // of possible scores with those reachable from the given score but does not
    // set the diagonal ranges for them.
    std::int32_t next_score(std::int32_t match_score) {
        // This score + a mismatch is a possible score.
        std::int32_t mismatch_score = match_score + this->mismatch;
        if (this->possible_scores.find(mismatch_score) == this->possible_scores.end()) {
            this->possible_scores[mismatch_score] = { 0, 0, false };
        }

        // This score + gap extend is a possible score reachable by a gap, if
        // this score was also reachable by a gap.
        auto match_iter = this->possible_scores.find(match_score);
        if (match_iter->second.reachable_with_gap) {
            std::int32_t extend_score = match_score + this->gap_extend;
            auto extend_iter = this->possible_scores.find(extend_score);
            if (extend_iter != this->possible_scores.end()) {
                extend_iter->second.reachable_with_gap = true;
            } else {
                this->possible_scores[extend_score] = { 0, 0, true };
            }
        }

        // This score + gap open + gap extend is a possible score reachable by a gap.
        std::int32_t open_score = match_score + this->gap_open + this->gap_extend;
        auto open_iter = this->possible_scores.find(open_score);
        if (open_iter != this->possible_scores.end()) {
            open_iter->second.reachable_with_gap = true;
        } else {
            this->possible_scores[open_score] = { 0, 0, true };
        }

        // We know that there are further values beyond match_score.
        ++match_iter;
        return match_iter->first;
    }

    // wf_next() in the paper, with a variant of wf_reduce().
    // If we reach the end of a node, we continue to the start of the next node even
    // if we do not use any characters in it.
    void next(std::int32_t score) {
        // diagonal_range is the potential range of diagonals for this score, based on the
        // diagonals of previous scores.
        std::pair<std::int32_t, std::int32_t> diagonal_range = this->get_diagonals(score);
        // actual_range is the range of diagonals with points that reach far enough
        // with this score.
        std::pair<std::int32_t, std::int32_t> actual_range = empty_diagonal_range();
        for (std::int32_t diagonal = diagonal_range.first; diagonal <= diagonal_range.second; diagonal++) {
            std::vector<std::uint32_t> leaves = this->get_leaves();
            // NOTE: We may do the same updates from multiple leaves.
            for (std::uint32_t leaf : leaves) {
                MatchPos ins = this->ins_predecessor(leaf, score, diagonal).first;
                if (!ins.empty()) {
                    ins.seq_offset++;
                    if (ins.distance(diagonal) >= this->min_distance) {
                        this->nodes[ins.node()].update(WFANode::INSERTIONS, score, diagonal, ins);
                        adjust_diagonal_range(actual_range, diagonal);
                    }
                }

                MatchPos del = this->del_predecessor(leaf, score, diagonal).first;
                if (!del.empty()) {
                    this->successor_offset(del);
                    if (del.distance(diagonal) >= this->min_distance) {
                        this->nodes[del.node()].update(WFANode::DELETIONS, score, diagonal, del);
                        adjust_diagonal_range(actual_range, diagonal);
                    }
                    this->expand_if_necessary(del);
                }

                MatchPos subst = this->find_pos(WFANode::MATCHES, leaf, score - this->mismatch, diagonal, true, true);
                if (!subst.empty()) {
                    subst.seq_offset++;
                    this->successor_offset(subst);
                    this->expand_if_necessary(subst);
                }

                // Determine the edit that reaches furthest on the diagonal.
                if (subst < ins) {
                    subst = std::move(ins);
                }
                if (subst < del) {
                    subst = std::move(del);
                }

                if (!subst.empty()) {
                    // If we reached the end position with the edit, we get a candidate
                    // alignment by assuming that the rest of the sequence is an insertion.
                    // If the edit is an insertion, we charge the gap open cost again, but
                    // we already got the same insertion without the extra cost from the
                    // match preceding the insertion.
                    WFANode& node = this->nodes[subst.node()];
                    if (subst.node_offset == node.target_offset) {
                        std::uint32_t gap_length = this->sequence.length() - subst.seq_offset;
                        std::int32_t gap_score = 0;
                        if (gap_length > 0) {
                            gap_score = this->gap_penalty(gap_length);
                        }
                        if (score + gap_score < this->candidate_point.score) {
                            this->candidate_point = { score + gap_score, diagonal, subst.seq_offset, subst.node_offset };
                            this->candidate_node = subst.node();
                        }
                    }
                    if (subst.distance(diagonal) >= this->min_distance) {
                        node.update(WFANode::MATCHES, score, diagonal, subst);
                        adjust_diagonal_range(actual_range, diagonal);
                    }
                }
            }
        }

        // Set the diagonal range for the score.
        auto iter = this->possible_scores.find(score);
        if (iter != this->possible_scores.end()) {
            iter->second.min_diagonal = actual_range.first;
            iter->second.max_diagonal = actual_range.second;
        }
    }

    // Returns the predecessor position for the furthest reaching insertion for
    // (score, diagonal) at the specified node or its ancestors, or an empty position
    // if it does not exist. Also returns the type of the predecessor.
    std::pair<MatchPos, WFAAlignment::Edit> ins_predecessor(uint32_t node, int32_t score, int32_t diagonal) const {
        MatchPos open = this->find_pos(WFANode::MATCHES, node, score - this->gap_open - this->gap_extend, diagonal - 1, true, false);
        MatchPos extend = this->find_pos(WFANode::INSERTIONS, node, score - this->gap_extend, diagonal - 1, true, false);
        return (open < extend ? std::make_pair(extend, WFAAlignment::insertion) : std::make_pair(open, WFAAlignment::match));
    }

    // Returns the predecessor position for the furthest reaching deletion for
    // (score, diagonal) at the specified node or its ancestors, or an empty position
    // if it does not exist. Also returns the type of the predecessor.
    std::pair<MatchPos, WFAAlignment::Edit> del_predecessor(uint32_t node, int32_t score, int32_t diagonal) const {
        MatchPos open = this->find_pos(WFANode::MATCHES, node, score - this->gap_open - this->gap_extend, diagonal + 1, false, true);
        MatchPos extend = this->find_pos(WFANode::DELETIONS, node, score - this->gap_extend, diagonal + 1, false, true);
        return (open < extend ? std::make_pair(extend, WFAAlignment::deletion) : std::make_pair(open, WFAAlignment::match));
    }

    // Returns the predecessor position for the furthest reaching run of matches
    // for (score, diagonal) at the specified node or its ancestors, or an empty
    // position if it does not exist. Also returns the type of the predecessor.
    std::pair<MatchPos, WFAAlignment::Edit> match_predecessor(uint32_t node, int32_t score, int32_t diagonal) const {
        MatchPos ins = this->find_pos(WFANode::INSERTIONS, node, score, diagonal, false, false);
        MatchPos del = this->find_pos(WFANode::DELETIONS, node, score, diagonal, false, false);
        MatchPos subst = this->find_pos(WFANode::MATCHES, node, score - this->mismatch, diagonal, false, false);
        if (!subst.empty()) {
            subst.seq_offset++;
            subst.node_offset++;
        }

        if (ins < del) {
            return (del < subst ? std::make_pair(subst, WFAAlignment::mismatch) : std::make_pair(del, WFAAlignment::deletion));
        } else {
            return (ins < subst ? std::make_pair(subst, WFAAlignment::mismatch) : std::make_pair(ins, WFAAlignment::insertion));
        }
    }

    // Move forward on the path corresponding to the position. If the position points
    // to the end of a node, we assume that it has not reached the end of the path.
    void successor_offset(MatchPos& pos) const {
        if (pos.node_offset >= this->nodes[pos.node()].length()) {
            pos.pop(); pos.node_offset = 0;
        }
        pos.node_offset++;
    }

    // Updates the node and an offset in it to the predecessor offset.
    void predecessor_offset(uint32_t& node, uint32_t& offset) const {
        if (offset > 0) {
            offset--;
        } else {
            node = this->parent(node);
            offset = this->nodes[node].length() - 1;
        }
    }

    // Returns true if the position is empty.
    static bool no_pos(pos_t pos) { return (id(pos) == 0); }

    // Replaces the candidate with the partial alignment with the highest alignment
    // score according to the aligner.
    void trim(const Aligner& aligner) {
        this->candidate_point = { 0, 0, 0, 0 };
        this->candidate_node = 0;
        std::int32_t best_score = 0;
        for (std::uint32_t node = 0; node < this->size(); node++) {
            for (auto& point_entry : this->nodes[node].wavefronts[WFANode::MATCHES]) {
                // Scan all stored points on the node.
                // TODO: Does iteration order matter?
                // Convert map entries to points.
                WFAPoint point = WFAPoint::from_map_entry(point_entry);
                std::int32_t alignment_score = point.alignment_score(aligner);
                if (alignment_score > best_score) {
                    // This is a new winner
                    this->candidate_point = point;
                    this->candidate_node = node;
                    best_score = alignment_score;
                }
            }
        }
    }

private:

    // wf_extend() on a specific diagonal for the set of (local) haplotypes corresponding to
    // the given list of leaves in the WFATree.
    void extend_over(std::int32_t score, std::int32_t diagonal, const std::vector<std::uint32_t>& leaves) {
        for (std::uint32_t leaf : leaves) {
            MatchPos pos = this->find_pos(WFANode::MATCHES, leaf, score, diagonal, false, false);
            if (pos.empty()) {
                continue; // An impossible score / diagonal combination.
            }
            while (true) {
                WFANode& node = this->nodes[pos.node()];
                bool may_reach_target = false;
                if (node.target_offset >= pos.node_offset && node.target_offset < node.length()) {
                    // This node covers the target, and we have not matched the offset yet.
                    may_reach_target = true;
                }
                node.match_forward(this->sequence, pos);

                // We got a match that reached the end or went past it.
                // Alternatively there is no end position and we have aligned the entire sequence.
                // This gives us a candidate where the rest of the sequence is an insertion.
                if (
                    (may_reach_target && pos.node_offset >= node.target_offset) ||
                    (no_pos(this->to) && pos.seq_offset >= this->sequence.length())
                ) {
                    // If we managed to match the target position, it is part of the overshoot.
                    std::uint32_t overshoot = (no_pos(to) ? 0 : pos.node_offset - node.target_offset);
                    std::uint32_t gap_length = (this->sequence.length() - pos.seq_offset) + overshoot;
                    std::int32_t gap_score = 0;
                    if (gap_length > 0) {
                        gap_score = this->gap_penalty(gap_length);
                    }

                    if (score + gap_score < this->candidate_point.score) {
                        this->candidate_point = {
                            score + gap_score, diagonal,
                            pos.seq_offset - overshoot, node.target_offset
                        };
                        this->candidate_node = pos.node();
                    }
                }

                // We may have matched some additional bases in the node.
                this->max_distance = std::max(this->max_distance, pos.distance(diagonal));
                node.update(WFANode::MATCHES, score, diagonal, pos);
                if (pos.node_offset < node.length()) {
                    break;
                }

                // We reached the end of the node. If the position does not specify the
                // path we should follow, we continue recursively in all children.
                this->expand_if_necessary(pos); // NOTE: Possibly invalidates `node`.
                if (pos.at_last_node()) {
                    assert(pos.node() == leaf);
                    // Create a copy of the child list, because further pushes to the
                    // node vector may invalidate the reference.
                    std::vector<std::uint32_t> new_leaves = this->nodes[pos.node()].children;
                    this->extend_over(score, diagonal, new_leaves);
                    break;
                }
                pos.pop();
                pos.node_offset = 0;
            }
        }
    }

    // TODO: Can we trim subtrees with no points far enough?
    std::vector<std::uint32_t> get_leaves() const {
        std::vector<std::uint32_t> leaves;
        for (std::uint32_t node = 0; node < this->size(); node++) {
            if (this->nodes[node].is_leaf()) {
                leaves.push_back(node);
            }
        }
        return leaves;
    }

    // Adjusts the diagonal range to include the given diagonal.
    static void adjust_diagonal_range(std::pair<std::int32_t, std::int32_t>& range, std::int32_t diagonal) {
        range.first = std::min(range.first, diagonal);
        range.second = std::max(range.second, diagonal);
    }

    // Updates the given diagonal range to include all diagonals covered by the
    // range for the given score.
    std::pair<std::int32_t, std::int32_t> update_diagonal_range(
        std::pair<std::int32_t, std::int32_t> range,
        std::int32_t score) const
    {
        if (score >= 0) {
            auto iter = this->possible_scores.find(score);
            if (iter != this->possible_scores.end()) {
                range.first = std::min(range.first, iter->second.min_diagonal);
                range.second = std::max(range.second, iter->second.max_diagonal);
            }
        }
        return range;
    }

    static std::pair<int32_t, int32_t> empty_diagonal_range() {
        return std::make_pair(std::numeric_limits<int32_t>::max(), std::numeric_limits<int32_t>::min());
    }

    // Determines the diagonal range for the given score.
    // Assumes that the score is valid.
    // Returns an empty range if the score is impossible.
    std::pair<int32_t, int32_t> get_diagonals(int32_t score) {
        // Determine the diagonal range for the given score.
        std::pair<int32_t, int32_t> range = empty_diagonal_range();
        range = this->update_diagonal_range(range, score - this->mismatch); // Mismatch.
        range = this->update_diagonal_range(range, score - this->gap_open - this->gap_extend); // New gap.
        range = this->update_diagonal_range(range, score - this->gap_extend); // Extend an existing gap.
        if (range.first > range.second) {
            return range;
        }
        range.first--; range.second++;
        return range;
    }

    // If we have reached the end of the current node, expand its children if necessary.
    // Call this whenever the alignment advances in the node.
    void expand_if_necessary(const MatchPos& pos) {
        // NOTE: Pushes to the node vector may invalidate references to the node.
        std::uint32_t node = pos.node();
        if (this->nodes[node].expanded() || pos.node_offset < this->nodes[node].length()) {
            return;
        }
        bool found = false;
        this->graph.follow_paths(this->nodes[node].state, [&](const gbwt::SearchState& child) -> bool {
            this->nodes[node].children.push_back(this->size());
            this->nodes.emplace_back(this->graph, child, this->to, node);
            found = true;
            return true;
        });
        if (!found) {
            this->nodes[node].dead_end = true;
        }
    }

    // WFATree::find_pos
    // Returns the furthest position in given WFA matrix for (score, diagonal) at the
    // specified node or its ancestors, or an empty position if it does not exist.
    // Returns an empty position if an extendable position is requested but the position
    // cannot be extended.
    MatchPos find_pos(size_t type, uint32_t node, int32_t score, int32_t diagonal, bool extendable_seq, bool extendable_graph) const {
        if (score < 0) {
            return MatchPos();
        }
        MatchPos::PathList path;
        while(true) {
            // Make the path a bit longer.
            path.push(node);
            // Find a position at this node.
            // The MatchPos will need to know the whole path, so we can return it.
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
            if (this->is_root(node)) {
                return MatchPos();
            }
            node = this->parent(node);
        }
    }

    // Assumes that the position is non-empty.
    bool at_dead_end(const MatchPos& pos) const {
        return (this->nodes[pos.node()].dead_end && pos.node_offset >= this->nodes[pos.node()].length());
    }
};

//------------------------------------------------------------------------------

//#define debug_connect

WFAAlignment WFAExtender::connect(std::string sequence, pos_t from, pos_t to) const {
    if (this->graph == nullptr || this->aligner == nullptr) {
#ifdef debug_connect
        std::cerr << "No graph or no aligner! Returning empty alignment!" << std::endl;
#endif
        return WFAAlignment();
    }
    if (!this->graph->has_node(id(from))) {
#ifdef debug_connect
        std::cerr << "No start node! Returning empty alignment!" << std::endl;
#endif
        return WFAAlignment();
    }

    this->mask(sequence);
    WFATree tree(*(this->graph), sequence, from, to, *(this->aligner), *(this->error_model));

    std::int32_t score = 0;
    while (true) {
#ifdef debug_connect
        std::cerr << "Extend for score " << score << std::endl;
#endif
        tree.extend(score);

        // Update the minimum distance for points we will still consider.
        std::int32_t distance_band = this->error_model->distance.evaluate(sequence.length());
        if (distance_band < tree.max_distance) {
            tree.min_distance = tree.max_distance - distance_band;
        }
#ifdef debug_connect
        std::cerr << "Max distance is " << tree.max_distance << ", min distance is " << tree.min_distance << std::endl;
#endif

        if (tree.candidate_point.score <= score) {
            break;
        }

        score = tree.next_score(score);
#ifdef debug_connect
        std::cerr << "Next score will be " << score << std::endl;
#endif

        if (score > tree.score_bound) {
            break;
        }

#ifdef debug_connect
        std::cerr << "Next for score " << score << std::endl;
#endif
        tree.next(score);
    }

    // If we do not have a full-length alignment within the score bound,
    // we find the best partial alignment if there was no destination or
    // return an empty alignment otherwise.
    bool full_length = true;
    std::uint32_t unaligned_tail = sequence.length() - tree.candidate_point.seq_offset;
    if (tree.candidate_point.score > tree.score_bound) {
#ifdef debug_connect
        std::cerr << "No alignment could be found under score bound of " << tree.score_bound << "; best found was " << tree.candidate_point.score << std::endl;
#endif
        unaligned_tail = 0;
        if (WFATree::no_pos(to)) {
            tree.trim(*(this->aligner));
            full_length = false;
        } else {
            return WFAAlignment();
        }
    }

    // Start building an alignment. Store the path first.
    // No need to convert the node offset because it is from the root state's node start.
    WFAAlignment result {
        {}, {}, std::uint32_t(offset(from) + 1), 0,
        tree.candidate_point.seq_offset + unaligned_tail,
        tree.candidate_point.alignment_score(*(this->aligner), unaligned_tail),
        true
    };
    std::uint32_t node = tree.candidate_node;
    while (true) {
        // Go back up the tree and compose the path in reverse order
        for (auto iter = tree.nodes[node].path.rbegin(); iter != tree.nodes[node].path.rend(); ++iter) {
            result.path.push_back(gbwtgraph::GBWTGraph::node_to_handle(*iter));
        }
        if (tree.is_root(node)) {
            // Stop when we reach the root
            break;
        }
        // Otherwise go to the parent
        node = tree.parent(node);
    }
    std::reverse(result.path.begin(), result.path.end());

    // We have a full-length alignment within the score bound with an implicit insertion at the end.
    WFAPoint point = tree.candidate_point;
    node = tree.candidate_node;
    if (unaligned_tail > 0) {
        std::uint32_t final_insertion = sequence.length() - tree.candidate_point.seq_offset;
        result.append(WFAAlignment::insertion, final_insertion);
        point.score -= tree.gap_penalty(unaligned_tail);
    }

    // Backtrace the edits.
    WFAAlignment::Edit edit = WFAAlignment::match;
    while (point.seq_offset > 0 || point.diagonal != 0) {
        std::pair<MatchPos, WFAAlignment::Edit> predecessor;
        switch (edit)
        {
        case WFAAlignment::match:
            predecessor = tree.match_predecessor(node, point.score, point.diagonal);
            result.append(WFAAlignment::match, point.seq_offset - predecessor.first.seq_offset);
            point.seq_offset = predecessor.first.seq_offset;
            point.node_offset = predecessor.first.node_offset;
            if (!predecessor.first.empty()) {
                node = predecessor.first.node();
            }
            edit = predecessor.second;
            break;
        case WFAAlignment::mismatch:
            result.append(WFAAlignment::mismatch, 1);
            point.seq_offset--;
            tree.predecessor_offset(node, point.node_offset);
            point.score -= tree.mismatch;
            edit = WFAAlignment::match;
            break;
        case WFAAlignment::insertion:
            predecessor = tree.ins_predecessor(node, point.score, point.diagonal);
            result.append(WFAAlignment::insertion, 1);
            point.seq_offset--;
            if (predecessor.second == WFAAlignment::insertion) {
                point.score -= tree.gap_extend;
            } else {
                point.score -= tree.gap_open + tree.gap_extend;
            }
            point.diagonal--;
            edit = predecessor.second;
            break;
        case WFAAlignment::deletion:
            predecessor = tree.del_predecessor(node, point.score, point.diagonal);
            result.append(WFAAlignment::deletion, 1);
            tree.predecessor_offset(node, point.node_offset);
            if (predecessor.second == WFAAlignment::deletion) {
                point.score -= tree.gap_extend;
            } else {
                point.score -= tree.gap_open + tree.gap_extend;
            }
            point.diagonal++;
            edit = predecessor.second;
            break;
        }
    }
    std::reverse(result.edits.begin(), result.edits.end());

    // We used "from + 1" as the starting position for the alignment. That could have
    // been a past-the-end position in the initial node. Once we have an actual path
    // instead of a tree of potential paths, we can remove the unused node.
    if (!result.path.empty() && result.node_offset >= this->graph->get_length(result.path.front())) {
        result.path.erase(result.path.begin());
        result.node_offset = 0;
    }

    // Due to the way we expand the tree of GBWT search states and store wavefront
    // information in the leaves, and how we coalesce graph nodes into WFANode
    // objects, we sometimes do not use any bases in trailing graph nodes.
    // We deal with this now to avoid facing the issue later.
    int64_t final_offset = result.final_offset(*(this->graph));

    while ((result.path.size() == 1 && final_offset == result.node_offset) || (result.path.size() > 1 && final_offset <= 0)) {
        // No bases actually used on the final node, so drop it and adjust the offset.
        result.path.pop_back();
        if (!result.path.empty()) {
            // Offset should now be relative to the start of this node, and not end of this node/start of node after
            final_offset += this->graph->get_length(result.path.back());

        } else {

        }
    }

#ifdef debug_connect
    std::cerr << "Found an alignment." << std::endl;
#endif
    return result;
}

WFAAlignment WFAExtender::suffix(const std::string& sequence, pos_t from) const {
    WFAAlignment result = this->connect(sequence, from, pos_t(0, false, 0));

    if (!result.edits.empty() && result.length == sequence.length() && (result.edits.back().first == WFAAlignment::match || result.edits.back().first == WFAAlignment::mismatch)) {
        // The alignment used all of the sequence and has a match/mismatch at the appropriate end
        result.score += this->aligner->full_length_bonus; 
    }

    return result;
}

WFAAlignment WFAExtender::prefix(const std::string& sequence, pos_t to) const {
    if (this->graph == nullptr) {
        return WFAAlignment();
    }

    // Flip the position, extend forward, and reverse the return value.
    to = reverse_base_pos(to, this->graph->get_length(this->graph->get_handle(id(to), is_rev(to))));
    WFAAlignment result = this->connect(reverse_complement(sequence), to, pos_t(0, false, 0));
    result.flip(*(this->graph), sequence);

    if (!result.edits.empty() && result.length == sequence.length() && (result.edits.front().first == WFAAlignment::match || result.edits.front().first == WFAAlignment::mismatch)) {
        result.score += this->aligner->full_length_bonus; 
    }

    return result;
}

//------------------------------------------------------------------------------

} // namespace vg
