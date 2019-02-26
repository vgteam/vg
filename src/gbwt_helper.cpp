#include "gbwt_helper.hpp"

#include <stack>

namespace vg {

//------------------------------------------------------------------------------

GBWTGraph::GBWTGraph(const gbwt::GBWT& gbwt_index, const HandleGraph& sequence_source) :
    index(gbwt_index), total_nodes(0) {

    // Sanity checks for the GBWT index.
    assert(this->index.bidirectional());

    // Determine the number of real nodes and the total length of the sequences.
    // Node n is real, if real_nodes[node_offset(n) / 2] is true.
    size_t total_length = 0, potential_nodes = this->index.sigma() - this->index.firstNode();
    this->real_nodes = std::vector<bool>(potential_nodes / 2, false);
    std::vector<handle_t> handle_cache(potential_nodes); // Getting handles from XG is slow.
    for (gbwt::node_type node = this->index.firstNode(); node < this->index.sigma(); node += 2) {
        if (this->index.empty(node)) {
            continue;
        }
        size_t offset = this->node_offset(node);
        this->real_nodes[offset / 2] = true;
        this->total_nodes++;
        handle_t source_handle = sequence_source.get_handle(gbwt::Node::id(node), false);
        total_length += sequence_source.get_length(source_handle);
        handle_cache[offset] = source_handle;
    }
    this->sequences.reserve(2 * total_length);
    this->offsets = sdsl::int_vector<0>(potential_nodes + 1, 0, gbwt::bit_length(this->sequences.capacity()));

    // Store the concatenated sequences and their offset ranges for both orientations of all nodes.
    // Given GBWT node n, the sequence is sequences[node_offset(n)] to sequences[node_offset(n + 1) - 1].
    for (gbwt::node_type node = this->index.firstNode(); node < this->index.sigma(); node += 2) {
        std::string seq;
        size_t offset = this->node_offset(node);
        if (this->real_nodes[offset / 2]) {
            seq = sequence_source.get_sequence(handle_cache[offset]);
        }
        this->sequences.insert(this->sequences.end(), seq.begin(), seq.end());
        this->offsets[offset + 1] = this->sequences.size();
        seq = reverse_complement(seq);
        this->sequences.insert(this->sequences.end(), seq.begin(), seq.end());
        this->offsets[offset + 2] = this->sequences.size();
    }
}

GBWTGraph::GBWTGraph(const GBWTGraph& source) :
    index(source.index) {
    this->sequences = source.sequences;
    this->offsets = source.offsets;
    this->real_nodes = source.real_nodes;
    this->total_nodes = source.total_nodes;
}

GBWTGraph::GBWTGraph(GBWTGraph&& source) :
    index(source.index) {
    this->sequences = std::move(source.sequences);
    this->offsets = std::move(source.offsets);
    this->real_nodes = std::move(source.real_nodes);
    this->total_nodes = std::move(total_nodes);
}

//------------------------------------------------------------------------------

bool GBWTGraph::has_node(id_t node_id) const {
    size_t offset = this->node_offset(gbwt::Node::encode(node_id, false)) / 2;
    return (offset < this->real_nodes.size() && this->real_nodes[offset]);
}

handle_t GBWTGraph::get_handle(const id_t& node_id, bool is_reverse) const {
    return as_handle(gbwt::Node::encode(node_id, is_reverse));
}

id_t GBWTGraph::get_id(const handle_t& handle) const {
    return gbwt::Node::id(as_integer(handle));
}

bool GBWTGraph::get_is_reverse(const handle_t& handle) const {
    return gbwt::Node::is_reverse(as_integer(handle));
}

handle_t GBWTGraph::flip(const handle_t& handle) const {
    return as_handle(gbwt::Node::reverse(as_integer(handle)));
}

size_t GBWTGraph::get_length(const handle_t& handle) const {
    size_t offset = this->node_offset(handle);
    return this->offsets[offset + 1] - this->offsets[offset];
}

std::string GBWTGraph::get_sequence(const handle_t& handle) const {
    size_t offset = this->node_offset(handle);
    return std::string(this->sequences.begin() + this->offsets[offset], this->sequences.begin() + this->offsets[offset + 1]);
}

// Using undocumented parts of the GBWT interface. --Jouni
bool GBWTGraph::follow_edges(const handle_t& handle, bool go_left, const std::function<bool(const handle_t&)>& iteratee) const {

    // Incoming edges correspond to the outgoing edges of the reverse node.
    gbwt::node_type curr = as_integer(handle);
    if (go_left) {
        curr = gbwt::Node::reverse(curr);
    }

    gbwt::CompressedRecord record = this->index.record(curr);
    for (gbwt::rank_type outrank = 0; outrank < record.outdegree(); outrank++) {
        gbwt::node_type next = record.successor(outrank);
        if (next == gbwt::ENDMARKER) {
            continue;
        }
        // If we started from the reverse node, we must reverse the successor nodes to get
        // the predecessor nodes of the original node.
        if (go_left) {
            next = gbwt::Node::reverse(next);
        }
        if (!iteratee(as_handle(next))) {
            return false;
        }
    }

    return true;
}

void GBWTGraph::for_each_handle(const std::function<bool(const handle_t&)>& iteratee, bool parallel) const {
    if (parallel) {
#pragma omp parallel for schedule(static)
        for (gbwt::node_type node = this->index.firstNode(); node < this->index.sigma(); node += 2) {
            if (!(this->real_nodes[this->node_offset(node) / 2])) {
                continue;
            }
            if (!iteratee(as_handle(node))) {
                // We should stop early but it's not worth the effort.
            }
        }
    } else {
        for (gbwt::node_type node = this->index.firstNode(); node < this->index.sigma(); node += 2) {
            if (!(this->real_nodes[this->node_offset(node) / 2])) {
                continue;
            }
            if (!iteratee(as_handle(node))) {
                return;
            }
        }
    }
}

size_t GBWTGraph::node_size() const {
    return total_nodes;
}

id_t GBWTGraph::min_node_id() const {
    return gbwt::Node::id(this->index.firstNode());
}

id_t GBWTGraph::max_node_id() const {
    id_t next_id = gbwt::Node::id(this->index.sigma());
    return next_id - 1;
}

//------------------------------------------------------------------------------

// Using undocumented parts of the GBWT interface. --Jouni
bool GBWTGraph::follow_edges(gbwt::SearchState state, const std::function<bool(const gbwt::SearchState&)>& iteratee) const {
    gbwt::CompressedRecord record = this->index.record(state.node);
    for (gbwt::rank_type outrank = 0; outrank < record.outdegree(); outrank++) {
        gbwt::node_type next_node = record.successor(outrank);
        if (next_node == gbwt::ENDMARKER) {
            continue;
        }
        gbwt::SearchState next_state(next_node, record.LF(state.range, next_node));
        if (!iteratee(next_state)) {
            return false;
        }
    }

    return true;
}

// Using undocumented parts of the GBWT interface. --Jouni
bool GBWTGraph::follow_edges(gbwt::BidirectionalState state, bool backward, const std::function<bool(const gbwt::BidirectionalState&)>& iteratee) const {
    if (backward) {
        state.flip();
    }

    gbwt::CompressedRecord record = this->index.record(state.forward.node);
    for (gbwt::rank_type outrank = 0; outrank < record.outdegree(); outrank++) {
        gbwt::node_type next_node = record.successor(outrank);
        if (next_node == gbwt::ENDMARKER) {
            continue;
        }
        gbwt::size_type reverse_offset = 0;
        gbwt::BidirectionalState next_state = state;
        next_state.forward.node = next_node;
        next_state.forward.range = record.bdLF(state.forward.range, next_node, reverse_offset);
        next_state.backward.range.first += reverse_offset;
        next_state.backward.range.second = next_state.backward.range.first + next_state.forward.size() - 1;
        if (backward) {
            next_state.flip();
        }
        if (!iteratee(next_state)) {
            return false;
        }
    }

    return true;
}

//------------------------------------------------------------------------------

// Stores haplotype-consistent traversal of the graph and the corresponding sequence.
struct GBWTTraversal {
    // The traversal as a sequence of (begin, length) pairs.
    std::vector<std::pair<pos_t, size_t>> traversal;
    // Length of the traversal.
    size_t length;
    // GBWT search state at the end of the traversal.
    gbwt::SearchState state;

    GBWTTraversal() : length(0) {}

    std::string get_sequence(const HandleGraph& graph) const {
        std::string result;
        result.reserve(this->length);
        for (auto& node : this->traversal) {
            handle_t handle = graph.get_handle(id(node.first), is_rev(node.first));
            result.append(graph.get_sequence(handle), offset(node.first), node.second);
        }
        return result;
    }
};

void extend_traversals(const GBWTGraph& graph,
                       std::stack<GBWTTraversal>& kmers, size_t target_length, size_t minimum_length,
                       const function<void(const std::vector<std::pair<pos_t, size_t>>&, const std::string&)>& lambda) {
    while (!kmers.empty()) {
        GBWTTraversal curr = kmers.top(); kmers.pop();
        // Report the full kmer.
        if (curr.length >= target_length) {
            lambda(curr.traversal, curr.get_sequence(graph));
            continue;
        }

        // Try to extend the kmer to all successor nodes.
        bool extend_success = false;
        graph.follow_edges(curr.state, [&](const gbwt::SearchState& next_state) -> bool {
            if (next_state.empty()) {
                return true;
            }
            size_t node_length = graph.get_length(GBWTGraph::node_to_handle(next_state.node));
            pos_t begin = make_pos_t(gbwt::Node::id(next_state.node), gbwt::Node::is_reverse(next_state.node), 0);
            size_t length = std::min(node_length, target_length - curr.length);
            GBWTTraversal next = curr;
            next.traversal.emplace_back(begin, length);
            next.length += length;
            next.state = next_state;
            kmers.push(next);
            extend_success = true;
            return true;
        });

        // Report sufficiently long kmers that cannot be extended.
        if (!extend_success && curr.length >= minimum_length) {
            lambda(curr.traversal, curr.get_sequence(graph));
        }
    }
}

void for_each_kmer(const GBWTGraph& graph, size_t k,
                   const function<void(const std::vector<std::pair<pos_t, size_t>>&, const std::string&)>& lambda,
                   bool parallel) {

    // Traverse all starting nodes in parallel.
    graph.for_each_handle([&](const handle_t& h) -> bool {
        // Initialize the stack with all starting positions in the current node.
        std::stack<GBWTTraversal> kmers;
        for (bool is_reverse : { false, true }) {
            handle_t handle = (is_reverse ? graph.flip(h) : h);
            gbwt::SearchState state = graph.get_state(handle);
            if (state.empty()) {
                continue;
            }
            id_t id = graph.get_id(handle);
            size_t node_length = graph.get_length(handle);
            for (size_t i = 0; i < node_length; i++) {
                pos_t begin = make_pos_t(id, is_reverse, i);
                size_t length = std::min(node_length - i, k);
                GBWTTraversal kmer;
                kmer.traversal.emplace_back(begin, length);
                kmer.length = length;
                kmer.state = state;
                kmers.push(kmer);
            }
        }

        // Extend with target length and minimum length k.
        extend_traversals(graph, kmers, k, k, lambda);
        return true;
    }, parallel);
}

void for_each_haplotype_window(const GBWTGraph& graph, size_t window_size,
                               const function<void(const std::vector<std::pair<pos_t, size_t>>&, const std::string&)>& lambda,
                               bool parallel) {

    // Traverse all starting nodes in parallel.
    graph.for_each_handle([&](const handle_t& h) -> bool {
        // Initialize the stack with both orientations.
        std::stack<GBWTTraversal> kmers;
        id_t id = graph.get_id(h);
        size_t node_length = graph.get_length(h);
        for (bool is_reverse : { false, true }) {
            handle_t handle = (is_reverse ? graph.flip(h) : h);
            gbwt::SearchState state = graph.get_state(handle);
            if (state.empty()) {
                continue;
            }
            GBWTTraversal kmer;
            kmer.traversal.emplace_back(make_pos_t(id, is_reverse, 0), node_length);
            kmer.length = graph.get_length(handle);
            kmer.state = state;
            kmers.push(kmer);
        }

        // Extend the windows.
        size_t target_length = node_length + window_size - 1;
        extend_traversals(graph, kmers, target_length, window_size, lambda);
        return true;
    }, parallel);
}

void for_each_window(const HandleGraph& graph, size_t window_size,
                     const function<void(const std::vector<std::pair<pos_t, size_t>>&, const std::string&)>& lambda,
                     bool parallel) {

    // Traverse all starting nodes in parallel.
    graph.for_each_handle([&](const handle_t& h) -> bool {
        // Initialize the stack with both orientations.
        std::stack<GBWTTraversal> windows;
        id_t id = graph.get_id(h);
        size_t node_length = graph.get_length(h);
        for (bool is_reverse : { false, true }) {
            handle_t handle = (is_reverse ? graph.flip(h) : h);
            GBWTTraversal window;
            window.traversal.emplace_back(std::make_pair(make_pos_t(id, is_reverse, 0), node_length));
            window.length = graph.get_length(handle);
            windows.push(window);
        }

        // Extend the windows.
        size_t target_length = node_length + window_size - 1;
        while (!windows.empty()) {
            GBWTTraversal window = windows.top(); windows.pop();
            // Report the full window.
            if (window.length >= target_length) {
                lambda(window.traversal, window.get_sequence(graph));
                continue;
            }

            // Extend the window to all successors.
            bool extend_success = false;
            handle_t curr = graph.get_handle(vg::id(window.traversal.back().first), is_rev(window.traversal.back().first));
            graph.follow_edges(curr, false, [&](const handle_t& next) {
                size_t node_length = graph.get_length(next);
                pos_t begin = make_pos_t(graph.get_id(next), graph.get_is_reverse(next), 0);
                size_t length = std::min(node_length, target_length - window.length);
                GBWTTraversal next_window = window;
                next_window.traversal.emplace_back(std::make_pair(begin, length));
                next_window.length += length;
                windows.push(next_window);
                extend_success = true;
            });

            // Report sufficiently long windows that cannot be extended.
            if (!extend_success && window.length >= window_size) {
                lambda(window.traversal, window.get_sequence(graph));
            }
        }

        return true;
    }, parallel);
}

gbwt::GBWT get_gbwt(const std::vector<gbwt::vector_type>& paths) {
    gbwt::size_type node_width = 1, total_length = 0;
    for (auto& path : paths) {
        for (auto node : path) {
            node_width = std::max(node_width, gbwt::bit_length(gbwt::Node::encode(node, true)));
        }
        total_length += 2 * (path.size() + 1);
    }

    gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
    gbwt::GBWTBuilder builder(node_width, total_length);
    for (auto& path : paths) {
        builder.insert(path, true);
    }
    builder.finish();

    std::string filename = temp_file::create("gbwt");
    sdsl::store_to_file(builder.index, filename);
    gbwt::GBWT gbwt_index;
    sdsl::load_from_file(gbwt_index, filename);
    temp_file::remove(filename);

    return gbwt_index;
}

//------------------------------------------------------------------------------

} // namespace vg
