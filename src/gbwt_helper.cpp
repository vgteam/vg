#include "gbwt_helper.hpp"

#include <stack>

namespace vg {

// Stores haplotype-consistent traversal of the graph and the corresponding sequence.
struct GBWTTraversal {
    // The traversal as a sequence of (begin, length) pairs.
    std::vector<std::pair<pos_t, size_t>> traversal;
    // The sequence.
    std::string seq;
    // GBWT search state at the end of the traversal.
    gbwt::SearchState state;
};

void extend_traversals(const HandleGraph& graph, const gbwt::GBWT& haplotypes,
                       std::stack<GBWTTraversal>& kmers, size_t target_length, size_t minimum_length,
                       const function<void(const std::vector<std::pair<pos_t, size_t>>&, const std::string&)>& lambda) {
    while (!kmers.empty()) {
        GBWTTraversal curr = kmers.top(); kmers.pop();
        // Report the full kmer.
        if (curr.seq.length() >= target_length) {
            lambda(curr.traversal, curr.seq);
            continue;
        }

        // Try to extend the kmer to all successor nodes.
        // We are using undocumented parts of the GBWT interface --Jouni
        bool extend_success = false;
        gbwt::CompressedRecord record = haplotypes.record(curr.state.node);
        for (gbwt::rank_type outrank = 0; outrank < record.outdegree(); outrank++) {
            gbwt::node_type next_node = record.successor(outrank);
            if (next_node == gbwt::ENDMARKER) {
                continue;
            }
            gbwt::range_type range = record.LF(curr.state.range, next_node);
            if (gbwt::Range::empty(range)) {
                continue;
            }
            id_t id = gbwt::Node::id(next_node);
            bool is_reverse = gbwt::Node::is_reverse(next_node);
            std::string seq = graph.get_sequence(gbwt_to_handle(graph, next_node));
            pos_t begin = make_pos_t(id, is_reverse, 0);
            size_t length = std::min(seq.length(), target_length - curr.seq.length());
            GBWTTraversal next = curr;
            next.traversal.emplace_back(std::make_pair(begin, length));
            next.seq.append(seq, offset(begin), length);
            next.state.node = next_node;
            next.state.range = range;
            kmers.push(next);
            extend_success = true;
        }

        // Report sufficiently long kmers that cannot be extended.
        if (!extend_success && curr.seq.length() >= minimum_length) {
            lambda(curr.traversal, curr.seq);
        }
    }
}

void for_each_kmer(const HandleGraph& graph, const gbwt::GBWT& haplotypes, size_t k,
                   const function<void(const std::vector<std::pair<pos_t, size_t>>&, const std::string&)>& lambda,
                   bool parallel) {

    // Traverse all starting nodes in parallel.
    graph.for_each_handle([&](const handle_t& h) {
        // Initialize the stack with all starting positions in the current node.
        std::stack<GBWTTraversal> kmers;
        for (bool is_reverse : { false, true }) {
            handle_t handle = (is_reverse ? graph.flip(h) : h);
            gbwt::SearchState state = haplotypes.find(handle_to_gbwt(graph, handle));
            if (state.empty()) {
                continue;
            }
            id_t id = graph.get_id(handle);
            std::string seq = graph.get_sequence(handle);
            for (size_t i = 0; i < seq.length(); i++) {
                pos_t begin = make_pos_t(id, is_reverse, i);
                size_t length = std::min(seq.length() - i, k);
                GBWTTraversal kmer;
                kmer.traversal.emplace_back(std::make_pair(begin, length));
                kmer.seq = seq.substr(offset(begin), length);
                kmer.state = state;
                kmers.push(kmer);
            }
        }

        // Extend with target length and minimum length k.
        extend_traversals(graph, haplotypes, kmers, k, k, lambda);
    }, parallel);
}

void for_each_window(const HandleGraph& graph, const gbwt::GBWT& haplotypes, size_t window_size,
                    const function<void(const std::vector<std::pair<pos_t, size_t>>&, const std::string&)>& lambda,
                    bool parallel) {

    // Traverse all starting nodes in parallel.
    graph.for_each_handle([&](const handle_t& h) {
        // Initialize the stack with both orientations.
        std::stack<GBWTTraversal> kmers;
        id_t id = graph.get_id(h);
        size_t node_length = graph.get_length(h);
        for (bool is_reverse : { false, true }) {
            handle_t handle = (is_reverse ? graph.flip(h) : h);
            gbwt::SearchState state = haplotypes.find(handle_to_gbwt(graph, handle));
            if (state.empty()) {
                continue;
            }
            GBWTTraversal kmer;
            kmer.traversal.emplace_back(std::make_pair(make_pos_t(id, is_reverse, 0), node_length));
            kmer.seq = graph.get_sequence(handle);
            kmer.state = state;
            kmers.push(kmer);
        }

        // Extend the windows.
        extend_traversals(graph, haplotypes, kmers, node_length + window_size - 1, window_size, lambda);
    }, parallel);
}

void for_each_window(const HandleGraph& graph, size_t window_size,
                     const function<void(const std::vector<std::pair<pos_t, size_t>>&, const std::string&)>& lambda,
                     bool parallel) {

    // Traverse all starting nodes in parallel.
    graph.for_each_handle([&](const handle_t& h) {
        // Initialize the stack with both orientations.
        std::stack<GBWTTraversal> windows;
        id_t id = graph.get_id(h);
        size_t node_length = graph.get_length(h);
        for (bool is_reverse : { false, true }) {
            handle_t handle = (is_reverse ? graph.flip(h) : h);
            GBWTTraversal window;
            window.traversal.emplace_back(std::make_pair(make_pos_t(id, is_reverse, 0), node_length));
            window.seq = graph.get_sequence(handle);
            windows.push(window);
        }

        // Extend the windows.
        size_t target_length = node_length + window_size - 1;
        while (!windows.empty()) {
            GBWTTraversal window = windows.top(); windows.pop();
            // Report the full window.
            if (window.seq.length() >= target_length) {
                lambda(window.traversal, window.seq);
                continue;
            }

            // Extend the window to all successors.
            bool extend_success = false;
            handle_t curr = graph.get_handle(vg::id(window.traversal.back().first), is_rev(window.traversal.back().first));
            graph.follow_edges(curr, false, [&](const handle_t& next) {
                std::string seq = graph.get_sequence(next);
                pos_t begin = make_pos_t(graph.get_id(next), graph.get_is_reverse(next), 0);
                size_t length = std::min(seq.length(), target_length - window.seq.length());
                GBWTTraversal next_window = window;
                next_window.traversal.emplace_back(std::make_pair(begin, length));
                next_window.seq.append(seq, offset(begin), length);
                windows.push(next_window);
                extend_success = true;
            });

            // Report sufficiently long windows that cannot be extended.
            if (!extend_success && window.seq.length() >= window_size) {
                lambda(window.traversal, window.seq);
            }
        }
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

} // namespace vg
