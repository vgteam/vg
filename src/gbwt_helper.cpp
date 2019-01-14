#include "gbwt_helper.hpp"

#include <stack>

namespace vg {

void for_each_kmer(const HandleGraph& graph, const gbwt::GBWT& haplotypes, size_t k,
                   const function<void(const GBWTTraversal&)>& lambda) {

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
                pos_t end = make_pos_t(id, is_reverse, std::min(seq.length(), i + k));
                GBWTTraversal kmer;
                kmer.traversal.emplace_back(std::make_pair(begin, end));
                kmer.seq = seq.substr(offset(begin), offset(end) - offset(begin));
                kmer.state = state;
                kmers.push(kmer);
            }
        }

        // Process all traversals.
        while (!kmers.empty()) {
            GBWTTraversal curr = kmers.top(); kmers.pop();
            // Report the full kmer.
            if (curr.seq.length() >= k) {
                lambda(curr);
                continue;
            }
            // Try to extend the kmer to all successor nodes.
            // We are using undocumented parts of the GBWT interface --Jouni
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
                pos_t end = make_pos_t(id, is_reverse, std::min(seq.length(), k - curr.seq.length()));
                GBWTTraversal next = curr;
                next.traversal.emplace_back(std::make_pair(begin, end));
                next.seq.append(seq, offset(begin), offset(end) - offset(begin));
                next.state.node = next_node;
                next.state.range = range;
                kmers.push(next);
            }
        }
    }, true);
}

} // namespace vg
