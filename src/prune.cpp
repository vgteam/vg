#include "prune.hpp"

#include <stack>

namespace vg {

constexpr size_t PRUNE_THREAD_BUFFER_SIZE = 128 * 1024;

pair_hash_set<edge_t> find_edges_to_prune(const HandleGraph& graph, size_t k, size_t edge_max) {

    // Each thread collects edges to be deleted into a separate buffer. When the buffer grows
    // large enough, flush it into a shared hash set.
    pair_hash_set<edge_t> result;
    auto flush_buffer = [&result](pair_hash_set<edge_t>& buffer) {
        #pragma omp critical (prune_flush)
        {
            for (const edge_t& edge : buffer) {
                result.insert(edge);
            }
        }
        buffer.clear();
    };

    // for each position on the forward and reverse of the graph
    std::vector<pair_hash_set<edge_t>> buffers(get_thread_count());
    graph.for_each_handle([&](const handle_t& h) {
        // for the forward and reverse of this handle
        // walk k bases from the end, so that any kmer starting on the node will be represented in the tree we build
        for (auto handle_is_rev : { false, true }) {
            handle_t handle = handle_is_rev ? graph.flip(h) : h;
            std::stack<walk_t> walks;
            // for each position in the node, set up a kmer with that start position and the node end or kmer length as the end position
            // determine next positions
            id_t handle_id = graph.get_id(handle);
            size_t handle_length = graph.get_length(handle);
            for (size_t i = 0; i < handle_length; i++) {
                pos_t begin = make_pos_t(handle_id, handle_is_rev, handle_length);
                pos_t end = make_pos_t(handle_id, handle_is_rev, std::min(handle_length, i + k));
                // We are only interested in walks that did not reach length k in the initial node.
                if (offset(end) - offset(begin) < k) {
                    size_t outdegree = graph.get_degree(handle, false);
                    graph.follow_edges(handle, false, [&](const handle_t& next) {
                        if (outdegree > 1 && edge_max == 0) { // our next step takes us over the max
                            int tid = omp_get_thread_num();
                            buffers[tid].insert(graph.edge_handle(handle, next));
                            if (buffers[tid].size() >= PRUNE_THREAD_BUFFER_SIZE) {
                                flush_buffer(buffers[tid]);
                            }
                        } else {
                            walk_t walk(offset(end) - offset(begin), begin, end, next, 0);
                            if (outdegree > 1) {
                                walk.forks++;
                            }
                            walks.push(walk);
                        }
                    });
                }
            }

            // Now expand the kmers until they reach k.
            while (!walks.empty()) {
                walk_t walk = walks.top();
                walks.pop();
                // Did we reach our target length?
                if (walk.length >= k) {
                    continue;
                }
                id_t curr_id = graph.get_id(walk.curr);
                size_t curr_length = graph.get_length(walk.curr);
                bool curr_is_rev = graph.get_is_reverse(walk.curr);
                size_t take = min(curr_length, k - walk.length);
                walk.end = make_pos_t(curr_id, curr_is_rev, take);
                walk.length += take;
                // Do we need to continue to the successor nodes?
                if (walk.length < k) {
                    size_t outdegree = graph.get_degree(walk.curr, false);
                    graph.follow_edges(walk.curr, false, [&](const handle_t& next) {
                        if (outdegree > 1 && edge_max == walk.forks) { // our next step takes us over the max
                            int tid = omp_get_thread_num();
                            buffers[tid].insert(graph.edge_handle(walk.curr, next));
                            if (buffers[tid].size() >= PRUNE_THREAD_BUFFER_SIZE) {
                                flush_buffer(buffers[tid]);
                            }
                        } else {
                            walk_t next_walk = walk;
                            next_walk.curr = next;
                            if (outdegree > 1) {
                                next_walk.forks++;
                            }
                            walks.push(next_walk);
                        }
                    });
                }
            }
        }
    }, true);

    // Flush the buffers and return the result.
    for (pair_hash_set<edge_t>& buffer : buffers) {
        flush_buffer(buffer);
    }
    return result;
}


}
