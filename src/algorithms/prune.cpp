#include "algorithms/prune.hpp"
#include "hash_map.hpp"
#include "position.hpp"
#include "source_sink_overlay.hpp"

#include <stack>

namespace vg {
namespace algorithms {

/// Record a <=k-length walk in the context of a graph.
struct walk_t {
    walk_t(uint16_t l,
           const pos_t& b,
           const pos_t& e,
           const handle_t& c,
           uint16_t f)
    : length(l), begin(b), end(e), curr(c), forks(f) { };
    /// our start position
    pos_t begin;
    pos_t end; /// one past the (current) end of the kmer
    handle_t curr; /// the next handle we extend into
    uint16_t forks; /// how many branching edge crossings we took to get here
    uint16_t length; /// how far we've been
};

constexpr size_t PRUNE_THREAD_BUFFER_SIZE = 1024 * 1024;

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

size_t prune_complex(DeletableHandleGraph& graph,
                     int path_length, int edge_max) {
    
    auto edges_to_destroy = find_edges_to_prune(graph, path_length, edge_max);
    for (auto& edge : edges_to_destroy) {
        graph.destroy_edge(edge);
    }
    return edges_to_destroy.size();
}

size_t prune_complex_with_head_tail(DeletableHandleGraph& graph,
                                    int path_length, int edge_max) {
    
    SourceSinkOverlay source_sink_graph(&graph, path_length);
    
    auto edges_to_destroy = find_edges_to_prune(source_sink_graph,
                                                path_length,
                                                edge_max);
    
    for (auto& edge : edges_to_destroy) {
        auto ss_handle_1 = source_sink_graph.forward(edge.first);
        auto ss_handle_2 = source_sink_graph.forward(edge.second);
        if (ss_handle_1 != source_sink_graph.get_source_handle()
            && ss_handle_1 != source_sink_graph.get_sink_handle()
            && ss_handle_2 != source_sink_graph.get_source_handle()
            && ss_handle_2 != source_sink_graph.get_sink_handle()) {
            // this is not an edge involving the artificial source/sink nodes
            graph.destroy_edge(source_sink_graph.get_underlying_handle(edge.first),
                               source_sink_graph.get_underlying_handle(edge.second));
            
        }
    }
    return edges_to_destroy.size();
}

size_t prune_short_subgraphs(DeletableHandleGraph& graph, int min_size) {
    
    unordered_set<handle_t> to_destroy;
    
    // DFS from all tips
    for (auto tip : handlealgs::find_tips(&graph)) {
        //cerr << "begin trav from " << graph.get_id(tip) << " " << graph.get_is_reverse(tip) << endl;
        auto start = graph.forward(tip);
        if (to_destroy.count(start)) {
            // we already found this subgraph from another tip
            //cerr << "skipping" << endl;
            continue;
        }
        vector<handle_t> stack(1, start);
        unordered_set<handle_t> seen{start};
        int size_seen = 0;
        // stop when we've seen a large enough subgraph
        while (!stack.empty() && size_seen < min_size) {
            handle_t handle = stack.back();
            stack.pop_back();
            size_seen += graph.get_length(handle);
            //cerr << "destack " << graph.get_id(handle) << ", update size seen to " << size_seen << endl;
            for (bool go_left : {true, false}) {
                graph.follow_edges(handle, go_left, [&](const handle_t& next) {
                    handle_t fwd_next = graph.forward(next);
                    if (!seen.count(fwd_next)) {
                        //cerr << "stack up " << graph.get_id(fwd_next) << ", update size seen to " << size_seen << endl;
                        stack.push_back(fwd_next);
                        seen.insert(fwd_next);
                    }
                });
            }
        }
        if (size_seen < min_size) {
            //cerr << "component is small enough to destroy" << endl;
            // this component is under the size limit, mark them for destruction
            for (auto handle : seen) {
                to_destroy.insert(handle);
            }
        }
    }
    
    // destroy all handles that we marked
    for (auto handle : to_destroy) {
        graph.destroy_handle(handle);
    }
    
    return to_destroy.size();
}


size_t remove_high_degree_nodes(DeletableHandleGraph& g, int max_degree) {
    vector<handle_t> to_remove;
    g.for_each_handle([&](const handle_t& h) {
        int edge_count = 0;
        g.follow_edges(h, false, [&](const handle_t& ignored) {
            ++edge_count;
        });
        g.follow_edges(h, true, [&](const handle_t& ignored) {
            ++edge_count;
        });
        if (edge_count > max_degree) {
            to_remove.push_back(h);
        }
    });
    // now destroy the high degree nodes
    for (auto& h : to_remove) {
        g.destroy_handle(h);
    }
    return to_remove.size();
}

}
}
