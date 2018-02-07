#include "prune.hpp"

namespace vg {

vector<edge_t> find_edges_to_prune(const HandleGraph& graph, size_t k, size_t edge_max) {
    // for each position on the forward and reverse of the graph
    //unordered_set<edge_t> edges_to_prune;
    vector<vector<edge_t> > edges_to_prune;
    edges_to_prune.resize(get_thread_count());
    graph.for_each_handle([&](const handle_t& h) {
            // for the forward and reverse of this handle
            // walk k bases from the end, so that any kmer starting on the node will be represented in the tree we build
            for (auto handle_is_rev : { false, true }) {
                //cerr << "###########################################" << endl;
                handle_t handle = handle_is_rev ? graph.flip(h) : h;
                list<walk_t> walks;
                // for each position in the node, set up a kmer with that start position and the node end or kmer length as the end position
                // determine next positions
                id_t handle_id = graph.get_id(handle);
                size_t handle_length = graph.get_length(handle);
                string handle_seq = graph.get_sequence(handle);
                for (size_t i = 0; i < handle_length;  ++i) {
                    pos_t begin = make_pos_t(handle_id, handle_is_rev, handle_length);
                    pos_t end = make_pos_t(handle_id, handle_is_rev, min(handle_length, i+k));
                    walk_t walk = walk_t(offset(end)-offset(begin), begin, end, handle, 0);
                    if (walk.length < k) {
                        // are we branching over more than one edge?
                        size_t next_count = 0;
                        graph.follow_edges(walk.curr, false, [&](const handle_t& next) { ++next_count; });
                        graph.follow_edges(walk.curr, false, [&](const handle_t& next) {
                                auto e = make_pair(walk.curr, next);
                                if (next_count > 1 && edge_max == walk.forks) { // our next step takes us over the max
                                    int tid = omp_get_thread_num();
                                    edges_to_prune[tid].push_back(make_pair(walk.curr, next));
                                } else {
                                    walks.push_back(walk);
                                    auto& todo = walks.back();
                                    todo.curr = next;
                                    if (next_count > 1) {
                                        ++todo.forks;
                                    }
                                }
                            });
                    } else {
                        walks.push_back(walk);
                    }
                }
                // now expand the kmers until they reach k
                while (!walks.empty()) {
                    // first we check which ones have reached length k in the current handle; for each of these we run lambda and remove them from our list
                    auto walks_end = walks.end();
                    for (list<walk_t>::iterator q = walks.begin(); q != walks_end; ++q) {
                        auto& walk = *q;
                        // did we reach our target length?
                        if (walk.length >= k) {
                            q = walks.erase(q);
                        } else {
                            id_t curr_id = graph.get_id(walk.curr);
                            size_t curr_length = graph.get_length(walk.curr);
                            bool curr_is_rev = graph.get_is_reverse(walk.curr);
                            size_t take = min(curr_length, k-walk.length);
                            walk.end = make_pos_t(curr_id, curr_is_rev, take);
                            walk.length += take;
                            if (walk.length < k) {
                                // if not, we need to expand through the node then follow on
                                size_t next_count = 0;
                                graph.follow_edges(walk.curr, false, [&](const handle_t& next) { ++next_count; });
                                graph.follow_edges(walk.curr, false, [&](const handle_t& next) {
                                        auto e = make_pair(walk.curr, next);
                                        if (next_count > 1 && edge_max == walk.forks) { // our next step takes us over the max
                                            int tid = omp_get_thread_num();
                                            edges_to_prune[tid].push_back(e);
                                        } else {
                                            walks.push_back(walk);
                                            auto& todo = walks.back();
                                            todo.curr = next;
                                            if (next_count > 1) {
                                                ++todo.forks;
                                            }
                                        }
                                    });
                                q = walks.erase(q);
                            } else {
                                // nothing, we'll remove it next time around
                            }
                        }
                    }
                }
            }
        }, true);
    uint64_t total_edges = 0;
    for (auto& v : edges_to_prune) total_edges += v.size();
    vector<edge_t> merged; merged.reserve(total_edges);
    for (auto& v : edges_to_prune) {
        merged.insert(merged.end(), v.begin(), v.end());
    }
    // duplicates are assumed to be dealt with externally
    return merged;
}


}
