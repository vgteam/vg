/**
 * \file k_widest_paths.cpp
 *
 * Implementation of Yen's Algorithm over the bidirected graph.
 */

#include "k_widest_paths.hpp"

#include <structures/updateable_priority_queue.hpp>

namespace vg {
namespace algorithms {

using namespace structures;

//#define debug_vg_algorithms

pair<double, vector<handle_t>> widest_dijkstra(const HandleGraph* g, handle_t source, handle_t sink,
                                               function<double(const handle_t&)> node_weight_callback,
                                               function<double(const edge_t&)> edge_weight_callback,
                                               function<bool(const handle_t&)> is_node_ignored_callback,
                                               function<bool(const edge_t&)> is_edge_ignored_callback) {

    // We keep a priority queue so we can visit the handle with the shortest
    // distance next. We put handles in here whenever we see them with shorter
    // distances (since STL priority queue can't update), so we also need to
    // make sure nodes coming out haven't been visited already.
    // (score, previous, current)
    using Record = tuple<double, handle_t, handle_t>;
    
    // We filter out handles that have already been visited.  And keep track of predecessors
    unordered_map<handle_t, pair<handle_t, double>> visited;

    // We need a custom ordering for the queue
    struct IsFirstGreater {
        inline bool operator()(const Record& a, const Record& b) {
            return get<0>(a) < get<0>(b);
        }
    };
    
    // We use a filtered priority queue for auto-Dijkstra
    UpdateablePriorityQueue<Record, handle_t, vector<Record>, IsFirstGreater> queue([](const Record& item) {
            return get<2>(item);
    });
    
    // We keep a current handle
    handle_t current;
    handle_t previous;
    // we don't include the score of the source
    // todo: should this be an option?
    double score = numeric_limits<double>::max();
    queue.push(make_tuple(score, source, source));
    
    while (!queue.empty()) {
        // While there are things in the queue, get the first.
        tie(score, previous, current) = queue.top();
        queue.pop();

#ifdef debug_vg_algorithms
        cerr << "Visit " << g->get_id(current) << " " << g->get_is_reverse(current) << " at width " << score << endl;
#endif    

        // Remember that we made it here.
        if (!visited.count(current) || score > visited[current].second) {
            visited[current] = make_pair(previous, score);
        }
                
        if (current != sink) {
            g->follow_edges(current, false, [&](const handle_t& next) {
                    // For each handle to the right of here
                    if (!visited.count(next) && 
                        !is_node_ignored_callback(next) &&
                        !is_edge_ignored_callback(g->edge_handle(current, next))) {

                        double next_score = score;

                        if (next != source && next != sink) {
                            // we don't include the source / sink
                            // todo: should we? should it be an option?
                            next_score = min(next_score, node_weight_callback(next));
                        }
                        next_score = min(next_score, edge_weight_callback(g->edge_handle(current, next)));

                        // New shortest distance. Will never happen after the handle comes out of the queue because of Dijkstra.
                        queue.push(make_tuple(next_score, current, next));
                        
#ifdef debug_vg_algorithms
                        cerr << "\tNew best path to " << g->get_id(next) << ":" << g->get_is_reverse(next)
                             << " at width " << next_score << endl;
#endif
                        
                    } else {
#ifdef debug_vg_algorithms
                        cerr << "\tDisregard path to " << g->get_id(next) << ":" << g->get_is_reverse(next)
                             << " at width " << score << " due to " << visited.count(next) << " || "
                             << is_node_ignored_callback(next) << " || " << is_edge_ignored_callback(g->edge_handle(current, next))
                             <<endl;
#endif
                    }
                });
        }
    }

    // trace our results back out of the visited table
    vector<handle_t> widest_path;
    double width = 0;
    if (visited.count(sink)) {
        width = visited[sink].second;
        for (handle_t tb = sink; widest_path.empty() || widest_path.back() != source; tb = visited[tb].first) {
            widest_path.push_back(tb);
        }
        std::reverse(widest_path.begin(), widest_path.end());
    }

#ifdef debug_vg_algorithms
    cerr << "Returning wideset path (w=" << width << "): ";
    for (auto h : widest_path) {
        cerr << g->get_id(h) << ":" << g->get_is_reverse(h) << ", ";
    }
    cerr << endl;
#endif
    
    return make_pair(width, widest_path);
}

// https://en.wikipedia.org/wiki/Yen%27s_algorithm
vector<pair<double, vector<handle_t>>> yens_k_widest_paths(const HandleGraph* g, handle_t source, handle_t sink,
                                                           size_t K,
                                                           function<double(const handle_t&)> node_weight_callback,
                                                           function<double(const edge_t&)> edge_weight_callback) {

    vector<pair<double, vector<handle_t>>> best_paths;
    best_paths.reserve(K);

    // get the widest path from dijkstra
    best_paths.push_back(widest_dijkstra(g, source, sink, node_weight_callback,
                                         edge_weight_callback, [](handle_t) {return false;},
                                         [](edge_t) {return false;}));
    
    // working path set, mapped to score
    // todo: make more efficient
    set<vector<handle_t>> B;
    // used to pull out the biggest element in B
    multimap<double, set<vector<handle_t>>::iterator> score_to_B;
    
    // start scanning for our k-1 next-widest paths
    for (size_t k = 1; k < K; ++k) {

        // we look for a "spur node" in the previous path.  the current path will be the previous path
        // up to that spur node, then a new path to the sink. (i is the index of the spur node in
        // the previous (k - 1) path
        vector<handle_t>& prev_path = best_paths[k - 1].second;
        for (size_t i = 0; i < prev_path.size() - 1; ++i) {
            
            handle_t spur_node = prev_path[i];
            // root path = prev_path[0 : i]

#ifdef debug_vg_algorithms
            cerr << "k=" << k << ": spur node=" << g->get_id(spur_node) << ":" << g->get_is_reverse(spur_node) << endl;
#endif
            unordered_set<edge_t> forgotten_edges;
            for (const auto& p_v : best_paths) {
                const vector<handle_t>& p = p_v.second;

                // check if the root path is a prefix of p
                bool is_common_root = true;
                for (size_t j = 0; j <= i && is_common_root; ++j) {
                    if (j >= prev_path.size() || p[j] != prev_path[j]) {
                        is_common_root = false;
                    }
                }

                // remove the links that are part of the previous shortest paths which share the same root path
                if (is_common_root) {
#ifdef debug_vg_algorithms
                  cerr << "forgetting edge " << g->get_id(p[i]) << ":" << g->get_is_reverse(p[i]) << " -- "
                       << g->get_id(p[i+1]) << ":" << g->get_is_reverse(p[i+1]) << endl;
#endif
                    forgotten_edges.insert(g->edge_handle(p[i], p[i+1]));
                }
            }
            
            // forget the root path too (except spur node)
            unordered_set<handle_t> forgotten_nodes;
            for (int j = 0; j < (int)i - 1; ++j) {
#ifdef debug_vg_algorithms
              cerr << "forgetting node " << g->get_id(prev_path[j]) << ":" << g->get_is_reverse(prev_path[j]) << endl;
#endif
                forgotten_nodes.insert(prev_path[j]);
            }

            // find our path from the the spur_node to the sink
            pair<double, vector<handle_t>> spur_path_v = widest_dijkstra(g, spur_node, sink, node_weight_callback, edge_weight_callback,
                                                                         [&](handle_t h) {return forgotten_nodes.count(h);},
                                                                         [&](edge_t e) {return forgotten_edges.count(e);});

            if (!spur_path_v.second.empty()) {
            
                // make the path by combining the root path and the spur path
                pair<double, vector<handle_t>> total_path;
                double total_width = numeric_limits<double>::max();
                for (size_t j = 0; j < i; ++j) {
                    total_path.second.push_back(prev_path[j]);
                    total_width = min(total_width, node_weight_callback(prev_path[j]));
                    if (j > 0) {
                        total_width = min(total_width, edge_weight_callback(g->edge_handle(prev_path[j-1], prev_path[j])));
                    }
                }
                if (!total_path.second.empty()) {
                    total_width = min(total_width, edge_weight_callback(g->edge_handle(total_path.second.back(), spur_path_v.second.front())));
                }

                // insert the path into our sorted set
                total_path.second.insert(total_path.second.end(), spur_path_v.second.begin(), spur_path_v.second.end());
                total_path.first = min(total_width, spur_path_v.first);
                pair<set<vector<handle_t>>::iterator, bool> ins = B.insert(total_path.second);
                if (ins.second == true) {
                    score_to_B.insert(make_pair(total_path.first, ins.first));
                } // todo: is there any reason we'd need to update the score of an existing entry in B?
            }
        }

        if (B.empty()) {
            break;
        }

        assert(score_to_B.size() == B.size());
        multimap<double, set<vector<handle_t>>::iterator>::iterator best_B_it = std::prev(score_to_B.end());
        
        // the best path gets put into our output list
        best_paths.push_back(make_pair(best_B_it->first, *best_B_it->second));
        B.erase(best_B_it->second);
        score_to_B.erase(best_B_it);

#ifdef debug_vg_algorithms
        cerr << "adding best path (w=" << best_paths.back().first << "): ";
        for (auto h : best_paths.back().second) {
            cerr << g->get_id(h) << ":" << g->get_is_reverse(h) << ", ";
        }
        cerr << endl;
#endif

    }

    return best_paths;
}


}
}
