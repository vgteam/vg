#include "eades_algorithm.hpp"

namespace vg {
namespace algorithms {

using namespace std;

    vector<handle_t> eades_algorithm(const HandleGraph* graph) {
        
        // decide which strand will be "forward" for each node
        vector<handle_t> canonical_orientation = single_stranded_orientation(graph);
        if (canonical_orientation.size() < graph->node_size()) {
            cerr << "error:[eades_algorithm] Eades' algorithm only valid on graphs with a single stranded orientation" << endl;
            exit(1);
        }
        
        // maps handles to records of ((in-degree, out-degree), bucket position)
        unordered_map<handle_t, pair<pair<int64_t, int64_t>, list<handle_t>::iterator>> degree_info;
        
        // buckets based on delta(u) among non-source, non-sink nodes (see paper)
        // buckets are numbered -n + 3, -n + 2, ... , n - 2, n - 3
        auto assign_bucket = [&](const int64_t& in_degree, const int64_t& out_degree) {
            return out_degree - in_degree + int64_t(canonical_orientation.size()) - 3;
        };
        vector<list<handle_t>> delta_buckets(2 * canonical_orientation.size() - 5);
        
        vector<handle_t> sources;
        vector<handle_t> sinks;
        
        for (const handle_t& handle : canonical_orientation) {
            
            // compute in- and out-degree
            int64_t in_degree = 0;
            int64_t out_degree = 0;
            graph->follow_edges(handle, true, [&](const handle_t& prev) {
                in_degree++;
            });
            graph->follow_edges(handle, false, [&](const handle_t& next) {
                out_degree++;
            });
            
            if (in_degree == 0) {
                // source
                sources.emplace_back(handle);
            }
            else if (out_degree == 0) {
                // sink
                sinks.emplace_back(handle);
            }
            else {
                // non-source, non-sink
                auto& bucket = delta_buckets[assign_bucket(in_degree, out_degree)];
                bucket.emplace_front(handle);
                degree_info[handle] = make_pair(make_pair(in_degree, out_degree), bucket.begin());
            }
        }
        
        // identify the highest non-empty bucket
        int64_t max_delta_bucket = delta_buckets.size() - 1;
        while (max_delta_bucket >= 0) {
            if (!delta_buckets[max_delta_bucket].empty()) {
                break;
            }
            max_delta_bucket--;
        }
        
        // init the layout to fill
        vector<handle_t> layout(canonical_orientation.size());
        
        // the next positions to add to in the layout (we fill from both sides)
        int64_t next_left_idx = 0;
        int64_t next_right_idx = layout.size() - 1;
        
        // update data structures to remove an edge into a node
        function<void(const handle_t&)> remove_inward_edge = [&](const handle_t& next) {
            auto iter = degree_info.find(next);
            if (iter != degree_info.end()) {
                // this node is in a delta bucket
                
                // remove it from the current bucket
                auto& degrees = iter->second.first;
                int64_t current_bucket_num = assign_bucket(degrees.first, degrees.second);
                auto& bucket = delta_buckets[current_bucket_num];
                bucket.erase(iter->second.second);
                
                // update the degrees to remove the inward edge
                degrees.first--;
                if (degrees.first == 0) {
                    // this is now a source
                    degree_info.erase(next);
                    sources.push_back(next);
                }
                else {
                    // this moves up one bucket
                    current_bucket_num++;
                    auto& new_bucket = delta_buckets[current_bucket_num];
                    new_bucket.emplace_front(next);
                    iter->second.second = new_bucket.begin();
                    
                    // if necessary, identify this as the new highest delta bucket
                    max_delta_bucket = max(max_delta_bucket, current_bucket_num);
                }
            }
        };
        
        // update data structures to remove an edge out of a node
        function<void(const handle_t&)> remove_outward_edge = [&](const handle_t& prev) {
            auto iter = degree_info.find(prev);
            if (iter != degree_info.end()) {
                // this node is in a delta bucket
                
                // remove it from the current bucket
                auto& degrees = iter->second.first;
                int64_t current_bucket_num = assign_bucket(degrees.first, degrees.second);
                auto& bucket = delta_buckets[current_bucket_num];
                bucket.erase(iter->second.second);
                
                // update the degrees to remove the outward edge
                degrees.second--;
                if (degrees.second == 0) {
                    // this is now a sink
                    degree_info.erase(prev);
                    sinks.push_back(prev);
                }
                else {
                    // this moves down one bucket
                    current_bucket_num--;
                    auto& new_bucket = delta_buckets[current_bucket_num];
                    new_bucket.emplace_front(prev);
                    iter->second.second = new_bucket.begin();
                }
            }
        };
        
        while (next_left_idx <= next_right_idx) {
            
            if (!sources.empty()) {
                // add it to the layout
                handle_t source = sources.back();
                sources.pop_back();
                layout[next_left_idx] = source;
                next_left_idx++;
                
                // remove a source from the graph
                graph->follow_edges(source, false, remove_inward_edge);
            }
            else if (!sinks.empty()) {
                // add it to the layout
                handle_t sink = sinks.back();
                sinks.pop_back();
                layout[next_right_idx] = sink;
                next_right_idx--;
                
                // remove a sink from the graph
                graph->follow_edges(sink, true, remove_outward_edge);
            }
            else {
                // remove a node in the highest delta bucket from the graph
                auto& bucket = delta_buckets[max_delta_bucket];
                handle_t next = bucket.back();
                bucket.pop_back();
                degree_info.erase(next);
                
                // add it to the layout
                layout[next_left_idx] = next;
                next_left_idx++;
                
                graph->follow_edges(next, false, remove_inward_edge);
                graph->follow_edges(next, true, remove_outward_edge);
            }
            
            // move the max bucket to the left if it has been emptied
            while (max_delta_bucket >= 0) {
                if (!delta_buckets[max_delta_bucket].empty()) {
                    break;
                }
                max_delta_bucket--;
            }
        }
        
        return layout;
    }
}
}
