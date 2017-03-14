//
//  multipath_mapper.cpp
//  
//
//

#include "multipath_mapper.hpp"

namespace vg {
    
    MultipathAligner::multipath_align(const Alignment& alignment,
                                      const vector<MaximalExactMatch>& mems){
        
        MultipathClusterer clusterer(alignment, mems, aligner, xgindex, full_length_bonus, max_strand_dist_probes,
                                     max_expected_dist_approx_error);
        
        vector<vector<pair<MaximalExactMatch* const, pos_t>>> clusters = clusterer.clusters();
        vector<VG*> cluster_graphs;
        cluster_graphs.reserve(clusters.size());
        
        unordered_map<id_t, size_t> node_id_to_cluster;
        
        for (size_t i = 0; i < clusters.size(); i++) {
            
            vector<pair<MaximalExactMatch* const, pos_t>>& cluster = clusters[i];
            vector<pos_t> positions;
            vector<size_t>& forward_max_dist;
            vector<size_t>& backward_max_dist;
            
            positions.reserve(cluster.size());
            forward_max_dist.reserve(cluster.size());
            backward_max_dist.reserve(cluster.size());
            
            for (auto mem_hit : cluster) {
                positions.push_back(mem_hit.second);
                forward_max_dist.push_back(aligner.longest_detectable_gap(alignment, mem_hit.first->end, full_length_bonus)
                                           + (mem_hit.first->end - mem_hit.first->begin));
                backward_max_dist.push_back(aligner.longest_detectable_gap(alignment, mem_hit.first->begin, full_length_bonus));
            }
            
            Graph graph;
            algorithms::extract_containing_graph(xgindex, graph, positions, forward_max_dist,
                                                 backward_max_dist, &node_cache);
            
            // check if this subgraph overlaps with any previous subgraph (indicates a probable clustering failure where
            // one cluster was split into multiple clusters)
            unordered_set<size_t> overlapping_graphs;
            for (size_t j = 0; j < graph.node_size(); j++) {
                id_t node_id = graph.node(j).id();
                if (node_id_to_cluster.count(node_id)) {
                    overlapping_graphs.insert(node_id_to_cluster[node_id]);
                }
                else {
                    node_id_to_cluster[node_id] = i;
                }
            }
            
            if (overlapping_graphs.empty()) {
                // hacky socketing into the VG constructor so I can use the move constructor and hopefully
                // avoid reallocating the Graph
                bool passed_graph = false;
                auto pass_graph = [&](Graph& into) {
                    if (passed_graph) {
                        return false;
                    }
                    else {
                        into = std::move(graph);
                        passed_graph = true;
                        return true;
                    }
                };
                cluster_graphs.push_back(new VG(pass_graph));
            }
            else if (overlapping_graphs.size() == 1) {
                // this subgraph overlaps with one other graph, so we merge it in rather than make a new one
                cluster_graphs[*overlapping_graphs.begin()]->extend(graph);
            }
            else {
                // this subgraph chains together multiple clusters, so we merge all of them into the cluster
                // with the smallest index
                
                size_t min_idx_cluster = std::numeric_limits<size_t>::max();
                for (size_t j : overlapping_graphs) {
                    if (j < min_idx_cluster) {
                        min_idx_cluster = j;
                    }
                }
                
                // merge in the new graph
                cluster_graphs[min_idx_cluster]->extend(graph);
                
                // merge in all the other graphs it connects to and remove them from the list
                overlapping_graphs.erase(min_idx_cluster);
                for (size_t j : overlapping_graphs) {
                    std::swap(cluster_graphs[j], cluster_graphs.back());
                    cluster_graphs[min_idx_cluster]->extend(cluster_graphs.back());
                    delete cluster_graphs.back();
                    cluster_graphs.pop_back();
                }
            }
        }
        
        
        for (VG* vg : cluster_graphs) {
            delete vg;
        }
    }
    
    
    MultipathClusterer::MultipathClusterer(const Alignment& alignment,
                                           const vector<MaximalExactMatch>& mems,
                                           const BaseAligner& aligner,
                                           xg::XG& xgindex,
                                           int8_t full_length_bonus,
                                           size_t max_strand_dist_probes,
                                           size_t max_expected_dist_approx_error) : aligner(aligner) {
        
        // the maximum graph distances to the right and left sides detectable from each node
        vector<pair<size_t, size_t>> maximum_detectable_gaps;
        
        // there generally will be at least as many nodes as MEMs, so we can speed up the reallocation
        nodes.reserve(mems.size());
        maximum_detectable_gaps.reserve(mems.size());
        
        for (const MaximalExactMatch& mem : mems) {
            
            // calculate the longest gaps we could detect to the left and right of this MEM
            pair<size_t, size_t> max_gaps(aligner.longest_detectable_gap(alignment, mem.begin, full_length_bonus),
                                          aligner.longest_detectable_gap(alignment, mem.end, full_length_bonus));
            
            for (gcsa::node_type mem_hit : mem.nodes) {
                nodes.emplace_back(mem, make_pos_t(mem_hit), mem_score);
                maximum_detectable_gaps.push_back(max_gaps);
            }
        }
        
        // now we use the distance approximation to cluster the MEM hits according to the strand
        // they fall on using the oriented distance estimation function
        
        // initialize with every hit in its own strand cluster
        vector<vector<size_t>> strand_clusters(nodes.size());
        for (size_t i = 0; i < nodes.size(); i++) {
            strand_clusters[i].push_back(i);
        }
        
        // for recording the distance of any pair that we check with a finite distance
        unordered_map<pair<size_t, size_t>, int64_t> recorded_finite_dists;
        // for recording the number of times elements of a strand cluster have been compared
        // and found an infinite distance
        map<pair<size_t, size_t>, size_t> num_infinite_dists;
        
        // when we merge strand clusters, we add a pointer to the strand they merged into
        // essentially, this become trees with pointers upward that we can traverse to
        // find the unmerged strand cluster any hit currently belongs to
        unordered_map<size_t, size_t> merged_strand;
        auto get_merged_strand = [&merged_strand](size_t i) {
            vector<size_t> path;
            while (merged_strand.count(i)) {
                path.push_back(i);
                i = merged_strand[i];
            }
            // path compression (similar to union find)
            for (size_t j : path) {
                merged_strand[j] = i;
            }
            return i;
        }
        
        // make a randomly shuffled list of pairs of nodes to compare distances in
        // TODO: is there a way to do this without generating all pairs?
        vector<pair<size_t, size_t>> comparisons;
        comparisons.reserve((nodes.size() * (nodes.size() - 1)) / 2);
        for (int64_t i = 1; i < nodes.size(); i++) {
            for (int64_t j = 0; j < i; j++) {
                comparisons.emplace_back(j, i):
            }
        }
        std::default_random_engine gen(std::random_device());
        std::shuffle(comparisons.begin(), comparisons.end(), gen);
        
        for (pair<size_t, size_t>& node_pair : comparisons) {
            
            // TODO: I don't love this component, it makes the algorithm O(n^2 log n) in number of MEM hits
            size_t strand_1 = get_merged_strand(node_pair.first);
            size_t strand_2 = get_merged_strand(node_pair.second);

            if (strand_1 == strand_2) {
                // these are already identified as on the same strand, don't need to do it again
                continue;
            }
            
            if (num_infinite_dists[make_pair(strand_1, strand_2)] >= max_strand_dist_probes) {
                // we've already checked multiple distances between these strand clusters and
                // none have returned a finite distance, so we conclude that they are in fact
                // on separate clusters and decline to check any more distances
                continue;
            }
            
            pos_t& pos_1 = nodes[node_pair.first];
            pos_t& pos_2 = nodes[node_pair.second];
            
            int64_t oriented_dist = xg_index.closest_shared_path_oriented_distance(id(pos_1), offset(pos_1), is_rev(pos_1),
                                                                                   id(pos_2), offset(pos_2), is_rev(pos_2));
            
            if (oriented_dist == std::numeric_limits<int64_t>::max()) {
                // distance is estimated at infinity, so these are either on different strands
                // or the path heuristic failed to find a shared path
                
                num_infinite_dists[make_pair(strand_1, strand_2)]++;
                num_infinite_dists[make_pair(strand_2, strand_1)]++;
            }
            else {
                // the distance is finite, so merge the strand clusters
                
                recorded_finite_dists[node_pair] = oriented_dist;
                
                // add the smaller cluster onto the larger one to minimize copies
                vector<size_t>* smaller_clust, larger_clust;
                
                // ensure that the first node in the pair belongs to the larger cluster
                if (strand_clusters[node_pair.first].size() < strand_clusters[node_pair.second].size()) {
                    std::swap(node_pair.first, node_pair.second);
                    std::swap(strand_1, strand_2);
                }
                
                auto& clust_1 = strand_clusters[node_pair.first];
                auto& clust_2 = strand_clusters[node_pair.second];
                
                clust_1.insert(clust_1.end(), clust_2.begin(), clust_2.end());
                
                // choose one of the strand clusters at random to remove (makes the merged_strand
                // tree have expected height O(log n))
                size_t strand_retaining, strand_removing;
                if (gen() % 2) {
                    merged_strand[node_pair.second] = node_pair.first;
                    std::swap(clust_2, strand_clusters.back());
                    strand_removing = strand_2;
                    strand_retaining = strand_1;
                }
                else {
                    merged_strand[node_pair.first] = node_pair.second;
                    clust_2 = std::move(clust_1);
                    std::swap(clust_1, strand_clusters.back());
                    strand_removing = strand_1;
                    strand_retaining = strand_2;
                }
                strand_clusters.pop_back();
                
                // collect the number of times this strand cluster has had an infinite distance to other strand
                vector<pair<size_t, size_t>> inf_counts;
                auto end = num_infinite_dists.upper_bound(make_pair(strand_removing, numeric_limits<size_t>::max()));
                auto begin = num_infinite_dists.lower_bound(make_pair(strand_removing, 0));
                for (auto iter = begin; iter != end; iter++) {
                    inf_counts.push_back((*iter).first.second, (*iter).second);
                }
                
                // add these counts to the other strand cluster and remove this one from the map
                for (const pair<size_t, size_t> inf_count : inf_counts) {
                    num_infinite_dists[make_pair(strand_retaining, inf_count.first)] += inf_count.second;
                    num_infinite_dists.erase(make_pair(strand_removing, inf_count.first));
                    num_infinite_dists.erase(make_pair(inf_count.first, strand_removing));
                }
            }
        }
        
        // build the graph of relative distances in adjacency list representation
        // by construction each strand cluster will be an undirected, unrooted tree
        unordered_map<size_t, vector<size_t>> strand_distance_tree;
        for (const auto& dist_record : recorded_finite_dists) {
            strand_distance_tree[dist_record.first.first].push_back(dist_record.first.second);
            strand_distance_tree[dist_record.first.second].push_back(dist_record.first.first);
        }
        
        // now approximate the relative positions along the strand by traversing each tree and
        // treating the distances we estimated as transitive
        vector<unordered_map<size_t, int64_t>> strand_relative_position;
        vector<bool> processed(nodes.size());
        for (const auto& adjacency_record : strand_distance_tree) {
            if (processed[adjacency_record.first]) {
                continue;
            }
            
            size_t strand = get_merged_strand(adjacency_record.first);
            unordered_map<size_t, int64_t>& relative_pos = strand_relative_position[strand];
            
            // arbitrarily make this node the 0 point
            relative_pos[adjacency_record.first] = 0;
            processed[adjacency_record.first] = true;
            
            // traverse the strand's tree with DFS
            list<size_t> queue{adjacency_record.first};
            while (!queue.empty()) {
                size_t curr = queue.back();
                queue.pop_back();
                
                int64_t curr_pos = relative_pos[curr];
                
                for (size_t next : strand_distance_tree[curr]) {
                    if (processed[next]) {
                        continue;
                    }
                    
                    // invert the sign of the distance if we originally measured it in the other order
                    int64_t dist = recorded_finite_dists.count(make_pair(curr, next)) ?
                                   recorded_finite_dists[make_pair(curr, next)] :
                                   -recorded_finite_dists[make_pair(next, furr)];
                    
                    // find the position relative to the previous node we just traversed
                    relative_pos[next] = curr_pos + dist;
                    processed[next] = true;
                    
                    queue.push_back(next);
                }
            }
        }
        
        // now we use the strand clusters and the estimated distances to make the DAG for the
        // approximate MEM alignment
        
        int64_t allowance = max_expected_dist_approx_error;
        for (const unordered_map<size_t, int64_t>& relative_pos : strand_relative_position) {
            
            // sort the nodes by relative position
            vector<pair<int64_t, size_t>> sorted_pos;
            for (const pair<size_t, int64_t>& pos_record : relative_pos) {
                sorted_pos.emplace_back(pos_record.second, pos_record.first);
            }
            std::sort(sorted_pos.begin(), sorted_pos.end());
            
            // find edges within each strand cluster by first identifying the interval of MEMs that meets
            // the graph distance constrant for each MEM and then checking for read colinearity and the
            // reverse distance constraint
            int64_t last_idx = sorted_pos.size() - 1;
            int64_t low = 0, hi = last_idx;
            for (int64_t i = 0; i < sorted_pos.size(); i++) {
                
                int64_t strand_pos = sorted_pos[i].first;
                size_t pivot_idx = sorted_pos[i].second;
                MPCNode& pivot = nodes[pivot_idx];
                int64_t pivot_length = pivot.mem->end - pivot.mem->begin;
                
                // the limits of how far away we might detect edges to add to the clustering graph
                int64_t target_low_pos = strand_pos - allowance;
                int64_t target_hi_pos = strand_pos + pivot_length + maximum_detectable_gaps[pivot_idx].second + allowance;
                
                // move the lower boundary of the search interval to the lowest value inside the
                // the target interval
                while (sorted_pos[low].first < target_low_pos) {
                    low++;
                }
                
                // move the upper boundary of the search interval to the highest value inside the
                // the target interval (this one can move in either direction because the maximum
                // detectable gap changes)
                if (sorted_pos[hi].first > target_hi_pos) {
                    while (sorted_pos[hi].first > target_hi_pos) {
                        hi--;
                    }
                }
                else {
                    while (hi == last_idx ? false : sorted_pos[hi + 1].first <= target_hi_pos) {
                        hi++;
                    }
                }
                
                for (int64_t j = low; j <= hi; j++) {
                    int64_t next_idx = sorted_pos[j].second;
                    MPCNode& next = nodes[next_idx];
                    
                    if (next.mem->begin <= pivot.mem->begin || next.mem->end <= pivot.mem->end) {
                        // the MEMs cannot be colinear along the read (also filters out j == i)
                        continue;
                    }
                    
                    // the length of the sequence in between the MEMs (can be negative if they overlap)
                    int64_t between_length = next.mem->begin - pivot.mem->end;
                    // the estimated distance between the end of the pivot and the start of the next MEM in the graph
                    int64_t graph_dist = max(0, sorted_pos[j].first - strand_pos - pivot_length);
                    // the discrepancy between the graph distance and the read distance
                    int64_t gap_length = abs(graph_dist - between_length);
                    
                    if (gap_length > maximum_detectable_gaps[next_idx].first + allowance) {
                        // the gap between the MEMs is too long to be believable from the next node
                        continue;
                    }
                    
                    int32_t edge_score;
                    if (between_length < 0) {
                        // the MEMs overlap, but this can occur in some insertions and deletions
                        // because the SMEM algorithm is "greedy" in taking up as much of the read
                        // as possible
                        // we can check if this happened with the SuffixTree, but it's expensive
                        // so for now we just give it the benefit of the doubt but adjust the edge
                        // score so that the matches don't get double counted
                        
                        int64_t extra_dist = max(0, gap_length);
                        
                        edge_score = -aligner.match * between_length
                                     + (extra_dist ? -(extra_dist - 1) * aligner.gap_extension - aligner.gap_open : 0);
                    }
                    else if (between_length > graph_dist) {
                        // the read length in between the MEMs is longer than the distance, suggesting a read insert
                        edge_score = -aligner.mismatch * graph_dist - (gap_length - 1) * aligner.gap_extension
                                     - aligner.gap_open;
                    }
                    else if (between_length < graph_dist) {
                        // the read length in between the MEMs is shorter than the distance, suggesting a read deletion
                        edge_score = -aligner.mismatch * between_length - (gap_length - 1) * aligner.gap_extension
                                     - aligner.gap_open;
                    }
                    else {
                        // the read length in between the MEMs is the same as the distance, suggesting a pure mismatch
                        edge_score = -aligner.mismatch * between_length;
                    }
                    
                    // add the edges in
                    pivot.edges_from.emplace_back(next_idx, edge_score);
                    next.edges_to.emplace_back(pivot_idx, edge_score);
                }
            }
        }
    }

    void MultipathClusterer::topological_order(vector<size_t>& order_out) {
        
        // initialize return value
        order_out.clear();
        order_out.resize(nodes.size());
        size_t order_idx = nodes.size() - 1;
        
        // initialize iteration structures
        vector<bool> enqueued = vector<bool>(nodes.size());
        vector<size_t> edge_index = vector<size_t>(nodes.size());
        vector<size_t> stack;
        
        // iterate through starting nodes
        for (size_t init_node_idx = 0; init_node_idx < nodes.size(); init_node_idx++) {
            if (enqueued[init_node_idx]) {
                continue;
            }
            // navigate through graph with DFS
            stack.push_back(init_node_idx);
            enqueued[init_node_idx] = true;
            while (!stack.empty()) {
                size_t node_idx = stack.back();
                if (edge_index[node_idx] < nodes[node_idx].edges_from.size()) {
                    int64_t target_idx = nodes[node_idx].edges_from[node_idx].to_idx;
                    if (enqueued[target_idx]) {
                        edge_index[node_idx]++;
                    }
                    else {
                        stack.push_back(target_idx);
                        enqueued[target_idx] = true;
                    }
                }
                else {
                    // add to topological order in reverse finishing order
                    stack.pop_back();
                    order_out[order_idx] = node_idx;
                    order_idx--;
                }
            }
        }
    }
    
    void MultipathMEMAligner::identify_sources_and_sinks(vector<size_t>& sources_out,
                                                         vector<size_t>& sinks_out) {
        
        sources_out.clear();
        sinks_out.clear();
        
        vector<bool> is_source(nodes.size(), true);
        
        for (size_t i = 0; i < nodes.size(); i++) {
            if (nodes[i].edges_from.empty()) {
                sinks_out.push_back(i);
            }
            
            for (MultipathMEMEdge& edge : nodes[i].edges_from) {
                is_source[edge.to_idx] = false;
            }
        }
        
        for (size_t i = 0; i < nodes.size(); i++) {
            if (is_source[i]) {
                sources_out.push_back(i);
            }
        }
    }
    
    
    
    void MultipathClusterer::connected_components(vector<vector<size_t>>& components_out) {
        
        components_out.clear();
        vector<bool> enqueued(nodes.size());
        
        // check each node in turn to find new components
        for (size_t dfs_start_idx = 0; dfs_start_idx < nodes.size(); dfs_start_idx++) {
            if (enqueued[dfs_start_idx]) {
                // we've already found this node from some component
                continue;
            }
            
            // this node belongs to a component we haven't found yet, use DFS to find the rest
            vector<size_t> stack {dfs_start_idx};
            enqueued[dfs_start_idx] = true;
            components_out.emplace_back(1, dfs_start_idx);
            
            while (!stack.empty()) {
                
                MultipathMEMNode& node = nodes[stack.back()];
                stack.pop_back();
                
                // search in both forward and backward directions
                
                for (MultipathMEMEdge& edge : node.edges_from) {
                    
                    if (!enqueued[edge.to_idx]) {
                        stack.push_back(edge.to_idx);
                        enqueued[edge.to_idx] = true;
                        components_out.back().push_back(edge.to_idx);
                    }
                }
                
                for (MultipathMEMEdge& edge : node.edges_to) {
                    
                    if (!enqueued[edge.to_idx]) {
                        stack.push_back(edge.to_idx);
                        enqueued[edge.to_idx] = true;
                        components_out.back().push_back(edge.to_idx);
                    }
                }
            }
        }
    }
    
    void MultipathClusterer::perform_dp() {
        
        // as in local alignment, minimum score is the score of node itself
        for (size_t i = 0; i < nodes.size(); i++) {
            nodes[i].dp_score = nodes[i].score;
        }
        
        vector<size_t> order;
        topological_order(order);
        
        for (size_t i : order) {
            MPCNode& node = nodes[i];
            
            // for each edge out of this node
            for (MPCEdge& edge : node.edges_from) {
                
                // check if the path through the node out of this edge increase score of target node
                MPCNode& target_node = nodes[edge.to_idx];
                int32_t extend_score = node.score + edge.weight + target_node.score;
                if (extend_score > target_node.dp_score) {
                    target_node.dp_score = extend_score;
                }
            }
        }
    }
    
    int64_t MultipathClusterer::trace_read_coverage(const vector<size_t>& trace) {
        
        MPCNode& node = nodes[trace[0]];
        auto curr_begin = node.mem->begin;
        auto curr_end = node.mem->end;
        
        int64_t total = 0;
        for (size_t i = 1; i < trace.size(); i++) {
            node = nodes[trace[i]];
            if (node.mem->begin >= curr_end) {
                total += (curr_end - curr_begin);
                curr_begin = node.mem->begin;
            }
            curr_end = node.mem->end;
        }
        return total + (curr_end - curr_begin);
    }
    
    vector<vector<pair<MaximalExactMatch* const, pos_t>>> MultipathClusterer::clusters(int32_t max_qual_score) {
        
        vector<vector<MaximalExactMatch* const, pos_t>> to_return;
        if (nodes.size() == 0) {
            // this should only happen if we have filtered out all MEMs, so there are none to cluster
            return to_return;
        }
        
        perform_dp();
        
        // find the weakly connected components, which should correspond to mappings
        vector<vector<size_t>> components;
        connected_components(components);
        
        // find the node with the highest DP score in each connected component
        vector<pair<int32_t, size_t>> component_traceback_ends(components.size(),
                                                               pair<int32_t, size_t>(numeric_limits<int32_t>::min(), 0));
        for (size_t i = 0; i < components.size(); i++) {
            vector<size_t>& component = components[i];
            pair<int32_t, size_t>& traceback_end = component_traceback_ends[i];
            for (size_t j = 0; j < component.size(); j++) {
                int32_t dp_score = nodes[component[i]].dp_score
                if (dp_score > traceback_end.first) {
                    traceback_end.first = dp_score;
                    traceback_end.second = j;
                }
            }
        }
        
        // sort indices in descending order by their highest traceback score
        vector<size_t> order = range_vector(0, components.size());
        std::sort(order.begin(), order.end() [&](const size_t i, const size_t j) {
            return component_traceback_ends[i].first > component_traceback_ends[j].first;
        });
        
        int32_t top_score = component_traceback_ends[order[0]].first;
        
        for (size_t i : order) {
            // get the component and the traceback end
            vector<size_t>& component = components[i];
            size_t trace_idx = component[component_traceback_ends[i].second];
            
            // if this cluster does not look like it even affect the mapping quality of the top scoring
            // cluster, don't bother forming it
            // TODO: this approximation could break down sometimes, need to look into it
            // TODO: is there a way to make the aligner do this? I don't like having this functionality outside of it
            if (4.3429448190325183 * aligner.log_base * (top_score - traceback_end.first) > max_qual_score ) {
                continue;
            }
            
            // traceback until hitting a node that has its own score (indicates beginning of a local alignment)
            vector<size_t> trace{trace_idx};
            while (nodes[trace_idx].dp_score > nodes[trace_idx].score) {
                for (MPCEdge& edge : nodes[trace_idx].edges_to) {
                    if (nodes[edge.to_idx].dp_score + edge.weight == nodes[trace_idx].dp_score) {
                        trace_idx = edge.to_idx;
                        trace.push_back(trace_idx);
                    }
                }
            }
            
            // make a cluster
            to_return.emplace_back();
            auto& cluster = to_return.back();
            for (auto iter = trace.rbegin(); iter != trace.rend(); iter++) {
                MPCNode& node = nodes[*iter];
                cluster.emplace_back(node.mem, node.start_pos);
            }
        }
        
        return to_return;
    }

}



