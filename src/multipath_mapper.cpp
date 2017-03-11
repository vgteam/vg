//
//  multipath_mapper.cpp
//  
//
//

#include "multipath_mapper.hpp"

namespace vg {
    MultipathClusterer::MultipathClusterer(const Alignment& alignment,
                                           const vector<MaximalExactMatch>& mems,
                                           const BaseAligner& aligner,
                                           xg::XG& xgindex,
                                           int8_t full_length_bonus,
                                           size_t num_pruning_tracebacks) :
    alignment(alignment), mems(mems), aligner(aligner), full_length_bonus(full_length_bonus)
    {
        
        init_mem_graph(xgindex);
        //    perform_dp();
        //    // TODO: is this the right approach? I might want to limit more than to tracebacks
        //    // overlap edges will almost always get chosen too because I'm assuming maximum score between them
        //    // maybe I need multiple scoring strategies? also worried about backward facing dangling overlaps once
        //    // I am doing tail alignment. Maybe would be best to keep track of overlap edges so I can treat them different
        //    //  - if the non-overlapping edge ends up being reachable in the graph, we probably never want the
        //    //    overlap edge because it's just going to get pulled into taking the same match
        //    prune_to_nodes_on_tracebacks(num_pruning_tracebacks);
        //    remerge_split_nodes
        //    query_node_matches(alignment, xgindex, node_cache);
        //    remove_snarls(alignment, snarl_manager, aligner);
        //    // past this point the edge weights don't really mean anything
        //    query_internal_edge_subgraphs(alignment, xgindex, node_cache);
    }
    
    void MultipathClusterer::init_mem_graph(const xg::XG& xgindex,
                                            size_t max_strand_dist_probes,
                                            size_t max_expected_dist_approx_error) {
        
        // TODO: handle inversions -- do I just need extra paths?
        // for now just not checking for orientation consistency and hoping it works out (but this
        // is only likely when inversion is small so that distance is not overestimated too much)
        auto distance = [&](pos_t& pos_1, pos_t& pos_2) {
            return xgindex.min_approx_path_distance(vector<string>(), id(pos_1), id(pos_2))
            - (int64_t) offset(pos_1) + (int64_t) offset(pos_2);
        };
        
        // the maximum graph distances to the right and left sides detectable from each node
        vector<pair<size_t, size_t>> maximum_detectable_gaps;
        
        // there generally will be at least as many nodes as MEMs, so we can speed up the reallocation
        nodes.reserve(mems.size());
        maximum_detectable_gaps.reserve(mems.size());
        
        for (const MaximalExactMatch& mem : mems) {
            
            // calculate the longest gaps we could detect to the left and right of this MEM
            pair<size_t, size_t> max_gaps(longest_detectable_gap(mem.begin),
                                          longest_detectable_gap(mem.end));
            
            for (gcsa::node_type mem_hit : mem.nodes) {
                nodes.emplace_back(mem.begin, mem.end, make_pos_t(mem_hit), 0, mem_score);
                maximum_detectable_gaps.push_back(max_gaps);
            }
        }
        
        // now we use the distance approximation to cluster the MEM hits according to the strand
        // they fall on using the oriented distance estimation function
        
        // initialize with every hit in its own strand cluster
        vector<vector<size_t>> strand_clusters(nodes.size());
        for (size_t i = 0; i < nodes.size(); i++) {
            strand_clusters[i].insert(i);
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
            while (merged_strand.count(i)) {
                i = merged_strand[i];
            }
            return i;
        }
        
        // make a randomly shuffled list of pairs of nodes to compare distances in
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
        // treating the distances we measured as transitive
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
        for (unordered_map<size_t, int64_t>& relative_pos : strand_relative_position) {
            
            // sort the nodes by relative position
            vector<pair<int64_t, size_t>> sorted_pos;
            for (const pair<size_t, int64_t>& pos_record : relative_pos) {
                sorted_pos.emplace_back(pos_record.second, pos_record.first );
            }
            std::sort(sorted_pos.begin(), sorted_pos.end());
            
            // first sweep: look for MEMs that are colinear along the read and have believable
            // distances between them
            int64_t last_idx = sorted_pos.size() - 1;
            int64_t low = 0, hi = last_idx;
            for (int64_t i = 0; i < sorted_pos.size(); i++) {
                
                int64_t strand_pos = sorted_pos[i].first;
                size_t pivot_idx = sorted_pos[i].second;
                MPCNode& pivot = nodes[pivot_idx];
                int64_t pivot_length = pivot.end - pivot.begin;
                
                int64_t target_low_pos = strand_pos - allowance;
                int64_t target_hi_pos = strand_pos + pivot_length + maximum_detectable_gaps[pivot_idx].second;
                
                while (sorted_pos[low].first < target_low_pos) {
                    low++;
                }
                
                // the upper bound of candidate edges can be in either direction since the maximum dectectable
                // gap varies by read position
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
                    if (next.begin <= pivot.begin || next.end <= pivot.end) {
                        // the nodes cannot be colinear along the read (also filters out j == i)
                        continue;
                    }
                    
                    int64_t graph_dist = max(0, sorted_pos[j].first - strand_pos - pivot_length);
                    
                    // is this distance believeable in the other direction?
                    if (graph_dist > maximum_detectable_gaps[next_idx].first) {
                        continue;
                    }
                    
                    int32_t edge_score;
                    int64_t between_length = next.begin - pivot.end;
                    if (between_length < 0) {
                        // the MEMs overlap, but this can occur in some insertions and deletions
                        // because the SMEM algorithm is "greedy" in taking up as much of the read
                        // as possible
                        // we can check if this happened with the SuffixTree, but it's expensive
                        // so for now we just give it the benefit of the doubt but adjust the edge
                        // score so that the matches don't get double counted
                        
                        int64_t extra_dist = max(0, graph_dist + between_length);
                        
                        edge_score = aligner.match * between_length
                                     + (extra_dist ? -(extra_dist - 1) * aligner.gap_extension - aligner.gap_open : 0);
                    }
                    else if (between_length > graph_dist) {
                        // the read length in between the MEMs is longer than the distance, suggesting a read insert
                        int64_t extra_dist = between_length - graph_dist;
                        
                        edge_score = -aligner.mismatch * graph_dist - (extra_dist - 1) * aligner.gap_extension
                                     - aligner.gap_open;
                    }
                    else if (between_length < graph_dist) {
                        // the read length in between the MEMs is shorter than the distance, suggesting a read deletion
                        int64_t extra_dist = graph_dist - between_length;
                        
                        edge_score = -aligner.mismatch * between_length - (extra_dist - 1) * aligner.gap_extension
                                     - aligner.gap_open;
                    }
                    else {
                        // the read length in between the MEMs is the same as the distance, suggesting a pure mismatch
                        edge_score = -aligner.mismatch * between_length;
                    }
                    
                    pivot.forward_edges_from.emplace_back(next_idx, edge_score);
                    next.forward_edges_to.emplace_back(pivot_idx, edge_score);
                }
            }
            
            // second sweep: look for MEMs that are anti-colinear along the read and have believable
            // distances between them (for the reverse complement alignment)
            // TODO: this code is pretty repetitive, is there an easy way unify it with the first sweep?
            low = 0;
            hi = last_idx;
            for (int64_t i = 0; i < sorted_pos.size(); i++) {
                
                int64_t strand_pos = sorted_pos[i].first;
                size_t pivot_idx = sorted_pos[i].second;
                
                MPCNode& pivot = nodes[pivot_idx];
                
                int64_t target_hi_pos = strand_pos + allowance;
                int64_t target_low_pos = strand_pos - pivot_length - (int64_t) maximum_detectable_gaps[pivot_idx].first;
                
                while (sorted_pos[hi].first > target_hi_pos) {
                    hi--;
                }
                
                // the lower bound of candidate edges can be in either direction since the maximum dectectable
                // gap varies by read position
                if (sorted_pos[low].first > target_low_pos) {
                    while (low == 0 ? false : sorted_pos[low - 1].first <= target_low_pos) {
                        low--;
                    }
                }
                else {
                    while (sorted_pos[low].first < target_low_pos) {
                        low++;
                    }
                }
                
                for (int64_t j = low; j <= hi; j++) {
                    int64_t next_idx = sorted_pos[j].second;
                    MPCNode& next = nodes[next_idx];
                    int64_t next_length = next.end - next.begin;
                    
                    if (next.begin >= pivot.begin || next.end >= pivot.end) {
                        // the nodes cannot be anti-colinear along the read (also filters out j == i)
                        continue;
                    }
                    
                    int64_t graph_dist = max(0, strand_pos - sorted_pos[j].first - next_length);
                    
                    // is this distance believeable in the other direction?
                    if (graph_dist > maximum_detectable_gaps[next_idx].second) {
                        continue;
                    }
                    
                    int32_t edge_score;
                    int64_t between_length = pivot.begin - next.end;
                    if (between_length < 0) {
                        // the MEMs overlap, but this can occur in some insertions and deletions
                        // because the SMEM algorithm is "greedy" in taking up as much of the read
                        // as possible
                        // we can check if this happened with the SuffixTree, but it's expensive
                        // so for now we just give it the benefit of the doubt but adjust the edge
                        // score so that the matches don't get double counted
                        
                        int64_t extra_dist = max(0, graph_dist + between_length);
                        
                        edge_score = aligner.match * between_length
                                     + (extra_dist ? -(extra_dist - 1) * aligner.gap_extension - aligner.gap_open : 0);
                    }
                    else if (between_length > graph_dist) {
                        // the read length in between the MEMs is longer than the distance, suggesting a read insert
                        int64_t extra_dist = between_length - graph_dist;
                        
                        edge_score = -aligner.mismatch * graph_dist - (extra_dist - 1) * aligner.gap_extension
                                     - aligner.gap_open;
                    }
                    else if (between_length < graph_dist) {
                        // the read length in between the MEMs is shorter than the distance, suggesting a read deletion
                        int64_t extra_dist = graph_dist - between_length;
                        
                        edge_score = -aligner.mismatch * between_length - (extra_dist - 1) * aligner.gap_extension
                                     - aligner.gap_open;
                    }
                    else {
                        // the read length in between the MEMs is the same as the distance, suggesting a pure mismatch
                        edge_score = -aligner.mismatch * between_length;
                    }
                    
                    pivot.reverse_edges_from.emplace_back(next_idx, edge_score);
                    next.reverse_edges_to.emplace_back(pivot_idx, edge_score);
                }
            }
        }
    }
    
    inline size_t MultipathMEMAligner::longest_detectable_gap(const string::const_iterator& read_pos) {
        
        // algebraic solution for when score is > 0 assuming perfect match other than gap
        return (min(read_pos - alignment.sequence().begin(), alignment.sequence().end() - read_pos)
                + full_length_bonus - aligner.gap_open) / aligner.gap_extension + 1;
        
    }
}



