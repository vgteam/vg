//
//  multipath_mem_aligner.cpp
//

#include "multipath_mem_aligner.hpp"

namespace vg {
    
    vector<size_t> range_vector(size_t begin, size_t end) {
        size_t len = end - begin;
        vector<size_t> range(len, begin);
        for (size_t i = 1; i < len; i++) {
            range[i] = begin + i;
        }
        return range;
    }
    
    MultipathMEMAligner::MultipathMEMAligner(const Alignment& alignment,
                                             const vector<MaximalExactMatch>& mems,
                                             const QualAdjAligner& aligner,
                                             xg::XG& xgindex,
                                             LRUCache<id_t, Node>& node_cache,
                                             SnarlManager& snarl_manager,
                                             int8_t full_length_bonus,
                                             size_t num_pruning_tracebacks) {
        
        init_first_pass_graph(alignment, mems, aligner, xgindex, full_length_bonus);
        perform_dp();
        // TODO: is this the right approach? I might want to limit more than to tracebacks
        // overlap edges will almost always get chosen too because I'm assuming maximum score between them
        // maybe I need multiple scoring strategies? also worried about backward facing dangling overlaps once
        // I am doing tail alignment. Maybe would be best to keep track of overlap edges so I can treat them different
        //  - if the non-overlapping edge ends up being reachable in the graph, we probably never want the
        //    overlap edge because it's just going to get pulled into taking the same match
        prune_to_nodes_on_tracebacks(num_pruning_tracebacks);
        remerge_split_nodes
        query_node_matches(alignment, xgindex, node_cache);
        remove_snarls(alignment, snarl_manager, aligner);
        // past this point the edge weights don't really mean anything
        query_internal_edge_subgraphs(alignment, xgindex, node_cache);
    }
    
    void MultipathMEMAligner::init_first_pass_graph(const Alignment& alignment,
                                                    const vector<MaximalExactMatch>& mems,
                                                    const QualAdjAligner& aligner,
                                                    const xg::XG& xgindex,
                                                    int8_t full_length_bonus) {
    
        // TODO: handle inversions -- do I just need extra paths?
        // for now just not checking for orientation consistency and hoping it works out (but this
        // is only likely when inversion is small so that distance is not overestimated too much)
        auto distance = [&](pos_t& pos_1, pos_t& pos_2) {
            return xgindex.min_approx_path_distance(vector<string>(), id(pos_1), id(pos_2))
                   - (int64_t) offset(pos_1) + (int64_t) offset(pos_2);
        };
        
        // the maximum graph distances to the right and left sides detectable from each node
        vector<pair<size_t, size_t>> maximum_detectable_gaps;
        
        int8_t match_score = aligner.adjusted_score_matrix[25 * aligner.max_qual_score];
        
        // there generally will be at least as many nodes as MEMs, so we can speed up the reallocation
        nodes.reserve(mems.size());
        maximum_detectable_gaps.reserve(mems.size());
        
        for (const MaximalExactMatch& mem : mems) {
            
            string::const_iterator base_qual_begin = alignment.quality().begin()
                                                     + (mem.begin - alignment.sequence().begin());
            
            int32_t mem_score = aligner.score_exact_match(mem.begin, mem.end, base_qual_begin);
            
            // calculate the longest gaps we could detect to the right and left of this MEM
            pair<size_t, size_t> max_gaps(longest_detectable_gap(alignment,
                                                                 mem,
                                                                 mem_score,
                                                                 match_score,
                                                                 true,
                                                                 aligner,
                                                                 full_length_bonus),
                                          longest_detectable_gap(alignment,
                                                                 mem,
                                                                 mem_score,
                                                                 match_score,
                                                                 false,
                                                                 aligner,
                                                                 full_length_bonus));
            
            for (gcsa::node_type mem_hit : mem.nodes) {
                nodes.emplace_back(mem.begin, mem.end, make_pos_t(mem_hit), 0, mem_score);
                maximum_detectable_gaps.push_back(max_gaps);
            }
        }
        
        // first pass: identify all pairs that could conceivably have edges, i.e. that have a
        // minimum insert of at most the maximum detectable gap and whos start positions are
        // sequential along the read. we do not enforce that the MEMs have no overlap, because
        // later we will check if we can get an edge by removing an overlap
        unordered_map<pair<int64_t, int64_t>, int64_t> pair_distances;
        for (int64_t i = 0, last = nodes.size() - 1; i < last; i++) {
            for (int64_t j = i + 1; j < nodes.size(); j++) {
                
                MultipathMEMNode& node_1 = nodes[i];
                MultipathMEMNode& node_2 = nodes[j];
                
                // if one MEM is contained in the other along the read, there is no way to make an edge
                // from the end of one to the interior of the other by chopping off an overlap
                if ((node_2.begin >= node_1.begin && node_2.end <= node_1.end)
                    || (node_1.begin >= node_2.begin && node_1.end <= node_2.end)) {
                    continue;
                }
                
                if (node_1.begin <= node_2.begin) {
                    int64_t dist = distance(node_1.start_pos, node_2.start_pos) - (node_1.end - node_1.begin);
                    // note: if the start of node_2 is contained inside the node_1 match, this will
                    // usually create a negative distance, which will pass the filter
                    
                    // the minimum of the max detectable gap is the max gap here
                    int64_t max_gap_length = (int64_t) min(maximum_detectable_gaps[i].first,
                                                           maximum_detectable_gaps[j].second);
                    
                    // record all pairs that are under the minimum
                    if (node_2.begin - node_1.end < max_gap_length && dist < max_gap_length) {
                        pair_distances[make_pair(i, j)] = dist;
                    }
                }
                else {
                    int64_t dist = distance(node_2.start_pos, node_1.start_pos) - (node_2.end - node_2.begin);
                    // note: if the start of node_1 is contained inside the node_2 match, this will
                    // usually create a negative distance, which will pass the filter
                    
                    // the minimum of the max detectable gap is the max gap here
                    int64_t max_gap_length = (int64_t) min(maximum_detectable_gaps[j].first,
                                                           maximum_detectable_gaps[i].second);
                    
                    // record all pairs that are under the minimum
                    if (node_1.begin - node_2.end < max_gap_length && dist < max_gap_length) {
                        pair_distances[make_pair(j, i)] = dist;
                    }
                }
            }
        }
        
        // as we make edges, record any overlaps that are necessary
        vector<pair<pair<int64_t, int64_t>, int64_t>> overlaps;
        
        for (const pair<pair<int64_t, int64_t>, size_t>& pair_distance : pair_distances) {
            
            MultipathMEMNode& node_1 = nodes[pair_distance.first.first];
            MultipathMEMNode& node_2 = nodes[pair_distance.first.second];
            
            // TODO: need to look up how approximate distance function handles strand, sign, etc.
            int64_t absolute_dist = abs(pair_distance.second);
            
            if (node_1.end <= node_2.begin) {
                // these matches are non-overlapping on the read
                
                // the amount of read between these MEMs
                int64_t between_length = node_2.begin - node_1.end;
                
                // TODO: assuming matches for intervening sequence will make it difficult to prune
                // off erroneous sub-MEMs and overlaps that we aren't going to filter out at this stage
                
                // find the maximum possible score (goal is sensitivity at this stage)
                int32_t weight;
                if (between_length == absolute_dist) {
                    // can account for gap with only mismatches
                    weight = between_length * match_score;
                }
                else if (between_length < absolute_dist) {
                    // need a deletion to account for gap
                    int64_t remaining = absolute_dist - between_length;
                    weight = between_length * match_score - aligner.scaled_gap_open - (remaining - 1) * aligner.scaled_gap_extension;
                }
                else {
                    // need an insertion to account for gap
                    int64_t extra = between_length - absolute_dist;
                    weight = absolute_dist * match_score - aligner.scaled_gap_open - (extra - 1) * aligner.scaled_gap_extension;
                }
                node_1.edges_from.emplace_back(pair_distance.first.second, weight);
            }
            else {
                // these matches overlap on the read, but maybe we can make an edge by trimming
                // off an overlap
                
                // figure out how much of the first MEM could potentially overlap with the second
                // TODO: some MEMs could represent the same read sequence, might be worth finding a
                // way to only build suffix tree once per sequence
                SuffixTree suffix_tree(node_1.begin, node_1.end);
                size_t overlap = suffix_tree.longest_overlap(node_2.begin, node_2.end);
                
                // record any non-trivial overlaps
                if (overlap && overlap != node_2.end - node_2.begin) {
                    overlaps.push_back(make_pair(pair_distance.first, overlap));
                }
            }
        }
        
        // sort overlaps in order of decreasing length (makes splitting nodes easier)
        sort(overlaps.begin(), overlaps.end(), [](pair<pair<int64_t, int64_t>, size_t> overlap_1,
                                                  pair<pair<int64_t, int64_t>, size_t> overlap_2) {
            return overlap_1.second > overlap_2.second;
        });
        
        for (const pair<pair<int64_t, int64_t>, size_t>& pair_overlap : overlaps) {
            
            MultipathMEMNode& from_node = nodes[pair_overlap.first.first];
            MultipathMEMNode& onto_node = nodes[pair_overlap.first.second];
            
            string::const_iterator split_point = onto_node.begin + pair_overlap.second;
            string::const_iterator base_qual_begin = alignment.quality().begin()
                                                     + (onto_node.begin - alignment.sequence().begin());
            
            // the score of the part of the match that overlaps
            int32_t overlap_score = aligner.score_exact_match(onto_node.begin, split_point,
                                                              base_qual_begin);
            
            // move the section of the node after the overlap onto a new node
            nodes.emplace_back(split_point, onto_node.end, onto_node.start_pos,
                               pair_overlap.second, onto_node.score - overlap_score);
            
            MultipathMEMNode& new_node = nodes.back();
            
            // move the edge list onto the new node
            new_node.edges_from = std::move(onto_node.edges_from);
            // the old edge list is now in an undefined state, so clear it to be sure it's empty
            onto_node.edges_from.clear();
            // add a 0 weight edge between old node and the new one
            onto_node.edges_from.emplace_back(nodes.size() - 1, 0);
            onto_node.end = split_point;
            onto_node.score = overlap_score;
            
            // add the edge from the other node onto the suffix node
            
            // we're guaranteed that this will be non-negative because the longest overlap
            // must at least contain any part of the read where they overlap
            int64_t between_length = new_node.begin - from_node.end;
            
            int64_t absolute_dist = abs(pair_distances[pair_overlap.first] - pair_overlap.second);
            
            // find the maximum possible score (goal is sensitivity at this stage)
            int32_t weight;
            if (between_length == absolute_dist) {
                // can account for gap with only mismatches
                weight = between_length * match_score;
            }
            else if (between_length < absolute_dist) {
                // need a deletion to account for gap
                int64_t remaining = absolute_dist - between_length;
                weight = between_length * match_score - aligner.scaled_gap_open - (remaining - 1) * aligner.scaled_gap_extension;
            }
            else {
                // need an insertion to account for gap
                int64_t extra = between_length - absolute_dist;
                weight = absolute_dist * match_score - aligner.scaled_gap_open - (extra - 1) * aligner.scaled_gap_extension;
            }
            
            // add an edge into the overlap
            from_node.edges_from.emplace_back(nodes.size() - 1, weight);
        }
        
        // forward edge lists are now stably on the same node, so construct the backward edge lists
        for (size_t i = 0; i < nodes.size(); i++) {
            // for each edge on the node
            for (MultipathMEMEdge& edge : nodes[i].edges_from) {
                // copy the edge in reverse
                nodes[edge.to_idx].edges_to.emplace_back(i, edge.weight);
            }
        }
        
        // TODO: remove transitively redundant edges (node-skipping edges)?
    }
    
    inline size_t MultipathMEMAligner::longest_detectable_gap(const Alignment& alignment,
                                                              const MaximalExactMatch& mem,
                                                              const int32_t mem_score,
                                                              const int8_t match_score,
                                                              const bool gap_to_right,
                                                              const QualAdjAligner& aligner,
                                                              const int8_t full_length_bonus) {
        
        // length of read on the same side of the gap as the MEM
        size_t same_side_remaining = gap_to_right ? mem.begin - alignment.sequence().begin()
                                                  : alignment.sequence().end() - mem.end;
        // length of read on the opposite side of the gap as the MEM
        size_t other_side_remaining = gap_to_right ? alignment.sequence().end() - mem.end
                                                   : mem.begin - alignment.sequence().begin();
        
        // algebraic solution for when score is > 0 assuming perfect match other than gap
        return min((same_side_remaining * match_score + mem_score + full_length_bonus
                    - aligner.scaled_gap_open) / aligner.scaled_gap_extension,
                   (other_side_remaining * match_score + full_length_bonus
                    - aligner.scaled_gap_open) / aligner.scaled_gap_extension) + 1;
        
    }
    
    inline size_t MultipathMEMAligner::longest_detectable_gap(const Alignment& alignment,
                                                              const string::const_iterator& gap_begin,
                                                              const int8_t match_score,
                                                              const QualAdjAligner& aligner,
                                                              const int8_t full_length_bonus) {
        
        // algebraic solution for when score is > 0 assuming perfect match other than gap
        return (min(gap_begin - alignment.sequence().begin(), alignment.sequence().end() - gap_begin) * match_score
                - aligner.scaled_gap_open) / aligner.scaled_gap_extension + 1;
        
    }
    
    void MultipathMEMAligner::topological_order(vector<size_t>& order_out) {
        
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
    
    void MultipathMEMAligner::perform_dp() {
        
        // as in local alignment, minimum score is the score of node itself
        for (size_t i = 0; i < nodes.size(); i++) {
            nodes[i].dp_score = nodes[i].score;
        }
        
        vector<size_t> order;
        topological_order(order);
        
        for (size_t i : order) {
            MultipathMEMNode& node = nodes[i];
            
            // for each edge out of this node
            for (MultipathMEMEdge& edge : node.edges_from) {
                
                // check if the path through the node out of this edge increase score of target node
                MultipathMEMNode& target_node = nodes[edge.to_idx];
                int32_t extend_score = node.score + edge.weight + target_node.score;
                if (extend_score > target_node.dp_score) {
                    target_node.dp_score = extend_score;
                }
            }
        }
    }
    
    void MultipathMEMAligner::prune_to_nodes_on_tracebacks(size_t num_tracebacks) {
        
        // initialize list of edges that are on tracebacks
        unordered_set<pair<size_t, size_t>> edge_idxs;
        
        // initialize traceback manager
        TracebackManager traceback_manager(nodes, num_tracebacks);
        
        // traceback manager will find the top scoring tracebacks in the course of
        // tracing back and keep track of them
        for (; traceback_manager.has_next(); traceback_manager.next()) {
            // where does the traceback start?
            size_t traceback_idx = traceback_manager.get_traceback_start();
            traceback_manager.mark_traced(traceback_idx);
            MultipathMEMNode& node = nodes[traceback_idx];
            
            // what is the traceback's score?
            int32_t traceback_score = traceback_manager.curr_traceback_score();
            
            // search backwards until reaching a node that has its initial score
            while (node.dp_score != node.score) {
                size_t prev_idx = traceback_idx;
                
                // does the traceback we are on deflect from the optimal traceback here?
                if (traceback_manager.at_next_deflection(traceback_idx)) {
                    // take the deflection
                    traceback_idx = traceback_manager.deflect();
                }
                else {
                    // take the optimal traceback
                    int32_t backward_score = node.dp_score - node.score;
                    MultipathMEMEdge* trace_edge = nullptr;
                    for (MultipathMEMEdge& edge : nodes[traceback_idx].edges_to) {
                        int32_t score_diff = backward_score - (nodes[edge.to_idx].dp_score + edge.weight);
                        if (score_diff == 0 && !trace_edge) {
                            // this is an optimal traceback
                            trace_edge = &edge;
                        }
                        else {
                            // this traceback is non optimal, but it might be part of a high scoring
                            // sub-optimal traceback
                            // TODO: detect whether this is an overlap edge and add it to the traced edges list
                            // without actually including a basically identical traceback?
                            traceback_manager.propose_deflection(traceback_idx, edge.to_idx,
                                                                 traceback_score - score_diff);
                        }
                    }
                    if (!trace_edge) {
                        cerr << "error:[MultipathMEMAligner] traceback error in MEM aligner" << endl;
                        exit(1);
                    }
                    
                    // follow the edge backwards
                    traceback_idx = trace_edge->to_idx;
                }
                
                // mark the edge as being involved in a traceback
                edge_idxs.insert(make_pair(traceback_idx, prev_idx));
                // get the next node
                traceback_manager.mark_traced(traceback_idx);
                node = nodes[traceback_idx];
            }
        }
        
        // prune back to the nodes and edges from these tracebacks
        prune_graph(traceback_manager.traced_node_idxs, edge_idxs);
    }
    
    void MultipathMEMAligner::connected_components(vector<vector<size_t>>& components_out) {
        
        components_out.clear();
        vector<bool> enqueued(nodes.size());
        
        for (size_t dfs_start_idx = 0; dfs_start_idx < nodes.size(); dfs_start_idx++) {
            if (enqueued[dfs_start_idx]) {
                continue;
            }
            
            vector<size_t> stack {dfs_start_idx};
            enqueued[dfs_start_idx] = true;
            components_out.emplace_back(1, dfs_start_idx);
            
            while (!stack.empty()) {
                
                MultipathMEMNode& node = nodes[stack.back()];
                stack.pop_back();
                
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
    
    void MultipathMEMAligner::prune_graph(const unordered_set<size_t>& node_idxs,
                                          const unordered_set<pair<size_t, size_t>>& edge_idxs) {
        
        // the number of nodes removed before this index
        vector<size_t> removed_before;
        removed_before.reserve(nodes.size());
        
        size_t nodes_removed_so_far = 0;
        
        for (size_t i = 0; i < nodes.size(); i++) {
            // keep track of how many nodes have been removed up to this point
            removed_before.push_back(nodes_removed_so_far);
            
            if (node_idxs.count(i)) {
                // keep this node
                
                vector<MultipathMEMEdge>& edges_from = nodes[i].edges_from;
                vector<MultipathMEMEdge>& edges_to = nodes[i].edges_to;
                
                // filter the edges out of this node
                size_t edges_removed_so_far = 0;
                for (size_t j = 0; j < edges_from.size(); j++) {
                    if (edge_idxs.count(make_pair(i, edges_from[j].to_idx))) {
                        swap(edges_from[j], edges_from[j - edges_removed_so_far]);
                    }
                    else {
                        edges_removed_so_far++;
                    }
                }
                edges_from.resize(edges_from.size() - edges_removed_so_far);
                
                // filter the edges into this node
                edges_removed_so_far = 0;
                for (size_t j = 0; j < edges_to.size(); j++) {
                    if (edge_idxs.count(make_pair(edges_to[j].to_idx, i))) {
                        swap(edges_from[j], edges_from[j - edges_removed_so_far]);
                    }
                    else {
                        edges_removed_so_far++;
                    }
                }
                edges_to.resize(edges_to.size() - edges_removed_so_far);
                
                // move the nodes we're going to keep into the front of the vector
                swap(nodes[i], nodes[i - nodes_removed_so_far]);
                
            }
            else {
                // skip this node and increase the count of the nodes we've removed
                nodes_removed_so_far++;
            }
        }
        
        // remove the nodes beyond those we kept
        nodes.resize(nodes.size() - nodes_removed_so_far);
        
        // renumber the indices in the edges
        for (MultipathMEMNode& node : nodes) {
            for (MultipathMEMEdge& edge : node.edges_from) {
                edge.to_idx -= removed_before[edge.to_idx];
            }
            for (MultipathMEMEdge& edge : node.edges_to) {
                edge.to_idx -= removed_before[edge.to_idx];
            }
        }
        
        // TODO: introduce logic for re-merging split overlap nodes?
    }
    
    void MultipathMEMAligner::query_node_matches(const Alignment& alignment, xg::XG& xgindex,
                                                 LRUCache<id_t, Node>& node_cache) {
        
        // group the nodes based on the MEM they originally came from
        // TODO: it would be better to just preserve this information somehow, rediscovering it is wasteful

        string::const_iterator read_begin = alignment.sequence().begin();
        
        unordered_map<pair<pos_t, int64_t>, vector<MultipathMEMNode*>> match_records;
        for (MultipathMEMNode& node : nodes) {
            /// make a record for this hit position, or get it if it already exists
            int64_t read_idx = (node.begin - read_begin) - node.offset;
            match_records[make_pair(node.start_pos, read_idx)].push_back(&node);
        }
        
        vector<pair<pair<pos_t, int64_t>, pair<string::const_iterator, string::const_iterator>>> matches;
        for (auto& match_record : match_records) {
            // put the nodes in order along the MEM
            vector<MultipathMEMNode*>& match_nodes = match_record.second;
            sort(match_nodes.begin(), match_nodes.end(),
                 [](const MultipathMEMNode* n1, const MultipathMEMNode* n2) {
                     return n1->offset < n2->offset;
                 });
            
            // combine the nodes into a single record of the match to query
            matches.emplace_back(match_record.first,
                                 make_pair(match_nodes.front()->begin - match_nodes.front()->offset,
                                           match_nodes.back()->end));
        }
        
        // sort the matches in descending order by length so that if there's a redundant
        // partial MEM we always find the containing MEM(s) first
        sort(matches.begin(), matches.end(),
             [](const pair<pos_t, pair<string::const_iterator, string::const_iterator>>& m1,
                const pair<pos_t, pair<string::const_iterator, string::const_iterator>>& m2) {
                 return m1.second.second - m1.second.first > m2.second.second - m2.second.first;
             });
        
        // map of node ids to the indices in the matches that contain them
        unordered_map<int64_t, vector<size_t>> node_matches;
        
        unordered_set<size_t> nonredundant_node_idxs;
        
        for (size_t i = 0; i < matches.size(); i++) {
            
            // the part of the read we're going to match
            string::const_iterator begin = matches[i].second.first;
            string::const_iterator end = matches[i].second.second;
            int64_t mem_length = end - begin;
            // the start of the hit
            pos_t hit_pos = matches[i].first.first;
            
            // check all MEMs that traversed this node to see if this is a sub-MEM
            bool is_partial_mem = false;
            for (size_t j : node_matches[id(hit_pos)]) {
                for (MultipathMEMNode* match_node : match_records[matches[j].first]) {
                    Path& path = match_node->path;
                    
                    int64_t relative_offset = begin - match_node->begin;
                    if (relative_offset < 0 || relative_offset >= match_node->end - match_node->begin) {
                        // the start of the hit does not fall on the same section of the read
                        // the match
                        continue;
                    }
                    
                    int64_t prefix_length = 0;
                    for (size_t k = 0; k < path.mapping_size(); k++) {
                        const Mapping& mapping = path.mapping(k);
                        // the length through this matpping
                        int64_t prefix_through_length = prefix_length + mapping.from_length();
                        if (relative_offset < prefix_length + mapping.from_length()) {
                            // the end of this mapping is past where this match would occur inside it
                            if (prefix_length <= relative_offset) {
                                // the possible start of this match is on this mapping
                                if (mapping.position().node_id() == id(hit_pos)
                                    && mapping.position().is_reverse() == is_rev(hit_pos)) {
                                    // we are on the node traversal we would expect if this is a false
                                    // partial match. technically this can happen when this isn't a false
                                    // partial match, but this will be extremely rare
                                    is_partial_mem = true;
                                    break;
                                }
                            }
                            else {
                                // we've gone past it
                                break;
                            }
                        }
                        prefix_length = prefix_through_length;
                    }
                    if (is_partial_mem) {
                        break;
                    }
                }
                if (is_partial_mem) {
                    break;
                }
            }
            
            // don't walk the match of false partial hits
            if (is_partial_mem) {
                continue;
            }
            
            // mark that this node passed the redundancy filter
            nonredundant_node_idxs.insert(i);
            
            // TODO: the following code is duplicative with Mapper::mem_positions_by_index
            
            // each record indicates the next edge index to traverse, the number of edges that
            // cannot reach a MEM end, and the positions along each edge out
            vector<pair<size_t, vector<pos_t>>> stack{make_pair(0, vector<pos_t>{hit_pos})};
            
            // it is possible for the MEM to have multiple exact hits from the same start position
            // but we will choose one greedily with DFS
            while (!stack.empty()) {
                
                size_t mem_idx = pos_stack.size() - 1;
                
                // which edge are we going to search out of this node next?
                size_t next_idx = pos_stack.back().first;
                
                // indicate which edge to check next
                pos_stack.back().first++;
                
                if (next_idx >= pos_stack.back().second.size()) {
                    // we have traversed all of the edges out of this position
                    
                    // backtrack to previous node
                    pos_stack.pop_back();
                    
                    // skip the forward search on this iteration
                    continue;
                }
                
                // increment to the next edge
                pos_stack.back().first.first++;
                
                pos_t graph_pos = pos_stack.back().second[next_idx];
                
                // does this graph position match the MEM?
                if (*(begin + mem_idx) != xg_cached_pos_char(graph_pos, xgindex, node_cache)) {
                    // increase the count of misses in this layer
                    pos_stack.back().first.second++;
                }
                else if (mem_idx < mem_length - 1) {
                    
                    // add a layer onto the stack for all of the edges out
                    pos_stack.push_back(make_pair(0, vector<pos_t>()));
                    
                    // fill the layer with the next positions
                    vector<pos_t>& nexts = pos_stack.back().second;
                    for (const pos_t& next_graph_pos : xg_cached_positions_bp_from(graph_pos, 1, false,
                                                                                   xgindex, node_cache)) {
                        nexts.push_back(next_graph_pos);
                    }
                }
                else {
                    break;
                }
            }
            
            assert(stack.size() == mem_length);
            
            // add pointers to this path to the index for node IDs
            int64_t node_id = 0;
            for (size_t j = 0; j < mem_length; j++) {
                int64_t id_here = id(stack[j].second[stack[j].first - 1]);
                if (id_here != node_id) {
                    node_matches[id_here].push_back(i);
                    node_id = id_here;
                }
            }
            
            vector<MultipathMEMNode*>& match_nodes = match_records[matches[i].first];
            for (MultipathMEMNode* match_node : match_nodes) {
                // get the path from the node
                Path& path = match_nodes->path;
                
                // initialize a the first mapping
                int64_t rank = 1;
                Mapping* mapping = path.add_mapping();
                mapping->set_rank(rank);
                Position* position = mapping->add_position();
                
                // get the first position of this match
                pos_t first_pos = stack[match_node->offset].second[stack[match_node->offset].first - 1];
                
                // keep track of how long the match is to one node
                int64_t node_match_len = 1;
                
                // transfer information to the Position
                position->set_offset(offset(first_pos));
                position->set_node_id(id(first_pos));
                position->set_is_reverse(is_rev(first_pos));
                
                // check the rest of the positions
                int64_t match_node_end = match_node->offset + (match_node->end - match_node->begin);
                for (int64_t j = match_node->offset + 1; j < match_node_end; i++) {
                    // get the next position
                    pos_t next_pos = stack[j].second[stack[j].first - 1];
                    
                    if (is_rev(next_pos) != position->is_reverse() || id(next_pos) != position->node_id()) {
                        // we've crossed an edge, finish the mapping
                        mapping->set_from_length(node_match_len);
                        mapping->set_to_length(node_match_len);
                        
                        // and start a new one
                        rank++;
                        mapping = path.add_mapping();
                        mapping->set_rank(rank);
                        position = mapping->add_position();
                        
                        // transfer the information to it
                        position->set_offset(offset(next_pos));
                        position->set_node_id(id(next_pos));
                        position->set_is_reverse(is_rev(next_pos));
                        
                        node_match_len = 1;
                    }
                    else {
                        // this is still matching the same node
                        node_match_len++;
                    }
                }
                // finish the last mapping
                mapping->set_from_length(node_match_len);
                mapping->set_to_length(node_match_len);
            }
        }
        
        // if walking out the matches uncovered any redundant sub-hits, prune
        // them from the graph
        if (nonredundant_node_idxs.size() < nodes.size()) {
            
            // gather all of the edges between
            unordered_set<pair<size_t, size_t>> nonredundant_edge_idxs;
            for (size_t i : nonredundant_node_idxs) {
                for (MultipathMEMEdge& edge : nodes[i].edges_from) {
                    if (nonredundant_node_idxs.count(edge.to_idx)) {
                        nonredundant_edge_idxs.emplace(i, edge.to_idx);
                    }
                }
            }
            
            prune_graph(nonredundant_node_idxs, nonredundant_edge_idxs);
        }
    }
    
    void MultipathMEMAligner::remove_snarls(const Alignment& alignment,
                                            const SnarlManager& snarl_manager,
                                            const QualAdjAligner& aligner) {
        // TODO: removing all Snarl interiors is going to be too aggressive, but I'm not sure what
        // a better principled alternative would be
        // for now I'm going to remove only the bottom level Snarl traversed by the path because that's
        // easy
        
        unordered_map<pair<int64_t, bool>, const Snarl*> snarl_start_index = snarl_manager.snarl_start_index();
        unordered_map<pair<int64_t, bool>, const Snarl*> snarl_end_index = snarl_manager.snarl_end_index();
        
        // we are going to add nodes, so we record the number of nodes at the beginning of the loop
        // and only iterate through those
        for (int64_t i = 0, starting_num_nodes = nodes.size(); i < starting_num_nodes; i++) {
            MultipathMEMNode& node = nodes[i];
            Path& path = node.path;
            
            // the start and end indices of the sections to cut out
            vector<pair<int64_t, int64_t>> cut_intervals;
            
            // the length of the read up to this point
            vector<int64_t> prefix_length;
            int64_t cumul_len = 0;
            
            // indicate that we can cut from the beginning
            int64_t cut begin = 0;
            
            for (size_t j = 0; j < path.mapping_size(); j++) {
                
                // update the prefix length index
                prefix_length.push_back(cumul_len);
                cumul_len += mapping_to_length(path.mapping(j));
                
                const Position& position = path.mapping(j).position();
                
                if (snarl_end_index.count(make_pair(position.node_id(), !position.is_reverse()))
                    snarl_start_index.count(make_pair(position.node_id(), !position.is_reverse()))
                    && cut_begin >= 0) {
                    // a cut ends here and we haven't made a cut yet in this branch of the tree
                    
                    // record the interval covered by this site
                    cut_intervals.emplace_back(cut_begin, j);
                    // mark that we made a cut in this branch
                    cut_begin = -1;
                }
                
                // we allow a nested site to reset the beginning of the cut
                if (snarl_start_index.count(make_pair(position.node_id(), position.is_reverse()))
                    snarl_end_index.count(make_pair(position.node_id(), position.is_reverse()))) {
                    // a cut starts here
                    cut_begin = j + 1;
                }
            }
            
            // add final entry to prefix length index
            prefix_length.push_back(cumul_len);
            
            // add final cut if it exists
            if (cut_begin > 0) {
                cut_intervals.emplace_back(cut_begin, path.mapping_size());
            }
            
            // create a node for each new cut section
            
            // the index of the last new node we created
            int64_t last_node_idx = i;
            // the beginning of the last cut we made
            int64_t last_cut_begin = path.mapping_size();
            for (auto iter = cut_intervals.rbegin(); iter != cut_intervals.rend(); iter++) {
                pair<int64_t, int64_t>& cut_interval = *iter;
                
                // make a new node with the section between the last cut, but don't make a new node
                // if we would be cutting all the way to the end
                if (cut_interval.second != path.mapping_size() && cut_interval.first != 0) {
                    string::const_iterator new_begin = node.begin + prefix_length[cut_interval.second];
                    string::const_iterator new_end = node.begin + prefix_length[last_cut_begin];
                    int32_t new_score = aligner.score_exact_match(new_begin, new_end,
                                                                  alignment.quality().begin()
                                                                  + (new_begin - alignment.sequence().begin()));
                    
                    nodes.emplace_back(new_begin,
                                       new_end,
                                       node.start_pos,
                                       node.offset + prefix_length[cut_interval.second],
                                       new_score);
                    
                    MultipathMEMNode& new_node = nodes.back();
                    
                    // transfer over the part of the path that corresponds to this inter-cut segment
                    for (int64_t j = cut_interval.second; j < last_cut_begin; j++) {
                        *(new_node.path.add_mapping()) = path.mapping(j);
                    }
                    
                    int64_t new_idx = nodes.size() - 1;
                    
                    // transfer over edges
                    new_node.edges_from = std::move(node.edges_from);
                    // update backwards references
                    for (MultipathMEMEdge& edge_from : new_node.edges_from) {
                        for (MultipathMEMEdge& edge_to : nodes[edge_from.to_idx].edges_to) {
                            if (edge_to.to_idx == i) {
                                edge_to.to_idx = new_idx;
                                break;
                            }
                        }
                    }
                    
                    // add an edge spanning the cut
                    new_node.edges_to.emplace_back(i, 0);
                    node.edges_from.emplace_back(new_idx, 0);
                    
                }
                
                // mark the end of the next node
                last_cut_begin = cut_interval.first;
            }
            
            if (!cut_intervals.empty()) {
                // make a new path to replace the one on the original node and update its information
                
                // the interval of mappings that correspond to what will be left on this node
                int64_t start_idx = 0;
                int64_t end_idx = cut_intervals.front().first;
                // if the cut goes all the way to the end of the node, move the interval down
                if (end_idx == 0) {
                    start_idx = cut_intervals.front().second;
                    end_idx = cut_intervals.size() > 1 ? cut_intervals[1].first : path.mapping_size();
                }
                
                // copy the segment of the path and replace the original
                Path temp_path;
                for (int64_t j = start_idx; j < end_idx; j++) {
                    *(temp_path.add_mapping()) = path.mapping(j);
                }
                path = std::move(temp_path);
                
                // update the other information on the node
                node.end = node.begin + prefix_length[end_idx];
                node.begin = node.begin + prefix_length[start_idx];
                node.score = aligner.score_exact_match(node.begin, node.end,
                                                       alignment.quality().begin()
                                                       + (node.end - alignment.sequence().begin()));
            }
        }
    }
    
    void MultipathMEMAligner::align_internal_edges(const Alignment& alignment,
                                                   const QualAdjAligner& aligner,
                                                   const int8_t full_length_bonus,
                                                   size_t max_alns_per_edge,
                                                   xg::XG& xgindex,
                                                   LRUCache<id_t, Node>& node_cache) {
        
        int8_t match_score = aligner.adjusted_score_matrix[25 * aligner.max_qual_score];
        
        unordered_set<size_t> prune_nodes;
        unordered_set<pair<size_t, size_t>> prune_edges;
        
        unordered_map<pair<size_t, size_t>, vector<Alignment>> edge_alns;
        
#pragma omp parallel for
        for (size_t i = 0; i < nodes.size(); i++) {
            
            MultipathMEMNode& node = nodes[i];
            
            pos_t left_cut_pos = final_position(node.path);
            
            size_t max_gap_len_left = longest_detectable_gap(alignment,
                                                             node.end,
                                                             match_score,
                                                             aligner,
                                                             full_length_bonus);
            
            for (MultipathMEMEdge& edge : node.edges_from) {
                
                size_t next_idx = edge.to_idx;
                MultipathMEMNode& next_node = nodes[next_idx];
                
                // don't extract subgraphs for split edges
                if (next.start_pos == next_node.start_pos &&
                    node.offset + (node.end - node.begin) == next_node.offset) {
                    continue;
                }
                
                pos_t right_cut_pos = initial_position(next_node.path);
                
                size_t max_gap_len_right = longest_detectable_gap(alignment,
                                                                  next_node.begin,
                                                                  match_score,
                                                                  aligner,
                                                                  full_length_bonus);
                
                int64_t read_segment_len = next_node.begin - node.end;
                size_t max_dist = min(max_gap_len_left, max_gap_len_right) + read_segment_len;
                
                
                // a heuristic to avoid querying large subgraphs looking for terminal cycles when the
                // two positions are on the same node and can reach each other (this tends to be a costly
                // case when the cycle doesn't actually exist). if there is a terminal cycle, it will
                // have to traverse the entire node, so if the read segment is not comparably sized to
                // the node sequence we don't look for cycles
                bool detect_terminal_cycles = true;
                if (id(left_cut_pos) == id(right_cut_pos)
                    && is_rev(left_cut_pos) == is_rev(right_cut_pos)
                    && offset(left_cut_pos) < offset(right_cut_pos)) {
                    
                    size_t node_len = node_length(id(left_cut_pos));
                    
                    detect_terminal_cycles = read_segment_len >= (3 * node_length) / 4
                }
                
                // TODO: unless i'm going to hang onto all of these ID translators, i'll
                // need to do the multi-alignment in this loop
                
                // get the graph between the two nodes
                Graph connecting_graph;
                auto id_trans = algorithms::extract_connecting_graph(xgindex, connecting_graph, max_dist,
                                                                     left_cut_pos, right_cut_pos, false,
                                                                     detect_terminal_cycles,
                                                                     true, true, true, &node_cache);
                
                // we could not find a path between the positions, but this might be because they
                // overlap each other unncessarily in the graph, so we look for overlap and requery
                if (!connecting_graph.node_size()) {
                    // TODO: another place that might make sense to do the longest overlap is
                    // after pruning with traceback (easier, but doesn't take into account the graph connectivity)
                    
                    // backtrack to reconstruct the MEM if it has been split
                    MultipathMEMNode& mem_start_node = node;
                    
                    while (mem_start_node.offset > 0) {
                        
                        MultipathMEMNode* prev = nullptr;
                        for (MultipathMEMEdge& mem_edge : mem_start_node.edges_to) {
                            MultipathMEMNode& from = nodes[edge.to_idx];
                            // is this part of the same MEM?
                            if (from.start_pos == mem_start_node.start_pos
                                && from.end + (mem_start_node.offset - from.offset) == mem_start_node.begin) {
                                prev = &from;
                                break;
                            }
                        }
                        
                        if (prev) {
                            mem_start_node = *prev;
                        }
                        else {
                            // the real start node may have been pruned off
                            break;
                        }
                    }
                    
                    // forward search to reconstruct next node's MEM if it has been split
                    MultipathMEMNode& mem_end_node = next_node;
                    while (!mem_end_node.edges_from.empty()) {
                        
                        MultipathMEMNode* next = nullptr;
                        for (MultipathMEMEdge& mem_edge : mem_end_node.edge_from) {
                            MultipathMEMNode& to = nodes[mem_edge.to_idx];
                            // is this part of the same MEM?
                            if (to.start_pos == mem_end_node.start_pos
                                && mem_end_node.begin + (to.offset - mem_end_node.offset) == to.begin) {
                                next = &from;
                                break;
                            }
                        }
                        
                        if (next) {
                            mem_end_node = *next;
                        }
                        else {
                            break;
                        }
                    }
                    
                    // find the longest overlap between the two MEMs
                    SuffixTree suffix_tree(mem_start_node.begin, node.end);
                    size_t overlap = suffix_tree.longest_overlap(next_node.begin, mem_end_node.end);
                    
                    // did we identify a non-trivial overlap that might rescue this unreachable edge?
                    if (overlap && overlap != mem_end_node.end - next_node.begin) {
                        // find the node that corresponds to this part of the split
                        size_t split_idx = edge.to_idx;
                        MultipathMEMNode& split_node = next_node
                        size_t beginning_offset = next_node.offset;
                        while (split_node.offset - beginning_offset + (split_node.end - split_node.begin) <= overlap) {
                            for (MultipathMEMEdge& split_edge : split_node.edges_from) {
                                split_idx = split_edge.to_idx;
                                MultipathMEMNode& to = nodes[split_idx];
                                // is this part of the same MEM?
                                if (to.start_pos == split_node.start_pos
                                    && split_node.begin + (to.offset - split_node.offset) == to.begin) {
                                    split_node = to;
                                    break;
                                }

                            }
                        }
                        
                        int64_t dist_on_node = overlap - (split_node.offset - beginning_offset);
                        string::const_iterator split_point = split_node.begin + dist_on_node;
                        
                        
                        int32_t section_score = aligner.score_exact_match(split_node.begin, split_point,
                                                                          alignment.quality().begin() + (split_node.begin -
                                                                                                         alignment.sequence().begin()));
                        
                        // make a new node
                        nodes.emplace_back(split_point, split_node.end, split_node.start_pos,
                                           split_node.offset + dist_on_node, split_node.score - section_score);
                        
                        // update old node information
                        split_node.end = split_point;
                        split_node.score = section_score;
                        
                        // split the path
                        
                        // find the split point
                        size_t split_mapping_idx = 0;
                        size_t split_mapping_offset = 0;
                        size_t length_so_far = 0;
                        for (size_t j = 0; j < split_node.path.mapping_size(); j++) {
                            size_t next_length_so_far = length_so_far + mapping_to_length(split_node.path.mapping(j));
                            if (next_length_so_far > dist_on_node) {
                                split_mapping_idx = j;
                                split_mapping_offset = dist_on_node - length_so_far;
                            }
                            length_so_far = next_length_so_far;
                        }
                        
                        // copy the end of the path past the split point
                        Mapping* mapping = new_node.path.add_mapping();
                        mapping->mutable_position()->set_offset(split_node.path.mapping(split_mapping_idx).position().offset()
                                                                + split_mapping_offset);
                        mapping->mutable_position()->set_node_id(split_node.path.mapping(split_mapping_idx).position().node_id());
                        mapping->set_from_length(split_node.path.mapping(split_mapping_idx).from_length() - split_mapping_offset);
                        mapping->set_to_length(mapping->from_length());
                        for (size_t j = split_mapping_idx + 1; j < split_node.path.mapping_size(); j++) {
                            mapping = new_node.path.add_mapping();
                            mapping->mutable_position()->set_node_id(split_node.path.mapping(j).position().node_id());
                            mapping->set_from_length(split_node.path.mapping(j).from_length());
                            mapping->set_to_length(mapping->from_length());
                        }
                        
                        // duplicate the path before the split point and replace the original
                        Path temp_path;
                        for (size_t j = 0; j < split_mapping_idx; j++) {
                            mapping = temp_path.add_mapping();
                            mapping->mutable_position()->set_node_id(split_node.path.mapping(j).position().node_id());
                            mapping->mutable_position()->set_offset(split_node.path.mapping(j).position().offset());
                            mapping->set_from_length(split_node.path.mapping(j).from_length());
                            mapping->set_to_length(mapping->from_length());
                        }
                        // don't copy the first last mapping if the entire thing was moved to the next path
                        if (split_mapping_offset > 0) {
                            mapping = temp_path.add_mapping();
                            mapping->mutable_position()->set_offset(split_node.path.mapping(split_mapping_idx).position().offset());
                            mapping->mutable_position()->set_node_id(split_node.path.mapping(split_mapping_idx).position().node_id());
                            mapping->set_from_length(split_mapping_offset);
                            mapping->set_to_length(mapping->from_length());
                        }
                        split_node.path = std::move(temp_path);
                        
                        // move forward edges onto the new node
                        nodes.back().edges_from = std::move(split_node.edges_from);
                        split_node.edges_from.clear();
                        
                        // create an edges between the split sides
                        size_t new_idx = nodes.size() - 1;
                        split_node.edges_from.emplace_back(new_idx, 0);
                        split_node.edges_to.emplace_back(split_idx, 0);
                        
                        // move the old edge that turned out to be unreachable
                        nodes.back().edges_to.emplace_back(i, edge.weight);
                        edge.to_idx = new_idx;
                        auto new_end = std::remove_if(next_node.edges_to.begin(), next_node.edges_to.end(),
                                                      auto [&](const MultipathMEMEdge& e) {return e.to_idx == i});
                        next_node.edges_to.resize(new_end - next_node.edges_to.begin());
                        
                        // flag any earlier part of the MEM that was only connected by this edge
                        if (next_node.edges_to.size() == 0) {
                            prune_nodes.insert(next_idx);
                            while (next_node.edges_from.size() == 1) {
                                next_idx = next_node.edges_from[0].to_idx;
                                next_node = nodes[next_idx];
                                if (next_node.edges_to.size() == 1) {
                                    prune_nodes.insert(next_idx);
                                }
                                else {
                                    break;
                                }
                            }
                        }
                        
                        // update the reference to the new node
                        next_node = new_node;
                        next_idx = new_idx;
                        
                        // requery the graph
                        right_cut_pos = initial_position(split_node.path);
                        id_trans = algorithms::extract_connecting_graph(xgindex, connecting_graph, max_dist,
                                                                        left_cut_pos, right_cut_pos, false,
                                                                        detect_terminal_cycles,
                                                                        true, true, true, &node_cache);
                        
                    }
                }
                
                // did the positions end up being connected under the max length?
                if (connecting_graph.node_size()) {
                    
                    // wrap the graph with a function to make a VG object
                    bool passed_graph = false;
                    auto pass_graph = [&](Graph& graph) {
                        if (passed_graph) {
                            return false;
                        }
                        else {
                            graph = std::move(connecting_graph);
                            passed_graph = true;
                            return true;
                        }
                    };
                    
                    VG vg(pass_graph);
                    
                    string::const_iterator base_qual_begin = alignment.quality().begin()
                                                             + (node.end - alignment.sequence().begin());
                    
                    // transfer in the
                    Alignment edge_aln;
                    edge_aln.set_sequence(string(node.end, next_node.begin));
                    edge_aln.set_quality(string(base_qual_begin, base_qual_begin + (next_node.begin - node.end)));
                    
                    vector<Alignment>& edge_alt_alns = edge_alns[make_pair(i, next_idx)];
                    
                    // TODO !!!!!!!!!! Dagify
                    
                    
                    aligner.align_global_banded_multi(edge_aln, edge_alt_alns,
                                                      connecting_graph, max_alns_per_edge,
                                                      1, true);
                    
                    // translate the node IDs in the alignments
                    for (Alignment& aln : edge_alt_alns) {
                        Path* path = aln.mutable_path();
                        for (size_t j = 0; j < path->mapping_size(); j++) {
                            Position* pos = path->mutable_mapping(j)->mutable_position();
                            pos->set_node_id(id_trans[pos->node_id()]);
                        }
                    }
                    
                }
                else {
                    // the edge does not represent a true connection
                    prune_edges.insert(make_pair(i, next_idx));
                }
                
            }
        }
        
    }
    
    MultipathMEMAligner::TracebackManager::TracebackManager(const vector<MultipathMEMNode>& nodes,
                                                            size_t max_num_tracebacks) :
        max_num_tracebacks(max_num_tracebacks), nodes(nodes),
        traceback_end_candidates(priority_queue<size_t,
                                                vector<size_t>,
                                                DPScoreComparator>(DPScoreComparator(nodes),
                                                                   range_vector(0, nodes.size())))
    {
        init_alt_traceback_stack();
    }
    
    inline void MultipathMEMAligner::TracebackManager::init_alt_traceback_stack() {
        size_t node_idx = traceback_end_candidates.top();
        traceback_end_candidates.pop();
        alt_tracebacks.push_back(make_pair(vector<Deflection>(), nodes[node_idx].dp_score));
        alt_tracebacks.back().first.push_back(Deflection(node_idx, node_idx));
    }
    
    inline void MultipathMEMAligner::TracebackManager::next() {
        if (curr_deflxn != (*curr_traceback).first.end()) {
            cerr << "warning:[MultipathMEMAligner] moving on to next alternate traceback without taking all deflections" << endl;
        }
        
        // advance to the next
        curr_traceback++;
        
        // find any non-redundant, non-branching alternate tracebacks
        find_local_traceback_end();
        
    }
    
    inline bool MultipathMEMAligner::TracebackManager::has_next() {
        return (curr_traceback != alt_tracebacks.end());
    }
    
    inline size_t MultipathMEMAligner::TracebackManager::get_traceback_start() {
        curr_deflxn = (*curr_traceback).first.begin();
        size_t node_idx = (*curr_deflxn).to_idx;
        curr_deflxn++;
        return node_idx;
    }
    
    inline void MultipathMEMAligner::TracebackManager::mark_traced(const size_t node_idx) {
        traced_node_idxs.insert(node_idx);
    }
    
    inline bool MultipathMEMAligner::TracebackManager::at_next_deflection(const size_t node_idx) {
        return curr_deflxn == (*curr_traceback).first.end() ? false : (*curr_deflxn).from_idx == node_idx;
    }
    
    inline int32_t MultipathMEMAligner::TracebackManager::curr_traceback_score() {
        return (*curr_traceback).second;
    }
    
    inline void MultipathMEMAligner::TracebackManager::find_local_traceback_end() {
        
        if (traceback_end_candidates.empty() ||
            (curr_traceback == alt_tracebacks.end() && alt_tracebacks.size() == max_num_tracebacks)) {
            // there are no more potential candidate local traceback starts or we are at the end
            // of a finished stack
            return;
        }
        else if (curr_traceback == alt_tracebacks.end() ? false :
                 nodes[traceback_end_candidates.top()].dp_score <= (*curr_traceback).second) {
            // our current traceback is higher scoring than any of the potential candidates
            return;
        }
        
        // search the nodes in reverse score order to find traceback starts that are not
        // just shorter versions of previously traversed tracebacks
        while (!traceback_end_candidates.empty()) {
            // get the next highest scoring unchecked node
            size_t candidate = traceback_end_candidates.top();
            int32_t candidate_score = nodes[candidate].dp_score;
            
            // we can stop looking if our current tracebook is already higher scoring than this one
            if (curr_traceback == alt_tracebacks.end() ? false : candidate_score <= (*curr_traceback).second) {
                break;
            }
            
            // this node is either going to be added to the stack or skipped because we've already
            // traced on it, so we can remove it from the list
            traceback_end_candidates.pop();
            
            // we've already seen the node inside another traceback, so it can only contribute a
            // shorter version of an already observed traceback
            if (!traced_node_idxs.count(candidate)) {
                continue;
            }
            
            // add a new traceback beginning here into the current position in the stack
            curr_traceback = alt_tracebacks.insert(curr_traceback,
                                                   make_pair(vector<Deflection>(), candidate_score));
            // set its starting point
            (*curr_traceback).first.push_back(Deflection(candidate, candidate));
            // remove any extra tracebacks if necessary
            trim_traceback_stack();
            // hold off searching the rest of the candidates until later so that all nodes from
            // this traceback will be marked as traced
            break;
        }
    }
    
    inline void MultipathMEMAligner::TracebackManager::propose_deflection(const size_t from,
                                                                          const size_t to,
                                                                          const int32_t score) {
        // only propose deflections if we're going through a new untraversed section of the traceback
        if (curr_deflxn != (*curr_traceback).first.end()) {
            return;
        }
        
        // is the score good enough to be put on the stack?
        if (score <= alt_tracebacks.back().second
            && alt_tracebacks.size() >= max_num_tracebacks) {
            return;
        }
        
        insert_traceback((*curr_traceback).first, from, to, score);
    }
    
    inline size_t MultipathMEMAligner::TracebackManager::deflect() {
        size_t to_idx = (*curr_deflxn).to_idx;
        curr_deflxn++;
        return to_idx;
    }
    
    inline void MultipathMEMAligner::TracebackManager::trim_traceback_stack() {
        if (alt_tracebacks.size() > max_num_tracebacks) {
            alt_tracebacks.pop_back();
        }
    }
    
    inline void MultipathMEMAligner::TracebackManager::insert_traceback(const vector<Deflection>& traceback_prefix,
                                                                        const size_t from, const size_t to,
                                                                        const int32_t score) {
        // find position in stack where this should go
        auto insert_after = alt_tracebacks.rbegin();
        while (score > (*insert_after).second) {
            insert_after++;
            if (insert_after == alt_tracebacks.rend()) {
                break;
            }
        }
        
        // insert if score is high enough or stack is not full yet
        if (insert_after != alt_tracebacks.rbegin()
            || alt_tracebacks.size() < max_num_tracebacks) {
            
            // create a new traceback here
            auto new_traceback = alt_tracebacks.emplace(insert_after.base(),
                                                        vector<Deflection>(traceback_prefix),
                                                        score);
            
            // add the final deflection
            (*new_traceback).first.emplace_back(from, to);
        }
        
        // remove lowest scoring tracebacks if stack is over capacity
        trim_traceback_stack();
    }
}
