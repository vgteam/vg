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
                                             int8_t full_length_bonus,
                                             size_t num_pruning_tracebacks) {
        
        init_first_pass_graph(alignment, mems, aligner, xgindex, full_length_bonus);
        perform_dp();
        prune_to_nodes_on_tracebacks(num_pruning_tracebacks);
        query_node_matches(alignment, xgindex, node_cache);
        // TODO: add logic to prune out Snarls
    }
    
    void MultipathMEMAligner::init_first_pass_graph(const Alignment& alignment,
                                                    const vector<MaximalExactMatch>& mems,
                                                    const QualAdjAligner& aligner,
                                                    const xg::XG& xgindex,
                                                    int8_t full_length_bonus){
    
        // TODO: handle inversions -- do I just need extra paths?
        // for now just not checking for orientation consistency and hoping it works out (but this
        // is only likely when inversion is small so that distance is not overestimated too much)
        auto distance = [&](pos_t& pos_1, pos_t& pos_2) {
            return xgindex.min_approx_path_distance(vector<string>(), id(pos_1), id(pos_2))
                   + (int64_t) offset(pos_1) - (int64_t) offset(pos_2);
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
            
            // TODO: should I include the full length alignment bonus in here? easy to include it
            // on the MEM itself, but hard to factor in the dangling end of the read near a mismatch
            
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
        
        // first pass: identify all pairs that could conceivably have edges
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
        
        vector<pair<pair<int64_t, int64_t>, int64_t>> overlaps;
        
        for (const pair<pair<int64_t, int64_t>, size_t>& pair_distance : pair_distances) {
            
            MultipathMEMNode& node_1 = nodes[pair_distance.first.first];
            MultipathMEMNode& node_2 = nodes[pair_distance.first.second];
            
            // TODO: need to look up how approximate distance function handles strand, sign, etc.
            int64_t absolute_dist = abs(pair_distance.second);
            
            // can we make an edge between these without trimming an overlap?
            if (node_1.end < node_2.begin) {
                // the amount of read between these MEMs
                int64_t between_length = node_2.begin - node_1.end;
                
                // TODO: assuming matches for intervening sequence will make it difficult to prune
                // off erroneous sub-MEMs that we aren't going to filter out at this stage
                
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
            
            // figure out how much of the first MEM could potentially overlap with the second
            // TODO: some of these could represent the same underlying MEM, might be worth finding a
            // way to only build suffix tree once per MEM
            SuffixTree suffix_tree(node_1.begin, node_1.end);
            size_t overlap = suffix_tree.longest_overlap(node_2.begin, node_2.end);
            
            // record any non-trivial overlaps
            if (overlap) {// && overlap != node_2.end - node_2.begin) {
                overlaps.push_back(make_pair(pair_distance.first, overlap));
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
                            // TODO: detect whether this is an overlap edge and add it to the edge list
                            // without actually included a basically identical traceback?
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
                edge_idxs.insert(make_pair(prev_idx, traceback_idx));
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
        
//        // sort the matches in descending order by length so that if there's a redundant
//        // sub-MEM we always find the containing MEM(s) first
//        sort(matches.begin(), matches.end(),
//             [](const pair<pos_t, pair<string::const_iterator, string::const_iterator>>& m1,
//                const pair<pos_t, pair<string::const_iterator, string::const_iterator>>& m2) {
//                 return m1.second.second - m1.second.first > m2.second.second - m2.second.first;
//             });
//        
//        // map of node ids to the indices in the matches that contain them
//        unordered_map<int64_t, vector<size_t>> node_matches;
//        
//        for (size_t i = 0; i < matches.size(); i++) {
//            
//        }
        
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
