//
//  multipath_mem_aligner.cpp
//

#include "multipath_mem_aligner.hpp"

namespace vg {
    
    MultipathMEMAligner::MultipathMEMAligner(const Alignment& alignment,
                                             const vector<MaximalExactMatch>& mems,
                                             const QualAdjAligner& aligner,
                                             const function<int64_t(pos_t&,pos_t&)>& distance,
                                             int8_t full_length_bonus) {
        
        
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
            
            for (gcsa::node_type mem_hit : mems.nodes) {
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
                    int64_t dist = distance(node_1.start_pos, node_2.start_pos);
                    
                    
                    // the minimum of the max detectable gap is the max gap here
                    size_t max_gap_length = min(maximum_detectable_gaps[i].first,
                                                maximum_detectable_gaps[j].second);
                    
                    // record all pairs that are under the minimum
                    if (node_2.begin - node_1.end < max_gap_length
                        && dist - (node_1.end - node_1.begin) < max_gap_length) {
                        pair_distances[make_pair(i, j)] = dist;
                    }
                }
                else {
                    int64_t dist = distance(node_2.start_pos, node_1.start_pos);
                    
                    
                    // the minimum of the max detectable gap is the max gap here
                    size_t max_gap_length = min(maximum_detectable_gaps[j].first,
                                                maximum_detectable_gaps[i].second);
                    
                    // record all pairs that are under the minimum
                    if (node_1.begin - node_2.end < max_gap_length
                        && dist - (node_2.end - node_2.begin) < max_gap_length) {
                        pair_distances[make_pair(j, i)] = dist;
                    }
                }
            }
        }
        
        vector<pair<pair<int64_t, int64_t>, int64_t>> overlaps;
        
        for (const <pair<int64_t, int64_t>, size_t>* pair_distance : pair_distances) {
            
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
                    weight = between_length * match_score - aligner.scaled_gap_open - (remaining - 1) * alginer.scaled_gap_extension;
                }
                else {
                    // need an insertion to account for gap
                    int64_t extra = between_length - absolute_dist;
                    weight = absolute_dist * match_score - aligner.scaled_gap_open - (extra - 1) * alginer.scaled_gap_extension;
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
        
        size_t num_original_nodes = nodes.size();
        
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
            
            int64_t absolute_dist = abs(pair_distance[pair_overlap.first] - pair_overlap.second);
            
            // find the maximum possible score (goal is sensitivity at this stage)
            int32_t weight;
            if (between_length == absolute_dist) {
                // can account for gap with only mismatches
                weight = between_length * match_score;
            }
            else if (between_length < absolute_dist) {
                // need a deletion to account for gap
                int64_t remaining = absolute_dist - between_length;
                weight = between_length * match_score - aligner.scaled_gap_open - (remaining - 1) * alginer.scaled_gap_extension;
            }
            else {
                // need an insertion to account for gap
                int64_t extra = between_length - absolute_dist;
                weight = absolute_dist * match_score - aligner.scaled_gap_open - (extra - 1) * alginer.scaled_gap_extension;
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
        
        int8_t max_match_score = aligner.adjusted_score_matrix[25 * aligner.max_qual_score];
        
        // algebraic solution for when score is > 0 assuming perfect match other than gap
        return min((same_side_remaining * max_match_score + mem_score + full_length_bonus
                    - aligner.scaled_gap_open) / aligner.scaled_gap_extension,
                   (other_side_remaining * max_match_score + full_length_bonus
                    - aligner.scaled_gap_open) / aligner.scaled_gap_extension) + 1;
        
    }
    
    void MultipathMEMAligner::topological_order(vector<size_t>& order_out) {
        
        // initialize return value
        order_out.clear();
        order_out.resize(nodes.size());
        size_t order_idx = nodes.size() - 1;
        
        // initialize iteration structures
        vector<bool> enqueued = vector<bool>(nodes.size());
        vector<int> edge_index = vector<int>(nodes.size());
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
    
    void MultipathMEMAligner::fill_dp() {
        
        // as in local alignment, minimum score is the score of node itself
        for (size_t i = 0; i < nodes.size(); i++) {
            nodes[i].dp_score = nodes[i].score;
        }
        
        vector<size_t> order;
        topological_order(order);
        
        for (size_t i : topological_order) {
            
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
    
    unordered_set<size_t> MultipathMEMAligner::find_nodes_on_tracebacks(size_t num_tracebacks) {
        
        TracebackManager traceback_manager(nodes, num_tracebacks);
        
        // locate the optimum traceback
        
        traceback_manager.find_local_traceback_end();
        
        for (; traceback_manager.has_next(); traceback_manager.next()) {
            size_t traceback_idx = traceback_manager.get_traceback_start();
            traceback_manager.mark_traced(traceback_idx);
            MultipathMEMNode& node = nodes[traceback_idx];
            int32_t traceback_score = traceback_manager.curr_traceback_score();
            
            while (node.dp_score != node.score) {
                
                if (traceback_manager.at_next_deflection(traceback_idx)) {
                    traceback_idx = traceback_manager.deflect();
                }
                else {
                    int32_t backward_score = node.dp_score - node.score;
                    bool found_trace = false;
                    MultipathMEMEdge& trace_edge;
                    for (MultipathMEMEdge& edge : nodes[traceback_idx].edges_to) {
                        int32_t score_diff = backward_score - (nodes[edge.to_idx].dp_score + edge.weight);
                        if (score_diff == 0 && !found_trace) {
                            found_trace = true;
                            trace_edge = edge;
                        }
                        else {
                            traceback_manager.propose_deflection(traceback_idx, edge.to_idx,
                                                                 traceback_score - score_diff);
                        }
                    }
                    if (!found_trace) {
                        cerr << "error:[MultipathMEMAligner] traceback error in MEM aligner" << endl;
                        exit(1);
                    }
                    
                    traceback_idx = trace_edge.to_idx;
                }
                
                traceback_manager.mark_traced(traceback_idx);
                node = nodes[traceback_idx];
            }
        }
        
        return std::move(traceback_manager.traced_node_idxs);
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
    
    MultipathMEMAligner::TracebackManager::TracebackManager(const vector<MultipathMEMNode>& nodes,
                                                              size_t max_num_tracebacks) :
        max_num_tracebacks(max_num_tracebacks), nodes(nodes)
    {
        DPScoreComparator comp = DPScoreComparator(nodes);
        pair<IncrementIter, IncrementIter> node_idxs = range(0, nodes.size());
        
        traceback_end_candidates = priority_queue<size_t,
                                                  vector<size_t>,
                                                  DPScoreComparator>(node_idxs.first
                                                                     node_idxs.second,
                                                                     comp);
        
        
    }
    
    inline void MultipathMEMAligner::TracebackManager::next() {
        if (curr_deflxn != (*curr_traceback).first.end()) {
            cerr << "warning:[MultipathMEMAligner] moving on to next alternate alignment without taking all deflections" << endl;
        }
        
        // look for tracebacks that don't branch out of previous tracebacks
        find_local_traceback_end();
        
        curr_traceback++;
    }
    
    inline void MultipathMEMAligner::TracebackManager::find_local_traceback_end() {
        
        // try to find a new end of a traceback if we haven't found any tracebacks yet or there
        // are potential traceback ends that are better the alternates we've already found
        bool attempt_to_find = true;
        if (!alt_tracebacks.empty()) {
            attempt_to_find = (nodes[traceback_end_candidates.top()].dp_score > alt_tracebacks.back().second);
        }
        
        if (attempt_to_find) {
            while (!traceback_end_candidates.empty()) {
                // get the next highest scoring unchecked node
                size_t candidate = traceback_end_candidates.top();
                traceback_end_candidates.pop();
                
                // if we've already found better alternates we can stop
                int32_t candidate_score = nodes[candidate].dp_score;
                if (candidate_score <= alt_tracebacks.back().second) {
                    break;
                }
                
                // we are not interested in shorter versions of alignments we've already traced,
                // so we don't consider a new traceback end if we've already seen the node inside
                // another traceback
                if (!traced_node_idxs.count(candidate)) {
                    continue;
                }
                
                // let the insert logic decide where this
                vector<Deflection> dummy_prefix;
                insert_traceback(dummy_prefix, candidate, candidate, candidate_score);
                break;
            }
        }
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
    
    inline void MultipathMEMAligner::TracebackManager::propose_deflection(const size_t from,
                                                                           const size_t to,
                                                                           const int32_t score) {
        // only propose deflections if we're going through a new untraversed section of the traceback
        if (curr_deflxn != (*curr_traceback).first.end()) {
            return;
        }
        
        // is the score good enough to be put on the stack?
        if (score <= alt_tracebacks.back().second && alt_tracebacks.size() >= max_multi_alns) {
            return;
        }
        
        insert_traceback((*curr_traceback).first, from, to, score);
    }
    
    inline size_t MultipathMEMAligner::deflect() {
        size_t to_idx = (*curr_deflxn).to_idx;
        curr_deflxn++;
        return to_idx;
    }
    
    inline void MultipathMEMAligner::insert_traceback(const vector<Deflection>& traceback_prefix,
                                                      const size_t from, const size_t to, const int32_t score) {
        // find position in stack where this should go
        auto insert_after = alt_tracebacks.rbegin();
        while (score > (*insert_after).second) {
            insert_after++;
            if (insert_after == alt_tracebacks.rend()) {
                break;
            }
        }
        
        // insert if score is high enough or stack is not full yet
        if (insert_after != alt_tracebacks.rbegin() || alt_tracebacks.size() < max_multi_alns) {
            
            // create a new traceback here
            auto new_traceback = alt_tracebacks.emplace(insert_after.base(), vector<Deflection>(), score);
            
            // add the deflections from the prefix
            (*new_traceback).first = traceback_prefix;
            
            // add the final deflection
            deflections.emplace_back(from, to);
        }
        
        // remove lowest scoring traceback if stack is over capacity
        if (alt_tracebacks.size() > max_multi_alns) {
            alt_tracebacks.pop_back();
        }
    }
    
    MultipathMEMAligner::MultipathMEMNode::MultipathMEMNode(string::const_iterator begin,
                                                            string::const_iterator end,
                                                            pos_t start_pos,
                                                            size_t offset,
                                                            int32_t score) :
        begin(begin), end(end), start_pos(start_pos), offset(offset), score(score)
    {
        // nothing to do
    }

    
}
