//
//  multipath_mapper.cpp
//  
//
//

#include "multipath_mapper.hpp"

namespace vg {
    MultipathClusterer::MultipathClusterer(const Alignment& alignment,
                                           const vector<MaximalExactMatch>& mems,
                                           const QualAdjAligner& aligner,
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
    
    void MultipathClusterer::init_mem_graph(const xg::XG& xgindex) {
        
        // TODO: handle inversions -- do I just need extra paths?
        // for now just not checking for orientation consistency and hoping it works out (but this
        // is only likely when inversion is small so that distance is not overestimated too much)
        auto distance = [&](pos_t& pos_1, pos_t& pos_2) {
            return xgindex.min_approx_path_distance(vector<string>(), id(pos_1), id(pos_2))
            - (int64_t) offset(pos_1) + (int64_t) offset(pos_2);
        };
        
        // the maximum graph distances to the right and left sides detectable from each node
        vector<pair<size_t, size_t>> maximum_detectable_gaps;
        
        int8_t match_score = aligner.match;
        
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
        
        unordered_map<pair<int64_t, int64_t>, int64_t> forward_pair_distances;
        unordered_map<pair<int64_t, int64_t>, int64_t> backward_pair_distances;
        
        for (int64_t i = 1; i < nodes.size(); i++) {
            for (int64_t j = 0; j < nodes.size(); j++) {
                
                MPCNode& node_1 = nodes[i];
                MPCNode& node_2 = nodes[j];
                
                // TODO: rethink this for a stranded graph
                
                // if one MEM is contained in the other along the read, there is no way to make an edge
                // from the end of one to the interior of the other by chopping off an overlap
                if ((node_2.begin >= node_1.begin && node_2.end <= node_1.end)
                    || (node_1.begin >= node_2.begin && node_1.end <= node_2.end)) {
                    continue;
                }
                
                // TODO: make a signed distance function in XG
                
                int64_t start_pos_dist = distance(node_1.start_pos, node_2.start_pos);
                // note: if the start of the first node is contained inside the second node match, this will
                // usually create a negative distance, which will pass the filter
                
                // look for pairs under the minimum distance in the forward orientation
                if (node_1.begin <= node_2.begin) {
                    int64_t graph_dist = start_pos_dist - (node_1.end - node_1.begin);
                    
                    // the minimum of the max detectable gap is the max gap here
                    int64_t max_gap_length = min(maximum_detectable_gaps[i].second,
                                                 maximum_detectable_gaps[j].first);
                    
                    // record all pairs that are under the maximum
                    if (node_2.begin - node_1.end <= max_gap_length && dist <= max_gap_length) {
                        forward_pair_distances[make_pair(i, j)] = dist;
                    }
                }
                else {
                    int64_t graph_dist = start_pos_dist - (node_2.end - node_2.begin);
                    
                    // the minimum of the max detectable gap is the max gap here
                    int64_t max_gap_length = min(maximum_detectable_gaps[i].first,
                                                 maximum_detectable_gaps[j].second);
                    
                    // record all pairs that are under the maximum
                    if (node_1.begin - node_2.end <= max_gap_length && graph_dist <= max_gap_length) {
                        forward_pair_distances[make_pair(j, i)] = graph_dist;
                    }
                }
                
                // look for pairs under the minimum distance in the reverse orientation
                if (node_1.end >= node_2.end) {
                    int64_t graph_dist = start_pos_dist - (node_1.end - node_1.begin);
                    
                    // the minimum of the max detectable gap is the max gap here
                    int64_t max_gap_length = min(maximum_detectable_gaps[i].second,
                                                 maximum_detectable_gaps[j].first);
                    
                    // record all pairs that are under the maximum
                    if (node_2.begin - node_1.end <= max_gap_length && dist <= max_gap_length) {
                        forward_pair_distances[make_pair(i, j)] = dist;
                    }
                }
                else {
                    
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



