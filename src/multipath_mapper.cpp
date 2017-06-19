//
//  multipath_mapper.cpp
//  
//
//

#include "multipath_mapper.hpp"

namespace vg {
    
    MultipathAligner::MultipathAligner(xg::XG& xgindex, LRUCache<id_t, Node>& node_cache,
                                       QualAdjAligner& qual_adj_aligner, SnarlManager* snarl_manager) :
                                       xgindex(xgindex), node_cache(node_cache), snarl_manager(snarl_manager),
                                       qual_adj_aligner(qual_adj_aligner) {
        // nothing to do
    }
    
    
    void MultipathAligner::multipath_align(const Alignment& alignment,
                                           const vector<MaximalExactMatch>& mems,
                                           list<MultipathAlignment>& multipath_alns_out,
                                           size_t max_alt_alns) {
        
        // cluster the MEMs
        MultipathClusterer clusterer(alignment, mems, qual_adj_aligner, xgindex, full_length_bonus,
                                     max_expected_dist_approx_error);
        vector<vector<pair<const MaximalExactMatch*, pos_t>>> clusters = clusterer.clusters();
        
        // subgraphs around each cluster
        vector<VG*> cluster_graphs;
        cluster_graphs.reserve(clusters.size());
        
        // we will ensure that nodes are in only one cluster, use this to record which one
        unordered_map<id_t, size_t> node_id_to_cluster;
        
        for (size_t i = 0; i < clusters.size(); i++) {
            
            // gather the parameters for subgraph extraction from the MEM hits
            
            vector<pair<const MaximalExactMatch*, pos_t>>& cluster = clusters[i];
            vector<pos_t> positions;
            vector<size_t> forward_max_dist;
            vector<size_t> backward_max_dist;
            
            positions.reserve(cluster.size());
            forward_max_dist.reserve(cluster.size());
            backward_max_dist.reserve(cluster.size());
            
            for (auto mem_hit : cluster) {
                // get the start position of the MEM
                positions.push_back(mem_hit.second);
                // search far enough away to get any hit detectable without soft clipping
                forward_max_dist.push_back(qual_adj_aligner.longest_detectable_gap(alignment, mem_hit.first->end, full_length_bonus)
                                           + (mem_hit.first->end - mem_hit.first->begin));
                backward_max_dist.push_back(qual_adj_aligner.longest_detectable_gap(alignment, mem_hit.first->begin, full_length_bonus));
            }
            
            // TODO: a progressive expansion of the subgraph if the MEM hit is already contained in
            // a cluster graph somewhere?
            
            // extract the subgraph within the search distance
            
            VG* cluster_graph = new VG();
            Graph& graph = cluster_graph->graph;
            
            // extract the protobuf Graph in place in the VG
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
                // there is no overlap with any other graph, suggesting a new unique hit
                
                cluster_graphs.push_back(cluster_graph);
                
                // now that we know we're going to save the graph, manually trigger the index building since
                // we circumvented the constructors
                cluster_graph->rebuild_indexes();
            }
            else {
                // this graph overlap at least one other graph, so we merge them into the (arbitrarily chosen) graph
                // with the minimum index in the vector
                size_t min_idx_cluster = *std::min_element(overlapping_graphs.begin(), overlapping_graphs.end());
                
                // merge in the new graph
                cluster_graphs[min_idx_cluster]->extend(graph);
                delete cluster_graph;
                
                // if this subgraph chains together multiple clusters, merge them and remove them from the list
                overlapping_graphs.erase(min_idx_cluster);
                for (size_t j : overlapping_graphs) {
                    std::swap(cluster_graphs[j], cluster_graphs.back());
                    cluster_graphs[min_idx_cluster]->extend(cluster_graphs.back()->graph);
                    delete cluster_graphs.back();
                    cluster_graphs.pop_back();
                }
                
                // relabel the cluster of any nodes that were not in the graph we merged into
                Graph& merged_graph = cluster_graphs[min_idx_cluster]->graph;
                for (size_t j = 0; j < merged_graph.node_size(); j++) {
                    node_id_to_cluster[merged_graph.node(j).id()] = min_idx_cluster;
                }
            }
        }
        
        // which MEMs are in play for which cluster?
        unordered_map<VG*, vector<pair<const MaximalExactMatch*, pos_t>>> cluster_graph_mems;
        for (const MaximalExactMatch& mem : mems) {
            for (gcsa::node_type hit : mem.nodes) {
                size_t cluster_idx = node_id_to_cluster[gcsa::Node::id(hit)];
                cluster_graph_mems[cluster_graphs[cluster_idx]].push_back(make_pair(&mem, make_pos_t(hit)));
            }
        }
        
        // record the sequence coverage of each cluster
        unordered_map<VG*, int64_t> cluster_graph_coverage;
        for (VG* vg : cluster_graphs) {
            cluster_graph_coverage[vg] = read_coverage(cluster_graph_mems[vg]);
        }
        
        // sort the cluster graphs descending by unique sequence coverage
        // TODO: figure out relationship between this and the clustering filter
        std::sort(cluster_graphs.begin(), cluster_graphs.end(), [&](VG* g1, VG* g2) {
            return cluster_graph_coverage[g1] > cluster_graph_coverage[g2];
        });
        
        // the longest path we could possibly align to (full gap and a full sequence)
        size_t target_length = qual_adj_aligner.longest_detectable_gap(alignment, full_length_bonus) + alignment.sequence().size();
        
        size_t num_alns = 0;
        for (VG* vg : cluster_graphs) {
            
            // if we have a cluster graph with small enough MEM coverage compared to the best one or we've made
            // the maximum number of alignmetns we stop producing alternate alignments
            if (cluster_graph_coverage[vg] < mem_coverage_min_ratio * cluster_graph_coverage[cluster_graphs[0]]
                || num_alns >= max_alt_alns) {
                break;
            }
            
            // convert from bidirected to directed
            unordered_map<id_t, pair<id_t, bool> > node_trans;
            VG align_graph = vg->split_strands(node_trans);
            // if necessary, convert from cyclic to acylic
            if (!align_graph.is_acyclic()) {
                unordered_map<id_t, pair<id_t, bool> > dagify_trans;
                align_graph = align_graph.dagify(target_length, // high enough that num SCCs is never a limiting factor
                                                 dagify_trans,
                                                 target_length,
                                                 0); // no maximum on size of component
                node_trans = align_graph.overlay_node_translations(dagify_trans, node_trans);
            }
            
            align_graph.sort();
            
            // create the injection translator, which maps a node in the original graph to every one of its occurrences
            // in the dagified graph
            unordered_multimap<id_t, pair<id_t, bool> > rev_trans;
            for (const auto& trans_record : node_trans) {
                rev_trans.insert(make_pair(trans_record.second.first,
                                           make_pair(trans_record.first, trans_record.second.second)));
            }
            
            vector<pair<const MaximalExactMatch*, pos_t>>& graph_mems = cluster_graph_mems[vg];
            
            // construct a graph that summarizes reachability between MEMs
            MultipathAlignmentGraph multi_aln_graph(align_graph, graph_mems, rev_trans, node_trans, snarl_manager, max_snarl_cut_size);
            
            vector<size_t> topological_order;
            multi_aln_graph.topological_sort(topological_order);
            
            // it's sometimes possible for transitive edges to survive the original construction algorithm, so remove them
            multi_aln_graph.remove_transitive_edges(topological_order);
            
            // prune this graph down the paths that have reasonably high likelihood
            multi_aln_graph.prune_to_high_scoring_paths(qual_adj_aligner, max_suboptimal_path_score_diff, topological_order);
            
            // create a new multipath alignment object and transfer over data from alignment
            multipath_alns_out.emplace_back();
            MultipathAlignment& multipath_aln = multipath_alns_out.back();
            multipath_aln.set_sequence(alignment.sequence());
            multipath_aln.set_quality(alignment.quality());
            multipath_aln.set_name(alignment.name());
            multipath_aln.set_sample_name(alignment.sample_name());
            multipath_aln.set_read_group(alignment.read_group());
            
            // add a subpath for each of the exact match nodes
            for (int64_t j = 0; j < multi_aln_graph.match_nodes.size(); j++) {
                ExactMatchNode& match_node = multi_aln_graph.match_nodes[j];
                Subpath* subpath = multipath_aln.add_subpath();
                *subpath->mutable_path() = match_node.path;
                int32_t match_score = qual_adj_aligner.score_exact_match(match_node.begin, match_node.end,
                                                                         alignment.quality().begin() + (match_node.begin - alignment.sequence().begin()));
                subpath->set_score(match_score + full_length_bonus * ((match_node.begin == alignment.sequence().begin()) +
                                                                      (match_node.end == alignment.sequence().end())));
            }
            
            // perform alignment in the intervening sections
            for (int64_t j = 0; j < multi_aln_graph.match_nodes.size(); j++) {
                ExactMatchNode& src_match_node = multi_aln_graph.match_nodes[j];
                Subpath* src_subpath = multipath_aln.mutable_subpath(j);
                
                const Path& path = src_subpath->path();
                const Mapping& final_mapping = path.mapping(path.mapping_size() - 1);
                // make a pos_t that points to the final base in the match
                pos_t src_pos = make_pos_t(final_mapping.position().node_id(),
                                           final_mapping.position().is_reverse(),
                                           mapping_from_length(final_mapping) - 1);
                
                // check whether this is a simple split MEM, in which case we can make the connection directly
                // and avoid making an intervening alignment
                if (src_match_node.edges.size() == 1) {
                    pair<size_t, size_t>& edge = src_match_node.edges[0];
                    // check whether the positions are right next to each other in both the read and the graph
                    //
                    // strictly speaking, there could be other paths between the match nodes if the ends of the paths
                    // are at node boundaries, but even then the best alignment would be from the empty sequence to
                    // the path of length 0
                    if (src_match_node.end == multi_aln_graph.match_nodes[edge.first].begin && edge.second == 0) {
                        // connect the supath directly without an intervening alignment
                        src_subpath->add_next(edge.first);
                        continue;
                    }
                }
                
                
                // the longest gap that could be detected at this position in the read
                size_t src_max_gap = qual_adj_aligner.longest_detectable_gap(alignment, src_match_node.end,
                                                                             full_length_bonus);
                                
                for (const pair<size_t, size_t>& edge : src_match_node.edges) {
                    ExactMatchNode& dest_match_node = multi_aln_graph.match_nodes[edge.first];
                    pos_t dest_pos = make_pos_t(multipath_aln.subpath(edge.first).path().mapping(0).position());
                    
                    size_t max_gap = std::min(src_max_gap,
                                              qual_adj_aligner.longest_detectable_gap(alignment, dest_match_node.begin,
                                                                                      full_length_bonus));
                    
                    // extract the graph between the matches
                    Graph connecting_graph;
                    unordered_map<id_t, id_t> connect_trans = algorithms::extract_connecting_graph(align_graph,      // DAG with split strands
                                                                                                   connecting_graph, // graph to extract into
                                                                                                   max_gap,          // longest distance necessary
                                                                                                   src_pos,          // end of earlier match
                                                                                                   dest_pos,         // beginning of later match
                                                                                                   false,            // do not extract the end positions in the matches
                                                                                                   false,            // do not bother finding all cycles (it's a DAG)
                                                                                                   true,             // remove tips
                                                                                                   true,             // only include nodes on connecting paths
                                                                                                   true);            // enforce max distance strictly
                    
                    
                    
                    
                    // transfer the substring between the matches to a new alignment
                    Alignment intervening_sequence;
                    intervening_sequence.set_sequence(alignment.sequence().substr(src_match_node.end - alignment.sequence().begin(),
                                                                                  dest_match_node.begin - src_match_node.end));
                    if (!alignment.quality().empty()) {
                        intervening_sequence.set_quality(alignment.quality().substr(src_match_node.end - alignment.sequence().begin(),
                                                                                    dest_match_node.begin - src_match_node.end));
                    }
                    
                    // TODO a better way of choosing the number of alternate alignments
                    // TODO alternate alignments restricted only to distinct node paths?
                    vector<Alignment> alts_alignments;
                    qual_adj_aligner.align_global_banded_multi(intervening_sequence, alts_alignments,
                                                               connecting_graph, num_alt_alns, band_padding, true);
                    
                    for (Alignment& connecting_alignment : alts_alignments) {
                        // create a subpath between the matches for this alignment
                        Subpath* connecting_subpath = multipath_aln.add_subpath();
                        *connecting_subpath->mutable_path() = connecting_alignment.path();
                        connecting_subpath->set_score(connecting_alignment.score());
                        
                        // add the appropriate connections
                        src_subpath->add_next(multipath_aln.subpath_size() - 1);
                        connecting_subpath->add_next(edge.first);
                        
                        // translate the path into the space of the main graph
                        translate_node_ids(*connecting_subpath->mutable_path(), connect_trans);
                        Mapping* first_mapping = connecting_subpath->mutable_path()->mutable_mapping(0);
                        if (first_mapping->position().node_id() == final_mapping.position().node_id()) {
                            first_mapping->mutable_position()->set_offset(mapping_from_length(final_mapping));
                        }
                    }
                }
            }
            
            vector<bool> is_source_node(multi_aln_graph.match_nodes.size(), true);
            for (size_t j = 0; j < multi_aln_graph.match_nodes.size(); j++) {
                ExactMatchNode& match_node = multi_aln_graph.match_nodes[j];
                if (match_node.edges.empty()) {
                    const Mapping& final_mapping = match_node.path.mapping(match_node.path.mapping_size() - 1);
                    if (match_node.end != alignment.sequence().end()) {
                        
                        Subpath* sink_subpath = multipath_aln.mutable_subpath(j);
                        
                        int64_t target_length = qual_adj_aligner.longest_detectable_gap(alignment, match_node.end) + (alignment.sequence().end() - match_node.end);
                        pos_t end_pos = final_position(match_node.path);
                        
                        Graph tail_graph;
                        unordered_map<id_t, id_t> tail_trans = algorithms::extract_extending_graph(align_graph,
                                                                                                   tail_graph,
                                                                                                   target_length,
                                                                                                   end_pos,
                                                                                                   false,         // search forward
                                                                                                   false);        // no need to preserve cycles (in a DAG)
                        
                        // get the sequence remaining in the right tail
                        Alignment right_tail_sequence;
                        right_tail_sequence.set_sequence(alignment.sequence().substr(match_node.end - alignment.sequence().begin(),
                                                                                     alignment.sequence().end() - match_node.end));
                        if (!alignment.quality().empty()) {
                            right_tail_sequence.set_quality(alignment.quality().substr(match_node.end - alignment.sequence().begin(),
                                                                                       alignment.sequence().end() - match_node.end));
                        }
                        
                        // align against the graph
                        vector<Alignment> alt_alignments;
                        qual_adj_aligner.align_pinned_multi(right_tail_sequence, alt_alignments, tail_graph, true, num_alt_alns, full_length_bonus);
                        
                        for (Alignment& tail_alignment : alt_alignments) {
                            Subpath* tail_subpath = multipath_aln.add_subpath();
                            *tail_subpath->mutable_path() = tail_alignment.path();
                            tail_subpath->set_score(tail_alignment.score());
                            
                            sink_subpath->add_next(multipath_aln.subpath_size() - 1);
                            
                            translate_node_ids(*tail_subpath->mutable_path(), tail_trans);
                            Mapping* first_mapping = tail_subpath->mutable_path()->mutable_mapping(0);
                            if (first_mapping->position().node_id() == final_mapping.position().node_id()) {
                                first_mapping->mutable_position()->set_offset(mapping_from_length(final_mapping));
                            }
                        }
                    }
                }
                else {
                    for (const pair<size_t, size_t>& edge : match_node.edges) {
                        is_source_node[edge.first] = false;
                    }
                }
            }
            
            for (size_t j = 0; j < multi_aln_graph.match_nodes.size(); j++) {
                if (is_source_node[j]) {
                    ExactMatchNode& match_node = multi_aln_graph.match_nodes[j];
                    if (match_node.begin != alignment.sequence().begin()) {
                        
                        int64_t target_length = qual_adj_aligner.longest_detectable_gap(alignment, match_node.begin) + (match_node.begin - alignment.sequence().begin());
                        pos_t begin_pos = initial_position(match_node.path);
                        
                        Graph tail_graph;
                        unordered_map<id_t, id_t> tail_trans = algorithms::extract_extending_graph(align_graph,
                                                                                                   tail_graph,
                                                                                                   target_length,
                                                                                                   begin_pos,
                                                                                                   true,          // search backward
                                                                                                   false);        // no need to preserve cycles (in a DAG)
                        
                        
                        Alignment left_tail_sequence;
                        left_tail_sequence.set_sequence(alignment.sequence().substr(0, match_node.begin - alignment.sequence().begin()));
                        if (!alignment.quality().empty()) {
                            left_tail_sequence.set_quality(alignment.quality().substr(0, match_node.begin - alignment.sequence().begin()));
                        }
                        
                        vector<Alignment> alt_alignments;
                        qual_adj_aligner.align_pinned_multi(left_tail_sequence, alt_alignments, tail_graph, false, num_alt_alns, full_length_bonus);
                        
                        for (Alignment& tail_alignment : alt_alignments) {
                            Subpath* tail_subpath = multipath_aln.add_subpath();
                            *tail_subpath->mutable_path() = tail_alignment.path();
                            tail_subpath->set_score(tail_alignment.score());
                            
                            tail_subpath->add_next(j);
                            multipath_aln.add_start(multipath_aln.subpath_size() - 1);
                            
                            translate_node_ids(*tail_subpath->mutable_path(), tail_trans);
                        }
                    }
                    else {
                        multipath_aln.add_start(j);
                    }
                }
            }
            
            for (size_t j = 0; j < multipath_aln.subpath_size(); j++) {
                translate_oriented_node_ids(*multipath_aln.mutable_subpath(j)->mutable_path(), node_trans);
            }
            
            num_alns++;
        }
        
        for (VG* vg : cluster_graphs) {
            delete vg;
        }
    }
    
    int64_t MultipathAligner::read_coverage(const vector<pair<const MaximalExactMatch*, pos_t>>& mem_hits) {
        if (mem_hits.empty()) {
            return 0;
        }
        
        vector<pair<string::const_iterator, string::const_iterator>> mem_read_segments;
        mem_read_segments.reserve(mem_hits.size());
        for (auto& mem_hit : mem_hits) {
            mem_read_segments.emplace_back(mem_hit.first->begin, mem_hit.first->end);
        }
        std::sort(mem_read_segments.begin(), mem_read_segments.end());
        auto curr_begin = mem_read_segments[0].first;
        auto curr_end = mem_read_segments[0].second;
        
        int64_t total = 0;
        for (size_t i = 1; i < mem_read_segments.size(); i++) {
            if (mem_read_segments[i].first >= curr_end) {
                total += (curr_end - curr_begin);
                curr_begin = mem_read_segments[i].first;
            }
            curr_end = mem_read_segments[i].second;
        }
        return total + (curr_end - curr_begin);
    }
    
    double MultipathAligner::read_coverage_z_score(int64_t coverage, const Alignment& alignment) {
        /* algebraically equivalent to
         *
         *      Coverage - ReadLen / 4
         *  -------------------------------
         *  sqrt(ReadLen * 1/4 * (1 - 1/4))
         * 
         * from the Normal approximation to a Binomal(ReadLen, 1/4)
         */
        double root_len = sqrt(alignment.sequence().size());
        return 0.5773502691896258 * (4.0 * coverage / root_len - root_len);
    }
    
    MultipathAlignmentGraph::MultipathAlignmentGraph(VG& vg, const vector<pair<const MaximalExactMatch*, pos_t>>& hits,
                                                     const unordered_multimap<id_t, pair<id_t, bool>>& injection_trans,
                                                     const unordered_map<id_t, pair<id_t, bool>>& projection_trans,
                                                     SnarlManager* cutting_snarls, int64_t max_snarl_cut_size) {

        
        // map of node ids in the dagified graph to the indices in the matches that contain them
        unordered_map<int64_t, vector<int64_t>> node_matches;
        
        // walk the matches and filter out redundant sub-MEMs
        for (int64_t i = 0; i < hits.size(); i++) {
            
            const pair<const MaximalExactMatch*, pos_t>& hit = hits[i];
            
            // the part of the read we're going to match
            string::const_iterator begin = hit.first->begin;
            string::const_iterator end = hit.first->end;
            int64_t mem_length = end - begin;
            // the start of the hit in the original graph
            const pos_t& hit_pos = hit.second;
            
            auto hit_range = injection_trans.equal_range(id(hit_pos));
            for (auto iter = hit_range.first; iter != hit_range.second; iter++) {
                // this graph is unrolled/dagified, so all orientations should match
                if ((*iter).second.second != is_rev(hit_pos)) {
                    continue;
                }
                
                // an id that corresponds to the original node
                id_t injected_id = (*iter).first;
                
                // check all MEMs that traversed this node to see if this is a redundant sub-MEM
                bool is_partial_mem = false;
                for (int64_t j : node_matches[id(hit_pos)]) {
                    ExactMatchNode& match_node = match_nodes[j];
                    
                    int64_t relative_offset = begin - match_node.begin;
                    if (relative_offset < 0 || relative_offset + (end - begin) >= match_node.end - match_node.begin) {
                        // the hit does not fall on the same section of the read as the other match, so
                        // it cannot be contained in it
                        continue;
                    }
                    
                    Path& path = match_node.path;
                    
                    // if this is a partial MEM, we should be able to predict its hit location by traversing the path
                    // of the parent MEM by a distance equal to the relative offset
                    
                    int64_t prefix_length = 0;
                    for (size_t k = 0; k < path.mapping_size(); k++) {
                        if (prefix_length > relative_offset) {
                            break;
                        }
                        const Mapping& mapping = path.mapping(k);
                        // the length through this mapping
                        int64_t prefix_through_length = prefix_length + mapping_from_length(mapping);
                        if (prefix_through_length > relative_offset) {
                            // we cross the relative offset on this node, so check if the path is in the predicted
                            // position for a redundant sub-MEM
                            id_t node_id_here = mapping.position().node_id();
                            is_partial_mem = is_partial_mem ||
                                             (injected_id == node_id_here
                                              && prefix_length + offset(hit_pos) == relative_offset
                                              && projection_trans.at(node_id_here).second == is_rev(hit_pos));
                        }
                        prefix_length = prefix_through_length;
                    }
                    if (is_partial_mem) {
                        break;
                    }
                }
                
                // don't walk the match of false partial hits
                if (is_partial_mem) {
                    continue;
                }
                
                // stack for DFS, each record contains tuples of (read begin, node offset, next node index, next node ids)
                vector<tuple<string::const_iterator, size_t, size_t, vector<NodeTraversal>>> stack;
                stack.emplace_back(begin, offset(hit_pos), 0,
                                   vector<NodeTraversal>{NodeTraversal(vg.get_node(injected_id))});
                
                while (!stack.empty()) {
                    auto& back = stack.back();
                    if (get<2>(back) == get<3>(back).size()) {
                        stack.pop_back();
                    }
                    NodeTraversal trav = get<3>(back)[get<2>(back)];
                    
                    const string& node_seq = trav.node->sequence();
                    size_t node_idx = get<1>(back);
                    string::const_iterator read_iter = get<0>(back);
                    
                    // look for a match along the entire node sequence
                    for (; node_idx < node_seq.size() && read_iter != end; node_idx++, read_iter++) {
                        if (node_seq[node_idx] != *read_iter) {
                            break;
                        }
                    }
                    
                    if (read_iter == end) {
                        // finished walking match
                        break;
                    }
                    else if (node_idx == node_seq.size()) {
                        // matched entire node
                        stack.emplace_back(read_iter, 0, 0, vector<NodeTraversal>());
                        vg.nodes_next(trav, get<3>(stack.back()));
                    }
                    else {
                        // match did not complete node or finish, this is a miss
                        get<2>(back)++;
                    }
                }
                
                if (!stack.empty()) {
                    // we left the trace in the stack, which means we found a complete match
                    int64_t match_node_idx = match_nodes.size();
                    match_nodes.emplace_back();
                    ExactMatchNode& match_node = match_nodes.back();
                    Path& path = match_node.path;
                    match_node.begin = begin;
                    match_node.end = end;
                    int64_t length_remaining = end - begin;
                    
                    // walk out the match
                    int32_t rank = 1;
                    for (auto search_record : stack) {
                        int64_t offset = get<1>(search_record);
                        Node* node = get<3>(search_record)[get<2>(search_record)].node;
                        int64_t length = std::min((int64_t) node->sequence().size() - offset, length_remaining);
                        
                        Mapping* mapping = path.add_mapping();
                        mapping->set_rank(rank);
                        
                        Edit* edit = mapping->add_edit();
                        edit->set_from_length(length);
                        edit->set_to_length(length);
                        
                        // note: the graph is dagified and unrolled, so all hits should be on the forward strand
                        Position* position = mapping->mutable_position();
                        position->set_node_id(node->id());
                        position->set_offset(offset);
                        
                        // record that each node occurs in this match so we can filter out sub-MEMs
                        node_matches[node->id()].push_back(match_node_idx);
                        
                        rank++;
                        length_remaining -= length;
                    }
                    
                    if (cutting_snarls) {
                        // we indicated a snarl manager that owns the snarls we want to cut out of exact matches
                        
                        // first compute the segments we want to cut out
                        
                        // this list holds the beginning of the current segment at each depth in the snarl hierarchy
                        // as we traverse the exact match, beginning is recorded in both sequence distance and node index
                        list<pair<size_t, size_t>> level_segment_begin;
                        level_segment_begin.emplace_back(0, 0);
                        
                        // we record which segments we are going to cut out of the match here
                        vector<pair<size_t, size_t>> cut_segments;
                        
                        auto curr_level = level_segment_begin.begin();
                        size_t prefix_length = 0;
                        for (size_t j = 0, last = path.mapping_size() - 1; j <= last; j++) {
                            
                            const Position& position = path.mapping(j).position();
                            
                            if (j > 0) {
                                // we have entered this node on this iteration
                                if (cutting_snarls->into_which_snarl(position.node_id(), !position.is_reverse())) {
                                    // as we enter this node, we are leaving the snarl we were in
                                    
                                    // since we're going up a level, we need to check whether we need to cut out the segment we've traversed
                                    if (prefix_length - (*curr_level).first <= max_snarl_cut_size) {
                                        cut_segments.emplace_back((*curr_level).second, j);
                                    }
                                    
                                    curr_level++;
                                    if (curr_level == level_segment_begin.end()) {
                                        // we were already at the highest level seen so far, so we need to add a new one
                                        // the entire previous part of the match is contained in this level, so we start
                                        // the segment from 0
                                        curr_level = level_segment_begin.insert(level_segment_begin.end(), make_pair(0, 0));
                                    }
                                }
                            }
                            
                            // cross to the other side of the node
                            prefix_length += mapping_from_length(path.mapping(j));
                            
                            if (j < last) {
                                // we are going to leave this node next iteration
                                if (cutting_snarls->into_which_snarl(position.node_id(), position.is_reverse())) {
                                    // as we leave this node, we are entering a new deeper snarl
                                    
                                    // the segment in the new level will begin at the end of the current node
                                    if (curr_level == level_segment_begin.begin()) {
                                        // we are already at the lowest level seen so far, so we need to add a new one
                                        level_segment_begin.emplace_front(prefix_length, j + 1);
                                        curr_level--;
                                    }
                                    else {
                                        // the lower level is in the record already, so we update its segment start
                                        curr_level--;
                                        *curr_level = make_pair(prefix_length, j + 1);
                                    }
                                }
                            }
                        }
                        
                        // check the final segment for a cut unless we're at the highest level in the match
                        auto last = level_segment_begin.end();
                        last--;
                        if (prefix_length - (*curr_level).first <= max_snarl_cut_size && curr_level != last) {
                            cut_segments.emplace_back((*curr_level).second, path.mapping_size());
                        }
                        
                        // did we cut out any segments?
                        if (!cut_segments.empty()) {
                            
                            // we may have decided to cut the segments of both a parent and child snarl, so now we
                            // collapse the list of intervals, which is sorted on the end index by construction
                            //
                            // snarl nesting properties guarantee that there will be at least one node between any
                            // cut segments that are not nested, so we don't need to deal with the case where the
                            // segments are partially overlapping (i.e. it's a bit easier than the general interval
                            // intersection problem)
                            vector<pair<size_t, size_t>> keep_segments;
                            size_t curr_keep_seg_end = path.mapping_size();
                            auto riter = cut_segments.rbegin();
                            if ((*riter).second == curr_keep_seg_end) {
                                // don't add an empty keep segment in the first position
                                curr_keep_seg_end = (*riter).first;
                                riter++;
                            }
                            for (; riter != cut_segments.rend(); riter++) {
                                if ((*riter).second < curr_keep_seg_end) {
                                    // this is a new interval
                                    keep_segments.emplace_back((*riter).second, curr_keep_seg_end);
                                    curr_keep_seg_end = (*riter).first;
                                }
                            }
                            if (curr_keep_seg_end > 0) {
                                // we are not cutting off the left tail, so add a keep segment for it
                                keep_segments.emplace_back(0, curr_keep_seg_end);
                            }
                            
                            // make a new node for all but one of the keep segments
                            auto last = keep_segments.end() - 1;
                            for (auto iter = keep_segments.begin(); iter != last; iter++) {
                                match_nodes.emplace_back();
                                ExactMatchNode& cut_node = match_nodes.back();
                                Path& cut_path = cut_node.path;
                                // transfer over the keep segment from the main path
                                int32_t rank = 1;
                                for (size_t j = (*iter).first; j < (*iter).second; j++, rank++) {
                                    Mapping* mapping = cut_path.add_mapping();
                                    *mapping = path.mapping(j);
                                    mapping->set_rank(rank);
                                }
                            }
                            // replace the path of the original node with the final keep segment
                            Path new_path;
                            int32_t rank = 1;
                            for (size_t j = (*last).first; j < (*last).second; j++, rank++) {
                                Mapping* mapping = new_path.add_mapping();
                                *mapping = path.mapping(j);
                                mapping->set_rank(rank);
                            }
                            path = new_path;
                        }
                    }
                }
            }
        }
        
        // now we calculate reachability between the walked matches so we know which ones
        // to connect with intervening alignments
        
        auto start_offset = [&](size_t idx) {
            return match_nodes[idx].path.mapping(0).position().offset();
        };
        
        auto end_offset = [&](size_t idx) {
            Path& path = match_nodes[idx].path;
            const Mapping& mapping = path.mapping(path.mapping_size() - 1);
            return mapping.position().offset() + mapping_from_length(mapping);
        };
        
        auto start_node_id = [&](size_t idx) {
            return match_nodes[idx].path.mapping(0).position().node_id();
        };
        
        auto end_node_id = [&](size_t idx) {
            Path& path = match_nodes[idx].path;
            return path.mapping(path.mapping_size() - 1).position().node_id();
        };

        // record the start and end node ids of every exact match
        unordered_map<id_t, vector<size_t>> exact_match_starts;
        unordered_map<id_t, vector<size_t>> exact_match_ends;
        for (size_t i = 0; i < match_nodes.size(); i++) {
            Path& path = match_nodes[i].path;
            exact_match_starts[path.mapping(0).position().node_id()].push_back(i);
            exact_match_ends[path.mapping(path.mapping_size() - 1).position().node_id()].push_back(i);
        }
        
        // sort the MEMs starting and ending on each node in node sequence order
        for (pair<const id_t, vector<size_t>>& node_match_starts : exact_match_starts) {
            std::sort(node_match_starts.second.begin(), node_match_starts.second.end(),
                      [&](const size_t idx_1, const size_t idx_2) {
                          return start_offset(idx_1) < start_offset(idx_2);
                      });
        }
        for (pair<const id_t, vector<size_t>>& node_match_ends : exact_match_ends) {
            std::sort(node_match_ends.second.begin(), node_match_ends.second.end(),
                      [&](const size_t idx_1, const size_t idx_2) {
                          return end_offset(idx_1) < end_offset(idx_2);
                      });
        }
        
        // some structures we will fill out with DP:
        
        // for each node, the starts and ends of MEMs that can reach this node without passing any other
        // MEM starts or ends
        unordered_map<id_t, unordered_map<size_t, size_t>> reachable_ends;
        unordered_map<id_t, unordered_map<size_t, size_t>> reachable_starts;
        
        // for the start of each MEM, the starts and ends of other MEMs that can reach it without passing any
        // other start or end
        unordered_map<size_t, vector<pair<size_t, size_t>>> reachable_starts_from_start;
        unordered_map<size_t, vector<pair<size_t, size_t>>> reachable_ends_from_start;
        
        // for the end of each MEM, the ends of other MEMs that can reach it without passing any
        // other start or end
        unordered_map<size_t, vector<pair<size_t, size_t>>> reachable_ends_from_end;
        unordered_map<size_t, vector<pair<size_t, size_t>>> reachable_starts_from_end;
        
        // note: graph has been sorted into topological order
        
        Graph& graph = vg.graph;
        for (int64_t i = 0; i < graph.node_size(); i++) {
            Node* node = graph.mutable_node(i);
            id_t node_id = node->id();
            
            size_t node_length = node->sequence().size();
            
            // do any MEMs start or end on this node?
            bool contains_starts = exact_match_starts.count(node_id);
            bool contains_ends = exact_match_ends.count(node_id);
            
            // we will use DP to carry reachability information forward onto the next nodes
            vector<NodeTraversal> nexts;
            vg.nodes_next(NodeTraversal(node), nexts);
            
            if (contains_starts && contains_ends) {
                // since there are both starts and ends on this node, we have to traverse both lists simultaneously
                // to assess reachability within the same node
                
                vector<size_t>& ends = exact_match_starts[node_id];
                vector<size_t>& starts = exact_match_starts[node_id];
                
                // find the range of starts and ends in the list with the same offset
                
                size_t start_range_begin = 0;
                size_t start_range_end = 0;
                size_t end_range_begin = 0;
                size_t end_range_end = 0;
                
                size_t curr_start_offset = start_offset(starts[start_range_begin]);
                size_t curr_end_offset = end_offset(ends[end_range_begin]);
                size_t prev_offset = 0;
                
                bool prev_is_end = curr_end_offset <= curr_start_offset;
                while (end_offset(starts[start_range_end]) == curr_end_offset) {
                    end_range_end++;
                    if (end_range_end == ends.size()) {
                        break;
                    }
                }
                while (start_offset(starts[start_range_end]) == curr_start_offset) {
                    start_range_end++;
                    if (start_range_end == starts.size()) {
                        break;
                    }
                }
                
                // connect the first range of starts or ends to the incoming starts and ends
                
                size_t prev_end_range_begin = end_range_begin;
                size_t prev_start_range_begin = start_range_begin;
                if (prev_is_end) {
                    for (size_t j = end_range_begin; j < end_range_end; j++) {
                        for (const pair<size_t, size_t>& incoming_end : reachable_ends[node_id]) {
                            reachable_ends_from_end[ends[j]].push_back(incoming_end);
                        }
                        for (const pair<size_t, size_t>& incoming_start : reachable_starts[node_id]) {
                            reachable_starts_from_end[ends[j]].push_back(incoming_start);
                        }
                    }
                    
                    end_range_begin = end_range_end;
                    prev_offset = curr_end_offset;
                    curr_end_offset = end_offset(ends[end_range_begin]);
                    if (end_range_begin != ends.size()) {
                        while (end_offset(starts[start_range_end]) == curr_end_offset) {
                            end_range_end++;
                            if (end_range_end == ends.size()) {
                                break;
                            }
                        }
                    }
                }
                else {
                    for (size_t j = start_range_begin; j < start_range_end; j++) {
                        for (const pair<size_t, size_t>& incoming_end : reachable_ends[node_id]) {
                            reachable_ends_from_start[starts[j]].push_back(incoming_end);
                        }
                        for (const pair<size_t, size_t>& incoming_start : reachable_starts[node_id]) {
                            reachable_starts_from_start[starts[j]].push_back(incoming_start);
                        }
                    }
                    
                    start_range_begin = start_range_end;
                    prev_offset = curr_start_offset;
                    curr_start_offset = start_offset(starts[start_range_begin]);
                    if (start_range_begin != starts.size()) {
                        while (start_offset(starts[start_range_end]) == curr_start_offset) {
                            start_range_end++;
                            if (start_range_end == starts.size()) {
                                break;
                            }
                        }
                    }
                }
                
                // iterate along ranges of starts or ends in order of their offsets
                
                while (start_range_begin < starts.size() && end_range_begin < ends.size()) {
                    if (curr_end_offset <= curr_start_offset) {
                        
                        size_t dist_between = curr_end_offset - prev_offset;
                        
                        // connect this range to the previous range
                        if (prev_is_end) {
                            for (size_t j = prev_end_range_begin; j < end_range_begin; j++) {
                                for (size_t k = end_range_begin; k < end_range_end; k++) {
                                    reachable_ends_from_end[ends[k]].push_back(make_pair(ends[j], dist_between));
                                }
                            }
                        }
                        else {
                            for (size_t j = prev_start_range_begin; j < start_range_begin; j++) {
                                for (size_t k = end_range_begin; k < end_range_end; k++) {
                                    reachable_starts_from_end[ends[k]].push_back(make_pair(starts[j], dist_between));
                                }
                            }
                        }
                        
                        // record the properties of this range
                        prev_end_range_begin = end_range_begin;
                        prev_is_end = true;
                        
                        // advance to the next range
                        end_range_begin = end_range_end;
                        prev_offset = curr_end_offset;
                        curr_end_offset = end_offset(ends[end_range_begin]);
                        if (end_range_begin != ends.size()) {
                            while (end_offset(starts[start_range_end]) == curr_end_offset) {
                                end_range_end++;
                                if (end_range_end == ends.size()) {
                                    break;
                                }
                            }
                        }
                    }
                    else {
                        
                        size_t dist_between = curr_end_offset - prev_offset;
                        
                        // connect this range to the previous range
                        if (prev_is_end) {
                            for (size_t j = prev_end_range_begin; j < end_range_begin; j++) {
                                for (size_t k = start_range_begin; k < start_range_end; k++) {
                                    reachable_ends_from_start[starts[k]].push_back(make_pair(ends[j], dist_between));
                                }
                            }
                        }
                        else {
                            for (size_t j = prev_start_range_begin; j < start_range_begin; j++) {
                                for (size_t k = start_range_begin; k < start_range_end; k++) {
                                    reachable_starts_from_start[starts[k]].push_back(make_pair(starts[j], dist_between));
                                }
                            }
                        }
                        
                        // record the properties of this range
                        prev_start_range_begin = start_range_begin;
                        prev_is_end = false;
                        
                        // advance to the next range
                        start_range_begin = start_range_end;
                        prev_offset = curr_start_offset;
                        curr_start_offset = start_offset(starts[start_range_begin]);
                        if (start_range_begin != starts.size()) {
                            while (start_offset(starts[start_range_end]) == curr_start_offset) {
                                start_range_end++;
                                if (start_range_end == starts.size()) {
                                    break;
                                }
                            }
                        }
                    }
                }
                
                // finish off the list of starts on this node
                
                while (start_range_begin < starts.size()) {
                    
                    size_t dist_between = curr_start_offset - prev_offset;
                    
                    if (prev_is_end) {
                        for (size_t j = prev_end_range_begin; j < end_range_begin; j++) {
                            for (size_t k = start_range_begin; k < start_range_end; k++) {
                                reachable_ends_from_start[starts[k]].push_back(make_pair(ends[j], dist_between));
                            }
                        }
                    }
                    else {
                        for (size_t j = prev_start_range_begin; j < start_range_begin; j++) {
                            for (size_t k = start_range_begin; k < start_range_end; k++) {
                                reachable_starts_from_start[starts[k]].push_back(make_pair(starts[j], dist_between));
                            }
                        }
                    }
                    
                    prev_start_range_begin = start_range_begin;
                    start_range_begin = start_range_end;
                    prev_offset = curr_start_offset;
                    curr_start_offset = start_offset(starts[start_range_begin]);
                    prev_is_end = false;
                    
                    if (start_range_begin != starts.size()) {
                        while (start_offset(starts[start_range_end]) == curr_start_offset) {
                            start_range_end++;
                            if (start_range_end == starts.size()) {
                                break;
                            }
                        }
                    }
                }
                
                // finish off the list of ends on this node
                
                while (end_range_begin < ends.size()) {
                    
                    size_t dist_between = curr_start_offset - prev_offset;
                    
                    if (prev_is_end) {
                        for (size_t j = prev_end_range_begin; j < end_range_begin; j++) {
                            for (size_t k = end_range_begin; k < end_range_end; k++) {
                                reachable_ends_from_end[ends[k]].push_back(make_pair(ends[j], dist_between));
                            }
                        }
                    }
                    else {
                        for (size_t j = prev_start_range_begin; j < start_range_begin; j++) {
                            for (size_t k = end_range_begin; k < end_range_end; k++) {
                                reachable_starts_from_end[ends[k]].push_back(make_pair(starts[j], dist_between));
                            }
                        }
                    }
                    
                    prev_end_range_begin = end_range_begin;
                    end_range_begin = end_range_end;
                    prev_offset = curr_end_offset;
                    curr_end_offset = end_offset(ends[end_range_begin]);
                    prev_is_end = true;
                    
                    if (end_range_begin != ends.size()) {
                        while (end_offset(starts[start_range_end]) == curr_end_offset) {
                            end_range_end++;
                            if (end_range_end == ends.size()) {
                                break;
                            }
                        }
                    }
                }
                
                // carry forward the reachability of the last range onto the next nodes
                size_t dist_thru = node_length - std::max(curr_end_offset, curr_start_offset);
                
                if (prev_is_end) {
                    for (NodeTraversal next : nexts) {
                        unordered_map<size_t, size_t> reachable_ends_next = reachable_ends[next.node->id()];
                        for (size_t j = prev_end_range_begin; j < ends.size(); j++) {
                            if (reachable_ends_next.count(ends[j])) {
                                reachable_ends_next[ends[j]] = std::min(reachable_ends_next[ends[j]], dist_thru);
                            }
                            else {
                                reachable_ends_next[ends[j]] = dist_thru;
                            }
                        }
                    }
                }
                else {
                    for (NodeTraversal next : nexts) {
                        unordered_map<size_t, size_t> reachable_starts_next = reachable_starts[next.node->id()];
                        for (size_t j = prev_start_range_begin; j < starts.size(); j++) {
                            if (reachable_starts_next.count(ends[j])) {
                                reachable_starts_next[starts[j]] = std::min(reachable_starts_next[starts[j]], dist_thru);
                            }
                            else {
                                reachable_starts_next[starts[j]] = dist_thru;
                            }
                        }
                    }
                }
            }
            else if (contains_starts) {
                
                // record the incoming starts/ends for the starts on this node
                
                vector<size_t>& starts = exact_match_starts[node_id];
                
                // the starts/ends coming into this node from outside
                size_t range_begin = 0;
                size_t range_end = 0;
                size_t curr_offset = start_offset(starts[range_begin]);
                size_t prev_offset = curr_offset;
                // find the range of starts that are at the first offset
                while (start_offset(starts[range_end]) == curr_offset) {
                    range_end++;
                    if (range_end == starts.size()) {
                        break;
                    }
                }
                // connect the range to the incoming starts/ends
                for (size_t j = range_begin; j < range_end; j++) {
                    for (const pair<size_t, size_t>& incoming_start : reachable_starts[node_id]) {
                        reachable_starts_from_start[starts[j]].push_back(incoming_start);
                    }
                    for (const pair<size_t, size_t>& incoming_end : reachable_ends[node_id]) {
                        reachable_ends_from_start[starts[j]].push_back(incoming_end);
                    }
                }
                
                // the reachable starts internal to this node
                size_t prev_range_begin = 0;
                range_begin = range_end;
                while (range_begin < starts.size()) {
                    // find the range of starts at this offset
                    prev_offset = curr_offset;
                    curr_offset = start_offset(starts[range_begin]);
                    while (start_offset(starts[range_end]) == curr_offset) {
                        range_end++;
                        if (range_end == starts.size()) {
                            break;
                        }
                    }
                    
                    size_t dist_between = curr_offset - prev_offset;
                    
                    // connect this range to the previous range
                    for (size_t j = range_begin; j < range_end; j++) {
                        for (size_t k = prev_range_begin; k < range_begin; k++) {
                            reachable_starts_from_start[starts[j]].push_back(make_pair(starts[k], dist_between));
                        }
                    }
                    prev_range_begin = range_begin;
                    range_begin = range_end;
                }
                
                // this node contains at least one start of a MEM, so carry forward the reachability of all
                // starts at the final offset onto the next nodes
                
                size_t dist_thru = node_length - curr_offset;
                
                for (NodeTraversal next : nexts) {
                    unordered_map<size_t, size_t> reachable_starts_next = reachable_starts[next.node->id()];
                    for (size_t j = prev_range_begin; j < starts.size(); j++) {
                        if (reachable_starts_next.count(starts[j])) {
                            reachable_starts_next[starts[j]] = std::min(reachable_starts_next[starts[j]], dist_thru);
                        }
                        else {
                            reachable_starts_next[starts[j]] = dist_thru;
                        }
                    }
                }
            }
            else if (contains_ends) {
                
                // record the incoming starts/ends for the ends on this node
                
                vector<size_t>& ends = exact_match_ends[node_id];
                
                // the ends coming into this node from outside
                size_t range_begin = 0;
                size_t range_end = 0;
                size_t curr_offset = end_offset(ends[range_begin]);
                size_t prev_offset = curr_offset;
                // find the range of ends that are at the first offset
                while (end_offset(ends[range_end]) == curr_offset) {
                    range_end++;
                    if (range_end == ends.size()) {
                        break;
                    }
                }
                // connect the range to the incoming ends
                for (size_t j = range_begin; j < range_end; j++) {
                    for (const pair<size_t, size_t>& incoming_end : reachable_ends[node_id]) {
                        reachable_ends_from_end[ends[j]].push_back(incoming_end);
                    }
                    for (const pair<size_t, size_t>& incoming_start : reachable_starts[node_id]) {
                        reachable_starts_from_end[ends[j]].push_back(incoming_start);
                    }
                }
                
                // the reachable ends internal to this node
                size_t prev_range_begin = range_begin;
                range_begin = range_end;
                while (range_begin < ends.size()) {
                    // find the range of ends at this offset
                    prev_offset = curr_offset;
                    curr_offset = end_offset(ends[range_begin]);
                    while (end_offset(ends[range_end]) == curr_offset) {
                        range_end++;
                        if (range_end == ends.size()) {
                            break;
                        }
                    }
                    
                    size_t dist_between = curr_offset - prev_offset;
                    
                    // connect this range to the previous range
                    for (size_t j = range_begin; j < range_end; j++) {
                        for (size_t k = prev_range_begin; k < range_begin; k++) {
                            reachable_ends_from_end[ends[j]].push_back(make_pair(ends[k], dist_between));
                        }
                    }
                    prev_range_begin = range_begin;
                    range_begin = range_end;
                }
                
                // this node contains at least one end of a MEM, so carry forward the reachability of all
                // ends at the final offset onto the next nodes
                
                size_t dist_thru = node_length - curr_offset;
                
                for (NodeTraversal next : nexts) {
                    unordered_map<size_t, size_t> reachable_ends_next = reachable_ends[next.node->id()];
                    for (size_t j = prev_range_begin; j < ends.size(); j++) {
                        if (reachable_ends_next.count(ends[j])) {
                            reachable_ends_next[ends[j]] = std::min(reachable_ends_next[ends[j]], dist_thru);
                        }
                        else {
                            reachable_ends_next[ends[j]] = dist_thru;
                        }
                    }
                }
            }
            else {
                // this doesn't contain the start or end of any MEM, so we carry forward the reachability
                // into this node onto the next nodes
                
                for (NodeTraversal next : nexts) {
                    unordered_map<size_t, size_t>& reachable_ends_next = reachable_ends[next.node->id()];
                    for (const pair<size_t, size_t>& reachable_end : reachable_ends[node_id]) {
                        size_t dist_thru = reachable_end.second + node_length;
                        if (reachable_ends_next.count(reachable_end.first)) {
                            reachable_ends_next[reachable_end.first] = std::min(reachable_ends_next[reachable_end.first],
                                                                                dist_thru);
                        }
                        else {
                            reachable_ends_next[reachable_end.first] = dist_thru;
                        }
                    }
                    
                    unordered_map<size_t, size_t>& reachable_starts_next = reachable_starts[next.node->id()];
                    for (const pair<size_t, size_t>& reachable_start : reachable_starts[node_id]) {
                        size_t dist_thru = reachable_start.second + node_length;
                        if (reachable_starts_next.count(reachable_start.first)) {
                            reachable_starts_next[reachable_start.first] = std::min(reachable_starts_next[reachable_start.first],
                                                                                    dist_thru);
                        }
                        else {
                            reachable_starts_next[reachable_start.first] = dist_thru;
                        }
                    }
                }
            }
        }
        
        // now we have the reachability information for the start and end of every MEM in the graph. we
        // will use this to navigate between the MEMs in a way that respects graph reachability so that this
        // phase of the algorithm only needs to pay attention to read colinearity and transitive reducibility
        
        vector<unordered_map<size_t, size_t>> noncolinear_shells(match_nodes.size());
        
        // tuples of (overlap size, index from, index onto, dist)
        vector<tuple<size_t, size_t, size_t, size_t>> confirmed_overlaps;
        
        for (size_t i = 0; i < graph.node_size(); i++) {
            id_t node_id = graph.node(i).id();
            
            if (!exact_match_starts.count(node_id)) {
                continue;
            }
            
            for (size_t start : exact_match_starts[node_id]) {
                // traverse all of the reachable starts to find the adjacent ends that might be colinear
                ExactMatchNode& start_node = match_nodes[start];
                
                unordered_map<size_t, size_t>& noncolinear_shell = noncolinear_shells[start];
                
                // pairs of (dist, index)
                priority_queue<pair<size_t, size_t>, vector<pair<size_t, size_t>>, std::greater<pair<size_t, size_t>>> start_queue;
                priority_queue<pair<size_t, size_t>, vector<pair<size_t, size_t>>, std::greater<pair<size_t, size_t>>> end_queue;
                start_queue.emplace(start, 0);
                
                unordered_set<size_t> traversed_start;
                
                while (!start_queue.empty()) {
                    pair<size_t, size_t> start_here = start_queue.top();
                    start_queue.pop();
                    if (traversed_start.count(start_here.first)) {
                        continue;
                    }
                    traversed_start.insert(start_here.second);
                    
                    // the minimum distance to each of the starts or ends this can reach is the sum of the min distance
                    // between them and the distance already traversed
                    for (const pair<size_t, size_t>& end : reachable_ends_from_start[start_here.second]) {
                        end_queue.emplace(start_here.first + end.second, end.first);
                    }
                    
                    for (const pair<size_t, size_t>& start_next : reachable_starts_from_start[start_here.second]) {
                        start_queue.emplace(start_here.first + start_next.second, start_here.second);
                    }
                }
                
                // now we've traversed all of the starts, we have the set of ends that can be reached
                // without passing another end
                
                unordered_set<size_t> traversed_end;
                
                while (!end_queue.empty()) {
                    size_t candidate_end = end_queue.top().second;
                    size_t candidate_dist = end_queue.top().first;
                    end_queue.pop();
                    if (traversed_end.count(candidate_end)) {
                        continue;
                    }
                    traversed_end.insert(candidate_end);
                    
                    ExactMatchNode& candidate_end_node = match_nodes[candidate_end];
                    
                    if (candidate_end_node.end <= start_node.begin) {
                        // these MEMs are read colinear and graph reachable, so connect them
                        candidate_end_node.edges.push_back(make_pair(start, candidate_dist));
                        
                        // skip to the predecessor's noncolinear shell, whose connections might not be blocked by
                        // this connection
                        for (const pair<size_t, size_t>& shell_pred : noncolinear_shells[candidate_end]) {
                            end_queue.emplace(candidate_dist + shell_pred.second, shell_pred.first);
                        }
                    }
                    else if (start_node.end > candidate_end_node.end && candidate_end_node.begin < start_node.begin) {
                        // the MEM can be made colinear by removing an overlap, which will not threaten reachability
                        size_t overlap = candidate_end_node.end - start_node.begin;
                        confirmed_overlaps.emplace_back(overlap, candidate_end, start, candidate_dist + overlap);
                        
                        // skip to the predecessor's noncolinear shell, whose connections might not be blocked by
                        // this connection
                        for (const pair<size_t, size_t>& shell_pred : noncolinear_shells[candidate_end]) {
                            end_queue.emplace(candidate_dist + shell_pred.second + overlap, shell_pred.first);
                        }
                    }
                    else {
                        // these MEMs are noncolinear, so add this predecessor to the noncolinear shell
                        if (noncolinear_shell.count(candidate_end)) {
                            noncolinear_shell[candidate_end] = std::min(candidate_dist, noncolinear_shell[candidate_end]);
                        }
                        else {
                            noncolinear_shell[candidate_end] = candidate_dist;
                        }
                        
                        // there is no connection to block further connections back, so any of this MEMs
                        // predecessors could still be colinear
                        
                        // find the ends that can reach it directly
                        for (const pair<size_t, size_t>& pred_end : reachable_ends_from_end[candidate_end]) {
                            end_queue.emplace(candidate_dist + pred_end.second, pred_end.first);
                        }
                        
                        // traverse backward through any starts to find more ends that can reach this MEM
                        priority_queue<pair<size_t, size_t>, vector<pair<size_t, size_t>>, std::greater<pair<size_t, size_t>>> pred_start_queue;
                        
                        // initialize the queue with the immediate start neighbors
                        for (const pair<size_t, size_t>& pred_start : reachable_starts_from_end[candidate_end]) {
                            pred_start_queue.emplace(candidate_dist + pred_start.second, pred_start.first);
                        }
                        
                        unordered_set<size_t> pred_traversed;
                        
                        // traverse backwards through starts, stopping at any ends
                        while (!pred_start_queue.empty()) {
                            size_t start_here = pred_start_queue.top().second;
                            size_t start_dist = pred_start_queue.top().first;
                            pred_start_queue.pop();
                            if (pred_traversed.count(start_here)) {
                                continue;
                            }
                            pred_traversed.insert(start_here);
                            
                            for (const pair<size_t, size_t>& pred_end : reachable_ends_from_start[start_here]) {
                                end_queue.emplace(start_dist + pred_end.second, pred_end.first);
                            }
                            for (const pair<size_t, size_t>& start_next : reachable_starts_from_start[start_here]) {
                                pred_start_queue.emplace(start_dist + start_next.second, start_next.first);
                            }
                        }
                    }
                }
                
                // record the node IDs that this MEM's path traverses
                unordered_set<id_t> match_path_nodes;
                Path& match_path = match_nodes[start].path;
                for (size_t j = 0; j < match_path.mapping_size(); j++) {
                    match_path_nodes.insert(match_path.mapping(j).position().node_id());
                }
                
                // queue records are (start/end idx, distance, is an end?)
                int64_t start_length = start_node.end - start_node.begin;
                list<tuple<size_t, int64_t, bool>> path_queue{make_tuple(start, -start_length, true)};
                unordered_set<pair<size_t, bool>> path_enqueued{make_pair(start, true)};
                
                unordered_map<size_t, int64_t> candidate_overlap_preds;
                
                while (!path_queue.empty()) {
                    size_t idx = get<0>(path_queue.front());
                    int64_t dist = get<1>(path_queue.front());
                    bool at_end = get<2>(path_queue.front());
                    path_queue.pop_front();
                    
                    // stop searching backwards at the start of the MEM
                    if (idx == start && !at_end) {
                        break;
                    }
                    
                    if (at_end) {
                        // this end is on the path of the MEM, it might be overlap-colinear
                        if (idx != start) {
                            candidate_overlap_preds[idx] = dist;
                        }
                        
                        // add the start and end predecessors that fall on the path of the MEM to the queue
                        // note: it is okay to not check whether the start/end is behind the start of the path on
                        // the same node because we are guaranteed to run into the start of the path before that
                        for (const pair<size_t, size_t>& start_pred : reachable_starts_from_end[idx]) {
                            if (!path_enqueued.count(make_pair(start_pred.first, false)) &&
                                match_path_nodes.count(start_node_id(start_pred.first))) {
                                
                                path_queue.emplace_back(start_pred.first, dist + (int64_t) start_pred.second, false);
                                path_enqueued.emplace(start_pred.first, false);
                            }
                        }
                        for (const pair<size_t, size_t>& end_pred : reachable_ends_from_end[idx]) {
                            if (!path_enqueued.count(make_pair(end_pred.first, true)) &&
                                match_path_nodes.count(end_node_id(end_pred.first))) {
                                
                                path_queue.emplace_back(end_pred.first, dist + (int64_t) end_pred.second, true);
                                path_enqueued.emplace(end_pred.first, true);
                            }
                        }
                    }
                    else {
                        if (candidate_overlap_preds.count(idx)) {
                            // this start is on the path of the MEM, it cannot be overlap-colinear even if
                            // the end is also on the path
                            ExactMatchNode& match_node = match_nodes[idx];
                            candidate_overlap_preds.erase(idx);
                            // it is also part of the noncolinear shell
                            noncolinear_shell[idx] = dist + (match_node.end - match_node.begin);
                        }
                        
                        // add the start and end predecessors that fall on the path of the MEM to the queue
                        for (const pair<size_t, size_t>& start_pred : reachable_starts_from_start[idx]) {
                            if (!path_enqueued.count(make_pair(start_pred.first, false)) &&
                                match_path_nodes.count(start_node_id(start_pred.first))) {
                                
                                path_queue.emplace_back(start_pred.first, dist + (int64_t) start_pred.second, false);
                                path_enqueued.emplace(start_pred.first, false);
                            }
                        }
                        for (const pair<size_t, size_t>& end_pred : reachable_ends_from_start[idx]) {
                            if (!path_enqueued.count(make_pair(end_pred.first, true)) &&
                                match_path_nodes.count(end_node_id(end_pred.first))) {
                                
                                path_queue.emplace_back(end_pred.first, dist + (int64_t) end_pred.second, true);
                                path_enqueued.emplace(end_pred.first, true);
                            }
                        }
                    }
                }
                
                for (const pair<size_t, int64_t>& overlap_candidate : candidate_overlap_preds) {
                    ExactMatchNode& overlap_node = match_nodes[overlap_candidate.first];
                    int64_t overlap_node_length = overlap_node.end - overlap_node.begin;
                    
                    // how much do the paths overlap (using negative distance)?
                    int64_t overlap = -overlap_candidate.second;
                    
                    // are the matches read colinear after removing the overlap?
                    if (overlap_node.end - overlap <= start_node.begin) {
                        confirmed_overlaps.emplace_back(overlap, overlap_candidate.first, start,
                                                        overlap_candidate.second + overlap);
                    }
                    else if (overlap_node.begin < start_node.begin) {
                        // there is still an even longer read overlap we need to remove
                        overlap = overlap_node.end - start_node.begin;
                        confirmed_overlaps.emplace_back(overlap, overlap_candidate.first, start,
                                                        overlap_candidate.second + overlap);
                    }
                    else {
                        // the overlapping node is still not reachable so it is in the noncolinear shell of this node
                        noncolinear_shell[overlap_candidate.first] = overlap_candidate.second + (overlap_node_length);
                    }
                }
            }
        }
        
        // now we've found all overlap edges, so we can add them into the graph in an order such that they don't
        // conflict (note that all overlap are from an earlier node onto a later one, so we don't need to worry
        // about overlaps coming in from both directions)
        
        // sort in descending order of overlap length
        std::sort(confirmed_overlaps.begin(), confirmed_overlaps.end(),
                  std::greater<tuple<size_t, size_t, size_t, size_t>>());
        
        // split up each node with an overlap edge onto it
        for (auto overlap_record : confirmed_overlaps) {
            ExactMatchNode& from_node = match_nodes[get<1>(overlap_record)];
            ExactMatchNode& onto_node = match_nodes[get<2>(overlap_record)];
            
            size_t suffix_idx = match_nodes.size();
            
            match_nodes.emplace_back();
            ExactMatchNode& suffix_node = match_nodes.back();
            
            // divide up the match on the read
            suffix_node.end = onto_node.end;
            suffix_node.begin = onto_node.begin + get<0>(overlap_record);
            onto_node.end = suffix_node.begin;
            
            // transfer the outgoing edges onto the new node
            suffix_node.edges = std::move(onto_node.edges);
            
            // clear the old edges and add a single edge to the suffix
            onto_node.edges.clear();
            onto_node.edges.emplace_back(suffix_idx, 0);
            
            // add the overlap edge
            from_node.edges.emplace_back(suffix_idx, get<3>(overlap_record));
            
            // store the full path and remove it from the node
            Path full_path = std::move(onto_node.path);
            onto_node.path.Clear();
            
            // add mappings from the path until reaching the overlap point
            int64_t remaining = get<0>(overlap_record);
            int64_t mapping_idx = 0;
            int64_t mapping_len = mapping_from_length(full_path.mapping(mapping_idx));
            while (remaining >= mapping_len) {
                *onto_node.path.add_mapping() = full_path.mapping(mapping_idx);
                mapping_idx++;
                mapping_len = mapping_from_length(full_path.mapping(mapping_idx));
            }
            
            if (remaining) {
                // the overlap point is in the middle of a node, need to split a mapping
                
                const Mapping& split_mapping = full_path.mapping(mapping_idx);
                
                // add the prefix of the mapping to the original node
                Mapping* prefix_split = onto_node.path.add_mapping();
                prefix_split->mutable_position()->set_node_id(split_mapping.position().node_id());
                prefix_split->mutable_position()->set_offset(split_mapping.position().offset());
                Edit* prefix_edit = prefix_split->add_edit();
                prefix_edit->set_from_length(remaining);
                prefix_edit->set_to_length(remaining);
                
                // add the suffix of the mapping to the new node
                Mapping* suffix_split = suffix_node.path.add_mapping();
                suffix_split->mutable_position()->set_node_id(split_mapping.position().node_id());
                suffix_split->mutable_position()->set_offset(remaining);
                Edit* suffix_edit = suffix_split->add_edit();
                suffix_edit->set_from_length(mapping_len - remaining);
                suffix_edit->set_to_length(mapping_len - remaining);
                
                mapping_idx++;
            }
            
            // add the remaining mappings to the suffix node
            for (; mapping_idx < full_path.mapping_size(); mapping_idx) {
                *suffix_node.path.add_mapping() = full_path.mapping(mapping_idx);
            }
        }
    }
    
    // Kahn's algorithm
    void MultipathAlignmentGraph::topological_sort(vector<size_t>& order_out) {
        order_out.resize(match_nodes.size());
        
        vector<size_t> in_degree(match_nodes.size());
        
        for (const ExactMatchNode& match_node : match_nodes) {
            for (const pair<size_t, size_t>& edge : match_node.edges) {
                in_degree[edge.first]++;
            }
        }
        
        list<size_t> source_queue;
        for (size_t i = 0; i < match_nodes.size(); i++) {
            if (in_degree[i] == 0) {
                source_queue.push_back(i);
            }
        }
        
        size_t next = 0;
        while (!source_queue.empty()) {
            size_t src = source_queue.front();
            source_queue.pop_front();
            
            for (const pair<size_t, size_t>& edge : match_nodes[src].edges) {
                in_degree[edge.first]--;
                if (in_degree[edge.first] == 0) {
                    source_queue.push_back(edge.first);
                }
            }
            
            order_out[next] = src;
            next++;
        }
    }
    
    void MultipathAlignmentGraph::remove_transitive_edges(const vector<size_t>& topological_order) {
        vector<size_t> topological_position(topological_order.size());
        for (size_t i = 0; i < topological_order.size(); i++) {
            topological_position[topological_order[i]] = i;
        }
        
        for (size_t i : topological_order) {
            vector<pair<size_t, size_t>>& edges = match_nodes[i].edges;
            std::sort(edges.begin(), edges.end(), [&](const pair<size_t, size_t>& edge_1, const pair<size_t, size_t>& edge_2) {
                return topological_position[edge_1.first] < topological_position[edge_2.first];
            });
            
            unordered_map<size_t, bool> keep_edge;
            for (size_t j = 0; j < edges.size(); j++) {
                keep_edge[edges[j].first] = true;
            }
            
            for (const pair<size_t, size_t>& edge : edges) {
                if (!keep_edge[edge.first]) {
                    continue;
                }
                
                vector<size_t> stack{edge.first};
                unordered_set<size_t> traversed{edge.first};
                while (!stack.empty()) {
                    size_t here = stack.back();
                    stack.pop_back();
                    
                    for (const pair<size_t, size_t>& next_edge : match_nodes[next_edge.first].edges) {
                        if (keep_edge.count(next_edge.first)) {
                            keep_edge[next_edge.first] = false;
                        }
                        
                        if (!traversed.count(next_edge.first)) {
                            traversed.insert(next_edge.first);
                            stack.push_back(next_edge.first);
                        }
                    }
                }
            }
            
            size_t j = 0;
            size_t end = edges.size();
            while (j < end) {
                if (keep_edge[edges[j].first]) {
                    j++;
                }
                else {
                    end--;
                    edges[j] = edges[end];
                }
            }
            edges.resize(end);
        }
        
    }
    
    void MultipathAlignmentGraph::prune_to_high_scoring_paths(const BaseAligner& aligner, int32_t max_suboptimal_score_diff,
                                                              const vector<size_t>& topological_order) {
        
        unordered_map<pair<size_t, size_t>, int32_t> edge_weights;

        vector<int32_t> node_weights(match_nodes.size());
    
        // compute the weight of edges and node matches
        for (size_t i = 0; i < match_nodes.size(); i++) {
            ExactMatchNode& from_node = match_nodes[i];
            node_weights[i] = aligner.match * (from_node.end - from_node.begin);
            
            for (const pair<size_t, size_t>& edge : from_node.edges) {
                ExactMatchNode& to_node = match_nodes[edge.first];
                
                int64_t graph_dist = edge.second;
                int64_t read_dist = to_node.begin - from_node.end;
                
                if (read_dist > graph_dist) {
                    // the read length in between the MEMs is longer than the distance, suggesting a read insert
                    int64_t gap_length = read_dist - graph_dist;
                    edge_weights[make_pair(i, edge.first)] = -aligner.mismatch * graph_dist - (gap_length - 1) * aligner.gap_extension
                                                             - aligner.gap_open;
                }
                else if (read_dist < graph_dist) {
                    // the read length in between the MEMs is shorter than the distance, suggesting a read deletion
                    int64_t gap_length = graph_dist - read_dist;
                    edge_weights[make_pair(i, edge.first)] = -aligner.mismatch * read_dist - (gap_length - 1) * aligner.gap_extension
                                                             - aligner.gap_open;
                }
                else {
                    // the read length in between the MEMs is the same as the distance, suggesting a pure mismatch
                    edge_weights[make_pair(i, edge.first)] = -aligner.mismatch * read_dist;
                }
            }
        }
        
        vector<int32_t> forward_scores = node_weights;
        vector<int32_t> backward_scores = node_weights;
        
        // forward DP
        for (int64_t i = 0; i < topological_order.size(); i++) {
            size_t idx = topological_order[i];
            ExactMatchNode& from_node = match_nodes[idx];
            for (const pair<size_t, size_t>& edge : from_node.edges) {
                forward_scores[edge.first] = std::max(forward_scores[edge.first],
                                                      forward_scores[idx] + edge_weights[make_pair(idx, edge.first)]);
            }
        }
        
        // backward DP
        for (int64_t i = topological_order.size() - 1; i >= 0; i--) {
            size_t idx = topological_order[i];
            ExactMatchNode& from_node = match_nodes[idx];
            for (const pair<size_t, size_t>& edge : from_node.edges) {
                backward_scores[idx] = std::max(backward_scores[idx],
                                                backward_scores[edge.first] + edge_weights[make_pair(idx, edge.first)]);
            }
        }
        
        // compute the minimum score we will require of a node or edge
        int32_t min_path_score = *std::max_element(forward_scores.begin(), forward_scores.end()) - max_suboptimal_score_diff;
        
        // use forward-backward to find nodes/edges on some path with a score above the minimum
        unordered_set<size_t> keep_nodes;
        unordered_set<pair<size_t, size_t>> keep_edges;
        for (size_t i = 0; i < match_nodes.size(); i++) {
            if (forward_scores[i] + backward_scores[i] - node_weights[i]) {
                keep_nodes.insert(i);
                for (const pair<size_t, size_t>& edge : match_nodes[i].edges) {
                    if (forward_scores[i] + backward_scores[edge.first] + edge_weights[make_pair(i, edge.first)] >= min_path_score) {
                        keep_edges.emplace(i, edge.first);
                    }
                }
            }
        }
        
        // prune down to these nodes and edges
        size_t next = 0;
        for (size_t i = 0; i < match_nodes.size(); i++) {
            if (keep_nodes.count(i)) {
                if (i != next) {
                    match_nodes[next] = std::move(match_nodes[i]);
                }
                ExactMatchNode& match_node = match_nodes[next];
                auto new_end = std::remove_if(match_node.edges.begin(), match_node.edges.end(),
                                              [&](const pair<size_t, size_t>& edge) { return keep_edges.count(make_pair(i, edge.first)); });
                match_node.edges.resize(new_end - match_node.edges.begin());
                next++;
            }
        }
        match_nodes.resize(next);
    }
    
    MultipathClusterer::MultipathClusterer(const Alignment& alignment,
                                           const vector<MaximalExactMatch>& mems,
                                           const QualAdjAligner& aligner,
                                           xg::XG& xgindex,
                                           int8_t full_length_bonus,
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
            
            
            int32_t mem_score = aligner.score_exact_match(mem.begin, mem.end,
                                                          alignment.quality().begin() + (mem.begin - alignment.sequence().begin()));
            
            for (gcsa::node_type mem_hit : mem.nodes) {
                nodes.emplace_back(mem, make_pos_t(mem_hit), mem_score);
                maximum_detectable_gaps.push_back(max_gaps);
            }
        }
        
        // now we use the distance approximation to cluster the MEM hits according to the strand
        // they fall on using the oriented distance estimation function
        
        // for recording the distance of any pair that we check with a finite distance
        unordered_map<pair<size_t, size_t>, int64_t> recorded_finite_dists;
        // for recording the number of times elements of a strand cluster have been compared
        // and found an infinite distance
        map<pair<size_t, size_t>, size_t> num_infinite_dists;
        
        // deterministically generate pseudo-shuffled pairs in constant time per pair
        // method taken from https://stackoverflow.com/questions/1866684/algorithm-to-print-out-a-shuffled-list-in-place-and-with-o1-memory
        size_t num_pairs = nodes.size() * nodes.size();
        assert(num_pairs <= 4294967296) // 2^32, else range_max can't get large enough
        size_t permutation_idx = 0;
        size_t range_max = 1;
        while (range_max < num_pairs) {
            range_max *= 2;
        }
        auto get_next_pair = [&]() {
            pair<size_t, size_t> to_return;
            do {
                size_t permuted = ((permutation_idx ^ (range_max - 1)) ^ (permutation_idx << 2) + 3) & (range_max - 1);
                to_return.first = current_pair_idx / nodes.size();
                to_return.second = current_pair_idx % nodes.size();
            } while (to_return.first >= to_return.second || permuted >= num_pairs);
            return to_return;
        }
        
        UnionFind union_find(nodes.size());
        
        size_t num_possible_merges_remaining = (nodes.size() * (nodes.size() - 1)) / 2;
        size_t pairs_checked = 0;
        
        // a simulated annealing parameter loosely inspired by the cutoff for an Erdos Renyi random graph
        // to be connected with probability approaching 1
        size_t current_max_num_probes = (size_t) (2.0 * ceil(log(nodes.size())));
        
        while (num_possible_merges_remaining > 0 && pairs_checked < num_pairs) {
            
            pair<size_t, size_t> node_pair = get_next_pair();
            
            pairs_checked++;
            
            size_t strand_1 = union_find.find(node_pair.first);
            size_t strand_2 = union_find.find(node_pair.second);

            if (strand_1 == strand_2) {
                // these are already identified as on the same strand, don't need to do it again
                continue;
            }
            
            if (num_infinite_dists[make_pair(strand_1, strand_2)] >= current_max_num_probes) {
                // we've already checked multiple distances between these strand clusters and
                // none have returned a finite distance, so we conclude that they are in fact
                // on separate clusters and decline to check any more distances
                continue;
            }
            
            pos_t& pos_1 = nodes[node_pair.first].start_pos;
            pos_t& pos_2 = nodes[node_pair.second].start_pos;
            
            int64_t oriented_dist = xgindex.closest_shared_path_oriented_distance(id(pos_1), offset(pos_1), is_rev(pos_1),
                                                                                  id(pos_2), offset(pos_2), is_rev(pos_2));
            
            if (oriented_dist == std::numeric_limits<int64_t>::max()) {
                // distance is estimated at infinity, so these are either on different strands
                // or the path heuristic failed to find a shared path
                
                num_infinite_dists[make_pair(strand_1, strand_2)]++;
                num_infinite_dists[make_pair(strand_2, strand_1)]++;
                
                // this infinite distance pushed the count over the maximum number of probes, so remove
                // these merges from the pool of potential merges remaining
                if (num_infinite_dists[make_pair(strand_1, strand_2)] >= current_max_num_probes) {
                    size_t strand_size_1 = union_find.group_size(strand_1);
                    size_t strand_size_2 = union_find.group_size(strand_2);
                    
                    num_possible_merges_remaining -= strand_size_1 * strand_size_2;
                }
            }
            else {
                // the distance is finite, so merge the strand clusters
                
                recorded_finite_dists[node_pair] = oriented_dist;
                
                size_t strand_size_1 = union_find.group_size(strand_1);
                size_t strand_size_2 = union_find.group_size(strand_2);
                
                union_find.union(node_pair.first, node_pair.second);
                
                // remove these from the pool of remaining merges
                num_possible_merges_remaining -= strand_size_1 * strand_size_2;
                
                size_t strand_retaining = union_find.find(node_pair.first);
                size_t strand_removing = strand_retaining == strand_1 ? strand_2 : strand_1;
                
                // collect the number of times the strand cluster thas is being removed has had an infinite distance
                // to other strands besides the one it's being merged into
                vector<pair<size_t, size_t>> inf_counts;
                auto end = num_infinite_dists.upper_bound(make_pair(strand_removing, numeric_limits<size_t>::max()));
                auto begin = num_infinite_dists.lower_bound(make_pair(strand_removing, 0));
                for (auto iter = begin; iter != end; iter++) {
                    if (*iter).first.second != strand_retaining) {
                        inf_counts.emplace_back((*iter).first.second, (*iter).second);
                    }
                }
                
                for (const pair<size_t, size_t> inf_count : inf_counts) {
                    size_t& curr_num_probes = num_infinite_dists[make_pair(strand_retaining, inf_count.first)];
                    bool already_blocked = curr_num_probes >= current_max_num_probes;
                    
                    // transfer these counts over to the strand cluster that is being retained
                    curr_num_probes += inf_count.second;
                    num_infinite_dists[make_pair(inf_count.first, strand_retaining)] += inf_count.second;
                    
                    // remove the strand from the infinite distance counter
                    num_infinite_dists.erase(make_pair(strand_removing, inf_count.first));
                    num_infinite_dists.erase(make_pair(inf_count.first, strand_removing));
                    
                    // if adding these counts pushed the cluster over the strand probe max, remove these merges from
                    // the pool remaining
                    if (curr_num_probes >= current_max_num_probes && !already_blocked) {
                        num_possible_merges_remaining -= (strand_size_1 + strand_size_2) * union_find.group_size(inf_count.first);
                    }
                }
            }
            
            // slowly lower the number of distances we need to check before we believe that two clusters are on
            // separate strands
            if (pairs_checked % nodes.size() == nodes.size() - 1) {
                current_max_num_probes--;
                
                for (const pair<pair<size_t, size_t>, size_t>& inf_dist_record : num_infinite_dists) {
                    // break symmetry so we don't repeat the operation twice
                    if (inf_dist_record.first.first < inf_dist_record.first.second && inf_dist_record.second == current_max_num_probes) {
                        // this merge just fell below the maximum number of distance probes
                        size_t strand_size_1 = union_find.group_size(inf_dist_record.first.first);
                        size_t strand_size_2 = union_find.group_size(inf_dist_record.first.second);
                        num_possible_merges_remaining -= strand_size_1 * strand_size_2;
                    }
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
                                   -recorded_finite_dists[make_pair(next, curr)];
                    
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
                    int64_t graph_dist = max(0ll, sorted_pos[j].first - strand_pos - pivot_length);
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
                        
                        int64_t extra_dist = max(0ll, gap_length);
                        
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
    
    void MultipathClusterer::identify_sources_and_sinks(vector<size_t>& sources_out,
                                                        vector<size_t>& sinks_out) {
        
        sources_out.clear();
        sinks_out.clear();
        
        vector<bool> is_source(nodes.size(), true);
        
        for (size_t i = 0; i < nodes.size(); i++) {
            if (nodes[i].edges_from.empty()) {
                sinks_out.push_back(i);
            }
            
            for (MPCEdge& edge : nodes[i].edges_from) {
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
                
                MPCNode& node = nodes[stack.back()];
                stack.pop_back();
                
                // search in both forward and backward directions
                
                for (MPCEdge& edge : node.edges_from) {
                    
                    if (!enqueued[edge.to_idx]) {
                        stack.push_back(edge.to_idx);
                        enqueued[edge.to_idx] = true;
                        components_out.back().push_back(edge.to_idx);
                    }
                }
                
                for (MPCEdge& edge : node.edges_to) {
                    
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
    
    vector<vector<pair<const MaximalExactMatch*, pos_t>>> MultipathClusterer::clusters(int32_t max_qual_score) {
        
        vector<vector<pair<const MaximalExactMatch*, pos_t>>> to_return;
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
                int32_t dp_score = nodes[component[i]].dp_score;
                if (dp_score > traceback_end.first) {
                    traceback_end.first = dp_score;
                    traceback_end.second = j;
                }
            }
        }
        
        // sort indices in descending order by their highest traceback score
        vector<size_t> order = range_vector(0, components.size());
        std::sort(order.begin(), order.end(), [&](const size_t i, const size_t j) {
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
            if (4.3429448190325183 * aligner.log_base * (top_score - component_traceback_ends[i].first) > max_qual_score ) {
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



