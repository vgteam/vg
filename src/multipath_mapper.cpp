//
//  multipath_mapper.cpp
//  
//
//

#define debug_multipath_mapper

#include "multipath_mapper.hpp"

namespace vg {
    
    MultipathMapper::MultipathMapper(xg::XG* xg_index, gcsa::GCSA* gcsa_index, gcsa::LCPArray* lcp_array,
                                     SnarlManager* snarl_manager) :
        BaseMapper(xg_index, gcsa_index, lcp_array),
        snarl_manager(snarl_manager)
    {
        // nothing to do
    }
    
    MultipathMapper::~MultipathMapper() {
        
    }
    
    void MultipathMapper::multipath_map(const Alignment& alignment,
                                        list<MultipathAlignment>& multipath_alns_out,
                                        size_t max_alt_alns) {
    
#ifdef debug_multipath_mapper
        cerr << "multipath mapping read " << pb2json(alignment) << endl;
        cerr << "querying MEMs..." << endl;
#endif
    
        // query MEMs using GCSA2
        double dummy;
        vector<MaximalExactMatch> mems = find_mems_deep(alignment.sequence().begin(), alignment.sequence().end(),
                                                        dummy, 0, min_mem_length, mem_reseed_length);
        
#ifdef debug_multipath_mapper
        cerr << "obtained MEMs:" << endl;
        for (MaximalExactMatch mem : mems) {
            cerr << "\t" << mem << endl;
        }
        cerr << "clustering MEMs..." << endl;
#endif
        
        // cluster the MEMs
        OrientedDistanceClusterer clusterer(alignment, mems, *qual_adj_aligner, xindex, max_expected_dist_approx_error);
        vector<vector<pair<const MaximalExactMatch*, pos_t>>> clusters = clusterer.clusters();
        
#ifdef debug_multipath_mapper
        cerr << "obtained clusters:" << endl;
        for (int i = 0; i < clusters.size(); i++) {
            cerr << "\tcluster " << i + 1 << endl;
            for (pair<const MaximalExactMatch*, pos_t>  hit : clusters[i]) {
                cerr << "\t\t" << hit.second << " " <<  hit.first->sequence() << endl;
            }
        }
        cerr << "extracting subgraphs..." << endl;
#endif
        
        // extract graphs around the clusters
        vector<VG*> cluster_graphs;
        unordered_map<VG*, vector<pair<const MaximalExactMatch*, pos_t>>> cluster_graph_mems;
        query_cluster_graphs(alignment, mems, clusters, cluster_graphs, cluster_graph_mems);
        
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
        
#ifdef debug_multipath_mapper
        cerr << "aligning to subgraphs..." << endl;
#endif
        
        // align to each cluster subgraph
        size_t num_alns = 0;
        for (VG* vg : cluster_graphs) {
            // if we have a cluster graph with small enough MEM coverage compared to the best one or we've made
            // the maximum number of alignments we stop producing alternate alignments
            if (cluster_graph_coverage[vg] < mem_coverage_min_ratio * cluster_graph_coverage[cluster_graphs[0]]
                || num_alns >= max_alt_alns) {
                break;
            }
            
#ifdef debug_multipath_mapper
            cerr << "performing alignment to subgraph " << pb2json(vg->graph) << endl;
#endif
            
            vector<pair<const MaximalExactMatch*, pos_t>>& graph_mems = cluster_graph_mems[vg];
            multipath_align(alignment, vg, graph_mems, multipath_alns_out);
            
            num_alt_alns++;
        }
        
        for (VG* vg : cluster_graphs) {
            delete vg;
        }
    }
    
    void MultipathMapper::query_cluster_graphs(const Alignment& alignment,
                                               const vector<MaximalExactMatch>& mems,
                                               const vector<vector<pair<const MaximalExactMatch*, pos_t>>>& clusters,
                                               vector<VG*>& cluster_graphs_out,
                                               unordered_map<VG*, vector<pair<const MaximalExactMatch*, pos_t>>>& cluster_graph_mems_out){
        
        // subgraphs around each cluster
        cluster_graphs_out.reserve(clusters.size());
        
        // we will ensure that nodes are in only one cluster, use this to record which one
        unordered_map<id_t, size_t> node_id_to_cluster;
        
        for (size_t i = 0; i < clusters.size(); i++) {
            
            // gather the parameters for subgraph extraction from the MEM hits
            
            const vector<pair<const MaximalExactMatch*, pos_t>>& cluster = clusters[i];
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
                forward_max_dist.push_back(qual_adj_aligner->longest_detectable_gap(alignment, mem_hit.first->end)
                                           + (mem_hit.first->end - mem_hit.first->begin));
                backward_max_dist.push_back(qual_adj_aligner->longest_detectable_gap(alignment, mem_hit.first->begin));
            }
            
            // TODO: a progressive expansion of the subgraph if the MEM hit is already contained in
            // a cluster graph somewhere?
            
            // extract the subgraph within the search distance
            
            VG* cluster_graph = new VG();
            Graph& graph = cluster_graph->graph;
            
            // extract the protobuf Graph in place in the VG
            algorithms::extract_containing_graph(*xindex, graph, positions, forward_max_dist,
                                                 backward_max_dist, &get_node_cache());
            
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
                
                cluster_graphs_out.push_back(cluster_graph);
                
                // now that we know we're going to save the graph, manually trigger the index building since
                // we circumvented the constructors
                cluster_graph->rebuild_indexes();
            }
            else {
                // this graph overlap at least one other graph, so we merge them into the (arbitrarily chosen) graph
                // with the minimum index in the vector
                size_t min_idx_cluster = *std::min_element(overlapping_graphs.begin(), overlapping_graphs.end());
                
                // merge in the new graph
                cluster_graphs_out[min_idx_cluster]->extend(graph);
                delete cluster_graph;
                
                // if this subgraph chains together multiple clusters, merge them and remove them from the list
                overlapping_graphs.erase(min_idx_cluster);
                for (size_t j : overlapping_graphs) {
                    std::swap(cluster_graphs_out[j], cluster_graphs_out.back());
                    cluster_graphs_out[min_idx_cluster]->extend(cluster_graphs_out.back()->graph);
                    delete cluster_graphs_out.back();
                    cluster_graphs_out.pop_back();
                }
                
                // relabel the cluster of any nodes that were not in the graph we merged into
                Graph& merged_graph = cluster_graphs_out[min_idx_cluster]->graph;
                for (size_t j = 0; j < merged_graph.node_size(); j++) {
                    node_id_to_cluster[merged_graph.node(j).id()] = min_idx_cluster;
                }
            }
        }
        
        // which MEMs are in play for which cluster?
        for (const MaximalExactMatch& mem : mems) {
            for (gcsa::node_type hit : mem.nodes) {
                id_t node_id = gcsa::Node::id(hit);
                if (node_id_to_cluster.count(node_id)) {
                    size_t cluster_idx = node_id_to_cluster[node_id];
                    cluster_graph_mems_out[cluster_graphs_out[cluster_idx]].push_back(make_pair(&mem, make_pos_t(hit)));
                }
            }
        }
    }
    
    void MultipathMapper::multipath_align(const Alignment& alignment, VG* vg,
                                          vector<pair<const MaximalExactMatch*, pos_t>>& graph_mems,
                                          list<MultipathAlignment>& multipath_alns_out) {
#ifdef debug_multipath_mapper
        cerr << "constructing alignment graph" << endl;
#endif
        
        // the longest path we could possibly align to (full gap and a full sequence)
        size_t target_length = qual_adj_aligner->longest_detectable_gap(alignment) + alignment.sequence().size();
        
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
            cerr << trans_record.second.first << "->" << trans_record.first << (trans_record.second.second ? "-" : "+") << endl;
            rev_trans.insert(make_pair(trans_record.second.first,
                                       make_pair(trans_record.first, trans_record.second.second)));
        }
#ifdef debug_multipath_mapper
        cerr << "making multipath alignment MEM graph" << endl;
#endif
        // construct a graph that summarizes reachability between MEMs
        MultipathAlignmentGraph multi_aln_graph(align_graph, graph_mems, rev_trans, node_trans, snarl_manager, max_snarl_cut_size);
        
        vector<size_t> topological_order;
        multi_aln_graph.topological_sort(topological_order);
        
        cerr << "computed topological sort" << endl;
        
        // it's sometimes possible for transitive edges to survive the original construction algorithm, so remove them
        multi_aln_graph.remove_transitive_edges(topological_order);
        
        cerr << "removed transitive edges" << endl;
        
        // prune this graph down the paths that have reasonably high likelihood
        multi_aln_graph.prune_to_high_scoring_paths(*qual_adj_aligner, max_suboptimal_path_score_diff, topological_order);
        
        cerr << "pruned to high scoring paths" << endl;
        
        // create a new multipath alignment object and transfer over data from alignment
        multipath_alns_out.emplace_back();
        MultipathAlignment& multipath_aln = multipath_alns_out.back();
        multipath_aln.set_sequence(alignment.sequence());
        multipath_aln.set_quality(alignment.quality());
        multipath_aln.set_name(alignment.name());
        multipath_aln.set_sample_name(alignment.sample_name());
        multipath_aln.set_read_group(alignment.read_group());
        
        cerr << "transferred over read information" << endl;
        
        // add a subpath for each of the exact match nodes
        for (int64_t j = 0; j < multi_aln_graph.match_nodes.size(); j++) {
            ExactMatchNode& match_node = multi_aln_graph.match_nodes[j];
            Subpath* subpath = multipath_aln.add_subpath();
            *subpath->mutable_path() = match_node.path;
            int32_t match_score = qual_adj_aligner->score_exact_match(match_node.begin, match_node.end,
                                                                     alignment.quality().begin() + (match_node.begin - alignment.sequence().begin()));
            subpath->set_score(match_score + qual_adj_aligner->full_length_bonus * ((match_node.begin == alignment.sequence().begin()) +
                                                                                   (match_node.end == alignment.sequence().end())));
        }
        
#ifdef debug_multipath_mapper
        cerr << "doing DP between MEMs" << endl;
#endif
        
        // perform alignment in the intervening sections
        for (int64_t j = 0; j < multi_aln_graph.match_nodes.size(); j++) {
#ifdef debug_multipath_mapper
            cerr << "checking for intervening alignments from match node " << j << endl;
#endif
            
            ExactMatchNode& src_match_node = multi_aln_graph.match_nodes[j];
            Subpath* src_subpath = multipath_aln.mutable_subpath(j);
            
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
#ifdef debug_multipath_mapper
                    cerr << "match node has a simple split connection onto node " << edge.first << endl;
#endif
                    // connect the supath directly without an intervening alignment
                    src_subpath->add_next(edge.first);
                    continue;
                }
            }
            
            const Path& path = src_subpath->path();
            const Mapping& final_mapping = path.mapping(path.mapping_size() - 1);
            const Position& final_mapping_position = final_mapping.position();
            // make a pos_t that points to the final base in the match
            pos_t src_pos = make_pos_t(final_mapping_position.node_id(),
                                       final_mapping_position.is_reverse(),
                                       final_mapping_position.offset() + mapping_from_length(final_mapping) - 1);
            
            // the longest gap that could be detected at this position in the read
            size_t src_max_gap = qual_adj_aligner->longest_detectable_gap(alignment, src_match_node.end);
            
            for (const pair<size_t, size_t>& edge : src_match_node.edges) {
                ExactMatchNode& dest_match_node = multi_aln_graph.match_nodes[edge.first];
                pos_t dest_pos = make_pos_t(multipath_aln.subpath(edge.first).path().mapping(0).position());
                
#ifdef debug_multipath_mapper
                cerr << "forming intervening alignment for edge to node " << edge.first << endl;
#endif
                
                size_t intervening_length = dest_match_node.begin - src_match_node.end;
                size_t max_dist = intervening_length + std::min(src_max_gap, qual_adj_aligner->longest_detectable_gap(alignment, dest_match_node.begin)) + 1;
                
                // extract the graph between the matches
                Graph connecting_graph;
                unordered_map<id_t, id_t> connect_trans = algorithms::extract_connecting_graph(align_graph,      // DAG with split strands
                                                                                               connecting_graph, // graph to extract into
                                                                                               max_dist,         // longest distance necessary
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
                
#ifdef debug_multipath_mapper
                cerr << "aligning to connecting graph: " << pb2json(connecting_graph) << endl;
#endif
                
                // TODO a better way of choosing the number of alternate alignments
                // TODO alternate alignments restricted only to distinct node paths?
                vector<Alignment> alt_alignments;
                qual_adj_aligner->align_global_banded_multi(intervening_sequence, alt_alignments,
                                                           connecting_graph, num_alt_alns, band_padding, true);
                                
                for (Alignment& connecting_alignment : alt_alignments) {
                    // create a subpath between the matches for this alignment
                    Subpath* connecting_subpath = multipath_aln.add_subpath();
                    connecting_subpath->set_score(connecting_alignment.score());
                    Path* subpath_path = connecting_subpath->mutable_path();
                    const Path& aligned_path = connecting_alignment.path();
                    const Mapping& first_mapping = aligned_path.mapping(0);
                    const Mapping& last_mapping = aligned_path.mapping(aligned_path.mapping_size() - 1);
                    int32_t rank = 1;
                    // check to make sure the first is not an empty anchoring mapping
                    if (mapping_from_length(first_mapping) != 0 || mapping_to_length(first_mapping) != 0) {
                        *subpath_path->add_mapping() = first_mapping;
                        rank++;
                    }
                    // add all mapping in between the ends
                    for (size_t j = 1; j < aligned_path.mapping_size() - 1; j++) {
                        Mapping* mapping = subpath_path->add_mapping();
                        *mapping = aligned_path.mapping(j);
                        mapping->set_rank(rank);
                        rank++;
                    }
                    // check to make sure the last is not an empty anchoring mapping or the same as the first
                    if ((mapping_from_length(last_mapping) != 0 || mapping_to_length(last_mapping) != 0)
                        && aligned_path.mapping_size() > 1) {
                        Mapping* mapping = subpath_path->add_mapping();
                        *mapping = last_mapping;
                        mapping->set_rank(rank);
                    }
                    
                    // add the appropriate connections
                    src_subpath->add_next(multipath_aln.subpath_size() - 1);
                    connecting_subpath->add_next(edge.first);
                    
                    // translate the path into the space of the main graph
                    translate_node_ids(*connecting_subpath->mutable_path(), connect_trans);
                    Mapping* firste_subpath_mapping = connecting_subpath->mutable_path()->mutable_mapping(0);
                    if (firste_subpath_mapping->position().node_id() == final_mapping.position().node_id()) {
                        firste_subpath_mapping->mutable_position()->set_offset(mapping_from_length(final_mapping));
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
                    
                    int64_t tail_length = alignment.sequence().end() - match_node.end;
                    int64_t target_length = tail_length + qual_adj_aligner->longest_detectable_gap(alignment, match_node.end) + (alignment.sequence().end() - match_node.end);
                    pos_t end_pos = final_position(match_node.path);
                    
                    Graph tail_graph;
                    unordered_map<id_t, id_t> tail_trans = algorithms::extract_extending_graph(align_graph,
                                                                                               tail_graph,
                                                                                               target_length,
                                                                                               end_pos,
                                                                                               false,         // search forward
                                                                                               false);        // no need to preserve cycles (in a DAG)
                    
                    // flip doubly reversing edges b/c gssw doesn't like them
                    for (size_t k = 0; k < tail_graph.edge_size(); k++) {
                        Edge* edge = tail_graph.mutable_edge(k);
                        if (edge->from_start() && edge->to_end()) {
                            id_t tmp = edge->from();
                            edge->set_from(edge->to());
                            edge->set_to(tmp);
                            edge->set_from_start(false);
                            edge->set_to_end(false);
                        }
                    }
                    
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
                    qual_adj_aligner->align_pinned_multi(right_tail_sequence, alt_alignments, tail_graph, true, num_alt_alns);
                    
#ifdef debug_multipath_mapper
                    cerr << "aligning sequence: " << right_tail_sequence.sequence() << endl << "to right tail graph: " << pb2json(tail_graph) << endl;
#endif
                    
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
                    
                    int64_t tail_length = match_node.begin - alignment.sequence().begin();
                    int64_t target_length = tail_length + qual_adj_aligner->longest_detectable_gap(alignment, match_node.begin) + (match_node.begin - alignment.sequence().begin());
                    pos_t begin_pos = initial_position(match_node.path);
                    
                    Graph tail_graph;
                    unordered_map<id_t, id_t> tail_trans = algorithms::extract_extending_graph(align_graph,
                                                                                               tail_graph,
                                                                                               target_length,
                                                                                               begin_pos,
                                                                                               true,          // search backward
                                                                                               false);        // no need to preserve cycles (in a DAG)
                    
                    // flip doubly reversing edges b/c gssw doesn't like them
                    for (size_t k = 0; k < tail_graph.edge_size(); k++) {
                        Edge* edge = tail_graph.mutable_edge(k);
                        if (edge->from_start() && edge->to_end()) {
                            id_t tmp = edge->from();
                            edge->set_from(edge->to());
                            edge->set_to(tmp);
                            edge->set_from_start(false);
                            edge->set_to_end(false);
                        }
                    }
                    
                    Alignment left_tail_sequence;
                    left_tail_sequence.set_sequence(alignment.sequence().substr(0, match_node.begin - alignment.sequence().begin()));
                    if (!alignment.quality().empty()) {
                        left_tail_sequence.set_quality(alignment.quality().substr(0, match_node.begin - alignment.sequence().begin()));
                    }
                    
#ifdef debug_multipath_mapper
                    cerr << "aligning sequence: " << left_tail_sequence.sequence() << endl << "to left tail graph: " << pb2json(tail_graph) << endl;
#endif
                    
                    vector<Alignment> alt_alignments;
                    qual_adj_aligner->align_pinned_multi(left_tail_sequence, alt_alignments, tail_graph, false, num_alt_alns);
                    
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
        
#ifdef debug_multipath_mapper
        cerr << "multipath alignment before translation: " << pb2json(multipath_aln) << endl;
#endif
        for (size_t j = 0; j < multipath_aln.subpath_size(); j++) {
            translate_oriented_node_ids(*multipath_aln.mutable_subpath(j)->mutable_path(), node_trans);
        }
#ifdef debug_multipath_mapper
        cerr << "completed multipath alignment: " << pb2json(multipath_aln) << endl;
#endif
        
    }
    
    int64_t MultipathMapper::read_coverage(const vector<pair<const MaximalExactMatch*, pos_t>>& mem_hits) {
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
    
    double MultipathMapper::read_coverage_z_score(int64_t coverage, const Alignment& alignment) {
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
#ifdef debug_multipath_mapper
        cerr << "walking out MEMs in graph" << endl;
#endif
        
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
            
#ifdef debug_multipath_mapper
            cerr << "walking MEM hit " << hit_pos << " " << hit.first->sequence() << endl;
#endif
            
            auto hit_range = injection_trans.equal_range(id(hit_pos));
            for (auto iter = hit_range.first; iter != hit_range.second; iter++) {
                // this graph is unrolled/dagified, so all orientations should match
                if ((*iter).second.second != is_rev(hit_pos)) {
                    continue;
                }
                
                // an id that corresponds to the original node
                id_t injected_id = (*iter).second.first;
                
#ifdef debug_multipath_mapper
                cerr << "hit node exists in graph as " << injected_id << endl;
#endif
                
                // check all MEMs that traversed this node to see if this is a redundant sub-MEM
                bool is_partial_mem = false;
                for (int64_t j : node_matches[id(hit_pos)]) {
                    ExactMatchNode& match_node = match_nodes[j];
                    
#ifdef debug_multipath_mapper
                    cerr << "there is a previous node that visited this hit, checking whether it is a parent MEM" << endl;
#endif
                    
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
                
#ifdef debug_multipath_mapper
                cerr << "performing DFS to walk out match" << endl;
#endif
                
                // stack for DFS, each record contains tuples of (read begin, node offset, next node index, next node ids)
                vector<tuple<string::const_iterator, size_t, size_t, vector<NodeTraversal>>> stack;
                stack.emplace_back(begin, offset(hit_pos), 0,
                                   vector<NodeTraversal>{NodeTraversal(vg.get_node(injected_id))});
                
                while (!stack.empty()) {
                    auto& back = stack.back();
                    if (get<2>(back) == get<3>(back).size()) {
#ifdef debug_multipath_mapper
                        cerr << "traversed all edges out of current traversal" << endl;
#endif
                        stack.pop_back();
                        continue;
                    }
                    NodeTraversal trav = get<3>(back)[get<2>(back)];
                    get<2>(back)++;
                    
#ifdef debug_multipath_mapper
                    cerr << "checking node " << trav.node->id() << endl;
#endif
                    
                    const string& node_seq = trav.node->sequence();
                    size_t node_idx = get<1>(back);
                    string::const_iterator read_iter = get<0>(back);
                    
                    // look for a match along the entire node sequence
                    for (; node_idx < node_seq.size() && read_iter != end; node_idx++, read_iter++) {
                        if (node_seq[node_idx] != *read_iter) {
#ifdef debug_multipath_mapper
                            cerr << "node sequence does not match read" << endl;
#endif
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
                }
                
                // we left the trace in the stack, which means we found a complete match
                if (stack.empty()) {
                    cerr << "error:[MultipathMapper] Couldn't find a match for MEM " << *hit.first << endl;
                    exit(1);
                }
                
#ifdef debug_multipath_mapper
                cerr << "converting into a Path" << endl;
#endif
                
                int64_t match_node_idx = match_nodes.size();
                match_nodes.emplace_back();
                ExactMatchNode* match_node = &match_nodes.back();
                Path* path = &match_node->path;
                match_node->begin = begin;
                match_node->end = end;
                int64_t length_remaining = end - begin;
                
                // walk out the match
                int32_t rank = 1;
                for (auto search_record : stack) {
                    int64_t offset = get<1>(search_record);
                    Node* node = get<3>(search_record)[get<2>(search_record) - 1].node;
                    int64_t length = std::min((int64_t) node->sequence().size() - offset, length_remaining);
                    
                    Mapping* mapping = path->add_mapping();
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
                
#ifdef debug_multipath_mapper
                cerr << pb2json(*path) << endl;
#endif
                
                if (cutting_snarls) {
#ifdef debug_multipath_mapper
                    cerr << "cutting with snarls" << endl;
#endif
                    // we indicated a snarl manager that owns the snarls we want to cut out of exact matches
                    
                    // first compute the segments we want to cut out
                    
                    // this list holds the beginning of the current segment at each depth in the snarl hierarchy
                    // as we traverse the exact match, the beginning is recorded in both sequence distance and node index
                    list<pair<size_t, size_t>> level_segment_begin;
                    level_segment_begin.emplace_back(0, 0);
                    
                    // we record which segments we are going to cut out of the match here
                    vector<pair<size_t, size_t>> cut_segments;
                    
                    auto curr_level = level_segment_begin.begin();
                    size_t prefix_length = 0;
                    for (size_t j = 0, last = path->mapping_size() - 1; j <= last; j++) {
                        const Position& position = path->mapping(j).position();
                        const auto& projection = projection_trans.at(position.node_id());
                        id_t projected_id = projection.first;
                        id_t projected_rev = projection.second != position.is_reverse();
                        
                        if (j > 0) {
                            // we have entered this node on this iteration
                            if (cutting_snarls->into_which_snarl(projected_id, !projected_rev)) {
                                // as we enter this node, we are leaving the snarl we were in
                                
                                // since we're going up a level, we need to check whether we need to cut out the segment we've traversed
                                if (prefix_length - curr_level->first <= max_snarl_cut_size) {
                                    cut_segments.emplace_back(curr_level->second, j);
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
                        prefix_length += mapping_from_length(path->mapping(j));
                        
                        if (j < last) {
                            // we are going to leave this node next iteration
                            if (cutting_snarls->into_which_snarl(projected_id, projected_rev)) {
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
                    if (prefix_length - curr_level->first <= max_snarl_cut_size && curr_level != last) {
                        cut_segments.emplace_back(curr_level->second, path->mapping_size());
                    }
                    
                    // did we cut out any segments?
                    if (!cut_segments.empty()) {
#ifdef debug_multipath_mapper
                        cerr << "found cut segments:" << endl;
                        for (auto seg : cut_segments) {
                            cerr << "\t" << seg.first << " " << seg.second << endl;
                        }
#endif
                        
                        // we may have decided to cut the segments of both a parent and child snarl, so now we
                        // collapse the list of intervals, which is sorted on the end index by construction
                        //
                        // snarl nesting properties guarantee that there will be at least one node between any
                        // cut segments that are not nested, so we don't need to deal with the case where the
                        // segments are partially overlapping (i.e. it's a bit easier than the general interval
                        // intersection problem)
                        vector<pair<size_t, size_t>> keep_segments;
                        size_t curr_keep_seg_end = path->mapping_size();
                        auto riter = cut_segments.rbegin();
                        if (riter->second == curr_keep_seg_end) {
                            // don't add an empty keep segment in the first position
                            curr_keep_seg_end = riter->first;
                            riter++;
                        }
                        for (; riter != cut_segments.rend(); riter++) {
                            if (riter->second < curr_keep_seg_end) {
                                // this is a new interval
                                keep_segments.emplace_back(riter->second, curr_keep_seg_end);
                                curr_keep_seg_end = riter->first;
                            }
                        }
                        if (curr_keep_seg_end > 0) {
                            // we are not cutting off the left tail, so add a keep segment for it
                            keep_segments.emplace_back(0, curr_keep_seg_end);
                        }

                        // make a new node for all but one of the keep segments
                        size_t prefix_length = 0;
                        size_t prefix_idx = 0;
                        for (auto iter = keep_segments.end() - 1; iter != keep_segments.begin(); iter--) {
#ifdef debug_multipath_mapper
                            cerr << "making path for keep segment " << iter->first << " " << iter->second << endl;
#endif
                            match_nodes.emplace_back();
                            
                            // update pointers in case the vector reallocates
                            match_node = &match_nodes[match_node_idx];
                            path = &match_node->path;
                            
                            // measure the length of path between this keep segment and the last one
                            while (prefix_idx < iter->first) {
                                prefix_length += mapping_from_length(path->mapping(prefix_idx));
                                prefix_idx++;
                            }
                            
                            size_t keep_segment_length = 0;
                            ExactMatchNode& cut_node = match_nodes.back();
                            Path& cut_path = cut_node.path;
                            // transfer over the keep segment from the main path and measure the length
                            int32_t rank = 1;
                            for (size_t j = iter->first; j < iter->second; j++, rank++) {
                                Mapping* mapping = cut_path.add_mapping();
                                *mapping = path->mapping(j);
                                mapping->set_rank(rank);
                                keep_segment_length += mapping_from_length(*mapping);
                            }
                            
                            // identify the substring of the MEM that stays on this node
                            cut_node.begin = match_node->begin + prefix_length;
                            cut_node.end = cut_node.begin + keep_segment_length;
                            
                            prefix_length += keep_segment_length;
                            prefix_idx = iter->second;
#ifdef debug_multipath_mapper
                            cerr << "new cut path: " << pb2json(cut_path) << endl;
#endif
                            
                        }
                        
                        while (prefix_idx < keep_segments.front().first) {
                            prefix_length += mapping_from_length(path->mapping(prefix_idx));
                            prefix_idx++;
                        }
                        
                        // replace the path of the original node with the final keep segment
                        size_t keep_segment_length = 0;
                        Path new_path;
                        int32_t rank = 1;
                        for (size_t j = keep_segments.front().first; j < keep_segments.front().second; j++, rank++) {
                            Mapping* mapping = new_path.add_mapping();
                            *mapping = path->mapping(j);
                            mapping->set_rank(rank);
                            keep_segment_length += mapping_from_length(*mapping);
                        }
                        *path = new_path;
                        
                        // update the substring of MEM to match the path
                        match_node->begin += prefix_length;
                        match_node->end = match_node->begin + keep_segment_length;
#ifdef debug_multipath_mapper
                        cerr << "new in place cut path: " << pb2json(new_path) << endl;
#endif
                    }
                }
            }
        }
        
#ifdef debug_multipath_mapper
        cerr << "computing reachability" << endl;
#endif
        
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
        
        auto endpoint_offset = [&](size_t idx, bool end) {
            return end ? end_offset(idx) : start_offset(idx);
        };
        
        auto endpoint_node_id = [&](size_t idx, bool end) {
            return end ? end_node_id(idx) : start_node_id(idx);
        };

        // record the start and end node ids of every exact match
        unordered_map<id_t, vector<size_t>> exact_match_starts;
        unordered_map<id_t, vector<size_t>> exact_match_ends;
        for (size_t i = 0; i < match_nodes.size(); i++) {
            Path& path = match_nodes[i].path;
            exact_match_starts[path.mapping(0).position().node_id()].push_back(i);
            exact_match_ends[path.mapping(path.mapping_size() - 1).position().node_id()].push_back(i);
        }
        
#ifdef debug_multipath_mapper
        cerr << "recorded starts: " << endl;
        for (const auto& rec : exact_match_starts) {
            cerr << "\t" << rec.first << ": ";
            for (auto l : rec.second) {
                cerr << l << " ";
            }
            cerr << endl;
        }
        
        cerr << "recorded ends: " << endl;
        for (const auto& rec : exact_match_ends) {
            cerr << "\t" << rec.first << ": ";
            for (auto l : rec.second) {
                cerr << l << " ";
            }
            cerr << endl;
        }
#endif

        
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
            
#ifdef debug_multipath_mapper
            cerr << "DP step for node " << node_id << endl;
#endif
            
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
                
#ifdef debug_multipath_mapper
                cerr << "node " << node_id << " contains both starts and ends of MEMs" << endl;
#endif
                
                vector<size_t>& ends = exact_match_ends[node_id];
                vector<size_t>& starts = exact_match_starts[node_id];
                
                
                // find the range of starts and ends in the list with the same offset
                
                size_t start_range_begin = 0;
                size_t start_range_end = 0;
                size_t end_range_begin = 0;
                size_t end_range_end = 0;
                
                size_t curr_start_offset = start_offset(starts[start_range_begin]);
                size_t curr_end_offset = end_offset(ends[end_range_begin]);
                size_t prev_offset = 0;
                
                while (end_offset(ends[end_range_end]) == curr_end_offset) {
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
                
#ifdef debug_multipath_mapper
                cerr << "first starts in range " << start_range_begin << ":" << start_range_end << " are located at initial offset " << curr_start_offset << endl;
                cerr << "first ends in range " << end_range_begin << ":" << end_range_end << " are located at initial offset " << curr_end_offset << endl;
#endif
                
                // connect the first range of starts or ends to the incoming starts and ends
                
                size_t prev_end_range_begin = end_range_begin;
                size_t prev_start_range_begin = start_range_begin;
                
                bool at_end = (curr_end_offset <= curr_start_offset);
                unordered_map<id_t, unordered_map<size_t, size_t>>* reachable_endpoints;
                unordered_map<size_t, vector<pair<size_t, size_t>>>* reachable_ends_from_endpoint;
                unordered_map<size_t, vector<pair<size_t, size_t>>>* reachable_starts_from_endpoint;
                vector<size_t>* endpoints;
                size_t* range_begin;
                size_t* range_end;
                size_t* prev_range_begin;
                size_t* curr_offset;
                if (at_end) {
                    reachable_endpoints = &reachable_ends;
                    reachable_starts_from_endpoint = &reachable_starts_from_end;
                    reachable_ends_from_endpoint = &reachable_ends_from_end;
                    endpoints = &ends;
                    range_begin = &end_range_begin;
                    range_end = &end_range_end;
                    prev_range_begin = &prev_end_range_begin;
                    curr_offset = &curr_end_offset;
                    
#ifdef debug_multipath_mapper
                    cerr << "first endpoint is an end" << endl;
#endif
                }
                else {
                    reachable_endpoints = &reachable_starts;
                    reachable_starts_from_endpoint = &reachable_starts_from_start;
                    reachable_ends_from_endpoint = &reachable_ends_from_start;
                    endpoints = &starts;
                    range_begin = &start_range_begin;
                    range_end = &start_range_end;
                    prev_range_begin = &prev_start_range_begin;
                    curr_offset = &curr_start_offset;
                    
#ifdef debug_multipath_mapper
                    cerr << "first endpoint is a start" << endl;
#endif
                }
                
                for (size_t j = *range_begin; j < *range_end; j++) {
                    for (const pair<size_t, size_t>& incoming_end : reachable_ends[node_id]) {
#ifdef debug_multipath_mapper
                        cerr << "identifying end " << incoming_end.first << " as reachable from endpoint " << endpoints->at(j) << endl;
#endif
                        (*reachable_ends_from_endpoint)[endpoints->at(j)].emplace_back(incoming_end.first, incoming_end.second + *curr_offset);
                    }
                    for (const pair<size_t, size_t>& incoming_start : reachable_starts[node_id]) {
#ifdef debug_multipath_mapper
                        cerr << "identifying start " << incoming_start.first << " as reachable from endpoint " << endpoints->at(j) << endl;
#endif
                        (*reachable_starts_from_endpoint)[endpoints->at(j)].emplace_back(incoming_start.first, incoming_start.second + *curr_offset);
                    }
                }
                
                
                bool prev_is_end = at_end;
                *range_begin = *range_end;
                prev_offset = *curr_offset;
                if (*range_begin != endpoints->size()) {
                    *curr_offset = endpoint_offset(endpoints->at(*range_begin), at_end);
                    while (endpoint_offset(endpoints->at(*range_end), at_end) == *curr_offset) {
                        (*range_end)++;
                        if (*range_end == endpoints->size()) {
                            break;
                        }
                    }
                }
                
#ifdef debug_multipath_mapper
                cerr << "next range is " << *range_begin << ":" << *range_end << " at offset " << *curr_offset << endl;
#endif
                
                // iterate along ranges of starts or ends in order of their offsets
                
                while (start_range_begin < starts.size() && end_range_begin < ends.size()) {
                    at_end = (curr_end_offset <= curr_start_offset);
                    if (at_end) {
                        reachable_endpoints = &reachable_ends;
                        reachable_starts_from_endpoint = &reachable_starts_from_end;
                        reachable_ends_from_endpoint = &reachable_ends_from_end;
                        endpoints = &ends;
                        range_begin = &end_range_begin;
                        range_end = &end_range_end;
                        prev_range_begin = &prev_end_range_begin;
                        curr_offset = &curr_end_offset;
                    }
                    else {
                        reachable_endpoints = &reachable_starts;
                        reachable_starts_from_endpoint = &reachable_starts_from_start;
                        reachable_ends_from_endpoint = &reachable_ends_from_start;
                        endpoints = &starts;
                        range_begin = &start_range_begin;
                        range_end = &start_range_end;
                        prev_range_begin = &prev_start_range_begin;
                        curr_offset = &curr_start_offset;
                    }
                    
                    size_t dist_between = *curr_offset - prev_offset;
                    
                    // connect this range to the previous range
                    if (prev_is_end) {
                        for (size_t j = prev_end_range_begin; j < end_range_begin; j++) {
                            for (size_t k = *range_begin; k < *range_end; k++) {
#ifdef debug_multipath_mapper
                                cerr << "identifying end " << ends[j] << " as reachable from endpoint " << endpoints->at(k) << endl;
#endif
                                (*reachable_ends_from_endpoint)[endpoints->at(k)].emplace_back(ends[j], dist_between);
                            }
                        }
                    }
                    else {
                        for (size_t j = prev_start_range_begin; j < start_range_begin; j++) {
                            for (size_t k = end_range_begin; k < end_range_end; k++) {
#ifdef debug_multipath_mapper
                                cerr << "identifying start " << starts[j] << " as reachable from endpoint " << endpoints->at(k) << endl;
#endif
                                (*reachable_starts_from_endpoint)[endpoints->at(k)].emplace_back(starts[j], dist_between);
                            }
                        }
                    }
                    
                    // record the properties of this range
                    *prev_range_begin = *range_begin;
                    prev_is_end = at_end;
                    
                    // advance to the next range
                    *range_begin = *range_end;
                    prev_offset = *curr_offset;
                    if (*range_begin != endpoints->size()) {
                        *curr_offset = endpoint_offset(endpoints->at(*range_begin), at_end);
                        while (endpoint_offset(endpoints->at(*range_end), at_end) == *curr_offset) {
                            (*range_end)++;
                            if (*range_end == endpoints->size()) {
                                break;
                            }
                        }
                    }
                    
#ifdef debug_multipath_mapper
                    cerr << "next range is " << *range_begin << ":" << *range_end << " at offset " << *curr_offset << endl;
#endif
                }
                
                // finish off the list of starts or ends on this node
                
                at_end = (end_range_begin < ends.size());
                if (at_end) {
                    reachable_endpoints = &reachable_ends;
                    reachable_starts_from_endpoint = &reachable_starts_from_end;
                    reachable_ends_from_endpoint = &reachable_ends_from_end;
                    endpoints = &ends;
                    range_begin = &end_range_begin;
                    range_end = &end_range_end;
                    prev_range_begin = &prev_end_range_begin;
                    curr_offset = &curr_end_offset;
                    
#ifdef debug_multipath_mapper
                    cerr << "final endpoint(s) are end(s)" << endl;
#endif
                }
                else {
                    reachable_endpoints = &reachable_starts;
                    reachable_starts_from_endpoint = &reachable_starts_from_start;
                    reachable_ends_from_endpoint = &reachable_ends_from_start;
                    endpoints = &starts;
                    range_begin = &start_range_begin;
                    range_end = &start_range_end;
                    prev_range_begin = &prev_start_range_begin;
                    curr_offset = &curr_start_offset;
                    
#ifdef debug_multipath_mapper
                    cerr << "final endpoint(s) are start(s)" << endl;
#endif
                }
                
                while (*range_begin < endpoints->size()) {
                    
                    size_t dist_between = *curr_offset - prev_offset;
                    
                    if (prev_is_end) {
                        for (size_t j = prev_end_range_begin; j < end_range_begin; j++) {
                            for (size_t k = *range_begin; k < *range_end; k++) {
#ifdef debug_multipath_mapper
                                cerr << "identifying end " << ends[j] << " as reachable from endpoint " << endpoints->at(k) << endl;
#endif
                                (*reachable_ends_from_endpoint)[endpoints->at(k)].push_back(make_pair(ends[j], dist_between));
                            }
                        }
                    }
                    else {
                        for (size_t j = prev_start_range_begin; j < start_range_begin; j++) {
                            for (size_t k = *range_begin; k < *range_end; k++) {
#ifdef debug_multipath_mapper
                                cerr << "identifying start " << starts[j] << " as reachable from endpoint " << endpoints->at(k) << endl;
#endif
                                (*reachable_starts_from_endpoint)[endpoints->at(k)].push_back(make_pair(starts[j], dist_between));
                            }
                        }
                    }
                    
#ifdef debug_multipath_mapper
                    cerr << "moving to next endpoint range" << endl;
#endif
                    
                    *prev_range_begin = *range_begin;
                    *range_begin = *range_end;
                    prev_offset = *curr_offset;
                    prev_is_end = at_end;
                    
                    if (*range_begin != endpoints->size()) {
                        *curr_offset = endpoint_offset(endpoints->at(*range_begin), at_end);
                        while (endpoint_offset(endpoints->at(*range_end), at_end) == *curr_offset) {
                            (*range_end)++;
                            if (*range_end == endpoints->size()) {
                                break;
                            }
                        }
                    }
                    
#ifdef debug_multipath_mapper
                    cerr << "next range is " << *range_begin << ":" << *range_end << " at offset " << *curr_offset << endl;
#endif
                }
                
                // carry forward the reachability of the last range onto the next nodes
                size_t dist_thru = node_length - *curr_offset;
                
#ifdef debug_multipath_mapper
                cerr << "carrying forward reachability onto next nodes at distance " << dist_thru << endl;
#endif
                
                for (NodeTraversal next : nexts) {
                    unordered_map<size_t, size_t>& reachable_endpoints_next = (*reachable_endpoints)[next.node->id()];
                    for (size_t j = *prev_range_begin; j < endpoints->size(); j++) {
                        if (reachable_endpoints_next.count(endpoints->at(j))) {
                            reachable_endpoints_next[endpoints->at(j)] = std::min(reachable_endpoints_next[endpoints->at(j)], dist_thru);
                        }
                        else {
                            reachable_endpoints_next[endpoints->at(j)] = dist_thru;
                        }
                    }
                }
            }
            else if (contains_starts || contains_ends) {
                // this nodes contains at least one start or end, but either all starts or all ends
                // record the incoming starts/ends for the starts/ends on this node
                
                unordered_map<id_t, unordered_map<size_t, size_t>>* reachable_endpoints;
                unordered_map<size_t, vector<pair<size_t, size_t>>>* reachable_ends_from_endpoint;
                unordered_map<size_t, vector<pair<size_t, size_t>>>* reachable_starts_from_endpoint;
                unordered_map<size_t, vector<pair<size_t, size_t>>>* reachable_endpoints_from_endpoint;
                vector<size_t>* endpoints;
                if (contains_ends) {
#ifdef debug_multipath_mapper
                    cerr << "node " << node_id << " contains only ends of MEMs" << endl;
#endif
                    reachable_endpoints = &reachable_ends;
                    reachable_starts_from_endpoint = &reachable_starts_from_end;
                    reachable_ends_from_endpoint = &reachable_ends_from_end;
                    reachable_endpoints_from_endpoint = &reachable_ends_from_end;
                    endpoints = &exact_match_ends[node_id];
                }
                else {
#ifdef debug_multipath_mapper
                    cerr << "node " << node_id << " contains only starts of MEMs" << endl;
#endif
                    reachable_endpoints = &reachable_starts;
                    reachable_starts_from_endpoint = &reachable_starts_from_start;
                    reachable_ends_from_endpoint = &reachable_ends_from_start;
                    reachable_endpoints_from_endpoint = &reachable_starts_from_start;
                    endpoints = &exact_match_starts[node_id];
                }
                
                // the starts/ends coming into this node from outside
                size_t range_begin = 0;
                size_t range_end = 0;
                size_t curr_offset = endpoint_offset(endpoints->at(range_begin), contains_ends);
                size_t prev_offset = curr_offset;
                // find the range of endpoints that are at the first offset
                while (endpoint_offset(endpoints->at(range_end), contains_ends)) {
                    range_end++;
                    if (range_end == endpoints->size()) {
                        break;
                    }
                }
                // connect the range to the incoming starts/ends
                for (size_t j = range_begin; j < range_end; j++) {
                    for (const pair<size_t, size_t>& incoming_start : reachable_starts[node_id]) {
#ifdef debug_multipath_mapper
                        cerr << "identifying start " << incoming_start.first << " as reachable from endpoint " << endpoints->at(j) << endl;
#endif
                        (*reachable_starts_from_endpoint)[endpoints->at(j)].emplace_back(incoming_start.first, incoming_start.second + curr_offset);
                    }
                    for (const pair<size_t, size_t>& incoming_end : reachable_ends[node_id]) {
#ifdef debug_multipath_mapper
                        cerr << "identifying end " << incoming_end.first << " as reachable from endpoint " << endpoints->at(j) << endl;
#endif
                        (*reachable_ends_from_endpoint)[endpoints->at(j)].emplace_back(incoming_end.first, incoming_end.second + curr_offset);
                    }
                }
                
                // the reachable endpoints internal to this node
                size_t prev_range_begin = range_begin;
                range_begin = range_end;
                while (range_begin < endpoints->size()) {
                    // find the range of endpoints at this offset
                    prev_offset = curr_offset;
                    curr_offset = endpoint_offset(endpoints->at(range_begin), contains_ends);
                    while (endpoint_offset(endpoints->at(range_begin), contains_ends) == curr_offset) {
                        range_end++;
                        if (range_end == endpoints->size()) {
                            break;
                        }
                    }
                    
                    size_t dist_between = curr_offset - prev_offset;
                    
                    // connect this range to the previous range
                    for (size_t j = range_begin; j < range_end; j++) {
                        for (size_t k = prev_range_begin; k < range_begin; k++) {
#ifdef debug_multipath_mapper
                            cerr << "identifying endpoint " << endpoints->at(k) << " as reachable from endpoint " << endpoints->at(j) << endl;
#endif
                            reachable_endpoints_from_endpoint->at(endpoints->at(j)).push_back(make_pair(endpoints->at(k), dist_between));
                        }
                    }
                    prev_range_begin = range_begin;
                    range_begin = range_end;
                }
                
                // this node contains at least one endpoint of a MEM, so carry forward the reachability of all
                // endpoints at the final offset onto the next nodes
                
                size_t dist_thru = node_length - curr_offset;
                
#ifdef debug_multipath_mapper
                cerr << "carrying forward reachability onto next nodes at distance " << dist_thru << endl;
#endif
                
                for (NodeTraversal next : nexts) {

                    unordered_map<size_t, size_t>& reachable_endpoints_next = (*reachable_endpoints)[next.node->id()];
                    for (size_t j = prev_range_begin; j < endpoints->size(); j++) {
                        if (reachable_endpoints_next.count(endpoints->at(j))) {
                            reachable_endpoints_next[endpoints->at(j)] = std::min(reachable_endpoints_next[endpoints->at(j)], dist_thru);
                        }
                        else {
                            reachable_endpoints_next[endpoints->at(j)] = dist_thru;
                        }
                    }
                }
            }
            else {
                // this node doesn't contain the start or end of any MEM, so we carry forward the reachability
                // into this node onto the next nodes
                
#ifdef debug_multipath_mapper
                cerr << "node " << node_id << " does not contain starts or ends of MEMs, carrying forward reachability" << endl;
#endif
                
                for (NodeTraversal next : nexts) {
                    unordered_map<size_t, size_t>& reachable_ends_next = reachable_ends[next.node->id()];
                    for (const pair<size_t, size_t>& reachable_end : reachable_ends[node_id]) {
                        size_t dist_thru = reachable_end.second + node_length;
#ifdef debug_multipath_mapper
                        cerr << "\tend " << reachable_end.first << " at dist " << dist_thru << " to node " << next.node->id() << endl;
#endif
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
#ifdef debug_multipath_mapper
                        cerr << "\tstart " << reachable_start.first << " at dist " << dist_thru << " to node " << next.node->id() << endl;
#endif
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
        
#ifdef debug_multipath_mapper
        cerr << "final reachability:" << endl;
        cerr << "\tstarts from starts" << endl;
        for (const auto& record : reachable_starts_from_start) {
            cerr << "\t\t" << record.first << ":";
            for (const auto& endpoint : record.second) {
                cerr << " " << endpoint.first << "(" << endpoint.second << ")";
            }
            cerr << endl;
        }
        cerr << "\tends from starts" << endl;
        for (const auto& record : reachable_ends_from_start) {
            cerr << "\t\t" << record.first << ":";
            for (const auto& endpoint : record.second) {
                cerr << " " << endpoint.first << "(" << endpoint.second << ")";
            }
            cerr << endl;
        }
        cerr << "\tstarts from ends" << endl;
        for (const auto& record : reachable_starts_from_end) {
            cerr << "\t\t" << record.first << ":";
            for (const auto& endpoint : record.second) {
                cerr << " " << endpoint.first << "(" << endpoint.second << ")";
            }
            cerr << endl;
        }
        cerr << "\tends from ends" << endl;
        for (const auto& record : reachable_ends_from_end) {
            cerr << "\t\t" << record.first << ":";
            for (const auto& endpoint : record.second) {
                cerr << " " << endpoint.first << "(" << endpoint.second << ")";
            }
            cerr << endl;
        }
        cerr << "setting up structure of MEM graph" << endl;
#endif
        
        // now we have the reachability information for the start and end of every MEM in the graph. we
        // will use this to navigate between the MEMs in a way that respects graph reachability so that this
        // phase of the algorithm only needs to pay attention to read colinearity and transitive reducibility
        
        vector<unordered_map<size_t, size_t>> noncolinear_shells(match_nodes.size());
        
        // tuples of (overlap size, index from, index onto, dist)
        vector<tuple<size_t, size_t, size_t, size_t>> confirmed_overlaps;
        
        for (size_t i = 0; i < graph.node_size(); i++) {
            id_t node_id = graph.node(i).id();
            
#ifdef debug_multipath_mapper
            cerr << "looking for edges for starts on node " << node_id << endl;
#endif
            
            if (!exact_match_starts.count(node_id)) {
#ifdef debug_multipath_mapper
                cerr << "there are no starts on this node" << endl;
#endif
                continue;
            }
            
            for (size_t start : exact_match_starts[node_id]) {
                // traverse all of the reachable starts to find the adjacent ends that might be colinear
#ifdef debug_multipath_mapper
                cerr << "at start " << start << endl;
#endif
                
                ExactMatchNode& start_node = match_nodes[start];
                unordered_map<size_t, size_t>& noncolinear_shell = noncolinear_shells[start];
                
                // pairs of (dist, index)
                priority_queue<pair<size_t, size_t>, vector<pair<size_t, size_t>>, std::greater<pair<size_t, size_t>>> start_queue;
                priority_queue<pair<size_t, size_t>, vector<pair<size_t, size_t>>, std::greater<pair<size_t, size_t>>> end_queue;
                start_queue.emplace(0, start);
                
                unordered_set<size_t> traversed_start;
                
                while (!start_queue.empty()) {
                    pair<size_t, size_t> start_here = start_queue.top();
                    start_queue.pop();
                    if (traversed_start.count(start_here.second)) {
                        continue;
                    }
                    traversed_start.insert(start_here.second);
                    
                    // the minimum distance to each of the starts or ends this can reach is the sum of the min distance
                    // between them and the distance already traversed
                    for (const pair<size_t, size_t>& end : reachable_ends_from_start[start_here.second]) {
                        end_queue.emplace(start_here.first + end.second, end.first);
#ifdef debug_multipath_mapper
                        cerr << "found reachable end " << end.first << " at distance " << start_here.first + end.second << endl;
#endif
                    }
                    
                    for (const pair<size_t, size_t>& start_next : reachable_starts_from_start[start_here.second]) {
                        start_queue.emplace(start_here.first + start_next.second, start_next.first);
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
#ifdef debug_multipath_mapper
                    cerr << "considering end " << candidate_end << " as candidate for edge of dist " << candidate_dist << endl;
#endif
                    traversed_end.insert(candidate_end);
                    
                    ExactMatchNode& candidate_end_node = match_nodes[candidate_end];
                    
                    if (candidate_end_node.end <= start_node.begin) {
                        // these MEMs are read colinear and graph reachable, so connect them
                        candidate_end_node.edges.push_back(make_pair(start, candidate_dist));
                        
#ifdef debug_multipath_mapper
                        cerr << "connection is read colinear, adding edge" << endl;
#endif
                        
                        // skip to the predecessor's noncolinear shell, whose connections might not be blocked by
                        // this connection
                        for (const pair<size_t, size_t>& shell_pred : noncolinear_shells[candidate_end]) {
#ifdef debug_multipath_mapper
                            cerr << "enqueueing " << shell_pred.first << " at dist " << shell_pred.second << " from noncolinear shell" << endl;
#endif
                            end_queue.emplace(candidate_dist + shell_pred.second, shell_pred.first);
                        }
                    }
                    else if (start_node.end > candidate_end_node.end && candidate_end_node.begin < start_node.begin) {
                        // the MEM can be made colinear by removing an overlap, which will not threaten reachability
                        size_t overlap = candidate_end_node.end - start_node.begin;
                        confirmed_overlaps.emplace_back(overlap, candidate_end, start, candidate_dist + overlap);
                        
#ifdef debug_multipath_mapper
                        cerr << "connection is overlap colinear, recording to add edge later" << endl;
#endif
                        
                        // skip to the predecessor's noncolinear shell, whose connections might not be blocked by
                        // this connection
                        for (const pair<size_t, size_t>& shell_pred : noncolinear_shells[candidate_end]) {
#ifdef debug_multipath_mapper
                            cerr << "enqueueing " << shell_pred.first << " at dist " << shell_pred.second << " from noncolinear shell" << endl;
#endif
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
                        
#ifdef debug_multipath_mapper
                        cerr << "connection is noncolinear, add to shell and continue to search backwards" << endl;
#endif
                        
                        // there is no connection to block further connections back, so any of this MEMs
                        // predecessors could still be colinear
                        
                        // find the ends that can reach it directly
                        for (const pair<size_t, size_t>& pred_end : reachable_ends_from_end[candidate_end]) {
                            end_queue.emplace(candidate_dist + pred_end.second, pred_end.first);
#ifdef debug_multipath_mapper
                            cerr << "found reachable end " << pred_end.first << " at distance " << candidate_dist + pred_end.second << endl;
#endif
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
#ifdef debug_multipath_mapper
                                cerr << "found reachable end " << pred_end.first << " at distance " << candidate_dist + pred_end.second << endl;
#endif
                            }
                            for (const pair<size_t, size_t>& start_next : reachable_starts_from_start[start_here]) {
                                pred_start_queue.emplace(start_dist + start_next.second, start_next.first);
#ifdef debug_multipath_mapper
                                cerr << "traverse intermediate start " << start_next.first << " at distance " << start_dist + start_next.second << endl;
#endif
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
        
#ifdef debug_multipath_mapper
        cerr << "breaking nodes at overlap edges" << endl;
#endif
        
        // now we've found all overlap edges, so we can add them into the graph in an order such that they don't
        // conflict (note that all overlap are from an earlier node onto a later one, so we don't need to worry
        // about overlaps coming in from both directions)
        
        // sort in descending order of overlap length
        std::sort(confirmed_overlaps.begin(), confirmed_overlaps.end(),
                  std::greater<tuple<size_t, size_t, size_t, size_t>>());
        
        // split up each node with an overlap edge onto it
        for (auto overlap_record : confirmed_overlaps) {
#ifdef debug_multipath_mapper
            cerr << "performing an overlap split from node " << get<1>(overlap_record) << " onto " << get<2>(overlap_record) << " of length " << get<0>(overlap_record) << " at distance " << get<3>(overlap_record) << endl;
#endif
            
            size_t suffix_idx = match_nodes.size();
            
            match_nodes.emplace_back();
            ExactMatchNode& from_node = match_nodes[get<1>(overlap_record)];
            ExactMatchNode& onto_node = match_nodes[get<2>(overlap_record)];
            ExactMatchNode& suffix_node = match_nodes.back();
            
#ifdef debug_multipath_mapper
            cerr << "before splitting:" << endl;
            cerr << "from node:" << endl << "\t";
            for (auto iter = from_node.begin; iter != from_node.end; iter++) {
                cerr << *iter;
            }
            cerr << endl;
            cerr << "\t" << pb2json(from_node.path) << endl;
            cerr << "onto node:" << endl << "\t";
            for (auto iter = onto_node.begin; iter != onto_node.end; iter++) {
                cerr << *iter;
            }
            cerr << endl << "\t" << pb2json(onto_node.path) << endl;
#endif
            
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
                remaining -= mapping_len;
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
                suffix_split->mutable_position()->set_offset(split_mapping.position().offset() + remaining);
                Edit* suffix_edit = suffix_split->add_edit();
                suffix_edit->set_from_length(mapping_len - remaining);
                suffix_edit->set_to_length(mapping_len - remaining);
                
                mapping_idx++;
            }
            
            // add the remaining mappings to the suffix node
            for (; mapping_idx < full_path.mapping_size(); mapping_idx++) {
                *suffix_node.path.add_mapping() = full_path.mapping(mapping_idx);
            }
            
#ifdef debug_multipath_mapper
            cerr << "after splitting:" << endl;
            cerr << "from node:" << endl << "\t";
            for (auto iter = from_node.begin; iter != from_node.end; iter++) {
                cerr << *iter;
            }
            cerr << endl << "\t" << pb2json(from_node.path) << endl;
            cerr << "onto node:" << endl << "\t";
            for (auto iter = onto_node.begin; iter != onto_node.end; iter++) {
                cerr << *iter;
            }
            cerr << endl << "\t" << pb2json(onto_node.path) << endl;
            cerr << "suffix node:" << endl << "\t";
            for (auto iter = suffix_node.begin; iter != suffix_node.end; iter++) {
                cerr << *iter;
            }
            cerr << endl << "\t" << pb2json(suffix_node.path) << endl;
#endif
        }
        
#ifdef debug_multipath_mapper
        cerr << "final mem graph:" << endl;
        for (size_t i = 0; i < match_nodes.size(); i++) {
            ExactMatchNode& match_node = match_nodes[i];
            cerr << i << " " << pb2json(match_node.path) << " ";
            for (auto iter = match_node.begin; iter != match_node.end; iter++) {
                cerr << *iter;
            }
            cerr << endl;
            cerr << "\t";
            for (auto edge : match_node.edges) {
                cerr << "(to:" << edge.first << ", dist:" << edge.second << ") ";
            }
            cerr << endl;
        }
#endif
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
    
    void MultipathAlignmentGraph::reorder_adjacency_lists(const vector<size_t>& order) {
        vector<vector<pair<size_t, size_t>>> reverse_graph(match_nodes.size());
        for (size_t i = 0; i < match_nodes.size(); i++) {
            for (const pair<size_t, size_t>& edge : match_nodes[i].edges) {
                reverse_graph[edge.first].emplace_back(i, edge.second);
            }
        }
        for (ExactMatchNode& match_node : match_nodes) {
            size_t out_degree = match_node.edges.size();
            match_node.edges.clear();
            match_node.edges.reserve(out_degree);
        }
        for (size_t i : order) {
            for (const pair<size_t, size_t>& edge : reverse_graph[i]) {
                match_nodes[edge.first].edges.emplace_back(i, edge.second);
            }
        }
    }
    
    void MultipathAlignmentGraph::remove_transitive_edges(const vector<size_t>& topological_order) {
        // algorithm assumes edges are also sorted in topological order, which guarantees that we will
        // traverse a path that reveals an edge as transitive before actually traversing the transitive edge
        reorder_adjacency_lists(topological_order);
        
        for (size_t i : topological_order) {
            vector<pair<size_t, size_t>>& edges = match_nodes[i].edges;
            
            // if there is only one edge out of a node, that edge can never be transitive
            // (this optimization covers most cases)
            if (edges.size() <= 1) {
                continue;
            }
            
            vector<bool> keep(edges.size(), true);
            unordered_set<size_t> traversed;
            
            for (size_t j = 0; j < edges.size(); j++) {
                const pair<size_t, size_t>& edge = edges[j];
                if (traversed.count(edge.first)) {
                    // we can reach the target of this edge by another path, so it is transitive
                    keep[j] = false;
                    continue;
                }
                
                // DFS to mark all reachable nodes from this edge
                vector<size_t> stack{edge.first};
                traversed.insert(edge.first);
                while (!stack.empty()) {
                    size_t idx = stack.back();
                    stack.pop_back();
                    for (const pair<size_t, size_t>& edge_from : match_nodes[idx].edges) {
                        if (!traversed.count(edge_from.first)) {
                            stack.push_back(edge_from.first);
                            traversed.insert(edge_from.first);
                        }
                    }
                }
            }
            
            // remove the transitive edges we found
            size_t next_idx = 0;
            for (size_t j = 0; j < edges.size(); j++) {
                if (keep[j] && j != next_idx) {
                    edges[next_idx] = edges[j];
                    next_idx++;
                }
                else if (keep[j]) {
                    next_idx++;
                }
            }
            edges.resize(next_idx);
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
            int32_t from_score = forward_scores[idx];
            for (const pair<size_t, size_t>& edge : match_nodes[idx].edges) {
                forward_scores[edge.first] = std::max(forward_scores[edge.first],
                                                      node_weights[edge.first] + from_score + edge_weights[make_pair(idx, edge.first)]);
            }
        }
        
        // backward DP
        for (int64_t i = topological_order.size() - 1; i >= 0; i--) {
            size_t idx = topological_order[i];
            int32_t score_here = node_weights[idx];
            for (const pair<size_t, size_t>& edge : match_nodes[idx].edges) {
                backward_scores[idx] = std::max(backward_scores[idx],
                                                score_here + backward_scores[edge.first] + edge_weights[make_pair(idx, edge.first)]);
            }
        }
        
        // compute the minimum score we will require of a node or edge
        int32_t min_path_score = *std::max_element(forward_scores.begin(), forward_scores.end()) - max_suboptimal_score_diff;
        
        // use forward-backward to find nodes/edges on some path with a score above the minimum
        unordered_set<size_t> keep_nodes;
        unordered_set<pair<size_t, size_t>> keep_edges;
        for (size_t i = 0; i < match_nodes.size(); i++) {
            if (forward_scores[i] + backward_scores[i] - node_weights[i] >= min_path_score) {
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
                                              [&](const pair<size_t, size_t>& edge) { return !keep_edges.count(make_pair(i, edge.first)); });
                match_node.edges.resize(new_end - match_node.edges.begin());
                next++;
            }
        }
        match_nodes.resize(next);
        
#ifdef debug_multipath_mapper
        cerr << "mem graph after pruning:" << endl;
        for (size_t i = 0; i < match_nodes.size(); i++) {
            ExactMatchNode& match_node = match_nodes[i];
            cerr << i << " " << pb2json(match_node.path) << " ";
            for (auto iter = match_node.begin; iter != match_node.end; iter++) {
                cerr << *iter;
            }
            cerr << endl;
            cerr << "\t";
            for (auto edge : match_node.edges) {
                cerr << "(to:" << edge.first << ", dist:" << edge.second << ") ";
            }
            cerr << endl;
        }
#endif
    }
}



