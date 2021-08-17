//
//  multipath_alignment_graph.cpp
//

#include "multipath_alignment_graph.hpp"
#include "sequence_complexity.hpp"

//#define debug_multipath_alignment

using namespace std;
namespace vg {
    
    unordered_multimap<id_t, pair<id_t, bool>> MultipathAlignmentGraph::create_injection_trans(const unordered_map<id_t, pair<id_t, bool>>& projection_trans) {
        // create the injection translator, which maps a node in the original graph to every one of its occurrences
        // in the dagified graph
        unordered_multimap<id_t, pair<id_t, bool> > injection_trans;
        for (const auto& trans_record : projection_trans) {
#ifdef debug_multipath_alignment
            cerr << trans_record.second.first << " -> " << trans_record.first << (trans_record.second.second ? "-" : "+") << endl;
#endif
            injection_trans.emplace(trans_record.second.first, make_pair(trans_record.first, trans_record.second.second));
        }
        
        return injection_trans;
    }

    function<pair<id_t, bool>(id_t)> MultipathAlignmentGraph::create_projector(const unordered_map<id_t, pair<id_t, bool>>& projection_trans) {
        return [&](id_t node_id) { return projection_trans.at(node_id); };
    }

    unordered_multimap<id_t, pair<id_t, bool>> MultipathAlignmentGraph::create_injection_trans(const HandleGraph& graph,
                                                                                               const function<pair<id_t, bool>(id_t)>& project) {
        unordered_multimap<id_t, pair<id_t, bool>> injection_trans;
        graph.for_each_handle([&](const handle_t& handle) {
            id_t node_id = graph.get_id(handle);
            auto proj = project(node_id);
            injection_trans.emplace(proj.first, make_pair(node_id, proj.second));
        });
        return injection_trans;
    }
    
    unordered_map<id_t, pair<id_t, bool>> MultipathAlignmentGraph::create_identity_projection_trans(const HandleGraph& graph) {
        unordered_map<id_t, pair<id_t, bool>> to_return;
        
        graph.for_each_handle([&](const handle_t& handle) {
            // Each node just projects from itself forward.
            to_return[graph.get_id(handle)] = make_pair(graph.get_id(handle), false);
        });
        
        return to_return;
    }
    
    MultipathAlignmentGraph::MultipathAlignmentGraph(const HandleGraph& graph,
                                                     const vector<pair<pair<string::const_iterator, string::const_iterator>, Path>>& path_chunks,
                                                     const Alignment& alignment, const function<pair<id_t, bool>(id_t)>& project,
                                                     const unordered_multimap<id_t, pair<id_t, bool>>& injection_trans, bool realign_Ns,
                                                     bool preserve_tail_anchors) {
        
        // Set up the initial multipath graph from the given path chunks.
        create_path_chunk_nodes(graph, path_chunks, alignment, project, injection_trans);
        
        // trim indels off of nodes to make the score dynamic programmable across nodes
        trim_hanging_indels(alignment, realign_Ns, preserve_tail_anchors);
        
        // compute reachability and add edges
        add_reachability_edges(graph, project, injection_trans);
        
    }
    
    MultipathAlignmentGraph::MultipathAlignmentGraph(const HandleGraph& graph,
                                                     const vector<pair<pair<string::const_iterator, string::const_iterator>, Path>>& path_chunks,
                                                     const Alignment& alignment, const function<pair<id_t, bool>(id_t)>& project, bool realign_Ns,
                                                     bool preserve_tail_anchors) :
                                                     MultipathAlignmentGraph(graph, path_chunks, alignment, project,
                                                                             create_injection_trans(graph, project), realign_Ns, preserve_tail_anchors) {
        // Nothing to do
        
    }

    MultipathAlignmentGraph::MultipathAlignmentGraph(const HandleGraph& graph,
                                                     const vector<pair<pair<string::const_iterator, string::const_iterator>, Path>>& path_chunks,
                                                     const Alignment& alignment, const unordered_map<id_t, pair<id_t, bool>>& projection_trans, bool realign_Ns,
                                                     bool preserve_tail_anchors) :
                                                     MultipathAlignmentGraph(graph, path_chunks, alignment, create_projector(projection_trans),
                                                                             create_injection_trans(projection_trans), realign_Ns, preserve_tail_anchors) {
        // Nothing to do
        
    }
    
    MultipathAlignmentGraph::MultipathAlignmentGraph(const HandleGraph& graph, MultipathMapper::memcluster_t& hits,
                                                     const function<pair<id_t, bool>(id_t)>& project,
                                                     const unordered_multimap<id_t, pair<id_t, bool>>& injection_trans,
                                                     vector<size_t>& path_node_provenance,
                                                     size_t max_branch_trim_length, gcsa::GCSA* gcsa,
                                                     const MultipathMapper::match_fanouts_t* fanout_breaks) {
        
        // initialize the match nodes
        create_match_nodes(graph, hits, project, injection_trans, path_node_provenance, max_branch_trim_length, fanout_breaks);
        
        if (gcsa) {
            // we indicated that these MEMs came from a GCSA, so there might be order-length MEMs that we can combine
            // TODO: this can lose some provenance information
            collapse_order_length_runs(graph, gcsa, path_node_provenance);
        }
        
        if (max_branch_trim_length) {
            // we indicated that we'd like to trim the path nodes to avoid ends that cause us to get locked into one
            // branch after a branch point
            trim_to_branch_points(&graph, max_branch_trim_length);
        }
        
#ifdef debug_multipath_alignment
        cerr << "nodes after adding, jittering, trimming, and collapsing:" << endl;
        for (size_t i = 0; i < path_nodes.size(); i++) {
            PathNode& path_node = path_nodes.at(i);
            cerr << i << " (hit " << path_node_provenance[i] << ") " << debug_string(path_node.path) << " ";
            for (auto iter = path_node.begin; iter != path_node.end; iter++) {
                cerr << *iter;
            }
            cerr << endl;
        }
#endif
        
        // compute reachability and add edges
        add_reachability_edges(graph, project, injection_trans, &path_node_provenance);
    }
    
    MultipathAlignmentGraph::MultipathAlignmentGraph(const HandleGraph& graph, MultipathMapper::memcluster_t& hits,
                                                     const function<pair<id_t, bool>(id_t)>& project,
                                                     vector<size_t>& path_node_provenance,
                                                     size_t max_branch_trim_length, gcsa::GCSA* gcsa,
                                                     const MultipathMapper::match_fanouts_t* fanout_breaks) :
                                                     MultipathAlignmentGraph(graph, hits, project,
                                                                             create_injection_trans(graph, project),
                                                                             path_node_provenance, max_branch_trim_length,
                                                                             gcsa, fanout_breaks) {
        // Nothing to do
        
    }

    MultipathAlignmentGraph::MultipathAlignmentGraph(const HandleGraph& graph, MultipathMapper::memcluster_t& hits,
                                                     const unordered_map<id_t, pair<id_t, bool>>& projection_trans,
                                                     vector<size_t>& path_node_provenance,
                                                     size_t max_branch_trim_length, gcsa::GCSA* gcsa,
                                                     const MultipathMapper::match_fanouts_t* fanout_breaks) :
                                                     MultipathAlignmentGraph(graph, hits, create_projector(projection_trans),
                                                                             create_injection_trans(projection_trans),
                                                                             path_node_provenance, max_branch_trim_length,
                                                                             gcsa, fanout_breaks) {
        // Nothing to do
        
    }
    
    MultipathAlignmentGraph::MultipathAlignmentGraph(const HandleGraph& graph, const Alignment& alignment, SnarlManager* snarl_manager,
                                                     MinimumDistanceIndex* dist_index, size_t max_snarl_cut_size,
                                                     const function<pair<id_t, bool>(id_t)>& project,
                                                     const unordered_multimap<id_t, pair<id_t, bool>>& injection_trans) {
        
        // this can only be done on aligned sequences
        if (!alignment.has_path() || alignment.path().mapping_size() == 0) {
            has_reachability_edges = true;
            return;
        }
        
        // shim the aligned path into the path chunks constructor to make a node for it
        vector<pair<pair<string::const_iterator, string::const_iterator>, Path>> path_holder;
        path_holder.emplace_back(make_pair(alignment.sequence().begin(), alignment.sequence().end()), alignment.path());
        create_path_chunk_nodes(graph, path_holder, alignment, project, injection_trans);
        
        // cut the snarls out of the aligned path so we can realign through them
        if (max_snarl_cut_size) {
            resect_snarls_from_paths(snarl_manager, dist_index, project, max_snarl_cut_size);
        }
        
        // the snarls algorithm adds edges where necessary
        has_reachability_edges = true;
        
        // trim indels from the end of path nodes so that scores will be dynamic programmable across subpaths
        trim_hanging_indels(alignment, true);
    }

    MultipathAlignmentGraph::MultipathAlignmentGraph(const HandleGraph& graph, const Alignment& alignment, SnarlManager* snarl_manager,
                                                     MinimumDistanceIndex* dist_index, size_t max_snarl_cut_size,
                                                     const unordered_map<id_t, pair<id_t, bool>>& projection_trans) :
                                                     MultipathAlignmentGraph(graph, alignment, snarl_manager, dist_index, max_snarl_cut_size,
                                                                             create_projector(projection_trans),
                                                                             create_injection_trans(projection_trans)) {
        // Nothing to do
    }

    MultipathAlignmentGraph::MultipathAlignmentGraph(const HandleGraph& graph, const Alignment& alignment, SnarlManager* snarl_manager,
                                                     MinimumDistanceIndex* dist_index, size_t max_snarl_cut_size,
                                                     const function<pair<id_t, bool>(id_t)>& project) :
                                                     MultipathAlignmentGraph(graph, alignment, snarl_manager, dist_index, max_snarl_cut_size,
                                                                             project, create_injection_trans(graph, project)) {
        // Nothing to do
    }
    
    void MultipathAlignmentGraph::create_path_chunk_nodes(const HandleGraph& graph, const vector<pair<pair<string::const_iterator, string::const_iterator>, Path>>& path_chunks,
                                                          const Alignment& alignment, const function<pair<id_t, bool>(id_t)>& project,
                                                          const unordered_multimap<id_t, pair<id_t, bool>>& injection_trans) {
        
        for (const auto& path_chunk : path_chunks) {
            
#ifdef debug_multipath_alignment
            cerr << "performing DFS to walk out path " << pb2json(path_chunk.second) << endl;
#endif
            
            const Path& path = path_chunk.second;
            
            auto range = injection_trans.equal_range(path.mapping(0).position().node_id());
            for (auto iter = range.first; iter != range.second; iter++) {
                
                id_t injected_id = iter->second.first;
                
                if (iter->second.second != path.mapping(0).position().is_reverse()) {
                    continue;
                }
                
                // stack for DFS, each record contains records of (next trav index, next traversals)
                vector<pair<size_t, vector<handle_t>>> stack;
                stack.emplace_back(0, vector<handle_t>{graph.get_handle(injected_id)});
                
                while (!stack.empty()) {
                    auto& back = stack.back();
                    if (back.first == back.second.size()) {
#ifdef debug_multipath_alignment
                        cerr << "traversed all edges out of current traversal" << endl;
#endif
                        stack.pop_back();
                        continue;
                    }
                    handle_t trav = back.second[back.first];
                    back.first++;
                    
#ifdef debug_multipath_alignment
                    cerr << "checking node " << graph.get_id(trav) << endl;
#endif
                    pair<id_t, bool> projected_trav = project(graph.get_id(trav));
                    
                    const Position& pos = path.mapping(stack.size() - 1).position();
                    if (projected_trav.first == pos.node_id() &&
                        projected_trav.second == (projected_trav.second != graph.get_is_reverse(trav))) {
                        
                        // position matched the path
                        
#ifdef debug_multipath_alignment
                        cerr << "chunk position " << pb2json(pos) << " matches traversal " << projected_trav.first << (projected_trav.second ? "-" : "+") << ", walked " << stack.size() << " of " << path.mapping_size() << " mappings" << endl;
#endif
                        
                        if (stack.size() == path.mapping_size()) {
#ifdef debug_multipath_alignment
                            cerr << "finished walking path" << endl;
#endif
                            break;
                        }
                        stack.emplace_back(0, vector<handle_t>());
                        graph.follow_edges(trav, false, [&](const handle_t& next) {
                            stack.back().second.emplace_back(next);
                        });
                    }
                }
                
                // did we successfully walk the path out?
                if (stack.empty()) {
#ifdef debug_multipath_alignment
                    cerr << "failed to successfully walk path, skipping" << endl;
#endif
                    continue;
                }
                
                // now we can make a node in the subpath graph
                path_nodes.emplace_back();
                PathNode& path_node = path_nodes.back();
                
                path_node.begin = path_chunk.first.first;
                path_node.end = path_chunk.first.second;
                
                // move over the path
                for (size_t i = 0; i < path.mapping_size(); i++) {
                    const Mapping& mapping = path.mapping(i);
                    const Position& position = mapping.position();
                    
                    auto& stack_record = stack[i];
                    handle_t& trav = stack_record.second[stack_record.first - 1];
                    
                    path_mapping_t* new_mapping = path_node.path.add_mapping();
                    position_t* new_position = new_mapping->mutable_position();
                                        
                    // use the node space that we walked out in
                    new_position->set_node_id(graph.get_id(trav));
                    new_position->set_is_reverse(graph.get_is_reverse(trav));
                    
                    new_position->set_offset(position.offset());
                    
                    for (int64_t j = 0; j < mapping.edit_size(); j++) {
                        from_proto_edit(mapping.edit(j), *new_mapping->add_edit());
                    }
                }
                
#ifdef debug_multipath_alignment
                cerr << "walked path: " << debug_string(path_node.path) << endl;
                cerr << "sequence: " << string(path_node.begin, path_node.end) << endl;
#endif
            }
        }
#ifdef debug_multipath_alignment
        cerr << "final path chunks:" << endl;
        for (size_t i = 0; i < path_nodes.size(); ++i) {
            cerr << i << ": " << string(path_nodes[i].begin, path_nodes[i].end) << " " << debug_string(path_nodes[i].path) << endl;
        }
#endif
    }
   
    
    bool MultipathAlignmentGraph::trim_and_check_for_empty(const Alignment& alignment, bool trim_Ns, PathNode& path_node,
                                                           bool preserve_tail_anchors, int64_t* removed_start_from_length,
                                                           int64_t* removed_end_from_length) {
        
#ifdef debug_multipath_alignment
        cerr << "trimming path node " << string(path_node.begin, path_node.end) << " " << debug_string(path_node.path) << endl;
#endif
        
        // Trim down the given PathNode of everything except softclips.
        // Return true if it all gets trimmed away and should be removed.
        path_t& path = path_node.path;
            
        int64_t mapping_start_idx = 0;
        int64_t mapping_last_idx = path.mapping_size() - 1;
        
        int64_t edit_start_idx = 0;
        int64_t edit_last_idx = path.mapping(mapping_last_idx).edit_size() - 1;
        
        // don't cut off softclips (we assume the entire softclip is in one edit and the next edit is aligned bases)
        bool softclip_start = (path.mapping(0).edit(0).from_length() == 0 &&
                               path.mapping(0).edit(0).to_length() > 0 &&
                               path_node.begin == alignment.sequence().begin());
        bool softclip_end = (path.mapping(mapping_last_idx).edit(edit_last_idx).from_length() == 0 &&
                             path.mapping(mapping_last_idx).edit(edit_last_idx).to_length() > 0 &&
                             path_node.end == alignment.sequence().end());
        
        // if indicated, we may want to preserve the location of tail anchors in spite of deletions
        bool ignore_deletion_start = (path.mapping(0).edit(0).from_length() > 0 &&
                                      path.mapping(0).edit(0).to_length() == 0 &&
                                      path_node.begin == alignment.sequence().begin() &&
                                      preserve_tail_anchors);
        bool ignore_deletion_end = (path.mapping(mapping_last_idx).edit(edit_last_idx).from_length() > 0 &&
                                    path.mapping(mapping_last_idx).edit(edit_last_idx).to_length() == 0 &&
                                    path_node.end == alignment.sequence().end() &&
                                    preserve_tail_anchors);
        
#ifdef debug_multipath_alignment
        cerr << "preserving softclips: begin? " << softclip_start << ", end? " << softclip_end << endl;
        cerr << "preserving deletion anchors: begin? " << ignore_deletion_start << ", end? " << ignore_deletion_end << endl;
#endif
        
        // Track how much we trim off each end
        int64_t removed_start_to_length = 0;
        int64_t removed_end_to_length = 0;
        if (removed_start_from_length != nullptr) {
            *removed_start_from_length = 0;
        }
        if (removed_end_from_length != nullptr) {
            *removed_end_from_length = 0;
        }
        
        int64_t removed_start_mapping_from_length = 0;
        
        // find the first aligned, non-N bases from the start of the path
        if (!softclip_start && !ignore_deletion_start) {
            bool found_start = false;
            for (; mapping_start_idx < path.mapping_size(); mapping_start_idx++) {
                const path_mapping_t& mapping = path.mapping(mapping_start_idx);
                removed_start_mapping_from_length = 0;
                for (edit_start_idx = 0; edit_start_idx < mapping.edit_size(); edit_start_idx++) {
                    const edit_t& edit = mapping.edit(edit_start_idx);
                    
                    if (edit.from_length() > 0 && edit.to_length() > 0 &&
                        (edit.sequence().empty() || !trim_Ns || any_of(edit.sequence().begin(), edit.sequence().end(), [](char c) {return c != 'N';}))) {
                        found_start = true;
                        break;
                    }
                    removed_start_to_length += edit.to_length();
                    removed_start_mapping_from_length += edit.from_length();
                }
                
                if (removed_start_from_length != nullptr) {
                    // Record how much we trimmed
                    *removed_start_from_length += removed_start_mapping_from_length;
                }
                
                if (found_start) {
                    break;
                }
            }
        }
        
        // find the first aligned bases from the end of the path
        if (!softclip_end && !ignore_deletion_end) {
            bool found_last = false;
            for (; mapping_last_idx >= 0; mapping_last_idx--) {
                const path_mapping_t& mapping = path.mapping(mapping_last_idx);
                for (edit_last_idx = mapping.edit_size() - 1; edit_last_idx >= 0; edit_last_idx--) {
                    const edit_t& edit = mapping.edit(edit_last_idx);
                    if (edit.from_length() > 0 && edit.to_length() > 0 &&
                        (edit.sequence().empty() || !trim_Ns || any_of(edit.sequence().begin(), edit.sequence().end(), [](char c) {return c != 'N';}))) {
                        found_last = true;
                        break;
                    }
                    removed_end_to_length += edit.to_length();
                    if (removed_end_from_length != nullptr) {
                        *removed_end_from_length += edit.to_length();
                    }
                }
                if (found_last) {
                    break;
                }
            }
        }
#ifdef debug_multipath_alignment
        cerr << "after removing non-softclip flanking indels and N matches, path goes from (" << mapping_start_idx << ", " << edit_start_idx << ") to (" << mapping_last_idx << ", " << edit_last_idx << ")" << endl;
#endif
        
        // did we find any indels?
        if (mapping_start_idx != 0 || mapping_last_idx + 1 != path.mapping_size() ||
            edit_start_idx != 0 || edit_last_idx + 1 != path.mapping(mapping_last_idx).edit_size()) {
            
            // would we need to trim the whole node?
            if (mapping_start_idx < mapping_last_idx ||
                (mapping_start_idx == mapping_last_idx && edit_start_idx <= edit_last_idx)) {
                
                // update the read interval
                path_node.begin += removed_start_to_length;
                path_node.end -= removed_end_to_length;
                
                // make a new path with the indels trimmed
                path_t trimmed_path;
                for (int64_t j = mapping_start_idx; j <= mapping_last_idx; j++) {
                    const path_mapping_t& mapping = path.mapping(j);
                    const position_t& position = mapping.position();
                    
                    path_mapping_t* new_mapping = trimmed_path.add_mapping();
                    position_t* new_position = new_mapping->mutable_position();
                    
                    new_position->set_node_id(position.node_id());
                    new_position->set_is_reverse(position.is_reverse());
                    new_position->set_offset(position.offset() + (j == mapping_start_idx ? removed_start_mapping_from_length : 0));
                    
                    size_t k_start = (j == mapping_start_idx ? edit_start_idx : 0);
                    size_t k_end = (j == mapping_last_idx ? edit_last_idx + 1 : mapping.edit_size());
                    for (size_t k = k_start; k < k_end; k++) {
                        *new_mapping->add_edit() = mapping.edit(k);
                    }
                }
                
                path_node.path = move(trimmed_path);
            }
            else if (ignore_deletion_start) {
                // we would need to remove the whole node, except we indicated that we want
                // to preserve the tail anchors to the start
                
                path_node.end = path_node.begin;
                
                pos_t start_pos = initial_position(path_node.path);
                path.clear_mapping();
                path_mapping_t* mapping = path.add_mapping();
                position_t* pos = mapping->mutable_position();
                pos->set_node_id(id(start_pos));
                pos->set_offset(offset(start_pos));
                pos->set_is_reverse(is_rev(start_pos));
                mapping->add_edit();
                
#ifdef debug_multipath_alignment
                cerr << "preserving start deletion path read[" << (path_node.begin - alignment.sequence().begin()) << "] " << debug_string(path_node.path) << endl;
#endif
            }
            else if (ignore_deletion_end) {
                // we would need to remove the whole node, except we indicated that we want
                // to preserve the tail anchors to the end
                
                path_node.begin = path_node.end;
                
                pos_t end_pos = final_position(path_node.path);
                path.clear_mapping();
                path_mapping_t* mapping = path.add_mapping();
                position_t* pos = mapping->mutable_position();
                pos->set_node_id(id(end_pos));
                pos->set_offset(offset(end_pos));
                pos->set_is_reverse(is_rev(end_pos));
                mapping->add_edit();
                
#ifdef debug_multipath_alignment
                cerr << "preserving end deletion path read[" << (path_node.begin - alignment.sequence().begin()) << "] " << debug_string(path_node.path) << endl;
#endif
            }
            else {
#ifdef debug_multipath_alignment
                cerr << "entire path node is trimmed" << endl;
#endif
                // We do need to remove the whole node; now it is empty.
                return true;
            }
        }
        
        // The node itself can stay
        return false;
    }
   
    void MultipathAlignmentGraph::trim_hanging_indels(const Alignment& alignment, bool trim_Ns,
                                                      bool preserve_tail_anchors) {
        
        // if the path begins or ends with any gaps we have to remove them to make the score
        // dynamic programmable across Subpaths
        
        unordered_set<size_t> to_remove;
        
        for (size_t i = 0; i < path_nodes.size(); i++) {
            
            PathNode& path_node = path_nodes.at(i);
            
            if (trim_and_check_for_empty(alignment, trim_Ns, path_node, preserve_tail_anchors)) {
                // We trimmed it and it all trimmed away
                to_remove.insert(i);
            }
        }
        
        if (!to_remove.empty()) {
            // we need remove the nodes that were completely trimmed from the graph
            
            // first we need to make edges that bypass the nodes we're going to remove
            for (size_t i = 0; i < path_nodes.size(); i++) {
                if (to_remove.count(i)) {
                    // we don't need to update edges on nodes we're going to remove
                    continue;
                }
                
                vector<pair<size_t, size_t>> new_edges;
                
                // records of (distance, index) in a queue for Dijkstra traversal
                priority_queue<pair<size_t, size_t>, vector<pair<size_t, size_t>>, greater<pair<size_t, size_t>>> edge_queue;
                // the indexes of edges we've already added
                unordered_set<size_t> added_edges;
                
                for (pair<size_t, size_t>& edge : path_nodes.at(i).edges) {
                    edge_queue.emplace(edge.second, edge.first);
                }
                
                while (!edge_queue.empty()) {
                    pair<size_t, size_t> traversed_edge = edge_queue.top();
                    edge_queue.pop();
                    
                    if (to_remove.count(traversed_edge.second)) {
                        // we're removing this path node, so traverse it and add connections along its edges
                        PathNode& removing_node = path_nodes.at(traversed_edge.second);
                        size_t through_length = traversed_edge.first + path_from_length(removing_node.path);
                        for (pair<size_t, size_t>& edge : removing_node.edges) {
                            edge_queue.emplace(through_length + edge.second, edge.first);
                        }
                    }
                    else if (!added_edges.count(traversed_edge.second)) {
                        // we've finished walking a shortest path to a non-removed node, switch the order back and add it
                        new_edges.emplace_back(traversed_edge.second, traversed_edge.first);
                        
                        // and mark it as added so we don't make duplicate edges
                        added_edges.insert(traversed_edge.second);
                    }
                }
                
                // replace the old edges with the new ones
                path_nodes.at(i).edges = new_edges;
            }
            
            // move the nodes we're going to keep into the prefix of the vector
            vector<size_t> removed_so_far(path_nodes.size(), 0);
            for (size_t i = 0; i < path_nodes.size(); i++) {
                if (i > 0) {
                    removed_so_far[i] = removed_so_far[i - 1];
                }
                
                if (to_remove.count(i)) {
                    removed_so_far[i]++;;
                }
                else if (removed_so_far[i]) {
                    path_nodes.at(i - removed_so_far[i]) = move(path_nodes.at(i));
                }
            }
            
            // actually remove the nodes
            path_nodes.resize(path_nodes.size() - to_remove.size());
            
            // update the indexes of the edges
            for (size_t i = 0; i < path_nodes.size(); i++) {
                for (pair<size_t, size_t>& edge : path_nodes.at(i).edges) {
                    edge.first -= removed_so_far[edge.first];
                }
            }
            
#ifdef debug_multipath_alignment
            cerr << "removed " << removed_so_far.back() << " complete nodes" << endl;
#endif
        }
        
#ifdef debug_multipath_alignment
        cerr << "finished trimming hanging indels" << endl;
#endif
    }
    
    void MultipathAlignmentGraph::create_match_nodes(const HandleGraph& graph, MultipathMapper::memcluster_t& hits,
                                                     const function<pair<id_t, bool>(id_t)>& project,
                                                     const unordered_multimap<id_t, pair<id_t, bool>>& injection_trans,
                                                     vector<size_t>& path_node_provenance,
                                                     int64_t max_branch_trim_length,
                                                     const MultipathMapper::match_fanouts_t* fanout_breaks) {
        
#ifdef debug_multipath_alignment
        cerr << "walking out MEMs in graph" << endl;
#endif
        
        // map of node ids in the dagified graph to the indices in the matches that contain them
        unordered_map<int64_t, vector<int64_t>> node_matches;
        
        // records an existing path node in this map
        auto record_node_matches = [&](const size_t i) {
            for (const auto& m : path_nodes[i].path.mapping()) {
                // record that each node occurs in this match so we can filter out sub-MEMs
                node_matches[m.position().node_id()].push_back(i);
#ifdef debug_multipath_alignment
                cerr << "associating node " << m.position().node_id() << " with a match at idx " << i << endl;
#endif
            }
        };
        
        // performs a check to see if a hit is a redundant sub-MEM
        auto is_redundant = [&](string::const_iterator begin, string::const_iterator end,
                                const pos_t& hit_pos, id_t injected_id) {
            
            // check all MEMs that traversed this node to see if this is a redundant sub-MEM
            if (node_matches.count(injected_id)) {
#ifdef debug_multipath_alignment
                cerr << "we need to check if this is a redundant sub MEM, there are previous paths that visited this hit" << endl;
#endif
                
                for (int64_t j : node_matches[injected_id]) {
                    PathNode& match_node = path_nodes[j];
                    
                    if (begin < match_node.begin || end > match_node.end) {
#ifdef debug_multipath_alignment
                        if (begin < match_node.begin) {
                            cerr << "this MEM is earlier in the read than path node " << j << " by " << (match_node.begin - begin) << ", so this is not redundant" << endl;
                        }
                        else if (end > match_node.end) {
                            cerr << "this MEM is later in the read than path node " << j << " by " << (end - match_node.end) << ", so this is not redundant" << endl;
                        }
#endif
                        // the hit does not fall on the same section of the read as the other match, so
                        // it cannot be contained in it
                        continue;
                    }
                    
                    int64_t relative_offset = begin - match_node.begin;
#ifdef debug_multipath_alignment
                    cerr << "the match on node " << j << " has an relative offset of " << relative_offset << " to the this MEM in the read" << endl;
#endif
                    
                    path_t& path = match_node.path;
                    
                    // if this is a partial MEM, we should be able to predict its hit location by traversing the path
                    // of the parent MEM by a distance equal to the relative offset
                    
#ifdef debug_multipath_alignment
                    cerr << "traversing putative parent MEM with path " << debug_string(path) << endl;
#endif
                    
                    int64_t prefix_length = 0;
                    for (size_t k = 0; k < path.mapping_size(); k++) {
                        if (prefix_length > relative_offset) {
#ifdef debug_multipath_alignment
                            cerr << "we have passed where the location would be, breaking out of loop" << endl;
#endif
                            break;
                        }
                        const path_mapping_t& mapping = path.mapping(k);
                        // the length through this mapping
                        int64_t prefix_through_length = prefix_length + mapping_to_length(mapping);
#ifdef debug_multipath_alignment
                        cerr << "after traversing the " << k << "-th step, we have covered a distance of " << prefix_through_length << endl;
#endif
                        if (prefix_through_length > relative_offset) {
                            // we cross the relative offset on this node, so check if the path is in the predicted
                            // position for a redundant sub-MEM
                            id_t node_id_here = mapping.position().node_id();
                            
#ifdef debug_multipath_alignment
                            cerr << "this mapping crosses where we would expect a child to be: " << node_id_here << (project(node_id_here).second ? "-" : "+") << ":" << mapping.position().offset() + relative_offset - prefix_length << endl;
                            cerr << "this MEM is actually at: " << injected_id << (is_rev(hit_pos) ? "-" : "+") << ":" << offset(hit_pos) << endl;
#endif
                            
                            // TODO: shouldn't everything be on the forward strand? i think i could remove the
                            // reverse checking so that only the offset of hit_pos need be communicated to this
                            // function
                            if (injected_id == node_id_here
                                && offset(hit_pos) == mapping.position().offset() + relative_offset - prefix_length
                                && project(node_id_here).second == is_rev(hit_pos)) {
                                
                                // this MEM is redundant with the earlier
                                return true;
                            }
                            

                            
                        }
                        prefix_length = prefix_through_length;
                    }
                }
            }
            return false;
        };
        
        
        
        // we can't filter sub-MEMs on the fly when there are fan-out breaks in this cluster
        // because it's too hard to make sure match nodes get created in descending size order
        bool filter_sub_mems_on_fly;
        if (fanout_breaks != nullptr) {
            if (fanout_breaks->empty()) {
                filter_sub_mems_on_fly = true;
            }
            else {
                bool found_any_breaks = false;
                for (size_t i = 0; i < hits.first.size() && !found_any_breaks; ++i) {
                    found_any_breaks = fanout_breaks->count(hits.first[i].first);
                }
                filter_sub_mems_on_fly = !found_any_breaks;
            }
        }
        else {
            filter_sub_mems_on_fly = true;
        }
        
        size_t num_failed_walks = 0;
        
        // walk the matches and filter out redundant sub-MEMs
        for (int64_t i = 0; i < hits.first.size(); i++) {
            
            pair<const MaximalExactMatch*, pos_t>& hit = hits.first[i];
            
            // the part of the read we're going to match
            string::const_iterator begin = hit.first->begin;
            string::const_iterator end = hit.first->end;
            int64_t mem_length = end - begin;
            // the start of the hit in the original graph
            const pos_t& hit_pos = hit.second;
            
#ifdef debug_multipath_alignment
            cerr << "walking MEM hit " << i << ": " << hit_pos << " " << hit.first->sequence() << endl;
            if (fanout_breaks && fanout_breaks->count(hit.first)) {
                cerr << "fan-out breaks:" << endl;
                for (auto fanout : fanout_breaks->at(hit.first)) {
                    cerr << "\t" << (fanout.first - hit.first->begin) << ": " << *fanout.first << " -> " << fanout.second << endl;
                }
            }
#endif
            bool walked_out_hit = false;
            auto hit_range = injection_trans.equal_range(id(hit_pos));
            for (auto iter = hit_range.first; iter != hit_range.second; iter++) {
                // this graph is unrolled/dagified, so all orientations should match
                if (iter->second.second != is_rev(hit_pos)) {
                    continue;
                }
                
                // an id that corresponds to the original node
                id_t injected_id = iter->second.first;
                
#ifdef debug_multipath_alignment
                cerr << "hit node exists in graph as " << injected_id << endl;
#endif
                if (filter_sub_mems_on_fly) {
                    // don't walk the match of redundant partial hits
                    if (is_redundant(begin, end, hit_pos, injected_id)) {
#ifdef debug_multipath_alignment
                        cerr << "this MEM is identified as a redundant sub-MEM, so we skip it" << endl;
#endif
                        continue;
                    }
                }
                
#ifdef debug_multipath_alignment
                cerr << "performing DFS to walk out match" << endl;
#endif
                
                // TODO: magic constant
                size_t matches_found = 0;
                size_t max_matches_walked = 32;
                
                // stack for DFS, each record contains tuples of
                // (read begin, node offset, next node index, next node handles, fan-out index)
                vector<tuple<string::const_iterator, size_t, size_t, vector<handle_t>, size_t>> stack;
                stack.emplace_back(begin, offset(hit_pos), 0,
                                   vector<handle_t>{graph.get_handle(injected_id)}, 0);
                size_t fanout_size = 0;
                if (fanout_breaks && fanout_breaks->count(hit.first)) {
                    fanout_size = fanout_breaks->at(hit.first).size();
                }
                while (!stack.empty() && matches_found < max_matches_walked) {
                    auto& back = stack.back();
                    if (get<2>(back) == get<3>(back).size()) {
#ifdef debug_multipath_alignment
                        cerr << "traversed all edges out of traversals coming from ";
                        if (stack.size() > 1) {
                            cerr << graph.get_id(get<3>(stack[stack.size() - 2])[get<2>(stack[stack.size() - 2]) - 1]);
                        }
                        else {
                            cerr << "start";
                        }
                        cerr << endl;
#endif
                        stack.pop_back();
                        continue;
                    }
                    
                    handle_t trav = get<3>(back)[get<2>(back)];
                    get<2>(back)++;
                    size_t fanout_idx = get<4>(back);
                    
#ifdef debug_multipath_alignment
                    cerr << "checking node " << graph.get_id(trav) << " at fanout idx " << get<4>(back) << " of " << fanout_size << endl;
#endif
                    
                    string node_seq = graph.get_sequence(trav);
                    size_t node_idx = get<1>(back);
                    string::const_iterator read_iter = get<0>(back);
                    
                    // look for a match along the entire node sequence
                    for (; node_idx < node_seq.size() && read_iter != end; node_idx++, read_iter++) {
                        char read_char;
                        if (fanout_idx < fanout_size
                            && fanout_breaks->at(hit.first)[fanout_idx].first == read_iter) {
                            // we're at the next place where we substituted the read character
                            // for a different one
                            read_char = fanout_breaks->at(hit.first)[fanout_idx].second;
#ifdef debug_multipath_alignment
                            cerr << "\tapplying fanout break to " << read_char << " instead of " << *read_iter << " at index " << (read_iter - begin) << " of MEM" << endl;
#endif
                            ++fanout_idx;
                        }
                        else {
                            // we just want to match the read
                            read_char = *read_iter;
                        }
                        if (node_seq[node_idx] != read_char) {
#ifdef debug_multipath_alignment
                            cerr << "node sequence does not match read" << endl;
#endif
                            break;
                        }
                    }
                    
                    if (read_iter == end) {
                        // finished walking match
#ifdef debug_multipath_alignment
                        cerr << "reached end of read sequence, converting into path node(s) starting at idx " << path_nodes.size() << endl;
#endif
                        assert(fanout_idx == fanout_size);
                        ++matches_found;
                        walked_out_hit = true;
                        
                        path_t path;
                        int64_t path_length = end - begin;
                        int64_t length_remaining = path_length;
                        size_t fanout_idx = 0;
                        int64_t length_until_fanout;
                        if (fanout_size) {
                            length_until_fanout = fanout_breaks->at(hit.first).front().first - begin;
                        }
                        else {
                            length_until_fanout = length_remaining;
                        }
                        auto curr_node_begin = begin;
                        
                        // walk out the match, breaking it at fan-out positions as necessary
                        int64_t offset = get<1>(stack.front());
                        for (size_t j = 0; curr_node_begin < end; ) {
                            auto& search_record = stack[j];
                            handle_t handle = get<3>(search_record)[get<2>(search_record) - 1];
                            int64_t length_on_node = min(int64_t(graph.get_length(handle)) - offset, length_remaining);
                            
                            path_mapping_t* mapping = path.add_mapping();
                            
                            edit_t* edit = mapping->add_edit();
                            
                            // note: the graph is dagified and unrolled, so all hits should be on the forward strand
                            position_t* position = mapping->mutable_position();
                            position->set_node_id(graph.get_id(handle));
                            position->set_offset(offset);
                            
                            if (length_on_node >= length_until_fanout) {
                                // we're either at a fan-out position, or at the end of the path, so
                                // now we want to emit a path node with the path we've just walked out
                                
                                edit->set_from_length(length_until_fanout);
                                edit->set_to_length(length_until_fanout);
                                
                                auto node_end = end - length_remaining + length_until_fanout;
                                
                                if (curr_node_begin < node_end) {
                                    // the node is non-empty
                                    
#ifdef debug_multipath_alignment
                                    cerr << "adding path node for walked match of sequence " << string(curr_node_begin, node_end) << endl;
                                    cerr << debug_string(path) << endl;
                                    cerr << "provenance of path traces to hit " << i << endl;
#endif
                                    
                                    // create a path node
                                    path_nodes.emplace_back();
                                    PathNode& match_node = path_nodes.back();
                                    match_node.path = move(path);
                                    match_node.begin = curr_node_begin;
                                    match_node.end = node_end;
                                    
                                    path_node_provenance.emplace_back(i - num_failed_walks);
                                    
                                    if (filter_sub_mems_on_fly) {
                                        record_node_matches(path_nodes.size() - 1);
                                    }
                                }
#ifdef debug_multipath_alignment
                                else {
                                    cerr << "skipping a walked path that has no sequence" << endl;
                                }
#endif
                                
                                // set up the next path walk
                                path = path_t();
                                
                                // we need to advance past the fan-out character
                                curr_node_begin = node_end + 1;
                                int64_t length_to_advance = length_until_fanout + 1;
                                
                                // walk the path until finding the corresponding position
                                length_remaining -= length_to_advance;
                                while (j < stack.size() &&
                                       offset + length_to_advance >= graph.get_length(get<3>(stack[j])[get<2>(stack[j]) - 1])) {
                                    length_to_advance -= graph.get_length(get<3>(stack[j])[get<2>(stack[j]) - 1]) - offset;
                                    ++j;
                                    offset = 0;
                                }
                                // manipulate the offset so that we start on the correct position
                                offset = length_to_advance;
                                
                                
                                // move the marker for the next fan-out ahead
                                ++fanout_idx;
                                if (fanout_idx < fanout_size) {
                                    length_until_fanout = fanout_breaks->at(hit.first)[fanout_idx].first - curr_node_begin;
                                }
                                else {
                                    length_until_fanout = length_remaining;
                                }
                            }
                            else {
                                edit->set_from_length(length_on_node);
                                edit->set_to_length(length_on_node);
                                
                                // advance to the next stack record
                                ++j;
                                offset = 0;
                                
                                // tick down the length trackers
                                length_remaining -= length_on_node;
                                length_until_fanout -=  length_on_node;

                            }
                        }
                    }
                    else if (node_idx == node_seq.size()) {
                        // matched entire node, move to next node(s)
                        stack.emplace_back(read_iter, 0, 0, vector<handle_t>(), fanout_idx);
                        graph.follow_edges(trav, false, [&](const handle_t& next) {
                            get<3>(stack.back()).push_back(next);
                        });
                    }
                }
            }
            
            // filter out failed walks so they are marked as unclustered
            if (!walked_out_hit) {
                ++num_failed_walks;
            }
            else if (num_failed_walks) {
                hits.first[i - num_failed_walks] = hit;
            }
        }
        
        // clear out the space that we allocated for walks that failed (if any)
        hits.first.resize(hits.first.size() - num_failed_walks);
        
        if (!filter_sub_mems_on_fly) {
            // we weren't removing redundant sub-MEMs as we made them, but now we can do it
            // by sorting the path nodes descending by length
            
            vector<size_t> order(path_nodes.size(), 0);
            for (size_t i = 1; i < path_nodes.size(); ++i) {
                order[i] = i;
            }
            
            stable_sort(order.begin(), order.end(), [&](size_t i, size_t j) {
                const PathNode& a = path_nodes[i];
                const PathNode& b = path_nodes[j];
                return a.end - a.begin > b.end - b.begin;
            });
            vector<size_t> index(order.size());
            for (size_t i = 0; i < order.size(); ++i) {
                index[order[i]] = i;
            }
            for (size_t i = 0; i < index.size(); ++i) {
                std::swap(path_nodes[index[i]], path_nodes[i]);
                std::swap(path_node_provenance[index[i]], path_node_provenance[i]);
                std::swap(index[index[i]], index[i]);
            }
            
            size_t removed_so_far = 0;
            for (size_t i = 0; i < path_nodes.size(); ++i) {
                const position_t& pos = path_nodes[i].path.mapping(0).position();
                // TODO: this seems kinda like overkill, why do we need anything more than the offset in hit_pos?
                auto proj = project(pos.node_id());
                pos_t hit_pos(proj.first, proj.second != pos.is_reverse(), pos.offset());
                
                if (is_redundant(path_nodes[i].begin, path_nodes[i].end, hit_pos, pos.node_id())) {
                    ++removed_so_far;
                }
                else {
                    if (removed_so_far > 0) {
                        path_nodes[i - removed_so_far] = move(path_nodes[i]);
                        path_node_provenance[i - removed_so_far] = path_node_provenance[i];
                    }
                    record_node_matches(i - removed_so_far);
                }
            }
            if (removed_so_far) {
                path_nodes.resize(path_nodes.size() - removed_so_far);
                path_node_provenance.resize(path_nodes.size());
            }
        }
        
        // merge the identical portion of any nodes that overlap exactly
        merge_partially_redundant_match_nodes(node_matches, path_node_provenance);
        
        // if the end of a node is a homopolymer, see if it can be jittered at all
        // and still make a reasonable alignment
        jitter_homopolymer_ends(graph, path_node_provenance, hits, max_branch_trim_length);
    }

    void MultipathAlignmentGraph::merge_partially_redundant_match_nodes(const unordered_map<int64_t, vector<int64_t>>& node_matches,
                                                                        vector<size_t>& path_node_provenance) {
        
#ifdef debug_multipath_alignment
        cerr << "looking for MEMs with partial redundancies to merge" << endl;
#endif
        
        if (path_nodes.size() <= 1) {
            return;
        }
        
        // find the groups that share at least one node
        structures::UnionFind union_find(path_nodes.size(), false);
        for (const auto& match_record : node_matches) {
            for (int64_t i = 1; i < match_record.second.size(); ++i) {
                union_find.union_groups(match_record.second.front(), match_record.second[i]);
            }
        }
        
        // records of (vector of (path node idx, mapping start idx), length)
        vector<pair<vector<pair<size_t, size_t>>, size_t>> identical_segments;
        
        // big loop to identify the identical segments
        for (const auto& overlapping_group : union_find.all_groups()) {
            
            if (overlapping_group.size() == 1) {
                // doesn't overlap with anything
                continue;
            }
            
#ifdef debug_multipath_alignment
            cerr << "looking for merges in overlapping group:" << endl;
            for (auto i : overlapping_group) {
                cerr << "\t" << i << endl;
            }
#endif
            
            // the amount of sequence we've walked for each node in the overlap group
            vector<int64_t> to_length(overlapping_group.size(), 0);
            
            // some helper functions that keep things a bit more succinct later
            
            // go from index within overlapping group to current sequence position
            auto seq_pos = [&](size_t i) {
                return path_nodes[overlapping_group[i]].begin + to_length[i];
            };
            // go from heap record to mapping
            auto next_mapping = [&](const pair<size_t, size_t>& a) {
                return path_nodes[overlapping_group[a.first]].path.mapping(a.second);
            };
            // go from index within overlapping group to path length in mappings
            auto path_length = [&](size_t i) {
                return path_nodes[overlapping_group[i]].path.mapping_size();
            };
            // reverse ordering
            auto heap_cmp = [&](const pair<size_t, size_t>& a, const pair<size_t, size_t>& b) {
                return seq_pos(a.first) > seq_pos(b.first);
            };
            
            // heap of (idx in group, mapping idx)
            vector<pair<size_t, size_t>> heap(overlapping_group.size(), pair<size_t, size_t>(0, 0));
            for (size_t i = 1; i < heap.size(); ++i) {
                heap[i].first = i;
            }
            make_heap(heap.begin(), heap.end(), heap_cmp);
            
            while (heap.size() > 1) {
                // move all of the path nodes whose next mapping occurs
                // at the next sequence position into the unheaped back
                auto heaped_end = heap.end();
                pop_heap(heap.begin(), heaped_end--, heap_cmp);
                while (heap.begin() != heaped_end &&
                       seq_pos(heap.front().first) == seq_pos(heap.back().first)) {
                    pop_heap(heap.begin(), heaped_end--, heap_cmp);
                }
                
#ifdef debug_multipath_alignment
                cerr << "next mappings are at uncentered sequence index " << seq_pos(heap.back().first) - path_nodes.front().begin << ":" << endl;
                for (auto it = heaped_end; it != heap.end(); ++it) {
                    cerr << "\t" << overlapping_group[it->first] << ": " << it->second << " " << debug_string(next_mapping(*it)) << endl;
                }
                cerr << "the other heaped mappings:" << endl;
                for (auto it = heap.begin(); it != heaped_end; ++it) {
                    cerr << "\t" << overlapping_group[it->first] << " (seq index " << seq_pos(it->first) - path_nodes.front().begin << "): " << it->second << " " << debug_string(next_mapping(*it)) << endl;
                }
#endif
                
                // split into groups that have identical next mappings
                vector<vector<size_t>> groups;
                for (size_t i = 0, n = (heap.end() - heaped_end); i < n; ++i) {
                    bool found_match = false;
                    for (auto& group : groups) {
                        if (next_mapping(*(heaped_end + i)) == next_mapping(*(heaped_end + group.front()))) {
                            group.push_back(i);
                            found_match = true;
                            break;
                        }
                    }
                    if (!found_match) {
                        // no matches, becomes its own group
                        groups.emplace_back(1, i);
                    }
                }
                
#ifdef debug_multipath_alignment
                cerr << "partitioned into match groups (by index within unheaped records):" << endl;
                for (auto& group : groups) {
                    cerr << "\t";
                    for (auto i : group) {
                        cerr << " " << i;
                    }
                    cerr << endl;
                }
#endif
                
                for (auto& group : groups) {
                    if (group.size() == 1) {
#ifdef debug_multipath_alignment
                        cerr << "skipping group of size 1 containing " << group.front() << endl;
#endif
                        
                        // no identical mapping to this one
                        auto& heap_rec = *(heaped_end + group.front());
                        to_length[heap_rec.first] += mapping_to_length(next_mapping(heap_rec));
                        ++heap_rec.second;
                    }
                    else {
#ifdef debug_multipath_alignment
                        cerr << "walking match for group with";
                        for (auto i : group) {
                            cerr << " " << i;
                        }
                        cerr << endl;
#endif
                        
                        // at least two mappings are identical
                        identical_segments.emplace_back();
                        auto& segment_rec = identical_segments.back();
                        segment_rec.second = 1;
                        for (auto i : group) {
                            auto& heap_rec = *(heaped_end + i);
                            to_length[heap_rec.first] += mapping_to_length(next_mapping(heap_rec));
                            segment_rec.first.emplace_back(overlapping_group[heap_rec.first], heap_rec.second++);
                        }
                        while (true) {
                            // only continue if all of the segments of this group match
                            bool all_match = true;
                            for (auto i : group) {
                                auto& heap_rec = *(heaped_end + i);
                                if (heap_rec.second == path_length(heap_rec.first)
                                    || next_mapping(heap_rec) != next_mapping(*(heaped_end + group.front()))) {
                                    // we hit the end of a path or found one that doesn't match
                                    all_match = false;
                                    break;
                                }
                            }
                            if (!all_match) {
                                break;
                            }
                            // extend the identical group by 1
                            ++identical_segments.back().second;
                            for (auto i : group) {
                                auto& heap_rec = *(heaped_end + i);
                                to_length[heap_rec.first] += mapping_to_length(next_mapping(heap_rec));
                                ++heap_rec.second;
                            }
                        }
#ifdef debug_multipath_alignment
                        cerr << "walked match of length " << identical_segments.back().second << endl;
#endif
                    }
                }
                
                // now remove any heap records that have reached the end of their path
                for (auto it = heaped_end; it != heap.end();) {
                    if (it->second == path_length(it->first)) {
#ifdef debug_multipath_alignment
                        cerr << "reached end of path node " << overlapping_group[it->first] << endl;
#endif
                        *it = heap.back();
                        heap.pop_back();
                    }
                    else {
#ifdef debug_multipath_alignment
                        cerr << "have not yet exhausted (path node index) " << overlapping_group[it->first] << ", keeping on heap with a new uncentered seq index of " << seq_pos(it->first) - path_nodes.front().begin << endl;
#endif
                        ++it;
                    }
                }
#ifdef debug_multipath_alignment
                cerr << "restoring heap" << endl;
#endif
                
                // and restore the heap
                while (heaped_end != heap.end()) {
                    push_heap(heap.begin(), ++heaped_end, heap_cmp);
                }
            }
        }
        
        if (!identical_segments.empty()) {
            
#ifdef debug_multipath_alignment
            cerr << "found identical segments:" << endl;
            for (auto& segment : identical_segments) {
                for (auto& rec : segment.first) {
                    cerr << "(" << rec.first << ", " << rec.second << ") ";
                }
                cerr << segment.second << endl;
            }
#endif
            
            
            vector<PathNode> merged_path_nodes;
            vector<size_t> merged_provenances;
            merged_path_nodes.reserve(path_nodes.size() + identical_segments.size());
            merged_provenances.reserve(path_nodes.size() + identical_segments.size());
            
            // the index of the last mapping copied over for each path node
            vector<size_t> last_copied(path_nodes.size(), 0);
            
            vector<size_t> copied_to_length(path_nodes.size(), 0);
            
            // function to add a path node for a segment, updating the tracking variables as necessary
            auto add_segment_path_node = [&](size_t orig_idx, size_t seg_begin, size_t length) {
                
                PathNode& orig_path_node = path_nodes[orig_idx];
                int64_t to_length_added = 0;
                if (seg_begin == 0 && length == orig_path_node.path.mapping_size()) {
                    to_length_added = orig_path_node.end - orig_path_node.begin;
                    merged_path_nodes.emplace_back(move(orig_path_node));
                }
                else {
                    merged_path_nodes.emplace_back();
                    PathNode& new_path_node = merged_path_nodes.back();
                    new_path_node.begin = orig_path_node.begin + copied_to_length[orig_idx];
                    new_path_node.path.mutable_mapping()->reserve(length);
                    for (size_t i = seg_begin, n = seg_begin + length; i < n; ++i) {
                        auto& mapping = *orig_path_node.path.mutable_mapping(i);
                        to_length_added += mapping_to_length(mapping);
                        *new_path_node.path.add_mapping() = move(mapping);
                    }
                    new_path_node.end = new_path_node.begin + to_length_added;
                }
                merged_provenances.push_back(path_node_provenance[orig_idx]);
#ifdef debug_multipath_alignment
                cerr << "made merged node from original node " << orig_idx << ", start " << seg_begin << ", len " << length << endl;
                cerr << string(merged_path_nodes.back().begin, merged_path_nodes.back().end) << endl;
                cerr << debug_string(merged_path_nodes.back().path) << endl;
#endif
                return to_length_added;
            };
            
            // by construction, the identical segments are ordered by the read position
            // of their start, so we can iterate in order safely
            for (auto& segment : identical_segments) {
                // copy over any intervening segments
                for (auto& path_node_start : segment.first) {
                    if (last_copied[path_node_start.first] != path_node_start.second) {
                        auto to_length_added = add_segment_path_node(path_node_start.first,
                                                                     last_copied[path_node_start.first],
                                                                     path_node_start.second - last_copied[path_node_start.first]);
                        last_copied[path_node_start.first] = path_node_start.second;
                        copied_to_length[path_node_start.first] += to_length_added;
                    }
                }
                // copy over the merge segments
                auto to_length_added = add_segment_path_node(segment.first.front().first,
                                                             last_copied[segment.first.front().first],
                                                             segment.second);
                for (auto& path_node_start : segment.first) {
                    last_copied[path_node_start.first] = path_node_start.second + segment.second;
                    copied_to_length[path_node_start.first] += to_length_added;
                }
            }
            
            // copy any remaining segments
            for (size_t i = 0; i < path_nodes.size(); ++i) {
                size_t path_length = path_nodes[i].path.mapping_size();
                if (last_copied[i] != path_length) {
                    add_segment_path_node(i, last_copied[i], path_length - last_copied[i]);
                }
            }
            
            path_nodes = move(merged_path_nodes);
            path_node_provenance = move(merged_provenances);
        }
        
#ifdef debug_multipath_alignment
        cerr << "nodes after merging partially redundant segments:" << endl;
        for (size_t i = 0; i < path_nodes.size(); ++i) {
            cerr << i << " (hit " << path_node_provenance[i] << "): " << debug_string(path_nodes[i].path) << " " << string(path_nodes[i].begin, path_nodes[i].end) << endl;
        }
#endif
    }

    void MultipathAlignmentGraph::jitter_homopolymer_ends(const HandleGraph& graph,
                                                          vector<size_t>& path_node_provenance,
                                                          const MultipathMapper::memcluster_t& hits,
                                                          int64_t max_branch_trim_length) {
        
#ifdef debug_multipath_alignment
        cerr << "checking for opportunities to jitter homopolymer anchors" << endl;
#endif
        
        // TODO: magic constants
        static const int64_t min_homopolymer_length = 6;
        // a homopolymer jitter will be accepted if it meets either of these criteria:
        static const int64_t max_jitter_diff = 2;
        static const int64_t min_jitter_length = 5;
        
        size_t num_original_path_nodes = path_nodes.size();
        for (size_t i = 0; i < num_original_path_nodes; ++i) {
            if (path_nodes[i].end - path_nodes[i].begin != hits.first[path_node_provenance[i]].first->length()
                || path_nodes[i].end - path_nodes[i].begin <= min_homopolymer_length) {
                // this node has already been merged with some other node, which gives an indication
                // that alternate exact matches have already handled the local alignment uncertainty
                // or alternatively this node is too short to jitter
                continue;
            }
            
            // TODO: put bookkeeping in place to remove this restriction
            // only try to jitter of one side of a path node at most (both is complicated because
            // we have to sever the source node)
            bool did_jitter = false;
            for (bool left_side : {true, false}) {
                if (did_jitter) {
                    break;
                }
                
                int64_t j_begin, incr;
                if (left_side) {
                    j_begin = 0;
                    incr = 1;
                }
                else {
                    j_begin = path_nodes[i].end - path_nodes[i].begin - 1;
                    incr = -1;
                }
                
                // find the length of homopolymer there is at the end of this match
                // TODO: technically these are homodimers now, but whatever
                int64_t homopolymer_length = 0;
                for (int64_t j = j_begin, n = path_nodes[i].end - path_nodes[i].begin; j >= 0 && j < n; j += incr) {
                    if (*(path_nodes[i].begin + j) == *(path_nodes[i].begin + j_begin + (abs(j - j_begin) % 2) * incr)) {
                        ++homopolymer_length;
                    }
                    else {
                        break;
                    }
                }
                
                if (homopolymer_length >= min_homopolymer_length) {
                    // this is a long enough homopolymer that we'll consider some jitter
                    
#ifdef debug_multipath_alignment
                    cerr << "found homopolymer of length " << homopolymer_length << " on path node " << i << " on left side? " << left_side << endl;
#endif
                                        
                    // walk until the furthest mapping that can be reached by peeling off
                    // the homopolymer
                    int64_t k = left_side ? 0 : path_nodes[i].path.mapping_size() - 1;
                    int64_t length_before = 0;
                    for (; k + incr >= 0 && k + incr < path_nodes[i].path.mapping_size(); k += incr) {
                        int64_t length_thru = length_before + mapping_to_length(path_nodes[i].path.mapping(k));
                        if (length_thru > homopolymer_length) {
                            break;
                        }
                        length_before = length_thru;
                    }
                                        
                    // walk backwards looking for jittered matches
                    handle_t adj_handle = graph.get_handle(path_nodes[i].path.mapping(k).position().node_id(),
                                                           path_nodes[i].path.mapping(k).position().is_reverse());
                    for (; k - incr >= 0 && k - incr < path_nodes[i].path.mapping_size() && !did_jitter; k -= incr) {
                        if (length_before <= max_branch_trim_length) {
                            // we won't bother trying to jitter over regions that will
                            // be caught branch point trimming
                            break;
                        }
#ifdef debug_multipath_alignment
                        cerr << "at mapping " << k << " (" << debug_string(path_nodes[i].path.mapping(k).position()) << ") having already walked " << length_before << endl;
#endif
                        
                        // TODO: contains some redundant code with create_match_nodes
                        
                        handle_t handle = adj_handle;
                        adj_handle = graph.get_handle(path_nodes[i].path.mapping(k - incr).position().node_id(),
                                                      path_nodes[i].path.mapping(k - incr).position().is_reverse());
                        graph.follow_edges(handle, left_side, [&](const handle_t& next) {
                            if (next != adj_handle){
                                // this is an adjacency that the current path doesn't take, so we can
                                // try to jitter down it
                                
#ifdef debug_multipath_alignment
                                cerr << "homopolymer can branch from " << graph.get_id(handle) << " " << graph.get_is_reverse(handle) << " to " << graph.get_id(next) << " " << graph.get_is_reverse(next) << " with " << length_before << " bases to jitter" << endl;
#endif
                                
                                // stack for DFS, each record contains tuples of
                                // (read begin, next node index, next node handles,
                                vector<tuple<string::const_iterator, size_t, vector<handle_t>>> stack;
                                auto riter = left_side ? path_nodes[i].begin + length_before - 1 : path_nodes[i].end - length_before;
                                stack.emplace_back(riter, 0, vector<handle_t>(1, next));
                                while (!stack.empty() && !did_jitter) {
                                    auto& back = stack.back();
                                    if (get<1>(back) == get<2>(back).size()) {
                                        stack.pop_back();
                                        continue;
                                    }
                                    
                                    handle_t trav = get<2>(back)[get<1>(back)];
                                    get<1>(back)++;
                                    
#ifdef debug_multipath_alignment
                                    cerr << "checking for matches node " << graph.get_id(trav) << endl;
#endif
                                    string node_seq = graph.get_sequence(trav);
                                    int64_t node_idx = left_side ? node_seq.size() - 1 : 0;
                                    string::const_iterator read_iter = get<0>(back);
                                    
                                    
                                    // look for a match along the entire node sequence
                                    for (; node_idx >= 0 && node_idx < node_seq.size()
                                         && read_iter >= path_nodes[i].begin && read_iter < path_nodes[i].end; node_idx -= incr, read_iter -= incr) {
#ifdef debug_multipath_alignment
                                        cerr << "comparing MEM index " << (read_iter - path_nodes[i].begin) << " " << *read_iter << " and node index " << node_idx << " " << node_seq[node_idx] << endl;
#endif
                                        if (node_seq[node_idx] != *read_iter) {
#ifdef debug_multipath_alignment
                                            cerr << "found mismatch" << endl;
#endif
                                            
                                            break;
                                        }
                                    }
                                    
                                    if ((node_idx < 0 || node_idx == node_seq.size())
                                        && read_iter >= path_nodes[i].begin && read_iter < path_nodes[i].end) {
                                        // we went off the end of the node without exhausting the MEM
                                        stack.emplace_back(read_iter, 0, vector<handle_t>());
                                        graph.follow_edges(trav, left_side, [&](const handle_t& next) {
                                            get<2>(stack.back()).emplace_back(next);
                                        });
                                    }
                                    else {
                                        int64_t length_diff = left_side ? read_iter - path_nodes[i].begin + 1 : path_nodes[i].end - read_iter;
                                        int64_t jitter_length = length_before - length_diff;
                                        if (length_diff <= max_jitter_diff || jitter_length >= min_jitter_length) {
                                            // we found a jittered anchor with nearly the same length, let's add
                                            // the jittered portion as an alternate anchor
                                            did_jitter = true;
                                            path_nodes.emplace_back();
                                            path_nodes.emplace_back();
                                            auto& split_node = path_nodes[path_nodes.size() - 2];
                                            auto& jitter_node = path_nodes.back();
                                            path_node_provenance.emplace_back(path_node_provenance[i]);
                                            path_node_provenance.emplace_back(path_node_provenance[i]);
                                            if (left_side) {
                                                jitter_node.begin = read_iter + 1;
                                                jitter_node.end = path_nodes[i].begin + length_before;
                                                if (node_idx < (int64_t) node_seq.size() - 1) {
                                                    auto mapping = jitter_node.path.add_mapping();
                                                    auto pos = mapping->mutable_position();
                                                    pos->set_node_id(graph.get_id(trav));
                                                    pos->set_is_reverse(graph.get_is_reverse(trav));
                                                    pos->set_offset(node_idx + 1);
                                                    auto edit = mapping->add_edit();
                                                    edit->set_from_length(node_seq.size() - node_idx - 1);
                                                    edit->set_to_length(edit->from_length());
                                                }
                                                for (int64_t l = stack.size() - 2; l >= 0; --l) {
                                                    handle_t h = get<2>(stack[l])[get<1>(stack[l]) - 1];
                                                    auto mapping = jitter_node.path.add_mapping();
                                                    auto pos = mapping->mutable_position();
                                                    pos->set_node_id(graph.get_id(h));
                                                    pos->set_is_reverse(graph.get_is_reverse(h));
                                                    pos->set_offset(0);
                                                    auto edit = mapping->add_edit();
                                                    edit->set_from_length(graph.get_length(h));
                                                    edit->set_to_length(edit->from_length());
                                                }
                                                // split up the node that we jittered so that it finds the edge to the jitter
                                                split_node.begin = path_nodes[i].begin + length_before;
                                                split_node.end = path_nodes[i].end;
                                                for (int64_t l = k; l < path_nodes[i].path.mapping_size(); ++l) {
                                                    *split_node.path.add_mapping() = move(*path_nodes[i].path.mutable_mapping(l));
                                                }
                                                path_nodes[i].end = split_node.begin;
                                                path_nodes[i].path.mutable_mapping()->resize(k);
                                            }
                                            else {
                                                jitter_node.begin = path_nodes[i].end - length_before;
                                                jitter_node.end = read_iter;
                                                for (int64_t l = 0; l + 1 < stack.size(); ++l) {
                                                    handle_t h = get<2>(stack[l])[get<1>(stack[l]) - 1];
                                                    auto mapping = jitter_node.path.add_mapping();
                                                    auto pos = mapping->mutable_position();
                                                    pos->set_node_id(graph.get_id(h));
                                                    pos->set_is_reverse(graph.get_is_reverse(h));
                                                    pos->set_offset(0);
                                                    auto edit = mapping->add_edit();
                                                    edit->set_from_length(graph.get_length(h));
                                                    edit->set_to_length(edit->from_length());
                                                }
                                                if (node_idx > 0) {
                                                    auto mapping = jitter_node.path.add_mapping();
                                                    auto pos = mapping->mutable_position();
                                                    pos->set_node_id(graph.get_id(trav));
                                                    pos->set_is_reverse(graph.get_is_reverse(trav));
                                                    pos->set_offset(0);
                                                    auto edit = mapping->add_edit();
                                                    edit->set_from_length(node_idx);
                                                    edit->set_to_length(edit->from_length());
                                                }
                                                // split up the node that we jittered so that it finds the edge to the jitter
                                                split_node.begin = path_nodes[i].end - length_before;
                                                split_node.end = path_nodes[i].end;
                                                for (int64_t l = k + 1; l < path_nodes[i].path.mapping_size(); ++l) {
                                                    *split_node.path.add_mapping() = move(*path_nodes[i].path.mutable_mapping(l));
                                                }
                                                path_nodes[i].end = split_node.begin;
                                                path_nodes[i].path.mutable_mapping()->resize(k + 1);
                                            }
                                            
#ifdef debug_multipath_alignment
                                            cerr << "jitter difference of " << length_diff << " was small enough or length of " << jitter_length << " was large enough to make a new jitter anchor: " << endl;
                                            cerr << string(jitter_node.begin, jitter_node.end) << " (" << (jitter_node.begin - path_nodes[i].begin) << ":" << (jitter_node.end - path_nodes[i].begin) << ") " << debug_string(jitter_node.path) << endl;
                                            cerr << "new split nodes in path node " << i << endl;
                                            cerr << string(path_nodes[i].begin, path_nodes[i].end) << " (" << (path_nodes[i].begin - path_nodes[i].begin) << ":" << (path_nodes[i].end - path_nodes[i].begin) << ") " << debug_string(path_nodes[i].path) << endl;
                                            cerr << string(split_node.begin, split_node.end) << " (" << (split_node.begin - path_nodes[i].begin) << ":" << (split_node.end - path_nodes[i].begin) << ") " << debug_string(split_node.path) << endl;
#endif
                                        }
                                    }
                                }
                            }
                            return !did_jitter;
                        });
                        
                        if (!did_jitter) {
                            length_before -= mapping_to_length(path_nodes[i].path.mapping(k - incr));
                        }
                    }
                }
            }
        }
    }
    
    void MultipathAlignmentGraph::collapse_order_length_runs(const HandleGraph& graph, gcsa::GCSA* gcsa,
                                                             vector<size_t>& path_node_provenance) {
        
#ifdef debug_multipath_alignment
        cerr << "looking for runs of order length MEMs to collapse with gcsa order "  << gcsa->order() << endl;
#endif
        
        vector<vector<size_t>> merge_groups;
        
        size_t num_order_length_mems = 0;
        for (size_t i = 0; i < path_nodes.size(); i++) {
            
            PathNode& match_node = path_nodes.at(i);
            
            if (match_node.end - match_node.begin < gcsa->order()) {
                // we have passed all of the order length MEMs, bail out of loop
#ifdef debug_multipath_alignment
                cerr << "found " << i << " order length MEMs" << endl;
#endif
                
                num_order_length_mems = i;
                break;
            }
        }
        
        // find the order of the order-length MEMs lexicographically along the read
        vector<size_t> order(num_order_length_mems, 0);
        for (size_t i = 1; i < order.size(); i++) {
            order[i] = i;
        }
        sort(order.begin(), order.end(), [&](size_t i, size_t j) {
            return path_nodes.at(i).begin < path_nodes.at(j).begin || (path_nodes.at(i).begin == path_nodes.at(j).begin &&
                                                                      path_nodes.at(i).end < path_nodes.at(j).end);
        });
        
        for (size_t i : order) {
            
            PathNode& match_node = path_nodes[i];
            
#ifdef debug_multipath_alignment
            cerr << "checking if MEM " << i << " can be an extension: " << endl;
            cerr << "\t" << string(match_node.begin, match_node.end) << "\t" << debug_string(match_node.path) << endl;
#endif
            
            // try to find any run of MEMs that could be merged with this MEM
            bool found_merge_group = false;
            for (size_t j = 0; j < merge_groups.size(); j++) {
                
                vector<size_t>& merge_group = merge_groups[j];
                
                // because of the sort order, the last node in this run should overlap the current MEM
                // if any of them can
                PathNode& last_run_node = path_nodes[merge_group[merge_group.size() - 1]];
                
#ifdef debug_multipath_alignment
                cerr << "checking against extending MEM " << merge_group[merge_group.size() - 1] << " in merge group " << j << ":" << endl;
                cerr << "\t";
                for (auto iter = last_run_node.begin; iter != last_run_node.end; iter++) {
                    cerr << *iter;
                }
                cerr << endl;
                cerr << "\t" << debug_string(last_run_node.path) << endl;
#endif
                
                // do they overhang an amount on the read that indicates they overlap and could be merged?
                int64_t overhang = last_run_node.end - match_node.begin;
                if (last_run_node.begin < match_node.begin && match_node.end > last_run_node.end && overhang >= 0) {
                    
#ifdef debug_multipath_alignment
                    cerr << "MEMs overlap on the read, checking for consistency with overhang on the path" << endl;
#endif
                    
                    // get the initial position of the node further to the right
                    pos_t match_node_initial_pos = make_pos_t(match_node.path.mapping(0).position());
                    
                    // get the position at the overhang back from the end of the node further to the left
                    int64_t remaining = last_run_node.end - last_run_node.begin;
                    pos_t last_run_node_internal_pos;
                    for (size_t k = 0; k < last_run_node.path.mapping_size(); k++) {
                        
                        int64_t mapping_length = mapping_from_length(last_run_node.path.mapping(k));
                        
                        if (remaining - mapping_length < overhang) {
                            // we will cross the position that should line up with the initial position on this mapping
                            
                            const position_t& overhang_position = last_run_node.path.mapping(k).position();
                            
                            get_id(last_run_node_internal_pos) = overhang_position.node_id();
                            get_is_rev(last_run_node_internal_pos) = overhang_position.is_reverse();
                            get_offset(last_run_node_internal_pos) = overhang_position.offset() + remaining - overhang;
                            
                            break;
                        }
                        
                        remaining -= mapping_length;
                    }
                    
                    // get the final position of the node further to the left
                    const path_mapping_t& final_mapping = last_run_node.path.mapping(last_run_node.path.mapping_size() - 1);
                    const position_t& final_mapping_position = final_mapping.position();
                    pos_t last_run_node_final_pos = make_pos_t(final_mapping_position.node_id(),
                                                               final_mapping_position.is_reverse(),
                                                               final_mapping_position.offset() + mapping_from_length(final_mapping));
                    
                    // get the position at the overhang into the node further to the right
                    remaining = match_node.end - match_node.begin;
                    pos_t match_node_internal_pos;
                    for (int64_t k = match_node.path.mapping_size() - 1; k >= 0; k--) {
                        
                        int64_t mapping_length = mapping_from_length(match_node.path.mapping(k));
                        remaining -= mapping_length;
                        
                        if (remaining < overhang) {
                            // we will cross the position that should line up with the initial position on this mapping
                            
                            const position_t& overhang_position = match_node.path.mapping(k).position();
                            
                            get_id(match_node_internal_pos) = overhang_position.node_id();
                            get_is_rev(match_node_internal_pos) = overhang_position.is_reverse();
                            get_offset(match_node_internal_pos) = overhang_position.offset() + overhang - remaining;
                            
                            break;
                        }
                    }
                    
#ifdef debug_multipath_alignment
                    cerr << "overhang segments on extension node " << match_node_initial_pos << ", " << match_node_internal_pos << " and extending node " << last_run_node_internal_pos << ", " << last_run_node_final_pos << endl;
#endif
                    
                    // do the positions match up as we would expect if these are actually part of the same match?
                    if (match_node_initial_pos == last_run_node_internal_pos && last_run_node_final_pos == match_node_internal_pos) {
                        
#ifdef debug_multipath_alignment
                        cerr << "found matching overhang" << endl;
#endif
                        
                        // add the match node to this merge group
                        merge_group.push_back(i);
                        found_merge_group = true;
                        
                        break;
                    }
                    else if (overhang == 0 && offset(match_node_initial_pos) == 0) {
                        // it could still be that these are two end-to-end matches that got assigned to the beginning
                        // and end of two nodes connected by an edge
                        
                        if (offset(last_run_node_final_pos) == graph.get_length(graph.get_handle(final_mapping_position.node_id()))) {
                            if (graph.has_edge(graph.get_handle(id(last_run_node_final_pos), is_rev(last_run_node_final_pos)),
                                               graph.get_handle(id(match_node_initial_pos), is_rev(match_node_initial_pos)))) {
#ifdef debug_multipath_alignment
                                cerr << "found end to end connection over an edge" << endl;
#endif
                                
                                // add the match node to this merge group
                                merge_group.push_back(i);
                                found_merge_group = true;
                                
                                break;
                            }
                        }
                    }
                }
            }
            
            if (!found_merge_group) {
                // make a new merge group consisting of only this MEM
                merge_groups.emplace_back(1, i);
            }
        }
        
#ifdef debug_multipath_alignment
        cerr << "merge groups among order length MEMs:" << endl;
        for (auto& group : merge_groups) {
            cerr << "\t";
            for (auto i : group) {
                cerr << i << " ";
            }
            cerr << endl;
        }
#endif
        
        if (merge_groups.size() != num_order_length_mems) {
            // we found at least one merge to do, now we need to actually do the merges
            
            unordered_set<size_t> to_remove;
            
            for (const vector<size_t>& merge_group : merge_groups) {
                
                // merge the paths into the first node in the group (arbitrarily)
                PathNode& merge_into_node = path_nodes[merge_group[0]];
                
                for (size_t i = 1; i < merge_group.size(); i++) {
                    
                    // mark the node we're merging from for removal
                    to_remove.insert(merge_group[i]);
                    
                    PathNode& merge_from_node = path_nodes.at(merge_group[i]);
                    
#ifdef debug_multipath_alignment
                    cerr << "merging into node " << merge_group[0] << " path " << debug_string(merge_into_node.path) << endl;
                    cerr << "from node " << merge_group[i] << " path " << debug_string(merge_from_node.path) << endl;
#endif
                    
                    // walk backwards until we find the first mapping to add
                    int64_t to_add_length = merge_from_node.end - merge_into_node.end;
                    int64_t remaining = to_add_length;
                    int64_t first_mapping_to_add_idx = 0;
                    for (int64_t j = merge_from_node.path.mapping_size() - 1; j >= 0; j--) {
                        remaining -= mapping_from_length(merge_from_node.path.mapping(j));
                        if (remaining <= 0) {
                            first_mapping_to_add_idx = j;
                            break;
                        }
                    }
                    
#ifdef debug_multipath_alignment
                    cerr << "the first mapping to merge is at index " << first_mapping_to_add_idx << " with " << -remaining << " remaining" << endl;
#endif
                    
                    // handle the first mapping we add as a special case
                    const path_mapping_t& first_mapping_to_add = merge_from_node.path.mapping(first_mapping_to_add_idx);
                    path_mapping_t* final_merging_mapping = merge_into_node.path.mutable_mapping(merge_into_node.path.mapping_size() - 1);
                    if (final_merging_mapping->position().node_id() == first_mapping_to_add.position().node_id() &&
                        final_merging_mapping->position().is_reverse() == first_mapping_to_add.position().is_reverse() &&
                        first_mapping_to_add.position().offset() - final_merging_mapping->position().offset() - remaining == mapping_from_length(*final_merging_mapping)) {
                        
                        // the mappings are on the same node, so they can be combined
                        int64_t mapping_to_add_length = mapping_from_length(first_mapping_to_add) + remaining;
                        edit_t* final_edit = final_merging_mapping->mutable_edit(final_merging_mapping->edit_size() - 1);
                        final_edit->set_from_length(final_edit->from_length() + mapping_to_add_length);
                        final_edit->set_to_length(final_edit->to_length() + mapping_to_add_length);
                        
#ifdef debug_multipath_alignment
                        cerr << "merged mapping is " << debug_string(*final_merging_mapping) << endl;
#endif
                    }
                    else {
                        // we need to add this as a new mapping
                        path_mapping_t* new_mapping = merge_into_node.path.add_mapping();
                        *new_mapping = first_mapping_to_add;
                        
#ifdef debug_multipath_alignment
                        cerr << "new adjacent mapping is " << debug_string(*new_mapping) << endl;
#endif
                    }
                    
                    // add the remaining mappings as new mappings
                    for (size_t j = first_mapping_to_add_idx + 1; j < merge_from_node.path.mapping_size(); j++) {
                        path_mapping_t* new_mapping = merge_into_node.path.add_mapping();
                        *new_mapping = merge_from_node.path.mapping(j);
                        
#ifdef debug_multipath_alignment
                        cerr << "new transfer mapping is " << debug_string(*new_mapping) << endl;
#endif
                    }
                    
                    // merge the substrings on the node
                    merge_into_node.end = merge_from_node.end;
                    
#ifdef debug_multipath_alignment
                    cerr << "merged path is " << debug_string(merge_into_node.path) << endl;
                    cerr << "merged substring is ";
                    for (auto iter = merge_into_node.begin; iter != merge_into_node.end; iter++) {
                        cerr << *iter;
                    }
                    cerr << endl;
#endif
                }
            }
            
            // remove all of the nodes we merged into other nodes
            size_t removed_so_far = 0;
            for (size_t i = 0; i < path_nodes.size(); i++) {
                if (to_remove.count(i)) {
                    removed_so_far++;
                }
                else if (removed_so_far > 0) {
                    path_nodes[i - removed_so_far] = move(path_nodes[i]);
                    path_node_provenance[i - removed_so_far] = path_node_provenance[i];
#ifdef debug_multipath_alignment
                    cerr << "moving path node " << i << " into index " << i - removed_so_far << endl;
#endif
                }
#ifdef debug_multipath_alignment
                else {
                    cerr << "moving path node " << i << " into index " << i << endl;
                }
#endif
            }
            
            path_nodes.resize(path_nodes.size() - to_remove.size());
            path_node_provenance.resize(path_nodes.size());
        }
#ifdef debug_multipath_alignment
        cerr << "done merging MEMs" << endl;
#endif
    }
    
    void MultipathAlignmentGraph::trim_to_branch_points(const HandleGraph* graph, size_t max_trim_length) {
        
        assert(!has_reachability_edges);
        
        for (PathNode& path_node : path_nodes) {
            
#ifdef debug_multipath_alignment
            cerr << "trimming to branch points within " << max_trim_length << " of ends in path " << debug_string(path_node.path) << endl;
#endif
            
            // find the mapping where we are first pass the maximum trim length coming inward
            // from the left
            int64_t from_length = 0;
            int64_t to_length = 0;
            int64_t prefix_idx = 0;
            for (; prefix_idx < path_node.path.mapping_size()
                 && from_length <= max_trim_length
                 && to_length <= max_trim_length; prefix_idx++) {
                
                from_length += mapping_from_length(path_node.path.mapping(prefix_idx));
                to_length += mapping_to_length(path_node.path.mapping(prefix_idx));
            }
            
            // walk backwards to see if we passed any leftward branch points
            for (prefix_idx--; prefix_idx > 0; prefix_idx--) {
                const position_t& pos = path_node.path.mapping(prefix_idx).position();
                if (graph->get_degree(graph->get_handle(pos.node_id(), pos.is_reverse()), true) > 1) {
                    // this is the inward most branch point within the trim length
                    
#ifdef debug_multipath_alignment
                    cerr << "found leftward branch point at " << prefix_idx << endl;
#endif
                    
                    break;
                }
            }
            
            // find the mapping where we are first pass the maximum trim length coming inward
            // from the right
            from_length = 0;
            to_length = 0;
            int64_t suffix_idx = path_node.path.mapping_size() - 1;
            for (; suffix_idx >= 0
                 && from_length <= max_trim_length
                 && to_length <= max_trim_length; suffix_idx--) {
                
                from_length += mapping_from_length(path_node.path.mapping(suffix_idx));
                to_length += mapping_to_length(path_node.path.mapping(suffix_idx));
            }
            
            // walk forward to see if we passed any rightwards branch points
            for (suffix_idx++; suffix_idx + 1 < path_node.path.mapping_size(); suffix_idx++) {
                const position_t& pos = path_node.path.mapping(suffix_idx).position();
                if (graph->get_degree(graph->get_handle(pos.node_id(), pos.is_reverse()), false) > 1) {
                    // this is the inward most branch point within the trim length
                    
#ifdef debug_multipath_alignment
                    cerr << "found right branch point at " << suffix_idx << endl;
#endif
                    break;
                }
            }
            
            // the prefix/suffix idx now indicate the place after which we can trim
            
            if (prefix_idx > suffix_idx) {
                // this is a weird case, we seem to want to trim the whole anchor, which suggests
                // that maybe the trim length was chosen to be too long. there are other explanations
                // but we're just going to ignore the trimming for now
                // TODO: is this the right approach?
#ifdef debug_multipath_alignment
                cerr << "seem to want to trim entire path, skipping" << endl;
#endif
                continue;
            }
            
            if (prefix_idx > 0 || suffix_idx + 1 < path_node.path.mapping_size()) {
                
                // compute the amount of read we're trimming off from each end
                int64_t trimmed_prefix_to_length = 0;
                int64_t trimmed_suffix_to_length = 0;
                for (int64_t i = 0; i < prefix_idx; i++) {
                    trimmed_prefix_to_length += mapping_to_length(path_node.path.mapping(i));
                }
                for (int64_t i = suffix_idx + 1; i < path_node.path.mapping_size(); i++) {
                    trimmed_suffix_to_length += mapping_to_length(path_node.path.mapping(i));
                }
                
                // replace the path with the portion that we didn't trim
                path_t new_path;
                for (int64_t i = prefix_idx; i <= suffix_idx; i++) {
                    path_mapping_t* mapping = new_path.add_mapping();
                    *mapping = path_node.path.mapping(i);
                }
                path_node.path = move(new_path);
#ifdef debug_multipath_alignment
                cerr << "trimmed path: " << debug_string(path_node.path) << endl;
#endif
                // update the read interval
                path_node.begin += trimmed_prefix_to_length;
                path_node.end -= trimmed_suffix_to_length;
            }
        }
    }

    vector<pair<size_t, size_t>> MultipathAlignmentGraph::get_cut_segments(path_t& path,
                                                                           SnarlManager* cutting_snarls,
                                                                           MinimumDistanceIndex* dist_index,
                                                                           const function<pair<id_t, bool>(id_t)>& project,
                                                                           int64_t max_snarl_cut_size) const {
        
        // this list holds the beginning of the current segment at each depth in the snarl hierarchy
        // as we traverse the exact match, the beginning is recorded in both sequence distance and node index
        list<pair<size_t, size_t>> level_segment_begin;
        level_segment_begin.emplace_back(0, 0);
        
        // we record which segments we are going to cut out of the match here
        vector<pair<size_t, size_t>> cut_segments;
        
        auto curr_level = level_segment_begin.begin();
        size_t prefix_length = 0;
        for (size_t j = 0, last = path.mapping_size() - 1; j <= last; j++) {
            const position_t& position = path.mapping(j).position();
            const auto& projection = project(position.node_id());
            id_t projected_id = projection.first;
            bool projected_rev = (projection.second != position.is_reverse());
            
            if (j > 0) {
                // we have entered this node on this iteration
                if (into_cutting_snarl(projected_id, !projected_rev, cutting_snarls, dist_index)) {
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
                        curr_level = level_segment_begin.insert(level_segment_begin.end(), pair<size_t, size_t>(0, 0));
                    }
                }
            }
            
            // cross to the other side of the node
            prefix_length += mapping_from_length(path.mapping(j));
            
            if (j < last) {
                // we are going to leave this node next iteration
                if (into_cutting_snarl(projected_id, projected_rev, cutting_snarls, dist_index)) {
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
        if ((prefix_length - curr_level->first <= max_snarl_cut_size) && curr_level != last) {
            cut_segments.emplace_back(curr_level->second, path.mapping_size());
        }
        
#ifdef debug_multipath_alignment
        cerr << "found " << cut_segments.size() << " cut segments:" << endl;
        for (auto seg : cut_segments) {
            cerr << "\t" << seg.first << ":" << seg.second << endl;
        }
#endif
        
        
        return cut_segments;
    }
    
    void MultipathAlignmentGraph::resect_snarls_from_paths(SnarlManager* cutting_snarls,
                                                           MinimumDistanceIndex* dist_index,
                                                           const function<pair<id_t, bool>(id_t)>& project,
                                                           int64_t max_snarl_cut_size) {
#ifdef debug_multipath_alignment
        cerr << "cutting with snarls" << endl;
#endif
        
        size_t num_original_path_nodes = path_nodes.size();
        
        // we'll need to keep track of which nodes we trim the front off of to update edge lengths later
        vector<size_t> trimmed_prefix_length(path_nodes.size());
        bool trimmed_any_prefix = false;
        
        for (size_t i = 0; i < num_original_path_nodes; i++) {
            
            // first compute the segments we want to cut out
            
            PathNode* path_node = &path_nodes.at(i);
            path_t* path = &path_node->path;
            
#ifdef debug_multipath_alignment
            cerr << "cutting node at index " << i << " with path " << debug_string(*path) << endl;
#endif
            
            auto cut_segments = get_cut_segments(*path, cutting_snarls, dist_index, project, max_snarl_cut_size);
            
            // did we cut out any segments?
            if (!cut_segments.empty()) {
                
                vector<pair<size_t, size_t>> keep_segments;
                
                // we may have decided to cut the segments of both a parent and child snarl, so now we
                // collapse the list of intervals, which is sorted on the end index by construction
                //
                // snarl nesting properties guarantee that there will be at least one node between any
                // cut segments that are not nested, so we don't need to deal with the case where the
                // segments are partially overlapping (i.e. it's a bit easier than the general interval
                // intersection problem)
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
                
                // the keep segments are now stored last-to-first, let's reverse them to their more natural ordering
                reverse(keep_segments.begin(), keep_segments.end());
                
                // record the data stored on the original path node
                path_t original_path = move(*path);
                string::const_iterator original_begin = path_node->begin;
                string::const_iterator original_end = path_node->end;
                vector<pair<size_t, size_t>> forward_edges = move(path_node->edges);
                
                // and reinitialize the node
                path_node->edges.clear();
                path->clear_mapping();
                
                size_t prefix_from_length = 0;
                size_t prefix_to_length = 0;
                size_t prefix_idx = 0;
                while (prefix_idx < keep_segments.front().first) {
                    prefix_from_length += mapping_from_length(original_path.mapping(prefix_idx));
                    prefix_to_length += mapping_to_length(original_path.mapping(prefix_idx));
                    prefix_idx++;
                }
                
                // keep track whether we trimmed a prefix off the left side of the priginal path
                trimmed_prefix_length[i] = prefix_from_length;
                trimmed_any_prefix = trimmed_any_prefix || (prefix_from_length > 0);
                
#ifdef debug_multipath_alignment
                cerr << "making path for initial keep segment " << keep_segments.front().first << ":" << keep_segments.front().second << " at idx " << i << endl;
#endif
                
                // place the first keep segment into the original node
                path_node->begin = original_begin + prefix_to_length;
                for (int32_t rank = 1; prefix_idx < keep_segments.front().second; prefix_idx++, rank++) {
                    path_mapping_t* mapping = path->add_mapping();
                    *mapping = original_path.mapping(prefix_idx);
                    prefix_from_length += mapping_from_length(*mapping);
                    prefix_to_length += mapping_to_length(*mapping);
                }
                path_node->end = original_begin + prefix_to_length;
                
                
#ifdef debug_multipath_alignment
                cerr << "new cut path: " << debug_string(path_node->path) << endl;
#endif
                
                // keep track of the index in the node vector of the previous segment
                size_t prev_segment_idx = i;
                
                for (size_t j = 1; j < keep_segments.size(); j++) {
                    
                    auto& keep_segment = keep_segments[j];
                    
#ifdef debug_multipath_alignment
                    cerr << "making path for next keep segment " << keep_segments[j].first << ":" << keep_segments[j].second << " at idx " << path_nodes.size() << endl;
#endif
                    
                    // record the start of the intersegment section of the read
                    size_t intersegment_start = prefix_from_length;
                    
                    // advance to the next keep segment
                    while (prefix_idx < keep_segment.first) {
                        prefix_from_length += mapping_from_length(original_path.mapping(prefix_idx));
                        prefix_to_length += mapping_to_length(original_path.mapping(prefix_idx));
                        prefix_idx++;
                    }
                    
                    // create a new node for this keep segment
                    path_nodes.emplace_back();
                    PathNode& cut_node = path_nodes.back();
                    path_t& cut_path = cut_node.path;
                    
                    // add a connecting edge from the last keep segment
                    path_nodes.at(prev_segment_idx).edges.emplace_back(path_nodes.size() - 1, prefix_from_length - intersegment_start);
                    
                    // transfer over the path and the read interval
                    cut_node.begin = original_begin + prefix_to_length;
                    for (int32_t rank = 1; prefix_idx < keep_segment.second; prefix_idx++, rank++) {
                        path_mapping_t* mapping = cut_path.add_mapping();
                        *mapping = original_path.mapping(prefix_idx);
                        prefix_from_length += mapping_from_length(*mapping);
                        prefix_to_length += mapping_to_length(*mapping);
                    }
                    cut_node.end = original_begin + prefix_to_length;
                    
#ifdef debug_multipath_alignment
                    cerr << "new cut path: " << debug_string(cut_path) << endl;
#endif
                    
                    prev_segment_idx = path_nodes.size() - 1;
                }
                
                // move the edges from the original node onto the last keep segment
                path_nodes.at(prev_segment_idx).edges = move(forward_edges);
                
                // add the length of the trimmed portion of the path to the edge length
                size_t trimmed_suffix_length = path_from_length(original_path) - prefix_from_length;
                if (trimmed_suffix_length) {
                    for (pair<size_t, size_t>& edge : path_nodes.at(prev_segment_idx).edges) {
                        edge.second += trimmed_suffix_length;
                    }
                }
            }
        }
        
        if (trimmed_any_prefix) {
            // we need to add the length of the prefixes we trimmed off to the length of edges
            for (PathNode& path_node : path_nodes) {
                for (pair<size_t, size_t>& edge : path_node.edges) {
                    if (edge.first < trimmed_prefix_length.size()) {
                        edge.second += trimmed_prefix_length[edge.first];
                    }
                }
            }
        }
    }
   
    void MultipathAlignmentGraph::synthesize_tail_anchors(const Alignment& alignment, const HandleGraph& align_graph, const GSSWAligner* aligner,
                                                          size_t min_anchor_length, size_t max_alt_alns, bool dynamic_alt_alns, size_t max_gap,
                                                          double pessimistic_tail_gap_multiplier) {
    
        // Align the tails, not collecting a set of source subpaths.
        // TODO: factor of 1/2 is arbitray, but i do think it should be fewer than the max
        auto tail_alignments = align_tails(alignment, align_graph, aligner, max<size_t>(1, max_alt_alns / 2),
                                           dynamic_alt_alns, max_gap, pessimistic_tail_gap_multiplier, max_alt_alns, nullptr);
        
        
        for (bool handling_right_tail : {false, true}) {
            // For each tail we are processing
        
            for (auto kv : tail_alignments[handling_right_tail]) {
                // For each node that has alignments off in that direction
                
                // Grab the PathNode we are attached to on the right or left side
                auto& attached_path_node_index = kv.first;
                // And all the alignments off of there
                auto& alns = kv.second;
                
                // only make anchors from alignments that have score equal to the optimal
                for (size_t i = 1; i < alns.size(); ++i) {
                    if (alns[i].score() < alns.front().score()) {
                        alns.resize(i);
                        break;
                    }
                }
                
#ifdef debug_multipath_alignment
                cerr << "Handling " << (handling_right_tail ? "right" : "left") << " tail off of PathNode "
                    << attached_path_node_index << " with path " << debug_string(path_nodes.at(attached_path_node_index).path) << endl;
#endif
                
                for (auto& aln : alns) {
                    
#ifdef debug_multipath_alignment
                    cerr << "Tail alignment: " << pb2json(aln) << endl;
#endif
                    
                    auto seq_begin = alignment.sequence().begin() + (handling_right_tail ? (alignment.sequence().size() - aln.sequence().size()) : 0);
                    
                    // how far have we traveled in the graph since the start of the alignment
                    size_t cumul_from_length = 0;
                    // how far have we traveled along the read since the start of the alignment
                    size_t cumul_to_length = 0;
                    
                    // we will keep track of where we are on the current node
                    size_t offset_on_curr_node = numeric_limits<size_t>::max();
                    
                    // when we find a match, we will keep track of
                    size_t curr_match_length = 0;
                    size_t match_start_mapping_idx = numeric_limits<size_t>::max();
                    size_t match_start_edit_idx = numeric_limits<size_t>::max();
                    size_t curr_match_start_offset = numeric_limits<size_t>::max();
                    
                    // if it's the right tail, we know the previous index, otherwise there is no previous index
                    size_t prev_anchor_path_node = handling_right_tail ? attached_path_node_index : numeric_limits<size_t>::max();
                    size_t prev_anchor_final_from_length = 0;
                    
                    const Path& path = aln.path();
                    
                    // a function that will create a new match node based on these trackers
                    auto create_synthetic_anchor_node = [&](const size_t& i, const size_t& j) {
                        
                        path_nodes.emplace_back();
                        PathNode& synth_path_node = path_nodes.back();
                                                
                        // copy the first mapping, paying attention to the initial position
                        path_mapping_t* new_mapping = synth_path_node.path.add_mapping();
                        position_t* new_position = new_mapping->mutable_position();
                        new_position->set_node_id(path.mapping(match_start_mapping_idx).position().node_id());
                        new_position->set_is_reverse(path.mapping(match_start_mapping_idx).position().is_reverse());
                        new_position->set_offset(curr_match_start_offset);
                        
                        // we should only be able to copy one edit over from this mapping, either because the next one
                        // is a mismatch or because it's on the next node's mapping
                        from_proto_edit(path.mapping(match_start_mapping_idx).edit(match_start_edit_idx),
                                        *new_mapping->add_edit());
                                                
                        // copy any whole mappings from the middle of the anchor path
                        for (size_t copy_i = match_start_mapping_idx + 1; copy_i < i; copy_i++) {
                            assert(path.mapping(copy_i).edit_size() == 1);
                            from_proto_mapping(path.mapping(copy_i), *synth_path_node.path.add_mapping());
                        }
                        
                        // on the final mapping we don't need to pay special attention to the initial position, but
                        // we can't copy over the whole mapping since there might be more edits after the match
                        if (i > match_start_mapping_idx && j > 0) {
                            // This condition is broken because N matches get split into separate edits
                            assert(j == 1);
                            new_mapping = synth_path_node.path.add_mapping();
                            position_t* pos = new_mapping->mutable_position();
                            const Position& pos_from = path.mapping(i).position();
                            pos->set_node_id(pos_from.node_id());
                            pos->set_offset(pos_from.offset());
                            pos->set_is_reverse(pos_from.is_reverse());
                            from_proto_edit(path.mapping(i).edit(0), *new_mapping->add_edit());
                        }
                        
                        synth_path_node.end = seq_begin + cumul_to_length;
                        synth_path_node.begin = synth_path_node.end - curr_match_length;
                        
#ifdef debug_multipath_alignment
                        cerr << "\tyielded anchor with path " << debug_string(synth_path_node.path) << " and seq ";
                        for (auto it = synth_path_node.begin; it != synth_path_node.end; ++it) {
                            cerr << *it;
                        }
                        cerr << endl;
                        
#endif
                        // make an edge from the previous synthetic anchor (if it exists)
                        if (prev_anchor_path_node != numeric_limits<size_t>::max()) {
                            path_nodes[prev_anchor_path_node].edges.emplace_back(path_nodes.size() - 1,
                                                                                 cumul_from_length - curr_match_length - prev_anchor_final_from_length);
#ifdef debug_multipath_alignment
                            cerr << "\talso making edge from path node at " << prev_anchor_path_node << " with length " << cumul_from_length - curr_match_length - prev_anchor_final_from_length << endl;
                            
#endif
                        }
                        
                        // mark this anchor as the new anchor
                        prev_anchor_path_node = path_nodes.size() - 1;
                        prev_anchor_final_from_length = cumul_from_length;
                    };
                    
                    // iterate over the path, updating the tracking variables as we go
                    for (size_t i = 0; i < path.mapping_size(); i++) {
                        const Mapping& mapping = path.mapping(i);
                        
                        // new node, new offset
                        offset_on_curr_node = mapping.position().offset();
                        
                        for (size_t j = 0; j < mapping.edit_size(); j++) {
                            
                            const Edit& edit = mapping.edit(j);
                            if (edit.from_length() != edit.to_length() || !edit.sequence().empty()) {
                                // we've found a non-match edit
                                
                                if (curr_match_length >= min_anchor_length) {
                                    // we are coming out of a match that is long enough for us to make a new anchor
                                    create_synthetic_anchor_node(i, j);
                                }
                                
                                // mark the trackers that indicate that we are not in a match
                                curr_match_start_offset = numeric_limits<size_t>::max();
                                curr_match_length = 0;
                            }
                            else {
                                // we've found a match
                                
                                if (curr_match_start_offset == numeric_limits<size_t>::max()) {
                                    // we're starting a new match, update the
                                    curr_match_start_offset = offset_on_curr_node;
                                    match_start_mapping_idx = i;
                                    match_start_edit_idx = j;
                                }
                                
                                // update the length of the match
                                curr_match_length += edit.from_length();
                            }
                            
                            // update our positional trackers
                            offset_on_curr_node += edit.from_length();
                            cumul_from_length += edit.from_length();
                            cumul_to_length += edit.to_length();
                        }
                    }
                    
                    if (curr_match_length > 0) {
                        // we were still in a match when we finished up, so we want to finish off the anchor
                        create_synthetic_anchor_node(path.mapping_size() - 1, path.mapping(path.mapping_size() - 1).edit_size());
                    }
                    
                    if (!handling_right_tail && prev_anchor_path_node != numeric_limits<size_t>::max()) {
                        // we need to make an edge from the final new anchor to the anchor we pinned to
                        path_nodes[prev_anchor_path_node].edges.emplace_back(attached_path_node_index,
                                                                             cumul_from_length - prev_anchor_final_from_length);
#ifdef debug_multipath_alignment
                        cerr << "adding final edge to " << attached_path_node_index << " with length " << cumul_from_length - prev_anchor_final_from_length << endl;
                        
#endif
                    }
                }
            }
        }
        
        // Now we've created new PathNodes for all the perfect matches in the tail alignments.
        // They can be resected out of snarls just like the original ones.
    }

    void MultipathAlignmentGraph::add_reachability_edges(const HandleGraph& graph,
                                                         const function<pair<id_t, bool>(id_t)>& project,
                                                         const unordered_multimap<id_t, pair<id_t, bool>>& injection_trans,
                                                         vector<size_t>* path_node_provenance) {
                                                         
        
        // We're going to make "reachability" edges, which connect MEMs (which
        // also may just be path segments, or trimmed MEMs) where both MEMs can
        // be part of the same alignment traceback, one after the other. For
        // that to be true, the MEMs have to be "colinear": the first comes
        // before the second in both the read and the graph. The MEMs also have
        // to not have any intervening MEMs where the intervening MEM is
        // reachable from the first MEM and the second MEM is reachable from
        // the intervening MEM.
        
        // We think in terms of "starts" (places in the graph and read where a
        // MEM begins, and "ends" (places in the graph and read where a MEM
        // stops). We always work in the read's local forward orientation.
        
        // MEM paths in the graph may visit graph nodes in any orientation. We
        // probably assume that the graph is dagified so everything flows in a
        // local forward orientation.
        
        // Our MEMs all live in path_nodes, and are identified by their indexes
        // there.
        
        
#ifdef debug_multipath_alignment
        cerr << "computing reachability" << endl;
#endif
        
        // Don't let people do this twice.
        assert(!has_reachability_edges);
        
        // optimization: we never add edges unless there are multiple nodes, and frequently there is only one
        // so we can skip traversing over the entire graph
        if (path_nodes.size() <= 1) {
#ifdef debug_multipath_alignment
            cerr << "skipping reachability computation because there are " << path_nodes.size() << " path nodes" << endl;
#endif
            has_reachability_edges = true;
            return;
        }
        
        // now we calculate reachability between the walked paths so we know which ones
        // to connect with intervening alignments
        
        /// Get the offset in the first visited graph node at which the given MEM starts.
        /// Does not account for orientation.
        auto start_offset = [&](size_t idx) {
            return path_nodes[idx].path.mapping(0).position().offset();
        };
        
        /// Get the offset in the first visited graph node at which the given MEM ends (i.e. the past-the-end offset).
        /// Does not account for orientation.
        auto end_offset = [&](size_t idx) {
            path_t& path = path_nodes[idx].path;
            const path_mapping_t& mapping = path.mapping(path.mapping_size() - 1);
            return mapping.position().offset() + mapping_from_length(mapping);
        };
        
        /// Get the ID of the first node visited in the graph along the path for a MEM.
        /// Does not account for orientation.
        auto start_node_id = [&](size_t idx) {
            return path_nodes[idx].path.mapping(0).position().node_id();
        };
        
        /// Get the ID of the last node visited in the graph along the path for a MEM.
        /// Does not account for orientation.
        auto end_node_id = [&](size_t idx) {
            path_t& path = path_nodes[idx].path;
            return path.mapping(path.mapping_size() - 1).position().node_id();
        };
        
        /// Get the offset in the read of either the start or past-the-end position of the given MEM, according to the end flag.
        auto endpoint_offset = [&](size_t idx, bool end) {
            return end ? end_offset(idx) : start_offset(idx);
        };
        
        /// Get the node ID in the graph of either the start or end position of the given MEM, according to the end flag.
        auto endpoint_node_id = [&](size_t idx, bool end) {
            return end ? end_node_id(idx) : start_node_id(idx);
        };
        
        // record the start and end node ids of every path
        // Maps from node ID to the list of MEM numbers that start on that node.
        unordered_map<id_t, vector<size_t>> path_starts;
        // Maps from node ID to the list of MEM numbers that end on that node.
        unordered_map<id_t, vector<size_t>> path_ends;
        for (size_t i = 0; i < path_nodes.size(); i++) {
            path_t& path = path_nodes[i].path;
            path_starts[path.mapping(0).position().node_id()].push_back(i);
            path_ends[path.mapping(path.mapping_size() - 1).position().node_id()].push_back(i);
        }
        
#ifdef debug_multipath_alignment
        cerr << "recorded starts: " << endl;
        for (const auto& rec : path_starts) {
            cerr << "\t" << "Node " << rec.first << ": ";
            for (auto l : rec.second) {
                cerr << "M" << l << " ";
            }
            cerr << endl;
        }
        
        cerr << "recorded ends: " << endl;
        for (const auto& rec : path_ends) {
            cerr << "\t" << "Node " << rec.first << ":  ";
            for (auto l : rec.second) {
                cerr << "M" << l << " ";
            }
            cerr << endl;
        }
#endif
        
        
        // Sort the MEMs starting and ending on each node in node sequence order.
        // MEMs that start/end earlier will appear earlier in the vector for the node they start/end on.
        for (pair<const id_t, vector<size_t>>& node_starts : path_starts) {
            std::sort(node_starts.second.begin(), node_starts.second.end(),
                      [&](const size_t idx_1, const size_t idx_2) {
                return start_offset(idx_1) < start_offset(idx_2);
            });
        }
        for (pair<const id_t, vector<size_t>>& node_ends : path_ends) {
            std::sort(node_ends.second.begin(), node_ends.second.end(),
                      [&](const size_t idx_1, const size_t idx_2) {
                return end_offset(idx_1) < end_offset(idx_2);
            });
        }
        
        // The "ranges" that are used below (range_start, range_end, etc.)
        // refer to intervals in these sorted per-node lists in path_starts and
        // path_ends corresponding to sets of MEMs that all start or end at the
        // same position on the same graph node.
        
        // We want to distinguish MEM numbers in path_nodes from indexes in path_starts and path_ends that we put ranges over.
        // So we prefix path_nodes-space MEM numbers with "M" in the debug output.
        // TODO: rename the variables so we can tell which size_ts have which semantics.
        // Or use a using to invent some semantic types.
        
        // some structures we will fill out with DP:
        
        // for each node, the starts and ends of MEMs that can reach this node
        // without passing any other MEM starts or ends. TODO: What do the
        // unordered_maps map from/to?
        unordered_map<id_t, unordered_map<size_t, size_t>> reachable_ends;
        unordered_map<id_t, unordered_map<size_t, size_t>> reachable_starts;
        
        // for the start of each MEM, the starts and ends of other MEMs that can reach it without passing any
        // other start or end
        // TODO: in what space is the MEM start size_t? Read space?
        // The pairs are pairs of MEM start and end positions (TODO: in read space?)
        // TODO: Do the other MEMs reach this MEM without passing other MEMs counting from their starts or their ends?
        unordered_map<size_t, vector<pair<size_t, size_t>>> reachable_starts_from_start;
        unordered_map<size_t, vector<pair<size_t, size_t>>> reachable_ends_from_start;
        
        // for the end of each MEM, the ends of other MEMs that can reach it without passing any
        // other start or end
        // TODO: in what space is the MEM end size_t? Read space?
        // The pairs are pairs of MEM start and end positions (TODO: in read space?)
        // TODO: Do the other MEMs reach this MEM without passing other MEMs counting from their starts or their ends?
        unordered_map<size_t, vector<pair<size_t, size_t>>> reachable_ends_from_end;
        unordered_map<size_t, vector<pair<size_t, size_t>>> reachable_starts_from_end;
        
        // get a topological order over the nodes in the graph to iterate over
        vector<handle_t> topological_order = handlealgs::lazier_topological_order(&graph);
        for (int64_t i = 0; i < topological_order.size(); i++) {
            id_t node_id = graph.get_id(topological_order[i]);
            
#ifdef debug_multipath_alignment
            cerr << "DP step for graph node " << node_id << endl;
#endif
            
            size_t node_length = graph.get_length(topological_order[i]);
            
            // do any MEMs start or end on this node?
            bool contains_starts = path_starts.count(node_id);
            bool contains_ends = path_ends.count(node_id);
            
            // we will use DP to carry reachability information forward onto the next nodes
            vector<handle_t> nexts;
            graph.follow_edges(topological_order[i], false, [&](const handle_t& next) {
                nexts.push_back(next);
            });
            
            if (contains_starts && contains_ends) {
                // since there are both starts and ends on this node, we have to traverse both lists simultaneously
                // to assess reachability within the same node
                
#ifdef debug_multipath_alignment
                cerr << "\tnode " << node_id << " contains both starts and ends of MEMs" << endl;
#endif
                
                vector<size_t>& ends = path_ends[node_id];
                vector<size_t>& starts = path_starts[node_id];
                
                
                // find the range of starts and ends in the list with the same offset
                
                size_t start_range_begin = 0;
                size_t start_range_end = 1;
                size_t end_range_begin = 0;
                size_t end_range_end = 1;
                
                size_t curr_start_offset = start_offset(starts[start_range_begin]);
                size_t curr_end_offset = end_offset(ends[end_range_begin]);
                size_t prev_offset = 0;
                
                while (end_range_end == ends.size() ? false : end_offset(ends[end_range_end]) == curr_end_offset) {
                    end_range_end++;
                }
                while (start_range_end == starts.size() ? false : start_offset(starts[start_range_end]) == curr_start_offset) {
                    start_range_end++;
                }
                
#ifdef debug_multipath_alignment
                cerr << "\tMEMs " << start_range_begin << ":" << start_range_end << " ordered by start position start at initial offset " << curr_start_offset << endl;
                cerr << "\tMEMs " << end_range_begin << ":" << end_range_end << " ordered by end position end at initial offset " << curr_end_offset << endl;
#endif
                
                // connect the first range of starts or ends to the incoming starts and ends
                
                size_t prev_end_range_begin = end_range_begin;
                size_t prev_start_range_begin = start_range_begin;
                
                bool at_end = (curr_end_offset <= curr_start_offset);
                // TODO: What exactly do these variables hold?
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
                    
#ifdef debug_multipath_alignment
                    cerr << "\tfirst endpoint is an end" << endl;
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
                    
#ifdef debug_multipath_alignment
                    cerr << "\tfirst endpoint is a start" << endl;
#endif
                }
                
                for (size_t j = *range_begin; j < *range_end; j++) {
                    for (const pair<size_t, size_t>& incoming_end : reachable_ends[node_id]) {
#ifdef debug_multipath_alignment
                        cerr << "\t\tidentifying end of M" << incoming_end.first << " as reachable from endpoint of M" << endpoints->at(j) << endl;
#endif
                        (*reachable_ends_from_endpoint)[endpoints->at(j)].emplace_back(incoming_end.first, incoming_end.second + *curr_offset);
                    }
                    for (const pair<size_t, size_t>& incoming_start : reachable_starts[node_id]) {
#ifdef debug_multipath_alignment
                        cerr << "\t\tidentifying start of M" << incoming_start.first << " as reachable from endpoint of M" << endpoints->at(j) << endl;
#endif
                        (*reachable_starts_from_endpoint)[endpoints->at(j)].emplace_back(incoming_start.first, incoming_start.second + *curr_offset);
                    }
                }
                
                
                bool prev_is_end = at_end;
                *range_begin = *range_end;
                prev_offset = *curr_offset;
                if (*range_begin != endpoints->size()) {
                    *curr_offset = endpoint_offset(endpoints->at(*range_begin), at_end);
                    while (*range_end == endpoints->size() ? false : endpoint_offset(endpoints->at(*range_end), at_end) == *curr_offset) {
                        (*range_end)++;
                    }
                }
                
#ifdef debug_multipath_alignment
                cerr << "\tnext look at MEMs " << *range_begin << ":" << *range_end << " ordered by start or end position, at offset " << *curr_offset << endl;
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
#ifdef debug_multipath_alignment
                    cerr << "\tat MEMs " << *range_begin << ":" << *range_end << " ordered by " << (at_end ? "end" : "start") << "s" << endl;
#endif
                    
                    size_t dist_between = *curr_offset - prev_offset;
                    
                    // connect this range to the previous range
                    if (prev_is_end) {
#ifdef debug_multipath_alignment
                        cerr << "\t\tlooking backwards to ends of MEMs " << prev_end_range_begin << ":" << end_range_begin << " ordered by ends" << endl;
#endif
                        for (size_t j = prev_end_range_begin; j < end_range_begin; j++) {
                            for (size_t k = *range_begin; k < *range_end; k++) {
#ifdef debug_multipath_alignment
                                cerr << "\t\tidentifying end of M" << ends[j] << " as reachable from " << (at_end ? "end" : "start") << " of M" << endpoints->at(k) << endl;
#endif
                                (*reachable_ends_from_endpoint)[endpoints->at(k)].emplace_back(ends[j], dist_between);
                            }
                        }
                    }
                    else {
#ifdef debug_multipath_alignment
                        cerr << "\t\tlooking backwards to starts in range " << prev_start_range_begin << ":" << start_range_begin << endl;
#endif
                        for (size_t j = prev_start_range_begin; j < start_range_begin; j++) {
                            for (size_t k = *range_begin; k < *range_end; k++) {
#ifdef debug_multipath_alignment
                                cerr << "\t\tidentifying start of M" << starts[j] << " as reachable from " << (at_end ? "end" : "start") << " of M" << endpoints->at(k) << endl;
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
                        while (*range_end == endpoints->size() ? false : endpoint_offset(endpoints->at(*range_end), at_end) == *curr_offset) {
                            (*range_end)++;
                        }
                    }
                    
#ifdef debug_multipath_alignment
                    cerr << "\tnext look at MEMS " << *range_begin << ":" << *range_end << " ordered by start or end position, at offset " << *curr_offset << endl;
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
                    
#ifdef debug_multipath_alignment
                    cerr << "\tfinal endpoint(s) are end(s)" << endl;
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
                    
#ifdef debug_multipath_alignment
                    cerr << "\tfinal endpoint(s) are start(s)" << endl;
#endif
                }
                
                while (*range_begin < endpoints->size()) {
                    
                    size_t dist_between = *curr_offset - prev_offset;
                    
                    if (prev_is_end) {
#ifdef debug_multipath_alignment
                        cerr << "\t\tlooking backwards to MEMs " << prev_end_range_begin << ":" << end_range_begin << " ordered by ends" << endl;
#endif
                        for (size_t j = prev_end_range_begin; j < end_range_begin; j++) {
                            for (size_t k = *range_begin; k < *range_end; k++) {
#ifdef debug_multipath_alignment
                                cerr << "\t\tidentifying end of M" << ends[j] << " as reachable from endpoint of M" << endpoints->at(k) << endl;
#endif
                                (*reachable_ends_from_endpoint)[endpoints->at(k)].push_back(make_pair(ends[j], dist_between));
                            }
                        }
                    }
                    else {
#ifdef debug_multipath_alignment
                        cerr << "\t\tlooking backwards to MEMs " << prev_start_range_begin << ":" << start_range_begin << " ordered by starts" << endl;
#endif
                        for (size_t j = prev_start_range_begin; j < start_range_begin; j++) {
                            for (size_t k = *range_begin; k < *range_end; k++) {
#ifdef debug_multipath_alignment
                                cerr << "\t\tidentifying start of M" << starts[j] << " as reachable from endpoint of M" << endpoints->at(k) << endl;
#endif
                                (*reachable_starts_from_endpoint)[endpoints->at(k)].push_back(make_pair(starts[j], dist_between));
                            }
                        }
                    }
                    
#ifdef debug_multipath_alignment
                    cerr << "\tmoving to next endpoint range" << endl;
#endif
                    
                    *prev_range_begin = *range_begin;
                    *range_begin = *range_end;
                    prev_offset = *curr_offset;
                    prev_is_end = at_end;
                    
                    if (*range_begin != endpoints->size()) {
                        *curr_offset = endpoint_offset(endpoints->at(*range_begin), at_end);
                        while (*range_end == endpoints->size() ? false : endpoint_offset(endpoints->at(*range_end), at_end) == *curr_offset) {
                            (*range_end)++;
                        }
                    }
                    
#ifdef debug_multipath_alignment
                    cerr << "\tnext look at MEMs " << *range_begin << ":" << *range_end << " ordered by start or end, at offset " << *curr_offset << endl;
#endif
                }
                
                // carry forward the reachability of the last range onto the next nodes
                size_t dist_thru = node_length - *curr_offset;
                
#ifdef debug_multipath_alignment
                cerr << "\tcarrying forward reachability onto next nodes at distance " << dist_thru << endl;
#endif
                
                for (const handle_t& next : nexts) {
                    unordered_map<size_t, size_t>& reachable_endpoints_next = (*reachable_endpoints)[graph.get_id(next)];
                    for (size_t j = *prev_range_begin; j < endpoints->size(); j++) {
                        if (reachable_endpoints_next.count(endpoints->at(j))) {
                            reachable_endpoints_next[endpoints->at(j)] = std::min(reachable_endpoints_next[endpoints->at(j)], dist_thru);
                        }
                        else {
                            reachable_endpoints_next[endpoints->at(j)] = dist_thru;
                        }
                        
#ifdef debug_multipath_alignment
                        cerr << "\t\t" << "endpoint of M" << endpoints->at(j) << " at dist " << reachable_endpoints_next[endpoints->at(j)] << " to node " << graph.get_id(next) << endl;
#endif
                        
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
#ifdef debug_multipath_alignment
                    cerr << "\tnode " << node_id << " contains only ends of MEMs" << endl;
#endif
                    reachable_endpoints = &reachable_ends;
                    reachable_starts_from_endpoint = &reachable_starts_from_end;
                    reachable_ends_from_endpoint = &reachable_ends_from_end;
                    reachable_endpoints_from_endpoint = &reachable_ends_from_end;
                    endpoints = &path_ends[node_id];
                }
                else {
#ifdef debug_multipath_alignment
                    cerr << "\tnode " << node_id << " contains only starts of MEMs" << endl;
#endif
                    reachable_endpoints = &reachable_starts;
                    reachable_starts_from_endpoint = &reachable_starts_from_start;
                    reachable_ends_from_endpoint = &reachable_ends_from_start;
                    reachable_endpoints_from_endpoint = &reachable_starts_from_start;
                    endpoints = &path_starts[node_id];
                }
                
                // the starts/ends coming into this node from outside
                size_t range_begin = 0;
                size_t range_end = 1;
                size_t curr_offset = endpoint_offset(endpoints->at(range_begin), contains_ends);
                size_t prev_offset = curr_offset;
                // find the range of endpoints that are at the first offset
                while (range_end < endpoints->size() && endpoint_offset(endpoints->at(range_end), contains_ends) == curr_offset) {
                    range_end++;
                }
                
#ifdef debug_multipath_alignment
                cerr << "\tMEMs " << range_begin << ":" << range_end << " ordered by start or end are at initial offset " << curr_offset << endl;
#endif
                
                // connect the range to the incoming starts/ends
                for (size_t j = range_begin; j < range_end; j++) {
                    for (const pair<size_t, size_t>& incoming_start : reachable_starts[node_id]) {
#ifdef debug_multipath_alignment
                        cerr << "\t\tidentifying start of M" << incoming_start.first << " as reachable from " << (contains_ends ? "end" : "start") << " of M" << endpoints->at(j) << endl;
#endif
                        (*reachable_starts_from_endpoint)[endpoints->at(j)].emplace_back(incoming_start.first, incoming_start.second + curr_offset);
                    }
                    for (const pair<size_t, size_t>& incoming_end : reachable_ends[node_id]) {
#ifdef debug_multipath_alignment
                        cerr << "\t\tidentifying end of M" << incoming_end.first << " as reachable from " << (contains_ends ? "end" : "start") << " of M" << endpoints->at(j) << endl;
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
                    while (range_end < endpoints->size() && endpoint_offset(endpoints->at(range_end), contains_ends) == curr_offset) {
                        range_end++;
                    }
                    
#ifdef debug_multipath_alignment
                    cerr << "\tnext look at MEMs " << range_begin << ":" << range_end << " ordered by start or end, at offset " << curr_offset << endl;
#endif
                    
                    size_t dist_between = curr_offset - prev_offset;
                    
                    // connect this range to the previous range
                    for (size_t j = range_begin; j < range_end; j++) {
                        for (size_t k = prev_range_begin; k < range_begin; k++) {
#ifdef debug_multipath_alignment
                            cerr << "\t\tidentifying " << (contains_ends ? "end" : "start") << " of M" << endpoints->at(k) << " as reachable from " << (contains_ends ? "end" : "start") << " of M" << endpoints->at(j) << endl;
#endif
                            (*reachable_endpoints_from_endpoint)[endpoints->at(j)].push_back(make_pair(endpoints->at(k), dist_between));
                        }
                    }
                    prev_range_begin = range_begin;
                    range_begin = range_end;
                }
                
                // this node contains at least one endpoint of a MEM, so carry forward the reachability of all
                // endpoints at the final offset onto the next nodes
                
                size_t dist_thru = node_length - curr_offset;
                
#ifdef debug_multipath_alignment
                cerr << "\tcarrying forward reachability onto next nodes at distance " << dist_thru << endl;
#endif
                
                for (const handle_t& next : nexts) {
                    
                    unordered_map<size_t, size_t>& reachable_endpoints_next = (*reachable_endpoints)[graph.get_id(next)];
                    for (size_t j = prev_range_begin; j < endpoints->size(); j++) {
                        if (reachable_endpoints_next.count(endpoints->at(j))) {
                            reachable_endpoints_next[endpoints->at(j)] = std::min(reachable_endpoints_next[endpoints->at(j)], dist_thru);
                        }
                        else {
                            reachable_endpoints_next[endpoints->at(j)] = dist_thru;
                        }
                        
#ifdef debug_multipath_alignment
                        cerr << "\t\t" << (contains_ends ? "end" : "start") << " of M" << endpoints->at(j) << " at dist " << reachable_endpoints_next[endpoints->at(j)] << " to node " << graph.get_id(next) << endl;
#endif
                    }
                }
            }
            else {
                // this node doesn't contain the start or end of any MEM, so we carry forward the reachability
                // into this node onto the next nodes
                
#ifdef debug_multipath_alignment
                cerr << "\tnode " << node_id << " does not contain starts or ends of MEMs, carrying forward reachability" << endl;
#endif
                
                for (const handle_t& next : nexts) {
                    unordered_map<size_t, size_t>& reachable_ends_next = reachable_ends[graph.get_id(next)];
                    for (const pair<size_t, size_t>& reachable_end : reachable_ends[node_id]) {
                        size_t dist_thru = reachable_end.second + node_length;
#ifdef debug_multipath_alignment
                        cerr << "\t\tend of M" << reachable_end.first << " at dist " << dist_thru << " to node " << graph.get_id(next) << endl;
#endif
                        if (reachable_ends_next.count(reachable_end.first)) {
                            reachable_ends_next[reachable_end.first] = std::min(reachable_ends_next[reachable_end.first],
                                                                                dist_thru);
                        }
                        else {
                            reachable_ends_next[reachable_end.first] = dist_thru;
                        }
                    }
                    
                    unordered_map<size_t, size_t>& reachable_starts_next = reachable_starts[graph.get_id(next)];
                    for (const pair<size_t, size_t>& reachable_start : reachable_starts[node_id]) {
                        size_t dist_thru = reachable_start.second + node_length;
#ifdef debug_multipath_alignment
                        cerr << "\t\tstart of M" << reachable_start.first << " at dist " << dist_thru << " to node " << graph.get_id(next) << endl;
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
        
#ifdef debug_multipath_alignment
        cerr << "final reachability:" << endl;
        cerr << "\tstarts from starts" << endl;
        for (const auto& record : reachable_starts_from_start) {
            cerr << "\t\tstart of M" << record.first << " can reach:" << endl;
            for (const auto& endpoint : record.second) {
                cerr << "\t\t\tstart of M" << endpoint.first << " (dist " << endpoint.second << ")" << endl;
            }
        }
        cerr << "\tends from starts" << endl;
        for (const auto& record : reachable_ends_from_start) {
            cerr << "\t\tstart of M" << record.first << " can reach:" << endl;
            for (const auto& endpoint : record.second) {
                cerr << "\t\t\tend of M" << endpoint.first << " (dist " << endpoint.second << ")" << endl;
            }
        }
        cerr << "\tstarts from ends" << endl;
        for (const auto& record : reachable_starts_from_end) {
            cerr << "\t\tend of M" << record.first << " can reach:" << endl;
            for (const auto& endpoint : record.second) {
                cerr << "\t\t\tstart of M" << endpoint.first << " (dist " << endpoint.second << ")" << endl;
            }
        }
        cerr << "\tends from ends" << endl;
        for (const auto& record : reachable_ends_from_end) {
            cerr << "\t\tend of M" << record.first << " can reach:" << endl;
            for (const auto& endpoint : record.second) {
                cerr << "\t\t\tend of M" << endpoint.first << " (dist " << endpoint.second << ")" << endl;
            }
        }
        cerr << "setting up structure of MEM graph" << endl;
#endif
        
        // now we have the reachability information for the start and end of every MEM in the graph. we
        // will use this to navigate between the MEMs in a way that respects graph reachability so that this
        // phase of the algorithm only needs to pay attention to read colinearity and transitive reducibility
        
        vector<unordered_map<size_t, size_t>> noncolinear_shells(path_nodes.size());
        
        // map from index_from to maps of index_onto to (overlap to length, overlap from length, dist)
        unordered_map<size_t, map<size_t, tuple<size_t, size_t, size_t>>> confirmed_overlaps;
        // map from path index to set of indexes whose start occurs on the path
        unordered_map<size_t, set<size_t>> path_starts_on_path;
        
        for (size_t i = 0; i < topological_order.size(); i++) {
            id_t node_id = graph.get_id(topological_order[i]);
            
#ifdef debug_multipath_alignment
            cerr << "looking for edges for starts on node " << node_id << endl;
#endif
            
            if (!path_starts.count(node_id) && !path_ends.count(node_id)) {
#ifdef debug_multipath_alignment
                cerr << "there are no starts or ends on this node" << endl;
#endif
                continue;
            }
            
            // keep track of the starts that are located at the same offset at earlier positions in the vector of
            // of starts (used later in the overlap finding step)
            vector<size_t> colocated_starts;
            
            vector<size_t>& starts = path_starts[node_id];
            vector<size_t>& ends = path_ends[node_id];
            size_t start_idx = 0, end_idx = 0;
            size_t curr_start_offset = numeric_limits<size_t>::max(), curr_end_offset = numeric_limits<size_t>::max();
            if (!starts.empty()) {
                curr_start_offset = start_offset(starts[0]);
            }
            if (!ends.empty()) {
                curr_end_offset = end_offset(ends[0]);
            }
            while (start_idx < starts.size() || end_idx < ends.size()) {
                
                // TODO: would it be better to combine these into one queue?
                
                // initialize queues for the next start and next end, prioritized by shortest distance
                structures::RankPairingHeap<size_t, size_t, std::greater<size_t>> start_queue, end_queue;
                
                if (curr_start_offset >= curr_end_offset) {
                    
                    // the next endpoint is an end, the point of searching backwards from these is to
                    // fill out the non-colinear shell of the current end with any path whose end is between
                    // this path's start and end but is not overlap colinear
                    
                    size_t end = ends[end_idx];
                    
#ifdef debug_multipath_alignment
                    cerr << "searching backward from end " << end << endl;
#endif
                    PathNode& end_node = path_nodes[end];
                    unordered_map<size_t, size_t>& noncolinear_shell = noncolinear_shells[end];
                    
                    for (const pair<size_t, size_t>& next_end : reachable_ends_from_end[end]) {
                        end_queue.push_or_reprioritize(next_end.first, next_end.second);
                    }
                    
                    for (const pair<size_t, size_t>& start_next : reachable_starts_from_end[end]) {
                        start_queue.push_or_reprioritize(start_next.first, start_next.second);
                    }
                    
                    while (!start_queue.empty() || !end_queue.empty()) {
                        // is the next item on the queues a start or an end?
                        if (start_queue.empty() ? false : (end_queue.empty() ? true : start_queue.top().second < end_queue.top().second)) {
                            
                            // the next closest endpoint is a start, traverse through it to find ends (which is what we really want)
                            
                            pair<size_t, size_t> start_here = start_queue.top();
                            start_queue.pop();
                            
                            // don't keep looking backward earlier than the start of the current path
                            if (start_here.first == end) {
                                continue;
                            }
                            
                            // the minimum distance to each of the starts or ends this can reach is (at most) the sum of the min distance
                            // between them and the distance already traversed
                            for (const pair<size_t, size_t>& next_end : reachable_ends_from_start[start_here.first]) {
                                end_queue.push_or_reprioritize(next_end.first, start_here.second + next_end.second);
                            }
                            
                            for (const pair<size_t, size_t>& start_next : reachable_starts_from_start[start_here.first]) {
                                start_queue.push_or_reprioritize(start_next.first, start_here.second + start_next.second);
                            }
                        }
                        else {
                            
                            // the next closest endpoint is an end, so we'll check whether we can
                            
                            pair<size_t, size_t> end_here = end_queue.top();
                            end_queue.pop();
                            
                            PathNode& next_end_node = path_nodes[end_here.first];
                            
                            // these are non-colinear, so add it to the non-colinear shell
                            if (next_end_node.begin >= end_node.begin || next_end_node.end >= end_node.end) {
                                if (noncolinear_shell.count(end_here.first)) {
                                    noncolinear_shell[end_here.first] = std::min(end_here.second, noncolinear_shell[end_here.first]);
                                }
                                else {
                                    noncolinear_shell[end_here.first] = end_here.second;
                                }
                                continue;
                            }
                            
                            // if we get this far, the two paths are colinear or overlap-colinear, so we won't add it to the
                            // non-colinear shell. now we need to decide whether to keep searching backward. we'll check a
                            // few conditions that will guarantee that the rest of the search is redundant
                            
                            // TODO: this actually isn't a full set of criteria, we don't just want to know if there is an
                            // edge, we want to know if it is reachable along any series of edges...
                            // at least this will only cause a few false positive edges that we can remove later with
                            // the transitive reduction
                            
                            // see if this node has an edge forward
                            bool has_edge_forward = false;
                            for (auto& edge : next_end_node.edges) {
                                if (edge.first == end) {
                                    has_edge_forward = true;
                                    break;
                                }
                            }
                            
                            if (has_edge_forward) { // already has an edge, can stop
                                continue;
                            }
                            
                            auto overlap_iter = confirmed_overlaps.find(end_here.first);
                            if (overlap_iter != confirmed_overlaps.end()) {
                                if (overlap_iter->second.count(end)) {
                                    has_edge_forward = true;
                                    break;
                                }
                            }
                            
                            if (has_edge_forward) { // already has an overlap, can stop
                                continue;
                            }
                            
                            // we can't easily guarantee that this is non-colinear or colinear, so we're just going to treat it
                            // as non-colinear and accept some risk of this creating redundant edges to later nodes
                            if (noncolinear_shell.count(end_here.first)) {
                                noncolinear_shell[end_here.first] = std::min(end_here.second, noncolinear_shell[end_here.first]);
                            }
                            else {
                                noncolinear_shell[end_here.first] = end_here.second;
                            }
                        }
                    }
                    
                    end_idx++;
                    curr_end_offset = (end_idx == ends.size() ? numeric_limits<size_t>::max() : end_offset(ends[end_idx]));
                }
                else {
                    
                    size_t start = starts[start_idx];
                    
#ifdef debug_multipath_alignment
                    cerr << "searching backward from start " << start << " at index " << start_idx << endl;
#endif
                    
                    PathNode& start_node = path_nodes[start];
                    unordered_map<size_t, size_t>& noncolinear_shell = noncolinear_shells[start];
                    // TODO: kinda ugly
                    // init this to 0, we'll actually compute it if we need it ever
                    size_t start_node_from_length = 0;

                    // we begin at the start we're searching backward from
                    start_queue.push_or_reprioritize(start, 0);
                    
                    while (!start_queue.empty() || !end_queue.empty()) {
                        // is the next item on the queues a start or an end?
                        if (!start_queue.empty() && (end_queue.empty() || start_queue.top().second < end_queue.top().second)) {
                            
                            // the next closest endpoint is a start, traverse through it to find ends (which is what we really want)
                            
                            pair<size_t, size_t> start_here = start_queue.top();
                            start_queue.pop();
                            
#ifdef debug_multipath_alignment
                            cerr << "traversing start " << start_here.first << " at distance " << start_here.second << endl;
#endif
                            
                            // the minimum distance to each of the starts or ends this can reach is the sum of the min distance
                            // between them and the distance already traversed
                            for (const pair<size_t, size_t>& end : reachable_ends_from_start[start_here.first]) {
                                end_queue.push_or_reprioritize(end.first, start_here.second + end.second);
#ifdef debug_multipath_alignment
                                cerr << "found reachable end " << end.first << " at distance " << start_here.second + end.second << endl;
#endif
                            }
                            
                            for (const pair<size_t, size_t>& start_next : reachable_starts_from_start[start_here.first]) {
                                start_queue.push_or_reprioritize(start_next.first, start_here.second + start_next.second);
                            }
                        }
                        else {
                            
                            // the next closest endpoint is an end, so we check if we can make a connection to the start
                            // that we're searching backward from
                            
                            size_t candidate_end, candidate_dist;
                            tie(candidate_end, candidate_dist) = end_queue.top();
                            end_queue.pop();
                            
#ifdef debug_multipath_alignment
                            cerr << "considering end " << candidate_end << " as candidate for edge of dist " << candidate_dist << endl;
#endif
                            
                            PathNode& candidate_end_node = path_nodes[candidate_end];
                            
                            if (candidate_end_node.end <= start_node.begin) {
                                // these MEMs are read colinear and graph reachable
                                if (start != candidate_end) {
                                    // and they are not the same path node, so add an edge
                                    // (this is almost always the case, but some code paths will add empty read sequences
                                    // to anchor to specific locations, which slightly confuses the reachability logic)
                                    candidate_end_node.edges.emplace_back(start, candidate_dist);
                                    
#ifdef debug_multipath_alignment
                                    cerr << "connection is read colinear, adding edge on " << candidate_end << " for total of " << candidate_end_node.edges.size() << " edges so far" << endl;
                                    for (auto& edge : candidate_end_node.edges) {
                                        cerr << "\t-> " << edge.first << " dist " << edge.second << endl;
                                    }
#endif
                                }
                                
                                // skip to the predecessor's noncolinear shell, whose connections might not be blocked by
                                // this connection
                                for (const pair<size_t, size_t>& shell_pred : noncolinear_shells[candidate_end]) {
#ifdef debug_multipath_alignment
                                    cerr << "enqueueing " << shell_pred.first << " at dist " << shell_pred.second + candidate_dist << " from noncolinear shell" << endl;
#endif
                                    end_queue.push_or_reprioritize(shell_pred.first, candidate_dist + shell_pred.second);
                                }
                            }
                            else if (candidate_end_node.begin < start_node.begin && candidate_end_node.end < start_node.end) {
                                // the MEM can be made colinear by removing an overlap, which will not threaten reachability
                                size_t read_overlap = candidate_end_node.end - start_node.begin;
                                size_t graph_overlap = corresponding_from_length(start_node.path, read_overlap, false);
                                confirmed_overlaps[start][candidate_end] = make_tuple(read_overlap, graph_overlap,
                                                                                      candidate_dist + graph_overlap);
                                
#ifdef debug_multipath_alignment
                                cerr << "connection is overlap colinear, recording to add edge later" << endl;
#endif
                                
                                // the end of this node might not actually block connections since it's going to intersect the middle of the node
                                // so we need to find predecessors to this end too
                                
                                // add any ends directly reachable from the end
                                for (const pair<size_t, size_t>& exposed_end : reachable_ends_from_end[candidate_end]) {
                                    end_queue.push_or_reprioritize(exposed_end.first, candidate_dist + exposed_end.second);
#ifdef debug_multipath_alignment
                                    cerr << "found reachable exposed end " << exposed_end.first << " at distance " << candidate_dist + exposed_end.second << endl;
#endif
                                }
                                
                                // add the directly reachable exposed starts to the queue
                                for (const pair<size_t, size_t>& exposed_start : reachable_starts_from_end[candidate_end]) {
#ifdef debug_multipath_alignment
                                    cerr << "adding exposed start traversal with " << exposed_start.first << " at distance " << candidate_dist + exposed_start.second << endl;
#endif
                                    start_queue.push_or_reprioritize(exposed_start.first, candidate_dist + exposed_start.second);
                                }
                                
                                // also skip to the predecessor's noncolinear shell, whose connections might not be blocked by
                                // this connection
                                for (const pair<size_t, size_t>& shell_pred : noncolinear_shells[candidate_end]) {
#ifdef debug_multipath_alignment
                                    cerr << "enqueueing " << shell_pred.first << " at dist " << candidate_dist + shell_pred.second << " from noncolinear shell" << endl;
#endif
                                    end_queue.push_or_reprioritize(shell_pred.first, candidate_dist + shell_pred.second);
                                }
                            }
                            else {
                                // these MEMs are noncolinear, so add this predecessor to the noncolinear shell
                                if (start_node_from_length == 0) {
                                    start_node_from_length = path_from_length(start_node.path);
                                }
                                if (noncolinear_shell.count(candidate_end)) {
                                    noncolinear_shell[candidate_end] = std::min(candidate_dist + start_node_from_length,
                                                                                noncolinear_shell[candidate_end]);
                                }
                                else {
                                    noncolinear_shell[candidate_end] = candidate_dist + start_node_from_length;
                                }
                                
#ifdef debug_multipath_alignment
                                cerr << "connection is noncolinear, add to shell at dist " << candidate_dist + start_node_from_length << " and continue to search backwards" << endl;
#endif
                                
                                // there is no connection to block further connections back, so any of this MEMs
                                // predecessors could still be colinear
                                
                                // find the ends that can reach it directly
                                for (const pair<size_t, size_t>& pred_end : reachable_ends_from_end[candidate_end]) {
                                    end_queue.push_or_reprioritize(pred_end.first, candidate_dist + pred_end.second);
#ifdef debug_multipath_alignment
                                    cerr << "found reachable end " << pred_end.first << " at distance " << candidate_dist + pred_end.second << endl;
#endif
                                }
                                
                                // set the start queue up with the immediate start neighbors
                                for (const pair<size_t, size_t>& pred_start : reachable_starts_from_end[candidate_end]) {
                                    start_queue.push_or_reprioritize(pred_start.first, candidate_dist + pred_start.second);
                                }
                            }
                        }
                    }
                    
#ifdef debug_multipath_alignment
                    cerr << "walking path to look for overlaps" << endl;
#endif
                    
                    path_t& path = start_node.path;
                    
                    // update the path starts index for the paths that start at the same position
                    for (size_t colocated_start : colocated_starts) {
                        path_starts_on_path[colocated_start].insert(start);
                        path_starts_on_path[start].insert(colocated_start);
                    }
                    // records of (node_idx, graph overlap length)
                    vector<pair<size_t, size_t>> overlap_candidates;
                    
                    if (path.mapping_size() == 1) {
                        // TODO: this edge case is a little duplicative, probably could merge
                        
#ifdef debug_multipath_alignment
                        cerr << "path is one mapping long" << endl;
#endif
                        
                        size_t final_offset = end_offset(start);
                        // record which starts are on the path on this node
                        for (size_t path_start_idx = start_idx + 1;
                             path_start_idx < starts.size() && start_offset(starts[path_start_idx]) < final_offset;
                             path_start_idx++) {
                            
                            path_starts_on_path[start].insert(starts[path_start_idx]);
                            
                        }
                        // record which ends are on the path on this node
                        for (size_t path_end_idx = end_idx; path_end_idx < ends.size(); path_end_idx++) {
                            size_t end_offset_here = end_offset(ends[path_end_idx]);
                            if (end_offset_here < final_offset) {
                                overlap_candidates.emplace_back(ends[path_end_idx], end_offset_here - curr_start_offset);
                            }
                            else {
                                break;
                            }
                        }
                    }
                    else {
#ifdef debug_multipath_alignment
                        cerr << "path is multiple mappings long" << endl;
#endif
                        
                        // record which starts are on the path on the first node
                        for (size_t path_start_idx = start_idx + 1; path_start_idx < starts.size(); path_start_idx++) {
                            path_starts_on_path[start].insert(starts[path_start_idx]);
                        }
                        // record which ends are on the path on the first node
                        for (size_t path_end_idx = end_idx; path_end_idx < ends.size(); path_end_idx++) {
                            overlap_candidates.emplace_back(ends[path_end_idx], end_offset(ends[path_end_idx]) - curr_start_offset);
                        }
                        size_t traversed_length = mapping_from_length(path.mapping(0));
                        
                        for (size_t j = 1; j + 1 < path.mapping_size(); j++) {
                            id_t path_node_id = path.mapping(j).position().node_id();
                            // record which starts are on the path on this node
                            for (size_t path_start : path_starts[path_node_id]) {
                                path_starts_on_path[start].insert(path_start);
                            }
                            // record which ends are on the path on this node
                            for (size_t path_end : path_ends[path_node_id]) {
                                overlap_candidates.emplace_back(path_end, end_offset(path_end) + traversed_length);
                            }
                            
                            traversed_length += mapping_from_length(path.mapping(j));
                        }
                        
                        id_t final_node_id = path.mapping(path.mapping_size() - 1).position().node_id();
                        vector<size_t>& final_starts = path_starts[final_node_id];
                        vector<size_t>& final_ends = path_ends[final_node_id];
                        
                        size_t final_offset = end_offset(start);
                        // record which starts are on the path on the last node
                        for (size_t path_start_idx = 0;
                             path_start_idx < final_starts.size() && start_offset(final_starts[path_start_idx]) < final_offset;
                             path_start_idx++) {
                            
                            path_starts_on_path[start].insert(final_starts[path_start_idx]);
                            
                        }
                        // record which ends are on the path on the last node
                        for (size_t path_end_idx = 0; path_end_idx < final_ends.size(); path_end_idx++) {
                            size_t end_offset_here = end_offset(final_ends[path_end_idx]);
                            if (end_offset_here < final_offset) {
                                overlap_candidates.emplace_back(final_ends[path_end_idx], end_offset_here + traversed_length);
                            }
                            else {
                                break;
                            }
                        }
                    }
                    
                    for (const pair<size_t, size_t>& overlap_candidate : overlap_candidates) {
#ifdef debug_multipath_alignment
                        cerr << "considering candidate overlap from " << overlap_candidate.first << " at dist " << overlap_candidate.second << endl;
#endif
                        
                        if (path_starts_on_path[start].count(overlap_candidate.first)) {
                            // the start of this MEM is also on the path, so this can't be an overhanging overlap
                            continue;
                        }
                        
                        if (!path_starts_on_path[overlap_candidate.first].count(start)) {
                            // the path we are walking doesn't start on the other path, so this can't be a full overlap
                            continue;
                        }
                        
                        PathNode& overlap_node = path_nodes[overlap_candidate.first];
                        
                        // how much do the paths overlap?
                        size_t graph_overlap = overlap_candidate.second;
                        size_t read_overlap = corresponding_to_length(path, graph_overlap, false);
                        
                        // are the paths read colinear after removing the overlap?
                        if (start_node.begin + read_overlap >= overlap_node.end) {
#ifdef debug_multipath_alignment
                            cerr << "confirmed overlap colinear with read overlap of " << read_overlap << ", graph overlap " << graph_overlap << endl;
#endif
                            confirmed_overlaps[start][overlap_candidate.first] = tuple<size_t, size_t, size_t>(read_overlap, graph_overlap, 0);
                        }
                        else if (overlap_node.begin < start_node.begin && overlap_node.end < start_node.end) {
#ifdef debug_multipath_alignment
                            cerr << "confirmed overlap colinear with longer read overlap of " << overlap_node.end - start_node.begin << endl;
#endif
                            // there is still an even longer read overlap we need to remove
                            size_t extended_read_overlap = overlap_node.end - start_node.begin;
                            size_t extended_graph_overlap = corresponding_from_length(path, extended_read_overlap, false);
                            confirmed_overlaps[start][overlap_candidate.first] = tuple<size_t, size_t, size_t>(extended_read_overlap,
                                                                                                               extended_graph_overlap,
                                                                                                               extended_graph_overlap - graph_overlap);
                        }
                        else {
#ifdef debug_multipath_alignment
                            cerr << "not colinear even with overlap, adding to non-colinear shell at distance " << overlap_candidate.second << endl;
#endif
                            // the overlapping node is still not reachable so it is in the noncolinear shell of this node
                            if (start_node_from_length == 0) {
                                start_node_from_length = path_from_length(start_node.path);
                            }
                            noncolinear_shell[overlap_candidate.first] = start_node_from_length - overlap_candidate.second;
                        }
                    }
                    
                    start_idx++;
                    size_t new_start_offset = (start_idx == starts.size() ? numeric_limits<size_t>::max() : start_offset(starts[start_idx]));
                    if (new_start_offset != curr_start_offset) {
                        colocated_starts.clear();
                    }
                    if (start_idx < starts.size()) {
                        colocated_starts.emplace_back(starts[start_idx]);
                    }
                    curr_start_offset = new_start_offset;
                }
            }
        }
        
#ifdef debug_multipath_alignment
        cerr << "breaking nodes at overlap edges" << endl;
#endif
        
        // now we've found all overlap edges, so we can add them into the graph in an order such that they don't
        // conflict (note that all overlap are from an earlier node onto a later one, so we don't need to worry
        // about overlaps coming in from both directions)
        
        // sort in descending order of overlap length and group by the node that is being cut among overlaps of same length
        // tuples of (read overlap, graph overlap, index onto, index from, distance)
        vector<tuple<size_t, size_t, size_t, size_t, size_t>> ordered_overlaps;
        for (const auto& path_overlaps : confirmed_overlaps) {
            for (const auto& overlap_record : path_overlaps.second) {
                ordered_overlaps.emplace_back(get<0>(overlap_record.second),
                                              get<1>(overlap_record.second),
                                              path_overlaps.first,
                                              overlap_record.first,
                                              get<2>(overlap_record.second));
            }
        }
        // because both from and to lengths are monotonic, we should never get into a situations where the
        // pair (read overlap, graph overlap) is incomparable, so this partial order is actually a total order
        // barring equal pairs
        sort(ordered_overlaps.begin(), ordered_overlaps.end(), greater<tuple<size_t, size_t, size_t, size_t, size_t>>());
        
        // keep track of whether another node is holding the suffix of one of the original nodes because of a split
        unordered_map<size_t, size_t> node_with_suffix;
        
        // split up each node with an overlap edge onto it
        auto iter = ordered_overlaps.begin();
        while (iter != ordered_overlaps.end()) {
            // find the range of overlaps that want to cut this node at the same place
            auto iter_range_end = iter;
            while (get<0>(*iter_range_end) == get<0>(*iter) && get<1>(*iter_range_end) == get<1>(*iter)
                   && get<2>(*iter_range_end) == get<2>(*iter)) {
                iter_range_end++;
                if (iter_range_end == ordered_overlaps.end()) {
                    break;
                }
            }
            
#ifdef debug_multipath_alignment
            cerr << "performing an overlap split onto " << get<2>(*iter) << " from " << get<3>(*iter);
            auto it = iter;
            ++it;
            for (; it != iter_range_end; ++it) {
                cerr << ", " << get<3>(*it);
            }
            cerr << " of read length " << get<0>(*iter) << " and graph length " << get<1>(*iter);
            if (path_node_provenance) {
                cerr << ", provenances " << path_node_provenance->at(get<2>(*iter)) << " and " << path_node_provenance->at(get<3>(*iter));
            }
            cerr << endl;
            
#endif
            
            
            PathNode* onto_node = &path_nodes[get<2>(*iter)];
            
#ifdef debug_multipath_alignment
            cerr << "before splitting:" << endl;
            cerr << "onto node:" << endl << "\t";
            for (auto node_iter = onto_node->begin; node_iter != onto_node->end; node_iter++) {
                cerr << *node_iter;
            }
            cerr << endl << "\t" << debug_string(onto_node->path) << endl;
#endif
            
            // TODO: there should be a way to do this in a single pass over mappings and edits
            // rather than traversing the whole mapping twice
            
            // store the full path and remove it from the node
            path_t full_path = std::move(onto_node->path);
            onto_node->path.clear_mapping();
            
            // keep track of how the read sequence should get split up
            size_t prefix_to_length = 0;
            
            // add mappings from the path until reaching the overlap point
            int64_t to_remaining = get<0>(*iter);
            int64_t from_remaining = get<1>(*iter);
            int64_t mapping_idx = 0;
            int64_t mapping_from_len = mapping_from_length(full_path.mapping(mapping_idx));
            int64_t mapping_to_len = mapping_to_length(full_path.mapping(mapping_idx));
            while (to_remaining >= mapping_to_len && from_remaining >= mapping_from_len) {
                *onto_node->path.add_mapping() = full_path.mapping(mapping_idx);
                prefix_to_length += mapping_to_len;
                to_remaining -= mapping_to_len;
                from_remaining -= mapping_from_len;
                mapping_idx++;
                if (mapping_idx == full_path.mapping_size()) {
                    break;
                }
                
                mapping_from_len = mapping_from_length(full_path.mapping(mapping_idx));
                mapping_to_len = mapping_to_length(full_path.mapping(mapping_idx));
            }
            
            if (mapping_idx == full_path.mapping_size() && !to_remaining && !from_remaining) {
                // TODO: isn't this case covered by taking the entire range of splits at the same place?
                
                // the overlap covered the path, so connect it to the onto node's successors
                // rather than splitting it into two nodes
                while (iter != iter_range_end) {
                    for (const pair<size_t, size_t> edge : onto_node->edges) {
                        path_nodes.at(get<3>(*iter)).edges.emplace_back(edge.first, edge.second + get<4>(*iter));
                    }
                    iter++;
                }
            }
            else {
                // the overlap was in the middle of the path, so split the onto node into a
                // prefix and suffix
                
                // make a new node to hold the suffix of the path
                size_t suffix_idx = path_nodes.size();
                path_nodes.emplace_back();
                if (path_node_provenance) {
                    path_node_provenance->emplace_back((*path_node_provenance)[get<2>(*iter)]);
                }
                PathNode& suffix_node = path_nodes.back();
                
                // get the pointer from the onto node back in case the vector reallocated
                onto_node = &path_nodes.at(get<2>(*iter));
                
                // transfer the outgoing edges onto the new node
                suffix_node.edges = std::move(onto_node->edges);
                
                // clear the old edges and add a single edge to the suffix
                onto_node->edges.clear();
                onto_node->edges.emplace_back(suffix_idx, 0);
                
                // keep track of the relationship of suffix nodes to original nodes
                if (!node_with_suffix.count(get<2>(*iter))) {
                    // since we take longest overlaps first, only the first split onto a node will change
                    // which node contains the suffix of the original node
                    node_with_suffix[get<2>(*iter)] = suffix_idx;
                }
                
                if (to_remaining || from_remaining) {
                    // the overlap point is in the middle of a node, need to split a mapping
                    
                    const path_mapping_t& split_mapping = full_path.mapping(mapping_idx);
                    
                    // add the prefix of the mapping to the original node
                    path_mapping_t* prefix_split = onto_node->path.add_mapping();
                    prefix_split->mutable_position()->set_node_id(split_mapping.position().node_id());
                    prefix_split->mutable_position()->set_offset(split_mapping.position().offset());
                    
                    // add the suffix of the mapping to the new node
                    path_mapping_t* suffix_split = suffix_node.path.add_mapping();
                    suffix_split->mutable_position()->set_node_id(split_mapping.position().node_id());
                    suffix_split->mutable_position()->set_offset(split_mapping.position().offset() + from_remaining);
                                        
                    // add the edits up to the point where the split occurs
                    size_t edit_idx = 0;
                    int64_t mapping_to_remaining = to_remaining;
                    int64_t mapping_from_remaining = from_remaining;
                    for (; edit_idx < split_mapping.edit_size()
                         && mapping_from_remaining >= split_mapping.edit(edit_idx).from_length()
                         && mapping_to_remaining >= split_mapping.edit(edit_idx).to_length(); edit_idx++) {
                        const edit_t& split_edit = split_mapping.edit(edit_idx);
                        mapping_from_remaining -= split_edit.from_length();
                        mapping_to_remaining -= split_edit.to_length();
                        prefix_to_length += split_edit.to_length();
                        *prefix_split->add_edit() = split_edit;
                    }
                    
                    // do we need to split in the middle of an edit?
                    if (mapping_from_remaining || mapping_to_remaining) {
                        
                        const edit_t& split_edit = split_mapping.edit(edit_idx);
                        // add an edit for either side of the split
                        edit_t* prefix_split_edit = prefix_split->add_edit();
                        prefix_split_edit->set_from_length(mapping_from_remaining);
                        prefix_split_edit->set_to_length(mapping_to_remaining);
                        
                        prefix_to_length += prefix_split_edit->to_length();
                        
                        edit_t* suffix_split_edit = suffix_split->add_edit();
                        suffix_split_edit->set_from_length(split_edit.from_length() - mapping_from_remaining);
                        suffix_split_edit->set_to_length(split_edit.to_length() - mapping_to_remaining);
                        
                        if (!split_edit.sequence().empty()) {
                            prefix_split_edit->set_sequence(split_edit.sequence().substr(0, mapping_to_remaining));
                            suffix_split_edit->set_sequence(split_edit.sequence().substr(mapping_to_remaining, string::npos));
                        }
                        
                        edit_idx++;
                    }
                    
                    // add the remaining edits after the split
                    for (; edit_idx < split_mapping.edit_size(); edit_idx++) {
                        *suffix_split->add_edit() = split_mapping.edit(edit_idx);
                    }
                    
                    mapping_idx++;
                }
                
                // add the remaining mappings to the suffix node
                for (; mapping_idx < full_path.mapping_size(); mapping_idx++) {
                    *suffix_node.path.add_mapping() = full_path.mapping(mapping_idx);
                }
                
                // divide up the read interval
                suffix_node.end = onto_node->end;
                suffix_node.begin = onto_node->begin + prefix_to_length;
                onto_node->end = suffix_node.begin;
                
#ifdef debug_multipath_alignment
                cerr << "after splitting:" << endl;
                cerr << "onto node:" << endl << "\t";
                for (auto node_iter = onto_node->begin; node_iter != onto_node->end; node_iter++) {
                    cerr << *node_iter;
                }
                cerr << endl << "\t" << debug_string(onto_node->path) << endl;
                cerr << "suffix node:" << endl << "\t";
                for (auto node_iter = suffix_node.begin; node_iter != suffix_node.end; node_iter++) {
                    cerr << *node_iter;
                }
                cerr << endl << "\t" << debug_string(suffix_node.path) << endl;
#endif
                
                while (iter != iter_range_end) {
                    // index of the node that contains the end of the original node we recorded the overlap from
                    size_t splitting_idx = node_with_suffix.count(get<3>(*iter)) ? node_with_suffix[get<3>(*iter)] : get<3>(*iter);
                    
#ifdef debug_multipath_alignment
                    cerr << "adding an overlap edge from node " << splitting_idx << " at distance " << get<4>(*iter) << endl;
                    cerr << "\t";
                    for (auto node_iter = path_nodes.at(splitting_idx).begin; node_iter != path_nodes.at(splitting_idx).end; node_iter++) {
                        cerr << *node_iter;
                    }
                    cerr << endl;
#endif
                    
                    // get the next node that overlaps onto the other node at this index and add the overlap edge
                    path_nodes.at(splitting_idx).edges.emplace_back(suffix_idx, get<4>(*iter));
                    
                    iter++;
                }
            }
        }
        
        // Go to the state where we know the reachability edges exist.
        has_reachability_edges = true;
        
#ifdef debug_multipath_alignment
        cerr << "final graph after adding reachability edges:" << endl;
        for (size_t i = 0; i < path_nodes.size(); i++) {
            PathNode& path_node = path_nodes.at(i);
            cerr << i;
            if (path_node_provenance) {
                cerr << " (hit " << path_node_provenance->at(i) << ")";
            }
            cerr << " " << debug_string(path_node.path) << " " << string(path_node.begin, path_node.end) << endl;
            cerr << "\t";
            for (auto edge : path_node.edges) {
                cerr << "(to:" << edge.first << ", graph dist:" << edge.second << ", read dist: " << (path_nodes.at(edge.first).begin - path_node.end) << ") ";
            }
            cerr << endl;
        }
#endif
        
    }
    
    void MultipathAlignmentGraph::clear_reachability_edges() {
        // Don't let people clear the edges if they are clear already.
        // That suggests that someone has gotten confused over whether they should exist or not.
        assert(has_reachability_edges);
    
        for (auto& node : path_nodes) {
            // Just clear all the edges from each node.
            // add_reachability_edges can rebuild them all.
            node.edges.clear();
        }
        
        // Go to the state where reachability edges don't exist
        has_reachability_edges = false;
        
    }
    
    // Kahn's algorithm
    void MultipathAlignmentGraph::topological_sort(vector<size_t>& order_out) {
        // Can only sort if edges are present.
        assert(has_reachability_edges);
        
        order_out.resize(path_nodes.size());
       
        vector<size_t> in_degree(path_nodes.size());
        
        for (const PathNode& path_node : path_nodes) {
            for (const pair<size_t, size_t>& edge : path_node.edges) {
                in_degree[edge.first]++;
            }
        }
        
        list<size_t> source_queue;
        for (size_t i = 0; i < path_nodes.size(); i++) {
            if (in_degree[i] == 0) {
                source_queue.push_back(i);
            }
        }
        
        size_t next = 0;
        while (!source_queue.empty()) {
            size_t src = source_queue.front();
            source_queue.pop_front();
            
            for (const pair<size_t, size_t>& edge : path_nodes.at(src).edges) {
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
        vector<vector<pair<size_t, size_t>>> reverse_graph(path_nodes.size());
        for (size_t i = 0; i < path_nodes.size(); i++) {
            for (const pair<size_t, size_t>& edge : path_nodes.at(i).edges) {
                reverse_graph[edge.first].emplace_back(i, edge.second);
            }
        }
        for (PathNode& path_node : path_nodes) {
            size_t out_degree = path_node.edges.size();
            path_node.edges.clear();
            path_node.edges.reserve(out_degree);
        }
        for (size_t i : order) {
            for (const pair<size_t, size_t>& edge : reverse_graph[i]) {
                path_nodes.at(edge.first).edges.emplace_back(i, edge.second);
            }
        }
    }
    
    void MultipathAlignmentGraph::remove_transitive_edges(const vector<size_t>& topological_order) {
        // We can only remove edges when the edges are present
        assert(has_reachability_edges);
        
        // algorithm assumes edges are also sorted in topological order, which guarantees that we will
        // traverse a path that reveals an edge as transitive before actually traversing the transitive edge
        reorder_adjacency_lists(topological_order);
        
        for (size_t i : topological_order) {
            vector<pair<size_t, size_t>>& edges = path_nodes[i].edges;
            
            // if there is only one edge out of a node, that edge can never be transitive
            // (this optimization covers most cases)
            if (edges.size() <= 1) {
                continue;
            }
            
            vector<bool> keep(edges.size(), true);
            unordered_set<size_t> traversed;
            
            for (size_t j = 0; j < edges.size(); j++) {
                const pair<size_t, size_t>& edge = edges[j];
                if (traversed.count(edge.first) && edge.second != 0 &&
                    path_nodes[i].end != path_nodes[edge.first].begin) {
                    // we can reach the target of this edge by another path, so it is transitive
                    // and the path nodes don't abut on either the read or graph
                    keep[j] = false;
                    continue;
                }
                
                // DFS to mark all reachable nodes from this edge
                vector<size_t> stack{edge.first};
                traversed.insert(edge.first);
                while (!stack.empty()) {
                    size_t idx = stack.back();
                    stack.pop_back();
                    for (const pair<size_t, size_t>& edge_from : path_nodes.at(idx).edges) {
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
        
        
#ifdef debug_multipath_alignment
        cerr << "removed transitive edges, topology is:" << endl;
        for (size_t i = 0; i < path_nodes.size(); i++) {
            cerr << "node " << i << ", " << debug_string(path_nodes.at(i).path.mapping(0).position()) << " ";
            for (auto iter = path_nodes.at(i).begin; iter != path_nodes.at(i).end; iter++) {
                cerr << *iter;
            }
            cerr << endl;
            for (pair<size_t, size_t> edge : path_nodes.at(i).edges) {
                cerr << "\tto " << edge.first << ", dist " << edge.second << endl;
            }
        }
#endif
    }
    
    MultipathAlignmentGraph::~MultipathAlignmentGraph() {
        
    }
    
    void MultipathAlignmentGraph::prune_to_high_scoring_paths(const Alignment& alignment, const GSSWAligner* aligner,
                                                              double max_suboptimal_score_ratio, const vector<size_t>& topological_order,
                                                              vector<size_t>& path_node_provenance) {
        
        // Can only prune when edges exist.
        assert(has_reachability_edges);
        
        if (path_nodes.empty()) {
            return;
        }
        
        unordered_map<pair<size_t, size_t>, int32_t> edge_weights;
        
        vector<int32_t> node_weights(path_nodes.size());
        
        // TODO: is the lower bound too strict?
        
        // compute the weight of edges and node matches
        for (size_t i = 0; i < path_nodes.size(); i++) {
            PathNode& from_node = path_nodes.at(i);
            node_weights[i] = (aligner->score_exact_match(from_node.begin, from_node.end,
                                                          alignment.quality().begin() + (from_node.begin - alignment.sequence().begin()))
                               + (from_node.begin == alignment.sequence().begin() ? aligner->score_full_length_bonus(true, alignment) : 0)
                               + (from_node.end == alignment.sequence().end() ? aligner->score_full_length_bonus(false, alignment) : 0));
                        
            for (const pair<size_t, size_t>& edge : from_node.edges) {
                PathNode& to_node = path_nodes.at(edge.first);
                
                int64_t graph_dist = edge.second;
                int64_t read_dist = to_node.begin - from_node.end;
                
                if (read_dist > graph_dist) {
                    // the read length in between the MEMs is longer than the distance, suggesting a read insert
                    // and potentially another mismatch on the other end
                    int64_t gap_length = read_dist - graph_dist;
                    edge_weights[make_pair(i, edge.first)] = (-(gap_length - 1) * aligner->gap_extension - aligner->gap_open
                                                              - (graph_dist > 0) * aligner->mismatch);
                }
                else if (read_dist < graph_dist) {
                    // the read length in between the MEMs is shorter than the distance, suggesting a read deletion
                    // and potentially another mismatch on the other end
                    int64_t gap_length = graph_dist - read_dist;
                    edge_weights[make_pair(i, edge.first)] = (-(gap_length - 1) * aligner->gap_extension - aligner->gap_open
                                                              - (read_dist > 0) * aligner->mismatch);
                }
                else {
                    // the read length in between the MEMs is the same as the distance, suggesting a pure mismatch
                    edge_weights[make_pair(i, edge.first)] = -((graph_dist > 0) + (graph_dist > 1)) * aligner->mismatch;
                }
            }
        }
        
        vector<int32_t> forward_scores = node_weights;
        vector<int32_t> backward_scores = node_weights;
        
        // forward DP
        for (int64_t i = 0; i < topological_order.size(); i++) {
            size_t idx = topological_order[i];
            int32_t from_score = forward_scores[idx];
            for (const pair<size_t, size_t>& edge : path_nodes.at(idx).edges) {
                forward_scores[edge.first] = std::max(forward_scores[edge.first],
                                                      node_weights[edge.first] + from_score + edge_weights[make_pair(idx, edge.first)]);
            }
        }
        
        // backward DP
        for (int64_t i = topological_order.size() - 1; i >= 0; i--) {
            size_t idx = topological_order[i];
            int32_t score_here = node_weights[idx];
            for (const pair<size_t, size_t>& edge : path_nodes.at(idx).edges) {
                backward_scores[idx] = std::max(backward_scores[idx],
                                                score_here + backward_scores[edge.first] + edge_weights[make_pair(idx, edge.first)]);
            }
        }
        
        // compute the minimum score we will require of a node or edge
        int32_t min_path_score = *std::max_element(forward_scores.begin(), forward_scores.end()) / max_suboptimal_score_ratio;
        
        // use forward-backward to find nodes on some path with a score above the minimum
        vector<size_t> removed_in_prefix(path_nodes.size() + 1, 0);
        for (size_t i = 0; i < path_nodes.size(); i++) {
            removed_in_prefix[i + 1] = (removed_in_prefix[i]
                                        + int(forward_scores[i] + backward_scores[i] - node_weights[i] < min_path_score));
        }
        
        // remove any nodes that failed to meet the threshold
        for (size_t i = 0; i < path_nodes.size(); ++i) {
            if (removed_in_prefix[i] == removed_in_prefix[i + 1]) {
                // we're keeping this node
                auto& path_node = path_nodes[i];
                
                // remove its edges that aren't on a sufficiently high-scoring path too
                size_t edges_removed = 0;
                for (size_t j = 0; j < path_node.edges.size(); ++j) {
                    auto& edge = path_node.edges[j];
                    if (forward_scores[i] + backward_scores[edge.first] + edge_weights[make_pair(i, edge.first)] < min_path_score) {
                        ++edges_removed;
                    }
                    else {
                        path_node.edges[j - edges_removed] = make_pair(edge.first - removed_in_prefix[edge.first], edge.second);
                    }
                }
                path_node.edges.resize(path_node.edges.size() - edges_removed);
                if (removed_in_prefix[i]) {
                    path_nodes[i - removed_in_prefix[i]] = move(path_node);
                    path_node_provenance[i - removed_in_prefix[i]] = path_node_provenance[i];
                }
            }
        }
        
        path_nodes.resize(path_nodes.size() - removed_in_prefix.back());
        path_node_provenance.resize(path_nodes.size());
        
#ifdef debug_multipath_alignment
        cerr << "pruned to high scoring paths, topology is:" << endl;
        for (size_t i = 0; i < path_nodes.size(); i++) {
            cerr << "node " << i << " (hit " << path_node_provenance[i] << "), " << debug_string(path_nodes.at(i).path.mapping(0).position()) << " ";
            for (auto iter = path_nodes.at(i).begin; iter != path_nodes.at(i).end; iter++) {
                cerr << *iter;
            }
            cerr << endl;
            for (pair<size_t, size_t> edge : path_nodes.at(i).edges) {
                cerr << "\tto " << edge.first << ", dist " << edge.second << endl;
            }
        }
#endif
    }
    
void MultipathAlignmentGraph::align(const Alignment& alignment, const HandleGraph& align_graph, const GSSWAligner* aligner,
                                    bool score_anchors_as_matches, size_t max_alt_alns, bool dynamic_alt_alns, size_t max_gap,
                                    double pessimistic_tail_gap_multiplier, bool simplify_topologies, size_t unmergeable_len,
                                    size_t band_padding, multipath_alignment_t& multipath_aln_out, SnarlManager* cutting_snarls,
                                    MinimumDistanceIndex* dist_index, const function<pair<id_t, bool>(id_t)>* project,
                                    bool allow_negative_scores) {
        
        // don't dynamically choose band padding, shim constant value into a function type
        function<size_t(const Alignment&,const HandleGraph&)> constant_padding = [&](const Alignment& seq, const HandleGraph& graph) {
            return band_padding;
        };
        align(alignment,
              align_graph,
              aligner,
              score_anchors_as_matches,
              max_alt_alns,
              dynamic_alt_alns,
              max_gap,
              pessimistic_tail_gap_multiplier,
              simplify_topologies,
              unmergeable_len,
              constant_padding,
              multipath_aln_out,
              cutting_snarls,
              dist_index,
              project,
              allow_negative_scores);
    }

    pair<path_t, int32_t> MultipathAlignmentGraph::zip_alignments(vector<pair<path_t, int32_t>>& alt_alns, bool from_left,
                                                                  const Alignment& alignment, const HandleGraph& align_graph,
                                                                  string::const_iterator begin, const GSSWAligner* aligner) {
#ifdef debug_multipath_alignment
        cerr << "attempting to zip alignments from " << (from_left ? "left" : "right") << endl;
        for (auto& alt_aln : alt_alns) {
            cerr << debug_string(alt_aln.first) << endl;
        }
#endif
        
        pair<path_t, int32_t> return_val;
        
        int64_t i, j, incr;
        auto& path_1 = alt_alns.front().first;
        if (from_left) {
            incr = 1;
            i = 0;
            j = 0;
        }
        else {
            incr = -1;
            i = path_1.mapping_size() - 1;
            if (i >= 0) {
                j = path_1.mapping(i).edit_size() - 1;
            }
            else {
                j = 0;
            }
        }
        // past-the-last index of the longest match that concludes in a match/mismatvch
        int64_t last_aligned_i, last_aligned_j;
        
        // walk the first path to check for full prefix/suffix matches
        bool found_mismatch = false, in_indel = false;
        while (i >= 0 && i < path_1.mapping_size()) {
            // check to make sure the positions of all of the alts match at this mapping
            const auto& mapping_1 = path_1.mapping(i);
            // start from the appropriate side of the mapping and check for matches of indels
            if (from_left) {
                j = 0;
            }
            else {
                j = mapping_1.edit_size() - 1;
            }
            const auto& pos_1 = mapping_1.position();
            auto offset_1 = from_left ? pos_1.offset() : pos_1.offset() + mapping_from_length(mapping_1);
            for (size_t k = 1; k < alt_alns.size(); ++k) {
                //cerr << "checking for matching position on aln " << k << endl;
                const auto& path_2 = alt_alns[k].first;
                int64_t i2;
                if (from_left) {
                    i2 = i;
                }
                else {
                    i2 = path_2.mapping_size() - path_1.mapping_size() + i;
                }
                if (i2 < 0 || i2 >= path_2.mapping_size()) {
                    //cerr << "no corresponding mapping in alt aln " << k << endl;
                    found_mismatch = true;
                    break;
                }
                const auto& mapping_2 = path_2.mapping(i2);
                const auto& pos_2 = mapping_2.position();
                if (pos_1.node_id() != pos_2.node_id()
                    || pos_2.is_reverse() != pos_2.is_reverse()
                    || (from_left ? pos_2.offset() : pos_2.offset() + mapping_from_length(mapping_2)) != offset_1) {
                    //cerr << "positions mismatch" << endl;
                    found_mismatch = true;
                    break;
                }
            }
            if (found_mismatch) {
                break;
            }
            while(j >= 0 && j < path_1.mapping(i).edit_size()) {
                const auto& edit_1 = mapping_1.edit(j);
                //cerr << "checking for matches to edit " << i << " " << j << ": " << debug_string(edit_1) << endl;
                if (!in_indel) {
                    // we've matched up to here and the most recent match wasn't an indel
                    last_aligned_i = i;
                    last_aligned_j = j;
                }
                // check whether all of the other paths match the first path at the next edit
                bool found_unmatchable = false;
                for (size_t k = 1; k < alt_alns.size() && !found_mismatch; ++k) {
                    auto& path_2 = alt_alns[k].first;
                    // find the corresponding mapping index
                    int64_t i2;
                    if (from_left) {
                        i2 = i;
                    }
                    else {
                        i2 = path_2.mapping_size() - path_1.mapping_size() + i;
                    }
                    
                    if (i2 < 0 || i2 >= path_2.mapping_size()) {
                        //cerr << "no corresponding mapping in alt aln " << k << endl;
                        found_mismatch = true;
                        break;
                    }
                    const auto& mapping_2 = path_2.mapping(i2);
                    // find the corresponding edit index
                    int64_t j2;
                    if (from_left) {
                        j2 = j;
                    }
                    else {
                        j2 = mapping_2.edit_size() - mapping_1.edit_size() + j;
                    }
                    if (j2 < 0 || j2 >= mapping_2.edit_size()) {
                        //cerr << "no corresponding edit in alt aln " << k << endl;
                        found_mismatch = true;
                        break;
                    }
                    const auto& edit_2 = mapping_2.edit(j2);
                    //cerr << "comparing to alt aln " << k << " edit " << i2 << " " << j2 << ": " << debug_string(edit_2) << endl;
                    found_mismatch = (edit_1 != edit_2);
                    if ((from_left && j + 1 == mapping_1.edit_size()) || (!from_left && j == 0)) {
                        // check if there will be edits in this mapping that we won't hit by using the indexes
                        // of edits on path 1
                        found_unmatchable = found_unmatchable || (mapping_1.edit_size() != mapping_2.edit_size());
                    }
                }
                if (found_mismatch) {
                    break;
                }
                else {
                    j += incr;
                    in_indel = (edit_1.from_length() == 0 || edit_1.to_length() == 0);
                    if (found_unmatchable) {
                        found_mismatch = true;
                        break;
                    }
                }
            }
            if (found_mismatch) {
                break;
            }
            else {
                i += incr;
                if (from_left || i < 0) {
                    j = 0;
                }
                else {
                    j = path_1.mapping(i).edit_size() - 1;
                }
            }
        }
        
        if (in_indel) {
            //cerr << "last match is in an indel" << endl;
            // the match concluded in an insertion or deletion, we need to make sure it doesn't
            // continue onto the next edit
            bool all_nexts_aligned = true;
            if (i >= 0 && i < path_1.mapping_size()) {
                // we didn't match the entire first path, so we need to actually check
                const auto& edit = path_1.mapping(i).edit(j);
                all_nexts_aligned = (edit.from_length() != 0 && edit.to_length() != 0);
            }
            for (size_t k = 1; k < alt_alns.size() && all_nexts_aligned; ++k) {
                auto& path_2 = alt_alns[k].first;
                int64_t i2 = path_2.mapping_size() - path_1.mapping_size() + i;
                if (i2 >= 0 && i2 < path_2.mapping_size()) {
                    // we didn't run through the full path
                    const auto& mapping_2 = path_2.mapping(i2);
                    int64_t j2 = mapping_2.edit_size() - path_1.mapping(i).edit_size() + j;
                    if (j2 >= 0 && j2 < mapping_2.edit_size()) {
                        // we didn't run through the full mapping
                        const auto& edit = mapping_2.edit(j2);
                        all_nexts_aligned = (edit.from_length() != 0 && edit.to_length() != 0);
                    }
                }
            }
            if (!all_nexts_aligned) {
                // we need to backtrack to the last aligned base
                i = last_aligned_i;
                j = last_aligned_j;
            }
        }
        
        
        if ((from_left && (i != 0 || j != 0)) ||
             (!from_left && (i != path_1.mapping_size() - 1 || j != path_1.mapping().back().edit_size() - 1))) {
            //cerr << "matched up to " << i << " " << j << endl;
            // we matched part of the path across all alternate alignments, we can zip it up
            
            // delete the corresponding parts from all the other paths
            for (int64_t k = 1; k < alt_alns.size(); ++k) {
                //cerr << "deleting from aln " << k << endl;
                auto& path_2 = alt_alns[k].first;
                if (from_left) {
                    // we have to remove part of the path from the beginning of the existing paths
                    int64_t num_mappings, num_edits;
                    if (j == 0) {
                        if (path_1.mapping(i - 1).edit_size() == path_2.mapping(i - 1).edit_size()) {
                            // this is the past-the-last on the mapping
                            num_mappings = i;
                            num_edits = 0;
                        }
                        else {
                            // the full mapping on path 1 didn't use up the final mapping on path 2
                            num_mappings = i - 1;
                            num_edits = path_1.mapping(i - 1).edit_size();
                        }
                    }
                    else {
                        if (i < path_2.mapping_size() && j == path_2.mapping(i).edit_size()) {
                            // past-the-last beyond the end of this mapping
                            num_mappings = i + 1;
                            num_edits = 0;
                        }
                        else {
                            // past-the-last in the middle of this mapping
                            num_mappings = i;
                            num_edits = j;
                        }
                    }
                    //cerr << "deleting " << num_mappings << " mappings and " << num_edits << " from the left" << endl;
                    
                    if (num_mappings < path_2.mapping_size() && num_edits != 0) {
                        // we have to remove part of the mapping
                        auto mapping = path_2.mutable_mapping(i);
                        int64_t from_len = 0;
                        for (int64_t l = 0; l < num_edits; ++l) {
                            from_len += mapping->edit(l).from_length();
                        }
                        mapping->mutable_edit()->erase(mapping->mutable_edit()->begin(),
                                                       mapping->mutable_edit()->begin() + num_edits);
                        mapping->mutable_position()->set_offset(mapping->position().offset() + from_len);
                    }
                    // remove any full mappings that are shared
                    path_2.mutable_mapping()->erase(path_2.mutable_mapping()->begin(),
                                                    path_2.mutable_mapping()->begin() + num_mappings);
                }
                else {
                    // we remove from the back, more complicated because of needing to translate indexes
                    
                    int64_t i2 = path_2.mapping_size() - path_1.mapping_size() + i;
                    int64_t num_mappings, num_edits;
                    if (i2 < 0) {
                        // we delete all of path 2
                        num_mappings = path_2.mapping_size();
                        num_edits = 0;
                    }
                    else if (i < 0) {
                        // we ran through path 1 but not path 2
                        if (path_1.mapping().front().edit_size() == path_2.mapping(i2 + 1).edit_size()) {
                            // we ran through the full mapping on path 2
                            num_mappings = path_1.mapping_size();
                            num_edits = 0;
                        }
                        else {
                            // we didn't run through the entire mapping on path 2
                            num_mappings = path_1.mapping_size() - 1;
                            num_edits = path_1.mapping().front().edit_size();
                        }
                    }
                    else {
                        // both of the ending mappings can be indexed
                        const auto& mapping_1 = path_1.mapping(i);
                        if (j == mapping_1.edit_size() - 1) {
                            // mismatch at the beginning of a mapping on path 1
                            if (path_1.mapping(i + 1).edit_size() == path_2.mapping(i2 + 1).edit_size()) {
                                // the previous mappings on both paths end at the same edit
                                num_mappings = path_1.mapping_size() - i - 1;
                                num_edits = 0;
                            }
                            else {
                                // we didn't run through the entire mapping on path 2
                                num_mappings = path_1.mapping_size() - i - 2;
                                num_edits = path_1.mapping(i + 1).edit_size();
                            }
                        }
                        else {
                            // the mismatch happened after at least one edit in this mapping on path 1
                            int64_t j2 = path_2.mapping(i2).edit_size() - mapping_1.edit_size() + j;
                            if (j2 < 0) {
                                // past-the-last beyond the end of this mapping
                                num_mappings = path_1.mapping_size() - i;
                                num_edits = 0;
                            }
                            else {
                                // past-the-last in the middle of this mapping
                                num_mappings = path_1.mapping_size() - i - 1;
                                num_edits = mapping_1.edit_size() - j - 1;
                            }
                        }
                    }
                    //cerr << "deleting " << num_mappings << " mappings and " << num_edits << " from the right" << endl;
                    
                    // delete the number of mappings and edits we decided
                    path_2.mutable_mapping()->resize(path_2.mapping_size() - num_mappings);
                    if (num_edits) {
                        auto& mapping_2 = path_2.mutable_mapping()->back();
                        mapping_2.mutable_edit()->resize(mapping_2.edit_size() - num_edits);
                    }
                }
            }
            
            // figure out how much of the primary to move over to the zipped alignment
            int64_t k, num_mappings, num_edits;
            if (from_left) {
                k = 0;
                num_mappings = i;
                num_edits = j;
                if (i < path_1.mapping_size() && j >= path_1.mapping(i).edit_size()) {
                    ++num_mappings;
                    num_edits = 0;
                }
            }
            else {
                k = path_1.mapping_size() - 1;
                num_mappings = k - i;
                num_edits = 0;
                if (j < 0) {
                    ++num_mappings;
                }
                else if (i >= 0) {
                    num_edits = path_1.mapping(i).edit_size() - j - 1;
                }
            }
            //cerr << "need to copy " << num_mappings << " full mappings and " << num_edits << " edits" << endl;
            // copy to the zipped alignment
            for (int64_t m = 0; m < num_mappings; k += incr, ++m) {
                //cerr << "copy mapping " << debug_string(*path_1.mutable_mapping(k)) << endl;
                *return_val.first.add_mapping() = move(*path_1.mutable_mapping(k));
            }
            if (num_edits != 0) {
                auto mapping = path_1.mutable_mapping(k);
                int64_t mapping_len = mapping_from_length(*mapping);
                auto pos = mapping->mutable_position();
                auto copy_mapping = return_val.first.add_mapping();
                auto copy_pos = copy_mapping->mutable_position();
                for (int64_t l = from_left ? 0 : mapping->edit_size() - 1, e = 0; e < num_edits; l += incr, ++e) {
                    auto edit = mapping->mutable_edit(l);
                    //cerr << "copy edit " << debug_string(*edit) << endl;
                    *copy_mapping->add_edit() = move(*edit);
                }
                if (from_left) {
                    mapping->mutable_edit()->erase(mapping->mutable_edit()->begin(),
                                                   mapping->mutable_edit()->begin() + num_edits);
                    *copy_pos = *pos;
                    pos->set_offset(pos->offset() + mapping_from_length(*copy_mapping));
                }
                else {
                    // the edits were added in reverse order, switch them back
                    reverse(copy_mapping->mutable_edit()->begin(), copy_mapping->mutable_edit()->end());
                    
                    mapping->mutable_edit()->resize(mapping->edit_size() - num_edits);
                    copy_pos->set_node_id(pos->node_id());
                    copy_pos->set_is_reverse(pos->is_reverse());
                    copy_pos->set_offset(pos->offset() + (mapping_len - mapping_from_length(*copy_mapping)));
                }
            }
            
            string::const_iterator alignment_begin;
            if (from_left) {
                path_1.mutable_mapping()->erase(path_1.mutable_mapping()->begin(),
                                                path_1.mutable_mapping()->begin() + num_mappings);
                alignment_begin = begin;
            }
            else {
                // the mappings were added in reverse order, switch them back
                reverse(return_val.first.mutable_mapping()->begin(), return_val.first.mutable_mapping()->end());
                
                path_1.mutable_mapping()->resize(path_1.mapping_size() - num_mappings);
                alignment_begin = begin + path_to_length(path_1);
            }
            
            // score the zipped alignment
            return_val.second = aligner->score_partial_alignment(alignment, align_graph, return_val.first, alignment_begin);
            
            // remove this part of the score from the alt alignments
            for (auto& alt_aln : alt_alns) {
                alt_aln.second -= return_val.second;
            }
#ifdef debug_multipath_alignment
            cerr << "successfully made zipped alignment with score " << return_val.second << endl;
            cerr << debug_string(return_val.first) << endl;
#endif
        }
        
        return return_val;
    }
    
    void MultipathAlignmentGraph::align(const Alignment& alignment, const HandleGraph& align_graph, const GSSWAligner* aligner,
                                        bool score_anchors_as_matches, size_t max_alt_alns, bool dynamic_alt_alns, size_t max_gap,
                                        double pessimistic_tail_gap_multiplier, bool simplify_topologies, size_t unmergeable_len,
                                        function<size_t(const Alignment&,const HandleGraph&)> band_padding_function,
                                        multipath_alignment_t& multipath_aln_out, SnarlManager* cutting_snarls,
                                        MinimumDistanceIndex* dist_index, const function<pair<id_t, bool>(id_t)>* project,
                                        bool allow_negative_scores) {
        
        // Can only align if edges are present.
        assert(has_reachability_edges);
        
        // transfer over data from alignment
        transfer_read_metadata(alignment, multipath_aln_out);
        
#ifdef debug_multipath_alignment
        cerr << "transferred over read information" << endl;
#endif
        // TODO: could I get away with moving the paths instead of copying them?
        
        // add a subpath for each of the exact match nodes
        if (score_anchors_as_matches) {
            for (int64_t j = 0; j < path_nodes.size(); j++) {
                PathNode& path_node = path_nodes.at(j);
                subpath_t* subpath = multipath_aln_out.add_subpath();
                *subpath->mutable_path() = path_node.path;
                int32_t match_score = aligner->score_exact_match(path_node.begin, path_node.end,
                                                                 alignment.quality().begin() + (path_node.begin - alignment.sequence().begin()));
                
                subpath->set_score(match_score
                                   + (path_node.begin == alignment.sequence().begin() ? aligner->score_full_length_bonus(true, alignment) : 0)
                                   + (path_node.end == alignment.sequence().end() ? aligner->score_full_length_bonus(false, alignment) : 0));
            }
        }
        else {
            for (size_t j = 0; j < path_nodes.size(); j++) {
                PathNode& path_node = path_nodes.at(j);
                subpath_t* subpath = multipath_aln_out.add_subpath();
                *subpath->mutable_path() = path_node.path;
                subpath->set_score(aligner->score_partial_alignment(alignment, align_graph, path_node.path, path_node.begin));
            }
        }
        
        // the indexes of subpaths that we will not allow to merge at non-branching paths
        unordered_set<size_t> prohibited_merges;
        
        // we use this function to remove alignments that follow the same path, keeping only the highest scoring one
        auto deduplicate_alt_alns = [](vector<Alignment>& alt_alns, bool leftward, bool rightward) -> vector<pair<path_t, int32_t>> {
            // init the deduplicated vector in STL types with move operators
            vector<pair<path_t, int32_t>> deduplicated(alt_alns.size());
            for (size_t i = 0; i < alt_alns.size(); ++i) {
                const auto& aln = alt_alns[i];
                deduplicated[i].second = aln.score();
                from_proto_path(aln.path(), deduplicated[i].first);
            }
            
            // we use stable sort to keep the original score-ordering among alignments
            // that take the same path, which is descending by score
            stable_sort(deduplicated.begin(), deduplicated.end(),
                        [&](const pair<path_t, int32_t>& aln_1, const pair<path_t, int32_t>& aln_2) {
                const auto& path_1 = aln_1.first, path_2 = aln_2.first;
                int64_t i, j, incr;
                if (leftward) {
                    i = path_1.mapping_size() - 1;
                    j = path_2.mapping_size() - 1;
                    incr = -1;
                }
                else {
                    i = 0;
                    j = 0;
                    incr = 1;
                }
                bool is_less = false;
                for (; i >= 0 && j >= 0 && i < path_1.mapping_size() && j < path_2.mapping_size(); i += incr, j += incr) {
                    const auto& pos_1 = path_1.mapping(i).position(), pos_2 = path_2.mapping(j).position();
                    if (pos_1.node_id() < pos_2.node_id() ||
                        (pos_1.node_id() == pos_2.node_id() && pos_1.is_reverse() < pos_2.is_reverse())) {
                        is_less = true;
                        break;
                    }
                    else if (pos_1.node_id() > pos_2.node_id() ||
                             (pos_1.node_id() == pos_2.node_id() && pos_1.is_reverse() > pos_2.is_reverse())) {
                        break;
                    }
                }
                // supersequences earlier in the case of paths of different lengths
                return (is_less || ((i < path_1.mapping_size() && i >= 0) && (j == path_2.mapping_size() || j < 0)));
            });
            
            // move alignments that have the same path to the end of the vector, keeping the
            // first one (which has the highest score)
            auto new_end = unique(deduplicated.begin(), deduplicated.end(),
                                  [&](const pair<path_t, int32_t>& aln_1, const pair<path_t, int32_t>& aln_2) {
                const auto& path_1 = aln_1.first, path_2 = aln_2.first;
                int64_t i, j, incr;
                if (leftward) {
                    i = path_1.mapping_size() - 1;
                    j = path_2.mapping_size() - 1;
                    incr = -1;
                }
                else {
                    i = 0;
                    j = 0;
                    incr = 1;
                }
                // if this is a tail alignment, we allow paths of different lengths to be "equal"
                // if one is a prefix of the other and lower-scoring
                bool is_equal = (path_1.mapping_size() == path_2.mapping_size() || leftward || rightward);
                for (; i >= 0 && j >= 0 && i < path_1.mapping_size() && j < path_2.mapping_size() && is_equal; i += incr, j += incr) {
                    const auto& pos_1 = path_1.mapping(i).position(), pos_2 = path_2.mapping(j).position();
                    is_equal = (pos_1.node_id() == pos_2.node_id() && pos_1.is_reverse() == pos_2.is_reverse());
                }
                // TODO: there has to be a more succinct way to check this condition
                return (is_equal &&
                        (path_1.mapping_size() == path_2.mapping_size() ||
                         (path_1.mapping_size() > path_2.mapping_size() && aln_1.second > aln_2.second) ||
                         (path_1.mapping_size() < path_2.mapping_size() && aln_1.second < aln_2.second)));
            });
            
            // remove the duplicates at the end
            deduplicated.erase(new_end, deduplicated.end());
            return deduplicated;
        };
        
        // TODO: i wonder if it would be cleaner/more general to use branches rather than snarls
        // as the condition here...
        
        // add the alignment as subpath(s), possibly cutting it at snarl boundaries and marking
        // those cuts as prohibited to merge. the subpaths are created in order at the end of
        // the subpath vector and have edges between them
        auto add_and_permanently_cut = [&](pair<path_t, int32_t>& aln, string::const_iterator begin) {
            
#ifdef debug_multipath_alignment
            cerr << "assessing need to permanently cut path" << endl;
            cerr << debug_string(aln.first) << endl;
#endif
            
            // the intervals of the path that are inside snarls
            auto cut_segments = get_cut_segments(aln.first, cutting_snarls,dist_index, *project,
                                                 unmergeable_len);
            
            // collect the (internal) indexes that follow cuts
            vector<size_t> segment_boundaries;
            segment_boundaries.reserve(cut_segments.size() * 2);
            for (auto& cut_segment : cut_segments) {
                if (cut_segment.first != 0) {
                    segment_boundaries.push_back(cut_segment.first);
                }
                if (cut_segment.second != cut_segment.first && cut_segment.second != aln.first.mapping_size()) {
                    segment_boundaries.push_back(cut_segment.second);
                }
            }
            // make sure their in order (only ever not in order if there are nested snarls)
            if (!is_sorted(segment_boundaries.begin(), segment_boundaries.end())) {
                sort(segment_boundaries.begin(), segment_boundaries.end());
            }
            // don't allow edge-spanning deletions to be broken up (breaks dynamic programmability of scores)
            auto end = remove_if(segment_boundaries.begin(), segment_boundaries.end(), [&](size_t i) {
                return (aln.first.mapping(i - 1).edit().back().to_length() == 0 &&
                        aln.first.mapping(i).edit().front().to_length() == 0);
            });
            
#ifdef debug_multipath_alignment
            cerr << "need to cut path before mappings:" << endl;
            for (auto it = segment_boundaries.begin(); it != end; ++it) {
                cerr << "\t" << *it << endl;
            }
#endif
            
            if (segment_boundaries.begin() == end) {
                // we don't actually want to cut up this alignment at all
                auto subpath = multipath_aln_out.add_subpath();
                *subpath->mutable_path() = move(aln.first);
                subpath->set_score(aln.second);
            }
            else {
                // make a subpath for the first segment, but don't do anything with it
                size_t first_idx = multipath_aln_out.subpath_size();
                multipath_aln_out.add_subpath();
                
                for (auto it = segment_boundaries.begin(), next = segment_boundaries.begin() + 1; it != end; ++it, ++next) {
                    // add an unmergeable link from previous subpath
                    multipath_aln_out.mutable_subpath()->back().add_next(multipath_aln_out.subpath_size());
                    prohibited_merges.insert(multipath_aln_out.subpath_size() - 1);
                    // move path into the subpath
                    auto subpath = multipath_aln_out.add_subpath();
                    for (size_t i = *it, n = (next == end ? aln.first.mapping_size() : *next); i < n; ++i) {
                        *subpath->mutable_path()->add_mapping() = move(*aln.first.mutable_mapping(i));
                    }
                }
                // get the subpath here in case the vector reallocates
                auto first_subpath = multipath_aln_out.mutable_subpath(first_idx);
                aln.first.mutable_mapping()->resize(segment_boundaries.front());
                *first_subpath->mutable_path() = move(aln.first);
                
                // score the individual segments
                auto b = begin;
                for (size_t i = first_idx; i < multipath_aln_out.subpath_size(); ++i) {
                    auto subpath = multipath_aln_out.mutable_subpath(i);
                    subpath->set_score(aligner->score_partial_alignment(alignment, align_graph, subpath->path(), b));
                    b += path_to_length(subpath->path());
                }
                
#ifdef debug_multipath_alignment
                for (auto i = first_idx; i < multipath_aln_out.subpath_size(); ++i) {
                    cerr << "added permanently cut subpath " <<  i << endl;
                    cerr << debug_string(multipath_aln_out.subpath(i)) << endl;
                }
#endif
            }
        };
        
#ifdef debug_multipath_alignment
        cerr << "doing DP between MEMs" << endl;
#endif
        
        // perform alignment in the intervening sections
        for (int64_t j = 0; j < path_nodes.size(); j++) {
#ifdef debug_multipath_alignment
            cerr << "checking for intervening alignments from match node " << j << " with path " << debug_string(path_nodes.at(j).path) << " and sequence ";
            for (auto iter = path_nodes.at(j).begin; iter != path_nodes.at(j).end; iter++) {
                cerr << *iter;
            }
            cerr << endl;
#endif
            
            PathNode& src_path_node = path_nodes.at(j);
            subpath_t* src_subpath = multipath_aln_out.mutable_subpath(j);
            
            const path_t& path = src_subpath->path();
            const path_mapping_t& final_mapping = path.mapping(path.mapping_size() - 1);
            const position_t& final_mapping_position = final_mapping.position();
            // make a pos_t that points to the final base in the match
            pos_t src_pos = make_pos_t(final_mapping_position.node_id(),
                                       final_mapping_position.is_reverse(),
                                       final_mapping_position.offset() + mapping_from_length(final_mapping));
            
            // the longest gap that could be detected at this position in the read
            size_t src_max_gap = aligner->longest_detectable_gap(alignment, src_path_node.end);
            
            // This holds edges that we remove, because we couldn't actually get an alignment across them with a positive score.
            unordered_set<pair<size_t, size_t>> edges_for_removal;
            
            for (const pair<size_t, size_t>& edge : src_path_node.edges) {
                PathNode& dest_path_node = path_nodes.at(edge.first);
                pos_t dest_pos = make_pos_t(multipath_aln_out.subpath(edge.first).path().mapping(0).position());
                
#ifdef debug_multipath_alignment
                cerr << "forming intervening alignment for edge to node " << edge.first << endl;
#endif
                
                size_t intervening_length = dest_path_node.begin - src_path_node.end;

                // if negative score is allowed set maximum distance to the length between path nodes
                // otherwise set it to the maximum gap length possible while retaining a positive score 
                size_t max_dist = allow_negative_scores ?
                    edge.second :
                    intervening_length + min(min(src_max_gap, aligner->longest_detectable_gap(alignment, dest_path_node.begin)), max_gap);
                
#ifdef debug_multipath_alignment
                cerr << "read dist: " << intervening_length << ", graph dist " << edge.second << " source max gap: " << src_max_gap << ", dest max gap " << aligner->longest_detectable_gap(alignment, dest_path_node.begin) << ", max allowed gap " << max_gap << endl;
#endif
                
                // extract the graph between the matches
                bdsg::HashGraph connecting_graph;
                auto connect_trans = algorithms::extract_connecting_graph(&align_graph,      // DAG with split strands
                                                                          &connecting_graph, // graph to extract into
                                                                          max_dist,          // longest distance necessary
                                                                          src_pos,           // end of earlier match
                                                                          dest_pos,          // beginning of later match
                                                                          false);            // do not enforce max distance strictly
                                
                if (connecting_graph.get_node_count() == 0) {
                    // the MEMs weren't connectable with a positive score after all, mark the edge for removal
#ifdef debug_multipath_alignment
                    cerr << "Remove edge " << j << " -> " << edge.first << " because we got no nodes in the connecting graph "
                        << src_pos << " to " << dest_pos << endl;
#endif
                    edges_for_removal.insert(edge);
                    continue;
                }
                
                size_t num_alt_alns = dynamic_alt_alns ? min(max_alt_alns, handlealgs::count_walks(&connecting_graph)) :
                                                         max_alt_alns;
                
                
                // transfer the substring between the matches to a new alignment
                Alignment intervening_sequence;
                intervening_sequence.set_sequence(alignment.sequence().substr(src_path_node.end - alignment.sequence().begin(),
                                                                              dest_path_node.begin - src_path_node.end));
                if (!alignment.quality().empty()) {
                    intervening_sequence.set_quality(alignment.quality().substr(src_path_node.end - alignment.sequence().begin(),
                                                                                dest_path_node.begin - src_path_node.end));
                }
                
                // if we're doing dynamic alt alignments, possibly expand the number of tracebacks until we get an
                // alignment to every path or hit the hard max
                vector<pair<path_t, int32_t>> deduplicated;
                if (num_alt_alns > 0) {
                    
                    size_t num_alns_iter = num_alt_alns;
                    while (deduplicated.size() < num_alt_alns) {
                        
                        intervening_sequence.clear_path();
                        
#ifdef debug_multipath_alignment
                        cerr << "making " << num_alns_iter << " alignments of sequence " << intervening_sequence.sequence() << " to connecting graph" << endl;
                        connecting_graph.for_each_handle([&](const handle_t& handle) {
                            cerr << connecting_graph.get_id(handle) << " " << connecting_graph.get_sequence(handle) << endl;
                            connecting_graph.follow_edges(handle, true, [&](const handle_t& prev) {
                                cerr << "\t" << connecting_graph.get_id(prev) << " <-" << endl;
                            });
                            connecting_graph.follow_edges(handle, false, [&](const handle_t& next) {
                                cerr << "\t-> " << connecting_graph.get_id(next) << endl;
                            });
                        });
#endif
                        
                        vector<Alignment> alt_alignments;
                        aligner->align_global_banded_multi(intervening_sequence, alt_alignments, connecting_graph, num_alns_iter,
                                                           band_padding_function(intervening_sequence, connecting_graph), true);
                        
                        // remove alignments with the same path
                        deduplicated = deduplicate_alt_alns(alt_alignments, false, false);
                        
                        if (num_alns_iter >= max_alt_alns || !dynamic_alt_alns) {
                            // we don't want to try again even if we didn't find every path yet
                            break;
                        }
                        else {
                            // if we didn't find every path, we'll try again with this many tracebacks
                            num_alns_iter = min(max_alt_alns, num_alns_iter * 2);
                        }
                    }
                }
                
                bool added_direct_connection = false;
                for (auto& connecting_alignment : deduplicated) {
#ifdef debug_multipath_alignment
                    cerr << "translating connecting alignment: " << debug_string(connecting_alignment.first) << ", score " << connecting_alignment.second << endl;
#endif
                    
                    auto& aligned_path = connecting_alignment.first;
                    const auto& first_mapping = aligned_path.mapping(0);
                    const auto& last_mapping = aligned_path.mapping(aligned_path.mapping_size() - 1);
                    
                    bool add_first_mapping = mapping_from_length(first_mapping) != 0 || mapping_to_length(first_mapping) != 0;
                    bool add_last_mapping = mapping_from_length(last_mapping) != 0 || mapping_to_length(last_mapping) != 0;
                    
                    if (!(add_first_mapping || add_last_mapping) && aligned_path.mapping_size() <= 2) {
                        if (!added_direct_connection) {
                            // edge case where there is a simple split but other non-simple edges intersect the target
                            // at the same place (so it passes the previous deduplicating filter)
                            // it actually doesn't need an alignment, just a connecting edge
                            src_subpath->add_next(edge.first);
                            added_direct_connection = true;
                        }
                        continue;
                    }
                    
                    // create a subpath between the matches for this alignment
                    subpath_t* connecting_subpath = multipath_aln_out.add_subpath();
                    connecting_subpath->set_score(connecting_alignment.second);
                    
                    // get a pointer to the subpath again in case the vector reallocated
                    src_subpath = multipath_aln_out.mutable_subpath(j);
                    
                    // get rid of the ends if they are empty
                    if (!add_last_mapping) {
                        aligned_path.mutable_mapping()->pop_back();
                    }
                    if (!add_first_mapping && !aligned_path.mapping().empty()) {
                        aligned_path.mutable_mapping()->erase(aligned_path.mutable_mapping()->begin());
                    }
                    
                    *connecting_subpath->mutable_path() = move(aligned_path);
                    
                    // add the appropriate connections
                    src_subpath->add_next(multipath_aln_out.subpath_size() - 1);
                    connecting_subpath->add_next(edge.first);
                    
                    // translate the path into the space of the main graph unless the path is null
                    if (connecting_subpath->path().mapping_size() != 0) {
                        translate_node_ids(*connecting_subpath->mutable_path(), connect_trans);
                        path_mapping_t* first_subpath_mapping = connecting_subpath->mutable_path()->mutable_mapping(0);
                        if (first_subpath_mapping->position().node_id() == final_mapping.position().node_id()) {
                            first_subpath_mapping->mutable_position()->set_offset(offset(src_pos));
                        }
                    }
                    
#ifdef debug_multipath_alignment
                    cerr << "subpath from " << j << " to " << edge.first << " at index " << multipath_aln_out.subpath_size() - 1 <<  ":" << endl;
                    cerr << debug_string(*connecting_subpath) << endl;
#endif
                }
            }
            
            if (!edges_for_removal.empty()) {
                auto new_end = std::remove_if(src_path_node.edges.begin(), src_path_node.edges.end(),
                                              [&](const pair<size_t, size_t>& edge) {
                                                  return edges_for_removal.count(edge);
                                              });
                src_path_node.edges.resize(new_end - src_path_node.edges.begin());
            }
        }
        
        
        // Now do the tails
        
        // We need to know what subpaths are real sources
        unordered_set<size_t> sources;
        
        // Actually align the tails
        auto tail_alignments = align_tails(alignment, align_graph, aligner, max_alt_alns, dynamic_alt_alns,
                                           max_gap, pessimistic_tail_gap_multiplier, 0, &sources);
                
        // TODO: merge and simplify the tail alignments? rescoring would be kind of a pain...
        
        // Handle the right tails
        for (auto& kv : tail_alignments[true]) {
            // For each sink subpath number with its alignments
            size_t j = kv.first;
            
            // remove alignments with the same path
            auto deduplicated = deduplicate_alt_alns(kv.second, false, true);
        
            PathNode& path_node = path_nodes.at(j);
            pos_t end_pos = final_position(path_node.path);
            
            size_t sink_idx = j;
            
            // zip together identical prefixes and suffixes of the tail alignment
            pair<path_t, int32_t> left_zip_aln, right_zip_aln;
            if (deduplicated.size() > 1) {
                right_zip_aln = zip_alignments(deduplicated, false, alignment, align_graph, path_node.end, aligner);
                left_zip_aln = zip_alignments(deduplicated, true, alignment, align_graph, path_node.end, aligner);
            }
            
            string::const_iterator tail_begin = path_node.end;
            if (!left_zip_aln.first.mapping().empty()) {
                // we were able to zip a prefix of the tails together
                
                multipath_aln_out.mutable_subpath(sink_idx)->add_next(multipath_aln_out.subpath_size());
                sink_idx = multipath_aln_out.subpath_size();
                
                auto subpath = multipath_aln_out.add_subpath();
                *subpath->mutable_path() = move(left_zip_aln.first);
                subpath->set_score(left_zip_aln.second);
                
                auto first_mapping = subpath->mutable_path()->mutable_mapping(0);
                if (first_mapping->position().node_id() == id(end_pos)) {
                    first_mapping->mutable_position()->set_offset(offset(end_pos));
                }
                
                // the zipped alignment is now the attachment point
                end_pos = final_position(subpath->path());
                tail_begin += path_to_length(subpath->path());
                
#ifdef debug_multipath_alignment
                cerr << "left zipped alignment from " << j << " to right tail" << " at index " << multipath_aln_out.subpath_size() - 1 <<  ":" << endl;
                cerr << debug_string(*subpath) << endl;
#endif
            }
            
            // add in the non-redundant middle of the alignments
            bool zip_to_zip_connection = false;
            size_t frayed_tips_begin = multipath_aln_out.subpath_size();
            for (auto& tail_alignment : deduplicated) {
                
                if (tail_alignment.first.mapping().empty()) {
                    // it's possible that an entire tail get's swallowed up in zipping
                    zip_to_zip_connection = true;
                    continue;
                }
                
                multipath_aln_out.mutable_subpath(sink_idx)->add_next(multipath_aln_out.subpath_size());
                
                // TODO: kind of inelegant that i'm relying on there being no zipped alignments
                // and therefore no path that relies on frayed_tips_begin:frayed_tips_end being a
                // a slice of subpaths that align the same sequence
                subpath_t* subpath = nullptr;
                if (deduplicated.size() == 1 && tail_alignment.first.mapping_size() > 1) {
                    // this tail might have soft-clip issues, check if we need to cut it
                    size_t first_idx = multipath_aln_out.subpath_size();
                    add_and_permanently_cut(tail_alignment, tail_begin);
                    subpath = multipath_aln_out.mutable_subpath(first_idx);
                }
                else {
                    subpath = multipath_aln_out.add_subpath();
                    *subpath->mutable_path() = move(tail_alignment.first);
                    subpath->set_score(tail_alignment.second);
                }

                auto first_mapping = subpath->mutable_path()->mutable_mapping(0);
                if (first_mapping->position().node_id() == id(end_pos)) {
                    first_mapping->mutable_position()->set_offset(offset(end_pos));
                }
                else if (right_zip_aln.first.mapping().empty()
                         && subpath->path().mapping_size() == 1 && first_mapping->edit_size() == 1
                         && first_mapping->edit(0).from_length() == 0 && first_mapping->edit(0).to_length() != 0
                         && first_mapping->position().node_id() != id(end_pos)) {
                    // this is a pure soft-clip on the beginning of the next node, we'll move it to the end
                    // of the match node to match invariants expected by other parts of the code base
                    position_t* pos = first_mapping->mutable_position();
                    pos->set_node_id(id(end_pos));
                    pos->set_is_reverse(is_rev(end_pos));
                    pos->set_offset(offset(end_pos));
                }
#ifdef debug_multipath_alignment
                cerr << "subpath from " << sink_idx << " to right tail at index " << multipath_aln_out.subpath_size() - 1 <<  ":" << endl;
                cerr << debug_string(*subpath) << endl;
#endif
            }
            size_t frayed_tips_end = multipath_aln_out.subpath_size();
            
            if (!right_zip_aln.first.mapping().empty()) {
                
                if (right_zip_aln.first.mapping_size() == 1
                    && right_zip_aln.first.mapping(0).edit_size() == 1
                    && right_zip_aln.first.mapping(0).edit(0).from_length() == 0
                    && right_zip_aln.first.mapping(0).edit(0).to_length() != 0) {
                    
                    // this subpath is a soft-clip, to match the expectations of other code we need to
                    // relocate this edit onto the other subpaths (more complicated than in the tail
                    // alignment because there might be multiple previous nodes)
                    
                    auto edit = right_zip_aln.first.mutable_mapping(0)->mutable_edit(0);
                    
                    // add an edit to the tails
                    for (size_t i = frayed_tips_begin; i < frayed_tips_end; ++i) {
                        *multipath_aln_out.mutable_subpath(i)->mutable_path()->mutable_mapping()->back().add_edit() = *edit;
                    }
                    
                    if (zip_to_zip_connection) {
                        // we can't add this edit to the left zipped alignment because it has sucessor
                        // subpaths, need to add it as a separate subpath
                        multipath_aln_out.mutable_subpath(sink_idx)->add_next(multipath_aln_out.subpath_size());
                        
                        auto tail_subpath = multipath_aln_out.add_subpath();
                        auto mapping = tail_subpath->mutable_path()->add_mapping();
                        *mapping->add_edit() = move(*edit);
                        
                        auto pos = mapping->mutable_position();
                        pos->set_node_id(id(end_pos));
                        pos->set_is_reverse(is_rev(end_pos));
                        pos->set_offset(offset(end_pos));
                        tail_subpath->set_score(0);
                        
#ifdef debug_multipath_alignment
                        cerr << "right zipped softclip subpath from " << sink_idx << " to right tail at index " << multipath_aln_out.subpath_size() - 1 <<  ":" << endl;
                        cerr << debug_string(*tail_subpath) << endl;
#endif
                    }
                }
                else {
                    // add edges from all of the alt alns to the zipped tail
                    for (size_t i = frayed_tips_begin; i < frayed_tips_end; ++i) {
                        multipath_aln_out.mutable_subpath(i)->add_next(multipath_aln_out.subpath_size());
                    }
                    if (zip_to_zip_connection) {
                        multipath_aln_out.mutable_subpath(sink_idx)->add_next(multipath_aln_out.subpath_size());
                    }
                    
                    // possibly cut up the final subpath to mitigate soft-clip issues
                    tail_begin += path_to_length(multipath_aln_out.subpath(frayed_tips_begin).path());
                    add_and_permanently_cut(right_zip_aln, tail_begin);
                    auto tail_subpath = &multipath_aln_out.mutable_subpath()->back();
                    
                    // an arbitrary choice of a predecessor (this only matters if
                    // the previous paths end mid-node, in which case all of them must
                    // end on that node)
                    end_pos = final_position(multipath_aln_out.subpath(frayed_tips_begin).path());
                    auto first_mapping = tail_subpath->mutable_path()->mutable_mapping(0);
                    if (first_mapping->position().node_id() == id(end_pos)) {
                        first_mapping->mutable_position()->set_offset(offset(end_pos));
                    }
                    // note: don't need to check for softclips because that is handled separately above
#ifdef debug_multipath_alignment
                    cerr << "right zipped alignment from " << sink_idx << " to right tail at index " << multipath_aln_out.subpath_size() - 1 <<  ":" << endl;
                    cerr << debug_string(*tail_subpath) << endl;
#endif
                }
            }
            
            // TODO: no option to merge if the match is long enough...
            if (deduplicated.size() == 1) {
                // the tail alignment could be merged with its attachment point
                const auto& attachment = multipath_aln_out.subpath(sink_idx).path().mapping().back();
                const auto& pos = attachment.position();
                if (pos.offset() + mapping_from_length(attachment) == align_graph.get_length(align_graph.get_handle(pos.node_id()))
                    && into_cutting_snarl(pos.node_id(), pos.is_reverse(), cutting_snarls, dist_index)) {
                    // the tail is the start of an alignment into a snarl
                    prohibited_merges.insert(sink_idx);
                }
            }
        }
                
        // Now handle the left tails.
        // We need to handle all sources, whether or not they got alignments
        for (auto& j : sources) {
            PathNode& path_node = path_nodes.at(j);
                
            if (path_node.begin != alignment.sequence().begin()) {
                
                // There should be some alignments
                // remove alignments with the same path
                auto deduplicated = deduplicate_alt_alns(tail_alignments[false][j], true, false);
                
                // collapse identical prefixes/suffixes of the alignments
                pair<path_t, int32_t> left_zip_aln, right_zip_aln;
                if (deduplicated.size() > 1) {
                    right_zip_aln = zip_alignments(deduplicated, false, alignment, align_graph, alignment.sequence().begin(), aligner);
                    left_zip_aln = zip_alignments(deduplicated, true, alignment, align_graph, alignment.sequence().begin(), aligner);
                }
                                
                size_t source_idx = j;
                auto first_mapping = &path_node.path.mapping().front();
                
                if (!right_zip_aln.first.mapping().empty()) {
                    // we were able to zip a prefix of the tails together
                    
                    source_idx = multipath_aln_out.subpath_size();
                    
                    auto zip_subpath = multipath_aln_out.add_subpath();
                    *zip_subpath->mutable_path() = move(right_zip_aln.first);
                    zip_subpath->set_score(right_zip_aln.second);
                    
                    first_mapping = &zip_subpath->path().mapping().front();
                    zip_subpath->add_next(j);
                    
#ifdef debug_multipath_alignment
                    cerr << "right zipped alignment from " << j << " to left tail at index " << multipath_aln_out.subpath_size() - 1 <<  ":" << endl;
                    cerr << debug_string(*zip_subpath) << endl;
#endif
                }
                
                bool zip_to_zip_connection = false;
                size_t frayed_tips_begin = multipath_aln_out.subpath_size();
                for (auto& tail_alignment : deduplicated) {
                    if (tail_alignment.first.mapping().empty()) {
                        // it's possible that an entire tail get's swallowed up in zipping
                        zip_to_zip_connection = true;
                        continue;
                    }
                    
                    if (left_zip_aln.first.mapping().empty()) {
                        // the left part of these alignments hasn't been zipped, so this
                        // is a start
                        multipath_aln_out.add_start(multipath_aln_out.subpath_size());
                    }
                    
                    // TODO: kind of inelegant that i'm relying on there being no zipped alignments
                    // and therefore no path that relies on frayed_tips_begin:frayed_tips_end being a
                    // a slice of subpaths that align the same sequence
                    if (deduplicated.size() == 1 && tail_alignment.first.mapping_size() > 1) {
                        // this tail might have soft-clip issues, check if we need to cut it
                        size_t first_idx = multipath_aln_out.subpath_size();
                        add_and_permanently_cut(tail_alignment, alignment.sequence().begin());
                    }
                    else {
                        auto subpath = multipath_aln_out.add_subpath();
                        *subpath->mutable_path() = move(tail_alignment.first);
                        subpath->set_score(tail_alignment.second);
                    }
                    
                    auto tail_subpath = &multipath_aln_out.mutable_subpath()->back();
                    
                    tail_subpath->add_next(source_idx);
                    
#ifdef debug_multipath_alignment
                    cerr << "subpath from " << source_idx << " to left tail at index " << multipath_aln_out.subpath_size() - 1 <<  ":" << endl;
                    cerr << debug_string(*tail_subpath) << endl;
#endif
                    path_mapping_t* final_mapping = tail_subpath->mutable_path()->mutable_mapping(tail_subpath->path().mapping_size() - 1);
                    if (tail_subpath->path().mapping_size() == 1 && final_mapping->edit_size() == 1
                        && final_mapping->edit(0).from_length() == 0 && final_mapping->edit(0).to_length() > 0
                        && final_mapping->position().node_id() != first_mapping->position().node_id()) {
                        // this is a pure soft clip on the end of the previous node, so we move it over to the
                        // beginning of the match node to match invariants in rest of code base
                        *final_mapping->mutable_position() = first_mapping->position();
                    }
                }
                size_t frayed_tips_end = multipath_aln_out.subpath_size();
                
                if (!left_zip_aln.first.mapping().empty()) {
                    if (left_zip_aln.first.mapping_size() == 1
                        && left_zip_aln.first.mapping(0).edit_size() == 1
                        && left_zip_aln.first.mapping(0).edit(0).from_length() == 0
                        && left_zip_aln.first.mapping(0).edit(0).to_length() != 0) {
                        
                        // this is a pure soft clip, we need to attach it to the adjacent nodes
                        // to maintain expected invariants (because they might start on different
                        // nodes)
                        
                        auto edit = &left_zip_aln.first.mutable_mapping()->front().mutable_edit()->front();
                        
                        // add edit to subpaths and add starts
                        for (size_t i = frayed_tips_begin; i < frayed_tips_end; ++i) {
                            auto edits = multipath_aln_out.mutable_subpath(i)->mutable_path()->mutable_mapping()->front().mutable_edit();
                            edits->insert(edits->begin(), *edit);
                            multipath_aln_out.add_start(i);
                        }
                        
                        if (zip_to_zip_connection) {
                            // add subpath for soft clip (can't be added to existing path because it
                            // has sucessor subpaths, need to add it as a separate subpath
                            
                            multipath_aln_out.add_start(multipath_aln_out.subpath_size());
                            
                            auto subpath = multipath_aln_out.add_subpath();
                            subpath->add_next(source_idx);
                            
                            auto mapping = subpath->mutable_path()->add_mapping();
                            *mapping->add_edit() = move(*edit);
                            
                            *mapping->mutable_position() = multipath_aln_out.subpath(source_idx).path().mapping().front().position();
                            subpath->set_score(0);
#ifdef debug_multipath_alignment
                            cerr << "left zipped soft clip subpath from " << source_idx << " to left tail at index " << multipath_aln_out.subpath_size() - 1 <<  ":" << endl;
                            cerr << debug_string(*subpath) << endl;
#endif
                        }
                    }
                    else {
                        
                        multipath_aln_out.add_start(multipath_aln_out.subpath_size());
                        
                        add_and_permanently_cut(left_zip_aln, alignment.sequence().begin());
                        
                        auto subpath = &multipath_aln_out.mutable_subpath()->back();
                        
                        for (size_t i = frayed_tips_begin; i < frayed_tips_end; ++i) {
                            subpath->add_next(i);
                        }
                        if (zip_to_zip_connection) {
                            subpath->add_next(source_idx);
                        }
#ifdef debug_multipath_alignment
                        cerr << "left zipped subpath from " << source_idx << " to left tail at index " << multipath_aln_out.subpath_size() - 1 <<  ":" << endl;
                        cerr << debug_string(*subpath) << endl;
#endif
                    }
                }
                
                // TODO: no option to merge if the match is long enough...
                if (deduplicated.size() == 1) {
                    // the tail alignment could be merged with its attachment point
                    const auto& attachment = multipath_aln_out.subpath(source_idx).path().mapping().front();
                    const auto& pos = attachment.position();
                    if (pos.offset() == 0 && into_cutting_snarl(pos.node_id(), !pos.is_reverse(), cutting_snarls, dist_index)) {
                        // the tail is the start of an alignment into a snarl
                        prohibited_merges.insert(frayed_tips_begin);
                    }
                }
            }
            else {
#ifdef debug_multipath_alignment
                cerr << "No tails left off of subpath " << j << "; it can be a start" << endl;
#endif
                multipath_aln_out.add_start(j);
            }
        }
        
        if (simplify_topologies) {
#ifdef debug_multipath_alignment
            cerr << "merging non-branching subpaths" << endl;
            size_t num_pre_merge = multipath_aln_out.subpath_size();
#endif
            
            merge_non_branching_subpaths(multipath_aln_out, &prohibited_merges);
            
#ifdef debug_multipath_alignment
            cerr << "reduce from " << num_pre_merge << " to " << multipath_aln_out.subpath_size() << " subpaths during merge" << endl;
#endif
        }
    }
    
    // just to control the memory usage to a small value
    const size_t MultipathAlignmentGraph::tail_gap_memo_max_size = 1000;

    // make the memo live in this .o file
    thread_local unordered_map<double, vector<int64_t>> MultipathAlignmentGraph::pessimistic_tail_gap_memo;

    int64_t MultipathAlignmentGraph::pessimistic_tail_gap(int64_t tail_length, double multiplier) {
        int64_t gap_length;
        if (tail_length >= tail_gap_memo_max_size) {
            gap_length = multiplier * sqrt(tail_length);
        }
        else {
            vector<int64_t>& memo = pessimistic_tail_gap_memo[multiplier];
            while (memo.size() <= tail_length) {
                memo.emplace_back(multiplier * sqrt(memo.size()));
            }
            gap_length = memo[tail_length];
        }
        return gap_length;
    }
    
    unordered_map<bool, unordered_map<size_t, vector<Alignment>>>
    MultipathAlignmentGraph::align_tails(const Alignment& alignment, const HandleGraph& align_graph, const GSSWAligner* aligner,
                                         size_t max_alt_alns, bool dynamic_alt_alns, size_t max_gap, double pessimistic_tail_gap_multiplier,
                                         size_t min_paths, unordered_set<size_t>* sources) {
        
#ifdef debug_multipath_alignment
        cerr << "doing tail alignments to:" << endl;
        to_dot(cerr);
#endif
        
        // TODO: magic number
        int64_t anchor_low_cmplx_len = 16;
        
        // multiplier to account for low complexity sequences
        auto low_complexity_multiplier = [](string::const_iterator begin, string::const_iterator end) {
            // TODO: magic numbers
            static const double seq_cmplx_alpha = 0.05;
            static const double max_multiplier = 32.0;
            SeqComplexity<2> complexity(begin, end);
            if (complexity.p_value(1) < seq_cmplx_alpha || complexity.p_value(2) < seq_cmplx_alpha) {
                double repetitiveness = max(complexity.repetitiveness(1), complexity.repetitiveness(2));
                double low_complexity_len = repetitiveness * (end - begin);
                // TODO: empirically-derived function, not very principled
                return max(1.0, min(max_multiplier, 0.05 * low_complexity_len * low_complexity_len));
            }
            else {
                return 1.0;
            }
        };
                
        // Make a structure to populate
        unordered_map<bool, unordered_map<size_t, vector<Alignment>>> to_return;
        auto& left_alignments = to_return[false];
        auto& right_alignments = to_return[true];
        
        vector<bool> is_source_node(path_nodes.size(), true);
        for (size_t j = 0; j < path_nodes.size(); j++) {
            PathNode& path_node = path_nodes.at(j);
#ifdef debug_multipath_alignment
            cerr << "Visit PathNode " << j << " with " << path_node.edges.size() << " outbound edges" << endl;
#endif
            if (!path_node.edges.empty()) {
                // We go places from here.
                for (const pair<size_t, size_t>& edge : path_node.edges) {
                    // Make everywhere we go as not a source
                    is_source_node[edge.first] = false;
#ifdef debug_multipath_alignment
                    cerr << "Edge " << j << " -> " << edge.first << " makes " << edge.first << " not a source" << endl;
#endif
                }
            }
            else if (path_node.end != alignment.sequence().end()) {

#ifdef debug_multipath_alignment
                cerr << "doing right end alignment from sink node " << j << " with path " << debug_string(path_node.path) << " and sequence ";
                for (auto iter = path_node.begin; iter != path_node.end; iter++) {
                    cerr << *iter;
                }
                cerr << endl;
#endif
                
                // figure out how long we need to try to align out to
                int64_t tail_length = alignment.sequence().end() - path_node.end;
                int64_t gap =  min(aligner->longest_detectable_gap(alignment, path_node.end), max_gap);
                if (pessimistic_tail_gap_multiplier) {
                    gap = min(gap, pessimistic_tail_gap(tail_length, pessimistic_tail_gap_multiplier));
                }
                int64_t target_length = tail_length + gap;
                
                
                
                pos_t end_pos = final_position(path_node.path);
                
                bdsg::HashGraph tail_graph;
                unordered_map<id_t, id_t> tail_trans = algorithms::extract_extending_graph(&align_graph,
                                                                                           &tail_graph,
                                                                                           target_length,
                                                                                           end_pos,
                                                                                           false,         // search forward
                                                                                           false);        // no need to preserve cycles (in a DAG)
                
                size_t num_alt_alns;
                if (dynamic_alt_alns) {
                    size_t num_paths = handlealgs::count_walks(&tail_graph);
                    if (num_paths < min_paths) {
                        continue;
                    }
                    num_alt_alns = min(max_alt_alns, num_paths);
                }
                else {
                    num_alt_alns = max_alt_alns;
                }
                
                if (num_alt_alns > 0) {
                    
                    // get the sequence remaining in the right tail
                    Alignment right_tail_sequence;
                    right_tail_sequence.set_sequence(alignment.sequence().substr(path_node.end - alignment.sequence().begin(),
                                                                                 alignment.sequence().end() - path_node.end));
                    if (!alignment.quality().empty()) {
                        right_tail_sequence.set_quality(alignment.quality().substr(path_node.end - alignment.sequence().begin(),
                                                                                   alignment.sequence().end() - path_node.end));
                    }
                    
#ifdef debug_multipath_alignment
                    cerr << "making " << num_alt_alns << " alignments of sequence: " << right_tail_sequence.sequence() << endl << "to right tail graph" << endl;
                    tail_graph.for_each_handle([&](const handle_t& handle) {
                        cerr << tail_graph.get_id(handle) << " " << tail_graph.get_sequence(handle) << endl;
                        tail_graph.follow_edges(handle, true, [&](const handle_t& prev) {
                            cerr << "\t" << tail_graph.get_id(prev) << " <-" << endl;
                        });
                        tail_graph.follow_edges(handle, false, [&](const handle_t& next) {
                            cerr << "\t-> " << tail_graph.get_id(next) << endl;
                        });
                    });
#endif
                    
                    // align against the graph
                    auto& alt_alignments = right_alignments[j];
                    if (num_alt_alns == 1) {
#ifdef debug_multipath_alignment
                        cerr << "align right with dozeu with gap " << gap << endl;
#endif
                        // we can speed things up by using the dozeu pinned alignment
                        alt_alignments.emplace_back(move(right_tail_sequence));
                        aligner->align_pinned(alt_alignments.back(), tail_graph, true, true, gap);
                    }
                    else {
                        
#ifdef debug_multipath_alignment
                        cerr << "align right with gssw" << endl;
#endif
                        
                        // TODO: it would be cleaner to handle these in one multiplier, but i already tuned
                        // this for the tails and i don't feel like doing that again...
                        double tail_multiplier = low_complexity_multiplier(right_tail_sequence.sequence().begin(),
                                                                           right_tail_sequence.sequence().end());
                        double anchor_multiplier = low_complexity_multiplier(max(path_node.end - anchor_low_cmplx_len,
                                                                                 alignment.sequence().begin()),
                                                                             path_node.end);
                        double multiplier = max(tail_multiplier, anchor_multiplier);
                        if (multiplier != 1.0) {
                            num_alt_alns = round(multiplier * num_alt_alns);
#ifdef debug_multipath_alignment
                            cerr << "increase num alns for low complexity sequence to " << num_alt_alns << endl;
#endif
                        }
                        aligner->align_pinned_multi(right_tail_sequence, alt_alignments, tail_graph, true, num_alt_alns);
                    }
                    
                    // Translate back into non-extracted graph.
                    // Make sure to account for having removed the left end of the cut node relative to end_pos
                    for (auto& aln : alt_alignments) {
                        // We always remove end_pos's offset, since we
                        // search forward from it to extract the subgraph,
                        // but that may be from the left or right end of
                        // its node depending on orientation.
                        translate_node_ids(*aln.mutable_path(), tail_trans, id(end_pos), offset(end_pos), is_rev(end_pos));
                    }
#ifdef debug_multipath_alignment
                    cerr << "made " << alt_alignments.size() << " tail alignments" << endl;
                    for (size_t i = 0; i < alt_alignments.size(); ++i) {
                        cerr << i << ": " << pb2json(alt_alignments[i]) << endl;
                    }
#endif
                }
            }
        }
        
        // Now we know all the sources so we can do them and the left tails.
        for (size_t j = 0; j < path_nodes.size(); j++) {
            if (is_source_node[j]) {
                // This is a source, whether or not it has a tail to the left.
                
#ifdef debug_multipath_alignment
                cerr << "PathNode " << j << " had no edges into it and is a source" << endl;
#endif
                
                if (sources != nullptr) {
                    // Remember it.
                    sources->insert(j);
                }
            
                PathNode& path_node = path_nodes[j];
                if (path_node.begin != alignment.sequence().begin()) {
                    
                    // figure out how far we need to try to align out to
                    int64_t tail_length = path_node.begin - alignment.sequence().begin();
                    int64_t gap =  min(aligner->longest_detectable_gap(alignment, path_node.begin), max_gap);
                    if (pessimistic_tail_gap_multiplier) {
                        gap = min(gap, pessimistic_tail_gap(tail_length, pessimistic_tail_gap_multiplier));
                    }
                    int64_t target_length = tail_length + gap;
                    
                    pos_t begin_pos = initial_position(path_node.path);
                    
                    
                    bdsg::HashGraph tail_graph;
                    unordered_map<id_t, id_t> tail_trans = algorithms::extract_extending_graph(&align_graph,
                                                                                               &tail_graph,
                                                                                               target_length,
                                                                                               begin_pos,
                                                                                               true,          // search backward
                                                                                               false);        // no need to preserve cycles (in a DAG)
                    
                    size_t num_alt_alns;
                    if (dynamic_alt_alns) {
                        size_t num_paths = handlealgs::count_walks(&tail_graph);
                        if (num_paths < min_paths) {
                            continue;
                        }
                        num_alt_alns = min(max_alt_alns, num_paths);
                    }
                    else {
                        num_alt_alns = max_alt_alns;
                    }
                            
                    if (num_alt_alns > 0) {
                        
                        Alignment left_tail_sequence;
                        left_tail_sequence.set_sequence(alignment.sequence().substr(0, path_node.begin - alignment.sequence().begin()));
                        if (!alignment.quality().empty()) {
                            left_tail_sequence.set_quality(alignment.quality().substr(0, path_node.begin - alignment.sequence().begin()));
                        }
                        
#ifdef debug_multipath_alignment
                        cerr << "making " << num_alt_alns << " alignments of sequence: " << left_tail_sequence.sequence() << endl << "to left tail graph" << endl;
                        tail_graph.for_each_handle([&](const handle_t& handle) {
                            cerr << tail_graph.get_id(handle) << " " << tail_graph.get_sequence(handle) << endl;
                            tail_graph.follow_edges(handle, false, [&](const handle_t& next) {
                                cerr << "\t-> " << tail_graph.get_id(next) << endl;
                            });
                            tail_graph.follow_edges(handle, true, [&](const handle_t& prev) {
                                cerr << "\t" << tail_graph.get_id(prev) << " <-" << endl;
                            });
                        });
#endif
                        
                        // align against the graph
                        auto& alt_alignments = left_alignments[j];
                        if (num_alt_alns == 1) {
#ifdef debug_multipath_alignment
                            cerr << "align left with dozeu using gap " << gap << endl;
#endif
                            // we can speed things up by using the dozeu pinned alignment
                            alt_alignments.emplace_back(move(left_tail_sequence));
                            aligner->align_pinned(alt_alignments.back(), tail_graph, false, true, gap);
                        }
                        else {
#ifdef debug_multipath_alignment
                            cerr << "align left with gssw" << endl;
#endif
                            double anchor_multiplier = low_complexity_multiplier(path_node.begin,
                                                                                 min(path_node.begin + anchor_low_cmplx_len,
                                                                                     alignment.sequence().end()));
                            double tail_multiplier = low_complexity_multiplier(left_tail_sequence.sequence().begin(),
                                                                               left_tail_sequence.sequence().end());
                            double multiplier = max(tail_multiplier, anchor_multiplier);
                            if (multiplier != 1.0) {
                                num_alt_alns = round(multiplier * num_alt_alns);
#ifdef debug_multipath_alignment
                                cerr << "increase num alns for low complexity sequence to " << num_alt_alns << endl;
#endif
                            }
                            
                            aligner->align_pinned_multi(left_tail_sequence, alt_alignments, tail_graph, false, num_alt_alns);
                        }
                        
                        // Translate back into non-extracted graph.
                        // Make sure to account for having removed the right end of the cut node relative to begin_pos
                        for (auto& aln : alt_alignments) {
                            // begin_pos's offset is how much we keep, so we removed the node's length minus that.
                            // And by default it comes off the right side of the node.
                            translate_node_ids(*aln.mutable_path(), tail_trans, 
                                id(begin_pos),
                                align_graph.get_length(align_graph.get_handle(id(begin_pos))) - offset(begin_pos),
                                !is_rev(begin_pos));
                        }
#ifdef debug_multipath_alignment
                        cerr << "made " << alt_alignments.size() << " tail alignments" << endl;
                        for (size_t i = 0; i < alt_alignments.size(); ++i) {
                            cerr << i << ": " << pb2json(alt_alignments[i]) << endl;
                        }
#endif
                    }
                }
            }
        }
        
        // GSSW does some weird things with N's that we want to normalize away
         if (find(alignment.sequence().begin(), alignment.sequence().end(), 'N') != alignment.sequence().end()) {
             for (bool side : {true, false}) {
                 for (pair<const size_t, vector<Alignment>>& tail_alignments : to_return[side]) {
                     for (Alignment& aln : tail_alignments.second) {
                         normalize_alignment(aln);
                     }
                 }
             }
         }
        
#ifdef debug_multipath_alignment
        cerr << "made alignments for " << to_return[false].size() << " source PathNodes and " << to_return[true].size() << " sink PathNodes" << endl;
        if (sources != nullptr) {
            cerr << "Identified " << sources->size() << " source PathNodes" << endl;
        }
#endif 
        
        // Return all the alignments, organized by tail and subpath
        return to_return;
    }
    
    bool MultipathAlignmentGraph::empty() const {
        return path_nodes.empty();
    }

    size_t MultipathAlignmentGraph::size() const {
        return path_nodes.size();
    }
    
    void MultipathAlignmentGraph::to_dot(ostream& out, const Alignment* alignment) const {
        // We track the VG graph nodes we talked about already.
        set<id_t> mentioned_nodes;
        set<pair<id_t, id_t>> mentioned_edges;
    
        out << "digraph graphname {" << endl;
        out << "rankdir=\"LR\";" << endl;
        for (size_t i = 0; i < path_nodes.size(); i++) {
            // For each node, say the node itself as a mapping node, annotated with match length
            out << "m" << i << " [label=\"" << i << "\" shape=circle tooltip=\""
                << (path_nodes.at(i).end - path_nodes.at(i).begin) << " bp";
            if (alignment != nullptr) {
                // Include info on where that falls in the query sequence
                out << " (" << (path_nodes.at(i).begin - alignment->sequence().begin()) << " - "
                    << (path_nodes.at(i).end - alignment->sequence().begin()) << ")";
            }
            out << "\"];" << endl;
            for (pair<size_t, size_t> edge : path_nodes.at(i).edges) {
                // For each edge from it, say where it goes and how far it skips
                out << "m" << i << " -> m" << edge.first << " [label=" << edge.second << "];" << endl;
            }
            auto& path = path_nodes.at(i).path;
            for (size_t j = 0; j < path.mapping_size(); j++) {
                // For each mapping in the path, show the vg node in the graph too
                auto node_id = path.mapping(j).position().node_id();
                
                if (!mentioned_nodes.count(node_id)) {
                    // This graph node eneds to be made
                    mentioned_nodes.insert(node_id);
                    out << "g" << node_id << " [label=\"" << node_id << "\" shape=box];" << endl;
                }
                
                // Attach the mapping to each involved graph node.
                out << "m" << i << " -> g" << node_id << " [dir=none color=blue];" << endl;
                
                if (j != 0) {
                    // We have a previous node in this segment of path. What is it?
                    auto prev_id = path.mapping(j-1).position().node_id();
                    pair<id_t, id_t> edge_pair{prev_id, node_id};
                    
                    if (!mentioned_edges.count(edge_pair)) {
                        // This graph edge needs to be made
                        mentioned_edges.insert(edge_pair);
                        
                        out << "g" << prev_id << " -> g" << node_id << ";" << endl;
                    }
                }
            }
        }
        out << "}" << endl;
    }

    bool MultipathAlignmentGraph::into_cutting_snarl(id_t node_id, bool is_rev,
                                                     SnarlManager* snarl_manager, MinimumDistanceIndex* dist_index) const {
        
        if (dist_index) {
            auto result = dist_index->into_which_snarl(node_id, is_rev);
            // points into a snarl and snarl is nontrivial
            return get<0>(result) != 0 && !get<2>(result);
        }
        else if (snarl_manager) {
            // no mechanism to check for triviality without traversing graph
            return snarl_manager->into_which_snarl(node_id, is_rev);
        }
        else {
            // no snarls provided
            return false;
        }
    }
    
    vector<vector<id_t>> MultipathAlignmentGraph::get_connected_components() const {
        // Brea all the path_nodes into components using this union-find
        structures::UnionFind unionfind(path_nodes.size());
        
        for (size_t i = 0; i < path_nodes.size(); i++) {
            // For each node
            for (auto& edge : path_nodes.at(i).edges) {
                // For each outgoing edge...
                
                // Union both ends into the same component.
                unionfind.union_groups(unionfind.find_group(i), unionfind.find_group(edge.first));
            }
        }
        
        // Assemble the vectors of ids
        vector<vector<id_t>> to_return;
        
        for (auto& component : unionfind.all_groups()) {
            // For each connected component
            
            // Make a set of vg graph IDs in it
            set<id_t> vg_ids;
            
            for (auto& path_node : component) {
                // For each path in the vg graph used by the component
                auto& path = path_nodes.at(path_node).path;
                
                for (auto& mapping : path.mapping()) {
                    // For each mapping in the path, remember its vg graph node
                    vg_ids.insert(mapping.position().node_id());
                }
            }
            
            // Copy the unique vg nodes used byt this component of PathNodes into the list of components to return
            to_return.emplace_back(vg_ids.begin(), vg_ids.end());
        }
        
        return to_return;
        
    }
}






