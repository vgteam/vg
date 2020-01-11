/**
 * \file
 * surjector.cpp: implements a class that surjects alignments onto paths
 */

#include "surjector.hpp"

//#define debug_spliced_surject
//#define debug_anchored_surject
//#define debug_validate_anchored_multipath_alignment

namespace vg {

using namespace std;
    
    Surjector::Surjector(const PathPositionHandleGraph* graph) : graph(graph) {
        if (!graph) {
            cerr << "error:[Surjector] Failed to provide an graph to the Surjector" << endl;
        }
    }
    
    Alignment Surjector::surject(const Alignment& source, const set<string>& path_names,
                                 bool allow_negative_scores, bool preserve_deletions) const {
    
        // Allocate the annotation info
        string path_name_out;
        int64_t path_pos_out;
        bool path_rev_out;
        
        // Do the surjection
        Alignment surjected = surject(source, path_names, path_name_out, path_pos_out, path_rev_out, allow_negative_scores, preserve_deletions);
        
        // Pack all the info into the refpos field
        surjected.clear_refpos();
        auto* pos = surjected.add_refpos();
        pos->set_name(path_name_out);
        pos->set_offset(path_pos_out);
        pos->set_is_reverse(path_rev_out);
    
        return surjected;
    }
    
    Alignment Surjector::surject(const Alignment& source, const set<string>& path_names,
                                 string& path_name_out, int64_t& path_pos_out, bool& path_rev_out,
                                 bool allow_negative_scores, bool preserve_deletions) const {

#ifdef debug_anchored_surject
        cerr << "surjecting alignment: " << pb2json(source) << " onto paths ";
        for (const string& path_name : path_names) {
            cerr << path_name << " ";
        }
        cerr << endl;
#endif
        
        if (source.path().mapping_size() != 0) {
            // The read is mapped. Check the input alignment for basic
            // consistency. If the sequence and the graph path don't agree
            // about the read length, something is very wrong with the input.
            size_t source_to_length = path_to_length(source.path());
            if (source.sequence().size() != source_to_length) {
                cerr << "error[Surjector::surject]: read " << source.name() << " has "
                << source.sequence().size() << " sequence bases but an input alignment that aligns "
                << source_to_length << " bases instead. This is invalid and uninterpretable; check your mapper." << endl;
                exit(1);
            }
        }
        
        // translate the path names and path handles
        unordered_set<path_handle_t> surjection_path_handles;
        for (const string& path_name : path_names) {
            surjection_path_handles.insert(graph->get_path_handle(path_name));
        }
        
        // make an overlay that will memoize the results of some expensive XG operations
        MemoizingGraph memoizing_graph(graph);
        
        // get the chunks of the aligned path that overlap the ref path
        auto path_overlapping_anchors = extract_overlapping_paths(&memoizing_graph, source, surjection_path_handles);
        
#ifdef debug_anchored_surject
        cerr << "got path overlapping segments" << endl;
        for (const auto& surjection_record : path_overlapping_anchors) {
            cerr << "path " << graph->get_path_name(surjection_record.first) << endl;
            for (auto& anchor : surjection_record.second.first) {
                cerr << "\t read[" << (anchor.first.first - source.sequence().begin()) << ":" << (anchor.first.second - source.sequence().begin()) << "] : ";
                for (auto iter = anchor.first.first; iter != anchor.first.second; iter++) {
                    cerr << *iter;
                }
                cerr << endl;
                cerr << "\t" << pb2json(anchor.second) << endl;
            }
        }
#endif
        // the surjected alignment for each path we overlapped
        unordered_map<path_handle_t, Alignment> path_surjections;
        for (pair<const path_handle_t, pair<vector<path_chunk_t>, vector<pair<step_handle_t, step_handle_t>>>>& surjection_record : path_overlapping_anchors) {
            if (!preserve_deletions) {
                path_surjections[surjection_record.first] = realigning_surject(&memoizing_graph, source, surjection_record.first,
                                                                               surjection_record.second.first, allow_negative_scores);
            }
            else {
                path_surjections[surjection_record.first] = spliced_surject(&memoizing_graph, source, surjection_record.first,
                                                                            surjection_record.second.first, surjection_record.second.second,
                                                                            allow_negative_scores);
            }
        }
        
        // in case we didn't overlap any paths, add a sentinel so the following code still executes correctly
        if (path_surjections.empty()) {
            path_surjections[handlegraph::as_path_handle(-1)] = make_null_alignment(source);
        }
        
        // choose which path surjection was best
        path_handle_t best_path_handle;
        int32_t score = numeric_limits<int32_t>::min();
        for (const auto& surjection : path_surjections) {
            if (surjection.second.score() >= score) {
                score = surjection.second.score();
                best_path_handle = surjection.first;
            }
        }
        Alignment& best_surjection = path_surjections[best_path_handle];
                
        // find the position along the path (do it greedily if we preserved deletions so that
        // it doesn't blow up when we don't match the path exactly)
        if (preserve_deletions) {
            set_path_position_inexact(&memoizing_graph, best_surjection, best_path_handle,
                                      path_name_out, path_pos_out, path_rev_out);
        }
        else {
            set_path_position(&memoizing_graph, best_surjection, best_path_handle,
                              path_name_out, path_pos_out, path_rev_out);
        }
        
        
#ifdef debug_anchored_surject
        cerr << "chose path " << path_name_out << " at position " << path_pos_out << (path_rev_out ? "-" : "+") << endl;
#endif
        return move(best_surjection);
        
    }

    vector<vector<size_t>> Surjector::reverse_adjacencies(const vector<vector<size_t>>& adj) const {
        // make a reverse adjacency list
        vector<vector<size_t>> rev_adj(adj.size());
        for (size_t i = 0; i < adj.size(); ++i) {
            for (size_t j : adj[i]) {
                rev_adj[j].push_back(i);
            }
        }
        return rev_adj;
    }

    vector<size_t> Surjector::connected_components(const vector<vector<size_t>>& adj, const vector<vector<size_t>>& rev_adj,
                                                   size_t* num_comps_out) const {
        
        // DFS to find connected components
        vector<bool> enqueued(adj.size(), false);
        vector<size_t> comps(adj.size());
        size_t curr_comp = 0;
        for (size_t i = 0; i < adj.size(); ++i) {
            if (!enqueued[i]) {
                vector<size_t> stack(1, i);
                enqueued[i] = true;
                while (!stack.empty()) {
                    size_t here = stack.back();
                    stack.pop_back();
                    comps[here] = curr_comp;
                    for (const vector<vector<size_t>>* adj_list : {&adj, &rev_adj}) {
                        for (size_t j : (*adj_list)[here]) {
                            if (!enqueued[j]) {
                                stack.push_back(j);
                                enqueued[j] = true;
                            }
                        }
                    }
                    
                }
                
                curr_comp += 1;
            }
        }
        
        if (num_comps_out) {
            *num_comps_out = curr_comp;
        }
        
        return comps;
    }

    vector<vector<size_t>> Surjector::transitive_reduction(const vector<vector<size_t>>& adj) const {
        
        // by construction the graph here has edges in topological order
        
        vector<vector<size_t>> reduction(adj.size());
        
        for (size_t i = 0; i < adj.size(); ++i) {
            
            const vector<size_t>& edges = adj[i];
            
            if (edges.size() == 1) {
                // optimization: a single edge out can never be transitive
                reduction[i].push_back(edges[0]);
                continue;
            }
            
            vector<bool> traversed(adj.size() - i, false);
            
            for (size_t j = 0; j < edges.size(); j++) {
                
                size_t edge = edges[j];
                if (traversed[edge - i]) {
                    // we can reach the target of this edge by another path, so it is transitive
                    continue;
                }
                
                // this edge reaches a target we haven't traversed to yet, so it can't be transitive
                reduction[i].push_back(edge);
                
                // DFS to mark all reachable nodes from this edge
                vector<size_t> stack(1, edge);
                traversed[edge - i] = true;
                while (!stack.empty()) {
                    size_t idx = stack.back();
                    stack.pop_back();
                    for (size_t k : adj[idx]) {
                        if (!traversed[k - i]) {
                            stack.push_back(k);
                            traversed[k - i] = true;
                        }
                    }
                }
            }
            
            
        }
        
        return reduction;
    }

    vector<size_t> Surjector::find_constriction_edges(const vector<vector<size_t>>& adj) const {
        
        auto rev_adj = reverse_adjacencies(adj);
        
        size_t num_comps = 0;
        auto comps = connected_components(adj, rev_adj, &num_comps);
        
        vector<size_t> fwd(adj.size(), 0), bwd(adj.size(), 0);
        vector<size_t> total_comp_paths(num_comps, 0);
        
        for (int64_t i = adj.size() - 1; i >= 0; --i) {
            if (adj[i].empty()) {
                // sink node
                bwd[i] = 1;
            }
            else {
                for (size_t j : adj[i]) {
                    bwd[i] += bwd[j];
                }
            }
            if (rev_adj[i].empty()) {
                // source node, add paths to total
                total_comp_paths[comps[i]] += bwd[i];
            }
        }
        
        vector<size_t> constrictions;
        
        for (int64_t i = 0; i < adj.size(); ++i) {
            if (rev_adj[i].empty()) {
                // source node
                fwd[i] = 1;
            }
            else {
                for (size_t j : rev_adj[i]) {
                    fwd[i] += fwd[j];
                }
            }
            
            if (adj[i].size() == 1) {
                if (fwd[i] * bwd[adj[i][0]] == total_comp_paths[comps[i]]) {
                    // all paths in this component traverse this edge, must be a constriction
                    
                    constrictions.push_back(i);
                }
            }
        }
        
        return constrictions;
    }

    Alignment Surjector::spliced_surject(const PathPositionHandleGraph* path_position_graph, const Alignment& source,
                                         const path_handle_t& path_handle, const vector<path_chunk_t>& path_chunks,
                                         const vector<pair<step_handle_t, step_handle_t>>& ref_chunks,
                                         bool allow_negative_scores) const {
        
        cerr << source.name() << endl;
        
        assert(path_chunks.size() == ref_chunks.size());
        
        auto get_strand = [&](size_t i) {
            bool strand_is_rev;
            if (ref_chunks[i].first == ref_chunks[i].second) {
                // this whole path is only on one node, so we can determine strand by checking the orientation
                // of the ref path step and the aligned path position
                strand_is_rev = (path_position_graph->get_is_reverse(path_position_graph->get_handle_of_step(ref_chunks[i].first))
                                 != path_chunks[i].second.mapping(0).position().is_reverse());
            }
            else {
                // the ref chunk is in the order of the path, so we can use step offsets to tell
                // whether the path is on the forward or reverse strand
                strand_is_rev = (path_position_graph->get_position_of_step(ref_chunks[i].first)
                                 > path_position_graph->get_position_of_step(ref_chunks[i].second));
            }
            return strand_is_rev;
        };
        
#ifdef debug_spliced_surject
        cerr << "making colinearity graph for " << path_chunks.size() << " path chunks" << endl;
#endif
        
        // by construction, the path chunks are ordered by initial index on aln path, but not necessarily second index
        vector<vector<size_t>> colinear_adj(path_chunks.size());
        
        for (size_t i = 0; i < path_chunks.size(); ++i) {
            
            bool strand = get_strand(i);
            
            for (size_t j = i + 1; j < path_chunks.size(); ++j) {
                                
                if (get_strand(j) != strand) {
                    // these two path chunks are not on the same strand of the path
                    continue;
                }
                
                if (path_chunks[j].first.first >= path_chunks[i].first.second &&
                    (graph->get_position_of_step(ref_chunks[j].first) > graph->get_position_of_step(ref_chunks[i].second)) != strand) {
                    // the second one is further along both the read and the path, so it is colinear
                    colinear_adj[i].push_back(j);
                }
            }
        }
        
#ifdef debug_spliced_surject
        cerr << "computing transitive reduction" << endl;
#endif
        
        // remove transitive edges
        vector<vector<size_t>> colinear_adj_red = transitive_reduction(colinear_adj);
        
#ifdef debug_spliced_surject
        cerr << "finding constriction edges" << endl;
#endif
        
        // find edges that constrict the colinearity graph
        vector<size_t> constrictions = find_constriction_edges(colinear_adj_red);
        
        // remove all the constrictions except those that are pure deletions
        size_t removed_so_far = 0;
        for (size_t i : constrictions) {
            if (path_chunks[i].first.second != path_chunks[colinear_adj_red[i][0]].first.first) {
                // not a deletion, does not directly abut on read
                ++removed_so_far;
            }
            else if (removed_so_far) {
                constrictions[i - removed_so_far] = constrictions[i];
            }
        }
        constrictions.resize(constrictions.size() - removed_so_far);
        
#ifdef debug_spliced_surject
        cerr << "found " << constrictions.size() << " deletion constrictions" << endl;
#endif
        
        // remove the constrictions and record their targets
        vector<size_t> constriction_targets(constrictions.size());
        for (size_t i = 0; i < constrictions.size(); ++i) {
            constriction_targets[i] = colinear_adj_red[constrictions[i]][0];
            colinear_adj_red[constrictions[i]].clear();
        }
        
#ifdef debug_spliced_surject
        cerr << "computing constriction components" << endl;
#endif
                
        size_t num_comps = 0;
        vector<size_t> constriction_comps = connected_components(colinear_adj_red,
                                                                 reverse_adjacencies(colinear_adj_red),
                                                                 &num_comps);
        
        vector<vector<size_t>> comp_groups(num_comps);
        for (size_t i = 0; i < constriction_comps.size(); ++i) {
            comp_groups[constriction_comps[i]].push_back(i);
        }
        
#ifdef debug_spliced_surject
        cerr << "surjecting " << comp_groups.size() << " constriction sections" << endl;
#endif
        
        vector<Alignment> sections;
        vector<pair<string::const_iterator, string::const_iterator>> read_ranges;
        vector<pair<step_handle_t, step_handle_t>> ref_ranges;
        
        for (size_t i = 0; i < comp_groups.size(); ++i) {
            pair<string::const_iterator, string::const_iterator> read_range;
            pair<step_handle_t, step_handle_t> ref_range;
            vector<path_chunk_t> section_path_chunks;
            
            
            
            bool strand = get_strand(comp_groups[i].front());
            
            vector<size_t>& group = comp_groups[i];
            for (size_t j = 0; j < group.size(); ++j) {
                section_path_chunks.push_back(path_chunks[group[j]]);
                if (j == 0 || read_range.first > path_chunks[group[j]].first.first) {
                    read_range.first = path_chunks[group[j]].first.first;
                }
                if (j == 0 || read_range.first < path_chunks[group[j]].first.second) {
                    read_range.second = path_chunks[group[j]].first.second;
                }
                if (j == 0 ||
                    (graph->get_position_of_step(ref_chunks[group[j]].first) < graph->get_position_of_step(ref_range.first)) != strand) {
                    ref_range.first = ref_chunks[group[j]].first;
                }
                if (j == 0 ||
                    (graph->get_position_of_step(ref_chunks[group[j]].second) > graph->get_position_of_step(ref_range.second)) != strand) {
                    ref_range.second = ref_chunks[group[j]].second;
                }
            }
            
            // make a dummy alignment with the relevant portion of the sequence
            Alignment section_source;
            *section_source.mutable_sequence() = string(read_range.first, read_range.second);
            if (!source.quality().empty()) {
                *section_source.mutable_quality() = string(source.quality().begin() + (read_range.first - source.sequence().begin()),
                                                           source.quality().begin() + (read_range.second - source.sequence().begin()));
            }
            
            // update the path chunk ranges to point into the dummy section read
            for (size_t i = 0; i < section_path_chunks.size(); ++i) {
                auto& range = section_path_chunks[i].first;
                range.first = section_source.sequence().begin() + (range.first - read_range.first);
                range.second = section_source.sequence().begin() + (range.second - read_range.first);
            }
            
#ifdef debug_spliced_surject
            cerr << "surjecting section " << i << ": " << pb2json(section_source) << endl;
#endif
            
            // perform a full length surjection within the section section
            sections.push_back(realigning_surject(graph, section_source, path_handle, section_path_chunks, false));
            read_ranges.push_back(read_range);
            ref_ranges.push_back(ref_range);
            
            // remove any extraneous full length bonuses
            // TODO: technically, this can give a non-optimal alignment because it's post hoc to the dynamic programming
            if (sections.back().path().mapping_size() > 0) {
                if (read_range.first != source.sequence().begin()) {
                    if (sections.back().path().mapping(0).edit(0).from_length() > 0) {
                        sections.back().set_score(sections.back().score() - get_aligner()->full_length_bonus);
                    }
                }
                if (read_range.second != source.sequence().end()) {
                    const Mapping& m = sections.back().path().mapping(0);
                    if (m.edit(m.edit_size() - 1).from_length() > 0) {
                        sections.back().set_score(sections.back().score() - get_aligner()->full_length_bonus);
                    }
                }
            }
            
#ifdef debug_spliced_surject
            cerr << "surjected section " << i << " after score adjustment: " << pb2json(sections.back()) << endl;
#endif
        }
        
        
#ifdef debug_spliced_surject
        cerr << "computing optimal combination of sections" << endl;
#endif
        
        // compute edge weights for the constriction edges between sections
        vector<int32_t> section_edge_scores(constrictions.size());
        for (size_t i = 0; i < constrictions.size(); ++i) {
            size_t from = constriction_comps[constrictions[i]];
            size_t to = constriction_comps[constriction_targets[i]];
            
            bool strand = get_strand(comp_groups[from].front());
            
#ifdef debug_spliced_surject
            cerr << "constriction edge ref chunks: " << endl;
            cerr << "\tfrom: " << endl;
            cerr << "\t\t" << graph->get_id(graph->get_handle_of_step(ref_ranges[from].first)) << " " << graph->get_is_reverse(graph->get_handle_of_step(ref_ranges[from].first)) << endl;
            cerr << "\t\t" << graph->get_id(graph->get_handle_of_step(ref_ranges[from].second)) << " " << graph->get_is_reverse(graph->get_handle_of_step(ref_ranges[from].second)) << endl;
            cerr << "\tto: " << endl;
            cerr << "\t\t" << graph->get_id(graph->get_handle_of_step(ref_ranges[to].first)) << " " << graph->get_is_reverse(graph->get_handle_of_step(ref_ranges[to].first)) << endl;
            cerr << "\t\t" << graph->get_id(graph->get_handle_of_step(ref_ranges[to].second)) << " " << graph->get_is_reverse(graph->get_handle_of_step(ref_ranges[to].second)) << endl;
#endif
            
            // how much of the path did we skip
            size_t deletion_length;
            if (strand) {
                deletion_length = (graph->get_position_of_step(ref_ranges[from].second)
                                   - graph->get_position_of_step(ref_ranges[to].first)
                                   - graph->get_length(graph->get_handle_of_step(ref_ranges[to].first)));
            }
            else {
                deletion_length = (graph->get_position_of_step(ref_ranges[to].first)
                                   - graph->get_position_of_step(ref_ranges[from].second)
                                   - graph->get_length(graph->get_handle_of_step(ref_ranges[from].second)));
                                   
            }
            
            int32_t edge_score;
            if (deletion_length >= min_splice_length) {
                // this deletion is very long, so we will interpret it as a splice and not penalize for it
                edge_score = 0;
            }
            else {
                // short deletion, probably just a polymorphism
                edge_score = -(get_aligner()->gap_open +
                               (deletion_length > 1 ? (deletion_length - 1) * get_aligner()->gap_extension : 0));
            }
            section_edge_scores[i] = edge_score;
        }
        
        // now we find use dynamic programming to find the best alignment across chunks
        
        vector<int64_t> backpointer(sections.size(), -1);
        vector<int32_t> score_dp(sections.size(), numeric_limits<int32_t>::min());
        
        // initialize the scores
        for (size_t i = 0; i < sections.size(); ++i) {
            score_dp[i] = sections[i].score();
        }
        if (!allow_negative_scores) {
            // remove the initialization if it's not a source so we don't allow subpath local alignments
            for (size_t i = 0; i < constriction_targets.size(); ++i) {
                score_dp[constriction_comps[constriction_targets[i]]] = numeric_limits<int32_t>::min();
            }
        }
        
        // do the dynamic programming
        for (size_t i = 0; i < constrictions.size(); ++i) {
            size_t from = constriction_comps[constrictions[i]];
            size_t to = constriction_comps[constriction_targets[i]];
            
            int32_t extended_score = score_dp[from] + section_edge_scores[i] + sections[to].score();
            if (extended_score > score_dp[to]) {
                score_dp[to] = extended_score;
                backpointer[to] = from;
            }
        }
        
        // find the maximum
        vector<size_t> traceback(1, -1);
        int32_t max_score = numeric_limits<int32_t>::min();
        for (size_t i = 0; i < score_dp.size(); ++i) {
            if (score_dp[i] > max_score) {
                max_score = score_dp[i];
                traceback[0] = i;
            }
        }
        
        // follow the back pointers
        while (backpointer[traceback.back()] != -1) {
            traceback.push_back(backpointer[traceback.back()]);
        }
        
        
#ifdef debug_spliced_surject
        cerr << "combining " << traceback.size() << " sections into surjected alignment" << endl;
#endif
        
        // make an alignment to build out the path in
        Alignment surjected;
        transfer_read_metadata(source, surjected);
        surjected.set_mapping_quality(source.mapping_quality());
        
        Path* surj_path = surjected.mutable_path();
        for (int64_t i = traceback.size() - 1; i >= 0; --i) {
            
            size_t section_idx = traceback[i];
            const Path& copy_path = sections[section_idx].path();
            
#ifdef debug_spliced_surject
            cerr << "appending path section " << pb2json(copy_path) << endl;
#endif
            
            // we have to have some special logic for the first mapping in each new path
            if (surj_path->mapping_size() > 0) {
                
                // merge the mapping if possible
                
                const Mapping& copy_mapping = copy_path.mapping(0);
                Mapping* final_mapping = surj_path->mutable_mapping(surj_path->mapping_size() - 1);
                
                if (final_mapping->position().node_id() == copy_mapping.position().node_id()
                    && final_mapping->position().is_reverse() == copy_mapping.position().is_reverse()
                    && final_mapping->position().offset() + mapping_from_length(*final_mapping) == copy_mapping.position().offset()) {
                    // these can be merged
                    for (size_t j = 0; j < copy_mapping.edit_size(); ++j) {
                        *final_mapping->add_edit() = copy_mapping.edit(j);
                    }
                }
                else {
                    // can't be merged, just copy
                    Mapping* surj_mapping = surj_path->add_mapping();
                    *surj_mapping = copy_mapping;
                    surj_mapping->set_rank(surj_path->mapping_size());
                }
            }
            else {
                // the first mapping in the surjected path
                
                Mapping* surj_mapping = surj_path->add_mapping();
                surj_mapping->set_rank(1);
                
                const Mapping& copy_mapping = copy_path.mapping(0);
                *surj_mapping->mutable_position() = copy_mapping.position();
                
                // handle any front end soft clips
                if (read_ranges[section_idx].first != source.sequence().begin()) {
                    Edit* e = surj_mapping->add_edit();
                    e->set_to_length(read_ranges[section_idx].first - source.sequence().begin());
                    e->set_sequence(string(source.sequence().begin(), read_ranges[section_idx].first));
                }
                
                for (size_t j = 0; j < copy_mapping.edit_size(); ++j) {
                    *surj_mapping->add_edit() = copy_mapping.edit(j);
                }
            }
            
            // copy the remaining mappings
            for (size_t j = 1; j < copy_path.mapping_size(); ++j) {
                Mapping* surj_mapping = surj_path->add_mapping();
                *surj_mapping = copy_path.mapping(j);
                surj_mapping->set_rank(surj_path->mapping_size());
            }
        }
        
        // handle any back end soft clips
        if (read_ranges[traceback[0]].second != source.sequence().end()) {
            Edit* e = surj_path->mutable_mapping(surj_path->mapping_size() - 1)->add_edit();
            e->set_to_length(source.sequence().end() - read_ranges[traceback[0]].second);
            e->set_sequence(string(read_ranges[traceback[0]].second, source.sequence().end()));
        }
        
        // set the score to the combined value
        surjected.set_score(max_score);
        
#ifdef debug_spliced_surject
        cerr << "final spliced surjection " << pb2json(surjected) << endl;
#endif
        
        return move(surjected);
    }

    Alignment Surjector::realigning_surject(const PathPositionHandleGraph* path_position_graph, const Alignment& source,
                                            const path_handle_t& path_handle, const vector<path_chunk_t>& path_chunks,
                                            bool allow_negative_scores) const {
        
#ifdef debug_anchored_surject
        cerr << "using overlap chunks on path " << graph->get_path_name(path_handle) << ", performing realigning surjection" << endl;
        cerr << "chunks:" << endl;
        for (size_t i = 0; i < path_chunks.size(); ++i) {
            cerr << "\t" << string(path_chunks[i].first.first, path_chunks[i].first.second) << ", " << pb2json(path_chunks[i].second) << endl;
        }
#endif
        
        // find the interval of the ref path we need to consider
        pair<size_t, size_t> ref_path_interval = compute_path_interval(path_position_graph, source, path_handle,
                                                                       path_chunks);
        
#ifdef debug_anchored_surject
        cerr << "final path interval is " << ref_path_interval.first << ":" << ref_path_interval.second << endl;
#endif
        
        // get the path graph corresponding to this interval
        bdsg::HashGraph path_graph;
        unordered_map<id_t, pair<id_t, bool>> path_trans = extract_linearized_path_graph(path_position_graph, &path_graph, path_handle,
                                                                                         ref_path_interval.first, ref_path_interval.second);
        
        // split it into a forward and reverse strand
        StrandSplitGraph split_path_graph(&path_graph);
        
        // make a translator down to the original graph
        unordered_map<id_t, pair<id_t, bool>> node_trans;
        split_path_graph.for_each_handle([&](const handle_t& handle) {
            handle_t underlying = split_path_graph.get_underlying_handle(handle);
            const pair<id_t, bool>& original = path_trans[path_graph.get_id(underlying)];
            node_trans[split_path_graph.get_id(handle)] = make_pair(original.first,
                                                                    original.second != path_graph.get_is_reverse(underlying));
        });
        
#ifdef debug_anchored_surject
        cerr << "made split, linearized path graph with " << split_path_graph.get_node_count() << " nodes" << endl;
#endif
        
        // compute the connectivity between the path chunks
        MultipathAlignmentGraph mp_aln_graph(split_path_graph, path_chunks, source, node_trans);
        
        // we don't overlap this reference path at all or we filtered out all of the path chunks, so just make a sentinel
        if (mp_aln_graph.empty()) {
            return move(make_null_alignment(source));
        }
        
        // TODO: is this necessary in a linear graph?
        vector<size_t> topological_order;
        mp_aln_graph.topological_sort(topological_order);
        mp_aln_graph.remove_transitive_edges(topological_order);
        
        // align the intervening segments and store the result in a multipath alignment
        MultipathAlignment mp_aln;
        mp_aln_graph.align(source, split_path_graph, get_aligner(), false, 1, false, 1, mp_aln, allow_negative_scores);
        topologically_order_subpaths(mp_aln);
        
        for (size_t i = 0; i < mp_aln.subpath_size(); i++) {
            // translate back into the original ID space
            translate_oriented_node_ids(*mp_aln.mutable_subpath(i)->mutable_path(), node_trans);
        }
        
        // identify the source subpaths (necessary for subpath-global optimal alignment algorithm)
        identify_start_subpaths(mp_aln);
        
#ifdef debug_anchored_surject
        cerr << "made multipath alignment " << pb2json(mp_aln) << endl;
#endif
        
#ifdef debug_validate_anchored_multipath_alignment
        if (!validate_multipath_alignment(mp_aln, *graph)) {
            cerr << "WARNING: multipath alignment for surjection of " << source.name() << " failed to validate" << endl;
        }
#endif
        
        // concatenate the subpaths either locally or globally, depending on whether we're
        // allowing negative scores
        Alignment surjected;
        optimal_alignment(mp_aln, surjected, allow_negative_scores);
        
        // transfer applicable metadata
        surjected.set_mapping_quality(source.mapping_quality());
        if (source.has_fragment_next()) {
            *surjected.mutable_fragment_next() = source.fragment_next();
        }
        if (source.has_fragment_prev()) {
            *surjected.mutable_fragment_prev() = source.fragment_prev();
        }
        
#ifdef debug_anchored_surject
        cerr << "concatenated and translated alignment " << pb2json(surjected) << endl;
#endif
        
        return surjected;
    }
    
    unordered_map<path_handle_t, pair<vector<Surjector::path_chunk_t>, vector<pair<step_handle_t, step_handle_t>>>>
    Surjector::extract_overlapping_paths(const PathPositionHandleGraph* graph, const Alignment& source,
                                         const unordered_set<path_handle_t>& surjection_paths) const {
        
        unordered_map<path_handle_t, pair<vector<path_chunk_t>, vector<pair<step_handle_t, step_handle_t>>>> to_return;
        
        const Path& path = source.path();
        
        // for each path that we're extending, the previous step and the strand we were at on it
        // mapped to the index of that path chunk in the path's vector
        unordered_map<pair<step_handle_t, bool>, size_t> extending_steps;
        int64_t through_to_length = 0;
        
        for (size_t i = 0; i < path.mapping_size(); i++) {
            
            int64_t before_to_length = through_to_length;
            through_to_length += mapping_to_length(path.mapping(i));
            
            const Position& pos = path.mapping(i).position();
            handle_t handle = graph->get_handle(pos.node_id(), pos.is_reverse());
            
#ifdef debug_anchored_surject
            cerr << "looking for paths on mapping " << i << " at position " << make_pos_t(pos) << endl;
#endif
            
            unordered_map<pair<step_handle_t, bool>, size_t> next_extending_steps;
            
            for (const step_handle_t& step : graph->steps_of_handle(handle)) {
                
#ifdef debug_anchored_surject
                cerr << "found a step on " << graph->get_path_name(graph->get_path_handle_of_step(step)) << endl;
#endif
                
                path_handle_t path_handle = graph->get_path_handle_of_step(step);
                if (!surjection_paths.count(path_handle)) {
                    // we are not surjecting onto this path
#ifdef debug_anchored_surject
                    cerr << "not surjecting to this path, skipping" << endl;
#endif
                    continue;
                }
                
                bool path_strand = graph->get_is_reverse(handle) != graph->get_is_reverse(graph->get_handle_of_step(step));
                
                step_handle_t prev_step = path_strand ? graph->get_next_step(step) : graph->get_previous_step(step);
                
#ifdef debug_anchored_surject
                cerr << "path strand is " << (path_strand ? "rev" : "fwd") << ", prev step is ";
                if (prev_step == graph->path_end(path_handle)) {
                    cerr << " path end";
                }
                else if (prev_step == graph->path_front_end(path_handle)) {
                    cerr << " path front end";
                }
                else {
                    cerr << graph->get_id(graph->get_handle_of_step(prev_step)) << (graph->get_is_reverse(graph->get_handle_of_step(prev_step)) ? "-" : "+");
                }
                cerr << endl;
                cerr << "possible extensions from: " << endl;
                for (const auto& record : extending_steps) {
                    cerr << "\t" << graph->get_id(graph->get_handle_of_step(record.first.first)) << (graph->get_is_reverse(graph->get_handle_of_step(record.first.first)) ? "-" : "+") << " on " << graph->get_path_name(graph->get_path_handle_of_step(record.first.first)) << " " << (record.first.second ? "rev" : "fwd") << endl;
                }
#endif
                
                if (extending_steps.count(make_pair(prev_step, path_strand))) {
                    // we are extending from the previous step, so we continue with the extension
                    
                    auto& path_chunks = to_return[graph->get_path_handle_of_step(step)];
                    size_t chunk_idx = extending_steps[make_pair(prev_step, path_strand)];
                    auto& aln_chunk = path_chunks.first[chunk_idx];
                    auto& ref_chunk = path_chunks.second[chunk_idx];
                    
                    // extend the range of the path on the reference
                    ref_chunk.second = step;
                    
                    // move the end of the sequence out
                    aln_chunk.first.second = source.sequence().begin() + through_to_length;
                    Mapping* mapping = aln_chunk.second.add_mapping();
                    // add this mapping
                    *mapping = path.mapping(i);
                    mapping->set_rank(aln_chunk.second.mapping(aln_chunk.second.mapping_size() - 2).rank() + 1);
                    
                    // in the next iteration, this step should point into the chunk it just extended
                    next_extending_steps[make_pair(step, path_strand)] = extending_steps[make_pair(prev_step, path_strand)];
                }
                else {
                    // this step does not extend a previous step, so we start a new chunk
                    auto& path_chunks = to_return[graph->get_path_handle_of_step(step)];
                    path_chunks.first.emplace_back();
                    path_chunks.second.emplace_back();
                    auto& aln_chunk = path_chunks.first.back();
                    auto& ref_chunk = path_chunks.second.back();
                    
                    // init the ref interval with the interval along the embedded path
                    ref_chunk.first = step;
                    ref_chunk.second = step;
                    
                    // init the new chunk with the sequence interval
                    aln_chunk.first.first = source.sequence().begin() + before_to_length;
                    aln_chunk.first.second = source.sequence().begin() + through_to_length;
                    
                    // and with the first mapping
                    Mapping* mapping = aln_chunk.second.add_mapping();
                    *mapping = path.mapping(i);
                    mapping->set_rank(1);
                    
                    // keep track of where this chunk is in the vector and which step it came from
                    // for the next iteration
                    next_extending_steps[make_pair(step, path_strand)] = path_chunks.first.size() - 1;
                }
            }
            
            // we've finished extending the steps from the previous mapping, so we replace them
            // with the steps we found in this iteration that we want to extend on the next one
            extending_steps = next_extending_steps;
        }
        
        return to_return;
    }
    
    pair<size_t, size_t>
    Surjector::compute_path_interval(const PathPositionHandleGraph* graph, const Alignment& source, path_handle_t path_handle,
                                     const vector<path_chunk_t>& path_chunks) const {
        
        pair<size_t, size_t> interval(numeric_limits<size_t>::max(), numeric_limits<size_t>::min());
        
        for (const auto& path_chunk : path_chunks) {
            
            size_t path_length = graph->get_path_length(path_handle);
            
            size_t left_overhang = (get_aligner()->longest_detectable_gap(source, path_chunk.first.first)
                                    + (path_chunk.first.first - source.sequence().begin()));
            
            size_t right_overhang = (get_aligner()->longest_detectable_gap(source, path_chunk.first.second)
                                     + (source.sequence().end() - path_chunk.first.second));
            
            const Position& first_pos = path_chunk.second.mapping(0).position();
            handle_t first_handle = graph->get_handle(first_pos.node_id(), first_pos.is_reverse());
            for (const step_handle_t& step : graph->steps_of_handle(first_handle)) {
                
                if (graph->get_path_handle_of_step(step) != path_handle) {
                    // this step isn't on the path we're considering
                    continue;
                }
                
                if (first_pos.is_reverse() != graph->get_is_reverse(graph->get_handle_of_step(step))) {
                    size_t path_offset = graph->get_position_of_step(step) + graph->get_length(first_handle) - first_pos.offset();
                    interval.second = max(interval.second, min(path_offset + left_overhang, path_length - 1));
                }
                else {
                    size_t path_offset = graph->get_position_of_step(step) + first_pos.offset();
                    if (left_overhang > path_offset) {
                        // avoid underflow
                        interval.first = 0;
                    }
                    else {
                        interval.first = min(interval.first, path_offset - left_overhang);
                    }
                }
            }
            
            const Mapping& final_mapping = path_chunk.second.mapping(path_chunk.second.mapping_size() - 1);
            const Position& final_pos = final_mapping.position();
            handle_t final_handle = graph->get_handle(final_pos.node_id(), final_pos.is_reverse());
            for (const step_handle_t& step : graph->steps_of_handle(final_handle)) {
                
                if (graph->get_path_handle_of_step(step) != path_handle) {
                    // this step isn't on the path we're considering
                    continue;
                }
                
                if (final_pos.is_reverse() != graph->get_is_reverse(graph->get_handle_of_step(step))) {
                    size_t path_offset = graph->get_position_of_step(step) + graph->get_length(final_handle) - final_pos.offset() - mapping_from_length(final_mapping);
                    if (right_overhang > path_offset) {
                        // avoid underflow
                        interval.first = 0;
                    }
                    else {
                        interval.first = min(interval.first, path_offset - right_overhang);
                    }
                }
                else {
                    size_t path_offset = graph->get_position_of_step(step) + first_pos.offset() + mapping_to_length(final_mapping);
                    interval.second = max(interval.second, min(path_offset + right_overhang, path_length - 1));
                }
            }
        }
        
        return interval;
    }
    
    unordered_map<id_t, pair<id_t, bool>>
    Surjector::extract_linearized_path_graph(const PathPositionHandleGraph* graph, MutableHandleGraph* into,
                                             path_handle_t path_handle, size_t first, size_t last) const {
        
        // TODO: we need better semantics than an unsigned interval for surjecting to circular paths
        
#ifdef debug_anchored_surject
        cerr << "extracting path graph for position interval " << first << ":" << last << " in path of length " << graph->get_path_length(path_handle) << endl;
#endif
        
        unordered_map<id_t, pair<id_t, bool>> node_trans;
        
        step_handle_t begin = graph->get_step_at_position(path_handle, first);
        step_handle_t end = graph->get_step_at_position(path_handle, last);
        
        if (graph->get_position_of_step(end) < last && end != graph->path_end(path_handle)) {
            // we actually want part of this step too, so we use the next one as the end iterator
            end = graph->get_next_step(end);
        }
        
        handle_t prev_node;
        for (step_handle_t step = begin; step != end; step = graph->get_next_step(step)) {
            // copy the node with the local orientation now forward
            handle_t copying = graph->get_handle_of_step(step);
            handle_t node_here = into->create_handle(graph->get_sequence(copying));
            
            if (step != begin) {
                // add an edge from the previous node
                into->create_edge(prev_node, node_here);
            }
            
            // record the translation
            node_trans[into->get_id(node_here)] = pair<id_t, bool>(graph->get_id(copying),
                                                                   graph->get_is_reverse(copying));
            
            prev_node = node_here;
        }
        
        return node_trans;
    }

    void Surjector::set_path_position_inexact(const PathPositionHandleGraph* graph, const Alignment& surjected,
                                              path_handle_t best_path_handle,
                                              string& path_name_out, int64_t& path_pos_out, bool& path_rev_out) const {
        
        const Path& path = surjected.path();
        
        if (path.mapping_size() == 0){
            // hack: we want a 0 position once we convert this to 1-based indexes
            path_pos_out = -1;
            path_rev_out = false;
            path_name_out = "";
            return;
        }
        
        int64_t fwd_pos = numeric_limits<int64_t>::max();
        int64_t rev_pos = numeric_limits<int64_t>::max();
        for (const Mapping& mapping : path.mapping()) {
            
            const Position& pos = mapping.position();
            for (const step_handle_t& step : graph->steps_of_handle(graph->get_handle(pos.node_id()))) {
                if (graph->get_path_handle_of_step(step) != best_path_handle) {
                    continue;
                }
                
                if (graph->get_is_reverse(graph->get_handle_of_step(step)) == pos.is_reverse()) {
                    // forward strand of the path
                    int64_t path_pos = graph->get_position_of_step(step) + pos.offset();
                    fwd_pos = min(path_pos, fwd_pos);
                }
                else {
                    int64_t path_pos = (graph->get_position_of_step(step)
                                        + graph->get_length(graph->get_handle_of_step(step))
                                        - pos.offset()
                                        - mapping_from_length(mapping));
                    rev_pos = min(path_pos, fwd_pos);
                }
            
            }
        }
        
        if (fwd_pos == numeric_limits<int64_t>::max() && rev_pos == numeric_limits<int64_t>::max()) {
            throw runtime_error("error:[Surjector] could not identify path position of surjected alignment " + surjected.name());
        }
        else if (fwd_pos <= rev_pos) {
            path_pos_out = fwd_pos;
            path_rev_out = false;
        }
        else {
            path_pos_out = rev_pos;
            path_rev_out = true;
        }
        path_name_out = graph->get_path_name(best_path_handle);
    }
    
    void Surjector::set_path_position(const PathPositionHandleGraph* graph, const Alignment& surjected,
                                      path_handle_t best_path_handle,
                                      string& path_name_out, int64_t& path_pos_out, bool& path_rev_out) const {

        const Path& path = surjected.path();
        
        if (path.mapping_size() == 0){
            // hack: we want a 0 position once we convert this to 1-based indexes
            path_pos_out = -1;
            path_rev_out = false;
            path_name_out = "";
            return;
        }
        
#ifdef debug_anchored_surject
        cerr << "choosing a path position for surjected alignment on path " << graph->get_path_name(best_path_handle) << endl;
#endif
        
        const Position& start_pos = path.mapping(0).position();
        handle_t start_handle = graph->get_handle(start_pos.node_id(), start_pos.is_reverse());
        // TODO: depending on path coverage, it could be inefficient to iterate over all steps on the handle
        for (const step_handle_t& step : graph->steps_of_handle(start_handle)) {
#ifdef debug_anchored_surject
            cerr << "found step on " << graph->get_path_name(graph->get_path_handle_of_step(step)) << endl;
#endif
            if (graph->get_path_handle_of_step(step) != best_path_handle) {
                // this is not the path we surjected onto
#ifdef debug_anchored_surject
                cerr << "wrong path, skipping" << endl;
#endif
                continue;
            }
            
            bool strand = graph->get_is_reverse(graph->get_handle_of_step(step)) != start_pos.is_reverse();
            
            // walk this whole thing out to ensure that this is the correct match
            bool match = true;
            step_handle_t path_step = strand ? graph->get_previous_step(step) : graph->get_next_step(step);
            for (size_t i = 1; i < path.mapping_size(); ++i) {
                
                if (path_step == graph->path_end(best_path_handle) ||
                    path_step == graph->path_front_end(best_path_handle)) {
                    // we've gone off the end of the path
#ifdef debug_anchored_surject
                    cerr << "encountered end of path at mapping index " << i << ", skipping" << endl;
#endif
                    match = false;
                    break;
                }
                
                const Position& pos = path.mapping(i).position();
#ifdef debug_anchored_surject
                cerr << "at mapping pos " << make_pos_t(pos) << " and path pos " << graph->get_id(graph->get_handle_of_step(path_step)) << (graph->get_is_reverse(graph->get_handle_of_step(path_step)) ? "-" : "+") << endl;
#endif
                
                if (pos.node_id() != graph->get_id(graph->get_handle_of_step(path_step)) ||
                    (pos.is_reverse() != graph->get_is_reverse(graph->get_handle_of_step(path_step)) != strand)) {
                    // the path here doesn't match the surjected path
#ifdef debug_anchored_surject
                    cerr << "encountered mismatch at mapping index " << i << ", skipping" << endl;
#endif
                    match = false;
                    break;
                }
                
                path_step = strand ? graph->get_previous_step(path_step) : graph->get_next_step(path_step);
            }
            
            if (match) {
                // we found a match, so set the position values and return
                path_name_out = graph->get_path_name(best_path_handle);
                path_rev_out = strand;
                size_t start_offset = graph->get_position_of_step(step);
                if (strand) {
                    path_pos_out = start_offset + graph->get_length(start_handle) - start_pos.offset() - path_from_length(path);
                }
                else {
                    path_pos_out = start_offset + start_pos.offset();
                }
                return;
            }
        }
        
        // we ran through all of the occurrences without finding a full match...
        throw runtime_error("error:[Surjector] could not identify path position of surjected alignment " + surjected.name());
    }
    
    Alignment Surjector::make_null_alignment(const Alignment& source) {
        Alignment null;
        null.set_name(source.name());
        null.set_sequence(source.sequence());
        null.set_quality(source.quality());
        null.set_read_group(source.read_group());
        null.set_sample_name(source.sample_name());
        null.set_is_secondary(source.is_secondary());
        if (source.has_fragment_next()) {
            *null.mutable_fragment_next() = source.fragment_next();
        }
        if (source.has_fragment_prev()) {
            *null.mutable_fragment_prev() = source.fragment_prev();
        }
        return null;
    }
}


