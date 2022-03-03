/**
 * \file
 * surjector.cpp: implements a class that surjects alignments onto paths
 */

#include "surjector.hpp"

//#define debug_spliced_surject
//#define debug_anchored_surject
//#define debug_multipath_surject
//#define debug_constrictions
//#define debug_prune_unconnectable
//#define debug_filter_paths
//#define debug_validate_anchored_multipath_alignment
//#define debug_always_warn_on_too_long

namespace vg {

using namespace std;
    
    Surjector::Surjector(const PathPositionHandleGraph* graph) : graph(graph) {
        if (!graph) {
            cerr << "error:[Surjector] Failed to provide an graph to the Surjector" << endl;
        }
    }
    
    Alignment Surjector::surject(const Alignment& source, const unordered_set<path_handle_t>& paths,
                                 bool allow_negative_scores, bool preserve_deletions) const {
    
        // Allocate the annotation info
        string path_name_out;
        int64_t path_pos_out;
        bool path_rev_out;
        
        // Do the surjection
        Alignment surjected = surject(source, paths, path_name_out, path_pos_out, path_rev_out, allow_negative_scores, preserve_deletions);
        
        // Pack all the info into the refpos field
        surjected.clear_refpos();
        auto* pos = surjected.add_refpos();
        pos->set_name(path_name_out);
        pos->set_offset(path_pos_out);
        pos->set_is_reverse(path_rev_out);
    
        return surjected;
    }

    Alignment Surjector::surject(const Alignment& source, const unordered_set<path_handle_t>& paths, string& path_name_out,
                                 int64_t& path_pos_out, bool& path_rev_out, bool allow_negative_scores,
                                 bool preserve_deletions) const {
        Alignment surjected;
        surject_internal(&source, nullptr, &surjected, nullptr, paths, path_name_out, path_pos_out,
                         path_rev_out, allow_negative_scores, preserve_deletions);
        return surjected;
    }

    multipath_alignment_t Surjector::surject(const multipath_alignment_t& source, const unordered_set<path_handle_t>& paths,
                                             string& path_name_out, int64_t& path_pos_out, bool& path_rev_out,
                                             bool allow_negative_scores, bool preserve_deletions) const {

        multipath_alignment_t surjected;
        surject_internal(nullptr, &source, nullptr, &surjected, paths, path_name_out, path_pos_out,
                         path_rev_out, allow_negative_scores, preserve_deletions);
        return surjected;
    }
    
    void Surjector::surject_internal(const Alignment* source_aln, const multipath_alignment_t* source_mp_aln,
                                     Alignment* aln_out, multipath_alignment_t* mp_aln_out,
                                     const unordered_set<path_handle_t>& paths,
                                     string& path_name_out, int64_t& path_pos_out, bool& path_rev_out,
                                     bool allow_negative_scores, bool preserve_deletions) const {

        
        // we need one and only one data type: Alignment or multipath_alignment_t
        assert(!(source_aln && source_mp_aln));
        assert((source_aln && aln_out) || (source_mp_aln && mp_aln_out));
        
#ifdef debug_anchored_surject
        cerr << "surjecting alignment: ";
        if (source_mp_aln) {
            cerr << debug_string(*source_mp_aln);
        }
        else {
            cerr << pb2json(*source_aln);
        }
        cerr << " onto paths ";
        for (const path_handle_t& path : paths) {
            cerr << graph->get_path_name(path) << " ";
        }
        cerr << endl;
#endif
                
        if (source_aln && source_aln->path().mapping_size() != 0) {
            // The read is mapped. Check the input alignment for basic
            // consistency. If the sequence and the graph path don't agree
            // about the read length, something is very wrong with the input.
            size_t source_to_length = path_to_length(source_aln->path());
            if (source_aln->sequence().size() != source_to_length) {
                cerr << "error[Surjector::surject]: read " << source_aln->name() << " has "
                << source_aln->sequence().size() << " sequence bases but an input alignment that aligns "
                << source_to_length << " bases instead. This is invalid and uninterpretable; check your mapper." << endl;
                cerr << "error[Surjector::surject]: offending alignemnt: " << pb2json(*source_aln) << endl; 
                exit(1);
            }
        }
        
        // make an overlay that will memoize the results of some expensive XG operations
        MemoizingGraph memoizing_graph(graph);
        
        // get the chunks of the aligned path that overlap the ref path
        unordered_map<pair<path_handle_t, bool>, vector<tuple<size_t, size_t, int32_t>>> connections;
        auto path_overlapping_anchors = source_aln ? extract_overlapping_paths(&memoizing_graph, *source_aln, paths)
                                                   : extract_overlapping_paths(&memoizing_graph, *source_mp_aln,
                                                                               paths, connections);
        
        if (source_mp_aln) {
            // the multipath alignment anchor algorithm can produce redundant paths if
            // the alignment's graph is not parsimonious, so we filter the shorter ones out
            for (pair<const pair<path_handle_t, bool>, pair<vector<path_chunk_t>, vector<pair<step_handle_t, step_handle_t>>>>& path_chunk_record : path_overlapping_anchors) {
                filter_redundant_path_chunks(path_chunk_record.second.first, path_chunk_record.second.second,
                                             connections[path_chunk_record.first]);
            }
        }
        
#ifdef debug_anchored_surject
        cerr << "got path overlapping segments" << endl;
        for (const auto& surjection_record : path_overlapping_anchors) {
            cerr << "path " << graph->get_path_name(surjection_record.first.first) << ", rev? " << surjection_record.first.second << endl;
            
            for (size_t i = 0; i < surjection_record.second.first.size(); ++i) {
                auto& anchor = surjection_record.second.first[i];
                if (source_aln) {
                    cerr << "\tread[" << (anchor.first.first - source_aln->sequence().begin()) << ":" << (anchor.first.second - source_aln->sequence().begin()) << "] : ";
                }
                else {
                    cerr << "\tread[" << (anchor.first.first - source_mp_aln->sequence().begin()) << ":" << (anchor.first.second - source_mp_aln->sequence().begin()) << "] : ";
                }
                for (auto iter = anchor.first.first; iter != anchor.first.second; iter++) {
                    cerr << *iter;
                }
                cerr << endl;
                cerr << "\tpath interval " << graph->get_position_of_step(surjection_record.second.second[i].first) << " - " << graph->get_position_of_step(surjection_record.second.second[i].second) << endl;
                cerr << "\t" << pb2json(anchor.second) << endl;
            }
            if (connections.count(surjection_record.first)) {
                cerr << "\tconnections" << endl;
                for (const auto& connection : connections[surjection_record.first]) {
                    cerr << "\t\t" << get<0>(connection) << " -> " << get<1>(connection) << " (" << get<2>(connection) << ")" << endl;
                }
            }
        }
#endif
        
        if (prune_suspicious_anchors) {
            // we want to remove anchors that can be error-prone: short anchors in the tails and anchors in
            // low complexity sequences
            for (auto it = path_overlapping_anchors.begin(); it != path_overlapping_anchors.end(); ++it) {
                auto& path_chunks = it->second.first;
                auto& step_ranges = it->second.second;
                vector<bool> keep(path_chunks.size(), true);
                for (int i = 0; i < path_chunks.size(); ++i) {
                    auto& chunk = path_chunks[i];
                    if (((i == 0 || i + 1 == path_chunks.size()) && path_chunks.size() != 1)
                        && path_from_length(chunk.second) <= max_tail_anchor_prune &&
                        chunk.first.second - chunk.first.first <= max_tail_anchor_prune) {
#ifdef debug_anchored_surject
                        cerr << "anchor " << i << " pruned for being a short tail" << endl;
#endif
                        // this is a short anchor on one of the tails
                        keep[i] = false;
                        continue;
                    }
                    SeqComplexity<6> complexity(chunk.first.first, chunk.first.second);
                    for (int order = 1; order <= 6; ++order) {
                        if (complexity.p_value(order) < low_complexity_p_value) {
#ifdef debug_anchored_surject
                            cerr << "anchor " << i << " pruned being low complexity at order " << order << " with p-value " << complexity.p_value(order) << " and repetitive fraction " << complexity.repetitiveness(order) << endl;
#endif
                            // the sequences is repetitive at this order
                            keep[i] = false;
                            break;
                        }
                    }
                }
                // make sure we didn't flag all of the anchors for removal
                bool keep_any = false;
                for (bool b : keep) {
                    keep_any = keep_any || b;
                }
                if (!keep_any) {
                    // we filtered out all of the anchors, choose the longest one to keep
                    // even though it failed the filter
                    int64_t max_idx = -1;
                    int64_t max_len = 0;
                    for (int i = 0; i < path_chunks.size(); ++i) {
                        int64_t len = path_chunks[i].first.second - path_chunks[i].first.first;
                        if (len > max_len) {
                            max_idx = i;
                            max_len = len;
                        }
                    }
                    if (max_idx >= 0) {
#ifdef debug_anchored_surject
                        cerr << "reversing decision to prune " << max_idx << endl;
#endif
                        keep[max_idx] = true;
                    }
                }
                // we're keeping at least one anchor, so we should be able to throw away the other ones
                int removed_so_far = 0;
                for (int i = 0; i < path_chunks.size(); ++i) {
                    if (!keep[i]) {
                        ++removed_so_far;
                    }
                    else if (removed_so_far) {
                        path_chunks[i - removed_so_far] = move(path_chunks[i]);
                        step_ranges[i - removed_so_far] = move(step_ranges[i]);
                    }
                }
                if (removed_so_far) {
                    path_chunks.resize(path_chunks.size() - removed_so_far);
                    step_ranges.resize(step_ranges.size() - removed_so_far);
                }
            }
        }
        
        
        // the surjected alignment for each path we overlapped
        unordered_map<pair<path_handle_t, bool>, pair<Alignment, pair<step_handle_t, step_handle_t>>> aln_surjections;
        unordered_map<pair<path_handle_t, bool>, pair<multipath_alignment_t, pair<step_handle_t, step_handle_t>>> mp_aln_surjections;
        for (pair<const pair<path_handle_t, bool>, pair<vector<path_chunk_t>, vector<pair<step_handle_t, step_handle_t>>>>& surj_record : path_overlapping_anchors) {
            
            // to hold the path interval that corresponds to the path we surject to
            pair<step_handle_t, step_handle_t> path_range;
            if (!preserve_deletions && source_aln) {
                // unspliced GAM -> GAM surjection
                auto surjection = realigning_surject(&memoizing_graph, *source_aln, surj_record.first.first, surj_record.first.second,
                                                     surj_record.second.first, path_range, allow_negative_scores,
                                                     false, false);
                if (surjection.path().mapping_size() != 0) {
                    aln_surjections[surj_record.first] = make_pair(move(surjection), path_range);
                }
            }
            else if (source_aln) {
                // spliced GAM -> GAM surjection
                auto surjection = spliced_surject(&memoizing_graph, source_aln->sequence(), source_aln->quality(),
                                                  source_aln->mapping_quality(), surj_record.first.first, surj_record.first.second,
                                                  surj_record.second.first, surj_record.second.second,
                                                  connections[surj_record.first], path_range,
                                                  allow_negative_scores, preserve_deletions);
                if (surjection.subpath_size() != 0) {
                    // this internal method is written for multipath alignments, so we need to convert to standard alignments
                    aln_surjections[surj_record.first] = make_pair(Alignment(), path_range);
                    auto& surjected_aln = aln_surjections[surj_record.first].first;
                    optimal_alignment(surjection, surjected_aln, allow_negative_scores);
                    transfer_read_metadata(*source_aln, surjected_aln);
                }
            }
            else {
                // surjecting a multipath alignment (they always use the spliced pathway even if not
                // doing spliced alignment)
                auto surjection = spliced_surject(&memoizing_graph, source_mp_aln->sequence(),
                                                  source_mp_aln->quality(), source_mp_aln->mapping_quality(),
                                                  surj_record.first.first, surj_record.first.second,
                                                  surj_record.second.first, surj_record.second.second,
                                                  connections[surj_record.first], path_range,
                                                  allow_negative_scores, preserve_deletions);
                if (surjection.subpath_size() != 0) {
                    // the surjection was a success
                    
                    // copy over annotations
                    // TODO: also redundantly copies over sequence and quality
                    transfer_read_metadata(*source_mp_aln, surjection);
                    
                    // record the result for this path
                    mp_aln_surjections[surj_record.first] = make_pair(move(surjection), path_range);
                }
            }
        }
        
        // in case we didn't overlap any paths, add a sentinel so the following code still executes correctly
        if (aln_surjections.empty() && mp_aln_surjections.empty()) {
            // this surjection didn't get aligned
            path_name_out = "";
            path_rev_out = false;
            path_pos_out = -1;
            if (source_mp_aln) {
                *mp_aln_out = make_null_mp_alignment(source_mp_aln->sequence(), source_mp_aln->quality());
            }
            else {
                *aln_out = make_null_alignment(*source_aln);
            }
            return;
        }
    
        // choose which path surjection was best
        pair<path_handle_t, bool> best_path_strand;
        int32_t score = numeric_limits<int32_t>::min();
        for (const auto& surjection : aln_surjections) {
            if (surjection.second.first.score() >= score) {
#ifdef debug_anchored_surject
                cerr << "surjection against path " << graph->get_path_name(surjection.first.first) << " strand " << surjection.first.second << " achieves highest score of " << surjection.second.first.score() << ": " << pb2json(surjection.second.first) << endl;
#endif
                score = surjection.second.first.score();
                best_path_strand = surjection.first;
            }
        }
        for (const auto& surjection : mp_aln_surjections) {

            int32_t surj_score = optimal_alignment_score(surjection.second.first, allow_negative_scores);
            if (surj_score >= score) {
#ifdef debug_anchored_surject
                cerr << "surjection against path " << graph->get_path_name(surjection.first.first) << " strand " << surjection.first.second << " achieves highest score of " << surj_score << ": " << debug_string(surjection.second.first) << endl;
#endif
                score = surj_score;
                best_path_strand = surjection.first;
            }
        }
                
        // find the position along the path
        
        // retrieve the first/last positions of the best alignment and the corresponding
        // path range
        pair<step_handle_t, step_handle_t> path_range;
        pos_t initial_pos, final_pos;
        if (aln_out) {
            auto& surjection = aln_surjections[best_path_strand];
            initial_pos = initial_position(surjection.first.path());
            final_pos = final_position(surjection.first.path());
            path_range = surjection.second;
            *aln_out = move(surjection.first);
        }
        else {
            auto& surjection = mp_aln_surjections[best_path_strand];
            initial_pos = initial_position(surjection.first.subpath().front().path());
            final_pos = final_position(surjection.first.subpath().back().path());
            path_range = surjection.second;
            *mp_aln_out = move(surjection.first);
        }
        
        // use this info to set the path position
        set_path_position(&memoizing_graph, initial_pos, final_pos, path_range.first, path_range.second,
                          best_path_strand.second, path_name_out, path_pos_out, path_rev_out);
        
        
#ifdef debug_anchored_surject
        cerr << "chose path " << path_name_out << " at position " << path_pos_out << (path_rev_out ? "-" : "+") << endl;
#endif
        
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

    vector<vector<size_t>> Surjector::remove_dominated_chunks(const string& src_sequence,
                                                              const vector<vector<size_t>>& adj,
                                                              vector<path_chunk_t>& path_chunks,
                                                              vector<pair<step_handle_t, step_handle_t>>& ref_chunks,
                                                              vector<tuple<size_t, size_t, int32_t>>& connections) const {
        
        // this is an easy way to ensure that all adjacency lists are ordered by index
        auto rev_adj = reverse_adjacencies(adj);
        auto fwd_adj = reverse_adjacencies(rev_adj);
        
        map<pair<vector<size_t>, vector<size_t>>, vector<size_t>> neighbor_groups;
        for (size_t i = 0; i < adj.size(); ++i) {
            neighbor_groups[make_pair(rev_adj[i], fwd_adj[i])].push_back(i);
        }
        
#ifdef debug_spliced_surject
        cerr << "neighbor groups:" << endl;
        for (const auto& group : neighbor_groups) {
            cerr << "(";
            for (size_t i = 0; i < group.first.first.size(); ++i) {
                if (i != 0){
                    cerr << ", ";
                }
                cerr << group.first.first[i];
            }
            cerr << ") (";
            for (size_t i = 0; i < group.first.second.size(); ++i) {
                if (i != 0){
                    cerr << ", ";
                }
                cerr << group.first.second[i];
            }
            cerr << ")" << endl;
            for (auto i : group.second) {
                cerr << "\t" << i << endl;
            }
        }
#endif
        
        vector<size_t> to_remove;
        for (const auto& group : neighbor_groups) {
            // only remove dominated chunks if they have the same, non-empty set of neighbors
            if (group.second.size() > 1 && (!group.first.first.empty() || !group.first.second.empty())) {
                vector<int64_t> total_lengths(group.second.size());
                int64_t max_total_length = 0;
                for (size_t i = 0; i < group.second.size(); ++i) {
                    auto& chunk = path_chunks[group.second[i]];
                    total_lengths[i] = path_from_length(chunk.second) + (chunk.first.second - chunk.first.first);
                    // don't count softclips
                    if (chunk.first.first == src_sequence.begin()) {
                        const auto& first_edit = *chunk.second.mapping().begin()->edit().begin();
                        if (first_edit.from_length() == 0) {
                            total_lengths[i] -= first_edit.to_length();
                        }
                    }
                    if (chunk.first.second == src_sequence.end()) {
                        const auto& last_edit = *chunk.second.mapping().rbegin()->edit().rbegin();
                        if (last_edit.from_length() == 0) {
                            total_lengths[i] -= last_edit.to_length();
                        }
                    }
                    max_total_length = max(total_lengths[i], max_total_length);
                }
                for (size_t i = 0; i < group.second.size(); ++i) {
                    if (total_lengths[i] < max_total_length - 2 * dominated_path_chunk_diff) {
                        to_remove.push_back(group.second[i]);
                    }
                }
            }
        }
        
        if (!to_remove.empty()) {
#ifdef debug_spliced_surject
            cerr << "marked for removal (unless has a connection):" << endl;
            for (auto i : to_remove) {
                cerr << "\t" << i << endl;
            }
#endif
            
            unordered_set<size_t> connected;
            for (const auto& connection : connections) {
                connected.insert(get<0>(connection));
                connected.insert(get<1>(connection));
            }
            
            vector<bool> do_remove(fwd_adj.size(), false);
            for (size_t i : to_remove) {
                if (!connected.count(i)) {
                    do_remove[i] = true;
                }
            }
            
            vector<size_t> removed_before(fwd_adj.size() + 1, 0);
            for (size_t i = 0; i < fwd_adj.size(); ++i) {
                if (do_remove[i]) {
                    removed_before[i + 1] = removed_before[i] + 1;
                }
                else {
                    if (removed_before[i]) {
                        fwd_adj[i - removed_before[i]] = move(fwd_adj[i]);
                        path_chunks[i - removed_before[i]] = move(path_chunks[i]);
                        ref_chunks[i - removed_before[i]] = move(ref_chunks[i]);
                    }
                    removed_before[i + 1] = removed_before[i];
                }
            }
            fwd_adj.resize(fwd_adj.size() - removed_before.back());
            path_chunks.resize(fwd_adj.size());
            ref_chunks.resize(fwd_adj.size());
            
            for (auto& adj_list : fwd_adj) {
                size_t removed_so_far = 0;
                for (size_t i = 0; i < adj_list.size(); ++i) {
                    if (do_remove[adj_list[i]]) {
                        ++removed_so_far;
                    }
                    else {
                        adj_list[i - removed_so_far] = adj_list[i] - removed_before[adj_list[i]];
                    }
                }
                adj_list.resize(adj_list.size() - removed_so_far);
            }
            
            for (auto& connection : connections) {
                get<0>(connection) -= removed_before[get<0>(connection)];
                get<1>(connection) -= removed_before[get<1>(connection)];
            }
        }
        return fwd_adj;
    }

    void Surjector::prune_unconnectable(vector<vector<size_t>>& adj,
                                        vector<vector<tuple<size_t, int32_t, bool>>>& splice_adj,
                                        vector<size_t>& component,
                                        vector<vector<size_t>>& comp_groups,
                                        vector<path_chunk_t>& path_chunks,
                                        vector<pair<step_handle_t, step_handle_t>>& ref_chunks) const {
        
        // record which path chunks and component groups have connection adjacencies
        vector<bool> has_connection_from(adj.size(), false), has_connection_to(adj.size(), false);
        vector<bool> comp_has_connection_from(comp_groups.size(), false), comp_has_connection_to(comp_groups.size(), false);
        for (size_t i = 0; i < splice_adj.size(); ++i) {
            for (auto& edge : splice_adj[i]) {
                if (get<2>(edge)) {
                    has_connection_from[i] = true;
                    comp_has_connection_from[component[i]] = true;
                    has_connection_to[get<0>(edge)] = true;
                    comp_has_connection_to[component[get<0>(edge)]] = true;
                }
            }
        }
        
        // record which path chunks are reachable from the connections adjacencies
        vector<bool> connects_forward(adj.size(), false), connects_backward(adj.size(), false);
        for (size_t i = 0; i < adj.size(); ++i) {
            size_t j = adj.size() - i - 1;
            connects_forward[i] = connects_forward[i]  || has_connection_to[i];
            connects_backward[j] = connects_backward[j] || has_connection_from[j];
            for (auto k : adj[i]) {
                connects_forward[k] = connects_forward[k] || connects_forward[i];
            }
            for (auto k : adj[j]) {
                connects_backward[j] = connects_backward[j] || connects_backward[k];
            }
        }
        
#ifdef debug_prune_unconnectable
        cerr << "connects forward:" << endl;
        for (size_t i = 0; i < adj.size(); ++i) {
            cerr << "\t" << i << ": " << connects_forward[i] << endl;
        }
        cerr << "connects backward:" << endl;
        for (size_t i = 0; i < adj.size(); ++i) {
            cerr << "\t" << i << ": " << connects_backward[i] << endl;
        }
#endif
        
        // mark path chunks for removal if a component has a connection but the
        // chunks don't occur on any interconnectino paths
        vector<unordered_set<size_t>> to_remove_by_group(comp_groups.size());
        for (size_t i = 0; i < adj.size(); ++i) {
            size_t grp = component[i];
            if ((comp_has_connection_to[grp] && !connects_forward[i])
                || (comp_has_connection_from[grp] && !connects_backward[i])) {
                to_remove_by_group[grp].insert(i);
#ifdef debug_prune_unconnectable
                cerr << "marking " << i << " for removal" << endl;
#endif
            }
        }
        
        // filter out the chunks
        vector<size_t> removed(adj.size() + 1, 0);
        for (size_t i = 0; i < adj.size(); ++i) {
            size_t grp = component[i];
            // remove if not on an inter-connection path, but don't remove an entire group
            // TODO: but how can we be sure to produce sensible results when an entire group
            // should be removed?
            if (to_remove_by_group[grp].count(i)
                && to_remove_by_group[grp].size() < comp_groups[grp].size()) {
#ifdef debug_prune_unconnectable
                cerr << "removing " << i << endl;
#endif
                removed[i + 1] = removed[i] + 1;
            }
            else {
                size_t removed_so_far = removed[i];
                removed[i + 1] = removed_so_far;
                if (removed_so_far) {
                    adj[i - removed_so_far] = move(adj[i]);
                    splice_adj[i - removed_so_far] = move(splice_adj[i]);
                    path_chunks[i - removed_so_far] = move(path_chunks[i]);
                    ref_chunks[i - removed_so_far] = move(ref_chunks[i]);
                    component[i - removed_so_far] = component[i];
                }
            }
        }
        
        if (removed.back()) {
            adj.resize(adj.size() - removed.back());
            splice_adj.resize(adj.size());
            path_chunks.resize(adj.size());
            ref_chunks.resize(adj.size());
            component.resize(adj.size());
            
            // rewire the adjacencies to the correct chunks
            for (auto& adj_list : adj) {
                size_t adj_removed = 0;
                for (size_t i = 0; i < adj_list.size(); ++i) {
                    if (removed[adj_list[i]] != removed[adj_list[i] + 1]) {
                        ++adj_removed;
                    }
                    else {
                        adj_list[i - adj_removed] = adj_list[i] - removed[adj_list[i]];
                    }
                }
                adj_list.resize(adj_list.size() - adj_removed);
            }
            
#ifdef debug_prune_unconnectable
            cerr << "rewired adj list:" << endl;
            for (size_t i = 0; i < adj.size(); ++i) {
                cerr << i << ":";
                for (auto j : adj[i]) {
                    cerr << " " << j;
                }
                cerr << endl;
            }
#endif
            
            // rewire the splice adjacencies to the correct chunks
            for (auto& splice_adj_list : splice_adj) {
                size_t adj_removed = 0;
                for (size_t i = 0; i < splice_adj_list.size(); ++i) {
                    auto& edge = splice_adj_list[i];
                    if (removed[get<0>(edge)] != removed[get<0>(edge) + 1]) {
                        ++adj_removed;
                    }
                    else {
                        splice_adj_list[i - adj_removed] = make_tuple(get<0>(edge) - removed[get<0>(edge)], get<1>(edge), get<2>(edge));
                    }
                }
                splice_adj_list.resize(splice_adj_list.size() - adj_removed);
            }
            
#ifdef debug_prune_unconnectable
            cerr << "rewired splice adj list:" << endl;
            for (size_t i = 0; i < splice_adj.size(); ++i) {
                cerr << i << ":";
                for (auto edge : splice_adj[i]) {
                    cerr << " (" << get<0>(edge) << " " << get<1>(edge) << " " << get<2>(edge) << ")";
                }
                cerr << endl;
            }
#endif
            
            // correct the indexes in the group
            for (auto& group : comp_groups) {
                size_t grp_removed = 0;
                for (size_t i = 0; i < group.size(); ++i) {
                    if (removed[group[i]] != removed[group[i] + 1]) {
                        ++grp_removed;
                    }
                    else {
                        group[i - grp_removed] = group[i] - removed[group[i]];
                    }
                }
                group.resize(group.size() - grp_removed);
            }
            
        }
    }

    vector<pair<vector<size_t>, vector<size_t>>> Surjector::find_constriction_bicliques(const vector<vector<size_t>>& adj,
                                                                                        const string& src_sequence,
                                                                                        const vector<path_chunk_t>& path_chunks,
                                                                                        const vector<tuple<size_t, size_t, int32_t>>& connections) const {
        
        auto connected_by_edge = [&](size_t i, size_t j) {
            const auto& final_mapping = *path_chunks[i].second.mapping().rbegin();
            const auto& final_position = final_mapping.position();
            const auto& initial_position = path_chunks[j].second.mapping().begin()->position();
            handle_t handle = graph->get_handle(final_position.node_id(),
                                                final_position.is_reverse());
            if (initial_position.offset() == 0
                && final_position.offset() + mapping_from_length(final_mapping) == graph->get_length(handle)
                && graph->has_edge(handle, graph->get_handle(initial_position.node_id(), initial_position.is_reverse()))) {
                return true;
            }
            else {
                return false;
            }
        };
        
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
        }
        
#ifdef debug_constrictions
        cerr << "forward counts" << endl;
        for (size_t i = 0; i < fwd.size(); ++i) {
            cerr << "\t" << i << ": " << fwd[i] << endl;
        }
        cerr << "backward counts" << endl;
        for (size_t i = 0; i < bwd.size(); ++i) {
            cerr << "\t" << i << ": " << bwd[i] << endl;
        }
#endif
        
        // divide up chunks by their component and their begin/end position on the read
        unordered_map<pair<size_t, int64_t>, vector<size_t>> chunks_by_begin, chunks_by_end;
        for (size_t i = 0; i < path_chunks.size(); ++i) {
            pair<size_t, int64_t> key_begin(comps[i], path_chunks[i].first.first - src_sequence.begin());
            pair<size_t, int64_t> key_end(comps[i], path_chunks[i].first.second - src_sequence.begin());
            chunks_by_begin[key_begin].push_back(i);
            chunks_by_end[key_end].push_back(i);
        }
        
        // reorganize the connections into an adjacency list
        unordered_map<size_t, unordered_set<size_t>> connection_adj;
        for (const auto& connection : connections) {
            connection_adj[get<0>(connection)].emplace(get<1>(connection));
        }
        
        // init return value
        vector<pair<vector<size_t>, vector<size_t>>> return_val;
        
        // subroutine to check if a candidate partition is a biclique deletion constriction
        auto test_candidate_biclique = [&](const unordered_set<size_t>& left_side,
                                           const unordered_set<size_t>& right_side,
                                           size_t comp_idx) {
            // record which pairs have a connection
            bool incompatible = false;
            unordered_set<size_t> left_connected, right_connected;
            for (auto left_it = left_side.begin(); left_it != left_side.end() && !incompatible; ++left_it) {
                auto adj_it = connection_adj.find(*left_it);
                if (adj_it != connection_adj.end()) {
                    for (auto right_it = adj_it->second.begin(); right_it != adj_it->second.end() && !incompatible; ++right_it) {
                        
#ifdef debug_constrictions
                        cerr << "looking at connection between " << *left_it << " and " << *right_it << endl;
#endif
                        if (right_side.count(*right_it)) {
                            left_connected.insert(*left_it);
                            right_connected.insert(*right_it);
                        }
                        else {
                            // the direction of this connection are not consistent with the left and right
                            // side of this iteration
#ifdef debug_constrictions
                            cerr << "connection is incompatible" << endl;
#endif
                            incompatible = true;
                        }
                    }
                }
            }
            
            if (incompatible) {
                // the division of deletions to the left and right side is not compatible with the
                // connections
                return;
            }
            
            // do the non-connected edges form a biclique?
            for (auto left_it = left_side.begin(); left_it != left_side.end() && !incompatible; ++left_it) {
                if (left_connected.count(*left_it)) {
                    continue;
                }
                size_t num_clique_edges = 0;
                for (auto i : adj[*left_it]) {
                    if (right_connected.count(i)) {
                        // we don't worry about it if the node has a connection, because it will lose
                        // all of its edges anyway
                        continue;
                    }
                    if (right_side.count(i) && connected_by_edge(*left_it, i)) {
                        // this looks like it could be a splice junction
                        ++num_clique_edges;
                    }
                    else {
#ifdef debug_constrictions
                        cerr << "adjacency " << *left_it << " -> " << i << " is " << (right_side.count(i) ? "not connected by a graph edge" : "missing") << endl;
#endif
                        incompatible = true;
                    }
                }
                incompatible = incompatible || (num_clique_edges != right_side.size() - right_connected.size());
            }
            
            if (incompatible) {
                // we have edges going to outside the biclique, or we have edges missing
                // from the biclique
#ifdef debug_constrictions
                cerr << "this left-right combination is incompatible" << endl;
#endif
                return;
            }
            
            // count up the walks through this biclique
            size_t walk_total = 0;
            for (auto left_it = left_side.begin(); left_it != left_side.end() && !incompatible; ++left_it) {
                for (auto j : adj[*left_it]) {
                    walk_total += fwd[*left_it] * bwd[j];
                }
            }
            
#ifdef debug_constrictions
            cerr << "biclique has a walk total of " << walk_total << " compared to component total " << total_comp_paths[comp_idx] << endl;
#endif
            
            if (walk_total == total_comp_paths[comp_idx]) {
                // all of the walks in this component go through this biclique, we've found a splice
                
#ifdef debug_constrictions
                cerr << "recording a constriction biclique" << endl;
#endif
                
                return_val.emplace_back(vector<size_t>(left_side.begin(), left_side.end()),
                                        vector<size_t>(right_side.begin(), right_side.end()));
                
                auto& biclique = return_val.back();
                sort(biclique.first.begin(), biclique.first.end());
                sort(biclique.second.begin(), biclique.second.end());
            }
        };
        
        bool found_deletion_bicliques = false;
        
        for (const auto& end_record : chunks_by_end) {
            
            if (!chunks_by_begin.count(end_record.first)) {
                // there are no chunks starting at the same read position
                continue;
            }
            
            
#ifdef debug_constrictions
            cerr << "looking for a constriction at position " <<  end_record.first.second << " among lefts:" << endl;
            for (auto i : end_record.second) {
                cerr << "\t" << i << endl;
            }
            cerr << "and rights:" << endl;
            for (auto i : chunks_by_begin[end_record.first]) {
                cerr << "\t" << i << endl;
            }
#endif
            
            // collect all of the path chunks that have no aligned read sequence
            vector<size_t> deletion_chunks;
            for (auto i : end_record.second) {
                if (path_chunks[i].first.first == path_chunks[i].first.second) {
                    deletion_chunks.push_back(i);
                }
            }
            
            // deletion chunks can go on either end, so we try all combinations
            // TODO: magic number to prevent explosion
            for (size_t iter = 0, end = (1 << min<size_t>(deletion_chunks.size(), 16)); iter < end; ++iter) {
                
                
#ifdef debug_constrictions
                cerr << "deletion combination iteration " << iter << " out of " << end << endl;
#endif
                
                // we will fill out the left and right side of this potential splice biclique
                unordered_set<size_t> left_side, right_side;
                
                // fill left side, handling pure deletions according to the iteration
                size_t deletion_chunk_idx = 0;
                for (auto i : end_record.second) {
                    if (deletion_chunk_idx < deletion_chunks.size() && i == deletion_chunks[deletion_chunk_idx]) {
                        if (iter & (1 << deletion_chunk_idx)) {
                            left_side.insert(i);
                        }
                        else {
                            right_side.insert(i);
                        }
                    }
                    else {
                        left_side.insert(i);
                    }
                }
                // fill right side
                for (auto i : chunks_by_begin[end_record.first]) {
                    if (!left_side.count(i)) {
                        right_side.insert(i);
                    }
                }
                
                if (!deletion_chunks.empty()) {
                    // TODO: this solution will catch any bicliques that have a chance to be a constriction
                    // (which requires that there is only one adjacent to either side) as long as there is only
                    // one total to either side. in theory, it is possible to have further constrictions between
                    // layers of deletion chunks
                    
                    // collect the left neighbors of left deletions and the right neighbors
                    // of right deletions, as well as which of these are also in the left
                    // or right partition
                    unordered_set<size_t> left_neighbors, right_neighbors, lefter_side, righter_side;
                    for (auto i : deletion_chunks) {
                        if (left_side.count(i)) {
                            for (auto j : rev_adj[i]) {
                                left_neighbors.insert(j);
                                if (left_side.count(j)) {
                                    lefter_side.insert(j);
                                }
                            }
                        }
                        else {
                            for (auto j : adj[i]) {
                                right_neighbors.insert(j);
                                if (right_side.count(j)) {
                                    righter_side.insert(j);
                                }
                            }
                        }
                    }
                    
                    if (!lefter_side.empty()) {
                        // there's a group of nodes to the left of the left-side deletions that
                        // are still within the partition
                        
                        // check if these nodes connect *ONLY* through left side deletions
                        bool all_fully_to_left = true;
                        for (auto i : lefter_side) {
                            bool fully_to_left = true;
                            for (auto j : adj[i]) {
                                if (!left_side.count(j) || path_chunks[j].first.first != path_chunks[j].first.second) {
                                    fully_to_left = false;
                                    break;
                                }
                            }
                            if (fully_to_left) {
                                // this node is cut off from the rest of the graph by the left-side
                                // deletions, we shouldn't look for it in the main biclique
                                left_side.erase(i);
                            }
                            all_fully_to_left = all_fully_to_left && fully_to_left;
                        }
                        
                        if (all_fully_to_left && lefter_side.size() == left_neighbors.size()) {
                            // it looks like the left neighbors could even be their own biclique
                            
                            found_deletion_bicliques = true;
                            
#ifdef debug_constrictions
                            cerr << "left-er of left deletions" << endl;
                            for (auto i : lefter_side) {
                                cerr << "\t" << i << endl;
                            }
                            cerr << "iter lefts:" << endl;
                            for (auto i : left_side) {
                                cerr << "\t" << i << endl;
                            }
#endif
                            // and test the potential biclique to the left
                            test_candidate_biclique(lefter_side, left_side, end_record.first.first);
                        }
                    }
                    
                    // TODO: repetitive (but not structured well to re-use code...)
                    if (!righter_side.empty()) {
                        // there's a group of nodes to the right of the right-side deletions that
                        // are still within the partition
                        
                        // check if these nodes connect *ONLY* through right side deletions
                        bool all_fully_to_right = true;
                        for (auto i : righter_side) {
                            bool fully_to_right = true;
                            for (auto j : rev_adj[i]) {
                                if (!right_side.count(j) || path_chunks[j].first.first != path_chunks[j].first.second) {
                                    fully_to_right = false;
                                    break;
                                }
                            }
                            if (fully_to_right) {
                                // this node is cut off from the rest of the graph by the right-side
                                // deletions, we shouldn't look for it in the main biclique
                                right_side.erase(i);
                            }
                            all_fully_to_right = all_fully_to_right && fully_to_right;
                        }
                        
                        if (all_fully_to_right && righter_side.size() == right_neighbors.size()) {
                            // it looks like the right neighbors could even be their own biclique
                            
                            found_deletion_bicliques = true;
                            
#ifdef debug_constrictions
                            cerr << "iter rights:" << endl;
                            for (auto i : right_side) {
                                cerr << "\t" << i << endl;
                            }
                            cerr << "right-er of right deletions" << endl;
                            for (auto i : righter_side) {
                                cerr << "\t" << i << endl;
                            }
#endif
                            // and test the potential biclique to the right
                            test_candidate_biclique(right_side, righter_side, end_record.first.first);
                        }
                    }
                }
                
#ifdef debug_constrictions
                cerr << "iter left:" << endl;
                for (auto i : left_side) {
                    cerr << "\t" << i << endl;
                }
                cerr << "iter rights:" << endl;
                for (auto i : right_side) {
                    cerr << "\t" << i << endl;
                }
#endif
                test_candidate_biclique(left_side, right_side, end_record.first.first);
            }
        }
        
        if (found_deletion_bicliques) {
            // it's possible with all the deletion iterations to re-discover bicliques
            // multiple times. if we did, get rid of them
            
            sort(return_val.begin(), return_val.end());
            auto new_end = unique(return_val.begin(), return_val.end());
            return_val.resize(new_end - return_val.begin());
        }
        
        return return_val;
    }

    void Surjector::cut_anchors(bool rev_strand, vector<path_chunk_t>& path_chunks,
                                vector<pair<step_handle_t, step_handle_t>>& ref_chunks,
                                vector<tuple<size_t, size_t, int32_t>>& connections) const {
        
        // TODO: this is very repetitive with the similar function in the main spliced surject
        
        // distance along the path from end of chunk 1 to some mapping on chunk 2
        auto path_distance = [&](size_t chunk_idx_1,
                                 size_t chunk_idx_2, size_t mapping_idx_2) {
            
            step_handle_t step_1 = ref_chunks[chunk_idx_1].second;
            // move right from the beginning of the second chunk if necessary
            step_handle_t step_2 = ref_chunks[chunk_idx_2].first;
            for (size_t i = 0; i < chunk_idx_2; ++i) {
                step_2 = rev_strand ? graph->get_previous_step(step_2) : graph->get_next_step(step_2);
            }
            
            // get the distance component that is on the two mappings
            const auto& mapping_1 = *path_chunks[chunk_idx_1].second.mapping().rbegin();
            const auto& mapping_2 = path_chunks[chunk_idx_2].second.mapping(mapping_idx_2);
            int64_t dist = mapping_2.position().offset() - mapping_1.position().offset() - mapping_from_length(mapping_1);
            
            // get the distance component that is along the path
            if (rev_strand) {
                dist += (graph->get_position_of_step(step_1)
                         + graph->get_length(graph->get_handle_of_step(step_1))
                         - graph->get_position_of_step(step_2)
                         - graph->get_length(graph->get_handle_of_step(step_2)));
            }
            else {
                dist += (graph->get_position_of_step(step_2)
                         - graph->get_position_of_step(step_1));
            }
            return dist;
        };
        
        // put the input in lexicographic order by read interval
        vector<size_t> order = range_vector(path_chunks.size());
        stable_sort(order.begin(), order.end(), [&](size_t i, size_t j) {
            const auto& interval_1 = path_chunks[i].first;
            const auto& interval_2 = path_chunks[j].first;
            return (interval_1.first < interval_2.first ||
                    (interval_1.first == interval_2.first && interval_1.second < interval_2.second));
        });
        vector<size_t> index(order.size());
        for (size_t i = 0; i < order.size(); i++) {
            index[order[i]] = i;
        }
        // update the connection indexes
        for (auto& connection : connections) {
            get<0>(connection) = index[get<0>(connection)];
            get<1>(connection) = index[get<1>(connection)];
        }

        // swap the order
        for (size_t i = 0; i < path_chunks.size(); ++i) {
            while (index[i] != i) {
#ifdef debug_spliced_surject
                cerr << "reordering chunks, swapping " << i << " and " << index[i] << endl;
#endif
                std::swap(path_chunks[index[i]], path_chunks[i]);
                std::swap(ref_chunks[index[i]], ref_chunks[i]);
                std::swap(index[index[i]], index[i]);
            }
        }
        
        // find any overlaps we want to break
        vector<pair<size_t, size_t>> overlaps;
        for (size_t i = 0; i < path_chunks.size(); ++i) {
            auto i_end = path_chunks[i].first.second;
            for (size_t j = i + 1; j < path_chunks.size() && path_chunks[j].first.first < i_end; ++j) {
                if (i_end < path_chunks[j].first.second) {
                    // the first chunk overlaps the second (note: this check depends on sorting order)
                    
                    // figure out how far we have to go down the second chunk to get past the overlap
                    int64_t to_walk = i_end - path_chunks[j].first.first;
                    int64_t walked = 0;
                    size_t k = 0;
                    for (; k < path_chunks[j].second.mapping_size() && walked < to_walk; ++k) {
                        walked += mapping_to_length(path_chunks[j].second.mapping(k));
                    }
                    if (k < path_chunks[j].second.mapping_size() && path_distance(i, j, k) >= 0) {
                        // we didn't walk off the end of the of second chunk and they're colinear along the path
                        // so we record an overlap
                        overlaps.emplace_back(j, k);
#ifdef debug_spliced_surject
                        cerr << "path chunk " << i << " overlaps " << j << " by " << to_walk << " on read, marking an overlap to split before mapping " << k << " at " << pb2json(path_chunks[j].second.mapping(k)) << endl;
#endif
                    }
                }
            }
        }
        
        if (!overlaps.empty()) {
            
            sort(overlaps.begin(), overlaps.end());
            overlaps.resize(unique(overlaps.begin(), overlaps.end()) - overlaps.begin());
            
#ifdef debug_spliced_surject
            cerr << "performing overlap splits: " << endl;
            for (auto overlap : overlaps) {
                cerr << "\t" << overlap.first << ", " << overlap.second << endl;
            }
#endif
            
            vector<path_chunk_t> split_path_chunks;
            split_path_chunks.reserve(path_chunks.size() + overlaps.size());
            vector<pair<step_handle_t, step_handle_t>> split_ref_chunks;
            split_ref_chunks.reserve(ref_chunks.size() + overlaps.size());
            
            vector<size_t> added_before(path_chunks.size(), 0);
            for (size_t i = 0, j = 0; i < path_chunks.size(); ++i) {
                if (i > 0) {
                    added_before[i] = added_before[i - 1];
                }
                
                if (j < overlaps.size() && overlaps[j].first == i) {
                    // find out how many overlap splits we need to perform
                    size_t n = 1;
                    while (j + n < overlaps.size() && overlaps[j + n].first == i) {
                        ++n;
                    }
#ifdef debug_spliced_surject
                    cerr << "path chunk " << i << " has " << n << " overlap splits" << endl;
#endif
                    // we'll walk along the ref path and the read intervals as we go
                    step_handle_t step = ref_chunks[i].first;
                    auto read_begin = path_chunks[i].first.first;
                    
                    for (size_t k = 0; k <= n; ++k) {
                        // figure out the bounds of mappings we'll move over
                        size_t begin_idx = (k == 0 ? 0 : overlaps[j + k - 1].second);
                        size_t end_idx = (k == n ? path_chunks[i].second.mapping_size() : overlaps[j + k].second);
                        
                        // add the mappings
                        split_path_chunks.emplace_back();
                        auto& path_chunk = split_path_chunks.back();
                        for (size_t l = begin_idx; l < end_idx; ++l) {
                            auto mapping = path_chunk.second.add_mapping();
                            *mapping = path_chunks[i].second.mapping(l);
                            mapping->set_rank(l - begin_idx + 1);
                        }
                        // identify the read interval
                        path_chunk.first.first = read_begin;
                        path_chunk.first.second = path_chunk.first.first + path_to_length(path_chunk.second);
                        read_begin = path_chunk.first.second;
                        
                        // walk the reference path steps
                        split_ref_chunks.emplace_back();
                        auto& ref_chunk = split_ref_chunks.back();
                        ref_chunk.first = step;
                        for (size_t l = begin_idx + 1; l < end_idx; ++l) {
                            step = rev_strand ? graph->get_previous_step(step) : graph->get_next_step(step);
                        }
                        ref_chunk.second = step;
                        // set up the step for the next iteration
                        step = rev_strand ? graph->get_previous_step(step) : graph->get_next_step(step);
#ifdef debug_spliced_surject
                        cerr << "next split for chunk " << i << " as " << split_path_chunks.size() - 1 << ", consisting of " << endl;
                        cerr << "\t" << string(path_chunk.first.first, path_chunk.first.second) << endl;
                        cerr << "\t" << pb2json(path_chunk.second) << endl;
                        cerr << "\t" << graph->get_position_of_step(ref_chunk.first) << " : " << graph->get_position_of_step(ref_chunk.second) << endl;
#endif
                    }
                    j += n;
                    added_before[i] += n;
                }
                else {
#ifdef debug_spliced_surject
                    cerr << "no splits on chunk " << i << ", add as " << split_path_chunks.size() << endl;
                    cerr << "\t" << string(path_chunks[i].first.first, path_chunks[i].first.second) << endl;
                    cerr << "\t" << pb2json(path_chunks[i].second) << endl;
                    cerr << "\t" << graph->get_position_of_step(ref_chunks[i].first) << " : " << graph->get_position_of_step(ref_chunks[i].second) << endl;
#endif
                    split_path_chunks.emplace_back(move(path_chunks[i]));
                    split_ref_chunks.emplace_back(move(ref_chunks[i]));
                    
                }
            }
            
            // replace the original path chunks and ref chunks with the split ones
            path_chunks = move(split_path_chunks);
            ref_chunks = move(split_ref_chunks);
            
            // and update the indexes of the connections
            for (auto& connection : connections) {
                // edges out should be updated for the splits added in that iteration because
                // the come out of the last split segment
                get<0>(connection) += added_before[get<0>(connection)];
                // edges in should only be updated for the splits that happened in earlier
                // iterations (also, these should never be in index 0, would violate
                // colinearity)
                get<1>(connection) += added_before[get<1>(connection) - 1];
            }
        }
    }

    multipath_alignment_t Surjector::spliced_surject(const PathPositionHandleGraph* path_position_graph,
                                                     const string& src_sequence, const string& src_quality,
                                                     const int32_t src_mapping_quality,
                                                     const path_handle_t& path_handle, bool rev_strand,
                                                     vector<path_chunk_t>& path_chunks,
                                                     vector<pair<step_handle_t, step_handle_t>>& ref_chunks,
                                                     vector<tuple<size_t, size_t, int32_t>>& connections,
                                                     pair<step_handle_t, step_handle_t>& path_range_out,
                                                     bool allow_negative_scores, bool deletions_as_splices) const {
                
        assert(path_chunks.size() == ref_chunks.size());
        
        function<int64_t(size_t,size_t)> path_distance = [&](size_t i, size_t j) {
            const auto& final_mapping = *path_chunks[i].second.mapping().rbegin();
            int64_t dist = (path_chunks[j].second.mapping(0).position().offset()
                            - final_mapping.position().offset()
                            - mapping_from_length(final_mapping));
            if (rev_strand) {
                dist += (graph->get_position_of_step(ref_chunks[i].second)
                         + graph->get_length(graph->get_handle_of_step(ref_chunks[i].second))
                         - graph->get_position_of_step(ref_chunks[j].first)
                         - graph->get_length(graph->get_handle_of_step(ref_chunks[j].first)));
            }
            else {
                dist += (graph->get_position_of_step(ref_chunks[j].first)
                         - graph->get_position_of_step(ref_chunks[i].second));
            }
            return dist;
        };
        
        multipath_alignment_t surjected;
        
        // checks whether the end of i is connected to the beginning of j by an edge
        
        
#ifdef debug_spliced_surject
        cerr << "removing any pure insertion path chunks" << endl;
#endif
        
        vector<size_t> insertions_removed(path_chunks.size() + 1, 0);
        for (size_t i = 0; i < path_chunks.size(); ++i) {
            const auto& chunk = path_chunks[i].second;
            bool has_aligned_bases = false;
            for (size_t j = 0; j < chunk.mapping_size() && !has_aligned_bases; ++j) {
                const auto& mapping = chunk.mapping(j);
                for (size_t k = 0; k < mapping.edit_size() && !has_aligned_bases; ++k) {
                    has_aligned_bases = mapping.edit(k).from_length() != 0;
                }
            }
            if (!has_aligned_bases) {
                insertions_removed[i + 1] = insertions_removed[i] + 1;
            }
            else {
                insertions_removed[i + 1] = insertions_removed[i];
                if (insertions_removed[i]) {
                    path_chunks[i - insertions_removed[i]] = move(path_chunks[i]);
                    ref_chunks[i - insertions_removed[i]] = move(ref_chunks[i]);
                }
            }
        }
        
        if (insertions_removed.back()) {
            path_chunks.resize(path_chunks.size() - insertions_removed.back());
            ref_chunks.resize(path_chunks.size());
            
            // update connections with the new indexes
            size_t connections_removed = 0;
            for (size_t i = 0; i < connections.size(); ++i) {
                auto& connection = connections[i];
                if (insertions_removed[get<0>(connection)] != insertions_removed[get<0>(connection) + 1]
                    || insertions_removed[get<1>(connection)] != insertions_removed[get<1>(connection) + 1]) {
                    ++connections_removed;
                }
                else {
                    get<0>(connection) -= insertions_removed[get<0>(connection)];
                    get<1>(connection) -= insertions_removed[get<1>(connection)];
                    connections[i - connections_removed] = connection;
                }
            }
            
            connections.resize(connections.size() - connections_removed);
        }
        
#ifdef debug_spliced_surject
        cerr << "removed " << insertions_removed.back() << " chunks" << endl;
        cerr << "checking for need to cut anchors" << endl;
#endif
        
        cut_anchors(rev_strand, path_chunks, ref_chunks, connections);
        
#ifdef debug_spliced_surject
        cerr << "making colinearity graph for " << path_chunks.size() << " path chunks" << endl;
#endif
         vector<vector<size_t>> colinear_adj(path_chunks.size());
        
        for (size_t i = 0; i < path_chunks.size(); ++i) {
            for (size_t j = i + 1; j < path_chunks.size(); ++j) {
                if (path_chunks[i].first.second <= path_chunks[j].first.first
                    && path_distance(i, j) >= 0) {
                    // the second one is further along both the read and the path, so it is colinear
                    colinear_adj[i].push_back(j);
                }
            }
        }
        
        // TODO: use tricks from multipath alignment graph to make a smaller initial chunk graph
        
#ifdef debug_spliced_surject
        cerr << "initial graph:" << endl;
        for (size_t i = 0; i < colinear_adj.size(); ++i) {
            cerr << i << ":";
            for (auto j : colinear_adj[i]) {
                cerr << " " << j;
            }
            cerr << endl;
        }
        cerr << "connections:" << endl;
        for (const auto& connection : connections) {
            cerr << get<0>(connection) << " -> " << get<1>(connection) << ", " << get<2>(connection) << endl;
        }
        
        cerr << "computing transitive reduction" << endl;
#endif
        
        // remove transitive edges
        vector<vector<size_t>> colinear_adj_red = transitive_reduction(colinear_adj);
        
#ifdef debug_spliced_surject
        cerr << "reduced graph:" << endl;
        for (size_t i = 0; i < colinear_adj_red.size(); ++i) {
            cerr << i << ":";
            for (auto j : colinear_adj_red[i]) {
                cerr << " " << j;
            }
            cerr << endl;
        }
        
        cerr << "removing dominated path chunks" << endl;
#endif
        
        colinear_adj_red = remove_dominated_chunks(src_sequence, colinear_adj_red, path_chunks, ref_chunks, connections);
        
#ifdef debug_spliced_surject
        cerr << "with dominated chunks removed:" << endl;
        for (size_t i = 0; i < colinear_adj_red.size(); ++i) {
            cerr << i << ":";
            for (auto j : colinear_adj_red[i]) {
                cerr << " " << j;
            }
            cerr << endl;
        }
        cerr << "connections:" << endl;
        for (const auto& connection : connections) {
            cerr << get<0>(connection) << " -> " << get<1>(connection) << ", " << get<2>(connection) << endl;
        }
#endif
        
        // records of (to idx, score, is a connection)
        vector<vector<tuple<size_t, int32_t, bool>>> splice_edges(path_chunks.size());
        
        vector<bool> has_inward_connection(path_chunks.size(), false);
        
        if (!connections.empty()) {
            
            // clear outward edges for chunks that send connections, and record
            // the scored edge
            
#ifdef debug_spliced_surject
            cerr << "handling any connections" << endl;
#endif
            
            unordered_set<pair<size_t, size_t>> connection_set;
            for (const auto& connection : connections) {
                connection_set.emplace(get<0>(connection), get<1>(connection));
            }
            
            for (const auto& connection : connections) {
                splice_edges[get<0>(connection)].emplace_back(get<1>(connection), get<2>(connection), true);
                has_inward_connection[get<1>(connection)] = true;
                // move direct adjacency edges out of this node to the splice edges (unless they correspond to the
                // the connection itself).
                for (auto target : colinear_adj_red[get<0>(connection)]) {
                    if (!connection_set.count(make_pair(get<0>(connection), target))
                        && path_chunks[get<0>(connection)].first.second == path_chunks[target].first.first
                        && path_distance(get<0>(connection), target) == 0) {
                        // TODO: why do we find directly abutting connections in the first place?
                        splice_edges[get<0>(connection)].emplace_back(target, 0, false);
                    }
                }
                colinear_adj_red[get<0>(connection)].clear();
            }
            
            
            // move inward exactly abutting edges for path chunks that receive connections into the splice edges
            // and clear the rest
            for (auto& adj : colinear_adj_red) {
                for (size_t i = 0; i < adj.size();) {
                    if (has_inward_connection[adj[i]]) {
                        if (!connection_set.count(make_pair(i, adj[i]))
                            && path_chunks[i].first.second == path_chunks[adj[i]].first.first
                            && path_distance(i, adj[i]) == 0) {
                            splice_edges[i].emplace_back(adj[i], 0, false);
                        }
                        adj[i] = adj.back();
                        adj.pop_back();
                    }
                    else {
                        ++i;
                    }
                }
            }
            
#ifdef debug_spliced_surject
            cerr << "after removing connections:" << endl;
            for (size_t i = 0; i < colinear_adj_red.size(); ++i) {
                cerr << i << ":";
                for (auto j : colinear_adj_red[i]) {
                    cerr << " " << j;
                }
                cerr << endl;
            }
            cerr << "splice graph:" << endl;
            for (size_t i = 0; i < splice_edges.size(); ++i) {
                cerr << i << ":";
                for (auto edge : splice_edges[i]) {
                    cerr << " (" << get<0>(edge) << ", " << get<1>(edge) << ", " << get<2>(edge) << ")";
                }
                cerr << endl;
            }
#endif
        }
        
        
        if (deletions_as_splices) {
            
            // look for constrictions and move them into the splice edges iteratively
            // (some edges that are not originally constrictions can become constrictions
            // once other constriction edges are removed, which separates the component)
            
            bool removed_edges = true;
            while (removed_edges) {
                removed_edges = false;
                
#ifdef debug_spliced_surject
                cerr << "finding constrictions" << endl;
#endif
                
                // find bicliques that constrict the colinearity graph
                auto constrictions = find_constriction_bicliques(colinear_adj_red, src_sequence,
                                                                 path_chunks, connections);
                
#ifdef debug_spliced_surject
                cerr << "found " << constrictions.size() << " constriction bicliques:" << endl;
                for (auto& constriction : constrictions) {
                    cerr << "left:" << endl;
                    for (auto i : constriction.first) {
                        cerr << "\t" << i << endl;
                    }
                    cerr << "right:" << endl;
                    for (auto i : constriction.second) {
                        cerr << "\t" << i << endl;
                    }
                }
#endif
                
                // if any constrictions correspond to pure deletions, remove them from the colineary
                // graph and record them as splice edges
                
                for (const auto& constriction : constrictions) {
                    
                    vector<tuple<size_t, size_t, int64_t>> new_edges;
                    bool includes_splice = false;
                    for (auto i : constriction.first) {
                        if (colinear_adj_red[i].empty()) {
                            // the edges have been cleared when incorporating a connection
                            continue;
                        }
                        for (auto j : constriction.second) {
                            if (has_inward_connection[j]) {
                                // backward edgs have been removed
                                continue;
                            }
                            int64_t dist = path_distance(i, j);
                            int64_t score;
                            if (dist >= min_splice_length) {
                                includes_splice = true;
                                score = 0;
                            }
                            else {
                                score = get_aligner(!src_quality.empty())->score_gap(dist);
                            }
                            
#ifdef debug_spliced_surject
                            cerr << "deletion of length " << dist << " from " << i << " to " << j << " is recorded as part of a splice biclique, and given score " << score << endl;
#endif
                            
                            new_edges.emplace_back(i, j, score);
                        }
                    }
                    if (includes_splice) {
                        removed_edges = true;
                        // remove the colinearity edges
                        for (auto i : constriction.first) {
                            colinear_adj_red[i].clear();
                        }
                        // transfer them to splice edges
                        for (const auto& edge : new_edges) {
                            splice_edges[get<0>(edge)].emplace_back(get<1>(edge), get<2>(edge), false);
                        }
                    }
                }
                
                
#ifdef debug_spliced_surject
                cerr << "after removing long constriction deletions:" << endl;
                for (size_t i = 0; i < colinear_adj_red.size(); ++i) {
                    cerr << i << ":";
                    for (auto j : colinear_adj_red[i]) {
                        cerr << " " << j;
                    }
                    cerr << endl;
                }
                cerr << "splice graph:" << endl;
                for (size_t i = 0; i < splice_edges.size(); ++i) {
                    cerr << i << ":";
                    for (auto edge : splice_edges[i]) {
                        cerr << " (" << get<0>(edge) << ", " << get<1>(edge) << ", " << get<2>(edge) << ")";
                    }
                    cerr << endl;
                }
#endif
                
            }
        }
        
#ifdef debug_spliced_surject
        cerr << "computing constriction components" << endl;
#endif
      
        // find the connected components in the graph with the splice edges removed
        size_t num_comps = 0;
        vector<size_t> constriction_comps = connected_components(colinear_adj_red,
                                                                 reverse_adjacencies(colinear_adj_red),
                                                                 &num_comps);
        vector<vector<size_t>> comp_groups(num_comps);
        for (size_t i = 0; i < constriction_comps.size(); ++i) {
            comp_groups[constriction_comps[i]].push_back(i);
        }
        
#ifdef debug_spliced_surject
        for (size_t i = 0; i < comp_groups.size(); ++i) {
            cerr << "group " << i << ":";
            for (auto j : comp_groups[i]) {
                cerr << " " << j;
            }
            cerr << endl;
        }
#endif
        
        prune_unconnectable(colinear_adj_red, splice_edges, constriction_comps, comp_groups,
                            path_chunks, ref_chunks);
        
#ifdef debug_spliced_surject
        cerr << "groups after pruning unconnectable" << endl;
        for (size_t i = 0; i < comp_groups.size(); ++i) {
            cerr << "group " << i << ":";
            for (auto j : comp_groups[i]) {
                cerr << " " << j;
            }
            cerr << endl;
        }
#endif
        
        
        // convert the splice edges into edges between the components and identify sources/sinks
        vector<bool> comp_is_source(comp_groups.size(), true);
        vector<vector<tuple<size_t, int32_t, bool>>> comp_group_edges(comp_groups.size());
        for (size_t i = 0; i < splice_edges.size(); ++i) {
            for (auto& edge : splice_edges[i]) {
                if (constriction_comps[i] < constriction_comps[get<0>(edge)]) {
                    comp_group_edges[constriction_comps[i]].emplace_back(constriction_comps[get<0>(edge)],
                                                                         get<1>(edge), get<2>(edge));
                    comp_is_source[constriction_comps[get<0>(edge)]] = false;
                }
            }
        }
        
        
#ifdef debug_spliced_surject
        cerr << "component group edges:" << endl;
        for (size_t i = 0; i < comp_group_edges.size(); ++i) {
            cerr << i << ":";
            for (auto edge : comp_group_edges[i]) {
                cerr << " (" << get<0>(edge) << ", " << get<1>(edge) << ", " << get<2>(edge) << ")";
            }
            cerr << endl;
        }
#endif
        
#ifdef debug_spliced_surject
        cerr << "surjecting " << comp_groups.size() << " constriction sections" << endl;
#endif
        
        vector<Alignment> sections;
        vector<pair<step_handle_t, step_handle_t>> section_path_ranges;
        
        for (size_t i = 0; i < comp_groups.size(); ++i) {
            pair<string::const_iterator, string::const_iterator> read_range;
            vector<path_chunk_t> section_path_chunks;
                                    
            vector<size_t>& group = comp_groups[i];
                        
            // the other end points are determine by how the portion of the
            for (size_t j = 0; j < group.size(); ++j) {
                
                section_path_chunks.push_back(path_chunks[group[j]]);
                
                if (j == 0 || read_range.first > path_chunks[group[j]].first.first) {
                    read_range.first = path_chunks[group[j]].first.first;
                }
                if (j == 0 || read_range.second < path_chunks[group[j]].first.second) {
                    read_range.second = path_chunks[group[j]].first.second;
                }
            }
            
            // sources/sinks align all the way to the end
            if (comp_is_source[i]) {
                read_range.first = src_sequence.begin();
            }
            if (comp_group_edges[i].empty()) {
                read_range.second = src_sequence.end();
            }
            
            // make a dummy alignment with the relevant portion of the sequence
            Alignment section_source;
            *section_source.mutable_sequence() = string(read_range.first, read_range.second);
            if (!src_quality.empty()) {
                *section_source.mutable_quality() = string(src_quality.begin() + (read_range.first - src_sequence.begin()),
                                                           src_quality.begin() + (read_range.second - src_sequence.begin()));
            }
#ifdef debug_always_warn_on_too_long
            // give it the full sequence as a name so we can see it later
            section_source.set_name(src_sequence);
#endif
            
            // update the path chunk ranges to point into the dummy section read
            for (size_t i = 0; i < section_path_chunks.size(); ++i) {
                auto& chunk_range = section_path_chunks[i].first;
                chunk_range.first = section_source.sequence().begin() + (chunk_range.first - read_range.first);
                chunk_range.second = section_source.sequence().begin() + (chunk_range.second - read_range.first);
            }
            
#ifdef debug_spliced_surject
            cerr << "surjecting section " << i << ": " << pb2json(section_source) << endl;
            cerr << "consists of " << section_path_chunks.size() << " path chunks" << endl;
#endif
            
            // perform a full length surjection within the section section
            section_path_ranges.emplace_back();
            sections.push_back(realigning_surject(graph, section_source, path_handle, rev_strand, section_path_chunks,
                                                  section_path_ranges.back(), true, true, true));
            
            // remove any extraneous full length bonuses
            // TODO: technically, this can give a non-optimal alignment because it's post hoc to the dynamic programming
            const auto& aligner = *get_aligner(!src_quality.empty());
            if (sections.back().path().mapping_size() != 0) {
                if (read_range.first != src_sequence.begin()) {
                    if (sections.back().path().mapping(0).edit(0).from_length() > 0) {
                        sections.back().set_score(sections.back().score()
                                                  - aligner.score_full_length_bonus(true, section_source));
                    }
                }
                if (read_range.second != src_sequence.end()) {
                    const Mapping& m = sections.back().path().mapping(0);
                    if (m.edit(m.edit_size() - 1).from_length() > 0) {
                        sections.back().set_score(sections.back().score()
                                                  - aligner.score_full_length_bonus(false, section_source));
                    }
                }
            }
            
#ifdef debug_spliced_surject
            cerr << "surjected section " << i << " after score adjustment: " << pb2json(sections.back()) << endl;
            if (sections.back().path().mapping_size() == 0) {
                cerr << "null alignment has no path range" << endl;
            }
            else {
                cerr << "path range: " << graph->get_id(graph->get_handle_of_step(section_path_ranges.back().first)) << " " << graph->get_is_reverse(graph->get_handle_of_step(section_path_ranges.back().first)) << " " << graph->get_position_of_step(section_path_ranges.back().first) << " : " << graph->get_id(graph->get_handle_of_step(section_path_ranges.back().second)) << " " << graph->get_is_reverse(graph->get_handle_of_step(section_path_ranges.back().second)) << " " << graph->get_position_of_step(section_path_ranges.back().second) << endl;
            }
#endif
        }
        
        // distance between the path ranges of two sections
        // assumes direct adjacency over an edge, but this may not be true in the case of a connection
        // TODO: repetitive with path_dist
        auto section_path_dist = [&](size_t i, size_t j) {
            if (rev_strand) {
                return (graph->get_position_of_step(section_path_ranges[i].second)
                        - graph->get_position_of_step(section_path_ranges[j].first)
                        - graph->get_length(graph->get_handle_of_step(section_path_ranges[j].first)));
                
            }
            else {
                return (graph->get_position_of_step(section_path_ranges[j].first)
                        - graph->get_position_of_step(section_path_ranges[i].second)
                        - graph->get_length(graph->get_handle_of_step(section_path_ranges[i].second)));
            }
        };
        
#ifdef debug_spliced_surject
        cerr << "computing optimal combination of sections" << endl;
#endif
        
        // now we find use dynamic programming to find the best alignment across chunks
        
        vector<int64_t> backpointer(sections.size(), -1);
        vector<int32_t> score_dp(sections.size(), numeric_limits<int32_t>::min());
        
        // initialize the scores at sources or at any section if we're doing subpath
        // local alignments (i.e. not allowing negative scores)
        for (size_t i = 0; i < sections.size(); ++i) {
            if (!allow_negative_scores || comp_is_source[i]) {
                score_dp[i] = sections[i].score();
            }
        }
        
        // do the dynamic programming
        for (size_t i = 0; i < comp_groups.size(); ++i) {
            
            for (auto& edge : comp_group_edges[i]) {
                
                int32_t extended_score = score_dp[i] + get<1>(edge) + sections[get<0>(edge)].score();
                
#ifdef debug_spliced_surject
                cerr << "extending from component " << i << " (DP score " << score_dp[i] << ") with score of " << extended_score << " to " << get<0>(edge) << " (DP score " << score_dp[get<0>(edge)] << ") dist " << section_path_dist(i, get<0>(edge)) << endl;
#endif
                
                if (extended_score > score_dp[get<0>(edge)]) {
                    score_dp[get<0>(edge)] = extended_score;
                    backpointer[get<0>(edge)] = i;
                }
                else if (extended_score == score_dp[get<0>(edge)]
                         && sections[i].path().mapping_size() != 0
                         && sections[get<0>(edge)].path().mapping_size() != 0
                         && backpointer[get<0>(edge)] >= 0
                         && section_path_dist(i, get<0>(edge)) < section_path_dist(backpointer[get<0>(edge)], get<0>(edge))) {
                    // break ties in favor of the closer exon
                    backpointer[get<0>(edge)] = i;
                }
            }
        }
        
        // find the maximum, subject to full length requirements
        vector<size_t> traceback(1, -1);
        int32_t max_score = numeric_limits<int32_t>::min();
        for (size_t i = 0; i < score_dp.size(); ++i) {
            
            if (score_dp[i] > max_score && (!allow_negative_scores || comp_group_edges[i].empty())) {
                max_score = score_dp[i];
                traceback[0] = i;
            }
            else if (score_dp[i] == max_score
                     && backpointer[i] != -1
                     && traceback[0] != -1
                     && backpointer[traceback[0]] != -1
                     && (!allow_negative_scores || comp_group_edges[i].empty())
                     && section_path_dist(backpointer[i], i) < section_path_dist(backpointer[traceback[0]], traceback[0])) {
                // break ties in favor exon with closer connection
                traceback[0] = i;
            }
        }
        
        // follow the back pointers
        while (backpointer[traceback.back()] != -1) {
            traceback.push_back(backpointer[traceback.back()]);
        }
        
        
#ifdef debug_spliced_surject
        cerr << "combining " << traceback.size() << " sections into surjected alignment" << endl;
        for (int64_t i = traceback.size() - 1; i >= 0; --i) {
            cerr << "\t" << traceback[i] << endl;
        }
#endif
        
        if (traceback.empty()) {
            // sentinel for unmapped
            path_range_out.first = path_range_out.second = graph->path_end(path_handle);
        }
        else {
            path_range_out.first = section_path_ranges[traceback.back()].first;
            path_range_out.second = section_path_ranges[traceback.front()].second;
        }
        
        // make an alignment to build out the path in
        surjected.set_sequence(src_sequence);
        surjected.set_quality(src_quality);
        surjected.set_mapping_quality(src_mapping_quality);
                
        subpath_t* prev_subpath = nullptr;
        for (int64_t i = traceback.size() - 1; i >= 0; --i) {
            
            size_t section_idx = traceback[i];
            const Path& copy_path = sections[section_idx].path();
            
            if (copy_path.mapping_size() == 0) {
                // the DP chose a segment that was unsurjectable
                path_range_out.first = path_range_out.second = graph->path_end(path_handle);
                surjected = make_null_mp_alignment(src_sequence, src_quality);
                return surjected;
            }
#ifdef debug_spliced_surject
            cerr << "appending path section " << pb2json(copy_path) << endl;
#endif
            if (copy_path.mapping_size() == 0) {
                // this happens if the surjected section is a pure deletion, we can just skip it
                continue;
            }
            
            if (i != traceback.size() - 1) {
                // make an edge back to the previous section
                size_t prev_idx = traceback[i + 1];
                // find the edge between these sections that the DP used
                for (auto& edge : comp_group_edges[prev_idx]) {
                    if (get<0>(edge) == section_idx &&
                        get<1>(edge) + score_dp[prev_idx] + sections[section_idx].score() == score_dp[section_idx]) {
                        if (get<2>(edge)) {
                            // this is from a connection
                            auto connection = prev_subpath->add_connection();
                            connection->set_next(surjected.subpath_size());
                            connection->set_score(get<1>(edge));
                        }
                        else {
                            // this is from a constriction or a preserved edge around a connection
                            prev_subpath->add_next(surjected.subpath_size());
                        }
                        break;
                    }
                }
            }
            // TODO: merge adjacent mappings across subpaths
            
            // make a new subpath to hold the this path section
            auto surj_subpath = surjected.add_subpath();
            surj_subpath->set_score(sections[section_idx].score());
            from_proto_path(copy_path, *surj_subpath->mutable_path());
            
            prev_subpath = surj_subpath;
        }
        
        // since the mp aln is a non-branching path, this is always the only start
        if (surjected.subpath_size() != 0) {
            surjected.add_start(0);
        }
        
#ifdef debug_spliced_surject
        cerr << "final spliced surjection " << debug_string(surjected) << endl;
#endif
        
        return surjected;
    }

    Alignment Surjector::realigning_surject(const PathPositionHandleGraph* path_position_graph, const Alignment& source,
                                            const path_handle_t& path_handle, bool rev_strand, const vector<path_chunk_t>& path_chunks,
                                            pair<step_handle_t, step_handle_t>& path_range_out, bool allow_negative_scores,
                                            bool preserve_N_alignments, bool preserve_tail_indel_anchors) const {
        
#ifdef debug_anchored_surject
        cerr << "using overlap chunks on path " << graph->get_path_name(path_handle) << " strand " << rev_strand << ", performing realigning surjection" << endl;
        cerr << "chunks:" << endl;
        for (size_t i = 0; i < path_chunks.size(); ++i) {
            cerr << "\t" << string(path_chunks[i].first.first, path_chunks[i].first.second) << ", " << pb2json(path_chunks[i].second) << endl;
        }
#endif
        
        // the alignment we will fill out
        Alignment surjected;
        
        // find the end-inclusive interval of the ref path we need to consider
        pair<size_t, size_t> ref_path_interval = compute_path_interval(path_position_graph, source, path_handle,
                                                                       rev_strand, path_chunks);
        if (ref_path_interval.first <= ref_path_interval.second) {
            // We actually got a nonempty range, so expand it.
            
            // having a buffer helps ensure that we get the correct anchoring position for some edge cases
            // of a full deletion that occurs on a node boundary
            if (ref_path_interval.first > 0) {
                --ref_path_interval.first;
            }
            if (ref_path_interval.second + 1 < path_position_graph->get_path_length(path_handle)) {
                ++ref_path_interval.second;
            }
        }
        
        if (path_chunks.size() == 0) {
#ifdef debug_anchored_surject
            cerr << "no path chunks provided, surjecting as unmapped" << endl;
#endif
            // Leave surjected path empty
        }
        else if (path_chunks.size() == 1
            && path_chunks.front().first.first == source.sequence().begin()
            && path_chunks.front().first.second == source.sequence().end()) {
#ifdef debug_anchored_surject
            cerr << "path chunk already constitutes a full alignment, skipping realignment" << endl;
#endif
            
            // just copy it over
            surjected.set_sequence(source.sequence());
            surjected.set_quality(source.quality());
            *surjected.mutable_path() = path_chunks.front().second;
            surjected.set_score(get_aligner(!source.quality().empty())->score_contiguous_alignment(surjected));
            
        }
        else {
            // we're going to have to realign some portions
            
#ifdef debug_anchored_surject
            cerr << "final path interval is " << ref_path_interval.first << ":" << ref_path_interval.second << " on path of length " << path_position_graph->get_path_length(path_handle) << endl;
#endif

            // If we put in path chunks we need to have ended up with a
            // nonempty path interval that they cover.
            assert(ref_path_interval.first <= ref_path_interval.second);
            
            // get the path graph corresponding to this interval
            bdsg::HashGraph path_graph;
            unordered_map<id_t, pair<id_t, bool>> path_trans = extract_linearized_path_graph(path_position_graph, &path_graph, path_handle,
                                                                                             ref_path_interval.first, ref_path_interval.second);
            
            // split it into a forward and reverse strand
            // TODO: we only should need one strand of the graph, but it might be different strands on
            // different nodes...
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

            size_t subgraph_bases = split_path_graph.get_total_length();
            if (subgraph_bases > max_subgraph_bases) {
#ifdef debug_always_warn_on_too_long
                cerr << "gave up on too long read " + source.name() + "\n";
#endif
                if (!warned_about_subgraph_size.test_and_set()) {
                    cerr << "warning[vg::Surjector]: Refusing to perform very large alignment against "
                        << subgraph_bases << " bp strand split subgraph for read " << source.name()
                        << "; suppressing further warnings." << endl;
                }
                return move(make_null_alignment(source)); 
            }
            
            // compute the connectivity between the path chunks
            MultipathAlignmentGraph mp_aln_graph(split_path_graph, path_chunks, source, node_trans, !preserve_N_alignments,
                                                 preserve_tail_indel_anchors);
            
            // we don't overlap this reference path at all or we filtered out all of the path chunks, so just make a sentinel
            if (mp_aln_graph.empty()) {
                return move(make_null_alignment(source));
            }
            
            // TODO: is this necessary in a linear graph?
            vector<size_t> topological_order;
            mp_aln_graph.topological_sort(topological_order);
            mp_aln_graph.remove_transitive_edges(topological_order);
            
            // align the intervening segments and store the result in a multipath alignment
            multipath_alignment_t mp_aln;
            mp_aln_graph.align(source, split_path_graph, get_aligner(),
                               false,                                    // anchors as matches
                               1,                                        // max alt alns
                               false,                                    // dynamic alt alns
                               numeric_limits<int64_t>::max(),           // max gap
                               0.0,                                      // pessimistic tail gap multiplier
                               false,                                    // simplify topologies
                               0,                                        // unmergeable len
                               1,                                        // band padding
                               mp_aln,                                   // output
                               nullptr,                                  // snarl manager
                               nullptr,                                  // distance index
                               nullptr,                                  // projector
                               allow_negative_scores);
            topologically_order_subpaths(mp_aln);
            
            if (preserve_tail_indel_anchors) {
                // this code path sometimes produces subpaths that have no aligned bases, which
                // sometimes play poorly with other parts of the code base
                remove_empty_alignment_sections(mp_aln);
            }
            
            for (size_t i = 0; i < mp_aln.subpath_size(); i++) {
                // translate back into the original ID space
                translate_oriented_node_ids(*mp_aln.mutable_subpath(i)->mutable_path(), node_trans);
            }
            
            // identify the source subpaths (necessary for subpath-global optimal alignment algorithm)
            identify_start_subpaths(mp_aln);
            
#ifdef debug_anchored_surject
            cerr << "made multipath alignment " << debug_string(mp_aln) << endl;
#endif
            
#ifdef debug_validate_anchored_multipath_alignment
            if (!validate_multipath_alignment(mp_aln, *graph)) {
                cerr << "WARNING: multipath alignment for surjection of " << source.name() << " failed to validate" << endl;
            }
#endif
            // concatenate the subpaths either locally or globally, depending on whether we're
            // allowing negative scores
            optimal_alignment(mp_aln, surjected, allow_negative_scores);
        }
        
        const auto& surj_path = surjected.path();
        if (surj_path.mapping_size() > 0) {
            // the surjection is mapped
            
#ifdef debug_anchored_surject
            cerr << "assigning a path range to surjected path: " << pb2json(surj_path) << endl;
#endif
            size_t mappings_matched = 0;
            
            // look in either the forward or reverse orientation along the path
            
#ifdef debug_anchored_surject
            cerr << "looking for path range on " << (rev_strand ? "reverse" : "forward") << " strand, for " << surj_path.mapping_size() << " mappings" << endl;
#endif
            step_handle_t step = rev_strand ? graph->get_step_at_position(path_handle, ref_path_interval.second)
                                            : graph->get_step_at_position(path_handle, ref_path_interval.first);
            step_handle_t end = rev_strand ? graph->get_previous_step(graph->get_step_at_position(path_handle, ref_path_interval.first))
                                           : graph->get_next_step(graph->get_step_at_position(path_handle, ref_path_interval.second));
            
            // walk the identified interval
            for (; step != end && mappings_matched != surj_path.mapping_size();
                 step = rev_strand ? graph->get_previous_step(step) : graph->get_next_step(step)) {
                const auto& pos = surj_path.mapping(mappings_matched).position();
                handle_t handle = graph->get_handle_of_step(step);
                if (graph->get_id(handle) == pos.node_id() &&
                    ((graph->get_is_reverse(handle) != pos.is_reverse()) == rev_strand)) {
                    // we found the next position we were expecting to
                    if (mappings_matched == 0) {
                        path_range_out.first = step;
                    }
                    path_range_out.second = step;
                    ++mappings_matched;
#ifdef debug_anchored_surject
                    cerr << "\tmatch at " << graph->get_id(handle) << " " << graph->get_is_reverse(handle) << endl;
#endif
                }
                else {
                    // we mismatched the path
                    if (mappings_matched) {
                        // return as if you hadn't matched at the start of this potential match
                        mappings_matched = 0;
                        // and go back to where we started on the path
                        // TODO: this is potentially quadratic, there are faster algorithms
                        step = path_range_out.first;
                    }
#ifdef debug_anchored_surject
                    cerr << "\tmismatch at " << graph->get_id(handle) << " " << graph->get_is_reverse(handle) << endl;
#endif
                }
            }
            
            if (mappings_matched != surj_path.mapping_size()) {
                cerr << "error: couldn't identify a path corresponding to surjected read " << source.name() << endl;
            }
        }
        else {
            // sentinel to indicate that surjection is unmapped
            path_range_out.first = path_range_out.second = graph->path_end(path_handle);
        }
        
        // transfer applicable metadata (including data that doesn't transit through multipath_alignment_t)
        surjected.set_name(source.name());
        surjected.set_read_group(source.read_group());
        surjected.set_sample_name(source.sample_name());
        surjected.set_mapping_quality(source.mapping_quality());
        if (source.has_fragment_next()) {
            *surjected.mutable_fragment_next() = source.fragment_next();
        }
        if (source.has_fragment_prev()) {
            *surjected.mutable_fragment_prev() = source.fragment_prev();
        }
        if (source.has_annotation()) {
            *surjected.mutable_annotation() = source.annotation();
        }
        
#ifdef debug_anchored_surject
        cerr << "concatenated and translated alignment " << pb2json(surjected) << endl;
#endif
        
        return surjected;
    }

    unordered_map<pair<path_handle_t, bool>, pair<vector<Surjector::path_chunk_t>, vector<pair<step_handle_t, step_handle_t>>>>
    Surjector::extract_overlapping_paths(const PathPositionHandleGraph* graph,
                                         const multipath_alignment_t& source,
                                         const unordered_set<path_handle_t>& surjection_paths,
                                         unordered_map<pair<path_handle_t, bool>, vector<tuple<size_t, size_t, int32_t>>>& connections_out) const {
        
        unordered_map<pair<path_handle_t, bool>, pair<vector<path_chunk_t>, vector<pair<step_handle_t, step_handle_t>>>> to_return;
        
        // reverse the connection edges for easy backwards lookup
        vector<vector<pair<size_t, int32_t>>> rev_connections(source.subpath_size());
        
        // compute the start of the read interval that corresponds to each mapping
        vector<vector<int64_t>> mapping_to_lengths(source.subpath_size());
        for (int64_t i = 0; i < source.subpath_size(); ++i) {
            mapping_to_lengths[i].resize(source.subpath(i).path().mapping_size(), 0);
            for (const auto& connection : source.subpath(i).connection()) {
                rev_connections[connection.next()].emplace_back(i, connection.score());
            }
        }
        for (int64_t i = 0; i < source.subpath_size(); ++i) {
            const auto& subpath = source.subpath(i);
            const auto& path = subpath.path();
            auto& subpath_to_length = mapping_to_lengths[i];
            int64_t thru_length = subpath_to_length.front();
            for (size_t j = 0; j < path.mapping_size(); ++j) {
                thru_length += mapping_to_length(path.mapping(j));
                if (j + 1 < path.mapping_size()) {
                    subpath_to_length[j + 1] = thru_length;
                }
            }
            for (auto n : subpath.next()) {
                mapping_to_lengths[n][0] = thru_length;
            }
            for (auto c : subpath.connection()) {
                mapping_to_lengths[c.next()][0] = thru_length;
            }
        }
        
        // map from (path, strand, subpath idx) to indexes among path chunks that have outgoing connections
        unordered_map<tuple<path_handle_t, bool, size_t>, vector<size_t>> connection_sources;
                
        // the mappings (subpath, mapping) that have already been associated with a step
        unordered_set<tuple<int64_t, int64_t, step_handle_t>> associated;
        for (int64_t i = 0; i < source.subpath_size(); ++i) {
            const auto& path = source.subpath(i).path();
            int64_t prefix_to_length = 0;
            for (int64_t j = 0; j < path.mapping_size(); ++j) {
                const auto& mapping = path.mapping(j);
                const auto& pos = mapping.position();
                handle_t handle = graph->get_handle(pos.node_id(), pos.is_reverse());
                graph->for_each_step_on_handle(handle, [&](const step_handle_t& step) {
                    
                    path_handle_t path_handle = graph->get_path_handle_of_step(step);
                    
                    if (!surjection_paths.count(path_handle) || associated.count(make_tuple(i, j, step))) {
                        // this is not on a path we're surjecting to, or we've already
                        // done it
                        return;
                    }
                    
#ifdef debug_multipath_surject
                    cerr << "starting new path DFS for subpath " << i << ", mapping " << j << ", path " << graph->get_path_name(path_handle) << ", step at " << graph->get_position_of_step(step) << endl;
#endif
                    
                    // do DFS starting from this mapping to find maximal path overlapping chunks
                    
                    // records of (subpath, mapping, next edge idx, step here)
                    // internal mappings are treated as having a single edge
                    vector<tuple<int64_t, int64_t, int64_t, step_handle_t>> stack;
                    stack.emplace_back(i, j, 0, step);
                    bool added_new_mappings = true;
                    while (!stack.empty()) {
                        
                        int64_t s_idx, m_idx, n_idx;
                        step_handle_t step_here;
                        tie(s_idx, m_idx, n_idx, step_here) = stack.back();
#ifdef debug_multipath_surject
                        cerr << "stack frame s " << s_idx << ", m " << m_idx << ", n " << n_idx << ", step at " << graph->get_position_of_step(step_here) << endl;
#endif
                        
                        const auto& subpath_here = source.subpath(s_idx);
                        const auto& path_here = subpath_here.path();
                        
                        if ((m_idx + 1 < path_here.mapping_size() && n_idx != 0)
                            || (m_idx + 1 == path_here.mapping_size() && (n_idx == subpath_here.next_size()
                                                                          || !subpath_here.connection().empty()))) {
                            // we've exhausted all of the outgoing adjacencies from this mapping
#ifdef debug_multipath_surject
                            cerr << "adjacencies exhausted or hit a connection" << endl;
#endif
                            if (added_new_mappings) {
                                
                                // a DFS traveresal has gone as far as possible, output the stack as a path
                                auto path_strand = make_pair(path_handle, handle != graph->get_handle_of_step(step));
                                auto& section_record = to_return[path_strand];
                                
                                if (m_idx + 1 == path_here.mapping_size() && !subpath_here.connection().empty()) {
                                    // record that connections leave this patch chunk
                                    connection_sources[make_tuple(path_handle, path_strand.second, s_idx)].push_back(section_record.first.size());
                                }
                                if (j == 0) {
                                    // translate connections into the indexes of their path chunks
                                    // note: if these are on different strands, they'll be ignored
                                    for (const auto& c : rev_connections[i]) {
                                        for (auto source : connection_sources[make_tuple(path_handle, path_strand.second, c.first)]) {
                                            connections_out[path_strand].emplace_back(source, section_record.first.size(),
                                                                                      c.second);
                                        }
                                    }
                                }
                                
                                // the interval of steps
                                section_record.second.emplace_back(get<3>(stack.front()), get<3>(stack.back()));
                                
                                section_record.first.emplace_back();
                                auto& chunk = section_record.first.back();
                                
                                // the aligned path
                                auto& path_chunk = chunk.second;
                                for (const auto& record : stack) {
                                    associated.emplace(get<0>(record), get<1>(record), get<3>(record));
                                    const auto& next_mapping = source.subpath(get<0>(record)).path().mapping(get<1>(record));
                                    const auto& next_pos = next_mapping.position();
                                    bool merged_mapping = false;
                                    if (path_chunk.mapping_size() != 0) {
                                        
                                        auto prev_mapping = path_chunk.mutable_mapping(path_chunk.mapping_size() - 1);
                                        const auto& prev_pos = prev_mapping->position();
                                        
                                        if (next_pos.node_id() == prev_pos.node_id() &&
                                            next_pos.is_reverse() == prev_pos.is_reverse() &&
                                            next_pos.offset() == prev_pos.offset() + mapping_from_length(*prev_mapping)) {
                                            // the next mapping is contiguous on a node with the previous one, we can merge
                                            // the two mappings into one
                                            
                                            auto prev_edit = prev_mapping->mutable_edit(prev_mapping->edit_size() - 1);
                                            const auto& next_edit = next_mapping.edit(0);
                                            if ((prev_edit->from_length() != 0) == (next_edit.from_length() != 0) &&
                                                (prev_edit->to_length() != 0) == (next_edit.to_length() != 0) &&
                                                prev_edit->sequence().empty() == next_edit.sequence().empty()) {
                                                
                                                prev_edit->set_from_length(prev_edit->from_length() + next_edit.from_length());
                                                prev_edit->set_to_length(prev_edit->to_length() + next_edit.to_length());
                                                prev_edit->set_sequence(prev_edit->sequence() + next_edit.sequence());
                                            }
                                            else {
                                                to_proto_edit(next_edit, *prev_mapping->add_edit());
                                            }
                                            for (size_t k = 1; k < next_mapping.edit_size(); ++k) {
                                                to_proto_edit(next_mapping.edit(k), *prev_mapping->add_edit());
                                            }
                                            
                                            merged_mapping = true;
                                        }
                                    }
                                    if (!merged_mapping) {
                                        // make a new mapping
                                        to_proto_mapping(next_mapping, *path_chunk.add_mapping());
                                    }
                                }
                                // the read interval
                                chunk.first.first = source.sequence().begin() + mapping_to_lengths[i][j];
                                chunk.first.second = chunk.first.first + path_to_length(path_chunk);
                                
                                // remember that we've already emitted all the mappings currently on the stack
                                added_new_mappings = false;
#ifdef debug_multipath_surject
                                cerr << "converted stack into path " << pb2json(path_chunk) << endl;
                                cerr << "read interval is " << (chunk.first.first - source.sequence().begin()) << ":" << (chunk.first.second - source.sequence().begin()) << " " << string(chunk.first.first, chunk.first.second) << endl;
#endif
                            }
                            
                            stack.pop_back();
                            continue;
                        }
                        
                        // mark that we have used this adjacency up
                        ++get<2>(stack.back());
                        
                        // get the indexes of the next mapping
                        const auto& mapping_here = path_here.mapping(m_idx);
                        int64_t next_s_idx, next_m_idx;
                        const path_mapping_t* next_mapping;
                        if (m_idx + 1 == path_here.mapping_size()) {
                            // mapping is at a subpath boundary
                            if (n_idx == 0) {
                                // check whether we branch to different segments of the same path
                                unordered_map<path_handle_t, vector<step_handle_t>> next_steps;
                                for (auto n : subpath_here.next()) {
                                    const auto& path = source.subpath(n).path();
                                    // search through path until we find an aligned base
                                    size_t p_idx = 0;
                                    while (p_idx + 1 < path.mapping_size() && mapping_from_length(path.mapping(p_idx)) == 0) {
                                        ++p_idx;
                                    }
                                    const auto& pos = path.mapping(p_idx).position();
                                    handle_t h = graph->get_handle(pos.node_id(), pos.is_reverse());
                                    graph->for_each_step_on_handle(h, [&](const step_handle_t& step) {
                                        next_steps[graph->get_path_handle_of_step(step)].push_back(step);
                                    });
                                }
                                bool branches_along_path = false;
                                for (pair<const path_handle_t, vector<step_handle_t>>& path_steps : next_steps) {
                                    sort(path_steps.second.begin(), path_steps.second.end());
                                    if (unique(path_steps.second.begin(), path_steps.second.end()) - path_steps.second.begin() > 1) {
                                        // the next subpaths along this branch point reach different places on the same
                                        // path, we have to avoid this so that the splicing logic will work
                                        branches_along_path = true;
                                        break;
                                    }
                                }
                                
                                if (branches_along_path) {
                                    // we'll prematurely end the DFS at this mapping
                                    get<2>(stack.back()) = subpath_here.next_size();
                                    continue;
                                }
                            }
                            next_s_idx = subpath_here.next(n_idx);
                            if (!rev_connections[next_s_idx].empty()) {
                                // we always break a path chunk at a connection
                                continue;
                            }
                            next_m_idx = 0;
                            next_mapping = &source.subpath(next_s_idx).path().mapping().front();
                        }
                        else {
                            // mapping is not at a subpath boundary
                            next_s_idx = s_idx;
                            next_m_idx = m_idx + 1;
                            next_mapping = &path_here.mapping(next_m_idx);
                        }
#ifdef debug_multipath_surject
                        cerr << "next s " << next_s_idx << ", m " << next_m_idx << endl;
#endif
                        
                        // check if the next position is consistent with the path we're walking
                        const auto& pos_here = mapping_here.position();
                        const auto& next_pos = next_mapping->position();
                        if (pos_here.node_id() == next_pos.node_id()
                            && pos_here.is_reverse() == next_pos.is_reverse()
                            && pos_here.offset() + mapping_from_length(mapping_here) == next_pos.offset()) {
                            // mappings are abutting within a node, don't leave the current step
                            stack.emplace_back(next_s_idx, next_m_idx, 0, step_here);
                            added_new_mappings = true;
                        }
                        else {
                            // mappings cross an edge in the graph
                            bool strand_rev = pos_here.is_reverse() != graph->get_is_reverse(graph->get_handle_of_step(step_here));
                            step_handle_t next_step = strand_rev ? graph->get_previous_step(step_here) : graph->get_next_step(step_here);
                            if (next_step != graph->path_end(path_handle) && next_step != graph->path_front_end(path_handle)) {
                                
                                handle_t next_handle = graph->get_handle_of_step(next_step);
                                if (graph->get_id(next_handle) == next_pos.node_id() &&
                                    (graph->get_is_reverse(next_handle) != next_pos.is_reverse()) == strand_rev) {
                                    // the next mapping is along the path how we would expect
                                    stack.emplace_back(next_s_idx, next_m_idx, 0, next_step);
                                    added_new_mappings = true;
                                }
                            }
                        }
                    }
                });
            }
        }
                
        return to_return;
    }
    
    unordered_map<pair<path_handle_t, bool>, pair<vector<Surjector::path_chunk_t>, vector<pair<step_handle_t, step_handle_t>>>>
    Surjector::extract_overlapping_paths(const PathPositionHandleGraph* graph, const Alignment& source,
                                         const unordered_set<path_handle_t>& surjection_paths) const {
        
        unordered_map<pair<path_handle_t, bool>, pair<vector<path_chunk_t>, vector<pair<step_handle_t, step_handle_t>>>> to_return;
        
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
                
                auto& path_chunks = to_return[make_pair(path_handle, path_strand)];
                
                if (extending_steps.count(make_pair(prev_step, path_strand))) {
                    // we are extending from the previous step, so we continue with the extension
                    
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

    void Surjector::filter_redundant_path_chunks(vector<path_chunk_t>& path_chunks,
                                                 vector<pair<step_handle_t, step_handle_t>>& ref_chunks,
                                                 vector<tuple<size_t, size_t, int32_t>>& connections) const {
        
        
#ifdef debug_filter_paths
        cerr << "filtering redundant path chunks" << endl;
#endif
        
        assert(path_chunks.size() == ref_chunks.size());
        vector<size_t> order(path_chunks.size(), 0);
        for (size_t i = 1; i < order.size(); ++i) {
            order[i] = i;
        }
        
        // convert connections to adjacency lists
        vector<vector<pair<size_t, int32_t>>> inward_connections(path_chunks.size()), outward_connections(path_chunks.size());
        for (const auto& connection : connections) {
            inward_connections[get<1>(connection)].emplace_back(get<0>(connection), get<2>(connection));
            outward_connections[get<0>(connection)].emplace_back(get<1>(connection), get<2>(connection));
        }
        
        // sort the adjacency lists for easy determination of subsets
        for (auto& adj : inward_connections) {
            sort(adj.begin(), adj.end());
        }
        for (auto& adj : outward_connections) {
            sort(adj.begin(), adj.end());
        }
        
#ifdef debug_filter_paths
        cerr << "original order for chunks" << endl;
        for (size_t i = 0; i < path_chunks.size(); ++i) {
            cerr << i << ": " << string(path_chunks[i].first.first, path_chunks[i].first.second) << " " << pb2json(path_chunks[i].second) << endl;
        }
        cerr << "connections" << endl;
        for (size_t i = 0; i < outward_connections.size(); ++i) {
            cerr << i << ":";
            for (const auto& c : outward_connections[i]) {
                cerr << " (" << c.first << " " << c.second << ")";
            }
            cerr << endl;
        }
#endif
        
        // test it one adjacency list entry is a subset of another (assumes sort)
        auto is_subset = [](const vector<pair<size_t, int32_t>>& sub, const vector<pair<size_t, int32_t>>& super) {
            size_t i = 0, j = 0;
            while (i < sub.size() && j < super.size()) {
                if (sub[i] == super[j]) {
                    ++i;
                    ++j;
                }
                else {
                    ++j;
                }
            }
            return (i == sub.size());
        };
        
        // order the path chunks by the left index of their read interval
        // and break ties in favor of longer intervals (so that the filteree
        // will come later in the vector as we expect)
        // among paths that are still tied, prioritize path chunks with more
        // connections
        stable_sort(order.begin(), order.end(), [&](size_t i, size_t j) {
            auto& chunk1 = path_chunks[i];
            auto& chunk2 = path_chunks[j];
            return (chunk1.first.first < chunk2.first.first ||
                    (chunk1.first.first == chunk2.first.first && chunk1.first.second > chunk2.first.second) ||
                    (chunk1.first.first == chunk2.first.first && chunk1.first.second == chunk2.first.second &&
                     (inward_connections[i].size() + outward_connections[i].size()
                      > inward_connections[j].size() + outward_connections[j].size())));
        });
        
#ifdef debug_filter_paths
        cerr << "sort order for chunks" << endl;
        for (auto i : order) {
            cerr << i << ": " << string(path_chunks[i].first.first, path_chunks[i].first.second) << " " << pb2json(path_chunks[i].second) << endl;
        }
#endif
        
        vector<bool> redundant(path_chunks.size(), false);
        
        // a heap where the top always points to the leftmost end of a read interval
        auto cmp = [&](int64_t i, int64_t j) {
            return path_chunks[i].first.second > path_chunks[j].first.second;
        };
        vector<int64_t> curr_chunks;
        
        for (int64_t i = 0; i < order.size(); ++i) {
            auto& chunk_here = path_chunks[order[i]];
#ifdef debug_filter_paths
            cerr << "looking for overlapping chunks for " << order[i] << endl;
            cerr << string(chunk_here.first.first, chunk_here.first.second) << " " << pb2json(chunk_here.second) << endl;
#endif
            // remove items from the heap if they are outside the window of this read interval
            while (!curr_chunks.empty() && path_chunks[curr_chunks.front()].first.second <= chunk_here.first.first) {
                pop_heap(curr_chunks.begin(), curr_chunks.end(), cmp);
                curr_chunks.pop_back();
            }
            
            for (auto j : curr_chunks) {
                if (path_chunks[j].first.first > chunk_here.first.first ||
                    path_chunks[j].first.second < chunk_here.first.second) {
                    // doesn't contain the right read interval
                    continue;
                }
                
                auto& chunk_over = path_chunks[j];
                auto remaining = chunk_here.first.first - chunk_over.first.first;
#ifdef debug_filter_paths
                cerr << "overlap candidate " << j << endl;
                cerr << string(chunk_over.first.first, chunk_over.first.second) << " " << pb2json(chunk_over.second) << endl;
                cerr << "at relative read offset " << remaining << endl;
#endif
                
                // walk the part of the overlapping path that comes before the path here
                int64_t m_over_idx = 0, e_over_idx = 0;
                while (m_over_idx < chunk_over.second.mapping_size()
                       && remaining >= chunk_over.second.mapping(m_over_idx).edit(e_over_idx).to_length()
                       && remaining > 0) {
                    
                    remaining -= chunk_over.second.mapping(m_over_idx).edit(e_over_idx).to_length();
                    
#ifdef debug_filter_paths
                    cerr << "walk down overlapper before match at " << m_over_idx << " " << e_over_idx << endl;
#endif
                    
                    ++e_over_idx;
                    if (e_over_idx == chunk_over.second.mapping(m_over_idx).edit_size()) {
                        ++m_over_idx;
                        e_over_idx = 0;
                    }
                }
                
                // we might need to walk another subpath of with to length of 0 to get to the start
                // of the short path
                while (m_over_idx < chunk_over.second.mapping_size()
                       && chunk_over.second.mapping(m_over_idx).edit(e_over_idx).to_length() == 0
                       && (chunk_over.second.mapping(m_over_idx).position().node_id()
                           != chunk_here.second.mapping(0).position().node_id())
                       && (chunk_over.second.mapping(m_over_idx).position().is_reverse()
                           != chunk_here.second.mapping(0).position().is_reverse())
                       && (chunk_over.second.mapping(m_over_idx).position().offset()
                           != chunk_here.second.mapping(0).position().offset())) {
#ifdef debug_filter_paths
                    cerr << "walking forward through a deletion to find match" << endl;
#endif
                    ++e_over_idx;
                    if (e_over_idx == chunk_over.second.mapping(m_over_idx).edit_size()) {
                        ++m_over_idx;
                        e_over_idx = 0;
                    }
                }
                
#ifdef debug_filter_paths
                cerr << "search for overlap begins at over idx " << m_over_idx << " " << e_over_idx << endl;
#endif
                
                // we'll only consider it to match the paths meet at a mapping boundary
                bool matches = (remaining == 0 && e_over_idx == 0
                                && m_over_idx < chunk_over.second.mapping_size());
                if (!matches) {
#ifdef debug_filter_paths
                    cerr << "shorter path not at an internal mapping boundary" << endl;
#endif
                    continue;
                }
                
                bool shares_start = m_over_idx == 0;
                
                // try to walk the part of the path where they overlap
                int64_t m_here_idx = 0, e_here_idx = 0;
                while (m_over_idx < chunk_over.second.mapping_size() &&
                       m_here_idx < chunk_here.second.mapping_size()) {
#ifdef debug_filter_paths
                    cerr << "looking for match at " << m_over_idx << " " << e_over_idx << ", " << m_here_idx << " " << e_here_idx << endl;
#endif
                    
                    if (e_here_idx == 0) {
                        const auto& pos_over = chunk_over.second.mapping(m_over_idx).position();
                        const auto& pos_here = chunk_here.second.mapping(m_here_idx).position();
                        if (pos_here.node_id() != pos_over.node_id()
                            || pos_here.is_reverse() != pos_over.is_reverse()
                            || pos_here.offset() != pos_over.offset()) {
#ifdef debug_filter_paths
                            cerr << "mappings not at the same position" << endl;
#endif
                            matches = false;
                            break;
                        }
                    }
                    const auto& edit_over = chunk_over.second.mapping(m_over_idx).edit(e_over_idx);
                    const auto& edit_here = chunk_here.second.mapping(m_here_idx).edit(e_here_idx);
                    if (edit_here.from_length() != edit_over.from_length()
                        || edit_here.to_length() != edit_here.to_length()) {
                        // note: we don't need to worry about edit sequence because we know we're at
                        // the same read interval
                        
                        // the edits don't match
#ifdef debug_filter_paths
                        cerr << "edits don't match" << endl;
#endif
                        matches = false;
                        break;
                    }
                    
                    ++e_over_idx;
                    ++e_here_idx;
                    if (e_over_idx == chunk_over.second.mapping(m_over_idx).edit_size()) {
                        if (e_here_idx != chunk_here.second.mapping(m_here_idx).edit_size()) {
#ifdef debug_filter_paths
                            cerr << "mapping boundaries don't occur at the same place" << endl;
#endif
                            matches = false;
                            break;
                        }
                        ++m_over_idx;
                        e_over_idx = 0;
                    }
                    if (e_here_idx == chunk_here.second.mapping(m_here_idx).edit_size()) {
                        ++m_here_idx;
                        e_here_idx = 0;
                    }
                }
                
                                
                if (matches && m_here_idx == chunk_here.second.mapping_size()) {
                    
                    if (shares_start) {
#ifdef debug_filter_paths
                        cerr << "checking shared inward connections" << endl;
#endif
                        if (!is_subset(inward_connections[order[i]], inward_connections[j])) {
#ifdef debug_filter_paths
                            cerr << "connections are non-redundant" << endl;
#endif
                            // the chunk has unique inward connections, we can't eliminate it
                            continue;
                        }
                    }
                    else if (!inward_connections[order[i]].empty()) {
#ifdef debug_filter_paths
                        cerr << "has internal inward connections" << endl;
#endif
                        // has inward connections to the interior of the chunk
                        continue;
                    }
                    if (m_over_idx == chunk_over.second.mapping_size()) {
#ifdef debug_filter_paths
                        cerr << "checking shared outward connections" << endl;
#endif
                        if (!is_subset(outward_connections[order[i]], outward_connections[j])) {
#ifdef debug_filter_paths
                            cerr << "connections are non-redundant" << endl;
#endif
                            // the chunk has unique outward connections, we can't eliminate it
                            continue;
                        }
                    }
                    else if (!outward_connections[order[i]].empty()) {
#ifdef debug_filter_paths
                        cerr << "has internal outward connections" << endl;
#endif
                        // has outward connections to the interior of the chunk
                        continue;
                    }
                    
                    // the whole path matches an earlier, longer path
                    redundant[order[i]] = true;
#ifdef debug_filter_paths
                    cerr << "marking path chunk " << order[i] << " redundant" << endl;
#endif
                    break;
                }
            }
            
            // we only need to look at nonredundant chunks on later iterations
            if (!redundant[order[i]]) {
                curr_chunks.push_back(order[i]);
                push_heap(curr_chunks.begin(), curr_chunks.end());
            }
        }
        
        // filter down to nonredundant paths
        vector<size_t> removed_before(path_chunks.size() + 1, 0);
        for (size_t i = 0; i < path_chunks.size(); ++i) {
            if (redundant[i]) {
#ifdef debug_filter_paths
                cerr << "filtering path chunk " << i << ": " << string(path_chunks[i].first.first, path_chunks[i].first.second) << " " << pb2json(path_chunks[i].second) << endl;
#endif
                ++removed_before[i];
            }
            else if (removed_before[i]) {
                path_chunks[i - removed_before[i]] = move(path_chunks[i]);
                ref_chunks[i - removed_before[i]] = move(ref_chunks[i]);
            }
            removed_before[i + 1] = removed_before[i];
        }
        if (removed_before.back() != 0) {
            path_chunks.resize(path_chunks.size() - removed_before.back());
            ref_chunks.resize(ref_chunks.size() - removed_before.back());
        }
        
        // update the indexes on connections and remove redundant ones
        size_t removed_so_far = 0;
        for (size_t i = 0; i < connections.size(); ++i) {
            auto& connection = connections[i];
            if (redundant[get<0>(connection)] || redundant[get<1>(connection)]) {
                ++removed_so_far;
            }
            else {
#ifdef debug_filter_paths
                cerr << "updating connection " << get<0>(connection) << " -> " << get<1>(connection) << endl;
#endif
                get<0>(connection) -= removed_before[get<0>(connection)];
                get<1>(connection) -= removed_before[get<1>(connection)];
                if (removed_so_far) {
                    connections[i - removed_so_far] = move(connection);
                }
            }
        }
        if (removed_so_far) {
            connections.resize(connections.size() - removed_so_far);
        }
        
        // sort the path chunks in lexicographic order like the downstream code expects
        
        if (!is_sorted(path_chunks.begin(), path_chunks.end(),
                       [&](path_chunk_t& a, path_chunk_t& b) { return a.first < b.first; })) {
            
#ifdef debug_filter_paths
            cerr << "putting path chunks in lexicographic order" << endl;
#endif
            
            // compute which index the chunks should end up in
            for (size_t i = 0; i < path_chunks.size(); ++i) {
                order[i] = i;
            }
            order.resize(path_chunks.size());
            stable_sort(order.begin(), order.end(), [&](size_t i, size_t j) {
                return path_chunks[i].first < path_chunks[j].first;
            });
            vector<size_t> index(order.size());
            for (size_t i = 0; i < order.size(); ++i) {
                index[order[i]] = i;
            }
            
            // update the indexes of the connections
            for (auto& connection : connections) {
                get<0>(connection) = index[get<0>(connection)];
                get<1>(connection) = index[get<1>(connection)];
            }
            
            // and co-sort the vectors into the computed indexes
            for (size_t i = 0; i < index.size(); ++i) {
                while (index[i] != i) {
                    std::swap(path_chunks[i], path_chunks[index[i]]);
                    std::swap(ref_chunks[i], ref_chunks[index[i]]);
                    std::swap(index[i], index[index[i]]);
                }
            }
        }
    }
    
    pair<size_t, size_t>
    Surjector::compute_path_interval(const PathPositionHandleGraph* graph, const Alignment& source, path_handle_t path_handle,
                                     bool rev_strand, const vector<path_chunk_t>& path_chunks) const {
        
        pair<size_t, size_t> interval(numeric_limits<size_t>::max(), numeric_limits<size_t>::min());
        
        size_t path_length = graph->get_path_length(path_handle);
        
        for (const auto& path_chunk : path_chunks) {
            
            size_t left_overhang = (get_aligner()->longest_detectable_gap(source, path_chunk.first.first)
                                    + (path_chunk.first.first - source.sequence().begin()));
            
            size_t right_overhang = (get_aligner()->longest_detectable_gap(source, path_chunk.first.second)
                                     + (source.sequence().end() - path_chunk.first.second));
            
            const Position& first_pos = path_chunk.second.mapping(0).position();
            handle_t first_handle = graph->get_handle(first_pos.node_id(), first_pos.is_reverse());
            for (const step_handle_t& step : graph->steps_of_handle(first_handle)) {
                
                if (graph->get_path_handle_of_step(step) != path_handle ||
                    (first_pos.is_reverse() != graph->get_is_reverse(graph->get_handle_of_step(step))) != rev_strand) {
                    // this step isn't on the path strand we're considering
                    continue;
                }
                
                if (rev_strand) {
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
                
                if (graph->get_path_handle_of_step(step) != path_handle ||
                    (final_pos.is_reverse() != graph->get_is_reverse(graph->get_handle_of_step(step))) != rev_strand) {
                    // this step isn't on the path we're considering
                    continue;
                }
                
                if (rev_strand) {
                    size_t path_offset = (graph->get_position_of_step(step)
                                          + graph->get_length(final_handle)
                                          - final_pos.offset()
                                          - mapping_from_length(final_mapping));
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
        
        if (graph->get_position_of_step(end) <= last && end != graph->path_end(path_handle)) {
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

    void Surjector::set_path_position(const PathPositionHandleGraph* graph, const pos_t& init_surj_pos, const pos_t& final_surj_pos,
                                      const step_handle_t& range_begin, const step_handle_t& range_end,
                                      bool rev_strand, string& path_name_out, int64_t& path_pos_out, bool& path_rev_out) const {

        
        assert(graph->get_path_handle_of_step(range_begin) == graph->get_path_handle_of_step(range_end));
        
        if (range_begin == graph->path_end(graph->get_path_handle_of_step(range_begin))
            && range_begin == graph->path_end(graph->get_path_handle_of_step(range_end))) {
            // sentinel for unmapped
            path_name_out = "";
            path_pos_out = -1;
            path_rev_out = false;
        }
        else {
#if defined(debug_anchored_surject) || defined(debug_spliced_surject)
            cerr << "setting position with initial position " << init_surj_pos << " and final position " << final_surj_pos << " based on range:" << endl;
            cerr << "\tbegin: id " << graph->get_id(graph->get_handle_of_step(range_begin)) << ", rev " << graph->get_is_reverse(graph->get_handle_of_step(range_begin)) << ", pos " << graph->get_position_of_step(range_begin) << endl;
            cerr << "\tend: id " << graph->get_id(graph->get_handle_of_step(range_end)) << ", rev " << graph->get_is_reverse(graph->get_handle_of_step(range_end)) << ", pos " << graph->get_position_of_step(range_end) << endl;
#endif
            
            // the path name
            path_name_out = graph->get_path_name(graph->get_path_handle_of_step(range_begin));
            path_rev_out = rev_strand;
            
            // are we on the reverse strand?
            size_t path_pos_begin = graph->get_position_of_step(range_begin);
            size_t path_pos_end = graph->get_position_of_step(range_end);
            
            // the path offset
            if (rev_strand) {
                path_pos_out = (path_pos_end + graph->get_length(graph->get_handle_of_step(range_end))
                                - offset(final_surj_pos));
            }
            else {
                path_pos_out = path_pos_begin + offset(init_surj_pos);
            }
        }
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

    multipath_alignment_t Surjector::make_null_mp_alignment(const string& src_sequence,
                                                            const string& src_quality) {
        multipath_alignment_t null;
        null.set_sequence(src_sequence);
        null.set_quality(src_quality);
        return null;
    }
}


