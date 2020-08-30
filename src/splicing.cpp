/**
 * \file splicing.cpp
 *
 * Implements SpliceRegion and some other splicing tools
 *
 */

#include "splicing.hpp"

//#define debug_splice_region
//#define debug_trimming

namespace vg {
   
SpliceRegion::SpliceRegion(const pos_t& seed_pos, bool search_left, int64_t search_dist,
                           const HandleGraph& graph,
                           const DinucleotideMachine& dinuc_machine,
                           const vector<string>& splice_motifs) {
    
#ifdef debug_splice_region
    cerr << "constructing splice region" << endl;
#endif
    
    // add a buffer of 2 bases for the dinucleotide itself
    search_dist += 2;
    
    // ensure that all of the motifs get an entry
    for (const auto& motif : splice_motifs) {
        motif_matches[motif];
    }
    
    // TODO: make a splice table into an independent object?
    
    vector<string> rev_motifs;
    auto motifs = &splice_motifs;
    if (search_left) {
        // we will actually need to look for the motifs in reverse
        rev_motifs.reserve(splice_motifs.size());
        for (const auto& motif : splice_motifs) {
            rev_motifs.emplace_back(motif.rbegin(), motif.rend());
        }
        motifs = &rev_motifs;
    }
    
    // extract the subgraph and initialize the DP structure
    
#ifdef debug_splice_region
    cerr << "init subgraph with pos " << seed_pos << endl;
#endif
    
    IncrementalSubgraph subgraph(graph, seed_pos, search_left, search_dist);
    
    vector<pair<handle_t, vector<uint32_t>>> dinuc_states;
    handle_t handle = subgraph.handle_at_order(0);
    dinuc_states.emplace_back(handle, vector<uint32_t>(subgraph.get_length(handle),
                                                       dinuc_machine.init_state()));
    
    while (subgraph.is_extendable()) {
        handle = subgraph.extend();
        dinuc_states.emplace_back(handle, vector<uint32_t>(subgraph.get_length(handle),
                                                           dinuc_machine.init_state()));
#ifdef debug_splice_region
        cerr << "extract " << graph.get_id(subgraph.get_underlying_handle(handle)) << " " << graph.get_is_reverse(subgraph.get_underlying_handle(handle)) << endl;
#endif
    }
    int64_t incr = search_left ? -1 : 1;
    
    // check if we match any motifs at this location and if so remember it
    auto record_motif_matches = [&](handle_t handle, int64_t j,
                                    const vector<uint32_t>& states) {
        for (size_t i = 0; i < splice_motifs.size(); ++i) {
            if (dinuc_machine.matches(states[j], (*motifs)[i])) {
                if ((j == 0 && !search_left) || (j + 1 == states.size() && search_left)) {
                    // we need to cross a node boundary to backtrack
                    subgraph.follow_edges(handle, !search_left, [&](const handle_t& prev) {
                        if (search_left) {
                            if (subgraph.get_base(prev, 0) == (*motifs)[i].front()) {
                                auto underlying = subgraph.get_underlying_handle(prev);
                                pos_t pos(graph.get_id(underlying), graph.get_is_reverse(underlying), 1);
                                int64_t trav_dist = subgraph.distance_from_start(prev) + subgraph.get_length(prev) - 1;
                                motif_matches[splice_motifs[i]].emplace_back(pos, trav_dist);
                            }
                        }
                        else {
                            size_t k = subgraph.get_length(prev) - 1;
                            if (subgraph.get_base(prev, k) == (*motifs)[i].front()) {
                                auto underlying = subgraph.get_underlying_handle(prev);
                                pos_t pos(graph.get_id(underlying), graph.get_is_reverse(underlying), k);
                                int64_t trav_dist = subgraph.distance_from_start(prev) + k;
                                motif_matches[splice_motifs[i]].emplace_back(pos, trav_dist);
                            }
                        }
                    });
                }
                else {
                    auto underlying = subgraph.get_underlying_handle(handle);
                    pos_t pos(graph.get_id(underlying), graph.get_is_reverse(underlying), j - 2 * incr + !search_left);
                    int64_t trav_dist = subgraph.distance_from_start(handle);
                    if (search_left) {
                        trav_dist += states.size() - j - 2;
                    }
                    else {
                        trav_dist += j - 1;
                    }
#ifdef debug_splice_region
                    cerr << "record match to " << splice_motifs[i] << " at " << pos << ", dist " << trav_dist << endl;
#endif
                    motif_matches[splice_motifs[i]].emplace_back(pos, trav_dist);
                }
            }
        }
    };
    
    
    // now actually do the DP
    for (size_t i = 0; i < dinuc_states.size(); ++i) {
        
        handle_t here = dinuc_states[i].first;
        vector<uint32_t>& states = dinuc_states[i].second;
        string seq = subgraph.get_sequence(here);
        
        // determine where we'll start iterating from
        int64_t j;
        if (i == 0) {
            j = search_left ? offset(seed_pos) - 1 : offset(seed_pos);
        }
        else {
            j = search_left ? seq.size() - 1 : 0;
        }
        
        // determine the bounds of the iteration
        int64_t prev_dist = subgraph.distance_from_start(here);
        int64_t left_end = 0;
        int64_t right_end = seq.size();
        if (prev_dist + seq.size() >= search_dist) {
            if (search_left) {
                left_end = prev_dist + seq.size() - search_dist;
            }
            else {
                right_end = search_dist - prev_dist;
            }
        }
        
#ifdef debug_splice_region
        cerr << "node number " << i << ", iteration bounds: j = " << j << ", incr = " << incr << ", left end = " << left_end << ", right end " << right_end << ", node len = " << seq.size() << endl;
#endif
        // are we starting at the boundary of a node?
        if ((j == 0 && !search_left) || (j == seq.size() - 1 && search_left)) {
            // merge all of the incoming transition states
            subgraph.follow_edges(here, !search_left, [&](const handle_t& prev) {
                vector<uint32_t>& incoming_states = dinuc_states[subgraph.order_of(prev)].second;
                uint32_t incoming = search_left ? incoming_states.front() : incoming_states.back();
                states[j] = dinuc_machine.merge_state(states[j], dinuc_machine.update_state(incoming, seq[j]));
            });
            record_motif_matches(here, j, states);
            j += incr;
        }
        
        // carry forward the transitions to the end of the node
        for (; j >= left_end && j < right_end; j += incr) {
            states[j] = dinuc_machine.update_state(states[j - incr], seq[j]);
            record_motif_matches(here, j, states);
        }
    }
}

const vector<pair<pos_t, int64_t>>& SpliceRegion::candidate_splice_sites(const string& motif) const {
    return motif_matches.at(motif);
}

tuple<pos_t, int64_t, int32_t> trimmed_end(const Alignment& aln, int64_t len, bool from_end,
                                           const HandleGraph& graph, const GSSWAligner& aligner) {
    
#ifdef debug_trimming
    cerr << "trimming alignment " << pb2json(aln) << " by " << len << ", from end? " << from_end << endl;
#endif
    
    const Path& path = aln.path();
    Path dummy_path;
    
    tuple<pos_t, int64_t, int32_t> return_val;
    
    bool copied_full_path = false;
    if (path.mapping_size()) {
        if (from_end) {
            const Mapping& final_mapping = path.mapping(path.mapping_size() - 1);
            if (final_mapping.edit(final_mapping.edit_size() - 1).from_length() == 0) {
                // we have to walk further to skip the softclips
                len += final_mapping.edit(final_mapping.edit_size() - 1).to_length();
            }
            int64_t i = path.mapping_size() - 1;
            while (i >= 0 && (len > mapping_to_length(path.mapping(i))
                              || mapping_from_length(path.mapping(i)) == 0)) {
                auto to_length = mapping_to_length(path.mapping(i));
                len = max<int64_t>(len - to_length, 0);
                get<1>(return_val) += to_length;
                
#ifdef debug_trimming
                cerr << "after mapping " << i << ", remaining length " << len << endl;
#endif
                
                Mapping* dummy_mapping = dummy_path.add_mapping();
                *dummy_mapping->mutable_edit() = path.mapping(i).edit();
                --i;
            }
            if (i < 0) {
#ifdef debug_trimming
                cerr << "walked entire path" << endl;
#endif
                get<0>(return_val) = initial_position(path);
                get<1>(return_val) = path_to_length(path);
                dummy_path = path;
                copied_full_path = true;
            }
            else {
                const Mapping& mapping = path.mapping(i);
                int64_t j = mapping.edit_size() - 1;
                int64_t from_length = 0;
                Mapping* dummy_mapping = nullptr;
                while (j >= 0 && (len > mapping.edit(j).to_length()
                                  || mapping.edit(j).from_length() == 0)) {
                    auto to_length = mapping.edit(j).to_length();
                    len = max<int64_t>(len - to_length, 0);
                    get<1>(return_val) += to_length;
                    from_length += mapping.edit(j).from_length();
                    
#ifdef debug_trimming
                    cerr << "after edit " << j << ", remaining length " << len << endl;
#endif
                    if (!dummy_mapping) {
                        dummy_mapping = dummy_path.add_mapping();
                    }
                    *dummy_mapping->add_edit() = mapping.edit(j);
                    --j;
                }
                if (j >= 0 && len > 0) {
                    auto last_from_length = (len * mapping.edit(j).from_length()) / mapping.edit(j).to_length();
                    get<1>(return_val) += len;
                    from_length += last_from_length;
                    
#ifdef debug_trimming
                    cerr << "handling final (partial) edit with to length " << len << ", from length " << last_from_length << endl;
#endif
                    
                    if (!dummy_mapping) {
                        dummy_mapping = dummy_path.add_mapping();
                    }
                    Edit* dummy_edit = dummy_mapping->add_edit();
                    dummy_edit->set_from_length(last_from_length);
                    dummy_edit->set_to_length(len);
                    if (!mapping.edit(j).sequence().empty()) {
                        dummy_edit->set_sequence(mapping.edit(j).sequence().substr(mapping.edit(j).to_length() - len, len));
                    }
                }
                const Position& position = mapping.position();
                get_id(get<0>(return_val)) = position.node_id();
                get_is_rev(get<0>(return_val)) = position.is_reverse();
                get_offset(get<0>(return_val)) = position.offset() + mapping_from_length(mapping) - from_length;
            }
        }
        else {
            if (path.mapping(0).edit(0).from_length() == 0) {
                // we have to walk further to skip the softclips
                len += path.mapping(0).edit(0).to_length();
            }
            int64_t i = 0;
            while (i < path.mapping_size() && (len > mapping_to_length(path.mapping(i))
                                               || mapping_from_length(path.mapping(i)) == 0 )) {
                
                auto to_length = mapping_to_length(path.mapping(i));
                len = max<int64_t>(len - to_length, 0);
                get<1>(return_val) += to_length;
                
#ifdef debug_trimming
                cerr << "after mapping " << i << ", remaining length " << len << endl;
#endif
                
                Mapping* dummy_mapping = dummy_path.add_mapping();
                *dummy_mapping->mutable_edit() = path.mapping(i).edit();
                ++i;
            }
            if (i == path.mapping_size()) {
#ifdef debug_trimming
                cerr << "walked entire path" << endl;
#endif
                get<0>(return_val) = final_position(path);
                get<1>(return_val) = path_to_length(path);
                dummy_path = path;
                copied_full_path = true;
            }
            else {
                const Mapping& mapping = path.mapping(i);
                int64_t j = 0;
                int64_t from_length = 0;
                Mapping* dummy_mapping = nullptr;
                while (j < mapping.edit_size() && (len > mapping.edit(j).to_length()
                                                   || mapping.edit(j).from_length() == 0)) {
                    auto to_length = mapping.edit(j).to_length();
                    len = max<int64_t>(len - to_length, 0);
                    get<1>(return_val) += to_length;
                    from_length += mapping.edit(j).from_length();
                    
#ifdef debug_trimming
                    cerr << "after edit " << j << ", remaining length " << len << endl;
#endif
                    
                    if (!dummy_mapping) {
                        dummy_mapping = dummy_path.add_mapping();
                    }
                    *dummy_mapping->add_edit() = mapping.edit(j);
                    ++j;
                }
                if (j != mapping.edit_size() && len > 0) {
                    auto last_from_length = (len * mapping.edit(j).from_length()) / mapping.edit(j).to_length();
                    get<1>(return_val) += len;
                    from_length += last_from_length;
                    
#ifdef debug_trimming
                    cerr << "handling final (partial) edit with to length " << len << ", from length " << last_from_length << endl;
#endif
                    
                    if (!dummy_mapping) {
                        dummy_mapping = dummy_path.add_mapping();
                    }
                    Edit* dummy_edit = dummy_mapping->add_edit();
                    dummy_edit->set_from_length(last_from_length);
                    dummy_edit->set_to_length(len);
                    if (!mapping.edit(j).sequence().empty()) {
                        dummy_edit->set_sequence(mapping.edit(j).sequence().substr(0, len));
                    }
                }
                const Position& position = mapping.position();
                get_id(get<0>(return_val)) = position.node_id();
                get_is_rev(get<0>(return_val)) = position.is_reverse();
                get_offset(get<0>(return_val)) = position.offset() + from_length;
            }
        }
    }
    
    string::const_iterator begin;
    if (from_end) {
        begin = aln.sequence().end() - get<1>(return_val);
        // TODO: kind of inelegant
        if (!copied_full_path) {
            // the path was built in reverse, flip around
            for (size_t i = 0, end = dummy_path.mapping_size() / 2; i < end; ++i) {
                swap(*dummy_path.mutable_mapping(i),
                     *dummy_path.mutable_mapping(dummy_path.mapping_size() - i - 1));
            }
            // the final mapping was also built in reverse
            Mapping* mapping = dummy_path.mutable_mapping(0);
            for (size_t i = 0, end = mapping->edit_size() / 2; i < end; ++i) {
                swap(*mapping->mutable_edit(i),
                     *mapping->mutable_edit(mapping->edit_size() - i - 1));
            }
        }
    }
    else {
        begin = aln.sequence().begin();
    }
#ifdef debug_trimming
    cerr << "scoring trimmed subpath " << pb2json(dummy_path) << ", with substring " << (begin - aln.sequence().begin()) << ":" << (begin - aln.sequence().begin()) + get<1>(return_val) << endl;
#endif
    
    // TODO: refactor so we can either use the subgraph or use score_mismatch
    // rather than needing get_handle calls
    get<2>(return_val) = aligner.score_partial_alignment(aln, graph, dummy_path, begin);
    
    return return_val;
}

}
