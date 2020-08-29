/**
 * \file splice_region.cpp
 *
 * Implements SpliceRegion
 *
 */

#include "splice_region.hpp"

//#define debug_splice_region

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

}
