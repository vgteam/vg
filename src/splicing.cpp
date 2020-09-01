/**
 * \file splicing.cpp
 *
 * Implements SpliceRegion and some other splicing tools
 *
 */

#include "splicing.hpp"

//#define debug_splice_region
//#define debug_trimming

#ifdef debug_splice_region
#import <bitset>
#endif

namespace vg {

SpliceMotifs::SpliceMotifs(const GSSWAligner& scorer) {
    
    vector<tuple<string, string, double>> default_motifs;
    default_motifs.emplace_back("GT", "AG", 0.9924);
    default_motifs.emplace_back("GC", "AG", 0.0069);
    default_motifs.emplace_back("AT", "AC", 0.0005);
    init(default_motifs, scorer);
}

SpliceMotifs::SpliceMotifs(const vector<tuple<string, string, double>>& motifs,
                           const GSSWAligner& scorer) {
    init(motifs, scorer);
}

size_t SpliceMotifs::size() const {
    return data.size();
}


const string& SpliceMotifs::oriented_motif(size_t motif_num, bool left_side) const {
    return left_side ? get<1>(data[motif_num]) : get<0>(data[motif_num]);
}

int32_t SpliceMotifs::score(size_t motif_num) const {
    return get<2>(data[motif_num]);
}

void SpliceMotifs::update_scoring(const GSSWAligner& scorer) {
    init(unaltered_data, scorer);
}

void SpliceMotifs::init(const vector<tuple<string, string, double>>& motifs,
                        const GSSWAligner& scorer) {
    
    // TODO: does this normalization to 1 make sense?
    double total_frequency = 0.0;
    for (const auto& record : motifs) {
        if (get<0>(record).size() != 2 || get<1>(record).size() != 2) {
             cerr << "error:[SpliceMotifs] Splice motif " << get<0>(record) << "-" << get<1>(record) << " is not a pair of dinucleotides." << endl;
        }
        if (get<2>(record) < 0.0 || get<2>(record) > 1.0) {
            cerr << "error:[SpliceMotifs] Frequency of splice motif " << get<0>(record) << "-" << get<1>(record) << " given as " << get<2>(record) << ". Must be a number between 0 and 1." << endl;
        }
        total_frequency += get<2>(record);
    }
    // a little slop for numerical imprecision
    if (total_frequency > 1.000001) {
        cerr << "error:[SpliceMotifs] Frequency of splice motifs sum to " << total_frequency << ". Must be a number between 0 and 1." << endl;
    }
    
    // in case we're resetting
    data.clear();
    unaltered_data = motifs;
    
#ifdef debug_splice_region
    cerr << "recording splice table" << endl;
#endif
    
    data.reserve(motifs.size());
    for (const auto& record : motifs) {
        data.emplace_back();
        get<0>(data.back()) = get<0>(record);
        // reverse the second string because it's encountered in reverse when going into
        // an intron
        get<1>(data.back()) = string(get<1>(record).rbegin(), get<1>(record).rend());
        // convert frequency to a log likelihood
        get<2>(data.back()) = int32_t(round(log(get<2>(record)) / scorer.log_base));
#ifdef debug_splice_region
        cerr << "\t" << get<0>(data.back()) << "\t" << get<1>(data.back()) << "\t" << get<2>(data.back()) << endl;
#endif
    }
}

SpliceRegion::SpliceRegion(const pos_t& seed_pos, bool search_left, int64_t search_dist,
                           const HandleGraph& graph,
                           const DinucleotideMachine& dinuc_machine,
                           const SpliceMotifs& splice_motifs)
    : subgraph(graph, seed_pos, search_left, search_dist + 2), motif_matches(splice_motifs.size())
{
    
#ifdef debug_splice_region
    cerr << "constructing splice region" << endl;
#endif
    
    
    // add a buffer of 2 bases for the dinucleotide itself
    // TODO: feels inelegant to do this here and in the initializer list
    search_dist += 2;
    
    // remember the starting location
    handle_t handle = subgraph.handle_at_order(0);
    seed = pair<handle_t, size_t>(handle, offset(seed_pos));
    
    // extract the subgraph and initialize the DP structure
    vector<pair<handle_t, vector<uint32_t>>> dinuc_states;
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
            if (dinuc_machine.matches(states[j], splice_motifs.oriented_motif(i, search_left))) {
                if ((j == 0 && !search_left) || (j + 1 == states.size() && search_left)) {
                    // we need to cross a node boundary to backtrack
                    subgraph.follow_edges(handle, !search_left, [&](const handle_t& prev) {
                        if (search_left) {
                            if (subgraph.get_base(prev, 0) == splice_motifs.oriented_motif(i, true).front()) {
                                int64_t trav_dist = subgraph.distance_from_start(prev) + subgraph.get_length(prev) - 1;
                                motif_matches[i].emplace_back(prev, 1, trav_dist);
#ifdef debug_splice_region
                                cerr << "record match to motif " << i << " at " << subgraph.order_of(prev) << ", dist " << trav_dist << endl;
#endif
                            }
                        }
                        else {
                            size_t k = subgraph.get_length(prev) - 1;
                            if (subgraph.get_base(prev, k) == splice_motifs.oriented_motif(i, false).front()) {
                                int64_t trav_dist = subgraph.distance_from_start(prev) + k;
                                motif_matches[i].emplace_back(prev, k, trav_dist);
#ifdef debug_splice_region
                                cerr << "record match to motif " << i << " at " << subgraph.order_of(prev) << ", dist " << trav_dist << endl;
#endif
                            }
                        }
                    });
                }
                else {
                    int64_t trav_dist = subgraph.distance_from_start(handle);
                    if (search_left) {
                        trav_dist += states.size() - j - 2;
                    }
                    else {
                        trav_dist += j - 1;
                    }
#ifdef debug_splice_region
                    cerr << "record match to motif " << i << " at " << pos << ", dist " << trav_dist << endl;
#endif
                    motif_matches[i].emplace_back(handle, j - 2 * incr + !search_left, trav_dist);
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

const IncrementalSubgraph& SpliceRegion::get_subgraph() const {
    return subgraph;
}

const pair<handle_t, size_t>& SpliceRegion::get_seed_pos() const {
    return seed;
}

const vector<tuple<handle_t, size_t, int64_t>>& SpliceRegion::candidate_splice_sites(size_t motif_num) const {
    return motif_matches[motif_num];
}

JoinedSpliceGraph::JoinedSpliceGraph(const HandleGraph& parent_graph,
                                     const IncrementalSubgraph& left_subgraph,
                                     handle_t left_splice_node, size_t left_splice_offset,
                                     const IncrementalSubgraph& right_subgraph,
                                     handle_t right_splice_node, size_t right_splice_offset)
    : parent_graph(&parent_graph), left_subgraph(&left_subgraph), right_subgraph(&right_subgraph),
      left_handle_trans(left_subgraph.get_node_count(), -1), right_handle_trans(right_subgraph.get_node_count(), -1),
      left_splice_offset(left_splice_offset), right_splice_offset(right_splice_offset)
{
    // TODO: use the handle translator as scratch instead of these temporary vectors?
    vector<bool> keep_left(left_subgraph.get_node_count(), false);
    vector<bool> keep_right(right_subgraph.get_node_count(), false);
    
    keep_left[left_subgraph.order_of(left_splice_node)] = true;
    vector<handle_t> stack(1, left_splice_node);
    while (!stack.empty()) {
        handle_t here = stack.back();
        stack.pop_back();
        left_subgraph.follow_edges(here, true, [&](const handle_t& prev) {
            if (!keep_left[left_subgraph.order_of(prev)]) {
                keep_left[left_subgraph.order_of(prev)] = true;
                stack.emplace_back(prev);
            }
        });
    }
    
    stack.emplace_back(right_splice_node);
    keep_right[right_subgraph.order_of(right_splice_node)] = true;
    // TODO: repetitive code
    while (!stack.empty()) {
        handle_t here = stack.back();
        stack.pop_back();
        right_subgraph.follow_edges(here, false, [&](const handle_t& prev) {
            if (!keep_right[right_subgraph.order_of(prev)]) {
                keep_right[right_subgraph.order_of(prev)] = true;
                stack.emplace_back(prev);
            }
        });
    }
    
    for (int64_t i = 0; i < left_subgraph.get_node_count(); ++i) {
        if (keep_left[i]) {
            left_handle_trans[i] = handle_idxs.size();
            handle_idxs.push_back(i);
        }
    }
    
    num_left_handles = handle_idxs.size();
    
    // in reverse order
    for (int64_t i = right_subgraph.get_node_count() - 1; i >= 0; --i) {
        if (keep_right[i]) {
            right_handle_trans[i] = handle_idxs.size();
            handle_idxs.push_back(i);
        }
    }
}

pair<size_t, size_t> JoinedSpliceGraph::translate_node_ids(Path& path) const {
    
    pair<size_t, size_t> splice_idxs(numeric_limits<size_t>::max(),
                                     numeric_limits<size_t>::max());
    
    for (size_t i = 0; i < path.mapping_size(); ++i) {
        
        Position* position = path.mutable_mapping(i)->mutable_position();
        
        // record any splice positions
        if (position->node_id() == get_id(left_splice_node())) {
            splice_idxs.first = i;
        }
        else if (position->node_id() == get_id(right_splice_node())) {
            splice_idxs.second = i;
        }
        
        // project down to the parent graph
        size_t j = position->node_id() - 1;
        auto subgraph = j < num_left_handles ? left_subgraph : right_subgraph;
        handle_t underlying = subgraph->get_underlying_handle(subgraph->handle_at_order(handle_idxs[j]));
        if (position->is_reverse()) {
            underlying = parent_graph->flip(underlying);
        }
        // adjust offsets and IDs in the position
        auto interval = underlying_interval(get_handle(position->node_id(), position->is_reverse()));
        position->set_node_id(parent_graph->get_id(underlying));
        position->set_is_reverse(parent_graph->get_is_reverse(underlying));
        position->set_offset(position->offset() + interval.first);
    }
    return splice_idxs;
}

handle_t JoinedSpliceGraph::left_seed_node() const {
    return handlegraph::number_bool_packing::pack(0, false);
}

handle_t JoinedSpliceGraph::right_seed_node() const {
    return handlegraph::number_bool_packing::pack(handle_idxs.size() - 1, false);
}

handle_t JoinedSpliceGraph::left_splice_node() const {
    return handlegraph::number_bool_packing::pack(num_left_handles - 1, false);
}

handle_t JoinedSpliceGraph::right_splice_node() const {
    return handlegraph::number_bool_packing::pack(num_left_handles, false);
}

bool JoinedSpliceGraph::has_node(id_t node_id) const {
    return node_id >= 0 && node_id <= handle_idxs.size();
}

handle_t JoinedSpliceGraph::get_handle(const id_t& node_id, bool is_reverse) const {
    return handlegraph::number_bool_packing::pack(node_id - 1, is_reverse);
}

id_t JoinedSpliceGraph::get_id(const handle_t& handle) const {
    return handlegraph::number_bool_packing::unpack_number(handle) + 1;
}

bool JoinedSpliceGraph::get_is_reverse(const handle_t& handle) const {
    return handlegraph::number_bool_packing::unpack_bit(handle);
}

handle_t JoinedSpliceGraph::flip(const handle_t& handle) const {
    return handlegraph::number_bool_packing::toggle_bit(handle);
}

size_t JoinedSpliceGraph::get_length(const handle_t& handle) const {
    auto interval = underlying_interval(handle);
    return interval.second - interval.first;
}

string JoinedSpliceGraph::get_sequence(const handle_t& handle) const {
    auto interval = underlying_interval(handle);
    size_t i = handlegraph::number_bool_packing::unpack_number(handle);
    const IncrementalSubgraph& subgraph = i < num_left_handles ? *left_subgraph : *right_subgraph;
    handle_t under = subgraph.handle_at_order(handle_idxs[i]);
    if (get_is_reverse(handle)) {
        under = subgraph.flip(under);
    }
    return subgraph.get_subsequence(under, interval.first, interval.second - interval.first);
}

bool JoinedSpliceGraph::follow_edges_impl(const handle_t& handle, bool go_left,
                                          const function<bool(const handle_t&)>& iteratee) const {
    
    bool using_left_edges = go_left != get_is_reverse(handle);
    size_t i = handlegraph::number_bool_packing::unpack_number(handle);
        
    auto traverse_within_subgraph = [&](const IncrementalSubgraph& subgraph,
                                        const vector<int64_t>& handle_trans) {
        // traverse within the subgraph
        handle_t under = subgraph.handle_at_order(handle_idxs[i]);
        if (get_is_reverse(handle)) {
            under = subgraph.flip(under);
        }
        return subgraph.follow_edges(under, go_left, [&](const handle_t& next) {
            // filter to only the handles that are included in the joined graph
            int64_t translated = handle_trans[subgraph.order_of(next)];
            if (translated != -1) {
                return iteratee(handlegraph::number_bool_packing::pack(translated,
                                                                       get_is_reverse(next)));
            }
            else {
                return true;
            }
        });
    };
    
    if (i + 1 < num_left_handles || (i + 1 == num_left_handles && using_left_edges)) {
        // internal to the left subgraph
        return traverse_within_subgraph(*left_subgraph, left_handle_trans);
    }
    else if (i > num_left_handles || (i == num_left_handles && !using_left_edges)) {
        // internal to the right subgraph
        return traverse_within_subgraph(*right_subgraph, right_handle_trans);
    }
    else if (i + 1 == num_left_handles) {
        // rightward across the splice join
        return iteratee(handlegraph::number_bool_packing::pack(num_left_handles,
                                                               get_is_reverse(handle)));
    }
    else {
        // leftward across the splice join
        return iteratee(handlegraph::number_bool_packing::pack(num_left_handles - 1,
                                                               get_is_reverse(handle)));
    }
}

bool JoinedSpliceGraph::for_each_handle_impl(const function<bool(const handle_t&)>& iteratee,
                                             bool parallel) const {
    bool keep_going = true;
    for (size_t i = 0; i < handle_idxs.size() && keep_going; ++i) {
        keep_going = iteratee(handlegraph::number_bool_packing::pack(i, false));
    }
    // not doing parallel, never expect to use it
    return keep_going;
}

size_t JoinedSpliceGraph::get_node_count() const {
    return handle_idxs.size();
}

id_t JoinedSpliceGraph::min_node_id() const {
    return 1;
}

id_t JoinedSpliceGraph::max_node_id() const {
    return handle_idxs.size();
}

char JoinedSpliceGraph::get_base(const handle_t& handle, size_t index) const {
    size_t i = handlegraph::number_bool_packing::unpack_number(handle);
    const IncrementalSubgraph& subgraph = i < num_left_handles ? *left_subgraph : *right_subgraph;
    handle_t under = subgraph.handle_at_order(handle_idxs[i]);
    if (get_is_reverse(handle)) {
        under = subgraph.flip(under);
    }
    auto interval = underlying_interval(handle);
    return subgraph.get_base(under, interval.first + index);
}

string JoinedSpliceGraph::get_subsequence(const handle_t& handle, size_t index, size_t size) const {
    size_t i = handlegraph::number_bool_packing::unpack_number(handle);
    const IncrementalSubgraph& subgraph = i < num_left_handles ? *left_subgraph : *right_subgraph;
    handle_t under = subgraph.handle_at_order(handle_idxs[i]);
    if (get_is_reverse(handle)) {
        under = subgraph.flip(under);
    }
    auto interval = underlying_interval(handle);
    index = min(interval.first + index, interval.second);
    size = min(size, interval.second - index);
    return subgraph.get_subsequence(under, index, size);
}

pair<size_t, size_t> JoinedSpliceGraph::underlying_interval(const handle_t& handle) const {
    size_t i = handlegraph::number_bool_packing::unpack_number(handle);
    const IncrementalSubgraph& subgraph = i < num_left_handles ? *left_subgraph : *right_subgraph;
    handle_t under = subgraph.handle_at_order(handle_idxs[i]);
    size_t begin, end;
    if (i == 0) {
        begin = -subgraph.distance_from_start(under);
    }
    else if (i == num_left_handles) {
        begin = right_splice_offset;
    }
    else {
        begin = 0;
    }
    if (i + 1 == handle_idxs.size()) {
        end = subgraph.get_length(under) + subgraph.distance_from_start(under);
    }
    else if (i + 1 == num_left_handles) {
        end = left_splice_offset;
    }
    else {
        end = subgraph.get_length(under);
    }
    pair<size_t, size_t> return_val(begin, end);
    if (get_is_reverse(handle)) {
        return_val.first = subgraph.get_length(under) - end;
        return_val.second = subgraph.get_length(under) - begin;
    }
    return return_val;
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
