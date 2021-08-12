/**
 * \file splicing.cpp
 *
 * Implements SpliceRegion and some other splicing tools
 *
 */

#include "splicing.hpp"

//#define debug_splice_region
//#define debug_trimming
//#define debug_fusing
//#define debug_from_hit
//#define debug_linker_split

#ifdef debug_splice_region
#include <bitset>
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

bool SpliceMotifs::motif_is_reverse(size_t motif_num) const {
    return motif_num % 2;
}

string SpliceMotifs::unoriented_motif(size_t motif_num, bool left_side) const {
    return left_side ? get<1>(unaltered_data[motif_num / 2]) : get<0>(unaltered_data[motif_num / 2]);
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
        int32_t score = round(log(get<2>(record)) / scorer.log_base);
        data.emplace_back();
        get<0>(data.back()) = get<0>(record);
        // reverse the second string because it's encountered in reverse when going into
        // an intron
        get<1>(data.back()) = string(get<1>(record).rbegin(), get<1>(record).rend());
        // convert frequency to a log likelihood
        get<2>(data.back()) = score;
        
        // now do the reverse complement
        data.emplace_back();
        get<0>(data.back()) = reverse_complement(get<1>(record));
        get<1>(data.back()) = reverse_complement(string(get<0>(record).rbegin(), get<0>(record).rend()));
        get<2>(data.back()) = score;
    }
#ifdef debug_splice_region
    for (const auto& record : data) {
        cerr << "\t" << get<0>(record) << "\t" << get<1>(record) << "\t" << get<2>(record) << endl;
    }
#endif
}

SpliceRegion::SpliceRegion(const pos_t& seed_pos, bool search_left, int64_t search_dist,
                           const HandleGraph& graph,
                           const DinucleotideMachine& dinuc_machine,
                           const SpliceMotifs& splice_motifs)
    : subgraph(graph, seed_pos, search_left, search_dist + 2), motif_matches(splice_motifs.size())
{
    
#ifdef debug_splice_region
    cerr << "constructing splice region starting from seed pos " << seed_pos << " in direction left? " << search_left << ", max dist " << search_dist << endl;
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
                                int64_t trav_dist = subgraph.min_distance_from_start(prev) + subgraph.get_length(prev) - 1;
                                motif_matches[i].emplace_back(prev, 1, trav_dist);
#ifdef debug_splice_region
                                cerr << "record match to motif " << i << " at " << subgraph.order_of(prev) << ", dist " << trav_dist << endl;
#endif
                            }
                        }
                        else {
                            size_t k = subgraph.get_length(prev) - 1;
                            if (subgraph.get_base(prev, k) == splice_motifs.oriented_motif(i, false).front()) {
                                int64_t trav_dist = subgraph.min_distance_from_start(prev) + k;
                                motif_matches[i].emplace_back(prev, k, trav_dist);
#ifdef debug_splice_region
                                cerr << "record match to motif " << i << " at " << subgraph.order_of(prev) << ", dist " << trav_dist << endl;
#endif
                            }
                        }
                    });
                }
                else {
                    int64_t trav_dist = subgraph.min_distance_from_start(handle);
                    if (search_left) {
                        trav_dist += states.size() - j - 2;
                    }
                    else {
                        trav_dist += j - 1;
                    }
#ifdef debug_splice_region
                    cerr << "record match to motif " << i << " at " << (j - 2 * incr + !search_left) << ", dist " << trav_dist << endl;
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
        int64_t prev_dist = subgraph.min_distance_from_start(here);
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

int64_t JoinedSpliceGraph::min_link_length() const {
    handle_t splice_left = left_subgraph->handle_at_order(handle_idxs[num_left_handles - 1]);
    handle_t splice_right = right_subgraph->handle_at_order(handle_idxs[num_left_handles]);
    return (left_subgraph->min_distance_from_start(splice_left)
            + right_subgraph->min_distance_from_start(splice_right)
            + left_splice_offset
            + right_subgraph->get_length(splice_right) - right_splice_offset);
}

int64_t JoinedSpliceGraph::max_link_length() const {
    handle_t splice_left = left_subgraph->handle_at_order(handle_idxs[num_left_handles - 1]);
    handle_t splice_right = right_subgraph->handle_at_order(handle_idxs[num_left_handles]);
    return (left_subgraph->max_distance_from_start(splice_left)
            + right_subgraph->max_distance_from_start(splice_right)
            + left_splice_offset
            + right_subgraph->get_length(splice_right) - right_splice_offset);
    
}

bool JoinedSpliceGraph::has_node(id_t node_id) const {
    return node_id > 0 && node_id <= handle_idxs.size();
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
        begin = -subgraph.min_distance_from_start(under);
    }
    else if (i == num_left_handles) {
        begin = right_splice_offset;
    }
    else {
        begin = 0;
    }
    if (i + 1 == handle_idxs.size()) {
        end = subgraph.get_length(under) + subgraph.min_distance_from_start(under);
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

multipath_alignment_t from_hit(const Alignment& alignment, const HandleGraph& graph,
                               const pos_t& hit_pos, const MaximalExactMatch& mem,
                               const GSSWAligner& scorer) {
    // TODO: mostly copied from multipath alignment graph
    
    multipath_alignment_t multipath_aln;
    transfer_read_metadata(alignment, multipath_aln);
    
    // stack for DFS, each record contains tuples of
    // (read begin, node offset, next node index, next node handles,
    vector<tuple<string::const_iterator, size_t, size_t, vector<handle_t>>> stack;
    stack.emplace_back(mem.begin, offset(hit_pos), 0,
                       vector<handle_t>{graph.get_handle(id(hit_pos), is_rev(hit_pos))});
    while (!stack.empty()) {
        auto& back = stack.back();
        if (get<2>(back) == get<3>(back).size()) {
            stack.pop_back();
            continue;
        }
        
        handle_t trav = get<3>(back)[get<2>(back)];
        get<2>(back)++;
        
#ifdef debug_from_hit
        cerr << "checking node " << graph.get_id(trav) << endl;
#endif
        
        string node_seq = graph.get_sequence(trav);
        size_t node_idx = get<1>(back);
        string::const_iterator read_iter = get<0>(back);
        
        // look for a match along the entire node sequence
        for (; node_idx < node_seq.size() && read_iter != mem.end; node_idx++, read_iter++) {
            if (node_seq[node_idx] != *read_iter) {
#ifdef debug_from_hit
                cerr << "node sequence does not match read" << endl;
#endif
                break;
            }
        }
        
        if (read_iter == mem.end) {
            break;
        }
        
        if (node_idx == node_seq.size()) {
            stack.emplace_back(read_iter, 0, 0, vector<handle_t>());
            graph.follow_edges(trav, false, [&](const handle_t& next) {
                get<3>(stack.back()).emplace_back(next);
            });
        }
    }
    
    subpath_t* subpath = multipath_aln.add_subpath();
    path_t* path = subpath->mutable_path();
    for (size_t i = 0; i < stack.size(); ++i) {
        path_mapping_t* mapping = path->add_mapping();
        if (i == 0 && mem.begin != alignment.sequence().begin()) {
            edit_t* edit = mapping->add_edit();
            edit->set_to_length(mem.begin - alignment.sequence().begin());
            edit->set_from_length(0);
            edit->set_sequence(string(alignment.sequence().begin(), mem.begin));
        }
        
        handle_t handle = get<3>(stack[i])[get<2>(stack[i]) - 1];
        
        position_t* position = mapping->mutable_position();
        position->set_node_id(graph.get_id(handle));
        position->set_is_reverse(graph.get_is_reverse(handle));
        position->set_offset(get<1>(stack[i]));
        
        edit_t* edit = mapping->add_edit();
        int64_t len = min<int64_t>(graph.get_length(handle) - get<1>(stack[i]),
                                   mem.end - get<0>(stack[i]));
        edit->set_to_length(len);
        edit->set_from_length(len);
        
        if (i + 1 == stack.size() && mem.end != alignment.sequence().end()) {
            edit_t* edit = mapping->add_edit();
            edit->set_to_length(alignment.sequence().end() - mem.end);
            edit->set_from_length(0);
            edit->set_sequence(string(mem.end, alignment.sequence().end()));
        }
    }
    subpath->set_score(scorer.score_partial_alignment(alignment, graph, *path,
                                                      alignment.sequence().begin()));
    
    identify_start_subpaths(multipath_aln);
    return multipath_aln;
}

tuple<pos_t, int64_t, int32_t> trimmed_end(const Alignment& aln, int64_t len, bool from_end,
                                           const HandleGraph& graph, const GSSWAligner& aligner) {
    
#ifdef debug_trimming
    cerr << "trimming alignment " << pb2json(aln) << " by " << len << ", from end? " << from_end << endl;
#endif
    
    const Path& path = aln.path();
    path_t dummy_path;
    
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
                
                from_proto_mapping(path.mapping(i), *dummy_path.add_mapping());
                --i;
            }
            if (i < 0) {
#ifdef debug_trimming
                cerr << "walked entire path" << endl;
#endif
                get<0>(return_val) = initial_position(path);
                get<1>(return_val) = path_to_length(path);
                dummy_path.clear_mapping();
                from_proto_path(path, dummy_path);
                copied_full_path = true;
            }
            else {
                const Mapping& mapping = path.mapping(i);
                int64_t j = mapping.edit_size() - 1;
                int64_t from_length = 0;
                path_mapping_t* dummy_mapping = nullptr;
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
                    from_proto_edit(mapping.edit(j), *dummy_mapping->add_edit());
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
                    auto* dummy_edit = dummy_mapping->add_edit();
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
                if (dummy_mapping) {
                    from_proto_position(mapping.position(), *dummy_mapping->mutable_position());
                }
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
                
                from_proto_mapping(path.mapping(i), *dummy_path.add_mapping());
                ++i;
            }
            if (i == path.mapping_size()) {
#ifdef debug_trimming
                cerr << "walked entire path" << endl;
#endif
                get<0>(return_val) = final_position(path);
                get<1>(return_val) = path_to_length(path);
                dummy_path.clear_mapping();
                from_proto_path(path, dummy_path);
                copied_full_path = true;
            }
            else {
                const Mapping& mapping = path.mapping(i);
                int64_t j = 0;
                int64_t from_length = 0;
                path_mapping_t* dummy_mapping = nullptr;
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
                        from_proto_position(mapping.position(), *dummy_mapping->mutable_position());
                    }
                    from_proto_edit(mapping.edit(j), *dummy_mapping->add_edit());
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
                        from_proto_position(mapping.position(), *dummy_mapping->mutable_position());
                    }
                    auto* dummy_edit = dummy_mapping->add_edit();
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
            auto* mapping = dummy_path.mutable_mapping(0);
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
    
    get<2>(return_val) = aligner.score_partial_alignment(aln, graph, dummy_path, begin);
    
    return return_val;
}

// TODO: this implementation ended up requiring a lot of duplicated code, i could probably clean it up
bool trim_path(path_t* path, bool from_left, int64_t mapping_idx, int64_t edit_idx, int64_t base_idx) {
    
    bool do_trim = ((from_left && (mapping_idx != 0 || edit_idx != 0 || base_idx != 0)) ||
                    (!from_left && mapping_idx != path->mapping_size()));
    
    if (edit_idx == 0 && base_idx == 0) {
        // position is past-the-last on a mapping
        if (from_left) {
            auto mappings = path->mutable_mapping();
            mappings->erase(mappings->begin(), mappings->begin() + mapping_idx);
        }
        else {
            path->mutable_mapping()->resize(mapping_idx);
        }
    }
    else {
        // position is inside a mapping
        auto mapping = path->mutable_mapping(mapping_idx);
        if (base_idx == 0) {
            // position is past-the-last on an edit
            if (from_left) {
                int64_t from_length_removed = 0;
                for (int64_t i = 0; i < edit_idx; ++i) {
                    from_length_removed += mapping->edit(i).from_length();
                }
                auto edits = mapping->mutable_edit();
                edits->erase(edits->begin(), edits->begin() + edit_idx);
                mapping->mutable_position()->set_offset(mapping->position().offset()
                                                        + from_length_removed);
                auto mappings = path->mutable_mapping();
                mappings->erase(mappings->begin(), mappings->begin() + mapping_idx);
            }
            else {
                mapping->mutable_edit()->resize(edit_idx);
                path->mutable_mapping()->resize(mapping_idx + 1);
            }
        }
        else {
            // position is inside an edit
            auto edit = mapping->mutable_edit(edit_idx);
            if (from_left) {
                int64_t from_length_removed = 0;
                if (base_idx > 0) {
                    if (edit->from_length() > 0) {
                        from_length_removed += base_idx;
                    }
                    edit->set_from_length(max<int64_t>(edit->from_length() - base_idx, 0));
                    edit->set_to_length(max<int64_t>(edit->to_length() - base_idx, 0));
                    if (!edit->sequence().empty()) {
                        edit->set_sequence(edit->sequence().substr(base_idx, edit->to_length()));
                    }
                    if (edit->from_length() == 0 && edit->to_length() == 0) {
                        ++edit_idx;
                    }
                }
                for (int64_t i = 0; i < edit_idx; ++i) {
                    from_length_removed += mapping->edit(i).from_length();
                }
                auto edits = mapping->mutable_edit();
                edits->erase(edits->begin(), edits->begin() + edit_idx);
                mapping->mutable_position()->set_offset(mapping->position().offset()
                                                        + from_length_removed);
                auto mappings = path->mutable_mapping();
                mappings->erase(mappings->begin(), mappings->begin() + mapping_idx);
            }
            else {
                if (base_idx < max(edit->from_length(), edit->to_length())) {
                    edit->set_from_length(min<int64_t>(edit->from_length(), base_idx));
                    edit->set_to_length(min<int64_t>(edit->to_length(), base_idx));
                    if (!edit->sequence().empty()) {
                        edit->set_sequence(edit->sequence().substr(0, base_idx));
                    }
                    if (edit->from_length() == 0 && edit->to_length() == 0) {
                        --edit_idx;
                    }
                }
                mapping->mutable_edit()->resize(edit_idx + 1);
                path->mutable_mapping()->resize(mapping_idx + 1);
            }
        }
    }
    return do_trim;
}

pair<pair<path_t, int32_t>, pair<path_t, int32_t>> split_splice_segment(const Alignment& splice_segment,
                                                                        const tuple<int64_t, int64_t, int64_t>& left_trace,
                                                                        const tuple<int64_t, int64_t, int64_t>& right_trace,
                                                                        int64_t splice_junction_idx,
                                                                        const GSSWAligner& scorer,
                                                                        const HandleGraph& graph) {
    
#ifdef debug_linker_split
    cerr << "splitting splice segment " << pb2json(splice_segment) << endl;
    cerr << "split is at index " << splice_junction_idx << endl;
    cerr << "left trace: " << get<0>(left_trace) << " " << get<1>(left_trace) << " " << get<2>(left_trace) << endl;
    cerr << "right trace: " << get<0>(right_trace) << " " << get<1>(right_trace) << " " << get<2>(right_trace) << endl;
#endif
    
    // TODO: make the scoring robust to hanging indels at the trace positions
    
    pair<pair<path_t, int32_t>, pair<path_t, int32_t>> return_val;
    auto& left_path = return_val.first.first;
    
    // walk the part of the splice segment before the trace on the left side
    size_t left_to_length = 0;
    size_t left_leading_to_length = 0;
    for (int64_t i = 0; i < get<0>(left_trace); ++i) {
        left_leading_to_length += mapping_to_length(splice_segment.path().mapping(i));
    }
    // special logic to handle the mapping with the traced location
    if (get<0>(left_trace) < splice_junction_idx) {
        path_mapping_t* post_trace_mapping = nullptr;
        size_t trace_leading_from_length = 0;
        const auto& trace_mapping = splice_segment.path().mapping(get<0>(left_trace));
        for (int64_t j = 0; j < get<1>(left_trace); ++j) {
            left_leading_to_length += trace_mapping.edit(j).to_length();
            trace_leading_from_length += trace_mapping.edit(j).from_length();
        }
        if (get<1>(left_trace) < trace_mapping.edit_size()) {
            const auto& trace_edit = trace_mapping.edit(get<1>(left_trace));
            if (trace_edit.to_length() != 0) {
                left_leading_to_length += get<2>(left_trace);
            }
            if (trace_edit.from_length() != 0) {
                trace_leading_from_length += get<2>(left_trace);
            }
            if (get<2>(left_trace) < max(trace_edit.from_length(), trace_edit.to_length())) {
                post_trace_mapping = left_path.add_mapping();
                auto edit = post_trace_mapping->add_edit();
                edit->set_from_length(max<int64_t>(0, trace_edit.from_length() - get<2>(left_trace)));
                edit->set_to_length(max<int64_t>(0, trace_edit.to_length() - get<2>(left_trace)));
                if (!trace_edit.sequence().empty()) {
                    edit->set_sequence(trace_edit.sequence().substr(get<2>(left_trace), string::npos));
                }
                left_to_length += edit->to_length();
            }
        }
        for (int64_t j = get<1>(left_trace) + 1; j < trace_mapping.edit_size(); ++j) {
            if (!post_trace_mapping) {
                post_trace_mapping = left_path.add_mapping();
            }
            from_proto_edit(trace_mapping.edit(j), *post_trace_mapping->add_edit());
        }
        if (post_trace_mapping) {
            const auto& trace_pos = trace_mapping.position();
            auto pos = post_trace_mapping->mutable_position();
            pos->set_node_id(trace_pos.node_id());
            pos->set_is_reverse(trace_pos.is_reverse());
            pos->set_offset(trace_pos.offset() + trace_leading_from_length);
        }
    }
    
    // pull out the left half of the alignment
    for (int64_t i = get<0>(left_trace) + 1; i < splice_junction_idx; ++i) {
        const auto& mapping = splice_segment.path().mapping(i);
        if (mapping_from_length(mapping) == 0 && mapping_to_length(mapping) == 0) {
            continue;
        }
        auto left_mapping = left_path.add_mapping();
        from_proto_position(mapping.position(), *left_mapping->mutable_position());
        for (size_t j = 0; j < mapping.edit_size(); ++j) {
            const auto& edit = mapping.edit(j);
            if (edit.from_length() != 0 || edit.to_length() != 0) {
                from_proto_edit(edit, *left_mapping->add_edit());
                left_to_length += edit.to_length();
            }
        }
    }
    
    // and then pull out the right half of the alignment up to the trace point
    auto& right_path = return_val.second.first;
    for (size_t i = splice_junction_idx; i < get<0>(right_trace); ++i) {
        const auto& mapping = splice_segment.path().mapping(i);
        if (mapping_from_length(mapping) == 0 && mapping_to_length(mapping) == 0) {
            continue;
        }
        auto right_mapping = right_path.add_mapping();
        from_proto_position(mapping.position(), *right_mapping->mutable_position());
        for (size_t j = 0; j < mapping.edit_size(); ++j) {
            const auto& edit = mapping.edit(j);
            if (edit.from_length() != 0 || edit.to_length() != 0) {
                from_proto_edit(edit, *right_mapping->add_edit());
            }
        }
    }
    
    // and special logic for the right trace mapping
    if (get<0>(right_trace) >= splice_junction_idx &&
        get<0>(right_trace) < splice_segment.path().mapping_size() &&
        (get<1>(right_trace) != 0 || get<2>(right_trace) != 0)) {
        auto pre_trace_mapping = right_path.add_mapping();
        const auto& trace_mapping = splice_segment.path().mapping(get<0>(right_trace));
        from_proto_position(trace_mapping.position(), *pre_trace_mapping->mutable_position());
        for (int64_t j = 0; j < get<1>(right_trace); ++j) {
            from_proto_edit(trace_mapping.edit(j), *pre_trace_mapping->add_edit());
        }
        if (get<1>(right_trace) < trace_mapping.edit_size() && get<2>(right_trace) != 0) {
            const auto& trace_edit = trace_mapping.edit(get<1>(right_trace));
            auto edit = pre_trace_mapping->add_edit();
            if (trace_edit.from_length() != 0) {
                edit->set_from_length(get<2>(right_trace));
            }
            if (trace_edit.to_length() != 0) {
                edit->set_to_length(get<2>(right_trace));
            }
            if (!trace_edit.sequence().empty()) {
                edit->set_sequence(trace_edit.sequence().substr(0, get<2>(right_trace)));
            }
        }
    }
    
    // score the two halves (but don't take the full length bonus, since this isn't actually
    // the end of the full read)
    return_val.first.second = scorer.score_partial_alignment(splice_segment, graph, left_path,
                                                             splice_segment.sequence().begin() + left_leading_to_length, true);
    return_val.second.second = scorer.score_partial_alignment(splice_segment, graph, right_path,
                                                              splice_segment.sequence().begin() + left_to_length + left_leading_to_length,
                                                              true);
    
#ifdef debug_linker_split
    cerr << "left partial score " << return_val.first.second << ", partial path " << pb2json(return_val.first.first) << endl;
    cerr << "right partial score " << return_val.second.second << ", partial path " << pb2json(return_val.second.first) << endl;
#endif
    
    // deletions can span the splice junction, in which case they will have been scored incorrectly
    // by taking the gap open penalty twice
    if (return_val.first.first.mapping_size() != 0 && return_val.second.first.mapping_size() != 0) {
        
        int64_t left_del_size = 0, right_del_size = 0;
        
        // measure left deletion at the end
        bool in_deletion = true;
        for (int64_t i = left_path.mapping_size() - 1; i >= 0 && in_deletion; --i) {
            const auto& mapping = left_path.mapping(i);
            for (int64_t j = mapping.edit_size() - 1; j >= 0 && in_deletion; --j) {
                const auto& edit = mapping.edit(j);
                if (edit.to_length() == 0) {
                    left_del_size += edit.from_length();
                }
                else {
                    in_deletion = false;
                }
            }
        }
        
        // measure right deletion at the beginning
        in_deletion = true;
        for (int64_t i = 0; i < right_path.mapping_size() && in_deletion; ++i) {
            const auto& mapping = right_path.mapping(i);
            for (int64_t j = 0; j < mapping.edit_size() && in_deletion; ++j) {
                const auto& edit = mapping.edit(j);
                if (edit.to_length() == 0) {
                    right_del_size += edit.from_length();
                }
                else {
                    in_deletion = false;
                }
            }
        }
        
#ifdef debug_linker_split
        cerr << "spanning split, left deletion size: " << left_del_size << ", right deletion size " << right_del_size << endl;
#endif
        
        if (left_del_size != 0 && right_del_size != 0) {
            // split the total gap score between the two (can break dynamic programmability
            // a little bit, but i think it's worth it to have a good alignment across the
            // splice junction)
            int32_t total_gap_score = scorer.score_gap(left_del_size + right_del_size);
            return_val.first.second += (total_gap_score / 2 - scorer.score_gap(left_del_size));
            return_val.second.second += (total_gap_score  - total_gap_score / 2 - scorer.score_gap(right_del_size));
            
#ifdef debug_linker_split
            cerr << "re-divided deletion score, now left score " << return_val.first.second << ", right score " << return_val.second.second << endl;
#endif
        }
    }
    
    return return_val;
}

multipath_alignment_t&& fuse_spliced_alignments(const Alignment& alignment,
                                                multipath_alignment_t&& left_mp_aln, multipath_alignment_t&& right_mp_aln,
                                                int64_t left_bridge_point, const Alignment& splice_segment,
                                                int64_t splice_junction_idx, int32_t splice_score, const GSSWAligner& scorer,
                                                const HandleGraph& graph) {
    
#ifdef debug_fusing
    cerr << "fusing spliced alignments" << endl;
    cerr << "left alignment" << endl;
    cerr << debug_string(left_mp_aln) << endl;
    cerr << "right alignment" << endl;
    cerr << debug_string(right_mp_aln) << endl;
    cerr << "linker" << endl;
    cerr << pb2json(splice_segment) << endl;
#endif
    
    // TODO: allow multiple splices to happen on the same multipath alignment?
    
    pos_t pos_left = initial_position(splice_segment.path());
    pos_t pos_right = final_position(splice_segment.path());
    
    int64_t right_bridge_point = left_bridge_point + path_to_length(splice_segment.path());
    
    auto left_locations = search_multipath_alignment(left_mp_aln, pos_left, left_bridge_point);
    auto right_locations = search_multipath_alignment(right_mp_aln, pos_right, right_bridge_point);
    
#ifdef debug_fusing
    cerr << "left splice locations:" << endl;
    for (auto loc : left_locations) {
        cerr << "\t" << get<0>(loc) << " " << get<1>(loc) << " " << get<2>(loc) << " " << get<3>(loc) << endl;
    }
    cerr << "right splice locations:" << endl;
    for (auto loc : right_locations) {
        cerr << "\t" << get<0>(loc) << " " << get<1>(loc) << " " << get<2>(loc) << " " << get<3>(loc) << endl;
    }
#endif
    
    if (left_locations.empty() || right_locations.empty()) {
        cerr << "error: splice segment could not be located on multipath alignment of read " << alignment.name() << endl;
        exit(1);
    }
    
    // trace the splice segment from first location on the list
    tuple<int64_t, int64_t, int64_t> left_path_trace, right_path_trace;
    vector<tuple<int64_t, int64_t, int64_t, int64_t>> left_mp_aln_trace, right_mp_aln_trace;
    tie(left_path_trace, left_mp_aln_trace) = trace_path(left_mp_aln, splice_segment.path(), get<0>(left_locations.front()),
                                                         get<1>(left_locations.front()), get<2>(left_locations.front()),
                                                         get<3>(left_locations.front()), false, splice_junction_idx);
    tie(right_path_trace, right_mp_aln_trace) = trace_path(right_mp_aln, splice_segment.path(), get<0>(right_locations.front()),
                                                           get<1>(right_locations.front()), get<2>(right_locations.front()),
                                                           get<3>(right_locations.front()), true, splice_junction_idx);
    
    // trace the splice segment on subsequent locations on the left
    for (size_t i = 1; i < left_locations.size(); ++i) {
        int64_t j, k, l, m;
        tie(j, k, l, m) = left_locations[i];
        auto trace = trace_path(left_mp_aln, splice_segment.path(), j, k, l, m, false, splice_junction_idx);
        if (trace.first > left_path_trace) {
            left_path_trace = trace.first;
            left_mp_aln_trace = move(trace.second);
        }
        else if (trace.first == left_path_trace) {
            for (auto& mp_aln_trace : trace.second) {
                left_mp_aln_trace.push_back(mp_aln_trace);
            }
        }
    }
    
    // trace the splice segment on subsequent locations on the right
    for (size_t i = 1; i < right_locations.size(); ++i) {
        int64_t j, k, l, m;
        tie(j, k, l, m) = right_locations[i];
        auto trace = trace_path(right_mp_aln, splice_segment.path(), j, k, l, m, true, splice_junction_idx);
        if (trace.first < right_path_trace) {
            right_path_trace = trace.first;
            right_mp_aln_trace = move(trace.second);
        }
        else if (trace.first == right_path_trace) {
            for (auto& mp_aln_trace : trace.second) {
                right_mp_aln_trace.push_back(mp_aln_trace);
            }
        }
    }
    
    // now also walk back the multipath locations to match
    vector<vector<int64_t>> left_rev_nexts;
    for (size_t i = 0; i < left_mp_aln_trace.size(); ++i) {
        int64_t s, m, e, b;
        tie(s, m, e, b) = left_mp_aln_trace[i];
    }
    
    
    // ensure that the connection locations are unique and in increasing order by subpath index
    sort(left_mp_aln_trace.begin(), left_mp_aln_trace.end());
    sort(right_mp_aln_trace.begin(), right_mp_aln_trace.end());
    auto new_left_end = unique(left_mp_aln_trace.begin(), left_mp_aln_trace.end());
    left_mp_aln_trace.erase(new_left_end, left_mp_aln_trace.end());
    auto new_right_end = unique(right_mp_aln_trace.begin(), right_mp_aln_trace.end());
    right_mp_aln_trace.erase(new_right_end, right_mp_aln_trace.end());
        
#ifdef debug_fusing
    cerr << "left traced path location:" << endl;
    cerr << "\t" << get<0>(left_path_trace) << " " << get<1>(left_path_trace) << " " << get<2>(left_path_trace) << endl;
    cerr << "left traced multipath locations:" << endl;
    for (auto loc : left_mp_aln_trace) {
        cerr << "\t" << get<0>(loc) << " " << get<1>(loc) << " " << get<2>(loc) << " " << get<3>(loc) << endl;
    }
    cerr << "right traced path location:" << endl;
    cerr << "\t" << get<0>(right_path_trace) << " " << get<1>(right_path_trace) << " " << get<2>(right_path_trace) << endl;
    cerr << "right traced multipath locations:" << endl;
    for (auto loc : right_mp_aln_trace) {
        cerr << "\t" << get<0>(loc) << " " << get<1>(loc) << " " << get<2>(loc) << " " << get<3>(loc) << endl;
    }
#endif
    
    vector<bool> to_keep_left(left_mp_aln.subpath_size(), false);
    vector<bool> is_bridge_left(left_mp_aln.subpath_size(), false);
    for (const auto& pos : left_mp_aln_trace) {
        is_bridge_left[get<0>(pos)] = true;
    }
    for (int64_t i = left_mp_aln.subpath_size() - 1; i >= 0; --i) {
        if (is_bridge_left[i]) {
            to_keep_left[i] = true;
        }
        else {
            for (auto j : left_mp_aln.subpath(i).next()) {
                to_keep_left[i] = to_keep_left[j] || to_keep_left[i];
            }
            for (const auto& connection : left_mp_aln.subpath(i).connection()) {
                to_keep_left[i] = to_keep_left[connection.next()] || to_keep_left[i];
            }
        }
    }
    
#ifdef debug_fusing
    cerr << "deciding what to remove on left:" << endl;
    for (size_t i = 0; i < to_keep_left.size(); ++i) {
        cerr << "\t" << i << ": keep? " << to_keep_left[i] << ", bridge? " << is_bridge_left[i] << endl;
    }
#endif
    
    vector<int64_t> left_removed_so_far(left_mp_aln.subpath_size() + 1, 0);
    int64_t left_loc_idx = 0;
    vector<int64_t> left_to_length(left_mp_aln.subpath_size(), 0);
    for (int64_t i = 0; i < left_mp_aln.subpath_size(); ++i) {
        // keep track of the read interval
        int64_t to_len = path_to_length(left_mp_aln.subpath(i).path());
        for (auto n : left_mp_aln.subpath(i).next()) {
            left_to_length[n] = left_to_length[i] + to_len;
        }
        for (const auto& connection : left_mp_aln.subpath(i).connection()) {
            left_to_length[connection.next()] = left_to_length[i] + to_len;
        }
        
        if (!to_keep_left[i]) {
            left_removed_so_far[i + 1] = left_removed_so_far[i] + 1;
            continue;
        }
        // TODO: what if there are multiple locations with hits on the same subpath...
        if (left_loc_idx < left_mp_aln_trace.size() && i == get<0>(left_mp_aln_trace[left_loc_idx])) {
            int64_t s, m, e, b;
            tie(s, m, e, b) = left_mp_aln_trace[left_loc_idx];
            ++left_loc_idx;
            
#ifdef debug_fusing
            cerr << "trimming subpath " << s << " at location " << m << " " << e << " " << b << endl;
#endif
            
            auto path = left_mp_aln.mutable_subpath(i)->mutable_path();
            bool trimmed = trim_path(path, false, m, e, b);
            if (trimmed) {
                int32_t new_score = scorer.score_partial_alignment(alignment, graph, *path,
                                                                   alignment.sequence().begin() + left_to_length[i]);
                left_mp_aln.mutable_subpath(i)->set_score(new_score);
            }
            left_mp_aln.mutable_subpath(i)->mutable_next()->clear();
            // TODO: i don't like doing this... need to revise for multiple cuts
            left_mp_aln.mutable_subpath(i)->mutable_connection()->clear();
        }
        
        if (left_removed_so_far[i]) {
            *left_mp_aln.mutable_subpath(i - left_removed_so_far[i]) = move(*left_mp_aln.mutable_subpath(i));
        }
        left_removed_so_far[i + 1] = left_removed_so_far[i];
    }
    if (left_removed_so_far.back() != 0) {
        left_mp_aln.mutable_subpath()->resize(left_mp_aln.subpath_size() - left_removed_so_far.back());
        for (size_t i = 0; i < left_mp_aln.subpath_size(); ++i)  {
            auto subpath = left_mp_aln.mutable_subpath(i);
            size_t nexts_removed_so_far = 0;
            for (size_t j = 0; j < subpath->next_size(); ++j) {
                if (to_keep_left[subpath->next(j)]) {
                    subpath->set_next(j - nexts_removed_so_far,
                                      subpath->next(j) - left_removed_so_far[subpath->next(j)]);
                }
                else {
                    ++nexts_removed_so_far;
                }
            }
            if (nexts_removed_so_far != 0) {
                subpath->mutable_next()->resize(subpath->next_size() - nexts_removed_so_far);
            }
            size_t connections_removed_so_far = 0;
            for (size_t j = 0; j < subpath->connection_size(); ++j) {
                auto connection = subpath->mutable_connection(j);
                if (to_keep_left[connection->next()]) {
                    connection->set_next(connection->next() - left_removed_so_far[connection->next()]);
                    if (connections_removed_so_far) {
                        *subpath->mutable_connection(j - connections_removed_so_far) = *connection;
                    }
                }
                else {
                    ++connections_removed_so_far;
                }
            }
            if (connections_removed_so_far != 0) {
                subpath->mutable_connection()->resize(subpath->connection_size() - connections_removed_so_far);
            }
        }
    }
    
    
#ifdef debug_fusing
    cerr << "after processing left:" << endl;
    cerr << debug_string(left_mp_aln) << endl;
#endif
    
    size_t left_subpaths_end = left_mp_aln.subpath_size();
    
    auto splice_segment_halves = split_splice_segment(splice_segment, left_path_trace, right_path_trace,
                                                      splice_junction_idx, scorer, graph);
    
#ifdef debug_fusing
    cerr << "split the linker:" << endl;
    cerr << "left score " << splice_segment_halves.first.second << ", path " << pb2json(splice_segment_halves.first.first) << endl;
    cerr << "right score " << splice_segment_halves.second.second << ", path " << pb2json(splice_segment_halves.second.first) << endl;
#endif
    
    if (splice_segment_halves.first.first.mapping_size() != 0) {
        for (const auto& left_loc : left_mp_aln_trace) {
            auto i = get<0>(left_loc) - left_removed_so_far[get<0>(left_loc)];
            if (i < left_mp_aln.subpath_size()) {
                left_mp_aln.mutable_subpath(i)->add_next(left_mp_aln.subpath_size());
            }
        }
        
        auto subpath = left_mp_aln.add_subpath();
        subpath->set_score(splice_segment_halves.first.second);
        *subpath->mutable_path() = move(splice_segment_halves.first.first);
    }
    
    if (splice_segment_halves.second.first.mapping_size() != 0) {
        if (splice_segment_halves.first.first.mapping_size() == 0) {
            // we skipped the left side of the splice, so connect to the left splice points
            for (const auto& left_loc : left_mp_aln_trace) {
                auto i = get<0>(left_loc) - left_removed_so_far[get<0>(left_loc)];
                auto connection = left_mp_aln.mutable_subpath(i)->add_connection();
                connection->set_next(left_mp_aln.subpath_size());
                connection->set_score(splice_score);
            }
        }
        else {
            // the left side also exists as a subpath, connect to it
            auto connection = left_mp_aln.mutable_subpath(left_mp_aln.subpath_size() - 1)->add_connection();
            connection->set_next(left_mp_aln.subpath_size());
            connection->set_score(splice_score);
        }
        
        auto subpath = left_mp_aln.add_subpath();
        subpath->set_score(splice_segment_halves.second.second);
        *subpath->mutable_path() = move(splice_segment_halves.second.first);
    }
    
#ifdef debug_fusing
    cerr << "after processing linker:" << endl;
    cerr << debug_string(left_mp_aln) << endl;
#endif
    
    size_t right_subpaths_begin = left_mp_aln.subpath_size();
    
    // figure out the read position of the subpaths before we start getting rid of edges
    vector<int64_t> right_to_length(right_mp_aln.subpath_size(), 0);
    for (size_t i = 0; i < right_mp_aln.subpath_size(); ++i) {
        // keep track of the read interval
        int64_t to_len = path_to_length(right_mp_aln.subpath(i).path());
        for (auto n : right_mp_aln.subpath(i).next()) {
            right_to_length[n] = right_to_length[i] + to_len;
        }
        for (const auto& connection : right_mp_aln.subpath(i).connection()) {
            right_to_length[connection.next()] = right_to_length[i] + to_len;
        }
    }
    
    vector<bool> to_keep_right(right_mp_aln.subpath_size(), false);
    vector<bool> is_bridge_right(right_mp_aln.subpath_size(), false);
    for (const auto& pos : right_mp_aln_trace) {
        is_bridge_right[get<0>(pos)] = true;
    }
    for (int64_t i = 0; i < right_mp_aln.subpath_size(); ++i) {
        to_keep_right[i] = to_keep_right[i] || is_bridge_right[i];
        for (auto j : right_mp_aln.subpath(i).next()) {
            to_keep_right[j] = to_keep_right[j] || to_keep_right[i];
        }
        for (const auto& connection : right_mp_aln.subpath(i).connection()) {
            to_keep_right[connection.next()] = to_keep_right[connection.next()] || to_keep_right[i];
        }
    }
    
#ifdef debug_fusing
    cerr << "deciding what to remove on right:" << endl;
    for (size_t i = 0; i < to_keep_right.size(); ++i) {
        cerr << "\t" << i << ": keep? " << to_keep_right[i] << ", bridge? " << is_bridge_right[i] << endl;
    }
#endif
    
    // transfer the subpaths from the right multipath alignment onto the left one
    
    vector<int64_t> right_removed_so_far(right_mp_aln.subpath_size() + 1, 0);
    int64_t right_loc_idx = 0;
    for (int64_t i = 0; i < right_mp_aln.subpath_size(); ++i) {
        if (!to_keep_right[i]) {
            right_removed_so_far[i + 1] = right_removed_so_far[i] + 1;
            continue;
        }
        if (right_loc_idx < right_mp_aln_trace.size() && i == get<0>(right_mp_aln_trace[right_loc_idx])) {
            // this is where the splice alignment began
            int64_t s, m, e, b;
            tie(s, m, e, b) = right_mp_aln_trace[right_loc_idx];
            ++right_loc_idx;
            
            
#ifdef debug_fusing
            cerr << "trimming subpath " << s << " at location " << m << " " << e << " " << b << endl;
#endif
            
            auto path = right_mp_aln.mutable_subpath(i)->mutable_path();
            int64_t to_len = path_to_length(*path);
            bool trimmed = trim_path(path, true, m, e, b);
            if (trimmed) {
                int64_t new_to_len = path_to_length(*path);
                int32_t new_score = scorer.score_partial_alignment(alignment, graph, *path,
                                                                   (alignment.sequence().begin() + right_to_length[i]
                                                                    + to_len - new_to_len));
                right_mp_aln.mutable_subpath(i)->set_score(new_score);
            }
            // add the edges into the right side
            if (right_subpaths_begin == left_subpaths_end) {
                // both of the splice segments were empty, connect directly to the left splice point
                for (const auto& left_loc : left_mp_aln_trace) {
                    auto i = get<0>(left_loc) - left_removed_so_far[get<0>(left_loc)];
                    if (i < left_subpaths_end) {
                        auto connection = left_mp_aln.mutable_subpath(i)->add_connection();
                        connection->set_next(left_mp_aln.subpath_size());
                        connection->set_score(splice_score);
                    }
                }
            }
            else if (splice_segment_halves.second.first.mapping_size() == 0) {
                // the splice segment on the right side is empty, make connection the the left side
                auto connection = left_mp_aln.mutable_subpath(right_subpaths_begin - 1)->add_connection();
                connection->set_next(left_mp_aln.subpath_size());
                connection->set_score(splice_score);
            }
            else {
                // add the edge from the rightmost splice segment
                left_mp_aln.mutable_subpath(right_subpaths_begin - 1)->add_next(left_mp_aln.subpath_size());
            }
        }
        
        *left_mp_aln.add_subpath() = move(*right_mp_aln.mutable_subpath(i));
        right_removed_so_far[i + 1] = right_removed_so_far[i];
    }
    // fix up the edges on the transferred subpaths from the right alignment
    for (size_t i = right_subpaths_begin; i < left_mp_aln.subpath_size(); ++i)  {
        
        auto subpath = left_mp_aln.mutable_subpath(i);
        size_t nexts_removed_so_far = 0;
        for (size_t j = 0; j < subpath->next_size(); ++j) {
            if (to_keep_right[subpath->next(j)]) {
                subpath->set_next(j - nexts_removed_so_far,
                                  subpath->next(j) - right_removed_so_far[subpath->next(j)] + right_subpaths_begin);
            }
            else {
                ++nexts_removed_so_far;
            }
        }
        if (nexts_removed_so_far != 0) {
            subpath->mutable_next()->resize(subpath->next_size() - nexts_removed_so_far);
        }
        size_t connections_removed_so_far = 0;
        for (size_t j = 0; j < subpath->connection_size(); ++j) {
            auto connection = subpath->mutable_connection(j);
            if (to_keep_right[connection->next()]) {
                connection->set_next(connection->next() - right_removed_so_far[connection->next()] + right_subpaths_begin);
                if (connections_removed_so_far) {
                    *subpath->mutable_connection(j - connections_removed_so_far) = *connection;
                }
            }
            else {
                ++connections_removed_so_far;
            }
        }
        if (connections_removed_so_far != 0) {
            subpath->mutable_connection()->resize(subpath->connection_size() - connections_removed_so_far);
        }
    }
    
#ifdef debug_fusing
    cerr << "fused alignment:" << endl;
    cerr << debug_string(left_mp_aln) << endl;
#endif
    
    // starts can change pretty drastically, so just clear and reidentify
    identify_start_subpaths(left_mp_aln);
    
    // and remove any empty bits
    remove_empty_alignment_sections(left_mp_aln);
    
#ifdef debug_fusing
    cerr << "final product after removing empty subpaths:" << endl;
    cerr << debug_string(left_mp_aln) << endl;
#endif
    
    // pass the left (where we collected everything) out without copying
    return move(left_mp_aln);
}


}
