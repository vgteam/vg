/**
 * \file chain_items.cpp
 * Non-template function implementations for chaining pieces of a read-to-graph alignment.
 */


#include "chain_items.hpp"
#include "crash.hpp"

#include <handlegraph/algorithms/dijkstra.hpp>
#include <structures/immutable_list.hpp>
#include <structures/min_max_heap.hpp>

#include <algorithm>

//#define debug_chaining
//#define debug_transition
//#define debug_dp

namespace vg {
namespace algorithms {

using namespace std;

ostream& operator<<(ostream& out, const Anchor& anchor) {
    // TODO: Just friend class to get these?
    size_t margin_left = anchor.read_start() - anchor.read_exclusion_start();
    size_t margin_right = anchor.read_exclusion_end() - anchor.read_end();
    if (margin_left) {
        out << "(" << margin_left << ")";
    }
    out << "{R:" << anchor.read_start() << "=G:" << anchor.graph_start() << "(+" << anchor.start_hint_offset() 
        << ")-"  << anchor.graph_end() << "(-" << anchor.end_hint_offset() << ")*" << anchor.length() << "}";
    if (margin_right) {
        out << "(" << margin_right << ")";
    }
    return out;
}

ostream& operator<<(ostream& out, const TracedScore& value) {
    if (value.source == TracedScore::nowhere()) {
        return out << value.score << " from nowhere";
    }
    return out << value.score << " from #" << value.source;
}


void TracedScore::max_in(const vector<TracedScore>& options, size_t option_number) {
    auto& option = options[option_number];
    if (option.score > this->score || this->source == nowhere()) {
        // This is the new winner.
        this->score = option.score;
        this->source = option_number;
        this->paths = option.paths;
        this->rec_num = option.rec_num;
    }
}

TracedScore TracedScore::score_from(const vector<TracedScore>& options, size_t option_number) {
    TracedScore got = options[option_number];
    got.source = option_number;
    got.paths = options[option_number].paths;
    got.rec_num = options[option_number].rec_num;
    return got;
}

TracedScore TracedScore::add_points(int adjustment) const {
    return {this->score + adjustment, this->source, this->paths, this->rec_num};
}

TracedScore TracedScore::set_shared_paths(const std::pair<size_t,size_t>& new_paths) const {
    size_t updated_paths;
    size_t rec_num = this->rec_num;
    if(new_paths.first == new_paths.second) {
       // if the paths are the same, there is no recombination inside the anchor. check if there is a recombination between anchors now
        if ((this->paths & new_paths.first) == 0) {
           // there is a recombination between anchors, so we "reset" the current paths
            updated_paths = new_paths.first;
            rec_num++;
        } else {
            // there is no recombination between anchors, so we update the current paths
            updated_paths = this->paths & new_paths.first;
        }
    } else {
        // Otherwise, we have a recombinant anchor, we don't care about the recombination inside the anchor, we just "reset" the current paths
        updated_paths = new_paths.second;
    }
    return {
        this->score,
        this->source,
        updated_paths,
        rec_num
    };
}

void sort_anchor_indexes(const std::vector<Anchor>& items, std::vector<size_t>& indexes) {
    // Sort the indexes by read start ascending, and read end descending
    std::sort(indexes.begin(), indexes.end(), [&](const size_t& a, const size_t& b) {
        auto& a_item = items[a];
        auto& b_item = items[b];
        auto a_start = a_item.read_start();
        auto b_start = b_item.read_start();
        // a should be first if it starts earlier, or starts atthe same place and ends later.
        return (a_start < b_start || (a_start == b_start && a_item.read_end() > b_item.read_end()));
    });
}

transition_iterator lookback_transition_iterator(size_t max_lookback_bases,
                                                 size_t min_lookback_items,
                                                 size_t lookback_item_hard_cap) {

    
    // Capture all the arguments by value into a lambda
    transition_iterator iterator = [max_lookback_bases,
                                    min_lookback_items,
                                    lookback_item_hard_cap](const VectorView<Anchor>& to_chain,
                                                            const SnarlDistanceIndex& distance_index,
                                                            const HandleGraph& graph,
                                                            size_t max_indel_bases,
                                                            const transition_iteratee& callback) {

    


        // We want to consider all the important transitions in the graph of what
        // items can come before what other items. We aren't allowing any
        // transitions between items that overlap in the read. We're going through
        // the destination items in order by read start, so we should also keep a
        // list of them in order by read end, and sweep a cursor over that, so we
        // always know the fisrt item that overlaps with or passes the current
        // destination item, in the read. Then when we look for possible
        // predecessors of the destination item, we can start just before there and
        // look left.
        vector<size_t> read_end_order = sort_permutation(to_chain.begin(), to_chain.end(),
        [&](const Anchor& a, const Anchor& b) {
            return a.read_end() < b.read_end();
        });
        // We use first overlapping instead of last non-overlapping because we can
        // just initialize first overlapping at the beginning and be right.
        auto first_overlapping_it = read_end_order.begin();

        for (size_t i = 0; i < to_chain.size(); i++) {
            // For each item
            auto& here = to_chain[i];
            
            if (i > 0 && to_chain[i-1].read_start() > here.read_start()) {
                // The items are not actually sorted by read start
                throw std::runtime_error("lookback_transition_iterator: items are not sorted by read start");
            }
            
            while (to_chain[*first_overlapping_it].read_end() <= here.read_start()) {
                // Scan ahead through non-overlapping items that past-end too soon,
                // to the first overlapping item that ends earliest.
                // Ordering physics *should* constrain the iterator to not run off the end.
                ++first_overlapping_it;
                crash_unless(first_overlapping_it != read_end_order.end());
            }
            
#ifdef debug_chaining
            cerr << "Look at transitions to #" << i
                << " at " << here;
            cerr << endl;
#endif

#ifdef debug_chaining
            cerr << "\tFirst item overlapping #" << i << " beginning at " << here.read_start()
                 << " is #" << *first_overlapping_it << " past-ending at "
                 << to_chain[*first_overlapping_it].read_end() << " so start before there." << std::endl;
#endif
            
            // Set up lookback control algorithm.
            // Until we have looked at a certain number of items, we keep going
            // even if we meet other stopping conditions.
            size_t items_considered = 0;
            
            // Start considering predecessors for this item.
            auto predecessor_index_it = first_overlapping_it;
            while (predecessor_index_it != read_end_order.begin()) {
                --predecessor_index_it;
                
                // How many items have we considered before this one?
                size_t item_number = items_considered++;
                
                // For each source that ended before here started, in reverse order by end position...
                auto& source = to_chain[*predecessor_index_it];
                
#ifdef debug_chaining
                cerr << "\tConsider transition from #" << *predecessor_index_it << ": " << source << endl;
#endif
                
                // How far do we go in the read?
                size_t read_distance = get_read_distance(source, here);
                
                if (item_number > lookback_item_hard_cap) {
                    // This would be too many
#ifdef debug_chaining
                    cerr << "\t\tDisregard due to hitting lookback item hard cap" << endl;
#endif
                    break;
                }
                if (item_number >= min_lookback_items) {
                    // We have looked at enough predecessors that we might consider stopping.
                    // See if we should look back this far.
                    if (read_distance > max_lookback_bases) {
                        // This is further in the read than the real hard limit.
#ifdef debug_chaining
                        cerr << "\t\tDisregard due to read distance " << read_distance 
                             << " over limit " << max_lookback_bases << endl;
#endif
                        break;
                    } 
                }
                
                // Now it's safe to make a distance query
                
                // How far do we go in the graph?
                // Don't bother finding out exactly if it is too much longer than in the read.
                size_t graph_distance = get_graph_distance(source, here, distance_index, graph, 
                                                           read_distance + max_indel_bases);
                
                std::pair<int, int> scores = {std::numeric_limits<int>::min(), std::numeric_limits<int>::min()};
                if (read_distance != numeric_limits<size_t>::max() && graph_distance != numeric_limits<size_t>::max()) {
                    // Transition seems possible, so yield it.
                    callback(transition_info(*predecessor_index_it, i, graph_distance, read_distance));
                }
            } 
        }
    };

    return iterator;
}

transition_iterator zip_tree_transition_iterator(const std::vector<SnarlDistanceIndexClusterer::Seed>& seeds,
                                                 const ZipCodeTree& zip_code_tree,
                                                 size_t max_graph_lookback_bases,
                                                 size_t max_read_lookback_bases) {
    
    // TODO: Remove seeds because we only bring it here for debugging and it complicates the dependency relationships
    return [&seeds, &zip_code_tree, max_graph_lookback_bases, max_read_lookback_bases](
        const VectorView<Anchor>& to_chain,
        const SnarlDistanceIndex& distance_index,
        const HandleGraph& graph,
        size_t max_indel_bases,
        const transition_iteratee& callback
    ) {

        // We need a way to map from the seeds that zip tree thinks about to the anchors that we think about.
        // So we need to index the anchors by leading/trailing seed.
        // TODO: Should we make someone else do the indexing so we can make the Anchor not need to remember the seed?
        std::unordered_map<size_t, size_t> seed_to_starting;
        std::unordered_map<size_t, size_t> seed_to_ending;
        for (size_t anchor_num = 0; anchor_num < to_chain.size(); anchor_num++) {
            seed_to_starting[to_chain[anchor_num].seed_start()] = anchor_num;
            seed_to_ending[to_chain[anchor_num].seed_end()] = anchor_num;
        }

        // If we find we are actually walking through the graph in opposition
        // to the read, we need to defer transitions from source on the read
        // forward strand to dest on the read forward strand, so we can go them
        // in order along the read forward strand.
        // This holds source, dest, and graph distance.
        // We will fill it all in and then sort it by destination read position.
        std::vector<transition_info> all_transitions = 
            generate_zip_tree_transitions(seeds, zip_code_tree, max_graph_lookback_bases,
                                          max_read_lookback_bases, to_chain,
                                          seed_to_starting, seed_to_ending);

        // Sort the transitions so we handle them in an allowed order for dynamic programming.
        std::sort(all_transitions.begin(), all_transitions.end(), 
            [&](const transition_info& a, const transition_info& b) {
            // Return true if a's destination seed is before b's in the read, and false otherwise.
            return to_chain[a.to_anchor].read_start() < to_chain[b.to_anchor].read_start();
        });

        for (auto& transition : all_transitions) {
            callback(transition); 
        }
    };
}

std::vector<transition_info> generate_zip_tree_transitions(
    const std::vector<SnarlDistanceIndexClusterer::Seed>& seeds,
    const ZipCodeTree& zip_code_tree,
    size_t max_graph_lookback_bases,
    size_t max_read_lookback_bases,
    const VectorView<Anchor>& to_chain,
    const std::unordered_map<size_t, size_t>& seed_to_starting, 
    const std::unordered_map<size_t, size_t>& seed_to_ending) {

    std::vector<transition_info> all_transitions;
    // Save hopefully enough space for the transitions
    all_transitions.reserve(zip_code_tree.get_tree_size());

    for (auto seed_itr = zip_code_tree.begin(); seed_itr != zip_code_tree.end(); ++seed_itr) {
        // For each destination seed left to right
        vector<ZipCodeTree::oriented_seed_t> dest = *seed_itr;

        // Might be the start of an anchor if forward relative to the read,
        // or the end of an anchor if reverse relative to the read
        unordered_map<ZipCodeTree::oriented_seed_t, size_t> dest_anchors;
        for (const auto& dest_seed : dest) {
            auto anchor = dest_seed.is_reversed ? seed_to_ending.find(dest_seed.seed)
                                                : seed_to_starting.find(dest_seed.seed);
            if (anchor != (dest_seed.is_reversed ? seed_to_ending.end() : seed_to_starting.end())) {
                dest_anchors[dest_seed] = anchor->second;
            }
        }

        if (dest_anchors.empty()) {
            // Only find transitions if we find an anchor for this seed
            continue;
        }

#ifdef debug_transition
        std::cerr << "Destination seed";
        if (dest_anchors.size() > 1) {
            std::cerr << "s";
        }
        pos_t cur_pos;
        string rev_string = "";
        for (const auto& cur_dest_anchor : dest_anchors) {
            auto dest_seed = cur_dest_anchor.first;
            std::cerr << " (S" << dest_seed.seed << "/anchor #" << cur_dest_anchor.second << ")";
            cur_pos = seeds[dest_seed.seed].pos;
            rev_string = dest_seed.is_reversed ? "rev" : "";
        }
        std::cerr  << " at " << cur_pos << rev_string << std::endl;
#endif

        for (auto source = zip_code_tree.find_distances(seed_itr, max_graph_lookback_bases); 
             !source.done(); ++source) {
            // For each source seed right to left
            ZipCodeTree::seed_result_t source_seed = *source;
            for (const auto& cur_dest_anchor : dest_anchors) {
#ifdef debug_transition
                std::cerr << "\tSource seed S" << source_seed.seed << " " << seeds[source_seed.seed].pos 
                          << (source_seed.is_reversed ? "rev" : "") << " at distance " << source_seed.distance
                            << "/" << max_graph_lookback_bases;
#endif
                if (!source_seed.is_reversed && !cur_dest_anchor.first.is_reversed) {
                    // Both were traversed in the same orientation as the read.
                    // They might not be at anchor borders though, so check.
                    auto found_source_anchor = seed_to_ending.find(source_seed.seed);
                    if (found_source_anchor != seed_to_ending.end()) {
                        // We can transition between these seeds
                        // without jumping to/from the middle of an anchor.
#ifdef debug_transition
                        std::cerr << " is anchor #" << found_source_anchor->second << std::endl;
                        std::cerr << "\t\tFound transition from #" << found_source_anchor->second 
                                  << " to #" << cur_dest_anchor.second << std::endl;
#endif
                        add_transition_if_legal(all_transitions, to_chain, max_read_lookback_bases, 
                                                found_source_anchor->second, cur_dest_anchor.second, 
                                                source_seed.distance);
                    } else {
#ifdef debug_transition
                        std::cerr << " does not represent an anchor." << std::endl;
#endif
                    }
                } else if (source_seed.is_reversed && cur_dest_anchor.first.is_reversed) {
                    // Both were traversed in the opposite orientation as the read.
                    // We need to save them flipped for later.
                    auto found_source_anchor = seed_to_starting.find(source_seed.seed);
                    if (found_source_anchor != seed_to_starting.end()) {
                        // We can transition between these seeds
                        // without jumping to/from the middle of an anchor.
                        // Queue them up, flipped
                            
#ifdef debug_transition
                        std::cerr << " is anchor #" << found_source_anchor->second << std::endl;
                        std::cerr << "\t\tFound backward transition from #" << cur_dest_anchor.second << " to #"
                                  << found_source_anchor->second << std::endl;
#endif
                        add_transition_if_legal(all_transitions, to_chain, max_read_lookback_bases, 
                                                cur_dest_anchor.second, found_source_anchor->second,
                                                source_seed.distance);
                    } else {
#ifdef debug_transition
                        std::cerr << " does not represent an anchor." << std::endl;
#endif
                    }
                } else {
                    // We have a transition between different orientations
                    // relative to the read. That shouldn't happen.
                    crash_unless(source_seed.is_reversed == cur_dest_anchor.first.is_reversed);
                }
            }
        }
    }

    return all_transitions;
}

void add_transition_if_legal(vector<transition_info>& transitions, 
                             const VectorView<Anchor>& to_chain, size_t max_read_lookback_bases,
                             size_t from_anchor, size_t to_anchor, size_t graph_distance) {
    auto& source_anchor = to_chain[from_anchor];
    auto& dest_anchor = to_chain[to_anchor];

#ifdef debug_transition
    std::cerr << "Handle transition #" << from_anchor << " " << source_anchor
              << " to #" << to_anchor << " " << dest_anchor << std::endl;
    assert(graph_distance != std::numeric_limits<size_t>::max());
#endif

    size_t read_distance = get_read_distance(source_anchor, dest_anchor);
    if (read_distance == std::numeric_limits<size_t>::max()) {
        // Not reachable in read
#ifdef debug_transition
        std::cerr << "\tNot reachable in read." << std::endl;
#endif
        return;
    }

    if (read_distance > max_read_lookback_bases) {
        // Too far in read to consider
#ifdef debug_transition
        std::cerr << "\tToo far apart in read (" << read_distance
                  << "/" << max_read_lookback_bases << ")." << std::endl;
#endif
        return;
    }

    if (source_anchor.read_exclusion_end() > dest_anchor.read_exclusion_start()) {
        // The actual core anchor part is reachable in the read,
        // but we cut these down from overlapping minimizers.
#ifdef debug_transition
        std::cerr << "\tOriginally overlapped in read." << std::endl;
#endif
        return;
    }

    // The zipcode tree is about point positions,
    // but we need distances between whole anchors.
    // The stored zipcode positions will be at distances
    // from the start/end of the associated anchor.
    
    // If the offset between the zip code point
    // and the start of the destination is 0,
    // and between the zip code point and the end of the source is 0,
    // we subtract 0 from the measured distance.
    // Otherwise we need to subtract something.
    size_t distance_to_remove = dest_anchor.start_hint_offset() + source_anchor.end_hint_offset();

#ifdef debug_transition
    std::cerr << "\tZip code tree sees " << graph_distance
              << " but we should back out " << distance_to_remove << std::endl;
#endif

    if (distance_to_remove > graph_distance) {
        // We actually end further along the graph path to the next
        // thing than where the next thing starts, so we can't actually
        // get there.
#ifdef debug_transition
        std::cerr << "\tBacked out too much" << std::endl;
#endif
        return;
    }
    // Consume the length. 
    graph_distance -= distance_to_remove;

#ifdef debug_transition
    std::cerr << "\tZip code tree sees " << source_anchor << " and "
              << dest_anchor << " as " << graph_distance << " apart" << std::endl;
#endif
    transitions.emplace_back(from_anchor, to_anchor, graph_distance, read_distance);
}

/// Compute a gap score like minimap2.
///
/// They say they use the average anchor length, but really we need to use the
/// minimizer/base seed length here. Otherwise gaps cost more as your fragments
/// that you are chaining get longer, and cost more at chaining than at
/// fragmenting.
///
/// Returns a positive value (gap penalty).
int score_chain_gap(size_t distance_difference, size_t base_seed_length) {
    if (distance_difference == 0) {
        // Do nothing and score 0
        return 0;
    } else {
        // Compute the penalty
        return 0.01 * base_seed_length * distance_difference + 0.5 * log2(distance_difference);
    }
}

multipath_alignment_t fill_in_mp_aln(const VectorView<Anchor>& to_chain,
                                     const SnarlDistanceIndex& distance_index,
                                     const HandleGraph& graph,
                                     const ChainScoringScheme& scheme,
                                     const transition_iterator& for_each_transition,
                                     size_t max_indel_bases,
                                     bool show_work) {

    DiagramExplainer diagram(show_work);
    TSVExplainer dump(show_work, "chaindump");
    if (diagram) {
        diagram.add_globals({{"rankdir", "LR"}});
    }
   
    if (show_work) {
        cerr << "Chaining group of " << to_chain.size() << " items" << endl;
    }

    crash_unless(scheme.recombination_penalty >= 0);
    crash_unless(scheme.consistency_bonus >= 0);

    // This is what we will build up
    multipath_alignment_t multipath_aln;

    // Compute a base seed average length.
    // TODO: Weight anchors differently?
    // TODO: Will this always be the same for all anchors in practice?
    size_t base_seed_length = 0;
    for (size_t i = 0; i < to_chain.size(); i++) {
        base_seed_length += to_chain[i].base_seed_length();
        // "subpaths" are just what we use to represent anchors
        subpath_t* subpath = multipath_aln.add_subpath();
        // How many points this anchor is worth to collect
        subpath->set_score(to_chain[i].score() + scheme.item_bonus);
        // Add a fake position to this subpath to store the anchor ID
        position_t* position = subpath->mutable_path()->add_mapping()->mutable_position();
        position->set_node_id(i);
        position->set_is_reverse(false);
        position->set_offset(0);
    }
    base_seed_length /= to_chain.size();

    // We will run this over every transition in a good DP order.
    auto iteratee = [&](const transition_info& transition) {
        if (show_work) {
#ifdef debug_dp
            cerr << "DP: " << transition.from_anchor << "->" << transition.to_anchor 
                 << " rd " << transition.read_distance << " gd " << transition.graph_distance << endl;
#endif
        }
        
        // For each item
        auto& here = to_chain[transition.to_anchor];
        
        // For each source we could come from
        auto& source = to_chain[transition.from_anchor];
            
        if (show_work) {
#ifdef debug_dp
            cerr << "\t\tCome from " << source << " to " << here << endl;
#endif
        }
            
        // Decide how much length changed
        size_t indel_length = (transition.read_distance > transition.graph_distance) ? transition.read_distance - transition.graph_distance 
                                                                                     : transition.graph_distance - transition.read_distance;
        
        if (show_work) {
#ifdef debug_dp
            cerr << "\t\t\tFor read distance " << transition.read_distance << " and graph distance " << transition.graph_distance
                 << " an indel of length " << indel_length
                 << ((transition.read_distance > transition.graph_distance) ? " seems plausible" : " would be required") << endl;
#endif
        }
            
        if (indel_length <= max_indel_bases) {
            // Assign points for the assumed matches in the transition, and charge for the indel.
            //
            // The Minimap2 paper
            // <https://doi.org/10.1093/bioinformatics/bty191> at 2.1.1 says
            // that we ought to assign "α(j,i)=min{min{yi−yj,xi−xj},wi} is the
            // number of matching bases between the two anchors", minus the gap
            // penalty. Here, i is the destination anchor and j is the
            // predecessor, and x and y are read and query positions of the
            // *final* base in the anchor, while w is anchor width.
            //
            // As written, the gloss isn't really true; the number of matching
            // bases between the two anchors isn't bounded below by the width
            // of the second anchor. It looks more like we are counting the
            // number of new matching bases in the destination anchor that are
            // not overlapping matching bases in the source anchor.
            //
            // Our distances are between the end of the previous anchor and the
            // start of this one (not the end as in Minimap2's formulation).
            // And our anchors also thus never overlap. So we can just always
            // use the length of the destination anchor.
            //
            // But we account for anchor length in the item points, so don't use it
            // here.
            int jump_points = -score_chain_gap(indel_length, base_seed_length) * scheme.gap_scale;

            // Penalize transitions which force a recombination
            if (source.anchor_end_paths() & here.anchor_start_paths() == 0) {
                jump_points -= scheme.recombination_penalty;
            }

            // Connect these anchors
            connection_t* connection = multipath_aln.mutable_subpath(transition.from_anchor)->add_connection();
            connection->set_score(jump_points);
            connection->set_next(transition.to_anchor);
            
            /// TODO: Consistency bonus
            
            if (dump) {
                // Dump a TSV of source anchor description, dest anchor description, and score of edge.
                dump.line();
                std::stringstream source_stream;
                source_stream << source;
                dump.field(source_stream.str());
                std::stringstream here_stream;
                here_stream << here;
                dump.field(here_stream.str());
                dump.field(jump_points);
            }
        } else {
            if (show_work) {
#ifdef debug_dp
                cerr << "\t\tTransition is impossible." << endl;
#endif
            }
        }
    };

    // Run our DP step over all the transitions.
    for_each_transition(to_chain,
                        distance_index,
                        graph,
                        max_indel_bases,
                        iteratee);

    return multipath_aln;
}

SubchainGroup split_up_subchains(const size_t& anchor_count,
                                 const vector<Alignment>& tracebacks,
                                 const vector<pair<uint32_t, uint32_t>>& alternatives) {
    SubchainGroup output;
    // For each anchor (by index), which other anchors can it connect to?
    vector<vector<size_t>> outgoing_edges(anchor_count);
    // For each anchor (by index), how many sources does it have?
    vector<size_t> source_count(anchor_count, 0);

    // Where to search for subchains (we must guarantee they start in read order)
    priority_queue<size_t, vector<size_t>, std::greater<size_t>> trace_from;
    
    // Extract raw lists of anchors from the tracebacks
    for (const auto& cur_trace : tracebacks) {
        // Loop over all "subpaths" (anchor IDs) in the traceback
        size_t num_ids = cur_trace.path().mapping().size();
        // This will start a chain trace
        trace_from.emplace(cur_trace.path().mapping(0).position().node_id());
#ifdef debug_chaining
        cerr << "Chain traceforwards may start from " << cur_trace.path().mapping(0).position().node_id() << endl;
#endif
        for (size_t i = 0; i < num_ids - 1; i++) {
            // Remember that this pair of anchors can connect
            size_t prev = cur_trace.path().mapping(i).position().node_id();
            size_t next = cur_trace.path().mapping(i+1).position().node_id();
            source_count[next]++;
            outgoing_edges[prev].emplace_back(next);
        }
    }

    // Now check in on the extra edges.
    for (const auto& extra_edge : alternatives) {
        if ((outgoing_edges[extra_edge.first].size() + source_count[extra_edge.first]) > 0
            && (outgoing_edges[extra_edge.second].size() + source_count[extra_edge.second]) > 0) {
            // Both sides of this edge are used, so it must connect two different subchains
            // Add this as another possible next
            outgoing_edges[extra_edge.first].emplace_back(extra_edge.second);
            source_count[extra_edge.second]++;
        }
    }

#ifdef debug_chaining
    for (size_t i = 0; i < anchor_count; i++) {
        cerr << i << " has " << source_count[i] << " sources, and outgoing edges to ";
        for (const auto& next : outgoing_edges[i]) {
            cerr << next << " ";
        }
        cerr << endl;
    }
#endif

    // Set up the subchains
    // For each anchor (by index), which subchain did it end up in?
    vector<size_t> subchain_id(anchor_count, std::numeric_limits<size_t>::max());
    while (!trace_from.empty()) {
        // Start trace for this subchain, until we hit into an endpoint
        size_t cur_anchor_id = trace_from.top();
        trace_from.pop();
        if (subchain_id[cur_anchor_id] != std::numeric_limits<size_t>::max()) {
            // This one was already put as part of an earlier subchain
            continue;
        }

        // Create a new subchain to trace into
        size_t cur_subchain_id = output.subchains.size();
        output.subchains.emplace_back();
        output.subchains.back().emplace_back(cur_anchor_id);
        subchain_id[cur_anchor_id] = cur_subchain_id;

#ifdef debug_chaining
        cerr << "Assign " << cur_anchor_id << " to subchain " << cur_subchain_id << endl;
#endif

        while (true) { 
            if (outgoing_edges[cur_anchor_id].empty()) {
                // This anchor can't trace outwards at all
                break;
            }           
            // If we've reached a decision point, then save all next edges
            // We have reached the end of one subchain
            if (outgoing_edges[cur_anchor_id].size() > 1 
                || outgoing_edges[outgoing_edges[cur_anchor_id].front()].size() > 1
                || source_count[outgoing_edges[cur_anchor_id].front()] > 1) {
                for (const auto& next : outgoing_edges[cur_anchor_id]) {
                    // We need to start a new trace from here
#ifdef debug_chaining
                    cerr << "Chain traceforwards may start from " << next << endl;
#endif
                    trace_from.emplace(next);
                }
                break;
            }
            
            // Otherwise, trace into the next edge
            cur_anchor_id = outgoing_edges[cur_anchor_id].front();
            output.subchains.back().emplace_back(cur_anchor_id);
            subchain_id[cur_anchor_id] = cur_subchain_id;
#ifdef debug_chaining
            cerr << "Assign " << cur_anchor_id << " to subchain " << cur_subchain_id << endl;
#endif
        }
    }

    // Make all inter-subchain connections necessary
    for (size_t i = 0; i < output.subchains.size(); i++) {
        // Loop over any way to exit this subchain
        for (const auto& next : outgoing_edges[output.subchains[i].back()]) {
            if (next != std::numeric_limits<size_t>::max()) {
                output.connections.emplace_back(i, subchain_id[next]);
#ifdef debug_chaining
                cerr << "Connect subchains " << i << " -> " << subchain_id[next] << endl;
#endif
            }
        }
    }

    return output;
}

SubchainGroup find_best_chains(const VectorView<Anchor>& to_chain,
                               const SnarlDistanceIndex& distance_index,
                               const HandleGraph& graph,
                               const ChainScoringScheme& scheme,
                               size_t max_chains,
                               const transition_iterator& for_each_transition,
                               size_t max_indel_bases,
                               bool show_work) {

    if (to_chain.empty()) {
        // Nothing to chain
        return SubchainGroup();
    }
        
    // First, we build the multipath alignment
    multipath_alignment_t mp_aln = fill_in_mp_aln(to_chain,
                                                  distance_index,
                                                  graph,
                                                  scheme,
                                                  for_each_transition,
                                                  max_indel_bases,
                                                  show_work);
    // Then do the tracebacks
    vector<pair<uint32_t, uint32_t>> alternatives;
    vector<Alignment> tracebacks = optimal_alignments_with_disjoint_subpaths(mp_aln, max_chains, &alternatives);
    
    if (tracebacks.empty()) {
        // Somehow we got nothing
        return SubchainGroup();
    }

    SubchainGroup output = split_up_subchains(to_chain.size(), tracebacks, alternatives);
    // Also remember its maximum score
    output.max_sparse_chain_score = tracebacks.front().score();
    return output;
}

SparseAnchorChain find_best_chain(const VectorView<Anchor>& to_chain,
                                  const SnarlDistanceIndex& distance_index,
                                  const HandleGraph& graph,
                                  const ChainScoringScheme& scheme,
                                  const transition_iterator& for_each_transition,
                                  size_t max_indel_bases) {
    SubchainGroup group = find_best_chains(to_chain,
                                           distance_index,
                                           graph,
                                           scheme,
                                           1,
                                           for_each_transition,
                                           max_indel_bases);
    if (group.subchains.empty()) {
        // We got nothing
        return SparseAnchorChain();
    }

    SparseAnchorChain output;
    output.anchors = group.subchains.front();
    output.chain_score = group.max_sparse_chain_score;
    
    return output;
}

//#define skip_zipcodes
//#define debug
//#define double_check_distances
//#define stop_on_mismatch
//#define replace_on_mismatch
size_t get_graph_distance(const Anchor& from, const Anchor& to, const SnarlDistanceIndex& distance_index,
                          const HandleGraph& graph, size_t distance_limit) {
    auto from_pos = from.graph_end();
    auto& to_pos = to.graph_start();
    
    auto* from_hint = from.end_hint();
    auto* to_hint = to.start_hint();
    
    size_t distance;
    
#ifdef skip_zipcodes
    if (false) {
#else
    if (from_hint && to_hint) {
#endif
#ifdef debug
        #pragma omp critical (cerr)
        {
            std::cerr << "Finding distance from " << from_pos << " to " << to_pos << " using hints ";
            from_hint->dump(std::cerr);
            std::cerr << " and ";
            to_hint->dump(std::cerr);
            std::cerr << std::endl;
        }
#endif
    
        // Can use zip code based oriented distance
        distance = ZipCode::minimum_distance_between(*from_hint, from_pos, 
                                                     *to_hint, to_pos,
                                                     distance_index,
                                                     distance_limit,
                                                     false, 
                                                     &graph);

#ifdef debug
        #pragma omp critical (cerr)
        std::cerr << "Zipcodes report " << distance << std::endl;
#endif

#ifdef double_check_distances
        // Make sure the minimizers aren't way off from the distance index.
        size_t check_distance = distance_index.minimum_distance(
            id(from_pos), is_rev(from_pos), offset(from_pos),
            id(to_pos), is_rev(to_pos), offset(to_pos),
            false, &graph);

        if (check_distance > distance) {
#ifdef debug
            #pragma omp critical (cerr)
            std::cerr << "Distance index reports " << check_distance << " instead" << std::endl;
#endif  
          
#ifdef stop_on_mismatch
            throw std::runtime_error("Zipcode distance mismatch");
#endif
#ifdef replace_on_mismatch
            distance = check_distance;
#endif
    }

#endif
    } else {
        // Query the distance index directly.
        distance = distance_index.minimum_distance(
            id(from_pos), is_rev(from_pos), offset(from_pos),
            id(to_pos), is_rev(to_pos), offset(to_pos),
            false, &graph);
    }
    if (distance > distance_limit) {
        // Zip code logic can have to compute a number over the limit,
        // and in that case will return it.
        // Cut it off here.
        distance = std::numeric_limits<size_t>::max();
    }
    return distance;
}

size_t get_read_distance(const Anchor& from, const Anchor& to) {
    if (to.read_start() < from.read_end()) {
        return std::numeric_limits<size_t>::max();
    }
    return to.read_start() - from.read_end();
}

}
}
