/**
 * \file chain_items.cpp
 * Non-template function implementations for chaining pieces of a read-to-graph alignment.
 */


#include "chain_items.hpp"

#include <handlegraph/algorithms/dijkstra.hpp>
#include <structures/immutable_list.hpp>
#include <structures/min_max_heap.hpp>

//#define debug_chaining

namespace vg {
namespace algorithms {

using namespace std;

ostream& operator<<(ostream& out, const Anchor& anchor) {
    return out << "{R:" << anchor.read_start() << "=G:" << anchor.graph_start() << "*" << anchor.length() << "}";
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
    }
}

TracedScore TracedScore::score_from(const vector<TracedScore>& options, size_t option_number) {
    TracedScore got = options[option_number];
    got.source = option_number;
    return got;
}

TracedScore TracedScore::add_points(int adjustment) const {
    return {this->score + adjustment, this->source};
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

void sort_and_shadow(const std::vector<Anchor>& items, std::vector<size_t>& indexes) {
    
    // Sort everything by read start ascending, and read end descending
    sort_anchor_indexes(items, indexes);
    
    // Keep a collection of the diagonals that are already represented,
    // and the read end position of the latest-ending item on those pairs that
    // we have taken. A diagonal is defined as a graph node ID, a graph strand,
    // and the difference between the graph offset and the read position. So we
    // can represent them with pos_t, and subtract the read position out of the
    // stored offset to make them.
    std::unordered_map<pos_t, size_t> diagonal_progress;
    
    // Scan through and make a new collection of indexes, keeping the first on
    // any pair of diagonals, which will thus be the one with the earliest
    // start, and within those the latest end. Since we need to keep items
    // which partially overlap but don't contain each other, we also keep an
    // item if it is the new latest-ending thing we've seen for a pair of
    // diagonals.
    std::vector<size_t> kept_indexes;
    kept_indexes.reserve(indexes.size());
    for (auto i : indexes) {
        // For each item we might keep
        auto& item = items[i];
        
        // Prepare the key of the diagonals it visits
        pos_t diagonal = item.graph_start();
        // Make the offsets store a difference between graph and read offset so
        // they really represent diagonals.
        get_offset(diagonal) -= item.read_start();
        
        auto& furthest_read_end = diagonal_progress[diagonal];
        if (furthest_read_end < item.read_end()) {
            // This is the first, or latest-ending, item seen on this diagonal.
            // If there was an earlier-ending item taken, we know it started before this one, because of iteration order.
            // So take this item.
            kept_indexes.push_back(i);
            // And record that we got out this far
            furthest_read_end = item.read_end();
#ifdef debug_chaining
            std::cerr << "Keep " << item << " which gets us to R" << furthest_read_end << " on diagonal " << diagonal << std::endl;
#endif
        } else {
#ifdef debug_chaining
            std::cerr << "Discard " << item << " as shadowed because we already got to R" << furthest_read_end << " on diagonal " << diagonal << std::endl;
#endif
        }
    }
    
    // Replace the indexes with the sorted and deduplicated ones.
    indexes = std::move(kept_indexes);
}

void sort_and_shadow(std::vector<Anchor>& items) {
    // Use the index-based implementation and then apply those indexes
    std::vector<size_t> indexes = range_vector(items.size());
    sort_and_shadow(items, indexes);
    std::vector<Anchor> kept_items;
    kept_items.reserve(indexes.size());
    for (auto& index : indexes) {
        kept_items.emplace_back(std::move(items[index]));
    }
    items = std::move(kept_items);
}

TracedScore chain_items_dp(vector<TracedScore>& chain_scores,
                           const VectorView<Anchor>& to_chain,
                           const SnarlDistanceIndex& distance_index,
                           const HandleGraph& graph,
                           int gap_open,
                           int gap_extension,
                           size_t max_lookback_bases,
                           size_t min_lookback_items,
                           size_t lookback_item_hard_cap,
                           size_t initial_lookback_threshold,
                           double lookback_scale_factor,
                           double min_good_transition_score_per_base,
                           int item_bonus,
                           size_t max_indel_bases) {
    
    DiagramExplainer diagram;
    diagram.add_globals({{"rankdir", "LR"}});
    
#ifdef debug_chaining
    cerr << "Chaining group of " << to_chain.size() << " items" << endl;
#endif
    
    // We want to consider all the important transitions in the graph of what
    // items can come before what other items. We aren't allowing any
    // transitions between items that overlap in the read. We're going through
    // the destination items in order by read start, so we should also keep a
    // list of them in order by read end, and sweep a cursor over that, so we
    // always know the fisrt item that overlaps with or passes the current
    // destination item, in the read. Then when we look for possible
    // predecessors of the destination item, we can start just before there and
    // look left.
    vector<size_t> read_end_order = sort_permutation(to_chain.begin(), to_chain.end(), [&](const Anchor& a, const Anchor& b) {
        return a.read_end() < b.read_end();
    });
    // We use first overlapping instead of last non-overlapping because we can
    // just initialize first overlapping at the beginning and be right.
    auto first_overlapping_it = read_end_order.begin();
    
    // Make our DP table big enough
    chain_scores.clear();
    chain_scores.resize(to_chain.size(), TracedScore::unset());
    
    // What's the winner so far?
    TracedScore best_score = TracedScore::unset();
    
    for (size_t i = 0; i < to_chain.size(); i++) {
        // For each item
        auto& here = to_chain[i];
        
        if (i > 0 && to_chain[i-1].read_start() > here.read_start()) {
            // The items are not actually sorted by read start
            throw std::runtime_error("chain_items_dp: items are not sorted by read start");
        }
        
        while (to_chain[*first_overlapping_it].read_end() <= here.read_start()) {
            // Scan ahead through non-overlapping items that past-end too soon,
            // to the first overlapping item that ends earliest.
            // Ordering physics *should* constrain the iterator to not run off the end.
            ++first_overlapping_it;
            assert(first_overlapping_it != read_end_order.end());
        }
        
        // How many points is it worth to collect?
        auto item_points = here.score() + item_bonus;
        
        std::string here_gvnode = "i" + std::to_string(i);
        
        // If we come from nowhere, we get those points.
        chain_scores[i] = std::max(chain_scores[i], {item_points, TracedScore::nowhere()});
        
#ifdef debug_chaining
        cerr << "Look at transitions to #" << i
            << " at " << here;
        cerr << endl;
#endif

#ifdef debug_chaining
        cerr << "\tFirst item overlapping #" << i << " beginning at " << here.read_start() << " is #" << *first_overlapping_it << " past-ending at " << to_chain[*first_overlapping_it].read_end() << " so start before there." << std::endl;
#endif
        
        // Set up lookback control algorithm.
        // Until we have looked at a certain number of items, we keep going
        // even if we meet other stopping conditions.
        size_t items_considered = 0;
        // If we are looking back further than this
        size_t lookback_threshold = initial_lookback_threshold;
        // And a gooid score has been found, stop
        bool good_score_found = false;
        // A good score will be positive and have a transition component that
        // looks good relative to how far we are looking back. The further we
        // look back the lower our transition score standards get, so remember
        // the best one we have seen so far in case the standard goes below it. 
        int best_transition_found = std::numeric_limits<int>::min();
        
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
                cerr << "\t\tDisregard due to read distance " << read_distance << " over limit " << max_lookback_bases << endl;
#endif
                    break;
                } else if (read_distance > lookback_threshold && good_score_found) {
                    // We already found something good enough.
#ifdef debug_chaining
                cerr << "\t\tDisregard due to read distance " << read_distance << " over threashold " << lookback_threshold << " and good score already found" << endl;
#endif
                    break;
                }
            }
            if (read_distance > lookback_threshold && !good_score_found) {
                // We still haven't found anything good, so raise the threshold.
                lookback_threshold *= lookback_scale_factor;
            }
            
            // Now it's safe to make a distance query
#ifdef debug_chaining
            cerr << "\t\tCome from score " << chain_scores[*predecessor_index_it]
                << " across " << source << " to " << here << endl;
#endif
            
            // We will actually evaluate the source.
            
            // How far do we go in the graph? Don't bother finding out exactly if it is too much longer than in the read.
            size_t graph_distance = get_graph_distance(source, here, distance_index, graph, read_distance + max_indel_bases);
            
            // How much does it pay (+) or cost (-) to make the jump from there
            // to here?
            // Don't allow the transition if it seems like we're going the long
            // way around an inversion and needing a huge indel.
            int jump_points;
            
            if (read_distance == numeric_limits<size_t>::max()) {
                // Overlap in read, so not allowed.
                jump_points = std::numeric_limits<int>::min();
            } else if (graph_distance == numeric_limits<size_t>::max()) {
                // No graph connection
                jump_points = std::numeric_limits<int>::min();
            } else {
                // Decide how much length changed
                size_t indel_length = (read_distance > graph_distance) ? read_distance - graph_distance : graph_distance - read_distance;
                
#ifdef debug_chaining
                cerr << "\t\t\tFor read distance " << read_distance << " and graph distance " << graph_distance << " an indel of length " << indel_length << " would be required" << endl;
#endif

                if (indel_length > max_indel_bases) {
                    // Don't allow an indel this long
                    jump_points = std::numeric_limits<int>::min();
                } else {
                    // Then charge for that indel
                    jump_points = score_gap(indel_length, gap_open, gap_extension);
                }
            }
            
            // And how much do we end up with overall coming from there.
            int achieved_score;
            
            if (jump_points != numeric_limits<int>::min()) {
                // Get the score we are coming from
                TracedScore source_score = TracedScore::score_from(chain_scores, *predecessor_index_it);
                
                // And the score with the transition and the points from the item
                TracedScore from_source_score = source_score.add_points(jump_points + item_points);
                
                // Remember that we could make this jump
                chain_scores[i] = std::max(chain_scores[i], from_source_score);
                                               
#ifdef debug_chaining
                cerr << "\t\tWe can reach #" << i << " with " << source_score << " + " << jump_points << " from transition + " << item_points << " from item = " << from_source_score << endl;
#endif
                if (from_source_score.score > 0) {
                    // Only explain edges that were actual candidates since we
                    // won't let local score go negative
                    
                    std::string source_gvnode = "i" + std::to_string(*predecessor_index_it);
                    // Suggest that we have an edge, where the edges that are the best routes here are the most likely to actually show up.
                    diagram.suggest_edge(source_gvnode, here_gvnode, here_gvnode, from_source_score.score, {
                        {"label", std::to_string(jump_points)},
                        {"weight", std::to_string(std::max<int>(1, from_source_score.score))}
                    });
                }
                
                achieved_score = from_source_score.score;
            } else {
#ifdef debug_chaining
                cerr << "\t\tTransition is impossible." << endl;
#endif
                achieved_score = std::numeric_limits<size_t>::min();
            }
            
            // Note that we checked out this transition and saw the observed scores and distances.
            best_transition_found = std::max(best_transition_found, jump_points);
            if (achieved_score > 0 && best_transition_found >= min_good_transition_score_per_base * std::max(read_distance, graph_distance)) {
                // We found a jump that looks plausible given how far we have searched, so we can stop searching way past here.
                good_score_found = true;
            }
        }
        
#ifdef debug_chaining
        cerr << "\tBest way to reach #" << i << " is " << chain_scores[i] << endl;
#endif
        
        std::stringstream label_stream;
        label_stream << "#" << i << " " << here << " = " << item_points << "/" << chain_scores[i].score;
        diagram.add_node(here_gvnode, {
            {"label", label_stream.str()}
        });
        auto graph_start = here.graph_start();
        std::string graph_gvnode = "n" + std::to_string(id(graph_start)) + (is_rev(graph_start) ? "r" : "f");
        diagram.ensure_node(graph_gvnode, {
            {"label", std::to_string(id(graph_start)) + (is_rev(graph_start) ? "-" : "+")},
            {"shape", "box"}
        });
        // Show the item as connected to its source graph node
        diagram.add_edge(here_gvnode, graph_gvnode, {{"color", "gray"}});
        // Make the next graph node along the same strand
        std::string graph_gvnode2 = "n" + std::to_string(id(graph_start) + (is_rev(graph_start) ? -1 : 1)) + (is_rev(graph_start) ? "r" : "f");
        diagram.ensure_node(graph_gvnode2, {
            {"label", std::to_string(id(graph_start) + (is_rev(graph_start) ? -1 : 1)) + (is_rev(graph_start) ? "-" : "+")},
            {"shape", "box"}
        });
        // And show them as connected. 
        diagram.ensure_edge(graph_gvnode, graph_gvnode2, {{"color", "gray"}});
        
        // See if this is the best overall
        best_score.max_in(chain_scores, i);
        
#ifdef debug_chaining
        cerr << "\tBest chain end so far: " << best_score << endl;
#endif
        
    }
    
    return best_score;
}

vector<pair<vector<size_t>, int>> chain_items_traceback(const vector<TracedScore>& chain_scores,
                                                        const VectorView<Anchor>& to_chain,
                                                        const TracedScore& best_past_ending_score_ever,
                                                        int item_bonus,
                                                        size_t max_tracebacks) {
    
    // We will fill this in with all the tracebacks, and then sort and truncate.
    vector<pair<vector<size_t>, int>> tracebacks;
    tracebacks.reserve(chain_scores.size());
    
    // Get all of the places to start tracebacks, in score order.
    std::vector<size_t> starts_in_score_order;
    starts_in_score_order.resize(chain_scores.size());
    for (size_t i = 0; i < starts_in_score_order.size(); i++) {
        starts_in_score_order[i] = i;
    }
    std::sort(starts_in_score_order.begin(), starts_in_score_order.end(), [&](const size_t& a, const size_t& b) {
        // Return true if item a has a better score than item b and should come first.
        return chain_scores[a] > chain_scores[b];
    });
    
    // To see if an item is used we have this bit vector.
    vector<bool> item_is_used(chain_scores.size(), false);
    
    for (auto& trace_from : starts_in_score_order) {
        if (item_is_used[trace_from]) {
            continue;
        }
        // For each unused item in score order, start a traceback stack (in reverse order)
        std::vector<size_t> traceback;
        traceback.push_back(trace_from);
        // Track the penalty we are off optimal for this traceback
        int penalty = best_past_ending_score_ever - chain_scores[trace_from];
        size_t here = trace_from;
        while (here != TracedScore::nowhere()) {
            // Mark here as used. Happens once per item, and so limits runtime.
            item_is_used[here] = true;
            size_t next = chain_scores[here].source;
            if (next != TracedScore::nowhere()) {
                if (item_is_used[next]) {
                    // We need to stop early and accrue an extra penalty.
                    // Take away all the points we got for coming from there and being ourselves.
                    penalty += chain_scores[here].score;
                    // But then re-add our score for just us
                    penalty -= (to_chain[here].score() + item_bonus);
                    // TODO: Score this more simply.
                    // TODO: find the dege to nowhere???
                    break;
                } else {
                    // Add to the traceback
                    traceback.push_back(next);
                }
            }
            here = next;
        }
        // Now put the traceback in the output list
        tracebacks.emplace_back();
        tracebacks.back().second = penalty;
        // Make sure to order the steps left to right, and not right to left as we generated them.
        std::copy(traceback.rbegin(), traceback.rend(), std::back_inserter(tracebacks.back().first));
    }
    
    // Sort the tracebacks by penalty, ascending
    std::sort(tracebacks.begin(), tracebacks.end(), [](const std::pair<std::vector<size_t>, int>& a, const std::pair<std::vector<size_t>, int>& b) {
        // Return true if a has the smaller penalty and belongs first
        return a.second < b.second;
    });
    
    if (tracebacks.size() > max_tracebacks) {
        // Limit to requested number
        tracebacks.resize(max_tracebacks);
    }

    return tracebacks;
}

vector<pair<int, vector<size_t>>> find_best_chains(const VectorView<Anchor>& to_chain,
                                                   const SnarlDistanceIndex& distance_index,
                                                   const HandleGraph& graph,
                                                   int gap_open,
                                                   int gap_extension,
                                                   size_t max_chains,
                                                   size_t max_lookback_bases,
                                                   size_t min_lookback_items,
                                                   size_t lookback_item_hard_cap,
                                                   size_t initial_lookback_threshold,
                                                   double lookback_scale_factor,
                                                   double min_good_transition_score_per_base,
                                                   int item_bonus,
                                                   size_t max_indel_bases) {
                                                                         
    if (to_chain.empty()) {
        return {{0, vector<size_t>()}};
    }
        
    // We actually need to do DP
    vector<TracedScore> chain_scores;
    TracedScore best_past_ending_score_ever = chain_items_dp(chain_scores,
                                                             to_chain,
                                                             distance_index,
                                                             graph,
                                                             gap_open,
                                                             gap_extension,
                                                             max_lookback_bases,
                                                             min_lookback_items,
                                                             lookback_item_hard_cap,
                                                             initial_lookback_threshold,
                                                             lookback_scale_factor,
                                                             min_good_transition_score_per_base,
                                                             item_bonus,
                                                             max_indel_bases);
    // Then do the tracebacks
    vector<pair<vector<size_t>, int>> tracebacks = chain_items_traceback(chain_scores, to_chain, best_past_ending_score_ever, item_bonus, max_chains);
    
    if (tracebacks.empty()) {
        // Somehow we got nothing
        return {{0, vector<size_t>()}};
    }
        
    // Convert form traceback and penalty to score and traceback.
    // Everything is already sorted.
    vector<pair<int, vector<size_t>>> to_return;
    to_return.reserve(tracebacks.size());
    for (auto& traceback : tracebacks) {
        // Move over the list of items and convert penalty to score
        to_return.emplace_back(best_past_ending_score_ever.score - traceback.second, std::move(traceback.first));
    }
    
    return to_return;
}

pair<int, vector<size_t>> find_best_chain(const VectorView<Anchor>& to_chain,
                                          const SnarlDistanceIndex& distance_index,
                                          const HandleGraph& graph,
                                          int gap_open,
                                          int gap_extension,
                                          size_t max_lookback_bases,
                                          size_t min_lookback_items,
                                          size_t lookback_item_hard_cap,
                                          size_t initial_lookback_threshold,
                                          double lookback_scale_factor,
                                          double min_good_transition_score_per_base,
                                          int item_bonus,
                                          size_t max_indel_bases) {
                                                                 
    return find_best_chains(
        to_chain,
        distance_index,
        graph,
        gap_open,
        gap_extension,
        1,
        max_lookback_bases,
        min_lookback_items,
        lookback_item_hard_cap,
        initial_lookback_threshold,
        lookback_scale_factor,
        min_good_transition_score_per_base,
        item_bonus,
        max_indel_bases
    ).front();
}

int score_best_chain(const VectorView<Anchor>& to_chain, const SnarlDistanceIndex& distance_index, const HandleGraph& graph, int gap_open, int gap_extension) {
    
    if (to_chain.empty()) {
        return 0;
    } else {
        // Do the DP but without the traceback.
        vector<TracedScore> chain_scores;
        TracedScore winner = algorithms::chain_items_dp(chain_scores, to_chain, distance_index, graph, gap_open, gap_extension);
        return winner.score;
    }
}

//#define skip_zipcodes
//#define debug
#define double_check_distances
//#define stop_on_mismatch
#define replace_on_mismatch
size_t get_graph_distance(const Anchor& from, const Anchor& to, const SnarlDistanceIndex& distance_index, const HandleGraph& graph, size_t distance_limit) {
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
        // Zip code logic can have to compute a number over the limit, and in that case will return it.
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
