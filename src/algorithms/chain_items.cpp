/**
 * \file chain_items.cpp
 * Non-template function implementations for chaining pieces of a read-to-graph alignment.
 */


#include "chain_items.hpp"

#include <handlegraph/algorithms/dijkstra.hpp>

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

void sort_and_shadow(const std::vector<Anchor>& items, std::vector<size_t>& indexes) {
    
    // Sort the indexes by read start ascending, and read end descending
    std::sort(indexes.begin(), indexes.end(), [&](const size_t& a, const size_t& b) {
        auto& a_item = items[a];
        auto& b_item = items[b];
        auto a_start = a_item.read_start();
        auto b_start = b_item.read_start();
        // a should be first if it starts earlier, or starts atthe same place and ends later.
        return (a_start < b_start || (a_start == b_start && a_item.read_end() > b_item.read_end()));
    });
    
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

TracedScore chain_items_dp(vector<TracedScore>& best_chain_score,
                           const VectorView<Anchor>& to_chain,
                           const SnarlDistanceIndex& distance_index,
                           const HandleGraph& graph,
                           int gap_open,
                           int gap_extension,
                           size_t max_lookback_bases,
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
    best_chain_score.resize(to_chain.size(), TracedScore::unset());
    
    // What's the winner so far?
    TracedScore best_score = TracedScore::unset();
    
    for (size_t i = 0; i < to_chain.size(); i++) {
        // For each item
        auto& here = to_chain[i];
        
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
        best_chain_score[i] = std::max(best_chain_score[i], {item_points, TracedScore::nowhere()});
        
#ifdef debug_chaining
        cerr << "Look at transitions to #" << i
            << " at " << here;
        cerr << endl;
#endif

#ifdef debug_chaining
        cerr << "\tFirst item overlapping #" << i << " beginning at " << here.read_start() << " is #" << *first_overlapping_it << " past-ending at " << to_chain[*first_overlapping_it].read_end() << " so start before there." << std::endl;
#endif
        
        // Set up lookback control algorithm
        // If we are looking back forther than this
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
            // For each source that ended before here started, in reverse order by end position...
            auto& source = to_chain[*predecessor_index_it];
            
#ifdef debug_chaining
            cerr << "\tConsider transition from #" << *predecessor_index_it << ": " << source << endl;
#endif

            // How far do we go in the read?
            size_t read_distance = get_read_distance(source, here);

            // See if we should look back this far.
            if (read_distance > max_lookback_bases) {
                // This is further in the read than the real hard limit.
                break;
            } else if (read_distance > lookback_threshold && good_score_found) {
                // We already found something good enough.
                break;
            }
            if (read_distance > lookback_threshold && !good_score_found) {
                // We still haven't found anything good, so raise the threshold.
                lookback_threshold *= lookback_scale_factor;
            }
            
            // Now it's safe to make a distance query
#ifdef debug_chaining
            cerr << "\t\tCome from score " << best_chain_score[*predecessor_index_it]
                << " across " << source << " to " << here << endl;
#endif
            
            // We will actually evaluate the source.
            
            // How far do we go in the graph?
            size_t graph_distance = get_graph_distance(source, here, distance_index, graph);
            
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
                TracedScore source_score = TracedScore::score_from(best_chain_score, *predecessor_index_it);
                
                // And the score with the transition and the points from the item
                TracedScore from_source_score = source_score.add_points(jump_points + item_points);
                
                // Remember that we could make this jump
                best_chain_score[i] = std::max(best_chain_score[i],
                                               from_source_score);
                                               
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
        cerr << "\tBest way to reach #" << i << " is " << best_chain_score[i] << endl;
#endif
        
        std::stringstream label_stream;
        label_stream << "#" << i << " " << here << " = " << item_points << "/" << best_chain_score[i].score;
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
        best_score.max_in(best_chain_score, i);
        
#ifdef debug_chaining
        cerr << "\tBest chain end so far: " << best_score << endl;
#endif
        
    }
    
    return best_score;
}

vector<size_t> chain_items_traceback(const vector<TracedScore>& best_chain_score,
                                     const VectorView<Anchor>& to_chain,
                                     const TracedScore& best_past_ending_score_ever) {
    
    // Now we need to trace back.
    vector<size_t> traceback;
    size_t here = best_past_ending_score_ever.source;
    if (here != TracedScore::nowhere()) {
#ifdef debug_chaining
        cerr << "Chain ends at #" << here << " " << to_chain[here]
            << " with score " << best_past_ending_score_ever << endl;
#endif
        while(here != TracedScore::nowhere()) {
            traceback.push_back(here);
#ifdef debug_chaining
            cerr << "Which gets score " << best_chain_score[here] << endl;
#endif
            here = best_chain_score[here].source;
#ifdef debug_chaining
            if (here != TracedScore::nowhere()) {
                cerr << "And comes after #" << here
                << " " << to_chain[here] << endl;
            } else {
                cerr << "And is first" << endl;
            }
#endif
        }
        // Flip it around front-ways
        std::reverse(traceback.begin(), traceback.end());
    }
    
#ifdef debug_chaining
    cerr << "Best score of chain overall: " << best_past_ending_score_ever << endl;
#endif

    return traceback;
}

pair<int, vector<size_t>> find_best_chain(const VectorView<Anchor>& to_chain,
                                          const SnarlDistanceIndex& distance_index,
                                          const HandleGraph& graph,
                                          int gap_open,
                                          int gap_extension,
                                          size_t max_lookback_bases,
                                          size_t initial_lookback_threshold,
                                          double lookback_scale_factor,
                                          double min_good_transition_score_per_base,
                                          int item_bonus,
                                          size_t max_indel_bases) {
                                                                 
    if (to_chain.empty()) {
        return std::make_pair(0, vector<size_t>());
    } else {
        
        // We actually need to do DP
        vector<TracedScore> best_chain_score;
        TracedScore best_past_ending_score_ever = chain_items_dp(best_chain_score,
                                                                 to_chain,
                                                                 distance_index,
                                                                 graph,
                                                                 gap_open,
                                                                 gap_extension,
                                                                 max_lookback_bases,
                                                                 initial_lookback_threshold,
                                                                 lookback_scale_factor,
                                                                 min_good_transition_score_per_base,
                                                                 item_bonus,
                                                                 max_indel_bases);
        // Then do the traceback and pair it up with the score.
        return std::make_pair(
            best_past_ending_score_ever.score,
            chain_items_traceback(best_chain_score, to_chain, best_past_ending_score_ever));
    }
}

int score_best_chain(const VectorView<Anchor>& to_chain, const SnarlDistanceIndex& distance_index, const HandleGraph& graph, int gap_open, int gap_extension) {
    
    if (to_chain.empty()) {
        return 0;
    } else {
        // Do the DP but without the traceback.
        vector<TracedScore> best_chain_score;
        TracedScore winner = algorithms::chain_items_dp(best_chain_score, to_chain, distance_index, graph, gap_open, gap_extension);
        return winner.score;
    }
}

size_t get_graph_distance(const Anchor& from, const Anchor& to, const SnarlDistanceIndex& distance_index, const HandleGraph& graph) {
    // TODO: hide something in the Anchors so we can use the minimizer cache information
    // For now just measure between the graph positions.
    
    auto from_pos = from.graph_end();
    auto& to_pos = to.graph_start();
    
    return distance_index.minimum_distance(
        id(from_pos), is_rev(from_pos), offset(from_pos),
        id(to_pos), is_rev(to_pos), offset(to_pos),
        false, &graph);  
}

size_t get_read_distance(const Anchor& from, const Anchor& to) {
    if (to.read_start() < from.read_end()) {
        return std::numeric_limits<size_t>::max();
    }
    return to.read_start() - from.read_end();
}

}
}
