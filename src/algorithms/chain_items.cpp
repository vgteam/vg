/**
 * \file chain_items.cpp
 * Non-template function implementations for chaining pieces of a read-to-graph alignment.
 */


#include "chain_items.hpp"
#include "crash.hpp"

#include <handlegraph/algorithms/dijkstra.hpp>
#include <structures/immutable_list.hpp>
#include <structures/min_max_heap.hpp>

#define debug_chaining
#define debug_transition

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
    out << "{R:" << anchor.read_start() << "=G:" << anchor.graph_start() << "(+" << anchor.start_hint_offset() << ")-"  << anchor.graph_end() << "(-" << anchor.end_hint_offset() << ")*" << anchor.length() << "}";
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
        vector<size_t> read_end_order = sort_permutation(to_chain.begin(), to_chain.end(), [&](const Anchor& a, const Anchor& b) {
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
            cerr << "\tFirst item overlapping #" << i << " beginning at " << here.read_start() << " is #" << *first_overlapping_it << " past-ending at " << to_chain[*first_overlapping_it].read_end() << " so start before there." << std::endl;
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
                    cerr << "\t\tDisregard due to read distance " << read_distance << " over limit " << max_lookback_bases << endl;
#endif
                        break;
                    } 
                }
                
                // Now it's safe to make a distance query
                
                // How far do we go in the graph? Don't bother finding out exactly if it is too much longer than in the read.
                size_t graph_distance = get_graph_distance(source, here, distance_index, graph, read_distance + max_indel_bases);
                
                std::pair<int, int> scores = {std::numeric_limits<int>::min(), std::numeric_limits<int>::min()};
                if (read_distance != numeric_limits<size_t>::max() && graph_distance != numeric_limits<size_t>::max()) {
                    // Transition seems possible, so yield it.
                    callback(*predecessor_index_it, i, read_distance, graph_distance);
                }
            } 
        }
    };

    return iterator;
}

transition_iterator zip_tree_transition_iterator(const std::vector<SnarlDistanceIndexClusterer::Seed>& seeds, const ZipCodeTree& zip_code_tree, size_t max_lookback_bases) {
    // TODO: Remove seeds because we only bring it here for debugging and it complicates the dependency relationships
    return [&seeds, &zip_code_tree, max_lookback_bases](const VectorView<Anchor>& to_chain,
                                                        const SnarlDistanceIndex& distance_index,
                                                        const HandleGraph& graph,
                                                        size_t max_indel_bases,
                                                        const transition_iteratee& callback) {
                            
        // We need a way to map from the seeds that zip tree thinks about to the anchors that we think about. So we need to index the anchors by leading/trailing seed.
        // TODO: Should we make someone else do the indexing so we can make the Anchor not need to remember the seed?
        std::unordered_map<size_t, size_t> seed_to_starting;
        std::unordered_map<size_t, size_t> seed_to_ending;
        for (size_t anchor_num = 0; anchor_num < to_chain.size(); anchor_num++) {
            seed_to_starting[to_chain[anchor_num].seed_start()] = anchor_num;
            seed_to_ending[to_chain[anchor_num].seed_end()] = anchor_num;
        }

        // Emit a transition between a source and destination anchor, or skip if actually unreachable.
        auto handle_transition = [&](size_t source_anchor_index, size_t dest_anchor_index, size_t graph_distance) {
            
            auto& source_anchor = to_chain[source_anchor_index];
            auto& dest_anchor = to_chain[dest_anchor_index];

#ifdef debug_transition
            std::cerr << "Handle transition #" << source_anchor_index << " " << source_anchor << " to #" << dest_anchor_index << " " << dest_anchor << std::endl;
#endif

            if (graph_distance == std::numeric_limits<size_t>::max()) {
                // Not reachable in graph (somehow)
                // TODO: Should never happen!
#ifdef debug_transition
                std::cerr << "\tNot reachable in graph!" << std::endl;
#endif
                return;
            }

            size_t read_distance = get_read_distance(source_anchor, dest_anchor);
            if (read_distance == std::numeric_limits<size_t>::max()) {
                // Not reachable in read
#ifdef debug_transition
                std::cerr << "\tNot reachable in read." << std::endl;
#endif
                return;
            }

            if (source_anchor.read_exclusion_end() > dest_anchor.read_exclusion_start()) {
                // The actual core anchor part is reachable in the read, but we cut these down from overlapping minimizers.
#ifdef debug_transition
                std::cerr << "\tOriginally overlapped in read." << std::endl;
#endif
                return;
            }

            // The zipcode tree is about point positions, but we need distances between whole anchors.
            // The stored zipcode positions will be at distances from the start/end of the associated anchor.
            
            // If the offset between the zip code point and the start of the destination is 0, and between the zip code point and the end of the source is 0, we subtract 0 from the measured distance. Otherwise we need to subtract something.
            size_t distance_to_remove = dest_anchor.start_hint_offset() + source_anchor.end_hint_offset();

#ifdef debug_transition
            std::cerr << "\tZip code tree sees " << graph_distance << " but we should back out " << distance_to_remove << std::endl;
#endif

            if (distance_to_remove > graph_distance) {
                // We actually end further along the graph path to the next
                // thing than where the next thing starts, so we can't actually
                // get there.
                return;
            }
            // Consume the length. 
            graph_distance -= distance_to_remove;

#ifdef debug_transition
            std::cerr << "\tZip code tree sees " << source_anchor << " and " << dest_anchor << " as " << graph_distance << " apart" << std::endl;
#endif

#ifdef double_check_distances

            auto from_pos = source_anchor.graph_end();
            auto to_pos = dest_anchor.graph_start();
            size_t check_distance = distance_index.minimum_distance(
                id(from_pos), is_rev(from_pos), offset(from_pos),
                id(to_pos), is_rev(to_pos), offset(to_pos),
                false, &graph);
            if (check_distance != graph_distance) {
                #pragma omp critical (cerr)
                std::cerr << "\tZip code tree sees " << source_anchor << " and " << dest_anchor << " as " << graph_distance << " apart but they are actually " << check_distance << " apart" << std::endl;
                crash_unless(check_distance == graph_distance);
            }

#endif

            // Send it along.
            callback(source_anchor_index, dest_anchor_index, read_distance, graph_distance); 
        };

        // If we find we are actually walking through the graph in opposition
        // to the read, we need to defer transitions from source on the read
        // forward strand to dest on the read forward strand, so we can go them
        // in order along the read forward strand.
        // This holds source, dest, and graph distance.
        // We will fill it all in and then sort it by destination read position.
        std::vector<std::tuple<size_t, size_t, size_t>> all_transitions;

        for (ZipCodeTree::iterator dest = zip_code_tree.begin(); dest != zip_code_tree.end(); ++dest) {
            // For each destination seed left to right
            ZipCodeTree::oriented_seed_t dest_seed = *dest;



            // Might be the start of an anchor if it is forward relative to the read, or the end of an anchor if it is reverse relative to the read
            std::unordered_map<size_t, size_t>::iterator found_dest_anchor = dest_seed.is_reverse ? seed_to_ending.find(dest_seed.seed) : seed_to_starting.find(dest_seed.seed);

            if (found_dest_anchor == (dest_seed.is_reverse ? seed_to_ending.end() : seed_to_starting.end())) {
                // We didn't find an anchor for this seed, maybe it lives in a different cluster. Skip it.
                continue;
            }

#ifdef debug_transition
            std::cerr << "Destination seed S" << dest_seed.seed << " " << seeds[dest_seed.seed].pos << (dest_seed.is_reverse ? "rev" : "") << " is anchor #" << found_dest_anchor->second << std::endl;
#endif

            for (ZipCodeTree::reverse_iterator source = zip_code_tree.look_back(dest, max_lookback_bases); source != zip_code_tree.rend(); ++source) {
                // For each source seed right to left
                ZipCodeTree::seed_result_t source_seed = *source;

                if (!source_seed.is_reverse && !dest_seed.is_reverse) {
                    // Both of these are in the same orientation relative to
                    // the read, and we're going through the graph in the
                    // read's forward orientation as assigned by these seeds.
                    // So we can just visit this transition.

                    // They might not be at anchor borders though, so check.
                    auto found_source_anchor = seed_to_ending.find(source_seed.seed);
                    if (found_source_anchor != seed_to_ending.end()) {
                        // We can transition between these seeds without jumping to/from the middle of an anchor.
#ifdef debug_transition
                        std::cerr << "\tSource seed S" << source_seed.seed << " " << seeds[source_seed.seed].pos << (source_seed.is_reverse ? "rev" : "") << " at distance " << source_seed.distance << "/" << max_lookback_bases << " is anchor #" << found_source_anchor->second << std::endl;
                        std::cerr << "\t\tFound transition from #" << found_source_anchor->second << " to #" << found_dest_anchor->second << std::endl;
#endif
                        all_transitions.emplace_back(found_source_anchor->second, found_dest_anchor->second, source_seed.distance);
                    } 
                } else if (source_seed.is_reverse && dest_seed.is_reverse) {
                    // Both of these are in the same orientation but it is opposite to the read.
                    // We need to find source as an anchor *started*, and thensave them flipped for later.
                    auto found_source_anchor = seed_to_starting.find(source_seed.seed);
                    if (found_source_anchor != seed_to_starting.end()) {
                        // We can transition between these seeds without jumping to/from the middle of an anchor.
                        // Queue them up, flipped
                        
#ifdef debug_transition
                        std::cerr << "\tSource seed S" << source_seed.seed << " " << seeds[source_seed.seed].pos << (source_seed.is_reverse ? "rev" : "") << " at distance " << source_seed.distance << "/" << max_lookback_bases << " is anchor #" << found_source_anchor->second << std::endl;
                        std::cerr << "\t\tFound backward transition from #" << found_dest_anchor->second << " to #" << found_source_anchor->second << std::endl;
#endif

                        all_transitions.emplace_back(found_dest_anchor->second, found_source_anchor->second, source_seed.distance);
                    } 
                } else {
                    // We have a transition between different orientations relative to the read. Don't show that.
                    continue;
                }
            }
        }

        // Sort the transitions so we handle them in akl allowed order for dynamic programming.
        std::sort(all_transitions.begin(), all_transitions.end(), [&](const std::tuple<size_t, size_t, size_t>& a, const std::tuple<size_t, size_t, size_t>& b) {
            // Return true if a's destination seed is before b's in the read, and false otherwise.
            return to_chain[get<1>(a)].read_start() < to_chain[get<1>(b)].read_start();
        });

        for (auto& transition : all_transitions) {
            // And handle all of them.
            // TODO: Inline this now-useless lambda that we call once.
            handle_transition(std::get<0>(transition), std::get<1>(transition), std::get<2>(transition));
        }
    };
}

/// Compute a gap score like minimap2.
///
/// They say they use the average anchor length, but really we need to use the
/// minimizer/base seed length here. Otherwise gaps cost more as your fragments
/// that you are chaining get longer, and cost more at chaining than at
/// fragmenting.
int score_chain_gap(size_t distance_difference, size_t base_seed_length) {
    if (distance_difference == 0) {
        return 0;
    } else {
        return 0.01 * base_seed_length * distance_difference + 0.5 * log2(distance_difference);
    }
}

TracedScore chain_items_dp(vector<TracedScore>& chain_scores,
                           const VectorView<Anchor>& to_chain,
                           const SnarlDistanceIndex& distance_index,
                           const HandleGraph& graph,
                           int gap_open,
                           int gap_extension,
                           const transition_iterator& for_each_transition,
                           int item_bonus,
                           int item_scale,
                           double gap_scale,
                           size_t max_indel_bases,
                           bool show_work) {
    
#ifdef debug_chaining
    DiagramExplainer diagram(show_work);
#else
    DiagramExplainer diagram(false);
#endif
    diagram.add_globals({{"rankdir", "LR"}});
   
#ifdef debug_chaining
    show_work = true;
#endif

    if (show_work) {
        cerr << "Chaining group of " << to_chain.size() << " items" << endl;
    }

    // Compute a base seed average length.
    // TODO: Weight anchors differently?
    // TODO: Will this always be the same for all anchors in practice?
    size_t base_seed_length = 0;
    for (auto& anchor : to_chain) {
        base_seed_length += anchor.base_seed_length();
    }
    base_seed_length /= to_chain.size();

    chain_scores.resize(to_chain.size());
    for (size_t i = 0; i < to_chain.size(); i++) {
        // Set up DP table so we can start anywhere with that item's score.
        chain_scores[i] = {to_chain[i].score(), TracedScore::nowhere()};
    }

    // We will run this over every transition in a good DP order.
    auto iteratee = [&](size_t from_anchor, size_t to_anchor, size_t read_distance, size_t graph_distance) {
        
        crash_unless(chain_scores.size() > to_anchor);
        crash_unless(chain_scores.size() > from_anchor);
        
        // For each item
        auto& here = to_chain[to_anchor];
        
        // How many points is it worth to collect?
        auto item_points = here.score() * item_scale + item_bonus;
        
        std::string here_gvnode = "i" + std::to_string(to_anchor);
        
        // If we come from nowhere, we get those points.
        chain_scores[to_anchor] = std::max(chain_scores[to_anchor], {item_points, TracedScore::nowhere()});
        
        // For each source we could come from
        auto& source = to_chain[from_anchor];
            
        if (show_work) {
            cerr << "\t\tCome from score " << chain_scores[from_anchor]
                << " across " << source << " to " << here << endl;
        }
            
        // How much does it pay (+) or cost (-) to make the jump from there
        // to here?
        // Don't allow the transition if it seems like we're going the long
        // way around an inversion and needing a huge indel.
        int jump_points;
            
        // Decide how much length changed
        size_t indel_length = (read_distance > graph_distance) ? read_distance - graph_distance : graph_distance - read_distance;
        
        if (show_work) {
            cerr << "\t\t\tFor read distance " << read_distance << " and graph distance " << graph_distance << " an indel of length " << indel_length << " would be required" << endl;
        }

        if (indel_length > max_indel_bases) {
            // Don't allow an indel this long
            jump_points = std::numeric_limits<int>::min();
        } else {
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
            jump_points = -score_chain_gap(indel_length, base_seed_length) * gap_scale;
        }
            
        if (jump_points != numeric_limits<int>::min()) {
            // Get the score we are coming from
            TracedScore source_score = TracedScore::score_from(chain_scores, from_anchor);
            
            // And the score with the transition and the points from the item
            TracedScore from_source_score = source_score.add_points(jump_points + item_points);
            
            // Remember that we could make this jump
            chain_scores[to_anchor] = std::max(chain_scores[to_anchor], from_source_score);
                                           
            if (show_work) {
                cerr << "\t\tWe can reach #" << to_anchor << " with " << source_score << " + " << jump_points << " from transition + " << item_points << " from item = " << from_source_score << endl;
            }
                    
            if (from_source_score.score > 0) {
                // Only explain edges that were actual candidates since we
                // won't let local score go negative
                
                std::string source_gvnode = "i" + std::to_string(from_anchor);
                // Suggest that we have an edge, where the edges that are the best routes here are the most likely to actually show up.
                diagram.suggest_edge(source_gvnode, here_gvnode, here_gvnode, from_source_score.score, {
                    {"label", std::to_string(jump_points)},
                    {"weight", std::to_string(std::max<int>(1, from_source_score.score))}
                });
            }
        } else {
            if (show_work) {
                cerr << "\t\tTransition is impossible." << endl;
            }
        }
    };

    // Run our DP step over all the transitions.
    for_each_transition(to_chain,
                        distance_index,
                        graph,
                        max_indel_bases,
                        iteratee);
        
   
    TracedScore best_score = TracedScore::unset();

    for (size_t to_anchor = 0; to_anchor < to_chain.size(); ++to_anchor) {
        // For each destination anchor, now that it is finished, see if it is the winner.
        auto& here = to_chain[to_anchor];
        auto item_points = here.score() * item_scale + item_bonus;

        if (show_work) {
            cerr << "\tBest way to reach #" << to_anchor  << " " << to_chain[to_anchor] << " is " << chain_scores[to_anchor] << endl;
        }
        
        // Draw the item in the diagram
        std::string here_gvnode = "i" + std::to_string(to_anchor);
        std::stringstream label_stream;
        label_stream << "#" << to_anchor << " " << here << " = " << item_points << "/" << chain_scores[to_anchor].score;
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
        best_score.max_in(chain_scores, to_anchor);
        
        if (show_work) {
            cerr << "\tBest chain end so far: " << best_score << endl;
        }
        
    }
    
    return best_score;
}

vector<pair<vector<size_t>, int>> chain_items_traceback(const vector<TracedScore>& chain_scores,
                                                        const VectorView<Anchor>& to_chain,
                                                        const TracedScore& best_past_ending_score_ever,
                                                        int item_bonus,
                                                        int item_scale,
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
                    penalty -= (to_chain[here].score() * item_scale + item_bonus);
                    // TODO: Score this more simply.
                    // TODO: find the edge to nowhere???
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
                                                   const transition_iterator& for_each_transition, 
                                                   int item_bonus,
                                                   int item_scale,
                                                   double gap_scale,
                                                   size_t max_indel_bases,
                                                   bool show_work) {
                                                                         
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
                                                             for_each_transition,
                                                             item_bonus,
                                                             item_scale,
                                                             gap_scale,
                                                             max_indel_bases,
                                                             show_work);
    // Then do the tracebacks
    vector<pair<vector<size_t>, int>> tracebacks = chain_items_traceback(chain_scores, to_chain, best_past_ending_score_ever, item_bonus, item_scale, max_chains);
    
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
                                          const transition_iterator& for_each_transition,
                                          int item_bonus,
                                          int item_scale,
                                          double gap_scale,
                                          size_t max_indel_bases) {
                                                                 
    return find_best_chains(
        to_chain,
        distance_index,
        graph,
        gap_open,
        gap_extension,
        1,
        for_each_transition,
        item_bonus,
        item_scale,
        gap_scale,
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
//#define double_check_distances
//#define stop_on_mismatch
//#define replace_on_mismatch
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
