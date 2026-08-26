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
        return out << value.score << " (eval +" << value.eval_bonus << ") from nowhere";
    }
    return out << value.score << " (eval +" << value.eval_bonus << ") from #" << value.source;
}


void TracedScore::max_in(const vector<vector<TracedScore>>& options, size_t option_number) {
    auto& option = options[option_number].front();
    if (option.score > this->score || this->source == nowhere()) {
        // This is the new winner.
        this->score = option.score;
        this->eval_bonus = option.eval_bonus;
        this->source = option_number;
        this->paths = option.paths;
        this->rec_num = option.rec_num;
    }
}

TracedScore TracedScore::score_from(const vector<vector<TracedScore>>& options, size_t option_number) {
    TracedScore got = options[option_number].front();
    got.source = option_number;
    return got;
}

TracedScore TracedScore::add_points(int adjustment) const {
    return {this->score + adjustment, this->eval_bonus, this->source, this->paths, this->rec_num};
}

void TracedScore::try_to_insert_self(vector<TracedScore>& current) const {
    // Figure out which index to insert this element at
    int insert_after_index = current.size() - 1;
    while (insert_after_index >= 0) {
        if (current[insert_after_index] > *this) {
            // This item beats me, so we need to insert after it
            break;
        }
        // We beat this item, so shift things over
        if (insert_after_index > 0) {
            current[insert_after_index] = current[insert_after_index-1];
        }
        insert_after_index--;
    }

    if (insert_after_index != current.size() - 1) {
        current[insert_after_index+1] = *this;
    }
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
        this->eval_bonus,
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
                                          max_read_lookback_bases, max_indel_bases, to_chain,
                                          seed_to_starting, seed_to_ending);

        // Sort the transitions so we handle them in an allowed order for dynamic programming.
        std::sort(all_transitions.begin(), all_transitions.end(), 
            [&](const transition_info& a, const transition_info& b) {
            // Return true if a's destination seed is before b's in the read, and false otherwise.
            // Break ties by indel size (smaller = better)
            return to_chain[a.to_anchor].read_start() < to_chain[b.to_anchor].read_start()
                   || (to_chain[a.to_anchor].read_start() == to_chain[b.to_anchor].read_start()
                       && a.indel_size < b.indel_size);
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
    size_t max_indel_bases,
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
                                                max_indel_bases, found_source_anchor->second,
                                                cur_dest_anchor.second, source_seed.distance);
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
                                                max_indel_bases, cur_dest_anchor.second,
                                                found_source_anchor->second, source_seed.distance);
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

void add_transition_if_legal(vector<transition_info>& transitions, const VectorView<Anchor>& to_chain,
                             size_t max_read_lookback_bases, size_t max_indel_bases,
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

    size_t indel_size = (read_distance > graph_distance) ? read_distance - graph_distance 
                                                         : graph_distance - read_distance;

    if (indel_size > max_indel_bases) {
#ifdef debug_transition
        std::cerr << "\tIndel size " << indel_size << " over max of " << max_indel_bases << std::endl;
#endif 
        return;
    }
#ifdef debug_transition
    std::cerr << "\tZip code tree sees " << source_anchor << " and "
              << dest_anchor << " as " << graph_distance << " apart; "
              << indel_size << " indel" << std::endl;
#endif
    transitions.emplace_back(from_anchor, to_anchor, indel_size);
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

/// If the current anchor does not share paths with the chain, pay a penalty.
/// Returns a positive value (gap penalty)
int check_recombination(const TracedScore& from, const Anchor& to) {
    if ((from.paths & to.anchor_start_paths()) == 0) {
        return 1;
    } else {
        return 0;
    }
}

void chain_items_dp(vector<vector<TracedScore>>& chain_scores,
                    const VectorView<Anchor>& to_chain,
                    const SnarlDistanceIndex& distance_index,
                    const HandleGraph& graph,
                    const transition_iterator& for_each_transition,
                    size_t max_predecessors,
                    const ChainScoringScheme& scheme,
                    size_t max_indel_bases,
                    bool show_work) {

#ifdef debug_chaining
    show_work = true;
#endif

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


    // Compute a base seed average length.
    // TODO: Weight anchors differently?
    // TODO: Will this always be the same for all anchors in practice?
    size_t base_seed_length = 0;
    for (auto& anchor : to_chain) {
        base_seed_length += anchor.base_seed_length();
    }
    base_seed_length /= to_chain.size();

    chain_scores.resize(to_chain.size());

    // We want to prefer to come from seeds where the transition preserves
    // access to matching haplotypes, because we don't want to back ourselves
    // into a corner where we need a recombination when we don't really have
    // to. So we cheat on the dynamic programming by adding an "evaluation
    // bonus" to the scores of the different DP options when comparing them. We
    // keep this bonus out of the actual recorded scores because we don't want
    // it raising the scores we actually get the more transitions we take.
    //
    // We store the bonus used to select the current winning predecessor for
    // each seed in this vector, which runs alongside the DP table.
    //
    // Starting from nowhere means full path conservation, so bonus = scheme.consistency_bonus.
    for (size_t i = 0; i < to_chain.size(); i++) {
        // Set up DP table so we can start anywhere with that item's score, with bonus applied.
        chain_scores[i] = std::vector<TracedScore>(max_predecessors,
            {to_chain[i].score() + scheme.item_bonus, scheme.consistency_bonus,
             TracedScore::nowhere(), to_chain[i].anchor_end_paths()});
    }

    // We will run this over every transition in a good DP order.
    auto iteratee = [&](const transition_info& transition) {
        if (show_work) {
#ifdef debug_dp
            cerr << "DP: " << transition.from_anchor << "->" << transition.to_anchor 
                 << " indel " << transition.indel_size << endl;
#endif
        }
        
        crash_unless(chain_scores.size() > transition.to_anchor);
        crash_unless(chain_scores.size() > transition.from_anchor);
        
        // For each item
        auto& here = to_chain[transition.to_anchor];
        
        // How many points is it worth to collect?
        auto item_points = here.score() + scheme.item_bonus;
        
        std::string here_gvnode;
        if (diagram) {
            here_gvnode = "i" + std::to_string(transition.to_anchor);
        }
        
        // For each source we could come from
        auto& source = to_chain[transition.from_anchor];
            
        if (show_work) {
#ifdef debug_dp
            cerr << "\t\tCome from score " << chain_scores[transition.from_anchor].front()
                 << " across " << source << " to " << here << endl;
#endif
        }

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
        int jump_points = -score_chain_gap(transition.indel_size, base_seed_length) * scheme.gap_scale;

        // add recombination penalty if necessary
        jump_points -= check_recombination(chain_scores[transition.from_anchor].front(), here) * scheme.recombination_penalty;
            
        // Get the score we are coming from
        TracedScore source_score = TracedScore::score_from(chain_scores, transition.from_anchor);
        
        // And the score with the transition and the points from the item
        TracedScore from_source_score = source_score.add_points(jump_points + item_points)
                                                    .set_shared_paths(here.anchor_paths());
        
        // Evaluate heuristic to preserve path flexibility without inflating actual scoring DP.
        // Bonus = fraction of conserved paths * scheme.consistency_bonus.
        // Bonus is 0 when recombination occurs (no shared paths).
        from_source_score.eval_bonus = 0;
        if (scheme.consistency_bonus > 0) {
            int pre_count = __builtin_popcountll(source_score.paths);
            if (pre_count > 0 && (source_score.paths & here.anchor_start_paths()) != 0) {
                // No recombination: bonus = fraction of paths conserved * penalty
                int post_count = __builtin_popcountll(from_source_score.paths);
                from_source_score.eval_bonus = (scheme.consistency_bonus * post_count) / pre_count;
            }
            // Recombination case (no shared paths): bonus stays 0
        }

        from_source_score.try_to_insert_self(chain_scores[transition.to_anchor]);
                                        
        if (show_work) {
#ifdef debug_dp
            cerr << "\t\tWe can reach #" << transition.to_anchor << " with " << source_score << " + " << jump_points
                 << " from transition + " << item_points << " from item = " << from_source_score << endl;
#endif
        }
        
        if (diagram) {
            if (from_source_score.score > 0) {
                // Only explain edges that were actual candidates since we
                // won't let local score go negative
                
                std::string source_gvnode = "i" + std::to_string(transition.from_anchor);
                // Suggest that we have an edge, where the edges that are
                // the best routes here are the most likely to actually show up.
                diagram.suggest_edge(source_gvnode, here_gvnode, here_gvnode, from_source_score.score, {
                    {"label", std::to_string(jump_points)},
                    {"weight", std::to_string(std::max<int>(1, from_source_score.score))}
                });
            }
        }
        if (dump) {
            // Dump a TSV of source anchor description, dest anchor description, and score achieved.
            dump.line();
            std::stringstream source_stream;
            source_stream << source;
            dump.field(source_stream.str());
            std::stringstream here_stream;
            here_stream << here;
            dump.field(here_stream.str());
            dump.field(from_source_score.score);
        }
    };

    // Run our DP step over all the transitions.
    for_each_transition(to_chain,
                        distance_index,
                        graph,
                        max_indel_bases,
                        iteratee);
        
    if (show_work) {
        TracedScore best_score = TracedScore::unset();

        for (size_t to_anchor = 0; to_anchor < to_chain.size(); ++to_anchor) {
            // For each destination anchor, now that it is finished, see if it is the winner.
            auto& here = to_chain[to_anchor];

            cerr << "\tBest way to reach #" << to_anchor  << " " << to_chain[to_anchor]
                 << " is " << chain_scores[to_anchor].front() << "\n"
                 << "\t\tbut you can also do: ";
            for (size_t alt_i = 1; alt_i < chain_scores[to_anchor].size(); alt_i++) {
                cerr << chain_scores[to_anchor][alt_i] << " ";
            }
            cerr << endl;
            
            if (diagram) {
                // Draw the item in the diagram
                auto item_points = here.score() + scheme.item_bonus;
                std::string here_gvnode = "i" + std::to_string(to_anchor);
                std::stringstream label_stream;
                label_stream << "#" << to_anchor << " " << here << " = " << item_points
                            << "/" << chain_scores[to_anchor].front().score;
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
                std::string graph_gvnode2 = ("n" + std::to_string(id(graph_start) 
                                            + (is_rev(graph_start) ? -1 : 1))
                                            + (is_rev(graph_start) ? "r" : "f"));
                diagram.ensure_node(graph_gvnode2, {
                    {"label", std::to_string(id(graph_start) + (is_rev(graph_start) ? -1 : 1)) 
                                            + (is_rev(graph_start) ? "-" : "+")},
                    {"shape", "box"}
                });
                // And show them as connected. 
                diagram.ensure_edge(graph_gvnode, graph_gvnode2, {{"color", "gray"}});
            }
            
            // See if this is the best overall
            best_score.max_in(chain_scores, to_anchor);

#ifdef debug_dp
            cerr << "\tBest chain end so far: " << best_score << endl; 
#endif
        }
#ifdef debug_chaining
        std::cerr << "[REC INFO] Recombination number for chain: "
                  << best_score.rec_num << "\tscore: "
                  << best_score.score << "\tpaths: "
                  << best_score.paths << std::endl;
#endif
    }
}

size_t count_total_recombinations(const SubchainGroup& group, const vector<size_t>& subchains_used) {
    size_t total_rec_count = group.subchains[subchains_used.front()].rec_count;
    path_flags_t cur_paths = group.subchains[subchains_used.front()].end_paths;
    for (size_t i = 1; i < subchains_used.size(); i++) {
        const Subchain& subchain = group.subchains[subchains_used[i]];
        if ((cur_paths & subchain.start_paths) == 0) {
            // This inter-subchain connection requires a recombination
            total_rec_count += 1;
            // Paths reset to whatever this next subchain wants
            cur_paths = subchain.end_paths;
        } else if (subchain.rec_count > 0) {
            // The subchain itself requires a recombination, so paths reset
            cur_paths = subchain.end_paths;
        } else {
            // No forced recombinations, so paths narrow down to what is OK with the end
            cur_paths &= subchain.end_paths;
        }
        total_rec_count += subchain.rec_count;
    }

    return total_rec_count;
}

void chain_items_traceback(const vector<vector<TracedScore>>& chain_scores,
                           const VectorView<Anchor>& to_chain,
                           vector<SparseAnchorChain>& tracebacks,
                           vector<AltEdge>& connections,
                           const ChainScoringScheme& scheme,
                           size_t max_tracebacks) {
    tracebacks.reserve(chain_scores.size());
    connections.reserve(chain_scores.size());
    
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
    
    // We don't want to use an item multiple times
    vector<bool> is_used(chain_scores.size(), false);
    
    for (auto& trace_from : starts_in_score_order) {
        // Save alt edges
        for (size_t alt_i = 1; alt_i < chain_scores[trace_from].size(); alt_i++) {
            if (chain_scores[trace_from][alt_i].source != TracedScore::nowhere()) {
                connections.emplace_back(chain_scores[trace_from][alt_i].source,
                                        trace_from,
                                        chain_scores[trace_from].front().score - chain_scores[trace_from][alt_i].score);
            }
        }
        if (is_used[trace_from]) {
            // We can't trace back from here
            continue;
        }
        // For each unused item in score order, start a traceback stack (in reverse order)
        tracebacks.emplace_back();
        tracebacks.back().chain_score = chain_scores[trace_from].front().score;
        vector<size_t> cur_anchors = {trace_from};
        size_t here = trace_from;
#ifdef debug_chaining
        std::cerr << "[REC INFO] Starting traceback at item #" << here
                  << " with recs: " << chain_scores[here].front().rec_num
                  << " score: " << chain_scores[here].front().score << std::endl;
#endif
        while (here != TracedScore::nowhere()) {
#ifdef debug_chaining
            std::cerr << "\trecs: " << chain_scores[here].front().rec_num
                      << " paths:\t" << std::bitset<64>(chain_scores[here].front().paths) << std::endl;
            std::cerr << "\t\tanchor #" << here << ": " << to_chain[here] << std::endl;
#endif
            // Mark here as used. Happens once per item, and so limits runtime.
            is_used[here] = true;
            size_t next = chain_scores[here].front().source;
            if (next != TracedScore::nowhere()) {
                if (is_used[next]) {
                    // Save this extra edge we tried to use
                    connections.emplace_back(next, here, 0);

                    // We need to stop early and accrue an extra penalty.
                    // Take away all the points we got for coming from there and being ourselves.
                    tracebacks.back().chain_score -= chain_scores[here].front().score;
                    // But then re-add our score for just us
                    tracebacks.back().chain_score += (to_chain[here].score() + scheme.item_bonus);
                    // TODO: Score this more simply.
                    // TODO: find the edge to nowhere???
                    break;
                } else {
                    // Add to the traceback
                    cur_anchors.push_back(next);
                }
            }
            here = next;
        }

        // Make sure to order the steps left to right, and not right to left as we generated them.
        std::copy(cur_anchors.rbegin(), cur_anchors.rend(),
                  std::back_inserter(tracebacks.back().anchors));
    }
    
    // Sort the tracebacks by score
    std::sort(tracebacks.begin(), tracebacks.end(), 
    [](const SparseAnchorChain& a, const SparseAnchorChain& b) {
        // Return true if a has the larger score and belongs first
        return a.chain_score > b.chain_score;
    });
    
    if (tracebacks.size() > max_tracebacks) {
        // Limit to requested number
        tracebacks.resize(max_tracebacks);
    }
}

/// TODO: break this function up it's getting really big
vector<SubchainGroup> split_up_subchains(const VectorView<Anchor>& to_chain,
                                         const vector<SparseAnchorChain>& original_tracebacks,
                                         const vector<AltEdge>& connections) {
    // For each anchor (by index), which traceback was it originally in?
    vector<size_t> home_trace(to_chain.size(), numeric_limits<size_t>::max());
    // For each anchor (by index), which other anchors can it connect to?
    vector<unordered_set<size_t>> outgoing_edges(to_chain.size());
    vector<unordered_set<size_t>> incoming_edges(to_chain.size());

    // Where to search for subchains (we must guarantee they start in read order)
    priority_queue<size_t, vector<size_t>, std::greater<size_t>> trace_from;
    // Best tie-ins for each tail end
    unordered_map<TailAnchor, AltEdge> tail_edges;

    // Where we will save results
    vector<SubchainGroup> output;

    // Remember which anchors are used
    for (size_t i = 0; i < original_tracebacks.size(); i++) {
        const vector<size_t>& cur_anchors = original_tracebacks[i].anchors;
        // Placeholder tie-ins that we want to fix later
        tail_edges.emplace(original_tracebacks[i].left_tail(), AltEdge(TracedScore::nowhere(), cur_anchors.front()));
        tail_edges.emplace(original_tracebacks[i].right_tail(), AltEdge(cur_anchors.back(), TracedScore::nowhere()));
#ifdef debug_chaining
        cerr << "Chain: ";
#endif
        for (const auto& anchor : cur_anchors) {
#ifdef debug_chaining
            cerr << anchor << " ";
#endif
            home_trace[anchor] = i;
        }
#ifdef debug_chaining
        cerr << endl;
#endif
    }

    // Now check in on the extra edges.
    for (const auto& edge : connections) {
        if (home_trace[edge.start_anchor] != numeric_limits<size_t>::max() 
            && home_trace[edge.end_anchor] != numeric_limits<size_t>::max()
            && home_trace[edge.start_anchor] != home_trace[edge.end_anchor]) {
#ifdef debug_chaining
            cerr << "Extra edge " << edge.start_anchor << " --> " << edge.end_anchor << " score diff " << edge.score_diff << endl;
#endif
            TailAnchor right_tie_in(edge.start_anchor, false);
            TailAnchor left_tie_in(edge.end_anchor, true);
            if (tail_edges.count(right_tie_in) && tail_edges.at(right_tie_in) > edge) {
                // Best right tie-in so far
                tail_edges[right_tie_in] = edge;
            }
            if (tail_edges.count(left_tie_in) && tail_edges.at(left_tie_in) > edge) {
                // Best left tie-in so far
                tail_edges[left_tie_in] = edge;
            }
        }
    }

    // Which of the tracebacks link up?
    vector<bool> is_joined_up(original_tracebacks.size(), false);
    for (const auto& tie_in : tail_edges) {
        const AltEdge& edge = tie_in.second;
        if (!edge.is_max_score_diff()) {
            // Use this extra edge
            outgoing_edges[edge.start_anchor].emplace(edge.end_anchor);
            incoming_edges[edge.end_anchor].emplace(edge.start_anchor);
            // Remember that these tracebacks are connected
            is_joined_up[home_trace[edge.start_anchor]] = true;
            is_joined_up[home_trace[edge.end_anchor]] = true;
        }
    }

    // Save edges that are part of an original traceback
    int max_joined_chain_score = 0;
    for (size_t i = 0; i < original_tracebacks.size(); i++) {
        if (!is_joined_up[i]) {
#ifdef debug_chaining
            cerr << "Saving traceback " << i << " as its own SubchainGroup" << endl;
#endif
            // This traceback is disjoint and should be returned separately
            output.emplace_back();
            output.back().subchains.emplace_back(original_tracebacks[i].anchors, true);
            output.back().max_sparse_chain_score = original_tracebacks[i].chain_score;
        } else {
            // We'll tie this traceback in to the rest
            max_joined_chain_score = std::max(max_joined_chain_score, original_tracebacks[i].chain_score);
            const vector<size_t>& cur_anchors = original_tracebacks[i].anchors;
            trace_from.emplace(cur_anchors.front());
            for (size_t i = 0; i < cur_anchors.size() - 1; i++) {
                // Remember that this pair of anchors can connect
                incoming_edges[cur_anchors[i+1]].emplace(cur_anchors[i]);
                outgoing_edges[cur_anchors[i]].emplace(cur_anchors[i+1]);
            }
        }
    }

    if (output.size() == original_tracebacks.size()) {
        // Everything was disjoint; return right away
        return output;
    }

#ifdef debug_chaining
    for (size_t i = 0; i < to_chain.size(); i++) {
        cerr << i << " has " << incoming_edges[i].size() << " sources, and outgoing edges to ";
        for (const auto& next : outgoing_edges[i]) {
            cerr << next << " ";
        }
        cerr << endl;
    }
#endif

    // From now on we're building our final connected SubchainGroup
    output.emplace_back();

    // Set up the subchains
    // For each anchor (by index), which subchain did it end up in?
    vector<size_t> subchain_id(to_chain.size(), std::numeric_limits<size_t>::max());
    while (!trace_from.empty()) {
        // Start trace for this subchain, until we hit into an endpoint
        size_t cur_anchor_id = trace_from.top();
        trace_from.pop();
        if (subchain_id[cur_anchor_id] != std::numeric_limits<size_t>::max()) {
            // This one was already put as part of an earlier subchain
            continue;
        }

        // Create a new subchain to trace into
        size_t cur_subchain_id = output.back().subchains.size();
        output.back().subchains.emplace_back(vector<size_t>{cur_anchor_id});
        subchain_id[cur_anchor_id] = cur_subchain_id;

#ifdef debug_chaining
        cerr << "Assign " << cur_anchor_id << " to subchain " << cur_subchain_id << endl;
#endif

        while (true) { 
            if (outgoing_edges[cur_anchor_id].empty()) {
                // This anchor can't trace outwards at all
                break;
            }           
            // If we've reached a decision point,
            // or the end of a traceback, then save all next edges
            // We have reached the end of one subchain
            size_t next_anchor_id = *outgoing_edges[cur_anchor_id].begin();
            if (outgoing_edges[cur_anchor_id].size() > 1 
                || incoming_edges[next_anchor_id].size() > 1
                || tail_edges.count(TailAnchor(cur_anchor_id, false))
                || tail_edges.count(TailAnchor(next_anchor_id, true))) {
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
            cur_anchor_id = next_anchor_id;
            output.back().subchains.back().anchors.emplace_back(cur_anchor_id);
            subchain_id[cur_anchor_id] = cur_subchain_id;
#ifdef debug_chaining
            cerr << "Assign " << cur_anchor_id << " to subchain " << cur_subchain_id << endl;
#endif
        }
    }

    // Make all inter-subchain connections necessary
    for (size_t i = 0; i < output.back().subchains.size(); i++) {
        // Loop over any way to exit this subchain
        Subchain& cur_subchain = output.back().subchains[i];
        for (const auto& next : outgoing_edges[cur_subchain.anchors.back()]) {
            if (next != std::numeric_limits<size_t>::max()) {
                output.back().connections.emplace_back(i, subchain_id[next]);
#ifdef debug_chaining
                cerr << "Connect subchains " << i << " -> " << subchain_id[next] << endl;
#endif
            }
        }

        // Mark which ones should get tails
        if (tail_edges.count(cur_subchain.left_tail())) {
            cur_subchain.add_left_tail = true;
        }
        if (tail_edges.count(cur_subchain.right_tail())) {
            cur_subchain.add_right_tail = true;
        }
    }

    // Assume first is highest score
    output.back().max_sparse_chain_score = max_joined_chain_score;
    return output;
}

vector<SubchainGroup> find_best_chains(const VectorView<Anchor>& to_chain,
                                       const SnarlDistanceIndex& distance_index,
                                       const HandleGraph& graph,
                                       const transition_iterator& for_each_transition,
                                       const ChainScoringScheme& scheme,
                                       size_t max_chains,
                                       size_t max_indel_bases,
                                       size_t min_chain_score,
                                       bool show_work) {

    if (to_chain.empty()) {
        // Nothing to chain
        return {SubchainGroup()};
    }
        
    // We actually need to do DP
    vector<vector<TracedScore>> chain_scores;
    chain_items_dp(chain_scores, to_chain, distance_index, graph, for_each_transition,
                   3, /// TODO: make into a param
                   scheme, max_indel_bases, show_work);
    
    // Then do the tracebacks
    vector<SparseAnchorChain> tracebacks;
    vector<AltEdge> connections;
    chain_items_traceback(chain_scores, to_chain, tracebacks, connections, scheme, max_chains);
    
    if (tracebacks.empty()) {
        // Somehow we got nothing
        return {SubchainGroup()};
    }

    // Get rid of tracebacks that are much, much worse than the best
    for (size_t i = 1; i < tracebacks.size(); i++) {
        if (tracebacks[i].chain_score < min_chain_score) {
#ifdef debug_chaining
            cerr << "Cutting down to " << i << " tracebacks because a further one has score "
                 << tracebacks[i].chain_score << " < " << min_chain_score << endl; 
#endif
            // Cut off at this point
            tracebacks.resize(i);
            break;
        }
    }

    vector<SubchainGroup> subchain_groups = split_up_subchains(to_chain, tracebacks, connections);

    for (SubchainGroup& group : subchain_groups) {
        for (Subchain& subchain : group.subchains) {
            // Compute the anchor indices in this chain that introduce an
            // inter-anchor recombination event. We simulate the path-bit
            // propagation along the chain using the same logic as
            // TracedScore::set_shared_paths (but without modifying the state).
            // Start with the endpoint paths of the first anchor.
            const Anchor& first_anchor = to_chain[subchain.anchors.front()];
            // Paths at the start of the subchain are for its first anchor
            subchain.start_paths = first_anchor.anchor_start_paths();
            // Is the first anchor internally recombinant?
            subchain.rec_count = first_anchor.anchor_start_paths() != first_anchor.anchor_end_paths();

            // Walk the chain from the second anchor onward and apply the
            // same recombination-detection rules used in set_shared_paths.
            path_flags_t current_paths = first_anchor.anchor_end_paths();
            for (size_t i = 1; i < subchain.anchors.size(); ++i) {
                auto new_paths = to_chain[subchain.anchors[i]].anchor_paths();
                // If the anchor's start and end paths are equal, it's not an
                // internally recombinant anchor; check inter-anchor overlap.
                if (new_paths.first == new_paths.second) {
                    if ((current_paths & new_paths.first) == 0) {
                        // No overlap -> inter-anchor recombination occurred here.
                        subchain.rec_count++;
                        // Reset current paths to the anchor's start paths.
                        current_paths = new_paths.first;
                    } else {
                        // Intersect supported paths and continue.
                        current_paths &= new_paths.first;
                    }
                } else {
                    // Recombinant anchor: do not count as inter-anchor
                    // recombination per original logic; reset paths to the
                    // anchor's end paths.
                    // TODO: since no fragmenting, probably unnecessary to track
                    subchain.rec_count++;
                    current_paths = new_paths.second;
                }
            }
            subchain.end_paths = current_paths;
        }
    }

    return subchain_groups;
}

SparseAnchorChain find_best_chain(const VectorView<Anchor>& to_chain,
                                  const SnarlDistanceIndex& distance_index,
                                  const HandleGraph& graph,
                                  const transition_iterator& for_each_transition,
                                  const ChainScoringScheme& scheme,
                                  size_t max_indel_bases,
                                  size_t min_chain_score) {
    vector<SubchainGroup> groups = find_best_chains(to_chain,
                                                   distance_index,
                                                   graph,
                                                   for_each_transition,
                                                   scheme,
                                                   1, // Only one chain needed!
                                                   max_indel_bases,
                                                   min_chain_score);
    if (groups.empty() || groups.front().subchains.empty()) {
        // We got nothing
        return SparseAnchorChain();
    }

    SparseAnchorChain output;
    output.anchors = groups.front().subchains.front().anchors;
    output.chain_score = groups.front().max_sparse_chain_score;
    
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
