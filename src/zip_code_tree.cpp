//#define DEBUG_ZIP_CODE_TREE
//#define DEBUG_ZIP_CODE_SORTING

#include "zip_code_tree.hpp"
#include "crash.hpp"

// Set for verbose logging from the zip code tree parsing logic
//#define debug_parse

// Set to compile in assertions to check the zipcode tree parsing logic
//#define check_parse

using namespace std;
namespace vg {

void ZipCodeTree::print_self(const vector<Seed>* seeds) const {
    for (const tree_item_t item : zip_code_tree) {
        if (item.get_type() == SEED) {
            cerr << seeds->at(item.get_value()).pos;
            if (item.get_is_reversed()) {
                cerr << "rev";
            }
            if (item.has_other_values()) {
                cerr << "<" << item.get_value();
                for (const auto& other_value : item.get_other_values()) {
                    cerr << "/" << other_value;
                }
                cerr << ">";
            }
        } else if (item.get_type() == DAG_SNARL_START) {
            cerr << "(";
        } else if (item.get_type() == DAG_SNARL_END) {
            cerr << ")";
        } else if (item.get_type() == CYCLIC_SNARL_START) {
            cerr << "{";
        } else if (item.get_type() == CYCLIC_SNARL_END) {
            cerr << "}";
        } else if (item.get_type() == CHAIN_START) {
            cerr << "[";
        } else if (item.get_type() == CHAIN_END) {
            cerr << "]";
        } else if (item.get_type() == EDGE) {
            cerr << " ";
            if (item.get_value() == std::numeric_limits<size_t>::max()) {
                cerr << "inf";
            } else {
                cerr << item.get_value();
            }
            cerr << " ";
        } else if (item.get_type() == CHAIN_COUNT) {
            cerr << item.get_value() << " ";
        } else {
            throw std::runtime_error("[zip tree]: Trying to print a zip tree item of unsupported type: " 
                                     + std::to_string(item.get_type()));
        }
    }
    cerr << endl;
}

void ZipCodeTree::add_close_bound(size_t start_index) {
    tree_item_type_t closing_type;
    if (zip_code_tree[start_index].get_type() == DAG_SNARL_START) {
        closing_type = DAG_SNARL_END;
    } else if (zip_code_tree[start_index].get_type() == CYCLIC_SNARL_START) {
        closing_type = CYCLIC_SNARL_END;
    } else if (zip_code_tree[start_index].get_type() == CHAIN_START) {
        closing_type = CHAIN_END;
    } else {
        throw std::runtime_error("[zip tree]: Attempting to close a zip tree item that is not a snarl or chain start");
    }
    zip_code_tree.emplace_back(closing_type, zip_code_tree[start_index].get_value(), start_index);
    bound_pair_indexes[start_index] = zip_code_tree.size() - 1;
    bound_pair_indexes[zip_code_tree.size() - 1] = start_index;
}

void ZipCodeForest::open_chain(forest_growing_state_t& forest_state, const interval_state_t& interval) {
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "\t\tOpen new chain at depth " << interval.depth << endl;
    // Top-level chains should be opened elsewhere
    assert(interval.depth > 0);
#endif
    size_t first_seed_index = forest_state.seed_sort_order[interval.interval_start];
    const Seed& first_seed = forest_state.seeds->at(first_seed_index);
    // If the current depth equals this seed's depth, we are at a node
    bool first_is_node = first_seed.zipcode.max_depth() == interval.depth;

    // Now record the start of this chain
    // Look up parent snarl ID; open_chain() is never used for top-level chains
    size_t snarl_id = forest_state.sibling_indices_at_depth[interval.depth-1][0].snarl_id;
    trees[forest_state.active_tree_index].zip_code_tree.emplace_back(ZipCodeTree::CHAIN_START, snarl_id);

    // Remember the chain start and its prefix sum value as a child of the chain
    forest_state.sibling_indices_at_depth[interval.depth].emplace_back(ZipCodeTree::CHAIN_START, 
        interval.is_reversed ? first_seed.zipcode.get_length(interval.depth, true) 
                             : 0);
    forest_state.sibling_indices_at_depth[interval.depth].back().chain_component = 
        (interval.is_reversed && !first_is_node) ? first_seed.zipcode.get_last_chain_component(interval.depth, true) 
                                                 : 0;

    // Remember the chain as a child of its parent snarl
    size_t chain_index = trees[forest_state.active_tree_index].zip_code_tree.size() - 1;
    vector<child_info_t>& snarl_children = forest_state.sibling_indices_at_depth[interval.depth-1];
    snarl_children.emplace_back(ZipCodeTree::CHAIN_START, chain_index);

    size_t chain_length = first_seed.zipcode.get_length(interval.depth, true);
    for (size_t i = interval.interval_start; i < interval.interval_end; i++) {
        size_t cur_seed_index = forest_state.seed_sort_order[i];
        const Seed& cur_seed = forest_state.seeds->at(cur_seed_index);
        bool cur_is_node = cur_seed.zipcode.max_depth() == interval.depth;
        // The distances in the snarl include the distances
        // from the first/last children in the chain to the ends of the chains
        //
        // Remember the distance to the start of each child in the chain
        if (interval.is_reversed) {
            // If the chain is reversed, find the distance to the end of the chain
            // from the seed's prefix sum and the length of the chain
            // If the length of the chain is infinite, then this is
            // not the last component of the chain and the distance is infinite
            // Otherwise, find the length of the chain/last component and
            // the length of the child, if it is a snarl
            
            if (chain_length == std::numeric_limits<size_t>::max()) {
                snarl_children.back().distances.first[cur_seed_index] = std::numeric_limits<size_t>::max();
            } else {
                size_t dist_val = forest_state.sort_values_by_seed[cur_seed_index].get_distance_value();
                bool child_is_node = cur_is_node || cur_seed.zipcode.get_code_type(interval.depth+1) == ZipCode::NODE;
                // Offset is distance plus possibly the child's length
                size_t offset = SnarlDistanceIndex::sum(dist_val, 
                    child_is_node ? 0
                                  : cur_seed.zipcode.get_length(interval.depth+1));
                // Distance to front of chain is chain's full length minus offset
                snarl_children.back().distances.first[cur_seed_index] = SnarlDistanceIndex::minus(chain_length, offset);
            }
        } else {
            // If the chain is forward, find the prefix sum of the first component
            if (!cur_is_node && cur_seed.zipcode.get_chain_component(interval.depth+1) != 0) {
                // If this isn't the first component, then it is infinite
                snarl_children.back().distances.first[cur_seed_index] = std::numeric_limits<size_t>::max();
            } else {
                // Otherwise, just the prefix sum
                snarl_children.back().distances.first[cur_seed_index] 
                    = forest_state.sort_values_by_seed[cur_seed_index].get_distance_value();
            }
        }
    }

    // Remember chain start and if first child was far enough to start a subtree
    forest_state.open_chains.emplace_back(chain_index, 
        snarl_children.back().distances.first[first_seed_index] > forest_state.distance_limit);
}

void ZipCodeTree::forget_past_index(size_t start_index,
                                    unordered_map<size_t, size_t>& removed_snarls,
                                    unordered_map<size_t, size_t>& removed_bounds) {
    // Associative-container erase idiom: https://stackoverflow.com/a/8234813
    for (auto it = snarl_start_indexes.cbegin(); it != snarl_start_indexes.cend();) {
        if (it->second >= start_index) {
            // Indexes must be shifted before being saved in removed_snarls
            removed_snarls[it->first] = it->second - start_index;
            it = snarl_start_indexes.erase(it);
        } else {
            ++it;
        }
    }

    for (auto it = bound_pair_indexes.cbegin(); it != bound_pair_indexes.cend();) {
        if (it->first >= start_index && it->second >= start_index) {
            // Indexes must be shifted before being saved in removed_bounds
            removed_bounds[it->first - start_index] = it->second - start_index;
            it = bound_pair_indexes.erase(it);
        } else if (it->first >= start_index != it->second >= start_index) {
            it = bound_pair_indexes.erase(it);
        } else {
            ++it;
        }
    }
}

bool ZipCodeForest::move_slice(forest_growing_state_t& forest_state, const size_t& depth) {
    // Chain IDs are reset since now there is no parent snarl
    tree_item_t new_start = tree_item_t(ZipCodeTree::CHAIN_START, std::numeric_limits<size_t>::max());
    tree_item_t new_end = tree_item_t(ZipCodeTree::CHAIN_END, std::numeric_limits<size_t>::max());
    ZipCodeTree& active_zip_tree = trees[forest_state.active_tree_index];

    // Might need to add an artificial end to the slice
    if (active_zip_tree.zip_code_tree.back().get_type() != ZipCodeTree::CHAIN_END) {
        active_zip_tree.zip_code_tree.emplace_back(new_end);
    }

    size_t slice_start_i = forest_state.open_chains.back().first;
    // Iterators to ends of slice, for insertion/erasure
    auto start_of_slice = active_zip_tree.zip_code_tree.begin() + slice_start_i;
    auto end_of_slice = active_zip_tree.zip_code_tree.end();
    // Are we moving a full chain or just a slice of one?
    bool move_full_chain = start_of_slice->get_type() == ZipCodeTree::CHAIN_START;

    // Pull out memory map of snarl start indexes to transfer to the new tree
    unordered_map<size_t, size_t> transferred_snarl_starts, transferred_pair_indexes;
    active_zip_tree.forget_past_index(slice_start_i, transferred_snarl_starts, transferred_pair_indexes);
    trees.emplace_back();
    ZipCodeTree& new_tree = trees.back();
    // Transfer index memory
    new_tree.snarl_start_indexes = transferred_snarl_starts;
    new_tree.bound_pair_indexes = transferred_pair_indexes;

    if (!move_full_chain) {
#ifdef DEBUG_ZIP_CODE_TREE
        cerr << "Copy a slice from the middle of the chain to the end" << endl;
        assert(start_of_slice->get_type() == ZipCodeTree::SEED || start_of_slice->is_snarl_start());
#endif
        // We're copying a slice of the chain from the middle to the end
        // Must add an artificial start to the slice
        new_tree.zip_code_tree.emplace_back(new_start);

        // Scoot all the indexes forward one
        new_tree.shift_past_index(0, 1);
    }

    // Copy everything in the slice into the new tree
    new_tree.zip_code_tree.insert(new_tree.zip_code_tree.end(),
            std::make_move_iterator(start_of_slice), std::make_move_iterator(end_of_slice));
    // Erase the slice
    trees[forest_state.active_tree_index].zip_code_tree.erase(start_of_slice, end_of_slice);

    // Force-reset chain ends to have no ID
    new_tree.zip_code_tree.front() = new_start;
    new_tree.zip_code_tree.pop_back(); // Remove chain end to reset its ID
    new_tree.add_close_bound(0);

#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "Validate slice" << endl;
    new_tree.print_self(forest_state.seeds);
    new_tree.validate_zip_tree(*forest_state.distance_index, forest_state.seeds, forest_state.distance_limit);
#endif
    return move_full_chain;
}

void ZipCodeForest::close_chain(forest_growing_state_t& forest_state, 
                                const size_t& depth, const Seed& last_seed, bool chain_is_reversed) {
    // Look up the ID of the parent snarl, unless this is a top-level chain
    auto snarl_id = depth > 0 ? forest_state.sibling_indices_at_depth[depth-1][0].snarl_id 
                              : std::numeric_limits<size_t>::max();
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "\t\tclose a chain in snarl " << snarl_id << " at depth " << depth << endl;
#endif
    if (trees[forest_state.active_tree_index].zip_code_tree.back().get_type() == ZipCodeTree::CHAIN_START) {
        // If the chain was empty (e.g. it only had a snarl which was removed)

        // Take out the CHAIN_START
        trees[forest_state.active_tree_index].zip_code_tree.pop_back();
        // Forget about the chain
        if (depth != 0) {
            forest_state.sibling_indices_at_depth[depth-1].pop_back();
            forest_state.open_chains.pop_back();
        }
    } else {
        // Otherwise, the chain wasn't empty so actually close it

        // Add the end of the chain to the zip code tree
        trees[forest_state.active_tree_index].add_close_bound(
            depth > 0 ? forest_state.sibling_indices_at_depth[depth-1].back().value
                      : 0);

        // If this was a top-level chain, we're done
        if (depth == 0) { return; }

        // For chains in snarls, we want to know the distance
        // from the last thing in the chain to the end of the chain
        // If further than the distance limit, we may slice off a subtree
        // If the chain remains in the snarl, remember the distance to the end
        // of the chain to add to the relevant distances in the parent snarl.
        // These distances will be stored in sibling_indices_at_depth

#ifdef DEBUG_ZIP_CODE_TREE
        assert(forest_state.sibling_indices_at_depth[depth-1].size() > 0);
        assert(forest_state.sibling_indices_at_depth[depth-1].back().type == ZipCodeTree::CHAIN_START);
#endif
        // Only add the distance for a non-root chain

        // If this is/isn't reversed, distance is to the start/end of the chain. 
        // The value in forest_state.sibling_indices_at_depth was the prefix sum
        // traversing the chain according to its tree orientation;
        // either way the distance is the length of the chain - the prefix sum
        size_t distance_to_chain_end = chain_is_reversed 
                                     ? forest_state.sibling_indices_at_depth[depth].back().value
                                     : SnarlDistanceIndex::minus(last_seed.zipcode.get_length(depth),
                                          forest_state.sibling_indices_at_depth[depth].back().value);
        
        // We might get rid of the chain
        bool chain_kept = true;
        if (distance_to_chain_end > forest_state.distance_limit && forest_state.open_chains.back().second) {
            // Distance to end & start of chain are both too far; splice out
            bool moved_full_chain = move_slice(forest_state, depth);

            if (moved_full_chain) {
                chain_kept = false;
                // The chain no longer exists in the snarl, so forget it
                forest_state.sibling_indices_at_depth[depth-1].pop_back();
            } else {
                // Take out the last edge, leading to the spliced-out part
                trees[forest_state.active_tree_index].zip_code_tree.pop_back();
                // Close the chain in the original active tree
                trees[forest_state.active_tree_index].add_close_bound(
                    forest_state.sibling_indices_at_depth[depth-1].back().value);
                // Treat distance as inf since it is over the distance limit
                distance_to_chain_end = std::numeric_limits<size_t>::max();
            }
        }

        // We didn't splice out the entire chain
        if (chain_kept) {
            // Remember chain's orientation and distance from last seed to edge
            forest_state.sibling_indices_at_depth[depth-1].back().is_reversed = chain_is_reversed;
            forest_state.sibling_indices_at_depth[depth-1].back().distances.second = distance_to_chain_end;
        }
        // We've closed a chain, so take out the latest open chain
        forest_state.open_chains.pop_back();
    }
}

void ZipCodeForest::add_child_to_chain(forest_growing_state_t& forest_state, const size_t& depth,
                                       const size_t& seed_index, bool child_is_reversed, bool chain_is_reversed) {
    const Seed& current_seed = forest_state.seeds->at(seed_index);
    bool current_is_reversed = (child_is_reversed != is_rev(current_seed.pos));

    // Would this seed be an exact duplicate of the previous one?
    tree_item_t& prev_item = trees[forest_state.active_tree_index].zip_code_tree.back();
    bool prev_is_duplicate = (prev_item.get_type() == ZipCodeTree::SEED 
                              && forest_state.seeds->at(prev_item.get_value()).pos == current_seed.pos
                              && prev_item.get_is_reversed() == current_is_reversed);
    if (prev_is_duplicate) {
#ifdef DEBUG_ZIP_CODE_TREE
        cerr << "Skipping duplicate seed at " << current_seed.pos << endl;
#endif
        prev_item.add_extra_value(seed_index);
        return;
    }

    // Possibly the last item is the same position but reversed,
    // and the seed before that is the duplicate
    if (trees[forest_state.active_tree_index].zip_code_tree.size() >= 3) {
        tree_item_t& three_back_item = trees[forest_state.active_tree_index].zip_code_tree[
            trees[forest_state.active_tree_index].zip_code_tree.size() - 3];
        bool three_back_is_duplicate = (three_back_item.get_type() == ZipCodeTree::SEED 
                                        && forest_state.seeds->at(three_back_item.get_value()).pos == current_seed.pos
                                        && three_back_item.get_is_reversed() == current_is_reversed);
        if (three_back_is_duplicate) {
#ifdef DEBUG_ZIP_CODE_TREE
            cerr << "Skipping duplicate seed at " << current_seed.pos << endl;
#endif
            three_back_item.add_extra_value(seed_index);
            return;
        }
    }

    // There will only be a relevant snarl ID if the chain is a child of a snarl
    auto snarl_id = depth >= 1 ? forest_state.sibling_indices_at_depth[depth-1][0].snarl_id 
                               : std::numeric_limits<size_t>::max();
    ZipCode::code_type_t current_type = current_seed.zipcode.get_code_type(depth);

    // Is this chain actually a node pretending to be a chain?
    bool is_trivial_chain = current_type == ZipCode::CHAIN && depth == current_seed.zipcode.max_depth();

    // For a root node or trivial chain, the "chain" is actually just the node,
    // so the depth of the chain is the same. Otherwise, depth is depth-1
    size_t chain_depth = is_trivial_chain || current_type == ZipCode::ROOT_NODE ? depth : depth-1;

    ///////////////// Get the offset in the parent chain (or node)
    size_t current_offset;

    // First, get the prefix sum in the chain + offset in the node
    if (current_type == ZipCode::ROOT_NODE || current_type == ZipCode::NODE || is_trivial_chain) {
        // For a node, this is still the distance used to sort on
        current_offset = forest_state.sort_values_by_seed[seed_index].get_distance_value();
    } else {
        // Otherwise, get the distance to the start or end of the chain
        current_offset = current_seed.zipcode.get_offset_in_chain(depth);
    }

    if (chain_is_reversed && !(current_type == ZipCode::NODE || current_type == ZipCode::ROOT_NODE 
                               || is_trivial_chain)) {
        // If we are adding a snarl and the chain is being traversed backwards,
        // then make sure the prefix sum is going to the right end of the snarl
        current_offset = SnarlDistanceIndex::sum(current_offset, current_seed.zipcode.get_length(depth));
    }

    /////////////////////// Get previous thing's offset in parent chain/node
    size_t previous_offset = forest_state.sibling_indices_at_depth[chain_depth][0].value;

#ifdef DEBUG_ZIP_CODE_TREE
    assert(forest_state.sibling_indices_at_depth[chain_depth].size() == 1);
#endif

    ///////////////////// Record distance from the previous thing in chain/node
    //       Or add a new tree if the distance is too far
    auto first_type = forest_state.sibling_indices_at_depth[chain_depth][0].type;
    if (chain_depth > 0 && first_type == ZipCodeTree::CHAIN_START) {
        // If this is the first thing in a non-root chain or node,
        // check whether it is connected to the front of the snarl
        const auto& flank_distances = forest_state.sibling_indices_at_depth[chain_depth-1].back().distances.first;
        forest_state.open_chains.back().second = flank_distances.at(seed_index) > forest_state.distance_limit;
    } else if (first_type != ZipCodeTree::CHAIN_START) {
        // For all except the first thing in a node/chain, we need to add an edge
        size_t distance_between = std::numeric_limits<size_t>::max();
        if (is_trivial_chain || current_type == ZipCode::ROOT_NODE || 
            forest_state.sibling_indices_at_depth[chain_depth][0].chain_component 
              == current_seed.zipcode.get_chain_component(depth)) {
            if (previous_offset != std::numeric_limits<size_t>::max()) {
                // In same & reachable chain component, so can find distance from last thing
                distance_between = (std::max(current_offset, previous_offset) 
                                    - std::min(current_offset, previous_offset));
            }
        }

        if (chain_depth == 0 && distance_between > forest_state.distance_limit) {
            // The next thing in the zip tree will be the first seed (or snarl)
            // in a top-level chain, so start a new tree
#ifdef DEBUG_ZIP_CODE_TREE
            cerr << "Start a new tree in the forest" << endl;
#endif
            // Close the previous chain
            trees[forest_state.active_tree_index].add_close_bound(0);

            if (forest_state.active_tree_index == std::numeric_limits<size_t>::max() 
                || trees[forest_state.active_tree_index].zip_code_tree.size() != 0) {
                // Add a new tree and make sure it is the new active tree
#ifdef DEBUG_ZIP_CODE_TREE
                // If we're starting a new tree then the last one must be valid
                if (forest_state.active_tree_index != std::numeric_limits<size_t>::max()) {
                    cerr << "Last tree: " << endl;
                    trees[forest_state.active_tree_index].print_self(forest_state.seeds);
                    trees[forest_state.active_tree_index].validate_zip_tree(
                        *forest_state.distance_index, forest_state.seeds, forest_state.distance_limit);
                }
#endif
                trees.emplace_back();
                forest_state.active_tree_index = trees.size()-1;
            }

            // Add the start of the new chain (top-level chains have ID=inf)
            trees[forest_state.active_tree_index].zip_code_tree.emplace_back(
                ZipCodeTree::CHAIN_START, std::numeric_limits<size_t>::max());

            // The first sibling is now the chain start, not the previous seed
            forest_state.sibling_indices_at_depth[chain_depth].pop_back();
            forest_state.sibling_indices_at_depth[chain_depth].emplace_back(ZipCodeTree::CHAIN_START, 
                chain_is_reversed ? current_seed.zipcode.get_length(chain_depth, true) : 0); 
            forest_state.sibling_indices_at_depth[chain_depth].back().chain_component = 
                !is_trivial_chain ? current_seed.zipcode.get_last_chain_component(chain_depth, true) : 0;
        } else if (distance_between > forest_state.distance_limit) { 
            // Too far from the previous thing, but inside a snarl
            if (forest_state.open_chains.back().second) {
                bool moved_full_chain = move_slice(forest_state, chain_depth);
                // Current chain slice was also too far away from prior thing
                if (moved_full_chain) {
                    // Since the full chain was moved, the back is a bound
                    // We can safely copy its ID
                    snarl_id = trees[forest_state.active_tree_index].zip_code_tree.back().get_value();
                    // Add back the start of the chain
                    trees[forest_state.active_tree_index].zip_code_tree.emplace_back(
                        ZipCodeTree::CHAIN_START, snarl_id);

                    // Update the chain as a child of the snarl
#ifdef DEBUG_ZIP_CODE_TREE
                    assert(forest_state.sibling_indices_at_depth[chain_depth-1].back().type \
                           == ZipCodeTree::CHAIN_START);   
                    // The value should be the index of the last seed,
                    // which is the first seed in the new tree
                    assert(forest_state.sibling_indices_at_depth[chain_depth-1].back().value 
                                == trees[forest_state.active_tree_index].zip_code_tree.size()-1);
                    assert(forest_state.open_chains.back().second);
#endif
                    // Don't need to update open_chains, since the next slice 
                    // will also start at the chain start
                } else {
                    // If the slice starts and ends in the middle of the chain
                    // replace with edge with infinite length, since
                    // it will be bigger than the distance limit anyway
#ifdef DEBUG_ZIP_CODE_TREE
                    assert(trees[forest_state.active_tree_index].zip_code_tree.back().get_type() == ZipCodeTree::EDGE);
#endif
                    trees[forest_state.active_tree_index].zip_code_tree.pop_back();
                    trees[forest_state.active_tree_index].zip_code_tree.emplace_back(
                        ZipCodeTree::EDGE, std::numeric_limits<size_t>::max());

                    // Remember the next seed or snarl that gets added
                    // as the start of a new chain slice
                    forest_state.open_chains.pop_back();
                    forest_state.open_chains.emplace_back(
                        trees[forest_state.active_tree_index].zip_code_tree.size(), true);
                }
            } else {
#ifdef DEBUG_ZIP_CODE_TREE
                 cerr << "The slice didn't get copied but maybe start a new slice here" << endl;
#endif
                // If the slice is still connected at the front,
                // add an edge and remember that it could start a new slice

                trees[forest_state.active_tree_index].zip_code_tree.emplace_back(ZipCodeTree::EDGE, distance_between);

                // Remember the next seed or snarl that gets added 
                // as the start of a new chain slice
                forest_state.open_chains.pop_back();
                forest_state.open_chains.emplace_back(
                    trees[forest_state.active_tree_index].zip_code_tree.size(), true);
            }
        } else {
            // If we didn't start a new tree, then add the edge
            trees[forest_state.active_tree_index].zip_code_tree.emplace_back(ZipCodeTree::EDGE, distance_between);
        }
    }

    /////////////////////////////Record this thing in the chain
    if (current_type == ZipCode::NODE || current_type == ZipCode::ROOT_NODE || is_trivial_chain) {
#ifdef DEBUG_ZIP_CODE_TREE
        cerr << "\t\tContinue node/chain with seed " << current_seed.pos << " at depth " << depth << endl;
#endif
        // If this was a node, just remember the seed
        trees[forest_state.active_tree_index].zip_code_tree.emplace_back(
            ZipCodeTree::SEED, seed_index, current_is_reversed);
    } else {
        open_snarl(forest_state, depth, current_seed.zipcode.get_code_type(depth) == ZipCode::CYCLIC_SNARL); 

        // For finding the distance to the next thing in the chain, the offset
        // is the offset of the end bound of the snarl, so add the snarl length
        current_offset = chain_is_reversed 
                       ? SnarlDistanceIndex::minus(current_offset, current_seed.zipcode.get_length(depth))
                       : SnarlDistanceIndex::sum(current_offset, current_seed.zipcode.get_length(depth));
    }

    // Remember this thing for the next sibling in the chain
    forest_state.sibling_indices_at_depth[chain_depth].pop_back();
    tree_item_t just_added = trees[forest_state.active_tree_index].zip_code_tree.back();
    if (just_added.get_type() == ZipCodeTree::SEED) {
        forest_state.sibling_indices_at_depth[chain_depth].emplace_back(ZipCodeTree::SEED, current_offset);
    } else {
#ifdef DEBUG_ZIP_CODE_TREE
        assert(just_added.is_snarl_start());
#endif
        forest_state.sibling_indices_at_depth[chain_depth].emplace_back(just_added.get_type(), current_offset); 
        forest_state.sibling_indices_at_depth[chain_depth].back().snarl_id = just_added.get_value();
    }

    if (!is_trivial_chain && current_type != ZipCode::ROOT_NODE) {
        // If needed, remember the chain component
        forest_state.sibling_indices_at_depth[chain_depth].back().chain_component 
            = current_seed.zipcode.get_chain_component(depth);
    }
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "Added sibling with type " << current_type << endl;
#endif
}

void ZipCodeForest::open_snarl(forest_growing_state_t& forest_state, const size_t& depth, bool is_cyclic_snarl) {
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "\t\tOpen new " << (is_cyclic_snarl ? "cyclic" : "DAG") << " snarl at depth " << depth << endl;
#endif
    trees[forest_state.active_tree_index].open_snarl(is_cyclic_snarl);
    tree_item_t snarl_start = trees[forest_state.active_tree_index].zip_code_tree.back();

    // Remember the start of the snarl for distances & ID
    forest_state.sibling_indices_at_depth[depth].emplace_back(
        snarl_start.get_type(), std::numeric_limits<size_t>::max());
    forest_state.sibling_indices_at_depth[depth].back().snarl_id = snarl_start.get_value();
}

void ZipCodeTree::open_snarl(bool is_cyclic_snarl) {
    size_t cur_id = next_snarl_id++;
    zip_code_tree.emplace_back(is_cyclic_snarl ? CYCLIC_SNARL_START : DAG_SNARL_START, cur_id);
    // Remember where the start was put
    snarl_start_indexes[cur_id] = zip_code_tree.size() - 1;
}

void ZipCodeForest::close_snarl(forest_growing_state_t& forest_state,
                                const size_t& depth, const Seed& last_seed, bool last_is_reversed) {
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "\t\tclose a snarl at depth " << depth << endl;
#endif
    child_info_t snarl_start = forest_state.sibling_indices_at_depth[depth][0];

    if (trees[forest_state.active_tree_index].zip_code_tree.size() == 1) {
        // If this would be an empty snarl, then just remove it
        trees.erase(trees.begin() + forest_state.active_tree_index);
        trees[forest_state.active_tree_index].snarl_start_indexes.erase(snarl_start.snarl_id);
    } else if (depth == 0) {
        // If this is a root snarl, we don't need distances so just close it
        // Root snarls are treated as DAG snarls, but they don't store distances
        trees[forest_state.active_tree_index].add_close_bound(0);
    } else if (forest_state.sibling_indices_at_depth[depth].size() == 1) {
        // Only "child" is snarl start since all chains were spliced out
        // Thus, this snarl is empty and we should remove it
#ifdef DEBUG_ZIP_CODE_TREE
        cerr << "\t\t\tThe snarl is actually empty so remove it" << endl;
        assert(trees[forest_state.active_tree_index].zip_code_tree.back().is_snarl_start());
#endif        
        // Pop the snarl start out
        trees[forest_state.active_tree_index].zip_code_tree.pop_back();
        trees[forest_state.active_tree_index].snarl_start_indexes.erase(snarl_start.snarl_id);

        // If this was the first thing in the chain, then we're done.
        // Otherwise, there was an edge to remove
        if (trees[forest_state.active_tree_index].zip_code_tree.back().get_type() == ZipCodeTree::EDGE) {
            // Take out chain edge and update previous thing in the chain

            // Distance from the last thing to the snarl start
            size_t previous_edge = trees[forest_state.active_tree_index].zip_code_tree.back().get_value();
            trees[forest_state.active_tree_index].zip_code_tree.pop_back();
            // Distance from the start of the chain to the snarl end
            size_t snarl_prefix_sum = forest_state.sibling_indices_at_depth[depth-1].back().value;
            forest_state.sibling_indices_at_depth[depth-1].pop_back();

            // Shift prefix sum to the other side of the snarl
            size_t snarl_length = last_seed.zipcode.get_length(depth);
            snarl_prefix_sum = last_is_reversed ? SnarlDistanceIndex::sum(snarl_prefix_sum, snarl_length)
                                                : SnarlDistanceIndex::minus(snarl_prefix_sum, snarl_length);
            tree_item_t last_item = trees[forest_state.active_tree_index].zip_code_tree.back();
            tree_item_type_t last_type = last_item.get_type();
            if (last_item.is_snarl_end()) {
                // Fix snarl ends to actually be snarl starts
                last_type = last_item.get_type() == ZipCodeTree::CYCLIC_SNARL_END 
                    ? ZipCodeTree::CYCLIC_SNARL_START : ZipCodeTree::DAG_SNARL_START;
            }
            size_t distance = last_is_reversed ? SnarlDistanceIndex::sum(snarl_prefix_sum, previous_edge)
                                               : SnarlDistanceIndex::minus(snarl_prefix_sum, previous_edge);
            // Now update sibling_indices_at_depth to be previous thing in chain
            forest_state.sibling_indices_at_depth[depth-1].emplace_back(last_type, distance);

            // If it was in the first component, then this is correct.
            // If a later component, then it's too far away anyway
            forest_state.sibling_indices_at_depth[depth-1].back().chain_component = 0;

            // At this point, the open_chain for the parent chain is either
            // before the removed snarl, the snarl itself, or after the snarl.
            // If the open_chain was before or at the snarl, nothing has changed.
            // If it is after the snarl, then the snarl
            // wasn't the start of a new slice so we back it up to the previous
            // child and say that it was not the start of a new slice.
            // TODO
            // If it was the snarl itself, then the next child added to the chain
            // will be the next open_chain, but I
            // haven't implemented this yet- it won't change the correctness
            bool chain_after_snarl = forest_state.open_chains.back().first
                   >= trees[forest_state.active_tree_index].zip_code_tree.size();
            if (depth > 0 && forest_state.open_chains.size() > 0 && chain_after_snarl) {
                // If a chain slice could have started at or after this snarl
#ifdef DEBUG_ZIP_CODE_TREE
                assert(forest_state.open_chains.back().second);
#endif
                // Find the start of the previous child
                size_t previous_index = trees[forest_state.active_tree_index].zip_code_tree.size() - 1;
                // Snarl nestedness
                size_t opened_snarls = 0;
                tree_item_t previous_item;

                // Hunt backwards in ziptree for previous child
                while (true) {
                    previous_item = trees[forest_state.active_tree_index].zip_code_tree.at(previous_index);
                    
                    if (opened_snarls == 0 && (previous_item.get_type() == ZipCodeTree::SEED
                                               || previous_item.is_snarl_start())) {
                        // Found a non-nested previous item (snarl or seed)
                        break;
                    } else if (previous_item.is_snarl_end()) {
                        // Entering nested snarl
                        opened_snarls++;
                    } else if (previous_item.is_snarl_start()) {
                        // Leaving nested snarl
                        opened_snarls--;
                    }
                    // Try again, one back
                    previous_index--;
                }

                previous_item = trees[forest_state.active_tree_index].zip_code_tree.at(previous_index);
#ifdef DEBUG_ZIP_CODE_TREE
                assert(previous_item.get_type() == ZipCodeTree::SEED || previous_item.is_snarl_start());
                cerr << "New start of previous open chain: " << previous_index << endl;;
#endif
                forest_state.open_chains.back().first = previous_index;
                forest_state.open_chains.back().second = false;      
            }
#ifdef DEBUG_ZIP_CODE_TREE
            assert(forest_state.sibling_indices_at_depth[depth-1].back().value >= 0);
#endif
        } else {
            // If this was the first thing in the chain,
            // update the previous sibling in the chain to be the chain start
#ifdef DEBUG_ZIP_CODE_TREE
            assert(trees[forest_state.active_tree_index].zip_code_tree.back().get_type() == ZipCodeTree::CHAIN_START);
#endif
            bool prev_reversed = forest_state.open_intervals[forest_state.open_intervals.size()-2].is_reversed;
            size_t offset = prev_reversed ? last_seed.zipcode.get_length(depth-1, true) : 0;
            forest_state.sibling_indices_at_depth[depth-1].pop_back();
            forest_state.sibling_indices_at_depth[depth-1].emplace_back(ZipCodeTree::CHAIN_START, offset);
            forest_state.sibling_indices_at_depth[depth-1].back().chain_component
                = last_seed.zipcode.get_last_chain_component(depth-1, true);
        }
    } else {
        // Snarl was not removed, so finish off with distance matrix & end bound
        add_distance_matrix(forest_state, depth, last_is_reversed);
        auto closing_type = snarl_start.type == ZipCodeTree::CYCLIC_SNARL_START ? ZipCodeTree::CYCLIC_SNARL_END
                                                                                : ZipCodeTree::DAG_SNARL_END;
        size_t snarl_start_i = trees[forest_state.active_tree_index].snarl_start_indexes.at(snarl_start.snarl_id);
        trees[forest_state.active_tree_index].add_close_bound(snarl_start_i);
     }
}

void ZipCodeTree::shift_past_index(size_t start_index, size_t shift_amount) {
    for (auto& snarl_start : snarl_start_indexes) {
        if (snarl_start.second >= start_index) {
            // Shift the snarl start index
            snarl_start.second += shift_amount;
        }
    }

    unordered_map<size_t, size_t> changed_pair_indexes;

    // Associative-container erase idiom: https://stackoverflow.com/a/8234813
    for (auto it = bound_pair_indexes.cbegin(); it != bound_pair_indexes.cend();) {
        if (it->second >= start_index) {
            // Store new pair and then delete the old
            changed_pair_indexes[it->first + shift_amount] = it->second + shift_amount;
            it = bound_pair_indexes.erase(it);
        } else {
            ++it;
        }
    }

    bound_pair_indexes.insert(changed_pair_indexes.begin(), changed_pair_indexes.end());
}

ZipCodeForest::seed_info_t::seed_info_t(size_t index, bool is_rev, size_t flank, bool right_side,
                                        size_t nested_snarl_offset,
                                        const size_t& depth, const forest_growing_state_t& forest_state) : 
                                        zipcode(forest_state.seeds->at(index).zipcode), flank_offset(flank),
                                        right_side(right_side), nested_snarl_offset(nested_snarl_offset) {
    const Seed& seed = forest_state.seeds->at(index);
    // Pull seed metadata
    rank = zipcode.get_rank_in_snarl(depth+1);
    node_length = zipcode.get_length(zipcode.max_depth());
    // Find position pointing out of chain
    seed_pos = seed.pos;
    if (is_rev) {
        seed_pos = reverse_seed();
    }
}

vector<ZipCodeForest::seed_info_t> ZipCodeForest::get_edge_seeds(const forest_growing_state_t& forest_state, 
                                                                 const size_t& depth) const {
    vector<seed_info_t> edge_seeds;
    // Convenience reference
    const auto& active_zip_tree = trees[forest_state.active_tree_index].zip_code_tree;
    const auto& chain_indices = forest_state.sibling_indices_at_depth[depth];
    // Get edge seeds for each chain
    for (size_t i = 1; i < chain_indices.size(); ++i) {
        child_info_t cur_chain = chain_indices[i];
        
        // Search for first seed, starting from right after chain start
        size_t seed_i = cur_chain.value + 1;
        size_t nested_snarl_offset = trees[forest_state.active_tree_index].get_offset_to_seed(seed_i, false);
        size_t seed_val = active_zip_tree[seed_i].get_value();
        edge_seeds.emplace_back(seed_val, !active_zip_tree[seed_i].get_is_reversed(),
            cur_chain.distances.first[seed_val], cur_chain.is_reversed, nested_snarl_offset, depth, forest_state);

        // Search for last seed, starting from right before chain end
        seed_i = i == chain_indices.size() - 1 ? active_zip_tree.size() - 2
                                               : chain_indices[i+1].value - 2;
        nested_snarl_offset = trees[forest_state.active_tree_index].get_offset_to_seed(seed_i, true);
        edge_seeds.emplace_back(active_zip_tree[seed_i].get_value(), active_zip_tree[seed_i].get_is_reversed(),
            cur_chain.distances.second, !cur_chain.is_reversed, nested_snarl_offset, depth, forest_state);
    }

    return edge_seeds;
}

size_t ZipCodeTree::get_offset_to_seed(size_t& i, bool right_to_left) const {
    if (zip_code_tree[i].get_type() == ZipCodeTree::SEED) {
        // If we're already a seed, then no need to find an offset
        return 0;
    } else {
#ifdef DEBUG_ZIP_CODE_TREE
        // Edge seed is in a nested snarl
        assert((right_to_left && zip_code_tree[i].is_snarl_end())
            || (!right_to_left && zip_code_tree[i].is_snarl_start()));
#endif
        size_t snarl_start = snarl_start_indexes.at(zip_code_tree[i].get_value());
        bool is_cyclic_snarl = zip_code_tree[snarl_start].get_type() == ZipCodeTree::CYCLIC_SNARL_START;
        // Chain count stored one after snarl start
        size_t chain_count = zip_code_tree[snarl_start+1].get_value();
        size_t offset = 0;

        if (right_to_left) {
            // How many rows in distance matrix?
            size_t row_count = is_cyclic_snarl ? (chain_count + 1) * 2
                                               : chain_count + 1;
            // nth triangular number
            size_t dist_matrix_size = row_count * (row_count + 1) / 2;
            // Last chain's right to snarl end is last in DAG dist matrix
            // or 2nd-to-last in cyclic dist matrix
            // +1for CHAIN_COUNT is cancelled out by -1 for 2nd-to-last
            size_t i_offset = is_cyclic_snarl ? 0 : 1;
            offset += zip_code_tree[snarl_start + dist_matrix_size + i_offset].get_value();

            // Skip snarl end & chain end
            i -= 2;
            // Recurse inside nested chain
            return offset + get_offset_to_seed(i, right_to_left);
        } else {
            // First's chain left to snarl start is 1st in DAG dist matrix
            // or 2nd in cyclic dist matrix
            size_t dist_matrix_i = is_cyclic_snarl ? 2 : 1;
            // Also a +1 for CHAIN_COUNT
            offset += zip_code_tree[snarl_start + 1 + dist_matrix_i].get_value();

            // Skip to first chain start
            while (zip_code_tree[i].get_type() != ZipCodeTree::CHAIN_START) {
                i++;
            }
            // Jump over chain start
            i++;
            // Recurse inside nested chain
            return offset + get_offset_to_seed(i, right_to_left);
        }
    }
}

void ZipCodeForest::add_edges_to_end(vector<tree_item_t>& dist_matrix,
                                     forest_growing_state_t& forest_state, 
                                     const size_t& depth, const vector<seed_info_t>& edge_seeds,
                                     bool snarl_is_reversed, bool is_cyclic_snarl) const {
    bool is_regular_snarl = edge_seeds[0].zipcode.get_code_type(depth) == ZipCode::REGULAR_SNARL;
    // start -> end is simply length of snarl
    dist_matrix.emplace_back(ZipCodeTree::EDGE, edge_seeds[0].zipcode.get_length(depth));

    // DAG snarls only have distances from chain ends to snarl end
    // seed 1 -> end, 3 -> end, etc.
    // Cyclic snarls have distances from all sides to the end
    size_t start_i = is_cyclic_snarl ? 0 : 1;
    size_t increment = is_cyclic_snarl ? 1 : 2;
    size_t edge_dist;
    for (size_t i = start_i; i < edge_seeds.size(); i += increment) {
        if (is_regular_snarl) {
            // Regular snarls just use the flank distances
            edge_dist = edge_seeds[i].flank_offset;
        } else {
            // Distance from the start of the snarl to the start of the chain
            size_t between_bounds_dist = edge_seeds[i].zipcode.get_distance_to_snarl_bound(
                depth+1, snarl_is_reversed, !edge_seeds[i].right_side);
            
            // Overall edge distance
            edge_dist = SnarlDistanceIndex::sum(between_bounds_dist, edge_seeds[i].flank_offset);
        }
        dist_matrix.emplace_back(ZipCodeTree::EDGE, edge_dist);
    }

    if (is_cyclic_snarl) {
        // end -> end may be possible, but is not optimal so we ignore it
        dist_matrix.emplace_back(ZipCodeTree::EDGE, std::numeric_limits<size_t>::max());
    }
}

void ZipCodeForest::add_edges_for_chains(vector<tree_item_t>& dist_matrix,
                                         forest_growing_state_t& forest_state, 
                                         const size_t& depth, const vector<seed_info_t>& edge_seeds,
                                         bool snarl_is_reversed, bool is_cyclic_snarl) const {
    // Grab snarl handle if necessary
    bool is_regular_snarl = edge_seeds[0].zipcode.get_code_type(depth) == ZipCode::REGULAR_SNARL;
    net_handle_t snarl_handle;
    if (!is_regular_snarl) {
        // Snarl handles are only used for irregular and cyclic snarls
        snarl_handle = edge_seeds[0].zipcode.get_net_handle(depth, forest_state.distance_index);
    }

    // Distance in between chain bounds
    size_t between_chain_dist;
    // Distance from chain bounds to the edge of the seeds
    size_t chain_flank_dist;
    // Final edge distance
    size_t edge_dist;

    // DAG snarls find distances FROM ends TO starts only
    // seed 0 ->, seed 2 ->, etc.
    size_t increment = is_cyclic_snarl ? 1 : 2;
    for (size_t to_i = 0; to_i < edge_seeds.size(); to_i += increment) {
        // Edge from start
        between_chain_dist = edge_seeds[to_i].zipcode.get_distance_to_snarl_bound(
            depth+1, !snarl_is_reversed, !edge_seeds[to_i].right_side);
        edge_dist = SnarlDistanceIndex::sum(between_chain_dist, edge_seeds[to_i].flank_offset);
        dist_matrix.emplace_back(ZipCodeTree::EDGE, edge_dist);

        // TO chain starts means -> seed 1, -> seed 3, etc.
        size_t start_from = is_cyclic_snarl ? 0 : 1;
        for (int64_t from_i = start_from; from_i <= to_i; from_i += increment) {

            if (to_i / 2 == from_i / 2) {
                // Intra-chain distance; flank distances are untrustworthy
                // if the true self-loop path stays within the chain

                // to_pos points "out", but we need to come "in" to it
                pos_t to_pos = edge_seeds[to_i].reverse_seed();

                edge_dist = minimum_nontrivial_distance(*forest_state.distance_index,
                    edge_seeds[from_i].seed_pos, to_pos, edge_seeds[to_i].node_length);
                // Subtract distance from inner seed to edge of inner snarl
                edge_dist = SnarlDistanceIndex::minus(edge_dist, edge_seeds[from_i].nested_snarl_offset);
                edge_dist = SnarlDistanceIndex::minus(edge_dist, edge_seeds[to_i].nested_snarl_offset);
            } else {
                // Inter-chain distance; assume we pass through chain bounds
                chain_flank_dist = SnarlDistanceIndex::sum(edge_seeds[from_i].flank_offset, 
                                                           edge_seeds[to_i].flank_offset);
                if (is_regular_snarl) {
                    // Distance between chains in a regular snarl is always inf
                    between_chain_dist = std::numeric_limits<size_t>::max();
                } else {
                    // Irregular/cyclic snarls have a snarl_handle we can use
                    between_chain_dist = forest_state.distance_index->distance_in_snarl(snarl_handle, 
                            edge_seeds[from_i].rank, edge_seeds[from_i].right_side, 
                            edge_seeds[to_i].rank, edge_seeds[to_i].right_side);
                }
                edge_dist = SnarlDistanceIndex::sum(between_chain_dist, chain_flank_dist);
            }

            dist_matrix.emplace_back(ZipCodeTree::EDGE, edge_dist);
        }
    }
}

void ZipCodeForest::add_distance_matrix(forest_growing_state_t& forest_state, 
                                        const size_t& depth, bool snarl_is_reversed) {
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "\t\tadd distances for snarl at depth " << depth << endl;
    print_self(forest_state.seeds);
#endif
    vector<seed_info_t> edge_seeds = get_edge_seeds(forest_state, depth);

    // Metadata about the snarl
    bool is_cyclic_snarl = forest_state.sibling_indices_at_depth[depth][0].type == ZipCodeTree::CYCLIC_SNARL_START;
    size_t sibling_count = forest_state.sibling_indices_at_depth[depth].size();

    // Set up distance matrix
    size_t num_edges = is_cyclic_snarl ? (sibling_count*2) * (sibling_count*2 + 1) / 2
                                       : sibling_count * (sibling_count + 1) / 2;
    vector<tree_item_t> dist_matrix;
    dist_matrix.reserve(num_edges + 1);

    // Child count (all siblings minus the snarl start) preceeds the matrix
    // This provides context to make reading the matrix easier
    dist_matrix.emplace_back(ZipCodeTree::CHAIN_COUNT, sibling_count - 1);

    if (is_cyclic_snarl) {
        // start -> start may be possible, but is not optimal so we ignore it
        dist_matrix.emplace_back(ZipCodeTree::EDGE, std::numeric_limits<size_t>::max());
    }

    add_edges_for_chains(dist_matrix, forest_state, depth, edge_seeds, snarl_is_reversed, is_cyclic_snarl);
    add_edges_to_end(dist_matrix, forest_state, depth, edge_seeds, snarl_is_reversed, is_cyclic_snarl);

#ifdef DEBUG_ZIP_CODE_TREE
    size_t num_chains = dist_matrix.front().get_value();
    cerr << "Distance matrix: (" << num_chains << " chain" << (num_chains == 1 ? "" : "s") << ")" << endl;
    size_t matrix_i = 1;
    size_t rows = is_cyclic_snarl ? sibling_count * 2 : sibling_count;
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j <= i; j++) {
            size_t dist = dist_matrix[matrix_i++].get_value();
            if (dist == std::numeric_limits<size_t>::max()) {
                cerr << "inf\t";
            } else {
                cerr << dist << "\t";
            }
        }
        cerr << endl;
    }
    assert(dist_matrix.size() == num_edges + 1);
#endif

    // Put matrix at the start of the snarl
    size_t first_child_index = forest_state.sibling_indices_at_depth[depth][1].value;
    auto& active_zip_tree = trees[forest_state.active_tree_index];
    active_zip_tree.zip_code_tree.insert(active_zip_tree.zip_code_tree.begin() + first_child_index, 
                                         dist_matrix.begin(), dist_matrix.end());
    // Shift memorized snarl start indexes
    size_t cur_snarl_index = active_zip_tree.snarl_start_indexes.at(
        forest_state.sibling_indices_at_depth[depth][0].snarl_id);
    active_zip_tree.shift_past_index(cur_snarl_index + 1, num_edges + 1);
}

std::pair<size_t, size_t> ZipCodeTree::dag_and_cyclic_snarl_count() const {
    size_t dag_count = 0;
    size_t cyclic_count = 0;

    for (const auto& item : zip_code_tree) {
        if (item.get_type() == ZipCodeTree::CYCLIC_SNARL_START) {
            cyclic_count++;
        } else if (item.get_type() == ZipCodeTree::DAG_SNARL_START) {
            dag_count++;
        }
    }

    return std::make_pair(dag_count, cyclic_count);
}

bool ZipCodeTree::seed_is_reversed_at_depth(const Seed& seed, size_t depth, const SnarlDistanceIndex& distance_index) {
    if (seed.zipcode.get_is_reversed_in_parent(depth)) {
        return true;
    } else if (depth > 0 && (seed.zipcode.get_code_type(depth-1) == ZipCode::IRREGULAR_SNARL
                             || seed.zipcode.get_code_type(depth-1) == ZipCode::CYCLIC_SNARL)) {
        // If parent is an irregular snarl, check orientation of child in snarl
        net_handle_t snarl_handle = seed.zipcode.get_net_handle(depth-1, &distance_index);
        size_t rank = seed.zipcode.get_rank_in_snarl(depth);
        if (distance_index.distance_in_snarl(snarl_handle, 0, false, rank, false)
                    == std::numeric_limits<size_t>::max()
            &&
            distance_index.distance_in_snarl(snarl_handle, 1, false, rank, true)
                    == std::numeric_limits<size_t>::max()) {
            // If the distance from the snarl start to the child start is inf
            // and the distance from the snarl end to the child end is inf
            // then we assume that this child is "reversed" in the parent snarl
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
}

bool ZipCodeTree::node_is_invalid(nid_t id, const SnarlDistanceIndex& distance_index, size_t distance_limit) const {
    bool is_invalid = false;
    net_handle_t net = distance_index.get_node_net_handle(id);
    // Go upwards in snarl tree to see if any level is invalid
    while (!distance_index.is_root(net)) {
        if (distance_index.is_looping_chain(net)) {
            is_invalid = true;
            break;
        } else if (distance_index.is_chain(distance_index.get_parent(net)) && 
                    !distance_index.is_trivial_chain(distance_index.get_parent(net))) {
            // Check if this net_handle_t could be involved in
            // a chain loop that is smaller than the distance limit
            size_t forward_loop = distance_index.is_node(net) 
                                ? distance_index.get_forward_loop_value(net)
                                : distance_index.get_forward_loop_value(
                                        distance_index.get_node_from_sentinel(
                                            distance_index.get_bound(net, true, false)));
            size_t reverse_loop = distance_index.is_node(net) 
                                ? distance_index.get_reverse_loop_value(net)
                                : distance_index.get_reverse_loop_value(    
                                     distance_index.get_node_from_sentinel(
                                        distance_index.get_bound(net, false, false)));
            if (forward_loop < distance_limit ||
                reverse_loop < distance_limit) {
                is_invalid = true;
                break;
            }
        }
        net = distance_index.get_parent(net);
    }
    if (distance_index.is_root_snarl(net)) {
        is_invalid = true;
    }
                
    return is_invalid;
}

void ZipCodeTree::validate_boundaries(const SnarlDistanceIndex& distance_index, 
                                      const vector<Seed>* seeds,
                                      size_t distance_limit) const {
#ifdef DEBUG_ZIP_CODE_TREE
    std::cerr << "Validating that zip code tree's boundaries match up" << std::endl;
#endif
    // We want to have at least one seed
    bool has_seed = false;
    // Open bounds in snarl tree (e.g. chain, below that a parent snarl, etc.)
    std::stack<tree_item_t> tree_stack;
    for (size_t i = 0 ; i < zip_code_tree.size() ; i++) {
        const tree_item_t& item = zip_code_tree[i];
        if (item.is_snarl_start()) {
            // If there is a top-level chain, top-level snarls are depth 1
            // If there's a root snarl, top-level snarls are depth 2
            if (tree_stack.size() == 1 || tree_stack.size() == 2) {
                // Also check snarl distances and child count 
                // for non-root top-level snarls (is recursive)
                vector<tree_item_t>::const_iterator cur_snarl_start = zip_code_tree.begin() + i;
                validate_snarl(cur_snarl_start, distance_index, seeds, distance_limit);
            }

            // Snarls are not allowed to have inf IDs
            assert(item.get_value() != std::numeric_limits<size_t>::max());
            tree_stack.push(item);
        } else if (item.get_type() == CHAIN_START) {
            if (tree_stack.empty()) {
                // Top-level chain should have a snarl ID of inf
                assert(item.get_value() == std::numeric_limits<size_t>::max());
            } else {
                // Child chains should have the same snarl ID as parent snarl
                assert(tree_stack.top().get_value() == item.get_value());
            }
            tree_stack.push(item);
        } else if (item.is_snarl_end()) {
            // Should have opened with the correct snarl type
            assert(tree_stack.top().get_type() == (item.get_type() == DAG_SNARL_END ? DAG_SNARL_START 
                                                                                    : CYCLIC_SNARL_START));
            // IDs should match
            assert(tree_stack.top().get_value() == item.get_value());
            
            tree_stack.pop();
            // Either this was a top-level snarl, or there's a parent chain
            assert(tree_stack.empty() || tree_stack.top().get_type() == CHAIN_START);
        } else if (item.get_type() == CHAIN_END) {
            // Snarl ID should match, except for top-level chain
            assert(tree_stack.size() == 1 || tree_stack.top().get_value() == item.get_value());

            assert(tree_stack.top().get_type() == CHAIN_START);
            tree_stack.pop();
            // Either this was a top-level chain, or there's a parent snarl
            assert(tree_stack.empty() || tree_stack.top().is_snarl_start());
        } else if (item.get_type() == SEED) {
            has_seed = true;
        }
    }
    assert(has_seed);

    for (const auto& index_pair : bound_pair_indexes) {
        // Check that all bound pairs point to valid snarl starts/ends
        assert(index_pair.first < zip_code_tree.size());
        assert(index_pair.second < zip_code_tree.size());
        // The pair indexes should be for the same kind of thing
        if (zip_code_tree[index_pair.first].get_type() == CHAIN_START) {
            assert(zip_code_tree[index_pair.second].get_type() == CHAIN_END);
        } else if (zip_code_tree[index_pair.first].get_type() == CHAIN_END) {
            assert(zip_code_tree[index_pair.second].get_type() == CHAIN_START);
        } else if (zip_code_tree[index_pair.first].get_type() == DAG_SNARL_START) {
            assert(zip_code_tree[index_pair.second].get_type() == DAG_SNARL_END);
        } else if (zip_code_tree[index_pair.first].get_type() == DAG_SNARL_END) {
            assert(zip_code_tree[index_pair.second].get_type() == DAG_SNARL_START);
        } else if (zip_code_tree[index_pair.first].get_type() == CYCLIC_SNARL_START) {
            assert(zip_code_tree[index_pair.second].get_type() == CYCLIC_SNARL_END);
        } else if (zip_code_tree[index_pair.first].get_type() == CYCLIC_SNARL_END) {
            assert(zip_code_tree[index_pair.second].get_type() == CYCLIC_SNARL_START);
        } else {
            assert(false);
        }
    }
}

void ZipCodeTree::validate_zip_tree_order(const SnarlDistanceIndex& distance_index, 
                                          const vector<Seed>* seeds) const {
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "Validate tree order" << endl;
#endif
    size_t previous_seed_index = std::numeric_limits<size_t>::max();
    bool previous_is_invalid = false;
    for (size_t i = 0 ; i < zip_code_tree.size() ; i++) {
        const tree_item_t& current_item = zip_code_tree[i];
        if (current_item.get_type() == SEED) {
            // Check if this is worth validating
            // Use a distance limit of 0 so it will ignore looping chains
            bool current_is_invalid = node_is_invalid(id(seeds->at(current_item.get_value()).pos), distance_index, 0);

            if (previous_seed_index != std::numeric_limits<size_t>::max() &&
                !current_is_invalid && !previous_is_invalid) {
                assert(previous_seed_index < seeds->size());
                assert(current_item.get_value() < seeds->size());
#ifdef DEBUG_ZIP_CODE_TREE
                cerr << "Comparing seeds " << seeds->at(previous_seed_index).pos << " and " 
                        << seeds->at(current_item.get_value()).pos << endl;
#endif

                // Comparator returning previous_seed_index < current_item.value
                size_t depth = 0;

                // Keep track of the orientation of each seed
                // Everything is sorted by orientation in top-level structure,
                // so if things are traversed backwards, reverse the orientation
                bool a_is_reversed = false;
                bool b_is_reversed = false;
                const Seed& previous_seed = seeds->at(previous_seed_index);
                const Seed& current_seed = seeds->at(current_item.get_value());
                while (depth < previous_seed.zipcode.max_depth()
                       && depth < current_seed.zipcode.max_depth()
                       && ZipCode::is_equal(previous_seed.zipcode, current_seed.zipcode, depth)) {
                    // Remember the orientation
                    if (ZipCodeTree::seed_is_reversed_at_depth(previous_seed, depth, distance_index)) { 
                        a_is_reversed = !a_is_reversed;
                    }
                    if (ZipCodeTree::seed_is_reversed_at_depth(current_seed, depth, distance_index)) {
                        b_is_reversed = !b_is_reversed;
                    }

                    depth++;
                }

                // Remember the orientation of the parent too
                size_t parent_of_a_is_reversed = a_is_reversed;

                // Check the orientations one last time
                if (ZipCodeTree::seed_is_reversed_at_depth(previous_seed, depth, distance_index)) { 
                    a_is_reversed = !a_is_reversed;
                }
                if (ZipCodeTree::seed_is_reversed_at_depth(current_seed, depth, distance_index)) {
                    b_is_reversed = !b_is_reversed;
                }
                
#ifdef DEBUG_ZIP_CODE_TREE
                cerr << "\t different at depth " << depth << endl;
#endif
                // Either depth is last thing in previous_seed or current_seed
                // or they are different at this depth

                if ( ZipCode::is_equal(previous_seed.zipcode, current_seed.zipcode, depth)) {
#ifdef DEBUG_ZIP_CODE_TREE
                    cerr << "\tthey are on the same node" << endl;
#endif
                    // If they are equal, then they must be on the same node

                    size_t offset1 = is_rev(previous_seed.pos)
                                      ? previous_seed.zipcode.get_length(depth) - offset(previous_seed.pos)
                                      : offset(previous_seed.pos);
                    size_t offset2 = is_rev(current_seed.pos) 
                                     ? current_seed.zipcode.get_length(depth) - offset(current_seed.pos)
                                     : offset(current_seed.pos);
                    if (!a_is_reversed) {
                        // If they are in previous_seed_index snarl or
                        // facing forward on a chain, order by offset in node
                        assert(offset1 <= offset2);
                    } else {
                        // Otherwise, the node is facing backwards in the chain, 
                        // so order backwards in node
                        assert(offset2 <= offset1);
                    }
                }  else if (depth == 0) {
#ifdef DEBUG_ZIP_CODE_TREE
                    cerr << "\tThey are on different connected components" << endl;
#endif
                    // If they are on different connected components, 
                    // sort by connected component
                    assert(previous_seed.zipcode.get_distance_index_address(0)
                           <= current_seed.zipcode.get_distance_index_address(0));
                }  else if (previous_seed.zipcode.get_code_type(depth-1) == ZipCode::CHAIN 
                            || previous_seed.zipcode.get_code_type(depth-1) == ZipCode::ROOT_CHAIN) {
#ifdef DEBUG_ZIP_CODE_TREE
                    cerr << "\t they are children of a common chain" << endl;
#endif
                    // If previous_seed and current_seed are children of a chain
                    size_t component_a = previous_seed.zipcode.get_chain_component(depth);
                    size_t component_b = current_seed.zipcode.get_chain_component(depth);
                    size_t offset_a = previous_seed.zipcode.get_offset_in_chain(depth);
                    size_t offset_b = current_seed.zipcode.get_offset_in_chain(depth);
                    if (component_a == component_b) {
                        if (offset_a == offset_b) {
                            // If they have the same prefix sum,
                            // then the snarl comes first
                            // Will never be on the same child at this depth
                            if (parent_of_a_is_reversed) {
                                assert(current_seed.zipcode.get_code_type(depth) != ZipCode::NODE);
                                assert(previous_seed.zipcode.get_code_type(depth) == ZipCode::NODE); 
                            } else {
                                assert(previous_seed.zipcode.get_code_type(depth) != ZipCode::NODE);
                                assert(current_seed.zipcode.get_code_type(depth) == ZipCode::NODE); 
                            }
                        } else {
                            // Check if the parent chain is reversed
                            // and if so, then the order should be reversed
                            // Could happen in irregular snarls
                            if (parent_of_a_is_reversed) {
                                assert(offset_b <= offset_a);
                            } else {
                                assert(offset_a <= offset_b);
                            }
                        }
                    } else {
                        if (parent_of_a_is_reversed) {
                            assert(component_b <= component_a);
                        } else {
                            assert(component_a <= component_b);
                        }
                    }
                } else if (previous_seed.zipcode.get_code_type(depth-1) == ZipCode::REGULAR_SNARL
                            || previous_seed.zipcode.get_code_type(depth-1) == ZipCode::IRREGULAR_SNARL) {
#ifdef DEBUG_ZIP_CODE_TREE
                    cerr << "\t they are children of a common dag snarl" << endl;
#endif
                    // Otherwise, they are children of a snarl
                    // Sort by a topological ordering from snarl start
                    // Ranks of children in snarls are in a topological order, 
                    // so sort on the ranks
                    assert(previous_seed.zipcode.get_rank_in_snarl(depth)
                            <= current_seed.zipcode.get_rank_in_snarl(depth));
                } 

            }
            previous_seed_index = current_item.get_value();
            previous_is_invalid = current_is_invalid;
        } else if (current_item.get_type() == ZipCodeTree::CHAIN_START) {
            // Chains can't start with edges
            assert(zip_code_tree[i+1].get_type() != ZipCodeTree::EDGE);
        } else if (current_item.get_type() == ZipCodeTree::CHAIN_END) {
            // And can't end with edges
            assert(zip_code_tree[i-1].get_type() != ZipCodeTree::EDGE);
        } else if (current_item.is_snarl_start()) {
            if (i != 0) {
                // Non-root snarls start with their node counts
                assert(zip_code_tree[i+1].get_type() == ZipCodeTree::CHAIN_COUNT);
            }
        }
    }
}

void ZipCodeTree::validate_seed_distances(const SnarlDistanceIndex& distance_index, 
                                          const vector<Seed>* seeds,
                                          size_t distance_limit) const {
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "Validate distances between seeds via iterator" << endl;
#endif
#ifdef debug_parse
    std::cerr << "Validating tree:";
    print_self(seeds);
#endif
    // Walk from the start of the zip tree, checking each pair of seeds
    ZipCodeTree::seed_iterator dest = begin();
    bool right_to_left = true;
    while (dest != end()) {
        // The minimum distance is the true minimum so we make sure,
        // for each other seed, at least one of the outputted distances works
        unordered_set<oriented_seed_t> passing_seeds;
        // seed -> (tree_distance, index_distance)
        unordered_map<oriented_seed_t, pair<size_t, size_t>> failing_seeds;

        // The first seed that the iterator points to
        // While there might be multiple, they all have the same position
        // so verifying distances for one of them is enough
        const Seed& start_seed = seeds->at((*dest).front().seed);
        // Which direction do we traverse the seed?
        const bool start_is_reversed = right_to_left ? (*dest).front().is_reversed
                                                     : !(*dest).front().is_reversed;
#ifdef debug_parse
        std::cerr << "-----------------------------------------" << std::endl;
        std::cerr << start_seed.pos << ": " << (right_to_left ? "Right to left" : "Left to right") << std::endl;
#endif

        // Walk through the tree starting from dest and check the distance
        for (auto dist_itr = find_distances(dest, right_to_left, distance_limit); !dist_itr.done(); ++dist_itr) {
            // While there might be multiple seeds at this position,
            // we only need to check the first one
            seed_result_t next_seed_result = (*dist_itr).front();
            // The seed we reached in our walk
            const Seed& next_seed = seeds->at(next_seed_result.seed);
            const bool next_is_reversed = next_seed_result.is_reversed;
            size_t tree_distance = next_seed_result.distance;

            // Pull metadata about these seeds
            net_handle_t start_handle = distance_index.get_node_net_handle(id(start_seed.pos));
            net_handle_t next_handle = distance_index.get_node_net_handle(id(next_seed.pos));
            size_t start_length = distance_index.minimum_length(start_handle);
            size_t next_length = distance_index.minimum_length(next_handle);

            // Reverse seed positions if we traversed them backwards
            pos_t start_pos = start_is_reversed ? reverse(start_seed.pos, start_length)
                                                : start_seed.pos;
            pos_t next_pos = next_is_reversed ? reverse(next_seed.pos, next_length)
                                              : next_seed.pos;

            // Calculate orientated distance next_pos -> start_pos
            size_t index_distance = minimum_nontrivial_distance(distance_index, next_pos, start_pos);

            if (index_distance != std::numeric_limits<size_t>::max() && (*dest).front().is_reversed == right_to_left) {
                // If the seed we ended at got reversed, then add 1
                index_distance += 1;
            }
            if (index_distance != std::numeric_limits<size_t>::max() && next_seed_result.is_reversed) {
                // If the seed we're starting from got reversed, then subtract 1
                index_distance -= 1;
            }

            bool distance_is_invalid = node_is_invalid(id(next_seed.pos), distance_index, distance_limit)
                                       ||  node_is_invalid(id(start_seed.pos), distance_index, distance_limit);
            if (!distance_is_invalid && (index_distance <= distance_limit 
                                        || (tree_distance <= distance_limit && tree_distance < index_distance))) {
                if (tree_distance != index_distance) {
#ifdef debug_parse
                    cerr << "\tWarning: distance mismatch found" << endl;
                    cerr << "Distance between " << next_seed.pos << (next_is_reversed ? "rev" : "") 
                            << " and " << start_seed.pos << (start_is_reversed ? "rev" : "") << endl;
                    cerr << "Tree distance: " << tree_distance << " index distance: " << index_distance << endl;
                    cerr << "With distance limit: " << distance_limit << endl;
#endif
                    failing_seeds[{next_seed_result.seed, next_is_reversed}] = make_pair(tree_distance, index_distance);
                } else {
                    passing_seeds.insert({next_seed_result.seed, next_is_reversed});
                }
            }
        }

        for (const auto& failure : failing_seeds) {
            // If the same as a pair of seeds that passed,
            // then there are simply multiple paths between these seeds
            if (!passing_seeds.count(failure.first)) {
                for (auto& seed : *seeds) {
                    cerr << seed.pos << endl;
                }

                pos_t start_pos = start_seed.pos;
                pos_t next_pos = seeds->at(failure.first.seed).pos;
                bool next_is_reversed = failure.first.is_reversed;
                size_t tree_distance = failure.second.first;
                size_t index_distance = failure.second.second;

                // This pair of seeds never managed a correct distance
                cerr << "Distance between " << next_pos << (next_is_reversed ? "rev" : "")
                    << " and " << start_pos << (start_is_reversed ? "rev" : "") << endl;
                cerr << "Tree distance: " << tree_distance << " index distance: " << index_distance << endl;
                cerr << "With distance limit: " << distance_limit << endl;
                assert(false);
            }
        }

        // Check next seed/direction along ziptree
        if (dest.in_cyclic_snarl()) {
            if (right_to_left) {
                // Seeds in cyclic snarls are checked in both directions
                right_to_left = false;
            } else {
                right_to_left = true;
                // After going in both directions, we can advance
                ++dest;
            }
        } else {
            // Seeds in non-cyclic snarls are only checked in one direction.
            ++dest;
        }
    }
}

void ZipCodeTree::validate_zip_tree(const SnarlDistanceIndex& distance_index, 
                                    const vector<Seed>* seeds,
                                    size_t distance_limit) const {
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "Validate tree with distance limit " << distance_limit << endl;
#endif
    assert(zip_code_tree.size() != 0);

    // Quickly verify memory of snarl start indexes
    for (const auto& snarl_start : snarl_start_indexes) {
#ifdef DEBUG_ZIP_CODE_TREE
        cerr << "Checking snarl ID " << snarl_start.first << " at index " << snarl_start.second << endl;
#endif
        assert(snarl_start.second < zip_code_tree.size());
        assert(zip_code_tree[snarl_start.second].is_snarl_start());
        assert(zip_code_tree[snarl_start.second].get_value() == snarl_start.first);
    }

    validate_boundaries(distance_index, seeds, distance_limit);
    validate_zip_tree_order(distance_index, seeds);
    validate_seed_distances(distance_index, seeds, distance_limit);
}

void ZipCodeForest::validate_zip_forest(const SnarlDistanceIndex& distance_index, 
                                        const vector<Seed>* seeds, size_t distance_limit) const {
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "Validating zip code forest with distance limit " << distance_limit << endl;
    print_self(seeds);
#endif
    vector<bool> has_seed(seeds->size(), false);
    for (const auto& tree : trees) {
        tree.validate_zip_tree(distance_index, seeds, distance_limit);
        for (const auto seed : tree.get_all_seeds()) {
            has_seed[seed.seed] = true;
        }
    }

    for (size_t i = 0 ; i < has_seed.size() ; i++) {
        bool x = has_seed[i];
        if (!x) { cerr << "Missing seed " << seeds->at(i).pos << endl;}
        assert(x);
    }
}

void ZipCodeTree::validate_snarl(std::vector<tree_item_t>::const_iterator& zip_iterator, 
                                 const SnarlDistanceIndex& distance_index, 
                                 const vector<Seed>* seeds,
                                 size_t distance_limit) const {
#ifdef DEBUG_ZIP_CODE_TREE
    std::cerr << "Validating snarl" << std::endl;
#endif
    // Snarl header
    assert(zip_iterator->is_snarl_start());
    bool is_cyclic_snarl = (zip_iterator->get_type() == ZipCodeTree::CYCLIC_SNARL_START);
    
    size_t snarl_id = zip_iterator->get_value();
    zip_iterator++;
    assert(zip_iterator->get_type() == ZipCodeTree::CHAIN_COUNT);
    // Nodes are children plus the snarl start
    size_t node_count = zip_iterator->get_value()+1;
    zip_iterator++;

    // Read distance matrix - can't use until we read the chains
    vector<size_t> dist_matrix;
    while (zip_iterator->get_type() == ZipCodeTree::EDGE) {
        dist_matrix.push_back(zip_iterator->get_value());
        zip_iterator++;
    }

    assert(dist_matrix.size() == is_cyclic_snarl ? (node_count*2) * (node_count*2 + 1) / 2
                                                 : node_count * (node_count + 1) / 2);

    size_t chains_seen = 0;
    while (!zip_iterator->is_snarl_end()) {
        // Check that this child chain is valid
        assert(zip_iterator->get_type() == ZipCodeTree::CHAIN_START);
        assert(zip_iterator->get_value() == snarl_id);

        // Validate nested chain
        validate_chain(zip_iterator, distance_index, seeds, distance_limit);
        chains_seen++;

        // Skip to next chain start
        zip_iterator++;
    }
    // Verify snarl end bound
    assert(zip_iterator->get_type() == is_cyclic_snarl ? ZipCodeTree::CYCLIC_SNARL_END 
                                                       : ZipCodeTree::DAG_SNARL_END);
    assert(zip_iterator->get_value() == snarl_id);

    // Was the CHAIN_COUNT accurate?
    assert(node_count == chains_seen + 1);
}

void ZipCodeTree::validate_chain(vector<tree_item_t>::const_iterator& zip_iterator, 
                                  const SnarlDistanceIndex& distance_index, 
                                  const vector<Seed>* seeds,
                                  size_t distance_limit) const {
#ifdef DEBUG_ZIP_CODE_TREE
    std::cerr << "Validating chain" << std::endl;
#endif
    // Remember the last seed we saw
    size_t last_seed_i = std::numeric_limits<size_t>::max();
    // Skip past chain start
    zip_iterator++;

    while (zip_iterator->get_type() != ZipCodeTree::CHAIN_END) {
        if (zip_iterator->get_type() == SEED) {
            last_seed_i = zip_iterator->get_value();
        } else if (zip_iterator->is_snarl_start()) {
            // Validate nested snarl
            validate_snarl(zip_iterator, distance_index, seeds, distance_limit);
        }
        zip_iterator++;
    }
}

ZipCodeTree::seed_iterator::seed_iterator(size_t start_index, const ZipCodeTree& ziptree)
    : index(start_index), zip_code_tree(ziptree.zip_code_tree),
    snarl_start_indexes(ziptree.snarl_start_indexes), bound_pair_indexes(ziptree.bound_pair_indexes),
    cyclic_snarl_nestedness(0), chain_numbers(std::stack<size_t>()) {
    
    // If we begin on a snarl, remember that before incrementing
    if (current_item().is_snarl_start()) {
        chain_numbers.push(0);
    }
    // Immediately advance to the first seed
    while (index < zip_code_tree.size() && current_item().get_type() != SEED) {
        ++(*this);
    }
}

auto ZipCodeTree::seed_iterator::operator++() -> seed_iterator& {
    ++index;
    while (index < zip_code_tree.size() && current_item().get_type() != SEED) {
        // cyclic_snarl_nestedness remembers if we're in a cyclic snarl
        if (current_item().get_type() == ZipCodeTree::CYCLIC_SNARL_START) {
            cyclic_snarl_nestedness++;
        } else if (current_item().get_type() == ZipCodeTree::CYCLIC_SNARL_END) {
            cyclic_snarl_nestedness--;
        }

        // chain_numbers remembers which chain we're in for each snarl
        if (current_item().is_snarl_start()) {
            // chain_numbers start with 1, so the snarl start is 0
            chain_numbers.push(0);
        } else if (current_item().is_snarl_end()) {
            chain_numbers.pop();
        } else if (current_item().get_type() == ZipCodeTree::CHAIN_START) {
            chain_numbers.top()++;
        }

        // Advance to the next seed, or the end.
        ++index;
    }
    return *this;
}

ZipCodeTree::distance_iterator::distance_iterator(size_t start_index,
                                                  const vector<tree_item_t>& zip_code_tree,
                                                  const unordered_map<size_t, size_t>& snarl_start_indexes,
                                                  const unordered_map<size_t, size_t>& bound_pair_indexes,
                                                  std::stack<size_t> chain_numbers,
                                                  bool right_to_left,
                                                  size_t distance_limit) :
    index(start_index), original_index(start_index), end_index(right_to_left ? zip_code_tree.size() : 0),
    zip_code_tree(zip_code_tree), snarl_start_indexes(snarl_start_indexes), bound_pair_indexes(bound_pair_indexes),
    chain_numbers(chain_numbers), right_to_left(right_to_left), original_right_to_left(right_to_left),
    distance_limit(distance_limit), stack_data(std::stack<size_t>()), current_state(S_START) {
    if (done()) {
        // We are an end iterator. Nothing else to do.
        return;
    }
#ifdef debug_parse
    std::cerr << "Able to do first initial tick." << std::endl;
#endif
    tick();
#ifdef debug_parse
    if (!done()) {
        std::cerr << "Able to do another initial tick." << std::endl;
    }
#endif
    // Skip to the first seed we actually want to yield, or to the end
    ++(*this);
    // As the end of the constructor, the iterator points to
    // a seed that has been ticked and yielded, or is rend.
#ifdef debug_parse
    if (done()) {
        std::cerr << "Tree iteration halted looking for first seed." << std::endl;
    }
#endif
}

auto ZipCodeTree::distance_iterator::operator++() -> distance_iterator& {
    // Invariant: the iterator points to
    // a seed that has been ticked and yielded, or to rend.
    if (!done()) {
#ifdef debug_parse
        std::cerr << "Skipping over a " << current_item().get_type()
                << " which we assume was handled already." << std::endl;
#endif
        move_index();
    }
    while (!done() && !tick()) {
        // Skip to the next seed we actually want to yield, or to the end
        move_index();
    }
#ifdef debug_parse
    if (done()) {
        std::cerr << "Tree iteration halted looking for next seed." << std::endl;
    }
#endif
    return *this;
}

auto ZipCodeTree::distance_iterator::operator*() const -> vector<seed_result_t> {
    // We are always at a seed, so show that seed
#ifdef check_parse
    crash_unless(!done());
    crash_unless(current_item().get_type() == SEED);
    crash_unless(!stack_data.empty());
#endif
    // We know the running distance to this seed will be at the top of the stack
    // If we're at the exact same position, the distance must be 0
    size_t distance = stack_data.top();
    bool is_reversed = right_to_left ? current_item().get_is_reversed()
                                     : !current_item().get_is_reversed();
    vector<seed_result_t> to_return;
    to_return.emplace_back(distance, current_item().get_value(), is_reversed);
    if (current_item().has_other_values()) {
        // If this seed has other values, add them too
        for (const auto& other_value : current_item().get_other_values()) {
            to_return.emplace_back(distance, other_value, is_reversed);
        }
    }
    return to_return;
}

auto ZipCodeTree::distance_iterator::push(size_t value) -> void {
    stack_data.push(value);
}

auto ZipCodeTree::distance_iterator::pop() -> size_t {
    size_t value = stack_data.top();
    stack_data.pop();
    return value;
}

auto ZipCodeTree::distance_iterator::top() -> size_t& {
#ifdef check_parse
    crash_unless(depth() > 0);
#endif
    return stack_data.top();
}

auto ZipCodeTree::distance_iterator::dup() -> void {
    push(stack_data.top());
}

auto ZipCodeTree::distance_iterator::depth() const -> size_t {
    return stack_data.size();
}

auto ZipCodeTree::distance_iterator::swap() -> void {
    // Grab the top item
    size_t temp = stack_data.top();
    stack_data.pop();
    // Swap it with what was under it
    std::swap(temp, stack_data.top());
    // And put that back on top
    stack_data.push(temp);
}

vector<size_t> ZipCodeTree::distance_iterator::get_distances_from_chain(size_t snarl_id, size_t chain_num, 
                                                                        bool right_side) const {
    // Read snarl header
    bool is_cyclic = snarl_is_cyclic(snarl_id);
    // SNARL_START, then CHAIN_COUNT, then distance matrix
    size_t dist_matrix_start = snarl_start_indexes.at(snarl_id) + 2;
    tree_item_t chain_count = zip_code_tree[dist_matrix_start - 1];
    if (chain_count.get_type() == ZipCodeTree::CHAIN_START) {
        // This must be a root snarl, without a distance matrix
        return vector<size_t>();
    }
    size_t num_chains = chain_count.get_value();
#ifdef debug_parse
    cerr << "Get distances for snarl " << snarl_id << " with " << num_chains << " chain(s); "
         << "stacking for chain " << chain_num << "'s "
         << (right_side ? "right" : "left") << " side" << endl;
    assert(chain_num <= num_chains || chain_num == std::numeric_limits<size_t>::max());
#endif

    // Distances in bottom-to-top order (i.e. last will be top)
    vector<size_t> distances;

    if (is_cyclic) {
        // Space for either side of each chain, plus snarl bound
        distances.reserve(num_chains * 2 + 1);
        // Which row corresponds to this chain side?
        size_t cur_row; 
        if (chain_num == 0) {
            // Chain 0 means the snarl start
            cur_row = 0;
        } else if (chain_num == std::numeric_limits<size_t>::max()) {
            // Chain max means the snarl end
            cur_row = num_chains * 2 + 1;
        } else {
            // Otherwise, we are in a chain
            cur_row = chain_num * 2 - (right_side ? 0 : 1);
        }

        // Very bottom of the stack is the distance to a snarl bound
        if (right_side) {
            // Heading towards the snarl end
            distances.push_back(get_matrix_value(dist_matrix_start, true, cur_row, num_chains * 2 + 1));
        } else {
            // Heading towards the snarl start
            distances.push_back(get_matrix_value(dist_matrix_start, true, cur_row, 0));
        }

        // Cyclic snarl always stacks all lefts on top, then all rights 
        // We process left to right then "bounce" right to left
        // c1_L is on top, and c1_R is on bottom, directly above the snarl bound
        for (size_t i = 1; i <= num_chains; i++) {
            // Right sides from c1 to cN
            distances.push_back(get_matrix_value(dist_matrix_start, true, cur_row, i * 2));
        }
        for (size_t i = num_chains; i >= 1; i--) {
            // Left sides from cN to c1
            distances.push_back(get_matrix_value(dist_matrix_start, true, cur_row, i * 2 - 1));
        }
    } else {
        if (chain_num == std::numeric_limits<size_t>::max()) {
            // Chain max means the snarl end
            chain_num = num_chains + 1;
        }

        if (right_side) {
            // Space for all chains to the right, plus snarl end
            distances.reserve(num_chains + 1 - chain_num);
            // DAG snarl, going left to right: get dists to all chains to right
            for (size_t i = num_chains + 1; i > chain_num; i--) {
                distances.push_back(get_matrix_value(dist_matrix_start, false, chain_num, i));
            }
        } else {
            // Space for all chains to the left, plus snarl start
            distances.reserve(chain_num);
            // DAG snarl, going right to left: get dists to all chains to left
            for (size_t i = 0; i < chain_num; i++) {
                distances.push_back(get_matrix_value(dist_matrix_start, false, chain_num, i));
            }
        }
    }

#ifdef debug_parse
    cerr << "Distances (bottom to top): ";
    for (const auto& dist : distances) {
        cerr << dist << " ";
    }
    cerr << endl;
#endif
    return distances;
}

size_t ZipCodeTree::distance_iterator::get_matrix_value(size_t matrix_start_i, bool has_main_diagonal, 
                                                        size_t row, size_t col) const {
    if (row < col) {
        // We only store the lower triangle, so swap the row and column
        size_t temp = row;
        row = col;
        col = temp;
    }

    // Triangular number of elements in previous rows, then offset by col
    size_t within_matrix_i = has_main_diagonal ? (row * (row + 1)) / 2 + col
                                               : (row * (row - 1)) / 2 + col;
    return zip_code_tree[matrix_start_i + within_matrix_i].get_value();
}

void ZipCodeTree::distance_iterator::skip_chain() {
    // Discard distance for this chain
    pop();
    // Jump to other bound
    index = pop();
#ifdef debug_parse
    std::cerr << "Skip chain, jump to index " << index << std::endl;
#endif
    // We skipped because we want to try other chains in this snarl
    continue_snarl();
}

void ZipCodeTree::distance_iterator::initialize_chain() {
    // Where *would* we jump, if we jumped?
    push(bound_pair_indexes.at(index));
    swap();
    if (top() > distance_limit || top() == std::numeric_limits<size_t>::max()) {
#ifdef debug_parse
        std::cerr << "Skip chain with running distance "
                << (top() == std::numeric_limits<size_t>::max() ? "inf" : std::to_string(top())) 
                << " over limit " << distance_limit << std::endl;
#endif
        // Running distance is already too high so skip the chain
        skip_chain();
    } else {
        // Do the chain
        state(S_SCAN_CHAIN);
    }
}

bool ZipCodeTree::distance_iterator::initialize_snarl(size_t chain_num) {
#ifdef debug_parse
    std::cerr << "Initialize snarl from chain number " << chain_num 
              << " going " << (right_to_left ? "right to left" : "left to right")
              << " with running distance " << top() << std::endl;
#endif
    // Grab distances for this snarl
    size_t snarl_id = current_item().get_value();
    bool is_cyclic = snarl_is_cyclic(snarl_id);
    vector<size_t> distances = get_distances_from_chain(snarl_id, chain_num, !right_to_left);
    if (distances.empty()) {
#ifdef debug_parse
        std::cerr << "Tried to stack up distances for a root-level snarl; halting now." << std::endl;
#endif
        halt();
        return false;
    }

    if (is_cyclic) {
        // Memorize previous direction
        push(right_to_left ? 1 : 0);
        // Put previous direction underneath the parent running distance
        swap();
        // Cyclic snarls always go left to right first
        right_to_left = false;
    }

    // Add distances to running distance
    for (const auto& distance : distances) {
        dup();
        // Add in the edge value to make a running distance for child
#ifdef debug_parse
        std::cerr << "\tTo current running parent distance " << top() << " add edge " 
                << (distance == std::numeric_limits<size_t>::max() ? "inf" : std::to_string(distance))
                << std::endl;
#endif
        top() = SnarlDistanceIndex::sum(top(), distance);
        // Flip top 2 elements, so now parent running distance is on top,
        // over edge running distance.
        swap();
    }
    // Remove parent running distance
    pop();
    continue_snarl();
    
    if (is_cyclic || chain_num == 0) {
#ifdef debug_parse
        std::cerr << "Start with leftmost chain" << std::endl;
#endif
        // We want to start with the leftmost chain
        // Jump to CHAIN_COUNT
        index = snarl_start_indexes.at(snarl_id) + 1;
        // Distance matrix has chains & bounds
        size_t node_count = current_item().get_value() + 1;
        // Skip past distance matrix
        index += is_cyclic ? (node_count*2 * (node_count*2 + 1)) / 2
                           : (node_count * (node_count + 1)) / 2;
#ifdef debug_parse
        std::cerr << "Jump to index " << index << " after distance matrix" << std::endl;
#endif
    }
    return true;
}

void ZipCodeTree::distance_iterator::continue_snarl() {
    // Different scanning states based on snarl type
    if (snarl_is_cyclic(current_item().get_value())) {
#ifdef debug_parse
        std::cerr << "Continuing cyclic snarl" << std::endl;
#endif
        state(S_SCAN_CYCLIC_SNARL);
    } else {
#ifdef debug_parse
        std::cerr << "Continuing DAG snarl" << std::endl;
#endif
        state(S_SCAN_DAG_SNARL);
    }
}

auto ZipCodeTree::distance_iterator::halt() -> void {
#ifdef debug_parse
    std::cerr << "Halt iteration!" << std::endl;
#endif
    index = end_index;
}

void ZipCodeTree::distance_iterator::unimplemented_error() {
    throw std::domain_error("Unimplemented symbol " + std::to_string(current_item().get_type()) 
                                    + " for state " + std::to_string(current_state));
}

auto ZipCodeTree::distance_iterator::tick() -> bool {
#ifdef debug_parse
    std::cerr << "Tick for state " << current_state << " on symbol " << current_item().get_type() 
              << " at entry " << index << std::endl;
    std::cerr << "Stack: ";
    std::stack<size_t> temp_stack = stack_data;
    while (!temp_stack.empty()) {
        std::cerr << (temp_stack.top() == std::numeric_limits<size_t>::max() ? "inf" 
                                                                             : std::to_string(temp_stack.top())) 
                  << " ";
        temp_stack.pop();
    }
    std::cerr << std::endl;
#endif
    switch (current_state) {
    case S_START:
        // Initial state.
        //
        // Stack is empty and we must be at a seed to start at.
        if (current_item().get_type() == SEED) {
#ifdef debug_parse
            std::cerr << "Skip over seed " << current_item().get_value() << std::endl;
#endif
            push(0);
            state(S_SCAN_CHAIN);
        } else {
            unimplemented_error(); 
        }
        break;
    case S_SCAN_CHAIN:
        // State where we are scanning through a chain
        //
        // Stack has at the top the running distance along the chain, and under
        // that running distances to use at the other chains in the snarl, and
        // under that running distances to use for the other chains in the
        // snarl's parent snarl, etc.
        if (current_item().get_type() == SEED) {
            // Calculate which direction this seed would be traversed in
            bool cur_is_rev = right_to_left ? current_item().get_is_reversed() 
                                            : !current_item().get_is_reversed();
            bool origin_is_rev = original_right_to_left ? zip_code_tree.at(original_index).get_is_reversed()
                                                        : !zip_code_tree.at(original_index).get_is_reversed();
            if (cur_is_rev == origin_is_rev) {
                // Emit seed here with distance at top of stack.
#ifdef check_parse
                crash_unless(depth() > 0);
#endif
#ifdef debug_parse
                std::cerr << "Yield seed " << current_item().get_value() << ", distance " << top() << std::endl;
#endif
                return true;
            } else {
                // This seed is in the opposite read orientation,
                // thus it could not be used anyhow: we skip it
#ifdef debug_parse
                std::cerr << "Skip seed " << current_item().get_value()
                          << " with incompatible read orientation" << std::endl;
#endif
            }
        } else if (entered_snarl()) {
            // Running distance along chain is on stack,
            // and will need to be added to all the stored distances.
            return !initialize_snarl(current_item().is_snarl_start() ? 0 
                                                                     : std::numeric_limits<size_t>::max());
        } else if (exited_chain()) {
            if (depth() == 1) {
                // We never entered the parent snarl of this chain, so stack up
                // the distances left of here as options added to the
                // distance along this chain.
                //
                // Running distance along chain is on stack, and will need to
                // be added to all the stored distances.
                // Note that there may be 0 stored distances
                // if we are below the top-level snarl.
                if (chain_numbers.empty()) {
                    // Finished top-level chain, so we are done.
#ifdef debug_parse
                    std::cerr << "Halt because we are at the top-level chain." << std::endl;
#endif
                    halt();
                    // When we halt we have to return true to show the position.
                    return true;
                }
#ifdef debug_parse
                std::cerr << "\tNever entered snarl; collect distances from zipcode tree." << std::endl;
#endif
                size_t chain_num = chain_numbers.top();
                chain_numbers.pop();
                return !initialize_snarl(chain_num);
            } else {
                // We did enter the parent snarl already.
                // Discard the running distance along this chain,
                // which no longer matters.
                pop();
                // Discard the jump index for this chain
                pop();
                // Running distance for next chain,
                // or running distance to cross the snarl, will be under it.
                continue_snarl();
            }
        } else if (current_item().get_type() == EDGE) {
            // Add value to running distance
            top() = SnarlDistanceIndex::sum(top(), current_item().get_value());
            if (top() > distance_limit || top() == std::numeric_limits<size_t>::max()) {
                // Skip over the rest of this chain
                if (depth() == 1) {
                    // We never entered the parent snarl of this chain.
                    // So if the distance along the chain is too much, there
                    // are not going to be any results with a smaller distance.
#ifdef debug_parse
                    std::cerr << "Halt because adding " << current_item().get_value() << " bp " 
                            << "gives distance " << top() << " > " << distance_limit << std::endl;
#endif
                    halt();
                    // When we halt we have to return true to show the position.
                    return true;
                } else {
                    // We need to try the next thing in the parent snarl,
                    // so skip the rest of the chain.
                    skip_chain();
                }
            }
        } else {
            unimplemented_error(); 
        }
        break;
    case S_SCAN_DAG_SNARL:
        // State where we are going through a DAG snarl having already
        // stacked up the running distances to use for each chain in the snarl.
        //
        // Stack has at the top running distances to use for each chain still
        // to be visited, and under those the same for the snarl above, etc.
        if (current_item().get_type() == EDGE) {
            // Hit the distance matrix; finished snarl
            // Top of the stack should be the parent running distance
            // Jump to the start of the snarl
            index = snarl_start_indexes.at(
                current_item(1).get_value() // Use snarl ID from previous chain
            );
            // Resume scanning chain
            state(S_SCAN_CHAIN);
        } else if (exited_snarl()) {
            // Stack holds running distance to use for the parent chain
            state(S_SCAN_CHAIN);
        } else if (entered_chain()) {
            // We've encountered a chain to look at, and the running distance
            // into the chain is already on the stack.
            initialize_chain();
        } else {
            unimplemented_error(); 
        }
        break;
    case S_SCAN_CYCLIC_SNARL:
        // State where we are going through a cyclic snarl having already
        // stacked up the running distances to use for each chain in the snarl.
        //
        // Cyclic snarls are traversed by "bouncing": entering each chain
        // from the left side, starting with the leftmost, and then once we
        // hit the other end of the snarl, entering each chain from the right.
        if (current_item().get_type() == CYCLIC_SNARL_END) {
            // Finished left sides, now doing right sides.
#ifdef check_parse
            crash_unless(!right_to_left);
#endif
            right_to_left = true;
        } else if (current_item().get_type() == EDGE) {
            // Hit distance matrix again; finished snarl
            // Put original direction over distance to bound
            swap();
#ifdef debug_parse
            std::cerr << "Finished cyclic snarl, resetting direction to "
                    << (top() == 1 ? "right to left" : "left to right" ) << std::endl;
#endif
            // The top of the stack will be the original direction we were going
            right_to_left = (pop() == 1);
            // Use snarl ID from previous chain
            size_t snarl_id = current_item(1).get_value();
            size_t snarl_start_i = snarl_start_indexes.at(snarl_id);
            if (right_to_left) {
                // Jump to snarl start
                index = snarl_start_i;
            } else {
                // Jump to snarl end
                index = bound_pair_indexes.at(snarl_start_i);
            }
#ifdef debug_parse
            std::cerr << "Jump to index " << index << endl;
#endif
            state(S_SCAN_CHAIN);
        } else if (entered_chain()) {
            // We've encountered a chain to look at, and the running distance
            // into the chain is already on the stack.
            initialize_chain();
        } else {
            unimplemented_error();
        }
        break;
    default:
        throw std::domain_error("Unimplemented state " + std::to_string(current_state)); 
    }
    // Unless we yield something, we don't want to pause the scan here.
    return false;
}

auto ZipCodeTree::find_distances(const seed_iterator& from, bool right_to_left,
                                 size_t distance_limit) const -> distance_iterator {
    return distance_iterator(from.get_index(), from.zip_code_tree, from.snarl_start_indexes,
                             from.bound_pair_indexes, from.chain_numbers, right_to_left, distance_limit);
}

std::ostream& operator<<(std::ostream& out, const ZipCodeTree::tree_item_type_t& type) {
    return out << std::to_string(type);
}
std::ostream& operator<<(std::ostream& out, const ZipCodeTree::distance_iterator::State& state) {
    return out << std::to_string(state);
}

void ZipCodeForest::sort_one_interval(forest_growing_state_t& forest_state, const interval_state_t& interval) const {
    vector<size_t>& zipcode_sort_order = forest_state.seed_sort_order; 
    vector<sort_value_t>& sort_values_by_seed = forest_state.sort_values_by_seed;
    const vector<Seed>* seeds = forest_state.seeds;

#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "Sort interval at depth " << interval.depth << (interval.is_reversed ? " reversed" : "") << endl;
#endif

    /*** First, fill in sort_values_by_seed for the relevant seeds ***/

    // This doesn't take into account the orientation,
    // except for nodes offsets in chains
    // Used for sorting at the given depth, so use values at depth depth+1

    // Get the minimum and maximum values that are used for sorting.
    // These will be used to determine if radix sort will be more efficient
    
    // This must be done even if the interval is already sorted,
    // because we need to fill in the sort values

    size_t max_sort_value = 0;
    size_t min_sort_value = std::numeric_limits<size_t>::max();

    // The min and max chain components
    size_t min_component = 0;
    size_t max_component = 0;

    for (size_t i = interval.interval_start ; i < interval.interval_end ; i++) {
        const Seed& seed = seeds->at(zipcode_sort_order[i]); 
#ifdef DEBUG_ZIP_CODE_SORTING
        cerr << "\tGet the sort value of seed " << seed.pos << " at depth " << interval.depth+1 
                << " with parent type " << interval.code_type << endl;
#endif
        if (interval.code_type == ZipCode::EMPTY) {
            // If we are sorting the root into connected components 
#ifdef DEBUG_ZIP_CODE_SORTING
            cerr << "\t\tThis is the root snarl so sort by connected component: " 
                    << seed.zipcode.get_distance_index_address(0) << endl;
#endif
            sort_values_by_seed[zipcode_sort_order[i]].set_sort_value( seed.zipcode.get_distance_index_address(0));
            sort_values_by_seed[zipcode_sort_order[i]].set_code_type(seed.zipcode.get_code_type(0));
        } else if (interval.code_type == ZipCode::NODE || interval.code_type == ZipCode::ROOT_NODE 
                    || seed.zipcode.max_depth() == interval.depth) {
#ifdef DEBUG_ZIP_CODE_SORTING
            cerr << "\t\t this is a node: offset: " << ( is_rev(seed.pos) 
                                                            ? seed.zipcode.get_length(interval.depth) - offset(seed.pos)
                                                            : offset(seed.pos)) << endl;;
#endif
            sort_values_by_seed[zipcode_sort_order[i]].set_sort_value(
                       is_rev(seed.pos) ? seed.zipcode.get_length(interval.depth) - offset(seed.pos)
                                                             : offset(seed.pos));
            sort_values_by_seed[zipcode_sort_order[i]].set_code_type(ZipCode::NODE);
        } else if (interval.code_type == ZipCode::CHAIN || interval.code_type == ZipCode::ROOT_CHAIN) {
#ifdef DEBUG_ZIP_CODE_SORTING
            cerr << "\t\t this is a chain:";
#endif
            // Get the prefix sum and chain order of the chain child.
            // The chain order is the value added to the prefix
            // sum to specify the order of children with the same prefix sum.
            // 1 will be added to snarls, and 2 will be added to the node
            // with offset in node of 0 (node 3 if chain is traversed forward)
            // See sort_value_t for more details

            size_t prefix_sum = seed.zipcode.get_offset_in_chain(interval.depth+1);

            ZipCode::code_type_t child_type = seed.zipcode.get_code_type(interval.depth+1);
            sort_values_by_seed[zipcode_sort_order[i]].set_code_type(child_type);
            size_t chain_component = seed.zipcode.get_chain_component(interval.depth+1);
            sort_values_by_seed[zipcode_sort_order[i]].set_chain_component(chain_component);
            min_component = std::min(min_component, chain_component);
            max_component = std::max(max_component, chain_component);

            if (child_type == ZipCode::REGULAR_SNARL 
                || child_type == ZipCode::IRREGULAR_SNARL
                || child_type == ZipCode::CYCLIC_SNARL) { 

                // For a snarl, the order is prefix_sum*3+1
                sort_values_by_seed[zipcode_sort_order[i]].set_sort_value(prefix_sum);
                sort_values_by_seed[zipcode_sort_order[i]].set_chain_order(1);
            } else {
                // Order depends on where the position falls in the node
                bool node_is_rev = seed.zipcode.get_is_reversed_in_parent(interval.depth+1) != is_rev(seed.pos);
                size_t node_offset = node_is_rev ? seed.zipcode.get_length(interval.depth+1) - offset(seed.pos)
                                                 : offset(seed.pos);

                sort_values_by_seed[zipcode_sort_order[i]].set_sort_value(
                    SnarlDistanceIndex::sum(prefix_sum, node_offset));
                sort_values_by_seed[zipcode_sort_order[i]].set_chain_order(node_offset == 0 ? 2 : 0);
            }
#ifdef DEBUG_ZIP_CODE_SORTING
            cerr << "Prefix sum " << sort_values_by_seed[zipcode_sort_order[i]].get_distance_value()
                 << " and sort value " << sort_values_by_seed[zipcode_sort_order[i]].get_sort_value()
                 << " and type " << child_type << endl;
#endif
        } else {
#ifdef DEBUG_ZIP_CODE_SORTING
            cerr << "\tThis is a snarl, so return the rank in the snarl: " 
                << seed.zipcode.get_rank_in_snarl(interval.depth+1) << endl;
#endif
            // The ranks of children in irregular snarls are in
            // a topological order, so sort on the ranks
            // The rank of children in a regular snarl is arbitrary
            // but it doesn't matter anyway
            sort_values_by_seed[zipcode_sort_order[i]].set_sort_value(seed.zipcode.get_rank_in_snarl(interval.depth+1));
            sort_values_by_seed[zipcode_sort_order[i]].set_code_type(seed.zipcode.get_code_type(interval.depth+1));
        }
        min_sort_value = std::min(min_sort_value, sort_values_by_seed[zipcode_sort_order[i]].get_sort_value());
        max_sort_value = std::max(max_sort_value, sort_values_by_seed[zipcode_sort_order[i]].get_sort_value());
    }

    /***** Figure out which sort method we should use ***/

    bool use_radix;
    if (interval.code_type  == ZipCode::ROOT_CHAIN) {
        // If this is a root chain, then use the default sort, 
        // because it's probably too big for radix and we can't tell
        // anyways because we don't store the length of a root-chain
        use_radix = false;
    } else {
        // The cost of default sort is nlog(n)
        size_t default_cost = (interval.interval_end - interval.interval_start) 
                            * std::log2(interval.interval_end - interval.interval_start);
        // The cost of radix sort is linear in the number of distinct values
        size_t radix_cost = max_sort_value - min_sort_value;
        use_radix = radix_cost <= default_cost;
    }

    /**** Sort *********/

    // If this is a multicomponent chain, then sort by component first
    vector<interval_state_t> sub_intervals;
    if (min_component != max_component) {
        sub_intervals.reserve(max_component-min_component);
#ifdef DEBUG_ZIP_CODE_TREE
        cerr << "Sort by chain component" << endl;
#endif
        // Sort by component using radix sort. I doubt that there will be 
        // enough components to make it more efficient to use the default sort
        radix_sort_zipcodes(zipcode_sort_order, sort_values_by_seed, interval, 
                            interval.is_reversed, min_component, max_component, true);

        // Now get the next intervals in sub_intervals
        size_t start = interval.interval_start;
        size_t previous_component = sort_values_by_seed[zipcode_sort_order[start]].get_chain_component();
        for (size_t i = interval.interval_start ; i < interval.interval_end ; i++) {
            size_t current_component = sort_values_by_seed[zipcode_sort_order[i]].get_chain_component();
            if (current_component != previous_component) {
                sub_intervals.emplace_back(interval);
                sub_intervals.back().interval_start = start;
                sub_intervals.back().interval_end = i;
                start = i;
                previous_component = current_component;
            }
        }
        sub_intervals.emplace_back(interval);
        sub_intervals.back().interval_start = start;
        sub_intervals.back().interval_end = interval.interval_end;
    } else {
        // Copy the current interval
        sub_intervals.emplace_back(interval);
    }

    for (auto& sub_interval : sub_intervals) {
        // Snarls are already sorted by a topological order of
        // the orientation of the zip tree, so don't reverse them
        // And don't reverse the sort if that has
        // already been taken into account in the value finding
        bool reverse_order = (sub_interval.code_type == ZipCode::REGULAR_SNARL 
                              || sub_interval.code_type == ZipCode::IRREGULAR_SNARL) 
                             ? false
                             : sub_interval.is_reversed; 
        // Sort the given interval using the value-getter and orientation
        if (use_radix) {
            radix_sort_zipcodes(zipcode_sort_order, sort_values_by_seed, sub_interval, 
                                reverse_order, min_sort_value, max_sort_value);
        } else {
            default_sort_zipcodes(zipcode_sort_order, sort_values_by_seed, sub_interval, reverse_order);
        }
    }
    return;
}

void ZipCodeForest::get_next_intervals(forest_growing_state_t& forest_state, const interval_state_t& interval,
                                       std::forward_list<ZipCodeForest::interval_state_t>& next_intervals) const {
    vector<size_t>& zipcode_sort_order = forest_state.seed_sort_order;
    vector<sort_value_t>& sort_values_by_seed = forest_state.sort_values_by_seed;
    const vector<Seed>* seeds = forest_state.seeds;
    const SnarlDistanceIndex* distance_index = forest_state.distance_index;
    
    // Add new intervals to the front of next_intervals in their sort order.
    // This means that the first interval found gets added to the front of list.
    // insert_itr points to the item in front of where to add the next interval
    // so emplace/insert_after instert_itr, and move it forward after inserting
    std::forward_list<ZipCodeForest::interval_state_t>::iterator insert_itr = next_intervals.before_begin();

    /*********   Check for new intervals of the children ****************/

#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "Finding intervals after sorting at depth " << interval.depth << endl;
#endif
    // After sorting, find runs of equivalent values for new_interval_to_sort
    // Everything gets put in a new interval, even if it has a unique value
    // Since nodes are really just seeds on the same chain,
    // runs of nodes facing the same direction get put together
    // Also need to check the orientation

    // max() is used for the root, when the child's depth should be 0
    size_t child_depth = interval.code_type == ZipCode::EMPTY ? 0 : interval.depth+1;


    if (interval.code_type != ZipCode::EMPTY && 
        seeds->at(zipcode_sort_order[interval.interval_start]).zipcode.max_depth() == interval.depth) {
        // If this is a trivial chain, then return the same interval as a node
#ifdef DEBUG_ZIP_CODE_TREE
        cerr << "\tthis was a trivial chain so just return the same interval as a node" << endl;
#endif
        next_intervals.emplace_after(insert_itr, interval.interval_start, interval.interval_end, 
                                     interval.is_reversed, ZipCode::NODE, child_depth);
        return;
    }


    // These get compared to see if the next seeds is in the same interval
    ZipCode::code_type_t first_type = sort_values_by_seed[zipcode_sort_order[interval.interval_start]].get_code_type();

    // This is only for nodes in chains, since anything on nodes in chains
    // are considered just children of the chain
    bool previous_is_node = first_type == ZipCode::NODE;

    // This only matters if it isn't a node
    size_t previous_sort_value = previous_is_node 
                                ? (ZipCodeTree::seed_is_reversed_at_depth(
                                    seeds->at(zipcode_sort_order[interval.interval_start]), 
                                    child_depth, *distance_index) 
                                    ? 1 : 0)
                                : sort_values_by_seed[zipcode_sort_order[interval.interval_start]].get_sort_value();

    // Start the first interval. 
    // The end value and is_reversed gets set when ending the interval
    next_intervals.emplace_after(insert_itr, interval.interval_start, interval.interval_start, 
                                 interval.is_reversed, first_type, child_depth);
    ++insert_itr;

    for (size_t i = interval.interval_start+1 ; i < interval.interval_end ; i++) {
        // If the current seed is a node and has nothing at depth+1
        // or is different from the previous seed at this depth
        ZipCode::code_type_t current_type = sort_values_by_seed[zipcode_sort_order[i]].get_code_type();
        bool is_node = current_type == ZipCode::NODE;
        size_t sort_value = is_node 
                          ? (ZipCodeTree::seed_is_reversed_at_depth(seeds->at(zipcode_sort_order[i]), 
                                                                    child_depth, *distance_index) ? 1 : 0)
                          : sort_values_by_seed[zipcode_sort_order[i]].get_sort_value();
        bool is_different_from_previous = (is_node != previous_is_node) ? true 
                                                                        : sort_value != previous_sort_value;
        previous_is_node = is_node;
        previous_sort_value = sort_value;

        if (is_different_from_previous) {
            // If at run end, close previous run; add end value and orientation

            insert_itr->interval_end = i;
            insert_itr->is_reversed = ZipCodeTree::seed_is_reversed_at_depth(seeds->at(zipcode_sort_order[i-1]), 
                                                                             child_depth, *distance_index) 
                                    ? !interval.is_reversed
                                    : interval.is_reversed;
            // Open a new run
            next_intervals.emplace_after(insert_itr, i, i, interval.is_reversed, 
                                         is_node ? ZipCode::NODE : current_type, 
                                         child_depth);
            ++insert_itr;
        }
    }

    // Close the last run
    insert_itr->interval_end = interval.interval_end;
    insert_itr->is_reversed = ZipCodeTree::seed_is_reversed_at_depth(
                                seeds->at(zipcode_sort_order[interval.interval_end-1]), child_depth, *distance_index) 
                             ? !interval.is_reversed
                             : interval.is_reversed;
    return;
}

void ZipCodeForest::radix_sort_zipcodes(vector<size_t>& zipcode_sort_order, 
                                        const vector<sort_value_t>& sort_values_by_seed,
                                        const interval_state_t& interval, bool reverse_order,
                                        size_t min_value, size_t max_value, bool sort_by_chain_component) const {
    // Radix sort the interval of zipcode_sort_order in the given interval
#ifdef DEBUG_ZIP_CODE_SORTING
    cerr << "\tradix sort" << endl;
#endif

    // Mostly copied from Jordan Eizenga

    // ount up occurrences of each rank
    std::vector<size_t> counts (max_value-min_value+2, 0);
    for (size_t i  = interval.interval_start ; i < interval.interval_end ; i++) {
        size_t sort_value = sort_by_chain_component ? sort_values_by_seed[zipcode_sort_order[i]].get_chain_component()
                                                    : sort_values_by_seed[zipcode_sort_order[i]].get_sort_value();
#ifdef DEBUG_ZIP_CODE_SORTING
        assert(sort_value >= min_value);
        assert(sort_value <= max_value);
        assert(counts.size() > sort_value - min_value + 1);
#endif
        size_t next_rank = sort_value - min_value + 1;
        ++counts[next_rank];
    }
    
    // Make this a count of the number of things before it
    for (size_t i = 1; i < counts.size(); ++i) {
        counts[i] += counts[i - 1];
    }

    // Get the sorted order
    std::vector<size_t> sorted(interval.interval_end - interval.interval_start);
    for (size_t i = interval.interval_start ; i < interval.interval_end ; i++) {
        size_t sort_value = sort_by_chain_component ? sort_values_by_seed[zipcode_sort_order[i]].get_chain_component()
                                                    : sort_values_by_seed[zipcode_sort_order[i]].get_sort_value();
        size_t rank = sort_value - min_value;
        sorted[counts[rank]++] = zipcode_sort_order[i];
    }
    
    // And place everything in the correct position
    for (size_t i = 0 ; i < sorted.size() ; i++) {
        if (reverse_order) {
            // If this snarl tree node is reversed, then reverse the sort order
            zipcode_sort_order[interval.interval_end - i - 1] = sorted[i];
        } else {
            zipcode_sort_order[i + interval.interval_start] = sorted[i];
        }
    }

}
void ZipCodeForest::default_sort_zipcodes(vector<size_t>& zipcode_sort_order, 
                                          const vector<sort_value_t>& sort_values_by_seed,
                                          const interval_state_t& interval, bool reverse_order) const { 
    // std::sort interval of zipcode_sort_order [interval_start, interval_end]
    
#ifdef DEBUG_ZIP_CODE_SORTING
    cerr << "\tdefault sort between " << interval.interval_start  << " and " << interval.interval_end  << endl;
    cerr << "\tis rev: " << reverse_order << endl;
#endif
    // Sort using std::sort 
    std::sort(zipcode_sort_order.begin() + interval.interval_start, 
              zipcode_sort_order.begin() + interval.interval_end, [&] (size_t a, size_t b) {
        // If this snarl tree node is reversed, then reverse the sort order
        return reverse_order ? sort_values_by_seed[a].get_sort_value() > sort_values_by_seed[b].get_sort_value()
                             : sort_values_by_seed[a].get_sort_value() < sort_values_by_seed[b].get_sort_value();
    });
}

void ZipCodeForest::fill_in_forest(const vector<Seed>& seeds,
                                   const SnarlDistanceIndex& distance_index, size_t distance_limit) {
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "Make a new forest with " << seeds.size() << " seeds with distance limit " << distance_limit << endl;
    for (auto& x : seeds) {
        cerr << x.pos << endl;
    }
    cerr << endl;
#endif
    if (seeds.size() == 0) {
        return;
    }

    /*
    The zip forest is made by sorting the seeds along chains/snarls,
    then adding each seed, snarl/chain boundary, and distance to zip_code_tree.

    Sorting and tree-making are done at the same time,
    in a depth-first traversal of the snarl tree.
    Sorting is done per node in the snarl tree.
    
    Intervals representing ranges of seeds corresponding to snarl tree
    structures are stored in a stack. The algorithm starts with an interval for
    each child of the root snarl. An interval is popped from the stack.
    Any incomplete snarls or chains that the interval is not a child of must be
    completed. Then, the snarl or chain that the interval represents is added to
    the zip tree, along with any relevant distances. Intervals representing the
    children of the snarl or chain are found and added to the stack.
    This repeats until the stack is empty.
    */

    // Start by initializing the state
    // The forest state keeps track of the sort order of seeds,
    // the intervals to be sorted, and which intervals are open and incomplete. 
    forest_growing_state_t forest_state(seeds, distance_index, distance_limit);

    // Start with root as an interval over seed_sort_order containing everything
    interval_state_t first_interval (0, seeds.size(), false, ZipCode::EMPTY, 0);

    // Sort and get the intervals of the connected components
    sort_one_interval(forest_state, first_interval);
    get_next_intervals(forest_state, first_interval, forest_state.intervals_to_process);

    while (!forest_state.intervals_to_process.empty()) {
#ifdef DEBUG_ZIP_CODE_TREE
        print_self(&seeds);
#endif
        // Process each unprocessed interval
        // First, check if anything needs to be closed,
        // i.e the interval's depth is <= to that of an open interval.
        // Add this interval's children in reverse order to intervals_to_process
        // Open the current interval's snarl/chain

        // Get the interval
        interval_state_t current_interval = std::move(forest_state.intervals_to_process.front());
        forest_state.intervals_to_process.pop_front();

        /********************

         * First, check if anything needs to be closed and close it

         ************************/

#ifdef DEBUG_ZIP_CODE_TREE
        cerr << "Process interval of type " << current_interval.code_type << " with range " 
             << current_interval.interval_start << "-" << current_interval.interval_end << endl;
        assert(current_interval.depth <= 
                seeds.at(forest_state.seed_sort_order[current_interval.interval_start]).zipcode.max_depth()+1);
        cerr << "Close anything open" << endl;
#endif
        while (!forest_state.open_intervals.empty()) {
            if (current_interval.depth <= forest_state.open_intervals.back().depth) {
                // Close open interval if the current interval is not a child
                // i.e. we have finished its children in the snarl tree
#ifdef DEBUG_ZIP_CODE_TREE
                cerr << "\tclose something at depth " << forest_state.open_intervals.size()-1 << endl;
#endif
                size_t depth = forest_state.open_intervals.size()-1;

                // The ancestor interval to close and its last seed
                const interval_state_t& ancestor_interval = forest_state.open_intervals.back();
                const Seed& last_seed = seeds.at(forest_state.seed_sort_order[ancestor_interval.interval_end-1]);

                if (ancestor_interval.code_type == ZipCode::CHAIN ||
                    ancestor_interval.code_type == ZipCode::NODE ||
                    ancestor_interval.code_type == ZipCode::ROOT_CHAIN ||
                    ancestor_interval.code_type == ZipCode::ROOT_NODE) {
                    close_chain(forest_state, depth, 
                                last_seed, ancestor_interval.is_reversed); 
                } else {
#ifdef DEBUG_ZIP_CODE_TREE
                    assert(ancestor_interval.code_type == ZipCode::REGULAR_SNARL ||
                        ancestor_interval.code_type == ZipCode::IRREGULAR_SNARL ||
                        ancestor_interval.code_type == ZipCode::CYCLIC_SNARL ||
                        ancestor_interval.code_type == ZipCode::ROOT_SNARL);
#endif
                    close_snarl(forest_state, depth, last_seed, ancestor_interval.is_reversed); 
                }

                // Take out this ancestor & forgot its list of children
                forest_state.sibling_indices_at_depth[depth].clear();
                forest_state.open_intervals.pop_back();
            } else {
                // If the current interval is contained in this open interval,
                // then it is also contained in all other ancestors so break
                break;
            }
        }

        /************ 
         * Now start processing the current interval
         *
         *
         * Sort this interval and add the child intervals 
         * in reverse order to intervals_to_process  
         ***********/
        
        if (current_interval.code_type != ZipCode::NODE ) {
            // Sort the current interval and get the intervals for its children
            sort_one_interval(forest_state, current_interval);
            get_next_intervals(forest_state, current_interval, forest_state.intervals_to_process);
        }
    
        /**********
         *
         * Open the current interval
         * If the current interval is a snarl and a child of a chain,
         * then add the preceding sibling seeds before the snarl
         *
         *******/

#ifdef DEBUG_ZIP_CODE_TREE
        cerr << "Open next interval or (if the interval is for nodes), add seeds" << endl;
#endif
        if (forest_state.open_intervals.size()+1 > forest_state.sibling_indices_at_depth.size()) {
            forest_state.sibling_indices_at_depth.emplace_back();
        }
        if (forest_state.open_intervals.empty()) {
            // If nothing open, this starts a new connected component; open it
#ifdef DEBUG_ZIP_CODE_TREE
            cerr << "Start a new connected component" << endl;
            assert(current_interval.code_type == ZipCode::ROOT_NODE ||
                    current_interval.code_type == ZipCode::NODE ||
                    current_interval.code_type == ZipCode::ROOT_CHAIN ||
                    current_interval.code_type == ZipCode::ROOT_SNARL);
#endif
            if (forest_state.active_tree_index == std::numeric_limits<size_t>::max() 
                || trees[forest_state.active_tree_index].zip_code_tree.size() != 0) {
#ifdef DEBUG_ZIP_CODE_TREE
                // If we're starting a new tree then the last one must be valid
                if (forest_state.active_tree_index != std::numeric_limits<size_t>::max()) {
                    cerr << "Last connected component: " << endl;
                    trees[forest_state.active_tree_index].print_self(forest_state.seeds);
                    trees[forest_state.active_tree_index].validate_zip_tree(
                        *forest_state.distance_index, forest_state.seeds, forest_state.distance_limit);
                }
#endif
                trees.emplace_back();
                forest_state.active_tree_index = trees.size()-1;
            }

            if (current_interval.code_type == ZipCode::ROOT_SNARL) {
                // Open the root snarl (treated as non-cyclic)
                open_snarl(forest_state, 0, false);
            } else if (current_interval.code_type == ZipCode::NODE) {
                // For a root node, just add it as a chain with all the seeds
                // Root chains have no parent snarl, so their ID is inf
                trees[forest_state.active_tree_index].zip_code_tree.emplace_back(
                    ZipCodeTree::CHAIN_START, std::numeric_limits<size_t>::max());

                // Remember the start of the chain
                forest_state.sibling_indices_at_depth[0].emplace_back(ZipCodeTree::CHAIN_START, 0);
                forest_state.sibling_indices_at_depth[0].back().chain_component = 0;

                // If this is a node, the interval contains everything in it
                // so add the seeds and close the chain here
                for (size_t i = current_interval.interval_start; i < current_interval.interval_end; i++) {
                    add_child_to_chain(forest_state, current_interval.depth, 
                                       forest_state.seed_sort_order[i], current_interval.is_reversed,
                                       current_interval.is_reversed); 
                }
                close_chain(forest_state, current_interval.depth, 
                            seeds.at(forest_state.seed_sort_order[current_interval.interval_end-1]), 
                            current_interval.is_reversed); 
            } else {
                // Open the root chain/node
                trees[forest_state.active_tree_index].zip_code_tree.emplace_back(ZipCodeTree::CHAIN_START, 
                                                                                 std::numeric_limits<size_t>::max());

                // Remember the start of the chain
                forest_state.sibling_indices_at_depth[0].emplace_back(ZipCodeTree::CHAIN_START, 0);
                forest_state.sibling_indices_at_depth[0].back().chain_component = 0;
            }                       
        } else if (forest_state.open_intervals.back().code_type == ZipCode::CHAIN || 
                   forest_state.open_intervals.back().code_type == ZipCode::ROOT_CHAIN ||
                   forest_state.open_intervals.back().code_type == ZipCode::ROOT_NODE) {
            // This is the child of a chain
            if (current_interval.code_type == ZipCode::NODE) {
                // This is a range of seeds that are on nodes on the chain,
                // not necessarily on the same node. Add each seed
                bool is_trivial_chain = current_interval.depth-1 == 
                        seeds.at(forest_state.seed_sort_order[current_interval.interval_start]).zipcode.max_depth();
                for (size_t i = current_interval.interval_start; i < current_interval.interval_end; i++) {
                    add_child_to_chain(forest_state, 
                                       is_trivial_chain ? current_interval.depth-1 : current_interval.depth, 
                                       forest_state.seed_sort_order[i], current_interval.is_reversed,
                                       forest_state.open_intervals.back().is_reversed); 
                }

            } else {
#ifdef DEBUG_ZIP_CODE_TREE
                assert(current_interval.code_type == ZipCode::REGULAR_SNARL
                       || current_interval.code_type == ZipCode::IRREGULAR_SNARL
                       || current_interval.code_type == ZipCode::CYCLIC_SNARL);
#endif
                // Add the snarl to the chain
                add_child_to_chain(forest_state, current_interval.depth,
                                   forest_state.seed_sort_order[current_interval.interval_start], 
                                   current_interval.is_reversed, forest_state.open_intervals.back().is_reversed);
            }
        } else {
            // If there's an open ancestor that isn't a chain, it must be a snarl
#ifdef DEBUG_ZIP_CODE_TREE
            assert(forest_state.open_intervals.back().code_type == ZipCode::REGULAR_SNARL ||
                    forest_state.open_intervals.back().code_type == ZipCode::IRREGULAR_SNARL ||
                    forest_state.open_intervals.back().code_type == ZipCode::CYCLIC_SNARL ||
                    forest_state.open_intervals.back().code_type == ZipCode::ROOT_SNARL);
#endif
            // Open the child chain
            open_chain(forest_state, current_interval);
            
        }

        if (current_interval.code_type != ZipCode::NODE) {
            // Add non-seed thing (i.e. snarl tree structure) to open_intervals
            forest_state.open_intervals.emplace_back(std::move(current_interval));
        }
    }
    // Finished adding all intervals: close anything that remained open
    while (!forest_state.open_intervals.empty()) {
        interval_state_t& ancestor_interval = forest_state.open_intervals.back();
        const Seed& last_seed = seeds.at(forest_state.seed_sort_order[ancestor_interval.interval_end-1]);

        if (ancestor_interval.code_type == ZipCode::CHAIN ||
            ancestor_interval.code_type == ZipCode::ROOT_CHAIN ||
            ancestor_interval.code_type == ZipCode::ROOT_NODE) {
            close_chain(forest_state, forest_state.open_intervals.size()-1, 
                        last_seed, ancestor_interval.is_reversed); 
        } else {
#ifdef DEBUG_ZIP_CODE_TREE
            assert(ancestor_interval.code_type == ZipCode::REGULAR_SNARL ||
                    ancestor_interval.code_type == ZipCode::IRREGULAR_SNARL ||
                    ancestor_interval.code_type == ZipCode::CYCLIC_SNARL ||
                    ancestor_interval.code_type == ZipCode::ROOT_SNARL);
#endif
            close_snarl(forest_state, forest_state.open_intervals.size()-1, 
                        last_seed, ancestor_interval.is_reversed); 
        }

        forest_state.open_intervals.pop_back();
    }

    // Get rid of empty tree
    if (trees[forest_state.active_tree_index].zip_code_tree.size() == 0) {
        trees.erase(trees.begin() + forest_state.active_tree_index);
    }
#ifdef DEBUG_ZIP_CODE_TREE
    validate_zip_forest(distance_index, &seeds, distance_limit);
    assert(forest_state.open_chains.empty());
    assert(forest_state.open_intervals.empty());
#endif
}

}

namespace std {

std::string to_string(const vg::ZipCodeTree::tree_item_type_t& type) {
    switch (type) {
    case vg::ZipCodeTree::SEED:
        return "SEED";
    case vg::ZipCodeTree::DAG_SNARL_START:
        return "DAG_SNARL_START";
    case vg::ZipCodeTree::DAG_SNARL_END:
        return "DAG_SNARL_END";
    case vg::ZipCodeTree::CYCLIC_SNARL_START:
        return "CYCLIC_SNARL_START";
    case vg::ZipCodeTree::CYCLIC_SNARL_END:
        return "CYCLIC_SNARL_END";
    case vg::ZipCodeTree::CHAIN_START:
        return "CHAIN_START";
    case vg::ZipCodeTree::CHAIN_END:
        return "CHAIN_END";
    case vg::ZipCodeTree::EDGE:
        return "EDGE";
    case vg::ZipCodeTree::CHAIN_COUNT:
        return "CHAIN_COUNT";
    default:
        throw std::runtime_error("Unimplemented zip code tree item type");
    }
}

std::string to_string(const vg::ZipCodeTree::distance_iterator::State& state) {
    switch (state) {
    case vg::ZipCodeTree::distance_iterator::S_START:
        return "S_START";
    case vg::ZipCodeTree::distance_iterator::S_SCAN_CHAIN:
        return "S_SCAN_CHAIN";
    case vg::ZipCodeTree::distance_iterator::S_SCAN_DAG_SNARL:
        return "S_SCAN_DAG_SNARL";
    case vg::ZipCodeTree::distance_iterator::S_SCAN_CYCLIC_SNARL:
        return "S_SCAN_CYCLIC_SNARL";
    default:
        throw std::runtime_error("Unimplemented zip code tree reverse iterator state");
    }
}

}