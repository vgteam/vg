//#define DEBUG_ZIP_CODE_TREE
//#define PRINT_NON_DAG_SNARLS
//#define DEBUG_ZIP_CODE_SORTING

#include "zip_code_tree.hpp"
#include <structures/union_find.hpp>
#include "crash.hpp"
#include "minimizer_mapper.hpp"

// Set for verbose logging from the zip code tree parsing logic
#define debug_parse

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
        } else if (item.get_type() == NODE_COUNT) {
            cerr << item.get_value() << " ";
        } else {
            throw std::runtime_error("[zip tree]: Trying to print a zip tree item of the wrong type");
        }
    }
    cerr << endl;
}

void ZipCodeForest::open_chain(forest_growing_state_t& forest_state, 
                               const size_t& depth, size_t seed_index, bool chain_is_reversed) {
    //If this is the start of a new chain
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "\t\tOpen new chain at depth " << depth << endl;
#endif
    const Seed& current_seed = forest_state.seeds->at(seed_index);

    bool is_node = current_seed.zipcode.max_depth() == depth;

    if (depth == 0) {
        //If this is the start of a new top-level chain, make a new tree,
        //which will be the new active tree
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "Add a new tree" << endl;
#endif
        if (forest_state.active_tree_index == std::numeric_limits<size_t>::max() 
            || trees[forest_state.active_tree_index].zip_code_tree.size() != 0) {
            //Don't add a new tree if the current one is empty
#ifdef DEBUG_ZIP_CODE_TREE
    //If we're starting a new tree then the last one must be valid
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
    }

    //Now record the start of this chain
    //Look up the ID of the parent snarl, unless this is a top-level chain
    size_t snarl_id = forest_state.sibling_indices_at_depth[depth-1][0].snarl_id;
    trees[forest_state.active_tree_index].zip_code_tree.emplace_back(ZipCodeTree::CHAIN_START, snarl_id);

    //Remember the chain start and its prefix sum value as a child of the chain
    forest_state.sibling_indices_at_depth[depth].push_back({ZipCodeTree::CHAIN_START, 
        chain_is_reversed ? current_seed.zipcode.get_length(depth, true) 
                          : 0});
    forest_state.sibling_indices_at_depth[depth].back().chain_component = 
        (chain_is_reversed && !is_node) ? current_seed.zipcode.get_last_chain_component(depth, true) 
                                        : 0;

    //And, if it is the child of a snarl,
    //then remember the chain as a child of the snarl
    if (depth != 0) {
        forest_state.sibling_indices_at_depth[depth-1].push_back({ZipCodeTree::CHAIN_START,
                                                     trees[forest_state.active_tree_index].zip_code_tree.size()-1});

        //The distances in the snarl include the distances
        //from the first/last children in the chain to the ends of the chains
        //
        //Remember the distance to the start of this child in the chain
        if (chain_is_reversed) {
            //If the chain is reversed, then we need to find
            //the distance to the end of the chain from the seed's prefix sum
            //and the length of the chain
            //If the length of the chain is infinite, then this is
            //not the last component of the chain and the distance is infinite
            //Otherwise, find the length of the chain/last component and
            //the length of the child, if it is a snarl
            size_t chain_length = current_seed.zipcode.get_length(depth, true);
            if (chain_length == std::numeric_limits<size_t>::max()) {
                forest_state.sibling_indices_at_depth[depth-1].back().distances.first 
                    = std::numeric_limits<size_t>::max();
            } else {
                size_t dist_val = forest_state.sort_values_by_seed[seed_index].get_distance_value();
                bool prev_is_node = is_node || current_seed.zipcode.get_code_type(depth+1) == ZipCode::NODE;
                forest_state.sibling_indices_at_depth[depth-1].back().distances.first 
                    = SnarlDistanceIndex::minus(chain_length, 
                        SnarlDistanceIndex::sum(dist_val, prev_is_node ? 0
                                                                       : current_seed.zipcode.get_length(depth+1)));
            }
        } else {
            //If the chain is traversed forward,
            //then the value is the prefix sum of the first component
            if (!is_node && current_seed.zipcode.get_chain_component(depth+1) != 0) {
                //If this isn't the first component, then it is infinite
                forest_state.sibling_indices_at_depth[depth-1].back().distances.first 
                    = std::numeric_limits<size_t>::max();
            } else {
                //Otherwise, just the prefix sum
                forest_state.sibling_indices_at_depth[depth-1].back().distances.first 
                    = forest_state.sort_values_by_seed[seed_index].get_distance_value();
                
            }
        }

        //Remember the opening of this chain, 
        //and if its first child was far enough from the start to 
        //start a new subtree
        forest_state.open_chains.emplace_back(trees[forest_state.active_tree_index].zip_code_tree.size()-1, 
                 forest_state.sibling_indices_at_depth[depth-1].back().distances.first > forest_state.distance_limit);
    }
}

unordered_map<size_t, size_t> ZipCodeTree::forget_snarls_past_index(size_t index) {
    unordered_map<size_t, size_t> removed_snarls;
    for (size_t i = 0; i < next_snarl_id; i++) {
        if (snarl_start_indexes.count(i) && snarl_start_indexes[i] >= index) {
            removed_snarls[i] = snarl_start_indexes[i] - index;
            snarl_start_indexes.erase(i);
        }
    }
    return removed_snarls;
}

bool ZipCodeForest::move_slice(forest_growing_state_t& forest_state, const size_t& depth) {
    // Chain IDs are reset since now there is no parent snarl
    tree_item_t new_start = tree_item_t(ZipCodeTree::CHAIN_START, std::numeric_limits<size_t>::max(), false);
    tree_item_t new_end = tree_item_t(ZipCodeTree::CHAIN_END, std::numeric_limits<size_t>::max(), false);

    // Might need to add an artificial end to the slice
    if (trees[forest_state.active_tree_index].zip_code_tree.back().get_type() != ZipCodeTree::CHAIN_END) {
        trees[forest_state.active_tree_index].zip_code_tree.emplace_back(new_end);
    }

    size_t slice_start_i = forest_state.open_chains.back().first;
    auto start_of_slice = trees[forest_state.active_tree_index].zip_code_tree.begin() + slice_start_i;
    auto end_of_slice = trees[forest_state.active_tree_index].zip_code_tree.end();
    bool move_full_chain = start_of_slice->get_type() == ZipCodeTree::CHAIN_START;

    // Pull out memory map of snarl start indexes to transfer to the new tree
    auto transferred_snarl_starts = trees[forest_state.active_tree_index].forget_snarls_past_index(slice_start_i);
    trees.emplace_back();

    if (!move_full_chain) {
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "Copy a slice from the middle of the chain to the end" << endl;
    assert(start_of_slice->get_type() == ZipCodeTree::SEED || start_of_slice->is_snarl_start());
#endif
        //We're copying a slice of the chain from the middle to the end
        //Start a new chain in the new subtree
        trees.back().zip_code_tree.emplace_back(new_start);

        //Scoot all the snarl start indexes forward one
        for (auto& snarl_start : transferred_snarl_starts) {
            snarl_start.second++;
        }
    }

    //Copy everything in the slice into the new tree
    trees.back().zip_code_tree.insert(trees.back().zip_code_tree.end(),
            std::make_move_iterator(start_of_slice), std::make_move_iterator(end_of_slice));
    //Erase the slice
    trees[forest_state.active_tree_index].zip_code_tree.erase(start_of_slice, end_of_slice);

    //Reset chain ends to have no ID
    trees.back().zip_code_tree.front() = new_start;
    trees.back().zip_code_tree.back() = new_end;
    // Transfer index memory
    trees.back().snarl_start_indexes = transferred_snarl_starts;

#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "Validate slice" << endl;
    trees.back().print_self(forest_state.seeds);
    trees.back().validate_zip_tree(*forest_state.distance_index, forest_state.seeds, forest_state.distance_limit);
#endif
    return move_full_chain;
}

void ZipCodeForest::close_chain(forest_growing_state_t& forest_state, 
                                const size_t& depth, const Seed& last_seed, bool chain_is_reversed) {
    bool parent_is_cyclic_snarl = depth > 0 && forest_state.sibling_indices_at_depth[depth-1].size() > 0
        && forest_state.sibling_indices_at_depth[depth-1][0].type == ZipCodeTree::CYCLIC_SNARL_START;
    // Look up the ID of the parent snarl, unless this is a top-level chain
    auto snarl_id = depth > 0 ? forest_state.sibling_indices_at_depth[depth-1][0].snarl_id 
                              : std::numeric_limits<size_t>::max();
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "\t\tclose a chain at depth " << depth  << " in ";
    cerr << (parent_is_cyclic_snarl ? "cyclic snarl" : "DAG snarl") << endl;
#endif
    if (trees[forest_state.active_tree_index].zip_code_tree.back().get_type() == ZipCodeTree::CHAIN_START) {
        //If the chain was empty.
        //This could happen if there was only a snarl in it and it got removed

        //Take out the CHAIN_START
        trees[forest_state.active_tree_index].zip_code_tree.pop_back();

        //Forget about this chain in its parent snarl
        if (trees[forest_state.active_tree_index].zip_code_tree.size() > 0
            && (trees[forest_state.active_tree_index].zip_code_tree.back().is_snarl_start()
                || trees[forest_state.active_tree_index].zip_code_tree.back().get_type() == ZipCodeTree::CHAIN_END)) {
            forest_state.sibling_indices_at_depth[depth-1].pop_back();
        }

        //Forget about the chain
        if (depth != 0) {
            forest_state.open_chains.pop_back();
        }

    } else {
        //Otherwise, the chain wasn't empty so actually close it


        //Add the end of the chain to the zip code tree
        trees[forest_state.active_tree_index].zip_code_tree.emplace_back(ZipCodeTree::CHAIN_END, snarl_id, false);

        if (depth == 0) {
            return;
        }

        // For chains in snarls, we want to know the distance
        // from the last thing in the chain to the end of the chain
        // If the distance is greater than the distance limit,
        // we may make a new tree for a slice of the chain.
        // If the chain remains in the snarl, 
        // we need to remember the distance to the end
        // of the chain to add to the relevant distances in the parent snarl.
        // These distances will be stored in sibling_indices_at_depth

#ifdef DEBUG_ZIP_CODE_TREE
    assert(forest_state.sibling_indices_at_depth[depth-1].size() > 0);
    assert(forest_state.sibling_indices_at_depth[depth-1].back().type == ZipCodeTree::CHAIN_START);
#endif
        //Only add the distance for a non-root chain

        //If this is reversed, then the distance should be
        //the distance to the start of the chain. 
        //Otherwise, the distance to the end
        //The value that got stored in forest_state.sibling_indices_at_depth was
        //the prefix sum traversing the chain according to its tree orientation,
        //so either way the distance is the length of the chain - the prefix sum
        size_t distance_to_chain_end = chain_is_reversed 
                                     ? forest_state.sibling_indices_at_depth[depth].back().value
                                     : SnarlDistanceIndex::minus(last_seed.zipcode.get_length(depth),
                                          forest_state.sibling_indices_at_depth[depth].back().value);
        bool chain_kept = true;
        if (distance_to_chain_end > forest_state.distance_limit && forest_state.open_chains.back().second) {
            //If the distance to the end is greater than the distance limit,
            //and there was something in the chain with
            //a large distance to the thing before it, then splice out a slice
            bool moved_full_chain = move_slice(forest_state, depth);

            if (moved_full_chain) {
                chain_kept = false;
                //The chain no longer exists in the snarl, so forget it
                forest_state.sibling_indices_at_depth[depth-1].pop_back();
             } else {
                //Take out the last edge
                size_t last_edge = trees[forest_state.active_tree_index].zip_code_tree.back().get_value();
                trees[forest_state.active_tree_index].zip_code_tree.pop_back();

                //Close the chain in the original active tree
                trees[forest_state.active_tree_index].zip_code_tree.emplace_back(ZipCodeTree::CHAIN_END, snarl_id);
                //Update the distance to the end of the chain 
                //to be the distance from the previous child 
                size_t last_length = depth == last_seed.zipcode.max_depth() 
                                   ? 0 
                                   : last_seed.zipcode.get_length(depth+1);

                distance_to_chain_end = SnarlDistanceIndex::sum(distance_to_chain_end, 
                                        SnarlDistanceIndex::sum(last_edge,
                                                                last_length));
            }
        }
        if (chain_kept) {
            // Remember chain's orientation and distance from last seed to edge
            forest_state.sibling_indices_at_depth[depth-1].back().is_reversed = chain_is_reversed;
            forest_state.sibling_indices_at_depth[depth-1].back().distances.second = distance_to_chain_end;
        }
        //We've closed a chain, so take out the latest open chain
        forest_state.open_chains.pop_back();
    }
}

void ZipCodeForest::add_child_to_chain(forest_growing_state_t& forest_state, 
                       const size_t& depth, const size_t& seed_index, bool child_is_reversed,
                       bool chain_is_reversed) {
    const Seed& current_seed = forest_state.seeds->at(seed_index);
    auto cur_id = depth > 1 ? forest_state.sibling_indices_at_depth[depth-2][0].snarl_id 
                            : std::numeric_limits<size_t>::max();

    ZipCode::code_type_t current_type = current_seed.zipcode.get_code_type(depth);

    //Is this chain actually a node pretending to be a chain
    bool is_trivial_chain = current_type == ZipCode::CHAIN && depth == current_seed.zipcode.max_depth();

    //For a root node or trivial chain, the "chain" is actually just the node,
    //so the depth  of the chain we're working on is the same depth.
    //Otherwise, the depth is depth-1
    size_t chain_depth = is_trivial_chain || current_type == ZipCode::ROOT_NODE ? depth : depth-1;

    ///////////////// Get the offset in the parent chain (or node)
    size_t current_offset;


    //First, get the prefix sum in the chain + offset in the node
    if (current_type == ZipCode::ROOT_NODE || current_type == ZipCode::NODE || is_trivial_chain) {
        //For a node, this is still the distance used to sort on
        current_offset = forest_state.sort_values_by_seed[seed_index].get_distance_value();
    } else {
        //Otherwise, get the distance to the start or end of the chain

        current_offset = current_seed.zipcode.get_offset_in_chain(depth);
    }
    if (chain_is_reversed && !(current_type == ZipCode::NODE || current_type == ZipCode::ROOT_NODE 
                               || is_trivial_chain)) {
        //If we are adding a snarl and the chain is being traversed backwards,
        //then make sure the prefix sum is going to the right end of the snarl
        current_offset =  SnarlDistanceIndex::sum(current_offset, current_seed.zipcode.get_length(depth));
    }


    /////////////////////// Get previous thing's offset in parent chain/node
    size_t previous_offset = forest_state.sibling_indices_at_depth[chain_depth][0].value;


#ifdef DEBUG_ZIP_CODE_TREE
    assert(forest_state.sibling_indices_at_depth[chain_depth].size() == 1);
#endif

    ///////////////////// Record distance from the previous thing in chain/node
    //       Or add a new tree if the distance is too far
    auto first_type = forest_state.sibling_indices_at_depth[chain_depth][0].type;
    // Assume we want to use this child
    bool skip_child = false;
    if (chain_depth > 0 && first_type == ZipCodeTree::CHAIN_START) {

        //If this is the first thing in a non-root chain or node,
        //remember the distance to the start of the chain/node.
        //This distance will be added to distances in the parent snarl
        forest_state.sibling_indices_at_depth[chain_depth-1][0].distances.first = chain_is_reversed 
            ? SnarlDistanceIndex::minus(current_seed.zipcode.get_length(chain_depth, true),
                SnarlDistanceIndex::sum(current_offset,
                    (is_trivial_chain || current_type == ZipCode::NODE 
                        ? 0 : current_seed.zipcode.get_length(chain_depth+1))))
            : current_offset;

        //Update the last chain opened
        forest_state.open_chains.back().second = std::max(current_offset, previous_offset) 
                                                 - std::min(current_offset, previous_offset) 
                                                 > forest_state.distance_limit;


    } else if (first_type != ZipCodeTree::CHAIN_START) {
        //for all except the first thing in a node/chain, we need to add an edge

        size_t distance_between;
        if (!is_trivial_chain && !current_type == ZipCode::ROOT_NODE 
            && forest_state.sibling_indices_at_depth[chain_depth][0].chain_component 
               != current_seed.zipcode.get_chain_component(depth)) {
            //If the parent is a multicomponent chain,
            //then they might be in different components
            distance_between = std::numeric_limits<size_t>::max();
        } else {
            distance_between = std::max(current_offset, previous_offset) - std::min(current_offset, previous_offset);
        }

        if (chain_depth == 0 && distance_between > forest_state.distance_limit) {
            //The next thing in the zip tree will be the first seed (or snarl)
            // in a top-level chain, so start a new tree
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "Start a new tree in the forest" << endl;
#endif
            //Close the previous chain
            trees[forest_state.active_tree_index].zip_code_tree.emplace_back(ZipCodeTree::CHAIN_END, cur_id);

            if (forest_state.active_tree_index == std::numeric_limits<size_t>::max() 
                || trees[forest_state.active_tree_index].zip_code_tree.size() != 0) {
                //Add a new tree and make sure it is the new active tree
#ifdef DEBUG_ZIP_CODE_TREE
    //If we're starting a new tree then the last one must be valid
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

            //Add the start of the new chain
            trees[forest_state.active_tree_index].zip_code_tree.emplace_back(
                ZipCodeTree::CHAIN_START, std::numeric_limits<size_t>::max());

            //The first sibling in the chain is now the chain start,
            //not the previous seed, so replace it
            forest_state.sibling_indices_at_depth[chain_depth].pop_back();
            forest_state.sibling_indices_at_depth[chain_depth].push_back({ZipCodeTree::CHAIN_START, 
                chain_is_reversed ? current_seed.zipcode.get_length(chain_depth, true) : 0}); 
            forest_state.sibling_indices_at_depth[chain_depth].back().chain_component = 
                !is_trivial_chain ? current_seed.zipcode.get_last_chain_component(chain_depth, true) : 0;

        } else if (distance_between > forest_state.distance_limit) { 
            //If this is too far from the previous thing, but inside a snarl

            if (forest_state.open_chains.back().second) {
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "\tMake a new slice of the chain at depth " << depth << endl;
#endif
                bool moved_full_chain = move_slice(forest_state, chain_depth);
                //If the current chain slice was also
                //too far away from the thing before it then copy the slice
                if (moved_full_chain) {
                    //If the slice starts at the start of the chain
                    //and ends at the previous seed;

                    //Add back the start of the chain
                    trees[forest_state.active_tree_index].zip_code_tree.emplace_back(ZipCodeTree::CHAIN_START, cur_id);

                    //Update the chain as a child of the snarl
#ifdef DEBUG_ZIP_CODE_TREE
    assert(forest_state.sibling_indices_at_depth[chain_depth-1].back().type == ZipCodeTree::CHAIN_START);   
    //The value should be the index of the last seed,
    //which is the first seed in the new tree
    assert(forest_state.sibling_indices_at_depth[chain_depth-1].back().value 
                == trees[forest_state.active_tree_index].zip_code_tree.size()-1);
    assert(forest_state.open_chains.back().second);
#endif
                    forest_state.sibling_indices_at_depth[chain_depth-1].back().distances.first = chain_is_reversed 
                        ? SnarlDistanceIndex::minus(current_seed.zipcode.get_length(chain_depth, true),
                            SnarlDistanceIndex::sum(current_offset,
                                (is_trivial_chain || current_type == ZipCode::NODE 
                                    ? 0 : current_seed.zipcode.get_length(chain_depth+1))))
                        : current_offset;

                    //Don't need to update open_chains, since the next slice 
                    //will also start at the chain start and be able to make 
                    //a new thing
                } else {
                    //If the slice starts and ends in the middle of the chain
                    //The original tree gets an edge with infinite length, since
                    //it will be bigger than the distance limit anyway
#ifdef DEBUG_ZIP_CODE_TREE
    assert(trees[forest_state.active_tree_index].zip_code_tree.back().get_type() == ZipCodeTree::EDGE);
#endif
                    trees[forest_state.active_tree_index].zip_code_tree.pop_back();
                    trees[forest_state.active_tree_index].zip_code_tree.emplace_back(
                        ZipCodeTree::EDGE, std::numeric_limits<size_t>::max());

                    //Remember the next seed or snarl that gets added
                    //as the start of a new chain slice
                    forest_state.open_chains.pop_back();
                    forest_state.open_chains.emplace_back(
                        trees[forest_state.active_tree_index].zip_code_tree.size(), true);
                }
            } else {
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "The slice didn't get copied but maybe start a new slice here" << endl;
#endif
                //If the slice doesn't get copied
                //because it is still connected at the front,
                //add the edge to the chain
                //and remember that it could start a new slice

                trees[forest_state.active_tree_index].zip_code_tree.emplace_back(ZipCodeTree::EDGE, distance_between);

                //Remember the next seed or snarl that gets added 
                //as the start of a new chain slice
                forest_state.open_chains.pop_back();
                forest_state.open_chains.emplace_back(
                    trees[forest_state.active_tree_index].zip_code_tree.size(), true);
            }

        } else {
            //If we didn't start a new tree, then add the edge
            trees[forest_state.active_tree_index].zip_code_tree.emplace_back(ZipCodeTree::EDGE, distance_between);
        }
    }

    /////////////////////////////Record this thing in the chain
    if (current_type == ZipCode::NODE || current_type == ZipCode::ROOT_NODE || is_trivial_chain) {
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "\t\tContinue node/chain with seed " << current_seed.pos << " at depth " << depth << endl;
#endif
        //If this was a node, just remember the seed
        if (!skip_child) {
            trees[forest_state.active_tree_index].zip_code_tree.emplace_back(
                ZipCodeTree::SEED, seed_index, child_is_reversed != is_rev(current_seed.pos));
        }
    } else {

        open_snarl(forest_state, depth, current_seed.zipcode.get_code_type(depth) == ZipCode::CYCLIC_SNARL); 

        //For finding the distance to the next thing in the chain, the offset
        //stored should be the offset of the end bound of the snarl, so add the 
        //length of the snarl
        current_offset = chain_is_reversed 
                       ? SnarlDistanceIndex::minus(current_offset, current_seed.zipcode.get_length(depth))
                       : SnarlDistanceIndex::sum(current_offset, current_seed.zipcode.get_length(depth));

    }

    if (!skip_child) {
        //Remember this thing for the next sibling in the chain
        forest_state.sibling_indices_at_depth[chain_depth].pop_back();
        forest_state.sibling_indices_at_depth[chain_depth].push_back({
            (current_type == ZipCode::NODE || current_type == ZipCode::ROOT_NODE) ? ZipCodeTree::SEED 
                                                                                  : ZipCodeTree::DAG_SNARL_START,
            current_offset}); 
        
        if (forest_state.sibling_indices_at_depth[chain_depth].back().type == ZipCodeTree::DAG_SNARL_START) {
            // If copying a snarl, remember its ID
            forest_state.sibling_indices_at_depth[chain_depth].back().snarl_id = 
                trees[forest_state.active_tree_index].zip_code_tree.back().get_value();
        }

        if (!is_trivial_chain && !current_type == ZipCode::ROOT_NODE) {
            // If needed, remember the chain component
            forest_state.sibling_indices_at_depth[chain_depth].back().chain_component 
                = current_seed.zipcode.get_chain_component(depth);
        }
    }
#ifdef DEBUG_ZIP_CODE_TREE
    if (!skip_child) {
        cerr << "Added sibling with type " << current_type << endl;
    } else {
        cerr << "Child not added" << endl;
    }
#endif

}

void ZipCodeForest::open_snarl(forest_growing_state_t& forest_state, const size_t& depth, bool is_cyclic_snarl) {
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "\t\tOpen new snarl at depth " << depth << endl;
#endif
    auto opening_type = is_cyclic_snarl ? ZipCodeTree::CYCLIC_SNARL_START : ZipCodeTree::DAG_SNARL_START;
    //If this was a snarl, record the start of the snarl
    size_t cur_id = trees[forest_state.active_tree_index].open_snarl(is_cyclic_snarl);

    
    if (depth != 0) {
        //Remember the start of the snarl to find distances later
        //Don't do this for a root snarl because technically there is
        //no start node so there are no distances to it
        forest_state.sibling_indices_at_depth[depth].push_back({opening_type, std::numeric_limits<size_t>::max()});
        forest_state.sibling_indices_at_depth[depth].back().snarl_id = cur_id;
    }
}

size_t ZipCodeTree::open_snarl(bool is_cyclic_snarl) {
    size_t cur_id = next_snarl_id++;
    zip_code_tree.emplace_back(is_cyclic_snarl ? CYCLIC_SNARL_START : DAG_SNARL_START, cur_id);
    // Remember the start was put
    snarl_start_indexes[cur_id] = zip_code_tree.size() - 1;
    return cur_id;
}

void ZipCodeForest::close_snarl(forest_growing_state_t& forest_state, 
                                const size_t& depth, const Seed& last_seed, bool last_is_reversed) {
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "\t\tclose a snarl at depth " << depth << endl;
#endif
    bool is_cyclic_snarl = forest_state.sibling_indices_at_depth[depth][0].type == ZipCodeTree::CYCLIC_SNARL_START;
    size_t cur_id = forest_state.sibling_indices_at_depth[depth][0].snarl_id;

    if (trees[forest_state.active_tree_index].zip_code_tree.size() == 1) {
        //If this would be an empty snarl, then just remove it
        trees.erase(trees.begin() + forest_state.active_tree_index);
        trees[forest_state.active_tree_index].snarl_start_indexes.erase(cur_id);
    } else if (depth == 0) {
        //If this is a root snarl, then we don't need distances so just close it
        //Root snarls are treated as DAG snarls; 
        //they don't store distances so it doesn't matter
        trees[forest_state.active_tree_index].zip_code_tree.emplace_back(ZipCodeTree::DAG_SNARL_END, 
                                                                         std::numeric_limits<size_t>::max());

    } else if (forest_state.sibling_indices_at_depth[depth].size() == 1) {
        //Since some of the children of the snarl
        //may have been removed to separate subtrees,
        //the snarl may actually be empty now
        //If there is only one "child" (the snarl start), 
        //then the snarl is actually empty, so delete it

#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "\t\t\tThe snarl is actually empty so remove it" << endl;
    assert(trees[forest_state.active_tree_index].zip_code_tree.back().is_snarl_start());
#endif        
        //Pop the snarl start out
        trees[forest_state.active_tree_index].zip_code_tree.pop_back();
        trees[forest_state.active_tree_index].snarl_start_indexes.erase(cur_id);

        //If this was the first thing in the chain, then we're done.
        //Otherwise, there was an edge to remove
        if (trees[forest_state.active_tree_index].zip_code_tree.back().get_type() == ZipCodeTree::EDGE) {
            //If the snarl was in the middle of a chain,
            //then we need to take out the edge and update
            //the previous thing in the chain with its prefix sum

            //This was the distance from the last thing to the snarl start
            size_t previous_edge = trees[forest_state.active_tree_index].zip_code_tree.back().get_value();
            trees[forest_state.active_tree_index].zip_code_tree.pop_back();

            //This is the distance from the start of the chain to the snarl end
            size_t snarl_prefix_sum = forest_state.sibling_indices_at_depth[depth-1].back().value;
            forest_state.sibling_indices_at_depth[depth-1].pop_back();

            //Snarl prefix sum is now distance from chain start to snarl start
            snarl_prefix_sum = SnarlDistanceIndex::minus(snarl_prefix_sum, last_seed.zipcode.get_length(depth));

            //Now update sibling_indices_at_depth to be previous thing in chain
            forest_state.sibling_indices_at_depth[depth-1].push_back({
                trees[forest_state.active_tree_index].zip_code_tree.back().get_type(),
                SnarlDistanceIndex::minus(snarl_prefix_sum, previous_edge)});
            //If it was in the first component, then this is correct.
            //If it was in a later component, then it was too 
            //far away anyway so it doesn't matter
            //TODO: I think this might cause problems if it was a looping chain
            forest_state.sibling_indices_at_depth[depth-1].back().chain_component = 0;


            //At this point, the open_chain for the parent chain is either
            //before the removed snarl, the snarl itself, or after the snarl.
            //If the open_chain was before or at the snarl, nothing has changed.
            //If it is after the snarl, then the snarl
            //wasn't the start of a new slice so we back it up to the previous 
            //child and say that it was not the start of a new slice.
            //TODO
            //If it was the snarl itself, then the next child added to the chain
            //will be the next open_chain, but I
            //haven't implemented this yet- it won't change the correctness
            if (depth > 0 && forest_state.open_chains.size() > 0 
                && forest_state.open_chains.back().first 
                   >= trees[forest_state.active_tree_index].zip_code_tree.size()) {
                //If a chain slice could have started at or after this snarl
#ifdef DEBUG_ZIP_CODE_TREE
    assert(forest_state.open_chains.back().second);
#endif
                //Find the start of the previous child
                size_t previous_index = trees[forest_state.active_tree_index].zip_code_tree.size() - 1;
                bool found_sibling = false;
                size_t opened_snarls = 0;
                tree_item_t previous_item;
                while (!found_sibling) {
                    previous_item = trees[forest_state.active_tree_index].zip_code_tree.at(previous_index);
                    if (opened_snarls == 0 && previous_item.get_type() == ZipCodeTree::SEED) {
                        found_sibling = true;
                    } else if (previous_item.is_snarl_end()) {
                        opened_snarls ++;
                        previous_index--;
                    } else if (previous_item.is_snarl_start() && opened_snarls == 0) {
                        found_sibling = true;
                    } else if (previous_item.is_snarl_start()) {
                        opened_snarls--;
                        previous_index--;
                    } else {
                        previous_index--;
                    }
                }
                previous_item = trees[forest_state.active_tree_index].zip_code_tree.at(previous_index);
                if (previous_index != 0 && previous_item.get_type() == ZipCodeTree::CHAIN_START) {
                    previous_index--;
                }
#ifdef DEBUG_ZIP_CODE_TREE
    previous_item = trees[forest_state.active_tree_index].zip_code_tree.at(previous_index);
    assert(previous_item.get_type() == ZipCodeTree::SEED 
            || previous_item.is_snarl_start() || previous_item.get_type() == ZipCodeTree::CHAIN_START);
    cerr << "New start of previous open chain: " << previous_index << endl;;
#endif
                forest_state.open_chains.back().first = previous_index;
                forest_state.open_chains.back().second = false;
                
            }
#ifdef DEBUG_ZIP_CODE_TREE
    assert(forest_state.sibling_indices_at_depth[depth-1].back().value >= 0);
#endif
        } else {
            //If this was the first thing in the chain, 
            //update the previous sibling in the chain to be the chain start
#ifdef DEBUG_ZIP_CODE_TREE
    assert(trees[forest_state.active_tree_index].zip_code_tree.back().get_type() 
           == ZipCodeTree::CHAIN_START);
#endif
            bool prev_reversed = forest_state.open_intervals[forest_state.open_intervals.size()-2].is_reversed;
            size_t offset = prev_reversed ? last_seed.zipcode.get_length(depth-1, true) : 0;
            forest_state.sibling_indices_at_depth[depth-1].pop_back();
            forest_state.sibling_indices_at_depth[depth-1].push_back({ZipCodeTree::CHAIN_START, offset});
            forest_state.sibling_indices_at_depth[depth-1].back().chain_component 
                = last_seed.zipcode.get_last_chain_component(depth-1, true);

        }
    } else {
        // Finish off snarl by adding in a distance matrix
        add_distance_matrix(forest_state, depth, last_is_reversed);
        auto closing_type = is_cyclic_snarl ? ZipCodeTree::CYCLIC_SNARL_END : ZipCodeTree::DAG_SNARL_END;
        trees[forest_state.active_tree_index].zip_code_tree.emplace_back(closing_type, cur_id);
     }
}

void ZipCodeTree::shift_snarls_forward(size_t start_snarl_id, size_t shift_amount) {
    size_t move_after_index = snarl_start_indexes[start_snarl_id];
    for (auto& snarl_start : snarl_start_indexes) {
        if (snarl_start.second > move_after_index) {
            // Shift the snarl start index
            snarl_start.second += shift_amount;
        }
    }
}

vector<ZipCodeForest::seed_info_t> ZipCodeForest::get_edge_seeds(const forest_growing_state_t& forest_state, 
                                                                 const size_t& depth) const {
    vector<seed_info_t> edge_seeds;
    // Get edge seeds for each chain
    for (size_t i = 1; i < forest_state.sibling_indices_at_depth[depth].size(); ++i) {
        child_info_t cur_chain = forest_state.sibling_indices_at_depth[depth][i];

        // Search for first seed
        size_t seed_i = cur_chain.value + 1;
        while (trees[forest_state.active_tree_index].zip_code_tree[seed_i].get_type() != ZipCodeTree::SEED) {
            seed_i++;
        }
        edge_seeds.emplace_back(trees[forest_state.active_tree_index].zip_code_tree[seed_i].get_value(),
            cur_chain.distances.first, cur_chain.is_reversed, forest_state);

        // Search for chain end
        while (trees[forest_state.active_tree_index].zip_code_tree[seed_i].get_type() != ZipCodeTree::CHAIN_END) {
            seed_i++;
        }
        // Back up to last seed
        while (trees[forest_state.active_tree_index].zip_code_tree[seed_i].get_type() != ZipCodeTree::SEED) {
            seed_i--;
        }
        edge_seeds.emplace_back(trees[forest_state.active_tree_index].zip_code_tree[seed_i].get_value(),
            cur_chain.distances.second, !cur_chain.is_reversed, forest_state);
    }

    return edge_seeds;
}

void ZipCodeForest::add_edges_to_end(vector<tree_item_t>& dist_matrix,
                                     forest_growing_state_t& forest_state, 
                                     const size_t& depth, const vector<seed_info_t>& edge_seeds,
                                     bool snarl_is_reversed, bool is_cyclic_snarl) const {
    // start -> end is simply length of snarl
    dist_matrix.emplace_back(ZipCodeTree::EDGE, edge_seeds[0].seed.zipcode.get_length(depth));

    // DAG snarls only have distances from chain ends to snarl end
    size_t start_i = is_cyclic_snarl ? 0 : 1;
    size_t increment = is_cyclic_snarl ? 1 : 2;
    for (size_t i = start_i; i < edge_seeds.size(); i += increment) {
        // Distance from the start of the snarl to the start of the chain
        size_t between_bounds_dist = edge_seeds[i].seed.zipcode.get_distance_to_snarl_bound(
            depth+1, snarl_is_reversed, !edge_seeds[i].right_side);
        // Overall edge distance
        size_t edge_dist = SnarlDistanceIndex::sum(between_bounds_dist, edge_seeds[i].flank_offset);
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
    bool is_regular_snarl = edge_seeds[0].seed.zipcode.get_code_type(depth) == ZipCode::REGULAR_SNARL;
    net_handle_t snarl_handle;
    if (!is_regular_snarl) {
        // Snarl handles are only used for irregular and cyclic snarls
        snarl_handle = edge_seeds[0].seed.zipcode.get_net_handle(depth, forest_state.distance_index);
    }

    // Distance in between chain bounds
    size_t between_chain_dist;
    // Distance from chain bounds to the edge of the seeds
    size_t chain_flank_dist;
    // Final edge distance
    size_t edge_dist;

    // DAG snarls find distances FROM ends TO starts only
    size_t increment = is_cyclic_snarl ? 1 : 2;
    for (size_t to_i = 0; to_i < edge_seeds.size(); to_i += increment) {
        // Edge from start
        between_chain_dist = edge_seeds[to_i].seed.zipcode.get_distance_to_snarl_bound(
            depth+1, !snarl_is_reversed, !edge_seeds[to_i].right_side);
        edge_dist = SnarlDistanceIndex::sum(between_chain_dist, edge_seeds[to_i].flank_offset);
        dist_matrix.emplace_back(ZipCodeTree::EDGE, edge_dist);

        // Current chain
        size_t to_rank = edge_seeds[to_i].seed.zipcode.get_rank_in_snarl(depth+1);

        size_t start_from = is_cyclic_snarl ? 0 : 1;
        // Start at previous chain end for DAGs,
        // but make sure to not overflow size_t
        size_t end_from = is_cyclic_snarl ? to_i : (to_i == 0 ? 0 : to_i - 1);
        for (int64_t from_i = start_from; from_i <= end_from; from_i += increment) {
            size_t from_rank = edge_seeds[from_i].seed.zipcode.get_rank_in_snarl(depth+1);

            if (is_regular_snarl) {
                // Distance between chains in a regular snarl is always inf
                between_chain_dist = std::numeric_limits<size_t>::max();
            } else {
                // Irregular/cyclic snarls have a snarl_handle we can use
                between_chain_dist = forest_state.distance_index->distance_in_snarl(snarl_handle, 
                        from_rank, edge_seeds[from_i].right_side, to_rank, edge_seeds[to_i].right_side);
            }

            chain_flank_dist = SnarlDistanceIndex::sum(edge_seeds[from_i].flank_offset, edge_seeds[to_i].flank_offset);
            if (from_i == to_i) {
                // Self loop; if the between_chain_distance is non infinite
                // then there must be a self-loop at the edge of the chain
                // thus we only need to double back over the chain flanks
                bool loop_is_possible = between_chain_dist == std::numeric_limits<size_t>::max();
                edge_dist = loop_is_possible ? std::numeric_limits<size_t>::max()
                                             : chain_flank_dist;
            } else {
                // Distance between chains
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
    trees[forest_state.active_tree_index].print_self(forest_state.seeds);
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
    dist_matrix.emplace_back(ZipCodeTree::NODE_COUNT, sibling_count - 1);

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
    trees[forest_state.active_tree_index].zip_code_tree.insert(
        trees[forest_state.active_tree_index].zip_code_tree.begin() + first_child_index,
        dist_matrix.begin(), dist_matrix.end());
    // Shift memorized snarl start indexes
    trees[forest_state.active_tree_index].shift_snarls_forward(
        forest_state.sibling_indices_at_depth[depth][0].snarl_id, num_edges + 1);
}

std::pair<size_t, size_t> ZipCodeTree::dag_and_non_dag_snarl_count(const vector<Seed>& seeds, 
                                                                   const SnarlDistanceIndex& distance_index) const {
    size_t dag_count = 0;
    size_t non_dag_count = 0;

    /* Walk through everything in the zip code tree 
       and at the first seed in each snarl, 
       check if it is a dag or not
    */

    //Keep track of the depth to check the zip codes
    size_t current_depth = 0;

    //When we encounter the start of a snarl, make a note of the depth. 
    //At the next seed, check the snarls at the depths recorded
    vector<size_t> snarl_depths;

    for (size_t i = 0 ; i < zip_code_tree.size() ; i++ ) {
        const tree_item_t& current_item = zip_code_tree[i];
        if (current_item.is_snarl_start()) {
            //For the start of a snarl, make a note of depth to check the seed
            snarl_depths.emplace_back(current_depth);

            //Increment the depth
            current_depth++;
        } else if (current_item.get_type() == ZipCodeTree::CHAIN_START) {
            //For the start of a chain, increment the depth
            current_depth++;
        } else if (current_item.get_type() == ZipCodeTree::CHAIN_END || current_item.is_snarl_end()) {
            //For the end of a snarl or chain, decrement the depth
            current_depth--;
        } else if (current_item.get_type() == ZipCodeTree::SEED) {
            ZipCode cur_zip = seeds[current_item.get_value()].zipcode;
            //If this is a seed, check the snarls we've seen previously
            for (const size_t& snarl_depth : snarl_depths) {
                if (cur_zip.get_code_type(snarl_depth) == ZipCode::REGULAR_SNARL) {
                    //If this is a regular snarl, then it must be a DAG too
                    dag_count++;
                } else {
                    //If this is an irregular snarl

                    //Check the snarl in the distance index
                    net_handle_t snarl_handle = cur_zip.get_net_handle(snarl_depth, &distance_index);
#ifdef DEBUG_ZIP_CODE_TREE
    assert(cur_zip.get_code_type(snarl_depth) == ZipCode::IRREGULAR_SNARL ||
           cur_zip.get_code_type(snarl_depth) == ZipCode::CYCLIC_SNARL ||
           cur_zip.get_code_type(snarl_depth) == ZipCode::ROOT_SNARL);
    assert(distance_index.is_snarl(snarl_handle));
#endif
                    if (distance_index.is_dag(snarl_handle)) {
                        dag_count++;
                    } else {
                        non_dag_count++;
#ifdef PRINT_NON_DAG_SNARLS
    size_t child_count = 0;
    distance_index.for_each_child(snarl_handle, [&](const net_handle_t& child) {
        child_count++;
    });
    cerr << distance_index.net_handle_as_string(snarl_handle) << "\t" << child_count << endl;
#endif
                    }
                }

            }
            //Clear the snarls
            snarl_depths.clear();
        }
    }

    return std::make_pair(dag_count, non_dag_count);
}

bool ZipCodeTree::seed_is_reversed_at_depth (const Seed& seed, size_t depth, const SnarlDistanceIndex& distance_index) {
    if (seed.zipcode.get_is_reversed_in_parent(depth)) {
        return true;
    } else if (depth > 0 && (seed.zipcode.get_code_type(depth-1) == ZipCode::IRREGULAR_SNARL
                             || seed.zipcode.get_code_type(depth-1) == ZipCode::CYCLIC_SNARL)) {
        //If parent is an irregular snarl, check orientation of child in snarl
        net_handle_t snarl_handle = seed.zipcode.get_net_handle(depth-1, &distance_index);
        size_t rank = seed.zipcode.get_rank_in_snarl(depth);
        if (distance_index.distance_in_snarl(snarl_handle, 0, false, rank, false)
                    == std::numeric_limits<size_t>::max()
            &&
            distance_index.distance_in_snarl(snarl_handle, 1, false, rank, true)
                    == std::numeric_limits<size_t>::max()) {
            //If the distance from the snarl start to the child start is inf
            //and the distance from the snarl end to the child end is inf
            //then we assume that this child is "reversed" in the parent snarl
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
    while (!distance_index.is_root(net)) {
        if (distance_index.is_looping_chain(net)) {
            is_invalid = true;
            break;
        } else if (distance_index.is_chain(distance_index.get_parent(net)) && 
                    !distance_index.is_trivial_chain(distance_index.get_parent(net))) {
            //Check if this net_handle_t could be involved in
            //a chain loop that is smaller than the distance limit
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

bool ZipCodeTree::node_is_in_cyclic_snarl(nid_t id, const SnarlDistanceIndex& distance_index) const {
    bool is_cyclic_snarl = false;
    net_handle_t net = distance_index.get_node_net_handle(id);
    while (!distance_index.is_root(net)) {
        if (distance_index.is_snarl(net) && !distance_index.is_dag(net)) {
            //If this is a cyclic snarl
            is_cyclic_snarl = true;;
            break;
        }
        net = distance_index.get_parent(net);
    }
    return is_cyclic_snarl;
}

void ZipCodeTree::validate_boundaries(const SnarlDistanceIndex& distance_index, 
                                      const vector<Seed>* seeds,
                                      size_t distance_limit) const {
#ifdef DEBUG_ZIP_CODE_TREE
    std::cerr << "Validating that zip code tree's boundaries match up" << std::endl;
#endif
    bool has_seed = false;
    vector<tree_item_type_t> tree_stack; 
    vector<size_t> snarl_id_stack;
    for (size_t i = 0 ; i < zip_code_tree.size() ; i++) {
        const tree_item_t& item = zip_code_tree[i];
        if (item.is_snarl_start()) {
            if (tree_stack.size() == 1) {
                //Also check snarl distances and child count 
                //for non-root top-level snarls (is recursive)
                vector<tree_item_t>::const_iterator cur_snarl_start = zip_code_tree.begin() + i;
                validate_snarl(cur_snarl_start, distance_index, seeds, distance_limit);
            }

            if (!snarl_id_stack.empty()) {
                // Snarl IDs must be in increasing order
                assert(snarl_id_stack.back() < item.get_value());
            }
            snarl_id_stack.push_back(item.get_value());
            tree_stack.push_back(item.get_type());
        } else if (item.get_type() == CHAIN_START) {
            // Snarl ID should match, except for top-level chain
            assert(tree_stack.empty() || snarl_id_stack.back() == item.get_value());
            tree_stack.push_back(item.get_type());
        } else if (item.is_snarl_end()) {
            // Should have opened with the correct snarl type
            assert(tree_stack.back() == (item.get_type() == DAG_SNARL_END ? DAG_SNARL_START : CYCLIC_SNARL_START));
            tree_stack.pop_back();
            
            // IDs should match
            assert(snarl_id_stack.back() == item.get_value());
            snarl_id_stack.pop_back();
        } else if (item.get_type() == CHAIN_END) {
            // Snarl ID should match, except for top-level chain
            assert(tree_stack.size() == 1 || snarl_id_stack.back() == item.get_value());

            assert(tree_stack.back() == CHAIN_START);
            tree_stack.pop_back();
            assert(tree_stack.empty() || tree_stack.back() == DAG_SNARL_START 
                   || tree_stack.back() == CYCLIC_SNARL_START);
        } else if (item.get_type() == SEED) {
            has_seed = true;
        }
    }
    assert(has_seed);
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
            //Check if this is worth validating
            //Use a distance limit of 0 so it will ignore looping chains
            bool current_is_invalid = node_is_invalid(id(seeds->at(current_item.get_value()).pos), distance_index, 0);
            bool current_is_in_cyclic_snarl = node_is_in_cyclic_snarl(id(seeds->at(current_item.get_value()).pos), 
                                                                      distance_index);

            if (previous_seed_index != std::numeric_limits<size_t>::max() &&
                !current_is_invalid && !previous_is_invalid) {
                assert(previous_seed_index < seeds->size());
                assert(current_item.get_value() < seeds->size());
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "Comparing seeds " << seeds->at(previous_seed_index).pos << " and " 
            << seeds->at(current_item.get_value()).pos << endl;
#endif

                //Comparator returning previous_seed_index < current_item.value
                size_t depth = 0;

                //Keep track of the orientation of each seed
                //Everything should be sorted according to
                //the orientation in the top-level structure,
                //so if things are traversed backwards, reverse the orientation
                bool a_is_reversed = false;
                bool b_is_reversed = false;
                const Seed& previous_seed = seeds->at(previous_seed_index);
                const Seed& current_seed = seeds->at(current_item.get_value());
                while (depth < previous_seed.zipcode.max_depth()
                       && depth < current_seed.zipcode.max_depth()
                       && ZipCode::is_equal(previous_seed.zipcode, current_seed.zipcode, depth)) {

                    //Remember the orientation
                    if (ZipCodeTree::seed_is_reversed_at_depth(previous_seed, depth, distance_index)) { 
                        a_is_reversed = !a_is_reversed;
                    }
                    if (ZipCodeTree::seed_is_reversed_at_depth(current_seed, depth, distance_index)) {
                        b_is_reversed = !b_is_reversed;
                    }

                    depth++;
                }

                //Remember the orientation of the parent too
                size_t parent_of_a_is_reversed = a_is_reversed;

                //Check the orientations one last time
                if (ZipCodeTree::seed_is_reversed_at_depth(previous_seed, depth, distance_index)) { 
                    a_is_reversed = !a_is_reversed;
                }
                if (ZipCodeTree::seed_is_reversed_at_depth(current_seed, depth, distance_index)) {
                    b_is_reversed = !b_is_reversed;
                }
                
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "\t different at depth " << depth << endl;
#endif
                //Either depth is last thing in previous_seed or current_seed
                //or they are different at this depth


                if ( ZipCode::is_equal(previous_seed.zipcode, current_seed.zipcode, depth)) {
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "\tthey are on the same node" << endl;
#endif
                    //If they are equal, then they must be on the same node

                    size_t offset1 = is_rev(previous_seed.pos)
                                      ? previous_seed.zipcode.get_length(depth) - offset(previous_seed.pos)
                                      : offset(previous_seed.pos);
                    size_t offset2 = is_rev(current_seed.pos) 
                                     ? current_seed.zipcode.get_length(depth) - offset(current_seed.pos)
                                     : offset(current_seed.pos);
                    if (!a_is_reversed) {
                        //If they are in previous_seed_index snarl 
                        //or they are facing forward on a chain, then order by
                        //the offset in the node
                        assert( offset1 <= offset2);
                    } else {
                        //Otherwise, the node is facing backwards in the chain, 
                        //so order backwards in node
                        assert( offset2 <= offset1);
                    }
                }  else if (depth == 0) {
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "\tThey are on different connected components" << endl;
#endif
                    //If they are on different connected components, 
                    //sort by connected component
                    assert(previous_seed.zipcode.get_distance_index_address(0)
                           <= current_seed.zipcode.get_distance_index_address(0));
                    
                }  else if (previous_seed.zipcode.get_code_type(depth-1) == ZipCode::CHAIN 
                            || previous_seed.zipcode.get_code_type(depth-1) == ZipCode::ROOT_CHAIN) {
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "\t they are children of a common chain" << endl;
#endif
                    //If previous_seed and current_seed are children of a chain
                    size_t component_a = previous_seed.zipcode.get_chain_component(depth);
                    size_t component_b = current_seed.zipcode.get_chain_component(depth);
                    size_t offset_a = previous_seed.zipcode.get_offset_in_chain(depth);
                    size_t offset_b = current_seed.zipcode.get_offset_in_chain(depth);
                    if (!current_is_in_cyclic_snarl) {

                        if (component_a == component_b) {
                            if ( offset_a == offset_b) {
                                //If they have the same prefix sum,
                                //then the snarl comes first
                                //Will never be on the same child at this depth
                                if (parent_of_a_is_reversed) {
                                    assert(current_seed.zipcode.get_code_type(depth) != ZipCode::NODE);
                                    assert(previous_seed.zipcode.get_code_type(depth) == ZipCode::NODE); 
                                } else {
                                    assert(previous_seed.zipcode.get_code_type(depth) != ZipCode::NODE);
                                    assert(current_seed.zipcode.get_code_type(depth) == ZipCode::NODE); 
                                }
                            } else {
                                //Check if the parent chain is reversed
                                //and if so, then the order should be reversed
                                //Could happen in irregular snarls
                                if (parent_of_a_is_reversed) {
                                    assert( offset_b <= offset_a);
                                } else {
                                    assert( offset_a <= offset_b);
                                }
                            }
                        } else {
                            if (parent_of_a_is_reversed) {
                                assert( component_b <= component_a);
                            } else {
                                assert( component_a <= component_b);
                            }
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
                    if (!current_is_in_cyclic_snarl) {
                        assert(previous_seed.zipcode.get_rank_in_snarl(depth)
                               <= current_seed.zipcode.get_rank_in_snarl(depth));
                    }
                } 

            }
            previous_seed_index = current_item.get_value();
            previous_is_invalid = current_is_invalid;
        } else if (current_item.get_type() == ZipCodeTree::CHAIN_START) {
            //Chains can't start with edges
            assert(zip_code_tree[i+1].get_type() != ZipCodeTree::EDGE);
        } else if (current_item.get_type() == ZipCodeTree::CHAIN_END) {
            //And can't end with edges
            assert(zip_code_tree[i-1].get_type() != ZipCodeTree::EDGE);
        } else if (current_item.is_snarl_start()) {
            //Sarls start with their node counts
            assert(zip_code_tree[i+1].get_type() == ZipCodeTree::NODE_COUNT);
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
    //Walk from the start of the zip tree, checking each pair of seeds
    std::stack<size_t> chain_numbers;
    ZipCodeTree::seed_iterator dest = begin();
    bool right_to_left = true;
    while (dest != end()) {
#ifdef debug_parse
    std::cerr << "-----------------------------------------" << std::endl;
    std::cerr << (right_to_left ? "Right to left" : "Left to right") << std::endl;
#endif
        //The seed that the iterator points to
        const Seed& start_seed = seeds->at((*dest).seed);

        //Do we want the distance going left in the node
        //This takes into account the position & orientation of tree traversal
        bool start_is_reversed = (*dest).is_reverse ? !is_rev(start_seed.pos) 
                                                    : is_rev(start_seed.pos);

        // The distance between any pair of seeds is the minimum distance

        //Walk through the tree starting from dest and check the distance
        distance_iterator distance_itr_start = find_distances(dest, right_to_left);
        for (distance_iterator tree_itr_left = distance_itr_start; !tree_itr_left.done(); ++tree_itr_left) {

            seed_result_t next_seed_result = *tree_itr_left;
            const Seed& next_seed = seeds->at(next_seed_result.seed);
            const bool next_is_reversed = next_seed_result.is_reverse ? !is_rev(next_seed.pos) 
                                                                      : is_rev(next_seed.pos);

            size_t tree_distance = next_seed_result.distance;

            net_handle_t start_handle = distance_index.get_node_net_handle(
                                            id(start_seed.pos),
                                            is_rev(start_seed.pos) != start_is_reversed);
            net_handle_t next_handle = distance_index.get_node_net_handle(
                                            id(next_seed.pos),  
                                            is_rev(next_seed.pos) != next_is_reversed);

            size_t index_distance = distance_index.minimum_distance(
                id(next_seed.pos), is_rev(next_seed.pos), offset(next_seed.pos),
                id(start_seed.pos), is_rev(start_seed.pos), offset(start_seed.pos), true);

            if (index_distance != std::numeric_limits<size_t>::max() && is_rev(next_seed.pos) != next_is_reversed) {
                //If the seed we're starting from got reversed, then subtract 1
                index_distance -= 1;
            }
            if (index_distance != std::numeric_limits<size_t>::max() && is_rev(start_seed.pos) != start_is_reversed) {
                //If the seed we ended at got reversed, then add 1
                index_distance += 1;
            }
            pos_t start_pos = is_rev(start_seed.pos) 
                            ? make_pos_t(id(start_seed.pos), 
                                         false, 
                                         distance_index.minimum_length(start_handle) - offset(start_seed.pos) )
                            : start_seed.pos;
            pos_t next_pos = is_rev(next_seed.pos) 
                           ? make_pos_t(id(next_seed.pos), 
                                        false, 
                                        distance_index.minimum_length(next_handle) - offset(next_seed.pos) )
                           : next_seed.pos;
            size_t start_length = distance_index.minimum_length(start_handle);
            size_t next_length = distance_index.minimum_length(next_handle);

            bool distance_is_invalid = node_is_invalid(id(next_seed.pos), distance_index, distance_limit) ||
                               node_is_invalid(id(start_seed.pos), distance_index, distance_limit);

            if (!distance_is_invalid && index_distance <= distance_limit) {
                if (start_pos == next_pos) {
                    if (tree_distance != 0 && tree_distance != index_distance) {
                        for (auto& seed : *seeds) {
                            cerr << seed.pos << endl;
                        }
                        cerr << "Distance between " << next_seed.pos << (next_is_reversed ? "rev" : "") 
                             << " and " << start_seed.pos << (start_is_reversed ? "rev" : "") << endl;
                        cerr << "Forward positions: " << start_pos << " " << next_pos 
                             << " and length " << start_length << endl;
                        cerr << "Tree distance: " << tree_distance << " index distance: " << index_distance << endl;
                        cerr << "With distance limit: " << distance_limit << endl;
                    }
                    //This could be off by one if one of the seeds is reversed,
                    //but I'm being lazy and just checking against the index
                    assert((tree_distance == 0 || tree_distance == index_distance));
                } else {
                    if (tree_distance != index_distance) {
                        for (auto& seed : *seeds) {
                            cerr << seed.pos << endl;
                        }
                        cerr << "Distance between " << next_seed.pos << (next_is_reversed ? "rev" : "") 
                             << " and " << start_seed.pos << (start_is_reversed ? "rev" : "") << endl;
                        cerr << "Forward positions: " << start_pos << " " << next_pos << " and lengths " 
                             << start_length << " " << next_length << endl;
                        cerr << "Tree distance: " << tree_distance << " index distance: " << index_distance << endl;
                        cerr << "With distance limit: " << distance_limit << endl;
                    }
                    assert(tree_distance == index_distance);
                }
            }

        }
        if (dest.cyclic_snarl_nested_depth > 0) {
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
    vector<bool> has_seed (seeds->size(), false);
    for (const auto& tree : trees) {
        tree.validate_zip_tree(distance_index, seeds, distance_limit);
        for (size_t i = 0 ; i < tree.zip_code_tree.size() ; i++) {
            const tree_item_t& item = tree.zip_code_tree[i];
            if (item.get_type() == ZipCodeTree::SEED) {
                has_seed[item.get_value()] = true;
            }
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
    assert(zip_iterator->get_type() == ZipCodeTree::NODE_COUNT);
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

    // Read chains
    vector<pos_t> positions;
    // Placeholder for snarl start
    positions.emplace_back(make_pos_t(0, false, 0));
    while (!zip_iterator->is_snarl_end()) {
        assert(zip_iterator->get_type() == ZipCodeTree::CHAIN_START);
        assert(zip_iterator->get_value() == snarl_id);

        // Only need to store first seeds for cyclic snarls
        if (is_cyclic_snarl) {
            // Skip forward to first child
            zip_iterator++;
            store_seed_position(*zip_iterator, distance_index, seeds, positions, true);
            // Back up to start of the chain
            zip_iterator--;
        }

        validate_chain(zip_iterator, distance_index, seeds, distance_limit);

        // Back up to last child
        zip_iterator--;
        store_seed_position(*zip_iterator, distance_index, seeds, positions, false);
        // Skip to next chain start
        zip_iterator += 2;
    }
    assert(zip_iterator->get_type() == is_cyclic_snarl ? ZipCodeTree::CYCLIC_SNARL_END 
                                                       : ZipCodeTree::DAG_SNARL_END);
    assert(zip_iterator->get_value() == snarl_id);

    // Placeholder for snarl end
    positions.emplace_back(make_pos_t(0, false, 0));
    assert(positions.size() == is_cyclic_snarl ? node_count * 2
                                               : node_count + 1);
    
    validate_distance_matrix(distance_index, dist_matrix, positions, is_cyclic_snarl, distance_limit);
}

void ZipCodeTree::store_seed_position(tree_item_t child, 
                                      const SnarlDistanceIndex& distance_index, 
                                      const vector<Seed>* seeds,
                                      std::vector<pos_t>& positions,
                                      bool reverse) const {
    // Is this actually a seed?
    if (child.get_type() == SEED) {
        pos_t seed_pos = seeds->at(child.get_value()).pos;

        // is reversed in the ziptree XOR is first seed of a cyclic snarl chain
        if (child.get_is_reversed() != reverse) {
            seed_pos = make_pos_t(id(seed_pos), !is_rev(seed_pos),
                                  distance_index.minimum_length(
                                    distance_index.get_node_net_handle(id(seed_pos))) - offset(seed_pos));
        }

        positions.push_back(seed_pos);
    } else {
        // Can't verify if the current child is a snarl
        positions.emplace_back(make_pos_t(0, false, 0));
    }
}

void ZipCodeTree::validate_distance_matrix(const SnarlDistanceIndex& distance_index,
                                           const std::vector<size_t>& dist_matrix,
                                           const std::vector<pos_t>& positions,
                                           bool has_self_loops,
                                           size_t distance_limit) const {
#ifdef DEBUG_ZIP_CODE_TREE
    std::cerr << "Validating distance matrix" << std::endl;
#endif
    // Check that the distance matrix is the right size
    if (has_self_loops){
        assert(dist_matrix.size() == (positions.size() * (positions.size() + 1)) / 2);
    } else {
        assert(dist_matrix.size() == (positions.size() * (positions.size() - 1)) / 2);
    }

#ifdef DEBUG_ZIP_CODE_TREE
    for (const auto pos : positions) {
        cerr << "pos: " << id(pos) << (is_rev(pos) ? " rev " : " ")  << offset(pos) << endl;
    }
#endif

    // Current positions
    pos_t from_pos, to_pos;
    size_t matrix_i = 0;
    // Handles for calculating distances
    net_handle_t to_handle, from_handle, parent_handle;
    // Distances being compared
    size_t matrix_distance, true_distance;

    // Check distances between all pairs of seeds
    for (size_t i = 0; i < positions.size(); i++) {
        from_pos = positions[i];
        if (id(from_pos) == 0) {
            // Skip the placeholder
            continue;
        }

        for (size_t j = 0; j < (has_self_loops ? i + 1 : i); j++) {
            to_pos = positions[j];
            if (id(to_pos) == 0) {
                // Skip the placeholder
                continue;
            }

            to_handle = distance_index.get_node_net_handle(id(to_pos), is_rev(to_pos));
            from_handle = distance_index.get_node_net_handle(id(from_pos), is_rev(from_pos));

            // Find appropriate shared parent
            if (has_self_loops) {
                // Assume parent is to_pos's chain
                parent_handle = distance_index.get_parent(to_handle);
                if ((i + 1) / 2 != (j + 1) / 2) {
                    // Inter-chain connection, so go up one level to the snarl
                    parent_handle = distance_index.get_parent(to_handle);
                }
            } else {
                // Non-cyclic snarls work on the level of chains, not seeds
                to_handle = distance_index.get_parent(to_handle);
                from_handle = distance_index.get_parent(from_handle);
                parent_handle = distance_index.get_parent(to_handle);
            }

            true_distance = distance_index.distance_in_parent(parent_handle, to_handle, from_handle);
            // Add distance within node for to_handle
            true_distance = SnarlDistanceIndex::sum(true_distance, 
                distance_index.minimum_length(to_handle) - offset(to_pos));
            // Add distance within node for from_handle
            true_distance = SnarlDistanceIndex::sum(true_distance, 
                distance_index.minimum_length(from_handle) - offset(from_pos));

            matrix_distance = dist_matrix[matrix_i++];

            if (true_distance != matrix_distance && true_distance < distance_limit) {
                cerr << "Distance mismatch between " << from_pos << " and " << to_pos << endl;
                cerr << "True distance: " << true_distance << ", Matrix distance: " << matrix_distance << endl;
                cerr << "With distance limit: " << distance_limit << endl;
                assert(true_distance == matrix_distance);
            }
        }
    }
}

void ZipCodeTree::validate_chain(vector<tree_item_t>::const_iterator& zip_iterator, 
                                  const SnarlDistanceIndex& distance_index, 
                                  const vector<Seed>* seeds,
                                  size_t distance_limit) const {
#ifdef DEBUG_ZIP_CODE_TREE
    std::cerr << "Validating chain" << std::endl;
#endif
    size_t last_seed_i = std::numeric_limits<size_t>::max();
    // Skip past chain start
    zip_iterator++;

    while (zip_iterator->get_type() != ZipCodeTree::CHAIN_END) {
        if (zip_iterator->get_type() == SEED) {
            last_seed_i = zip_iterator->get_value();
        } else if (zip_iterator->get_type() == EDGE) {
            // Edges within chains are validated by the zip iterator checks
        } else if (zip_iterator->is_snarl_start()) {
            // Validate nested snarl
            validate_snarl(zip_iterator, distance_index, seeds, distance_limit);
        }
        zip_iterator++;
    }
}

ZipCodeTree::seed_iterator::seed_iterator(vector<tree_item_t>::const_iterator begin, 
    vector<tree_item_t>::const_iterator end, const vector<tree_item_t>& zip_code_tree,
    const unordered_map<size_t, size_t>& snarl_start_indexes) : it(begin), end(end), zip_code_tree(zip_code_tree),
    snarl_start_indexes(snarl_start_indexes), cyclic_snarl_nested_depth(0), chain_numbers(std::stack<size_t>()) {
    
    // Immediately advance to the first seed
    while (this->it != this->end && this->it->get_type() != SEED) {
        ++(*this);
    }
}

auto ZipCodeTree::seed_iterator::operator++() -> seed_iterator& {
    ++it;
    while (it != end && it->get_type() != SEED) {
        // cyclic_snarl_nested_depth remembers if we're in a cyclic snarl
        if (it->get_type() == ZipCodeTree::CYCLIC_SNARL_START) {
            cyclic_snarl_nested_depth++;
        } else if (it->get_type() == ZipCodeTree::CYCLIC_SNARL_END) {
            cyclic_snarl_nested_depth--;
        }

        // chain_numbers remembers which chain we're in for each snarl
        if (it->is_snarl_start()) {
            chain_numbers.push(0);
        } else if (it->is_snarl_end()) {
            chain_numbers.pop();
        } else if (it->get_type() == ZipCodeTree::CHAIN_START) {
            chain_numbers.top()++;
        }

        // Advance to the next seed, or the end.
        ++it;
    }
    return *this;
}

auto ZipCodeTree::seed_iterator::operator==(const seed_iterator& other) const -> bool {
    // Ends don't matter for comparison.
    return it == other.it;
}
    
auto ZipCodeTree::seed_iterator::operator*() const -> oriented_seed_t {
    return {it->get_value(), it->get_is_reversed()};
}

auto ZipCodeTree::seed_iterator::remaining_tree() const -> size_t {
    size_t to_return = end - it - 1;
#ifdef debug_parse
    std::cerr << "From " << &*it << " there are " << to_return << " slots after" << std::endl;
#endif
    return to_return;
}

auto ZipCodeTree::begin() const -> seed_iterator {
    return seed_iterator(zip_code_tree.begin(), zip_code_tree.end(), zip_code_tree, snarl_start_indexes);
}

auto ZipCodeTree::end() const -> seed_iterator {
    return seed_iterator(zip_code_tree.end(), zip_code_tree.end(), zip_code_tree, snarl_start_indexes);
}

ZipCodeTree::distance_iterator::distance_iterator(vector<tree_item_t>::const_reverse_iterator rbegin, 
    vector<tree_item_t>::const_reverse_iterator rend, const vector<tree_item_t>& zip_code_tree,
    const unordered_map<size_t, size_t>& snarl_start_indexes, std::stack<size_t> chain_numbers, 
    bool right_to_left, size_t distance_limit) : 
    it(rbegin), rend(rend), origin(rbegin), zip_code_tree(zip_code_tree), snarl_start_indexes(snarl_start_indexes),
    chain_numbers(chain_numbers), right_to_left(right_to_left),
    distance_limit(distance_limit), stack_data(nullptr), current_state(S_START) {
#ifdef debug_parse
    if (this->it != rend) {
        std::cerr << "Able to do first initial tick." << std::endl;
    }
#endif
    if (this->it == rend) {
        // We are an end iterator. Nothing else to do.
        return;
    }
    while (this->it != rend && !tick()) {
        // Skip ahead to the first seed we actually want to yield,
        //or to the end of the data.
        ++(*this);
#ifdef debug_parse
        if (this->it != rend) {
            std::cerr << "Able to do another initial tick." << std::endl;
        }
#endif
    }
    // As the end of the constructor, the iterator points to
    // a seed that has been ticked and yielded, or is rend.
#ifdef debug_parse
    if (this->it == rend) {
        std::cerr << "Tree iteration halted looking for first seed." << std::endl;
    }
#endif
}

ZipCodeTree::distance_iterator::distance_iterator(const distance_iterator& other) : 
    it(other.it), rend(other.rend), origin(other.origin), distance_limit(other.distance_limit), chain_numbers(std::move(other.chain_numbers)),
    stack_data(other.stack_data ? new std::stack<size_t>(*other.stack_data) : nullptr), 
    right_to_left(other.right_to_left), zip_code_tree(other.zip_code_tree),
    snarl_start_indexes(other.snarl_start_indexes), current_state(other.current_state) {
    // Nothing to do!
}

ZipCodeTree::distance_iterator::distance_iterator(distance_iterator&& other) : 
    it(std::move(other.it)), rend(std::move(other.rend)), origin(other.origin), chain_numbers(std::move(other.chain_numbers)),
    distance_limit(std::move(other.distance_limit)), stack_data(std::move(other.stack_data)),
    right_to_left(other.right_to_left), zip_code_tree(std::move(other.zip_code_tree)),
    snarl_start_indexes(std::move(other.snarl_start_indexes)), current_state(std::move(other.current_state)) {
    // Nothing to do!
}

auto ZipCodeTree::distance_iterator::operator=(const distance_iterator& other) -> distance_iterator& {
    it = other.it;
    rend = other.rend;
    origin = other.origin;
    chain_numbers = other.chain_numbers;
    distance_limit = other.distance_limit;
    stack_data.reset(other.stack_data ? new std::stack<size_t>(*other.stack_data) : nullptr);
    right_to_left = other.right_to_left;
    current_state = other.current_state;
    return *this;
}

auto ZipCodeTree::distance_iterator::operator=(distance_iterator&& other) -> distance_iterator& {
    it = std::move(other.it);
    rend = std::move(other.rend);
    origin = std::move(other.origin);
    chain_numbers = std::move(other.chain_numbers);
    distance_limit = std::move(other.distance_limit);
    stack_data = std::move(other.stack_data);
    right_to_left = std::move(other.right_to_left);
    current_state = std::move(other.current_state);
    return *this;
}

auto ZipCodeTree::distance_iterator::operator++() -> distance_iterator& {
    // Invariant: the iterator points to
    // a seed that has been ticked and yielded, or to rend.
    if (it != rend) {
#ifdef debug_parse
    std::cerr << "Skipping over a " << it->get_type() << " which we assume was handled already." << std::endl;
#endif
        if (right_to_left) {
            ++it;
        } else {
            --it;
        }

    }
    while (it != rend && !tick()) {
        // Skip ahead to the next seed we actually want to yield,
        // or to the end of the data.
        if (right_to_left) {
            ++it;
        } else {
            --it;
        }
    }
#ifdef debug_parse
    if (it == rend) {
        std::cerr << "Tree iteration halted looking for next seed." << std::endl;
    }
#endif
    return *this;
}

auto ZipCodeTree::distance_iterator::operator==(const distance_iterator& other) const -> bool {
    // Ends and other state don't matter for comparison.
    return it == other.it;
}

auto ZipCodeTree::distance_iterator::operator*() const -> seed_result_t {
    // We are always at a seed, so show that seed
#ifdef check_parse
    crash_unless(it != rend);
    crash_unless(it->get_type() == SEED);
    crash_unless(stack_data);
    crash_unless(!stack_data->empty());
#endif
    // We know the running distance to this seed will be at the top of the stack
    seed_result_t to_return;
    to_return.seed = it->get_value();
    cerr << (right_to_left ? "Right to left" : "Left to right") << " distance to seed " 
         << to_return.seed << " is " << stack_data->top() << endl;
    to_return.is_reverse = right_to_left ? it->get_is_reversed()
                                         : !it->get_is_reversed();
    to_return.distance = (it == origin) ? 0 : stack_data->top();
    return to_return;
}

auto ZipCodeTree::distance_iterator::push(size_t value) -> void {
    stack().push(value);
}

auto ZipCodeTree::distance_iterator::pop() -> size_t {
    size_t value = stack().top();
    stack().pop();
    return value;
}

auto ZipCodeTree::distance_iterator::top() -> size_t& {
#ifdef check_parse
    crash_unless(depth() > 0);
#endif
    return stack().top();
}

auto ZipCodeTree::distance_iterator::dup() -> void {
    push(stack().top());
}

auto ZipCodeTree::distance_iterator::depth() const -> size_t {
    if (!stack_data) {
        return 0;
    } else {
        return stack_data->size();
    }
}

auto ZipCodeTree::distance_iterator::swap() -> void {
    // Grab the top item
    size_t temp = stack().top();
    stack().pop();
    // Swap it with what was under it
    std::swap(temp, stack().top());
    // And put that back on top
    stack().push(temp);
}

vector<size_t> ZipCodeTree::distance_iterator::get_distances_from_chain(size_t snarl_id, size_t chain_num, 
                                                                        bool right_side) const {
    // Read snarl header
    bool is_cyclic = snarl_is_cyclic(snarl_id);
    size_t dist_matrix_start = snarl_start_indexes.at(snarl_id) + 2;
    tree_item_t node_count = zip_code_tree[dist_matrix_start - 1];
    if (node_count.get_type() == ZipCodeTree::CHAIN_START) {
        // This must be a root snarl, without a distance matrix
        return vector<size_t>();
    }
    size_t num_chains = node_count.get_value();
#ifdef debug_parse
    cerr << "Get distances for snarl " << snarl_id << " with " << num_chains << " chain(s); "
         << "stacking for chain " << chain_num << "'s "
         << (right_side ? "right" : "left") << " side" << endl;
    assert(chain_num <= num_chains || chain_num == std::numeric_limits<size_t>::max());
#endif

    vector<size_t> distances;

    if (is_cyclic) {
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

        if (right_side) {
            // Heading towards the snarl end
            distances.push_back(get_matrix_value(dist_matrix_start, true, cur_row, num_chains * 2 + 1));
        } else {
            // Heading towards the snarl start
            distances.push_back(get_matrix_value(dist_matrix_start, true, cur_row, 0));
        }

        // Cyclic snarl always stacks all lefts on top of then all rights 
        // Due to "bouncing", c1_L is on top, and c1_R is on bottom
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
            // DAG snarl, going left to right: get dists to all chains to right
            for (size_t i = num_chains + 1; i > chain_num; i--) {
                distances.push_back(get_matrix_value(dist_matrix_start, false, chain_num, i));
            }
        } else {
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


auto ZipCodeTree::distance_iterator::state(State new_state) -> void {
    current_state = new_state;
}

bool ZipCodeTree::distance_iterator::entered_snarl() const {
    return (right_to_left && it->is_snarl_end())
           || (!right_to_left && it->is_snarl_start());
}

bool ZipCodeTree::distance_iterator::exited_snarl() const {
    return (right_to_left && it->is_snarl_start())
           || (!right_to_left && it->is_snarl_end());
}

bool ZipCodeTree::distance_iterator::entered_chain() const {
    return (right_to_left && it->get_type() == ZipCodeTree::CHAIN_END)
           || (!right_to_left && it->get_type() == ZipCodeTree::CHAIN_START);
}

bool ZipCodeTree::distance_iterator::exited_chain() const {
    return (right_to_left && it->get_type() == ZipCodeTree::CHAIN_START)
           || (!right_to_left && it->get_type() == ZipCodeTree::CHAIN_END);
}

void ZipCodeTree::distance_iterator::skip_chain() {
    // Tracker for how nested we are in snarls
    // We are only allowed to stop skipping when at the same nestedness
    push(0);
    state(S_SKIP_CHAIN);
}

void ZipCodeTree::distance_iterator::initialize_chain() {
    if (top() > distance_limit || top() == std::numeric_limits<size_t>::max()) {
#ifdef debug_parse
    std::cerr << "Skip chain " << it->get_value() << " with running distance "
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
    // By default, after stacking up distances we need to find the first chain
    State final_state = S_FIND_CHAIN;
    bool original_right_to_left = right_to_left;
    size_t chain_nestedness; 
    if (snarl_is_cyclic(it->get_value())) {
        // Memorize previous direction
        push(right_to_left ? 1 : 0);
        if (right_to_left) {
            // Cyclic snarl R->L needs to skip the chains to get
            // to starting spot; thus, scan to find distance matrix
            final_state = S_FIND_DIST_MATRIX;
            chain_nestedness = 0;
        } else if (!right_to_left && chain_num != 0) {
            // Cyclic snarl L->R needs to turn around to get
            // back to the distance matrix
            final_state = S_FIND_DIST_MATRIX;
            right_to_left = true;
            // We're inside a chain already
            chain_nestedness = 1;
        }
    }

    vector<size_t> distances = get_distances_from_chain(it->get_value(), chain_num, !original_right_to_left);
    if (distances.empty()) {
#ifdef debug_parse
    std::cerr << "Tried to stack up distances for a root-level snarl; halting now." << std::endl;
#endif
        halt();
        return false;
    }

    // Add distances to chain's running distance
    for (const auto& distance : distances) {
        dup();
        // Add in the edge value to make a running distance
        // for the thing this edge is for.
        // Account for if the edge is actually unreachable.
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

    state(final_state);
    if (final_state == S_FIND_DIST_MATRIX) {
        // Chain nestedness to emerge from
        push(chain_nestedness);
    }
    return true;
}

void ZipCodeTree::distance_iterator::continue_snarl() {
    // Different scanning states based on snarl type
    if (snarl_is_cyclic(it->get_value())) {
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
    it = rend;
}

void ZipCodeTree::distance_iterator::unimplemented_error() {
    throw std::domain_error("Unimplemented symbol " + std::to_string(it->get_type()) 
                                    + " for state " + std::to_string(current_state));
}

auto ZipCodeTree::distance_iterator::tick() -> bool {
#ifdef debug_parse
    std::cerr << "Tick for state " << current_state << " on symbol " << it->get_type() 
              << " at entry " << (rend - it + 1) << std::endl;
    std::cerr << "Stack: ";
    if (stack_data) {
        std::stack<size_t> temp_stack = *stack_data;
        while (!temp_stack.empty()) {
            std::cerr << (temp_stack.top() == std::numeric_limits<size_t>::max() ? "inf" 
                                                                                 : std::to_string(temp_stack.top())) 
                      << " ";
            temp_stack.pop();
        }
    } else {
        std::cerr << "empty";
    }
    std::cerr << std::endl;
#endif
    switch (current_state) {
    case S_START:
        // Initial state.
        //
        // Stack is empty and we must be at a seed to start at.
        if (it->get_type() == SEED) {
#ifdef debug_parse
    std::cerr << "Skip over seed " << it->get_value() << std::endl;
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
        if (it->get_type() == SEED) {
            // Emit seed here with distance at top of stack.
#ifdef check_parse
    crash_unless(depth() > 0);
#endif
#ifdef debug_parse
    std::cerr << "Yield seed " << it->get_value() << ", distance " << top() << std::endl;
#endif
            return true;
        } else if (entered_snarl()) {
            // Running distance along chain is on stack,
            // and will need to be added to all the stored distances.
            return !initialize_snarl(it->is_snarl_start() ? 0 : std::numeric_limits<size_t>::max());
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
                // Running distance for next chain,
                // or running distance to cross the snarl, will be under it.
                continue_snarl();
            }
        } else if (it->get_type() == EDGE) {
            // Distance between things in a chain.
            // Add value into running distance, maxing it if value is max.
            top() = SnarlDistanceIndex::sum(top(), it->get_value());
            if (top() > distance_limit || top() == std::numeric_limits<size_t>::max()) {
                // Skip over the rest of this chain
                if (depth() == 1) {
                    // We never entered the parent snarl of this chain.
                    // So if the distance along the chain is too much, there
                    // are not going to be any results with a smaller distance.
#ifdef debug_parse
    std::cerr << "Halt because adding " << it->get_value() << " bp " 
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
    case S_FIND_CHAIN:
        // State where we are finding the first chain in a snarl.
        // Distances have already been stacked up by initialize_snarl()
        //
        // Stack has the running distance along the parent chain, and over that
        // that the stacked running distances for items in the snarl.
        if (it->get_type() == EDGE || it->get_type() == NODE_COUNT) {
            // Skip over distance matrix; we have already used it
        } else if (entered_chain()) {
            initialize_chain();
        } else if (exited_snarl()) {
            // We didn't hit another chain in the snarl, we hit the start of
            // the snarl. We stacked one distance, to the snarl end.
            if (depth() > 1) {
                // We have a distance for the parent chain, so toss out the
                // inner chain's distance
                pop();
            }
            state(S_SCAN_CHAIN);
        } else {
            unimplemented_error(); 
        }
        break;
    case S_SCAN_DAG_SNARL:
        // State where we are going through a DAG snarl having already
        // stacked up the running distances to use for each chain in the snarl.
        //
        // Stack has at the top running distances to use for each chain still
        // to be visited in the snarl, and under those the same for the snarl
        // above that, etc.
        if (it->get_type() == EDGE || it->get_type() == NODE_COUNT) {
            // Skip over distance matrix; we have already used it
        } else if (exited_snarl()) {
            // Stack holds running distance along parent chain plus edge
            // distance to cross the snarl, or running distance out of chain we
            // started in plus distance to exit the snarl.
            //
            // This is the running distance to use for the parent chain now.
            // So go back to scanning the parent chain.
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
        if (it->get_type() == CYCLIC_SNARL_END) {
            // Finished left sides, now doing right sides.
            right_to_left = !right_to_left;
        } else if (it->get_type() == EDGE) {
#ifdef debug_parse
    std::cerr << "Finished cyclic snarl, resetting direction to "
              << (top() == 1 ? "right to left" : "left to right" ) << std::endl;
#endif
            // Hit distance matrix again; finished snarl
            // The top of the stack will be the original direction we were going
            right_to_left = (pop() == 1);
            // Current chain nesting level, to remember what we're skipping over
            push(0);
            state(S_SKIP_SNARL);
        } else if (entered_chain()) {
            // We've encountered a chain to look at, and the running distance
            // into the chain is already on the stack.
            initialize_chain();
        } else {
            unimplemented_error();
        }
        break;
    case S_SKIP_CHAIN:
        // State where we are skipping over the rest of a chain because we hit
        // the distance limit, but we might need to do other chains in a parent
        // snarl.
        //
        // Stack has the nesting level of child snarls we are reading over
        // until we get back to the level we want to skip past the chain start.
        //
        // Under that is the running distance along the chain being skipped.
        // And under that it has the running distance for ther next thing in
        // the snarl, which had better exist or we shouldn't be trying to skip
        // the chain, we should have halted.
        if (it->get_type() == SEED || it->get_type() == EDGE || it->get_type() == NODE_COUNT) {
            // We don't emit seeds until the chain is over,
            // and we ignore edges/distance matrices
        } else if (exited_snarl()) {
            // We might now be able to match chains again
            top() -= 1;
        } else if (entered_snarl()) {
            // We can't match chains until we leave the snarl
            top() += 1;
        } else if (exited_chain()) {
            if (top() == 0) {
                // This is the other end of the chain we were wanting to skip.
                // Discard the snarl nestedness level
                pop();
#ifdef check_parse
    crash_unless(depth() >= 1);
#endif
                // Discard the running distance along this chain
                pop();
                // Running distance for next chain,
                // or running distance to cross the snarl, will be under it.
                continue_snarl();
            }
            // Otherwise this is a chain inside a child snarl and we ignore it.
        } else if (entered_chain()) {
            // Nested chain that we're skipping over
        } else {
            unimplemented_error(); 
        }
        break;
    case S_SKIP_SNARL:
        // State where we are skipping over a cyclic snarl to get to the
        // correct exit, after processing all chains
        //
        // Stack has the nesting level of child chains we are reading over
        // until we get back to the level we want to skip past the snarl start.
        //
        // Under that is the running distance along the parent chain.
        if (it->get_type() == SEED || it->get_type() == EDGE || it->get_type() == NODE_COUNT) {
            // We don't emit seeds until the snarl is over,
            // and we ignore edges/distance matrices
        } else if (exited_chain()) {
            // We might now be able to match snarls again
            top() -= 1;
        } else if (entered_chain()) {
            // We can't match snarls until we leave the chain
            top() += 1;
        } else if (exited_snarl()) {
            if (top() == 0) {
                // This is the start of the snarl we were wanting to skip.
                // Discard the chain nestedness level
                pop();
#ifdef check_parse
    crash_unless(depth() >= 1);
#endif
                // Running distance for next chain,
                // or running distance to cross the snarl, will be under it.
                state(S_SCAN_CHAIN);
            }
            // Otherwise this is a snarl inside a child chain and we ignore it.
        } else if (entered_snarl()) {
            // Nested snarl that we're skipping over
        } else {
            unimplemented_error(); 
        }
        break;
    case S_FIND_DIST_MATRIX:
        // State where we are finding the distance matrix for a snarl.
        // We are in a cyclic snarl, and we need to find the distance matrix
        // to start scanning the chains.
        //
        // We have already stacked up the running distances to use for each
        // chain if we start on the left, which is where the matrix is.
        if (it->get_type() == SEED || it->get_type() == NODE_COUNT) {
            // We don't emit seeds until the snarl is over
        } else if (entered_snarl() || exited_snarl()) {
            // Child snarl we're skipping over
        } else if (exited_chain()) {
            // We might now be able to match snarls again
            top() -= 1;
        } else if (entered_chain()) {
            // We can't match snarls until we leave the chain
            top() += 1;
        } else if (it->get_type() == EDGE) {
            if (top() == 0) {
#ifdef debug_parse
    std::cerr << "Found distance matrix for cyclic snarl; turning around now" << std::endl;
#endif
                // This is the distance matrix we were wanting to skip to.
                // Discard the chain nestedness level
                pop();
#ifdef check_parse
    crash_unless(depth() >= 1);
#endif
                // Turn around and find the first chain
                right_to_left = false;
                state(S_FIND_CHAIN);
            }
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
    return distance_iterator(zip_code_tree.rbegin() + from.remaining_tree(), 
                             right_to_left ? zip_code_tree.rend() : zip_code_tree.rbegin(), 
                             from.zip_code_tree, from.snarl_start_indexes,
                             from.chain_numbers, right_to_left, distance_limit);
}
auto ZipCodeTree::rend() const -> distance_iterator {
    return distance_iterator(zip_code_tree.rend(), zip_code_tree.rend(), zip_code_tree,
                             snarl_start_indexes, std::stack<size_t>(), 0);
}


std::ostream& operator<<(std::ostream& out, const ZipCodeTree::tree_item_type_t& type) {
    return out << std::to_string(type);
}

std::ostream& operator<<(std::ostream& out, const ZipCodeTree::distance_iterator::State& state) {
    return out << std::to_string(state);
}

void ZipCodeForest::sort_one_interval(forest_growing_state_t& forest_state, 
    const interval_state_t& interval) const {

    vector<size_t>& zipcode_sort_order = forest_state.seed_sort_order; 
    vector<sort_value_t>& sort_values_by_seed = forest_state.sort_values_by_seed;
    const vector<Seed>* seeds = forest_state.seeds;

#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "Sort interval at depth " << interval.depth << (interval.is_reversed ? " reversed" : "") << endl;
#endif



    /*** First, fill in sort_values_by_seed for the relevant seeds ***/

    //This doesn't take into account the orientation,
    //except for nodes offsets in chains
    //Used for sorting at the given depth, so use values at depth depth+1

    //Get the minimum and maximum values that are used for sorting.
    //These will be used to determine if radix sort will be more efficient
    
    //This must be done even if the interval is already sorted,
    //because we need to fill in the sort values

    size_t max_sort_value = 0;
    size_t min_sort_value = std::numeric_limits<size_t>::max();

    //The min and max chain components
    size_t min_component = 0;
    size_t max_component = 0;

    for (size_t i = interval.interval_start ; i < interval.interval_end ; i++) {
        const Seed& seed = seeds->at(zipcode_sort_order[i]); 
#ifdef DEBUG_ZIP_CODE_SORTING
    cerr << "\tGet the sort value of seed " << seed.pos << " at depth " << interval.depth+1 
            << " with parent type " << interval.code_type << endl;
#endif
        if (interval.code_type == ZipCode::EMPTY) {
            // If we are sorting the root int connected components 

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

                //For a snarl, the order is prefix_sum*3+1
                sort_values_by_seed[zipcode_sort_order[i]].set_sort_value(prefix_sum);
                sort_values_by_seed[zipcode_sort_order[i]].set_chain_order(1);
            } else {
                //If this is a node, then the order
                //depends on where the position falls in the node
                bool node_is_rev = seed.zipcode.get_is_reversed_in_parent(interval.depth+1) != is_rev(seed.pos);
                node_is_rev = node_is_rev;
                size_t node_offset = node_is_rev ? seed.zipcode.get_length(interval.depth+1) - offset(seed.pos)
                                                 : offset(seed.pos);

                sort_values_by_seed[zipcode_sort_order[i]].set_sort_value(
                    SnarlDistanceIndex::sum(prefix_sum, node_offset));
                if (node_offset == 0) {
                    sort_values_by_seed[zipcode_sort_order[i]].set_chain_order(2);
                } else {
                    sort_values_by_seed[zipcode_sort_order[i]].set_chain_order(0);
                }
            }
#ifdef DEBUG_ZIP_CODE_SORTING
    cerr << "Prefix sum " << sort_values_by_seed[zipcode_sort_order[i]].get_distance_value() << " and sort value " 
            << sort_values_by_seed[zipcode_sort_order[i]].get_sort_value() << " and type " << child_type << endl;
#endif
        } else {
#ifdef DEBUG_ZIP_CODE_SORTING
    cerr << "\tThis is snarl, so return the rank in the snarl: " 
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
        //If this is a root chain, then use the default sort, 
        //because it's probably too big for radix and we can't tell
        //anyways because we don't store the length of a root-chain
        use_radix = false;
    } else {
        //The cost of default sort is nlog(n) 
        //where n is the number of things to sort
        size_t default_cost = (interval.interval_end - interval.interval_start) 
                            * std::log2(interval.interval_end - interval.interval_start);
        //The cost of radix sort is linear in the number of distinct values 
        //(since we will subtract the minimum)
        size_t radix_cost = max_sort_value - min_sort_value;
        use_radix = radix_cost <= default_cost;
    }

    /**** Sort *********/

    //If this is a multicomponent chain, then sort by component first
    vector<interval_state_t> sub_intervals;
    if (min_component != max_component) {
        sub_intervals.reserve(max_component-min_component);
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "Sort by chain component" << endl;
#endif
        //Sort by component using radix sort. I doubt that there will be 
        //enough components to make it more efficient to use the default sort
        radix_sort_zipcodes(zipcode_sort_order, sort_values_by_seed, interval, 
                            interval.is_reversed, min_component, max_component, true);

        //Now get the next intervals in sub_intervals
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
        //Copy the current interval
        sub_intervals.emplace_back(interval);
    }


    for (auto& sub_interval : sub_intervals) {
        //Snarls are already sorted by a topological order of
        //the orientation of the zip tree, so don't reverse them
        //And don't reverse the sort if that has
        //already been taken into account in the value finding
        bool reverse_order = (sub_interval.code_type == ZipCode::REGULAR_SNARL 
                              || sub_interval.code_type == ZipCode::IRREGULAR_SNARL) 
                             ? false
                             : sub_interval.is_reversed; 
        if (use_radix) {
            //Sort the given interval using the value-getter and orientation
            radix_sort_zipcodes(zipcode_sort_order, sort_values_by_seed, sub_interval, 
                                reverse_order, min_sort_value, max_sort_value);
        } else {
            //Sort the given interval using the value-getter and orientation
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
    

    //New intervals get added to the front of next intervals,
    //in the sort order that they are found in.
    //This means that the first interval found gets added to the front of list,
    //then the next one gets added after that one.
    //insert_itr will always point to the item in front of wherever
    //the next interval should be added,
    //so always emplace/insert_after the instert_itr,
    //and move it forward after inserting
    std::forward_list<ZipCodeForest::interval_state_t>::iterator insert_itr = next_intervals.before_begin();



    /*********   Check for new intervals of the children ****************/

#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "Finding intervals after sorting at depth " << interval.depth << endl;
#endif
    //After sorting, find runs of equivalent values for new_interval_to_sort
    //Everything gets put into a new interval,
    //even if it is the only thing with that partitioning value
    //Since nodes are really just seeds on the same chain,
    //runs of nodes get put together even if they're actually on different nodes
    //as long as the nodes are facing in the same direction
    //Also need to check the orientation

    //max() is used for the root, when the child's depth should be 0
    size_t child_depth = interval.code_type == ZipCode::EMPTY ? 0 : interval.depth+1;


    if (interval.code_type != ZipCode::EMPTY && 
        seeds->at(zipcode_sort_order[interval.interval_start]).zipcode.max_depth() == interval.depth ) {
        //If this is a trivial chain, then return the same interval as a node
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "\tthis was a trivial chain so just return the same interval as a node" << endl;
#endif
        next_intervals.emplace_after(insert_itr, interval.interval_start, interval.interval_end, 
                                     interval.is_reversed, ZipCode::NODE, child_depth);
        return;
    }


    //These get compared to see if the next seeds is in the same interval
    ZipCode::code_type_t first_type = sort_values_by_seed[zipcode_sort_order[interval.interval_start]].get_code_type();

    //This is only for nodes in chains, since anything on nodes in chains
    //are considered just children of the chain
    bool previous_is_node = first_type == ZipCode::NODE;

    //This only matters if it isn't a node
    size_t previous_sort_value = previous_is_node 
                                ? (ZipCodeTree::seed_is_reversed_at_depth(
                                    seeds->at(zipcode_sort_order[interval.interval_start]), 
                                    child_depth, *distance_index) 
                                    ? 1 : 0)
                                : sort_values_by_seed[zipcode_sort_order[interval.interval_start]].get_sort_value();

    //Start the first interval. 
    //The end value and is_reversed gets set when ending the interval
    next_intervals.emplace_after(insert_itr, interval.interval_start, interval.interval_start, interval.is_reversed, 
                               first_type, child_depth);
    ++insert_itr;

    for (size_t i = interval.interval_start+1 ; i < interval.interval_end ; i++) {
        
        //If the current seed is a node and has nothing at depth+1
        //or is different from the previous seed at this depth
        ZipCode::code_type_t current_type = sort_values_by_seed[zipcode_sort_order[i]].get_code_type();
        bool is_node = current_type == ZipCode::NODE;
        size_t sort_value = is_node 
                          ? (ZipCodeTree::seed_is_reversed_at_depth(seeds->at(zipcode_sort_order[i]), 
                                                                    child_depth, *distance_index) ? 1 : 0)
                          : sort_values_by_seed[zipcode_sort_order[i]].get_sort_value();
        bool is_different_from_previous = is_node != previous_is_node ? true : sort_value != previous_sort_value;
        previous_is_node = is_node;
        previous_sort_value = sort_value;

        if (is_different_from_previous) {
            //If this is the end of a run, close the previous run
            //Add its end value and orientation

            insert_itr->interval_end = i;

            
            insert_itr->is_reversed = ZipCodeTree::seed_is_reversed_at_depth(seeds->at(zipcode_sort_order[i-1]), 
                                                                                      child_depth, *distance_index) 
                                    ? !interval.is_reversed
                                    : interval.is_reversed;

            //Open a new run
            next_intervals.emplace_after(insert_itr, i, i, interval.is_reversed, 
                                         is_node ? ZipCode::NODE : current_type, 
                                         child_depth);
            ++insert_itr;
        }
    }

    //Close the last run
    insert_itr->interval_end = interval.interval_end;

    insert_itr->is_reversed = ZipCodeTree::seed_is_reversed_at_depth(
                                seeds->at(zipcode_sort_order[interval.interval_end-1]), child_depth, *distance_index) 
                             ? !interval.is_reversed
                             : interval.is_reversed;
#ifdef DEBUG_ZIP_CODE_SORTING
    cerr << "New sort order " << endl;
    for (auto& interval : new_intervals) {
        for (size_t i = interval.interval_start ; i < interval.interval_end ; i++) {
            cerr << seeds->at(zipcode_sort_order[i]).pos << ", ";
        }
        cerr << "|";
    }
    cerr << endl;
#endif
    return;
}

void ZipCodeForest::radix_sort_zipcodes(vector<size_t>& zipcode_sort_order, 
                                        const vector<sort_value_t>& sort_values_by_seed,
                                        const interval_state_t& interval, bool reverse_order,
                                        size_t min_value, size_t max_value, bool sort_by_chain_component) const {
    //Radix sort the interval of zipcode_sort_order in the given interval
#ifdef DEBUG_ZIP_CODE_SORTING
    cerr << "\tradix sort" << endl;
#endif

    //Mostly copied from Jordan Eizenga

    // count up occurrences of each rank
    std::vector<size_t> counts (max_value-min_value+2, 0);
    for (size_t i  = interval.interval_start ; i < interval.interval_end ; i++) {
        size_t sort_value = sort_by_chain_component ?  sort_values_by_seed[zipcode_sort_order[i]].get_chain_component()
                                                    : sort_values_by_seed[zipcode_sort_order[i]].get_sort_value();
#ifdef DEBUG_ZIP_CODE_SORTING
    assert(sort_value >= min_value);
    assert(sort_value <= max_value);
    cerr << "Sort value for seed " << seeds->at(zipcode_sort_order[i]).pos << ": " 
         << sort_value << endl;
    assert(counts.size() > sort_value - min_value + 1);
#endif
        size_t next_rank = sort_value - min_value + 1;

        ++counts[next_rank];
    }
    
    //Make this a count of the number of things before it
    for (size_t i = 1; i < counts.size(); ++i) {
        counts[i] += counts[i - 1];
    }

    //Get the sorted order
    std::vector<size_t> sorted(interval.interval_end - interval.interval_start);
    for (size_t i = interval.interval_start ; i < interval.interval_end ; i++) {
        size_t sort_value = sort_by_chain_component ?  sort_values_by_seed[zipcode_sort_order[i]].get_chain_component()
                                                    : sort_values_by_seed[zipcode_sort_order[i]].get_sort_value();
        size_t rank = sort_value - min_value;
        sorted[counts[rank]++] = zipcode_sort_order[i];
    }
    
    //And place everything in the correct position
    for (size_t i = 0 ; i < sorted.size() ; i++) {


        //If this is reversed in the top-level chain,
        //then the order should be backwards
        if (reverse_order) {
            zipcode_sort_order[interval.interval_end - i - 1] = sorted[i];
        } else {
            zipcode_sort_order[i + interval.interval_start] = sorted[i];
        }
    }

}
void ZipCodeForest::default_sort_zipcodes(vector<size_t>& zipcode_sort_order, 
                                          const vector<sort_value_t>& sort_values_by_seed,
                                          const interval_state_t& interval, bool reverse_order) const { 
    //std::sort the interval of zipcode_sort_order
    //between interval_start and interval_end
    
#ifdef DEBUG_ZIP_CODE_SORTING
    cerr << "\tdefault sort between " << interval.interval_start  << " and " << interval.interval_end  << endl;
    cerr << "\tis rev: " << reverse_order << endl;
#endif
    //Sort using std::sort 
    std::sort(zipcode_sort_order.begin() + interval.interval_start, 
              zipcode_sort_order.begin() + interval.interval_end, [&] (size_t a, size_t b) {
        //If this snarl tree node is reversed, then reverse the sort order
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
    the zip  tree, along with any relevant distances. Intervals representing the
    children of the snarl or  chain are found and added to the stack.
    This repeats until the stack is empty.

    */

    //Start by initializing the state
    //The forest state keeps track of the sort order of seeds,
    //the intervals that need to be sorted,
    //and which intervals are open and incomplete. 
    forest_growing_state_t forest_state(seeds, distance_index, distance_limit);

    //Start with root as the interval over seed_sort_order containing everything
    interval_state_t first_interval (0, seeds.size(), false, ZipCode::EMPTY, 0);

    //Sort and get the intervals of the connected components
    sort_one_interval(forest_state, first_interval);
    get_next_intervals(forest_state, first_interval, forest_state.intervals_to_process);


    while (!forest_state.intervals_to_process.empty()) {
#ifdef DEBUG_ZIP_CODE_TREE
    print_self(&seeds);
#endif
        // For each unprocessed interval, process it
        // First, check if anything needs to be closed,
        //   which will happen if the interval's depth is 
        //   greater than or equal to that of an open interval.
        // Get the intervals of this interval's children and add them in
        //   reverse order to the stack intervals_to_process
        // Open the current interval's snarl/chain


        //Get the interval
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
                //If the current interval is not a child of the open interval
                //close the last thing in open_intervals
                //There will be an interval for every ancestor in the snarl tree
                //so this can just check depth

#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "\tclose something at depth " << forest_state.open_intervals.size()-1 << endl;
#endif

                size_t depth = forest_state.open_intervals.size()-1;

                //The ancestor interval to close and its last seed
                const interval_state_t& ancestor_interval = forest_state.open_intervals.back();
                const Seed& last_seed = seeds.at(forest_state.seed_sort_order[ancestor_interval.interval_end-1]);

                if (ancestor_interval.code_type == ZipCode::CHAIN ||
                    ancestor_interval.code_type == ZipCode::NODE ||
                    ancestor_interval.code_type == ZipCode::ROOT_CHAIN ||
                    ancestor_interval.code_type == ZipCode::ROOT_NODE) {
                    //Close a chain

                    close_chain(forest_state, depth, 
                                last_seed, ancestor_interval.is_reversed); 
                } else {
#ifdef DEBUG_ZIP_CODE_TREE
    assert(ancestor_interval.code_type == ZipCode::REGULAR_SNARL ||
            ancestor_interval.code_type == ZipCode::IRREGULAR_SNARL ||
            ancestor_interval.code_type == ZipCode::CYCLIC_SNARL ||
            ancestor_interval.code_type == ZipCode::ROOT_SNARL);
#endif
                    //Close a snarl
                    close_snarl(forest_state, depth, last_seed, ancestor_interval.is_reversed); 
                }

                //Clear list of children of snarl tree structure at this level
                forest_state.sibling_indices_at_depth[depth].clear();

                //Take out this ancestor
                forest_state.open_intervals.pop_back();
            } else {
                //If the current interval is contained in this open interval,
                //then it is also contained in all other ancestors so break
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
            //Sort the current interval and get the intervals for its children
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
            // If there is nothing open, this starts a new connected component
            // Just open it

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
    //If we're starting a new tree then the last one must be valid
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
                // Open the root snarl (never cyclic, right?)
                open_snarl(forest_state, 0, false);
            } else if (current_interval.code_type == ZipCode::NODE) {
                //For a root node, just add it as a chain with all the seeds

                trees[forest_state.active_tree_index].zip_code_tree.emplace_back(ZipCodeTree::CHAIN_START, 
                                                                                 std::numeric_limits<size_t>::max());

                //Remember the start of the chain
                forest_state.sibling_indices_at_depth[0].push_back({ZipCodeTree::CHAIN_START, 0});
                forest_state.sibling_indices_at_depth[0].back().chain_component = 0;

                //If this is a node, then the interval contains everything in it
                //so add the seeds and close the chain here
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

                //Remember the start of the chain
                forest_state.sibling_indices_at_depth[0].push_back({ZipCodeTree::CHAIN_START, 0});
                forest_state.sibling_indices_at_depth[0].back().chain_component = 0;

            }                       
        } else if (forest_state.open_intervals.back().code_type == ZipCode::CHAIN || 
                   forest_state.open_intervals.back().code_type == ZipCode::ROOT_CHAIN ||
                   forest_state.open_intervals.back().code_type == ZipCode::ROOT_NODE) {
            // This is the child of a chain

            if (current_interval.code_type == ZipCode::NODE) {
                // If the type of this interval is NODE, then this is
                // a range of seeds that are on nodes on the chain,
                // not necessarily on the same node
                // Add each seed

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

                //Add the snarl to the chain
                add_child_to_chain(forest_state, current_interval.depth,
                                   forest_state.seed_sort_order[current_interval.interval_start], 
                                   current_interval.is_reversed, forest_state.open_intervals.back().is_reversed);
            }
            

        } else {
        //If there's an open ancestor that isn't a chain, it must be a snarl
#ifdef DEBUG_ZIP_CODE_TREE
    assert(forest_state.open_intervals.back().code_type == ZipCode::REGULAR_SNARL ||
            forest_state.open_intervals.back().code_type == ZipCode::IRREGULAR_SNARL ||
            forest_state.open_intervals.back().code_type == ZipCode::CYCLIC_SNARL ||
            forest_state.open_intervals.back().code_type == ZipCode::ROOT_SNARL);
#endif

            //Open the child chain
            open_chain(forest_state, forest_state.open_intervals.size(), 
                       forest_state.seed_sort_order[current_interval.interval_start], current_interval.is_reversed);
            
        }

        if (current_interval.code_type != ZipCode::NODE) {
            // Add to open_intervals
            forest_state.open_intervals.emplace_back(std::move(current_interval));
        }
    }
    //Finished adding all intervals
    

    //Now close anything that remained open
    while (!forest_state.open_intervals.empty()) {
        interval_state_t& ancestor_interval = forest_state.open_intervals.back();
        const Seed& last_seed = seeds.at(forest_state.seed_sort_order[ancestor_interval.interval_end-1]);

        if (ancestor_interval.code_type == ZipCode::CHAIN ||
            ancestor_interval.code_type == ZipCode::ROOT_CHAIN ||
            ancestor_interval.code_type == ZipCode::ROOT_NODE) {
            //Close a chain

            close_chain(forest_state, forest_state.open_intervals.size()-1, 
                        last_seed, ancestor_interval.is_reversed); 
        } else {
#ifdef DEBUG_ZIP_CODE_TREE
    assert(ancestor_interval.code_type == ZipCode::REGULAR_SNARL ||
            ancestor_interval.code_type == ZipCode::IRREGULAR_SNARL ||
            ancestor_interval.code_type == ZipCode::CYCLIC_SNARL ||
            ancestor_interval.code_type == ZipCode::ROOT_SNARL);
#endif
            //Close a snarl
            close_snarl(forest_state, forest_state.open_intervals.size()-1, 
                        last_seed, ancestor_interval.is_reversed); 
        }

        forest_state.open_intervals.pop_back();
    }

    if (trees[forest_state.active_tree_index].zip_code_tree.size() == 0) {
        trees.erase(trees.begin() + forest_state.active_tree_index);
    }
#ifdef DEBUG_ZIP_CODE_TREE
    print_self(&seeds);
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
    case vg::ZipCodeTree::NODE_COUNT:
        return "NODE_COUNT";
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
    case vg::ZipCodeTree::distance_iterator::S_FIND_CHAIN:
        return "S_FIND_CHAIN";
    case vg::ZipCodeTree::distance_iterator::S_SCAN_DAG_SNARL:
        return "S_SCAN_DAG_SNARL";
    case vg::ZipCodeTree::distance_iterator::S_SCAN_CYCLIC_SNARL:
        return "S_SCAN_CYCLIC_SNARL";
    case vg::ZipCodeTree::distance_iterator::S_SKIP_CHAIN:
        return "S_SKIP_CHAIN";
    case vg::ZipCodeTree::distance_iterator::S_SKIP_SNARL:
        return "S_SKIP_SNARL";
    case vg::ZipCodeTree::distance_iterator::S_FIND_DIST_MATRIX:
        return "S_FIND_DIST_MATRIX";
    default:
        throw std::runtime_error("Unimplemented zip code tree reverse iterator state");
    }
}



}


