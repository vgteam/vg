//#define DEBUG_ZIP_CODE_TREE
//#define PRINT_NON_DAG_SNARLS

#include "zip_code_tree.hpp"

#include "crash.hpp"

//#define debug_parse

using namespace std;
namespace vg {

void ZipCodeForest::fill_in_forest(vector<Seed>& all_seeds, const SnarlDistanceIndex& distance_index,
                               size_t distance_limit) {
    if (all_seeds.size() == 0) {
        return;
    }
    seeds = &all_seeds;

    /*
    Make a ZipCodeForest
    Takes a vector of seeds and fills in the forest

    Forest making is done by first sorting the seeds along chains/snarls
    Then, adding each seed, snarl/chain boundary, and distance to zip_code_tree
    A new tree is added to the forest for each connected component, and for any
    slice of a chain that is farther than the given distance_limit from anything
    on either side
    */

    //////////////////// Sort the seeds


    //Sort the seeds roughly linearly along top-level chains
    vector<size_t> seed_indices = sort_seeds_by_zipcode(distance_index);

#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "Sorted positions:" << endl;
    for (const size_t& i : seed_indices) {
        cerr << seeds->at(i).pos << endl;
    }
#endif

    //seed_indices is now sorted roughly along snarls and chains


    ///////////////////// Build the tree
    forest_growing_state_t forest_state;
    forest_state.active_zip_tree = std::numeric_limits<size_t>::max();

    /* The tree will hold all seeds and the bounds of snarls and chains
       For each chain, there must be a distance between each element of the chain (seeds and snarls)
       For each snarl, each element (chain or boundary) is preceded by the distances to everything
         before it in the snarl.
    */

    for (size_t i = 0 ; i < seed_indices.size() ; i++) {
#ifdef DEBUG_ZIP_CODE_TREE
        cerr << "At " << i << "st/nd/th seed: " << seeds->at(seed_indices[i]).pos << endl;
#endif

        //1. First, find the lowest common ancestor with the previous seed.
        //2. To finish the ancestors of the previous seed that are different from this one,
        //   walk up the snarl tree from the previous max depth and mark the end of the ancestor,
        //   adding distances for snarl ends 
        //3. To start anything for this seed, start from the first ancestor that is different
        //   and walk down the snarl tree, adding distances for each ancestor

        Seed& current_seed = seeds->at(seed_indices[i]);

        size_t current_max_depth = current_seed.zipcode_decoder->max_depth();
        //Make sure forest_state.sibling_indices_at_depth has enough spaces for this zipcode
        while (forest_state.sibling_indices_at_depth.size() < current_max_depth+1) {
            forest_state.sibling_indices_at_depth.emplace_back();
        }

        //Get the previous seed (if this isn't the first one)
        Seed& previous_seed = i == 0 ? current_seed : seeds->at(seed_indices[i-1]);
        //And the previous max depth
        size_t previous_max_depth = i == 0 ? 0 : previous_seed.zipcode_decoder->max_depth();

        //Remember the orientation for the seeds at the current depth
        //We start the first traversal (2) from previous_max_depth
        //The second traversal (3) starts from first_different_ancestor_depth 
        //This one is for the first traversal, so it will be for previous_max_depth
        bool previous_is_reversed = false;
        //This is for the second traversal, find it when finding first_different_ancestor_depth
        bool current_is_reversed = false;


        //Find the depth at which the two seeds are on different snarl tree nodes
        size_t first_different_ancestor_depth = 0;
        bool same_node = false;
        size_t max_depth = std::min(current_max_depth, previous_max_depth);
        size_t max_depth_checked = max_depth;

        for (size_t depth = 0 ; depth <= max_depth ; depth++) {
            first_different_ancestor_depth = depth;
            
            if (ZipCodeTree::seed_is_reversed_at_depth(current_seed, depth, distance_index)) {

                current_is_reversed = !current_is_reversed;

#ifdef DEBUG_ZIP_CODE_TREE
                cerr << "\tcurrent is reversed at depth " << depth << endl;
#endif
            }
            if (i != 0 && ZipCodeTree::seed_is_reversed_at_depth(previous_seed, depth, distance_index)) {

                previous_is_reversed = !previous_is_reversed;

#ifdef DEBUG_ZIP_CODE_TREE
                cerr << "\tprevious is reversed at depth " << depth << endl;
#endif
            }
            if (!ZipCodeDecoder::is_equal(*current_seed.zipcode_decoder, 
                        *previous_seed.zipcode_decoder, depth)) {
                max_depth_checked = depth;
                break;
            } else if (depth == max_depth) {
                same_node = true;
            }
        }
        if (previous_max_depth > max_depth_checked) {
            //We might need to update previous_is_reversed
            for (size_t depth = max_depth_checked+1 ; depth <= previous_max_depth ; depth++) {
                
                if (ZipCodeTree::seed_is_reversed_at_depth(previous_seed, depth, distance_index)) {
                    previous_is_reversed = !previous_is_reversed;

#ifdef DEBUG_ZIP_CODE_TREE
                    cerr << "\tprevious is reversed at depth " << depth << endl;
#endif
                }
            }
        }
        if (i == 0) { 
            same_node = false;
        }
#ifdef DEBUG_ZIP_CODE_TREE
        cerr << "\tthe depth of the first ancestor different than the previous seed is " << first_different_ancestor_depth << endl;
        cerr << "\tWalk up the snarl tree from depth " << previous_max_depth << " and close any snarl/chains" << endl;
#endif

        //Now, close anything that ended at the previous seed, starting from the leaf of the previous seed
        //If there was no previous seed, then the loop is never entered
        for (int depth = previous_max_depth ; !same_node && i!=0 && depth >= first_different_ancestor_depth && depth >= 0 ; depth--) {
            ZipCode::code_type_t previous_type = previous_seed.zipcode_decoder->get_code_type(depth);
            if (previous_type == ZipCode::CHAIN || previous_type == ZipCode::ROOT_CHAIN || previous_type == ZipCode::ROOT_NODE) {

                close_chain(forest_state, distance_index, distance_limit, depth, 
                            previous_seed, previous_is_reversed );

            } else if (previous_type == ZipCode::REGULAR_SNARL || previous_type == ZipCode::IRREGULAR_SNARL) { 

                close_snarl(forest_state, distance_index, depth, previous_seed, previous_is_reversed);
     
            }
            //Update previous_is_reversed to the one before this
            if (ZipCodeTree::seed_is_reversed_at_depth(previous_seed, depth, distance_index)) {
                previous_is_reversed = !previous_is_reversed;
            }

            //Clear the list of children of the thing at this level
            forest_state.sibling_indices_at_depth[depth].clear();
        }
#ifdef DEBUG_ZIP_CODE_TREE
        cerr << "\tWalk down the snarl tree from depth " << first_different_ancestor_depth << " to " << current_max_depth  << " and open any snarl/chains" << endl;
#endif

        //Now go through everything that started a new snarl tree node going down the snarl tree
        //For each new snarl or seed in a chain, add the distance to the thing preceding it in the chain
        //For each new chain in a snarl, add the distance to everything preceding it in the snarl
        //If this is the same node as the previous, then first_different_ancestor_depth is the depth 
        //of the node
        for (size_t depth = first_different_ancestor_depth ; depth <= current_max_depth ; depth++) {
            ZipCode::code_type_t current_type = current_seed.zipcode_decoder->get_code_type(depth);

            if (current_type == ZipCode::NODE || current_type == ZipCode::REGULAR_SNARL || current_type == ZipCode::IRREGULAR_SNARL
                || current_type == ZipCode::ROOT_NODE) {

                if (current_type == ZipCode::ROOT_NODE && forest_state.sibling_indices_at_depth[depth].empty()) {
                    //If this is a root-level node and the first time we've seen it,
                    //then open the node
#ifdef DEBUG_ZIP_CODE_TREE
                    cerr << "Add new root node as new tree" << endl;
#endif

                    //First, add this as a new connected component
                    trees.emplace_back(seeds);
                    forest_state.active_zip_tree = 0;

                    //Start the new tree
                    trees[forest_state.active_zip_tree].zip_code_tree.push_back({ZipCodeTree::CHAIN_START, std::numeric_limits<size_t>::max(), false});
                    forest_state.sibling_indices_at_depth[depth].push_back({ZipCodeTree::CHAIN_START, 0});
                }

                //Add the seed to its chain
                add_child_to_chain(forest_state, distance_index, distance_limit, depth, seed_indices[i], current_seed, current_is_reversed );
            } else if (current_type == ZipCode::ROOT_SNARL) {
                //If this is a root snarl, then just add the start of the snarl
                if (forest_state.sibling_indices_at_depth[depth].size() == 0) {
                    //IF this is the start of a new root snarl
#ifdef DEBUG_ZIP_CODE_TREE
                    cerr << "\t\tOpen new root snarl at depth " << depth << endl;
#endif

                    //Add a new subtree for the connected component
                    trees.emplace_back(seeds);
                    forest_state.active_zip_tree = trees.size()-1;

                    //Now record the start of this snarl
                    open_snarl(forest_state, 0);

                }
            } else {

                //Otherwise, this is a chain or root chain
                //If it is a chain, then it is the child of a snarl, so we need to find distances
                //to everything preceding it in the snarl
                assert(current_type == ZipCode::CHAIN || current_type == ZipCode::ROOT_CHAIN);

                //If this is the first time seeing the chain, then open it
                if (forest_state.sibling_indices_at_depth[depth].size() == 0) {
                    open_chain(forest_state, distance_index, distance_limit, depth, current_seed, current_is_reversed);
                }

                if (depth == current_max_depth) {
                    //If this is a trivial chain, then also add the seed and the distance to the 
                    //thing before it
                    add_child_to_chain(forest_state, distance_index, distance_limit, depth, seed_indices[i], current_seed, current_is_reversed);                    
                }
            }
            
            //Finished with this depth, so update current_is_reversed to be for the next ancestor
            if (depth < current_max_depth && ZipCodeTree::seed_is_reversed_at_depth(current_seed, depth+1, distance_index)) {
                current_is_reversed = !current_is_reversed;
            }
        }


    }
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "Close any snarls or chains that remained open" << endl;
#endif

    // Now close anything that remained open
    const Seed& last_seed = seeds->at(seed_indices.back());
    size_t last_max_depth = last_seed.zipcode_decoder->max_depth();

    //Find out if this seed is reversed at the leaf of the snarl tree (the node)
    bool last_is_reversed = false;
    for (size_t depth = 0 ; depth <= last_max_depth ; depth++) {
        if (ZipCodeTree::seed_is_reversed_at_depth(last_seed, depth, distance_index)) {
            last_is_reversed = !last_is_reversed;
        }
    }
    for (int depth = last_max_depth ; depth >= 0 ; depth--) {
        if (forest_state.sibling_indices_at_depth[depth].size() > 0) {
            ZipCode::code_type_t last_type = last_seed.zipcode_decoder->get_code_type(depth);
            if (last_type == ZipCode::CHAIN || last_type == ZipCode::ROOT_CHAIN || last_type == ZipCode::ROOT_NODE) {
                close_chain(forest_state, distance_index, distance_limit, depth, 
                            last_seed, last_is_reversed );

            } else if (last_type == ZipCode::REGULAR_SNARL || last_type == ZipCode::IRREGULAR_SNARL 
                        || last_type == ZipCode::ROOT_SNARL) { 

                close_snarl(forest_state, distance_index, depth, last_seed, last_is_reversed);

            }
        }
        //Update last_is_reversed to the one before this
        if (depth > 0 && ZipCodeTree::seed_is_reversed_at_depth(last_seed, depth-1, distance_index)) {
            last_is_reversed = !last_is_reversed;
        }
    }
#ifdef DEBUG_ZIP_CODE_TREE
    assert(forest_state.open_chains.empty());
#endif

}

void ZipCodeForest::open_chain(forest_growing_state_t& forest_state, const SnarlDistanceIndex& distance_index,
                       const size_t& distance_limit, const size_t& depth, Seed& current_seed, bool current_is_reversed) {
    //If this is the start of a new chain
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "\t\tOpen new chain at depth " << depth << endl;
#endif

    size_t current_max_depth = current_seed.zipcode_decoder->max_depth();

    if (depth == 0) {
        //If this is the start of a new top-level chain, make a new tree, which will be the new active tree
#ifdef DEBUG_ZIP_CODE_TREE
        cerr << "Add a new tree" << endl;
#endif
        trees.emplace_back(seeds);
        forest_state.active_zip_tree = trees.size()-1;
    } else  {
        //If this is the start of a non-root chain, then it is the child of a snarl and 
        //we need to find the distances to the previous things in the snarl
        //The distances will be filled in when the chain is closed, since parts of the
        //chain may be removed, and the distance to the start of the chain may change
        for (size_t i = 0 ; i < forest_state.sibling_indices_at_depth[depth-1].size() ; i++) {
            trees[forest_state.active_zip_tree].zip_code_tree.push_back({ZipCodeTree::EDGE, 
                  std::numeric_limits<size_t>::max(),
                  false});
        }

    }

    //Now record the start of this chain
    trees[forest_state.active_zip_tree].zip_code_tree.push_back({ZipCodeTree::CHAIN_START, std::numeric_limits<size_t>::max(), false});

    //Remember the start of the chain, with the prefix sum value
    forest_state.sibling_indices_at_depth[depth].push_back({ZipCodeTree::CHAIN_START, 0});

    //And, if it is the child of a snarl, then remember the chain as a child of the snarl
    if (depth != 0) {
        forest_state.sibling_indices_at_depth[depth-1].push_back({ZipCodeTree::CHAIN_START,
                                                     trees[forest_state.active_zip_tree].zip_code_tree.size()-1});

        //The distances in the snarl include the distances from the first/last children in the
        //chain to the ends of the chains
        //
        //Remember the distance to the start of this child in the chain
        if (depth == current_max_depth) {
            //If this is really a node, then get the distance to the start of the node
            forest_state.sibling_indices_at_depth[depth-1].back().distances.first =
                current_is_reversed != is_rev(current_seed.pos)
                    ? current_seed.zipcode_decoder->get_length(depth) - offset(current_seed.pos)
                    : offset(current_seed.pos);
        } else {
            //Otherwise, this is really a chain, so get the prefix sum in the chain

            forest_state.sibling_indices_at_depth[depth-1].back().distances.first = current_is_reversed 
                ? SnarlDistanceIndex::minus(current_seed.zipcode_decoder->get_length(depth) ,
                                            SnarlDistanceIndex::sum(
                                                current_seed.zipcode_decoder->get_offset_in_chain(depth+1),
                                                current_seed.zipcode_decoder->get_length(depth+1))) 
                : current_seed.zipcode_decoder->get_offset_in_chain(depth+1);

            if (depth+1 == current_max_depth) {
                //If this is a node, then add the offset of the position in the node
                bool child_is_reversed = ZipCodeTree::seed_is_reversed_at_depth(current_seed, depth+1, distance_index) 
                    ? !current_is_reversed : current_is_reversed;
                forest_state.sibling_indices_at_depth[depth-1].back().distances.first = 
                    SnarlDistanceIndex::sum(forest_state.sibling_indices_at_depth[depth-1].back().distances.first, 
                      child_is_reversed != is_rev(current_seed.pos)
                          ? current_seed.zipcode_decoder->get_length(depth+1) - offset(current_seed.pos)
                          : offset(current_seed.pos));
            }
        }

        //Remember the opening of this chain, and if its first child was far enough from the start to 
        //start a new subtree
        forest_state.open_chains.emplace_back(trees[forest_state.active_zip_tree].zip_code_tree.size()-1, 
                                              depth == 0 ? false : forest_state.sibling_indices_at_depth[depth-1].back().distances.first 
                                                                   > distance_limit);
    }
}

void ZipCodeForest::close_chain(forest_growing_state_t& forest_state, const SnarlDistanceIndex& distance_index,
                       const size_t& distance_limit, const size_t& depth, const Seed& last_seed, bool last_is_reversed) {

#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "\t\tclose a chain at depth " << depth << endl;
#endif
    if (trees[forest_state.active_zip_tree].zip_code_tree.back().type == ZipCodeTree::CHAIN_START) {
        //If the chain was empty.
        //This could happen if there was only a snarl in it and it got removed

        //Take out the CHAIN_START
        trees[forest_state.active_zip_tree].zip_code_tree.pop_back();

        //Forget about this chain in its parent snarl
        if (trees[forest_state.active_zip_tree].zip_code_tree.back().type == ZipCodeTree::EDGE) {
            forest_state.sibling_indices_at_depth[depth-1].pop_back();
        }

        //If the chain was part of a snarl, then take out the edges
        while (trees[forest_state.active_zip_tree].zip_code_tree.back().type == ZipCodeTree::EDGE) {
            trees[forest_state.active_zip_tree].zip_code_tree.pop_back();
        }

        //Forget about the chain
        forest_state.open_chains.pop_back();

    } else {
        //Add the end of the chain to the zip code tree
        trees[forest_state.active_zip_tree].zip_code_tree.push_back({ZipCodeTree::CHAIN_END, std::numeric_limits<size_t>::max(), false});


        // For chains in snarls, we want to know the distance from the last thing
        // in the chain to the end of the chain
        // If the distance is greater than the distance limit, we may make a new tree
        // for a slice of the chain.
        // If the chain remains in the snarl, we need to remember the distance to the end
        // of the chain to add to the relevant distances in the parent snarl.
        // These distances will be stored in forest_state.sibling_indices_at_depth

        if ( depth != 0 ) {
#ifdef DEBUG_ZIP_CODE_TREE
            assert(forest_state.sibling_indices_at_depth[depth-1].size() > 0);
            assert(forest_state.sibling_indices_at_depth[depth-1].back().type == ZipCodeTree::CHAIN_START);
#endif
            //Only add the distance for a non-root chain

            //If this is reversed, then the distance should be the distance to the start of 
            //the chain. Otherwise, the distance to the end
            //The value that got stored in forest_state.sibling_indices_at_depth was the prefix sum
            //traversing the chain according to its orientation in the tree, so either way
            //the distance is the length of the chain - the prefix sum
            size_t distance_to_chain_end = SnarlDistanceIndex::minus(last_seed.zipcode_decoder->get_length(depth),
                                          forest_state.sibling_indices_at_depth[depth].back().value);
            bool add_distances = true;
            if (distance_to_chain_end > distance_limit && forest_state.open_chains.back().second) {
                //If the distance to the end is greater than the distance limit, and there was something
                // in the chain with a large distance to the thing before it, then splice out a chain slice

                //Add a new tree
                trees.emplace_back(seeds);

                if (trees[forest_state.active_zip_tree].zip_code_tree.at(forest_state.open_chains.back().first).type 
                        == ZipCodeTree::CHAIN_START) {
                    //If we're copying the entire chain child of a snarl
#ifdef DEBUG_ZIP_CODE_TREE
                    cerr << "Copy the entire chain to a new subtree" << endl;
#endif

                    //Copy everything in the child chain into the new tree
                    trees.back().zip_code_tree.insert(trees.back().zip_code_tree.end(),
                        std::make_move_iterator(trees[forest_state.active_zip_tree].zip_code_tree.begin() 
                            + forest_state.open_chains.back().first),
                        std::make_move_iterator(trees[forest_state.active_zip_tree].zip_code_tree.end()));

                    //Remove the child chain from the active tree
                    trees[forest_state.active_zip_tree].zip_code_tree.erase(
                           trees[forest_state.active_zip_tree].zip_code_tree.begin() + forest_state.open_chains.back().first,
                           trees[forest_state.active_zip_tree].zip_code_tree.end());

                    //The chain no longer exists in the snarl, so forget that it exists
                    forest_state.sibling_indices_at_depth[depth-1].pop_back();

                    //And remove all the edges
                    while (trees[forest_state.active_zip_tree].zip_code_tree.back().type == ZipCodeTree::EDGE) {
                        trees[forest_state.active_zip_tree].zip_code_tree.pop_back();
                    }
#ifdef DEBUG_ZIP_COD
                    assert((trees[forest_state.active_zip_tree].zip_code_tree.back().type == ZipCodeTree::CHAIN_END ||
                            trees[forest_state.active_zip_tree].zip_code_tree.back().type == ZipCodeTree::SNARL_START));
#endif
                    // Since we took out the whole chain, we shouldn't add the distances later
                    add_distances = false;
                } else {
#ifdef DEBUG_ZIP_CODE_TREE
                    cerr << "Copy a slice from the middle of the chain to the end" << endl;
                    assert((trees[forest_state.active_zip_tree].zip_code_tree.at(forest_state.open_chains.back().first).type == ZipCodeTree::SEED || 
                            trees[forest_state.active_zip_tree].zip_code_tree.at(forest_state.open_chains.back().first).type == ZipCodeTree::SNARL_START));
#endif
                    //We're copying a slice of the chain from the middle to the end
                    //Start a new chain in the new subtree
                    trees.back().zip_code_tree.push_back({ZipCodeTree::CHAIN_START, 
                                                          std::numeric_limits<size_t>::max(), false});

                    //Copy everything in the slice into the new tree
                    trees.back().zip_code_tree.insert(trees.back().zip_code_tree.end(),
                        std::make_move_iterator(trees[forest_state.active_zip_tree].zip_code_tree.begin() + forest_state.open_chains.back().first),
                        std::make_move_iterator(trees[forest_state.active_zip_tree].zip_code_tree.end()));
                    //Erase the slice
                    trees[forest_state.active_zip_tree].zip_code_tree.erase(
                            trees[forest_state.active_zip_tree].zip_code_tree.begin() + forest_state.open_chains.back().first,
                            trees[forest_state.active_zip_tree].zip_code_tree.end());


                    //Take out the last edge
                    size_t last_edge = trees[forest_state.active_zip_tree].zip_code_tree.back().value;
                    trees[forest_state.active_zip_tree].zip_code_tree.pop_back();

                    //Close the chain in the original active tree
                    trees[forest_state.active_zip_tree].zip_code_tree.push_back({ZipCodeTree::CHAIN_END, 
                                                std::numeric_limits<size_t>::max(), false});

                    //Update the distance to the end of the chain to be the distance from the previous child 
                    distance_to_chain_end = SnarlDistanceIndex::sum(distance_to_chain_end, 
                                            SnarlDistanceIndex::sum(last_edge,
                                                                    last_seed.zipcode_decoder->get_length(depth+1)));
                }
            }
            if (add_distances) {
                // If this chain (or chain slice) remains in the snarl, then add the distances
                // in the snarl

                //remember the distance to the end to be used in snarl distances
                forest_state.sibling_indices_at_depth[depth-1].back().distances.second = distance_to_chain_end;

                add_snarl_distances(forest_state, distance_index, depth-1, last_seed, last_is_reversed, false);
            }
            //We've closed a chain, so take out the latest open chain
            forest_state.open_chains.pop_back();
        }
    }
}

void ZipCodeForest::add_child_to_chain(forest_growing_state_t& forest_state, const SnarlDistanceIndex& distance_index,
                       const size_t& distance_limit, const size_t& depth, const size_t& seed_index, Seed& current_seed, bool current_is_reversed) {
    //For these things, we need to remember the offset in the node/chain

    ZipCode::code_type_t current_type = current_seed.zipcode_decoder->get_code_type(depth);

    //Is this chain actually a node pretending to be a chain
    bool is_trivial_chain = current_type == ZipCode::CHAIN && depth == current_seed.zipcode_decoder->max_depth();

    //For a root node or trivial chain, the "chain" is actually just the node, so the depth 
    // of the chain we're working on is the same depth. Otherwise, the depth is depth-1
    size_t chain_depth = is_trivial_chain || current_type == ZipCode::ROOT_NODE ? depth : depth-1;

    ///////////////// Get the offset in the parent chain (or node)
    size_t current_offset;


    //First, get the prefix sum in the chain
    if (current_type == ZipCode::ROOT_NODE || is_trivial_chain) {
        //Which is 0 if this is just a node
        current_offset = 0;
    } else {
        //And the distance to the start or end of the chain if it's a node/snarl in a chain

        //If we're traversing this chain backwards, then the offset is the offset from the end
        bool chain_is_reversed = ZipCodeTree::seed_is_reversed_at_depth(current_seed, depth, distance_index) 
            ? !current_is_reversed : current_is_reversed;

        current_offset = chain_is_reversed 
                ? SnarlDistanceIndex::minus(current_seed.zipcode_decoder->get_length(chain_depth) ,
                                            SnarlDistanceIndex::sum(
                                                current_seed.zipcode_decoder->get_offset_in_chain(depth),
                                                current_seed.zipcode_decoder->get_length(depth))) 
                : current_seed.zipcode_decoder->get_offset_in_chain(depth);
    }

    if (depth == current_seed.zipcode_decoder->max_depth()) {
        //If this is a node, then add the offset of the seed in the node
        current_offset = SnarlDistanceIndex::sum(current_offset, 
            current_is_reversed != is_rev(current_seed.pos)
                ? current_seed.zipcode_decoder->get_length(depth) - offset(current_seed.pos)
                : offset(current_seed.pos));

    }

    /////////////////////// Get the offset of the previous thing in the parent chain/node
    size_t previous_offset = forest_state.sibling_indices_at_depth[chain_depth][0].value;


#ifdef DEBUG_ZIP_CODE_TREE
    assert(forest_state.sibling_indices_at_depth[chain_depth].size() == 1);
#endif

    ///////////////////// Record the distance from the previous thing in the chain/node
    //       Or add a new tree if the distance is too far
    if (chain_depth > 0 && forest_state.sibling_indices_at_depth[chain_depth][0].type == ZipCodeTree::CHAIN_START){
        //If this is the first thing in a non-root chain or node, remember the distance to the 
        //start of the chain/node.
        //This distance will be added to distances in the parent snarl
        forest_state.sibling_indices_at_depth[chain_depth-1][0].distances.first = current_offset;

        //Also update the last chain opened
        forest_state.open_chains.back().second = current_offset > distance_limit;


    } else if (forest_state.sibling_indices_at_depth[chain_depth][0].type != ZipCodeTree::CHAIN_START) {
        //for everything except the first thing in a node/chain, we need to add the edge

        size_t distance_between;
        if (previous_offset > current_offset) {
            //If the parent is a multicomponent chain, then they might be in different components
            //TODO: This won't catch all cases of different components in the chain
            distance_between = std::numeric_limits<size_t>::max();
        } else {
            distance_between = current_offset - previous_offset;
        }

        if (chain_depth == 0 && distance_between > distance_limit) {
            //The next thing in the zip tree will be the first seed (or snarl) in a top-level chain, 
            // so start a new tree
#ifdef DEBUG_ZIP_CODE_TREE
            cerr << "Start a new tree in the forest" << endl;
#endif
            //Add the end of the first chain
            trees[forest_state.active_zip_tree].zip_code_tree.push_back({ZipCodeTree::CHAIN_END, 
                                                                         std::numeric_limits<size_t>::max(), 
                                                                         false});

            //Add a new tree and make sure it is the new active tree
            trees.emplace_back(seeds);
            forest_state.active_zip_tree = trees.size()-1;

            //Add the start of the new chain
            trees[forest_state.active_zip_tree].zip_code_tree.push_back({ZipCodeTree::CHAIN_START, 
                                                                         std::numeric_limits<size_t>::max(), 
                                                                         false});

            //The first sibling in the chain is now the chain start, not the previous seed, so replace it
            forest_state.sibling_indices_at_depth[chain_depth].pop_back();
            forest_state.sibling_indices_at_depth[chain_depth].push_back({ZipCodeTree::CHAIN_START, 0}); 

        } else if (distance_between > distance_limit) { 
            //If this is too far from the previous thing in a nested chain

            if (forest_state.open_chains.back().second) {
#ifdef DEBUG_ZIP_CODE_TREE
                cerr << "\tMake a new slice of the chain at depth " << depth << endl;
#endif
                //If the current chain slice was also too far away from the thing before it
                // then copy the slice
                if (trees[forest_state.active_zip_tree].zip_code_tree.at(forest_state.open_chains.back().first).type 
                        == ZipCodeTree::CHAIN_START) {
                    //If the slice starts at the start of the chain and ends at the previous seed

                    //Copy everything in the slice to the end of a new tree
                    trees.emplace_back(seeds);
                    trees.back().zip_code_tree.insert(trees.back().zip_code_tree.end(),
                        std::make_move_iterator(trees[forest_state.active_zip_tree].zip_code_tree.begin() + forest_state.open_chains.back().first),
                        std::make_move_iterator(trees[forest_state.active_zip_tree].zip_code_tree.end()));

                    //Erase the slice from the active tree
                    trees[forest_state.active_zip_tree].zip_code_tree.erase(
                        trees[forest_state.active_zip_tree].zip_code_tree.begin() + forest_state.open_chains.back().first,
                        trees[forest_state.active_zip_tree].zip_code_tree.end());

                    //Add the end of the chain to the new slice
                    trees.back().zip_code_tree.push_back({ZipCodeTree::CHAIN_END, 
                                                          std::numeric_limits<size_t>::max(), 
                                                          false});

                    //Add back the start of the chain
                    trees[forest_state.active_zip_tree].zip_code_tree.push_back({ZipCodeTree::CHAIN_START, 
                                                                                 std::numeric_limits<size_t>::max(), 
                                                                                 false});

                    //Update the chain as a child of the snarl
#ifdef DEBUG_ZIP_CODE_TREE
                    assert(forest_state.sibling_indices_at_depth[chain_depth-1].back().type == ZipCodeTree::CHAIN_START);
                    //The value should be the index of the last seed, which is the first seed in the new tree
                    assert(forest_state.sibling_indices_at_depth[chain_depth-1].back().value 
                                == trees[forest_state.active_zip_tree].zip_code_tree.size()-1);
                    assert(forest_state.open_chains.back().second);

#endif
                    forest_state.sibling_indices_at_depth[chain_depth-1].back().distances.first = current_offset;

                    //Don't need to update open_chains, since the next slice will also start at the chain start and be able to make 
                    //a new thing

                } else {
#ifdef DEBUG_ZIP_CODE_TREE
                    assert((trees[forest_state.active_zip_tree].zip_code_tree.at(forest_state.open_chains.back().first).type 
                               == ZipCodeTree::SEED ||
                           trees[forest_state.active_zip_tree].zip_code_tree.at(forest_state.open_chains.back().first).type 
                                == ZipCodeTree::SNARL_START));
#endif
                    //If the slice starts and ends in the middle of the chain

                    //Copy everything in the slice to a new chain in a new tree
                    trees.emplace_back(seeds);
                    trees.back().zip_code_tree.push_back({ZipCodeTree::CHAIN_START, 
                                                          std::numeric_limits<size_t>::max(), 
                                                          false});
                    trees.back().zip_code_tree.insert(trees.back().zip_code_tree.end(),
                        std::make_move_iterator(trees[forest_state.active_zip_tree].zip_code_tree.begin() 
                                                 + forest_state.open_chains.back().first),
                        std::make_move_iterator(trees[forest_state.active_zip_tree].zip_code_tree.end()));

                    //Erase the slice from the active tree
                    trees[forest_state.active_zip_tree].zip_code_tree.erase(
                        trees[forest_state.active_zip_tree].zip_code_tree.begin() + forest_state.open_chains.back().first,
                        trees[forest_state.active_zip_tree].zip_code_tree.end());
                    //Add the end of the chain to the new slice
                    trees.back().zip_code_tree.push_back({ZipCodeTree::CHAIN_END, 
                                                          std::numeric_limits<size_t>::max(), 
                                                          false});
                    //The original tree gets an edge with infinite length, since it will be bigger than the distance limit anyway
#ifdef DEBUG_ZIP_CODE_TREE
                    assert(trees[forest_state.active_zip_tree].zip_code_tree.back().type == ZipCodeTree::EDGE);
#endif
                    trees[forest_state.active_zip_tree].zip_code_tree.pop_back();
                    trees[forest_state.active_zip_tree].zip_code_tree.push_back({ZipCodeTree::EDGE, 
                                                                                 std::numeric_limits<size_t>::max(),
                                                                                 false});

                    //Remember the next seed or snarl that gets added as the start of a new chain slice
                    forest_state.open_chains.pop_back();
                    forest_state.open_chains.emplace_back(trees[forest_state.active_zip_tree].zip_code_tree.size(), true);
                }
            } else {
#ifdef DEBUG_ZIP_CODE_TREE
                cerr << "The slice didn't get copied but maybe start a new slice here" << endl;
#endif
                //If the slice doesn't get copied because it is still connected at the front,
                //add the edge anyway

                trees[forest_state.active_zip_tree].zip_code_tree.push_back({ZipCodeTree::EDGE, distance_between, false});

                //Remember the next seed or snarl that gets added as the start of a new chain slice
                forest_state.open_chains.pop_back();
                forest_state.open_chains.emplace_back(trees[forest_state.active_zip_tree].zip_code_tree.size(), true);
            }

        } else {
            //If we didn't start a new tree, then remember the edge
            trees[forest_state.active_zip_tree].zip_code_tree.push_back({ZipCodeTree::EDGE, distance_between, false});
        }
    }

    /////////////////////////////Record this thing in the chain
    if (current_type == ZipCode::NODE || current_type == ZipCode::ROOT_NODE || is_trivial_chain) {
#ifdef DEBUG_ZIP_CODE_TREE
        cerr << "\t\tContinue node/chain with seed " << current_seed.pos << " at depth " << depth << endl;
#endif
        //If this was a node, just remember the seed
        trees[forest_state.active_zip_tree].zip_code_tree.push_back({ZipCodeTree::SEED, seed_index, current_is_reversed != is_rev(current_seed.pos)});
    } else {

        open_snarl(forest_state, depth); 

        //For finding the distance to the next thing in the chain, the offset
        //stored should be the offset of the end bound of the snarl, so add the 
        //length of the snarl
        current_offset = SnarlDistanceIndex::sum(current_offset,
            current_seed.zipcode_decoder->get_length(depth));

    }

    //Remember this thing for the next sibling in the chain
    forest_state.sibling_indices_at_depth[chain_depth].pop_back();
    forest_state.sibling_indices_at_depth[chain_depth].push_back({(
         current_type == ZipCode::NODE || current_type == ZipCode::ROOT_NODE) ? ZipCodeTree::SEED 
                                                                              : ZipCodeTree::SNARL_START,
         current_offset}); 
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "Add sibling with type " << current_type << endl;
#endif

}

void ZipCodeForest::open_snarl(forest_growing_state_t& forest_state, const size_t& depth) {
#ifdef DEBUG_ZIP_CODE_TREE
        cerr << "\t\tOpen new snarl at depth " << depth << endl;
#endif
    //If this was a snarl, record the start of the snarl
    trees[forest_state.active_zip_tree].zip_code_tree.push_back({ZipCodeTree::SNARL_START, std::numeric_limits<size_t>::max(), false});
    
    if (depth != 0) {
        //Remember the start of the snarl to find distances later
        //Don't do this for a root snarl because technically there is no start node so there are no distances to it
        forest_state.sibling_indices_at_depth[depth].push_back({ZipCodeTree::SNARL_START, std::numeric_limits<size_t>::max()});
    }
}

void ZipCodeForest::close_snarl(forest_growing_state_t& forest_state, const SnarlDistanceIndex& distance_index, 
                                const size_t& depth, const Seed& last_seed, bool last_is_reversed) {
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "\t\tclose a snarl at depth " << depth << endl;
#endif

    if (depth == 0) {
        //If this is a root snarl, then we don't need distances so just close it 
        trees[forest_state.active_zip_tree].zip_code_tree.push_back({ZipCodeTree::SNARL_END, std::numeric_limits<size_t>::max(), false});

    } else if (forest_state.sibling_indices_at_depth[depth].size() == 1) {
        //Since some of the children of the snarl may have been removed to separate subtrees,
        //the snarl may actually be empty now
        //If there is only one "child" (the snarl start), then the snarl is actually empty, so delete it

#ifdef DEBUG_ZIP_CODE_TREE
        cerr << "\t\t\tThe snarl is actually empty so remove it" << endl;
#endif
        //Take out the edges
        while (trees[forest_state.active_zip_tree].zip_code_tree.back().type == ZipCodeTree::EDGE) {
            trees[forest_state.active_zip_tree].zip_code_tree.pop_back();
        }
#ifdef DEBUG_ZIP_CODE_TREE
        assert(trees[forest_state.active_zip_tree].zip_code_tree.back().type == ZipCodeTree::SNARL_START);
#endif        //Pop the snarl start out
        trees[forest_state.active_zip_tree].zip_code_tree.pop_back();

        //If this was the first thing in the chain, then we're done. Otherwise, there was an edge to remove
        if (trees[forest_state.active_zip_tree].zip_code_tree.back().type == ZipCodeTree::EDGE) {
            //If the snarl was in the middle of a chain, then we need to take out the edge and update
            //the previous thing in the chain with its prefix sum

            //This was the distance from the last thing to the start of this snarl
            size_t previous_edge = trees[forest_state.active_zip_tree].zip_code_tree.back().value;
            trees[forest_state.active_zip_tree].zip_code_tree.pop_back();

            //This is the distance from the start of the chain to the end of the snarl
            size_t snarl_prefix_sum = forest_state.sibling_indices_at_depth[depth-1].back().value;
            forest_state.sibling_indices_at_depth[depth-1].pop_back();

            //Snarl prefix sum is now the distance from the start of the chain to the start of the snarl
            snarl_prefix_sum = SnarlDistanceIndex::minus(snarl_prefix_sum, last_seed.zipcode_decoder->get_length(depth));

            //Now update forest_state.sibling_indices_at_depth to be the previous thing in the chain
            forest_state.sibling_indices_at_depth[depth-1].push_back({
                trees[forest_state.active_zip_tree].zip_code_tree.back().type == ZipCodeTree::SEED ? ZipCodeTree::SEED 
                                                                                                   : ZipCodeTree::SNARL_START,
                SnarlDistanceIndex::minus(snarl_prefix_sum, previous_edge)});
#ifdef DEBUG_ZIP_CODE_TREE
            assert(forest_state.sibling_indices_at_depth[depth-1].back().value >= 0);
#endif
        } else {
            //If this was the first thing in the chain, update the previous sibling in the chain to be the start of the chain
#ifdef DEBUG_ZIP_CODE_TREE
            assert(trees[forest_state.active_zip_tree].zip_code_tree.back().type == ZipCodeTree::CHAIN_START);
#endif
            forest_state.sibling_indices_at_depth[depth-1].pop_back();
            forest_state.sibling_indices_at_depth[depth-1].push_back({ ZipCodeTree::CHAIN_START, 0});
        }
    } else {

        //If this is the end of the snarl that still has children, then we need to save the distances to 
        //all previous children of the snarl
        trees[forest_state.active_zip_tree].zip_code_tree.resize(trees[forest_state.active_zip_tree].zip_code_tree.size() 
                                                                    + forest_state.sibling_indices_at_depth[depth].size());
       
        add_snarl_distances(forest_state, distance_index, depth, last_seed, last_is_reversed, true);

        //Note the count of children and the end of the snarl
        trees[forest_state.active_zip_tree].zip_code_tree.push_back({ZipCodeTree::NODE_COUNT, 
                                                                     forest_state.sibling_indices_at_depth[depth].size()-1, 
                                                                     false});
        trees[forest_state.active_zip_tree].zip_code_tree.push_back({ZipCodeTree::SNARL_END, 
                                                                     std::numeric_limits<size_t>::max(), 
                                                                     false});
     }
}

void ZipCodeForest::add_snarl_distances(forest_growing_state_t& forest_state, const SnarlDistanceIndex& distance_index, 
                                const size_t& depth, const Seed& seed, bool is_reversed, bool to_snarl_end) {

    // This adds distances from everything in the snarl to the last thing in the snarl, which is either the snarl end
    // or a chain child of the snarl

    //Distances from this child to add 
    size_t distance_to_chain_end = to_snarl_end ? 0 : forest_state.sibling_indices_at_depth[depth].back().distances.second;
    size_t distance_to_chain_start = to_snarl_end ? 0 : forest_state.sibling_indices_at_depth[depth].back().distances.first;
    
    // This is the index of the thing in the snarl right before the distances start. Used to figure out
    // where to put the distances
    size_t last_child_index = to_snarl_end ? trees[forest_state.active_zip_tree].zip_code_tree.size()
                                           : forest_state.sibling_indices_at_depth[depth].back().value;
    
    //Now add the distances from the start of the chain to everything before it in the snarl
    
    
    //If the parent snarl is reversed
    bool snarl_is_reversed = to_snarl_end ? is_reversed 
                                          : (ZipCodeTree::seed_is_reversed_at_depth(seed, depth+1, distance_index) 
                                                ? !is_reversed : is_reversed);
    
    
    // If this is to the end bound, get the distance to all siblings. If it is to the last child, don't get
    // the distance to itself
    size_t sibling_count = to_snarl_end ? forest_state.sibling_indices_at_depth[depth].size() 
                                        : forest_state.sibling_indices_at_depth[depth].size()-1;
    for ( size_t sibling_i = 0 ; sibling_i < sibling_count ; sibling_i++) {
        const auto& sibling = forest_state.sibling_indices_at_depth[depth][sibling_i];

        if (sibling.type == ZipCodeTree::SNARL_START) {
            //Get the distance to the start (or end if it's reversed) of the snarl

            //If we're getting the distance to the end of the snarl, then this is the length of the snarl
            // otherwise, it is the distance from the seed to the start (or end) of the snarl
            size_t snarl_distance = to_snarl_end ? seed.zipcode_decoder->get_length(depth)
                                                 : SnarlDistanceIndex::sum (distance_to_chain_start,
                                                    snarl_is_reversed
                                                        ? seed.zipcode_decoder->get_distance_to_snarl_end(depth+1)
                                                        : seed.zipcode_decoder->get_distance_to_snarl_start(depth+1));

            //Add the edge
            trees[forest_state.active_zip_tree].zip_code_tree.at(last_child_index - 1 - sibling_i) = 
                {ZipCodeTree::EDGE, snarl_distance, false};

        } else {
            //Otherwise, the previous thing was another child of the snarl
            //and we need to record the distance between these two
            //TODO: This can be improved for simple snarls
            size_t distance;
            if (seed.zipcode_decoder->get_code_type(depth) == ZipCode::REGULAR_SNARL) {
                //If this is the child of a regular snarl, then the distance between
                //any two chains is inf, and the distance to any bound is 0
                distance = to_snarl_end ? sibling.distances.second : std::numeric_limits<size_t>::max();
            } else {
                size_t seed_i = sibling.value+1;
                while (trees[forest_state.active_zip_tree].zip_code_tree[seed_i].type != ZipCodeTree::SEED) {
                    seed_i++;
                }
                auto& sibling_seed = seeds->at(trees[forest_state.active_zip_tree].zip_code_tree[seed_i].value);

                if (to_snarl_end) {
                    distance = SnarlDistanceIndex::sum( sibling.distances.second,
                                                        is_reversed ? sibling_seed.zipcode_decoder->get_distance_to_snarl_start(depth+1)
                                                                    : sibling_seed.zipcode_decoder->get_distance_to_snarl_end(depth+1));
                } else {
                    size_t rank2 = seed.zipcode_decoder->get_rank_in_snarl(depth+1);
                    size_t rank1 = sibling_seed.zipcode_decoder->get_rank_in_snarl(depth+1);
                    bool rev2 = is_reversed;
                    bool rev1 = ZipCodeTree::seed_is_reversed_at_depth(sibling_seed, depth+1, distance_index);

                    size_t distance_to_end_of_last_child = sibling.type == ZipCodeTree::SNARL_START ? 0
                                                     : sibling.distances.second;
                    //TODO: idk about this distance- I think the orientations need to change
                    //The bools for this are true if the distance is to/from the right side of the child
                    //We want the right side of 1 (which comes first in the dag ordering) to the left side of 2
                    //relative to the orientation of the snarl
                    net_handle_t snarl_handle = seed.zipcode_decoder->get_net_handle(depth, &distance_index);
                    distance = SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(
                        distance_index.distance_in_snarl(snarl_handle, rank1, !rev1, rank2, rev2),
                        distance_to_chain_start),
                        distance_to_end_of_last_child);
                }
            }
            trees[forest_state.active_zip_tree].zip_code_tree.at(last_child_index - 1 - sibling_i) = {ZipCodeTree::EDGE, distance, false};
        }
    
    }
}

std::pair<size_t, size_t> ZipCodeTree::dag_and_non_dag_snarl_count(vector<Seed>& seeds, const SnarlDistanceIndex& distance_index) const {
    size_t dag_count = 0;
    size_t non_dag_count = 0;

    /* Walk through everything in the zip code tree and at the first seed in each snarl, 
       check if it is a dag or not
    */

    //Keep track of the depth to check the zip codes
    size_t current_depth = 0;

    //When we encounter the start of a snarl, make a note of the depth. At the next seed,
    //check the snarls at the depths recorded
    vector<size_t> snarl_depths;

    for (size_t i = 0 ; i < zip_code_tree.size() ; i++ ) {
        const tree_item_t& current_item = zip_code_tree[i];
        if (current_item.type == ZipCodeTree::SNARL_START) {
            //For the start of a snarl, make a note of the depth to check the next seed
            snarl_depths.emplace_back(current_depth);

            //Increment the depth
            current_depth++;
        } else if (current_item.type == ZipCodeTree::CHAIN_START) {
            //For the start of a chain, increment the depth
            current_depth++;
        } else if (current_item.type == ZipCodeTree::CHAIN_END || current_item.type == ZipCodeTree::SNARL_END) {
            //For the end of a snarl or chain, decrement the depth
            current_depth--;
        } else if (current_item.type == ZipCodeTree::SEED) {
            //If this is a seed, check the snarls we've seen previously
            for (const size_t& snarl_depth : snarl_depths) {
                if (seeds[current_item.value].zipcode_decoder->get_code_type(snarl_depth) == ZipCode::REGULAR_SNARL) {
                    //If this is a regular snarl, then it must be a DAG too
                    dag_count++;
                } else {
                    //If this is an irregular snarl

                    //Check the snarl in the distance index
                    net_handle_t snarl_handle = seeds[current_item.value].zipcode_decoder->get_net_handle(snarl_depth, &distance_index);
#ifdef DEBUG_ZIP_CODE_TREE
                    assert(seeds[current_item.value].zipcode_decoder->get_code_type(snarl_depth) == ZipCode::IRREGULAR_SNARL ||
                           seeds[current_item.value].zipcode_decoder->get_code_type(snarl_depth) == ZipCode::ROOT_SNARL);
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

void ZipCodeTree::print_self() const {
    for (const tree_item_t item : zip_code_tree) {
        if (item.type == SEED) {
            cerr << seeds->at(item.value).pos;
            if (item.is_reversed) {
                cerr << "rev";
            }
        } else if (item.type == SNARL_START) {
            cerr << "(";
        } else if (item.type == SNARL_END) {
            cerr << ")";
        } else if (item.type == CHAIN_START) {
            cerr << "[";
        } else if (item.type == CHAIN_END) {
            cerr << "]";
        } else if (item.type == EDGE) {
            cerr << " " << item.value << " ";
        } else if (item.type == NODE_COUNT) {
            cerr << " " << item.value;
        } else {
            throw std::runtime_error("[zip tree]: Trying to print a zip tree item of the wrong type");
        }
    }
    cerr << endl;
}

void ZipCodeTree::validate_zip_tree(const SnarlDistanceIndex& distance_index, size_t distance_limit) const {

    //Make sure that everything is in a valid order
    size_t previous_seed_index = std::numeric_limits<size_t>::max();
    bool previous_is_valid = true;
    for (size_t i = 0 ; i < zip_code_tree.size() ; i++) {
        const tree_item_t& current_item = zip_code_tree[i];
        if (current_item.type == SEED) {
            bool current_is_valid = true;
            //Check if this is worth validating
            //TODO: For now, ignore anything with non-dag snarls, multicomponent or looping chains
            net_handle_t net = distance_index.get_node_net_handle(id(seeds->at(current_item.value).pos));
            while (!distance_index.is_root(net)) {
                if ((distance_index.is_snarl(net) && !distance_index.is_dag(net)) ||
                    distance_index.is_multicomponent_chain(net) || distance_index.is_looping_chain(net)) {
                    //If this is something that we haven't handled
                    current_is_valid = false;
                    cerr << "warning: validating a zip tree with a non-dag snarl, multicomponent chain, or looping chain" << endl;
                    break;
                }
                net = distance_index.get_parent(net);
            }
            if (previous_seed_index != std::numeric_limits<size_t>::max() &&
                current_is_valid && previous_is_valid) {
                assert(previous_seed_index < seeds->size());
                assert(current_item.value < seeds->size());
#ifdef DEBUG_ZIP_CODE_TREE
                cerr << "Comparing seeds " << seeds->at(previous_seed_index).pos << " and " << seeds->at(current_item.value).pos << endl;
#endif

                //Comparator returning previous_seed_index < current_item.value
                size_t depth = 0;

                //Keep track of the orientation of each seed
                //Everything should be sorted according to the orientation in the top-level structure,
                //so if things are traversed backwards, reverse the orientation
                bool a_is_reversed = false;
                bool b_is_reversed = false;
                while (depth < seeds->at(previous_seed_index).zipcode_decoder->max_depth() &&
                       depth < seeds->at(current_item.value).zipcode_decoder->max_depth() &&
                       ZipCodeDecoder::is_equal(*seeds->at(previous_seed_index).zipcode_decoder, *seeds->at(current_item.value).zipcode_decoder, depth)) {

                    //Remember the orientation
                    if (ZipCodeTree::seed_is_reversed_at_depth(seeds->at(previous_seed_index), depth, distance_index)) { 
                        a_is_reversed = !a_is_reversed;
                    }
                    if (ZipCodeTree::seed_is_reversed_at_depth(seeds->at(current_item.value), depth, distance_index)) {
                        b_is_reversed = !b_is_reversed;
                    }

                    depth++;
                }

                //Remember the orientation of the parent too
                size_t parent_of_a_is_reversed = a_is_reversed;

                //Check the orientations one last time
                if (ZipCodeTree::seed_is_reversed_at_depth(seeds->at(previous_seed_index), depth, distance_index)) { 
                    a_is_reversed = !a_is_reversed;
                }
                if (ZipCodeTree::seed_is_reversed_at_depth(seeds->at(current_item.value), depth, distance_index)) {
                    b_is_reversed = !b_is_reversed;
                }
                
#ifdef DEBUG_ZIP_CODE_TREE
                cerr << "\t different at depth " << depth << endl;
#endif
                //Either depth is the last thing in previous_seed_index or current_item.value, or they are different at this depth


                if ( ZipCodeDecoder::is_equal(*seeds->at(previous_seed_index).zipcode_decoder, *seeds->at(current_item.value).zipcode_decoder, depth)) {
#ifdef DEBUG_ZIP_CODE_TREE
                    cerr << "\tthey are on the same node" << endl;
#endif
                    //If they are equal, then they must be on the same node

                    size_t offset1 = is_rev(seeds->at(previous_seed_index).pos)
                                   ? seeds->at(previous_seed_index).zipcode_decoder->get_length(depth) - offset(seeds->at(previous_seed_index).pos) - 1
                                   : offset(seeds->at(previous_seed_index).pos);
                    size_t offset2 = is_rev(seeds->at(current_item.value).pos)
                                   ? seeds->at(current_item.value).zipcode_decoder->get_length(depth) - offset(seeds->at(current_item.value).pos) - 1
                                   : offset(seeds->at(current_item.value).pos);
                    if (!a_is_reversed) {
                        //If they are in previous_seed_index snarl or they are facing forward on a chain, then order by
                        //the offset in the node
                        assert( offset1 <= offset2);
                    } else {
                        //Otherwise, the node is facing backwards in the chain, so order backwards in node
                        assert( offset2 <= offset1);
                    }
                }  else if (depth == 0) {
#ifdef DEBUG_ZIP_CODE_TREE
                    cerr << "\tThey are on different connected components" << endl;
#endif
                    //If they are on different connected components, sort by connected component
                    assert( seeds->at(previous_seed_index).zipcode_decoder->get_distance_index_address(0) <= 
                            seeds->at(current_item.value).zipcode_decoder->get_distance_index_address(0));
                    
                }  else if (seeds->at(previous_seed_index).zipcode_decoder->get_code_type(depth-1) == ZipCode::CHAIN 
                            || seeds->at(previous_seed_index).zipcode_decoder->get_code_type(depth-1) == ZipCode::ROOT_CHAIN) {
#ifdef DEBUG_ZIP_CODE_TREE
                    cerr << "\t they are children of a common chain" << endl;
#endif
                    //If previous_seed_index and current_item.value are both children of a chain
                    size_t offset_a = seeds->at(previous_seed_index).zipcode_decoder->get_offset_in_chain(depth);
                    size_t offset_b = seeds->at(current_item.value).zipcode_decoder->get_offset_in_chain(depth);

                    if ( offset_a == offset_b) {
                        //If they have the same prefix sum, then the snarl comes first
                        //They will never be on the same child at this depth
                        if (parent_of_a_is_reversed) {
                            assert(seeds->at(current_item.value).zipcode_decoder->get_code_type(depth) != ZipCode::NODE && 
                                   seeds->at(previous_seed_index).zipcode_decoder->get_code_type(depth) == ZipCode::NODE); 
                        } else {
                            assert( seeds->at(previous_seed_index).zipcode_decoder->get_code_type(depth) != ZipCode::NODE && 
                                    seeds->at(current_item.value).zipcode_decoder->get_code_type(depth) == ZipCode::NODE); 
                        }
                    } else {
                        //Check if the parent chain is reversed and if so, then the order should be reversed
                        //The parent could be reversed if it is in an irregular snarl and the 
                        if (parent_of_a_is_reversed) {
                            assert( offset_b <= offset_a);
                        } else {
                            assert( offset_a <= offset_b);
                        }
                    }
                } else {
#ifdef DEBUG_ZIP_CODE_TREE
                    cerr << "\t they are children of a common snarl" << endl;
#endif
                    // Otherwise, they are children of a snarl
                    // Sort by a topological ordering from the start of the snarl
                    // The ranks of children in snarls are in a topological order, so 
                    // sort on the ranks
                    assert( seeds->at(previous_seed_index).zipcode_decoder->get_rank_in_snarl(depth) <=
                            seeds->at(current_item.value).zipcode_decoder->get_rank_in_snarl(depth));
                } 

            }
            previous_seed_index = current_item.value;
            previous_is_valid = current_is_valid;
        } else if (current_item.type == CHAIN_START) {
            //Chains can't start with edges
            assert(zip_code_tree[i+1].type != EDGE);
        } else if (current_item.type == CHAIN_END) {
            //And can't end with edges
            assert(zip_code_tree[i-1].type != EDGE);
        }
    }

    //Make sure that all snarls/chains are opened and closed in a valid order
    vector<tree_item_type_t> snarl_stack; 
    for (const tree_item_t& item : zip_code_tree) {
        if (item.type == SNARL_START) {
            snarl_stack.push_back(SNARL_START);
        } else if (item.type == CHAIN_START) {
            snarl_stack.push_back(CHAIN_START);
        } else if (item.type == SNARL_END) {
            assert(snarl_stack.back() == SNARL_START);
            snarl_stack.pop_back();
        } else if (item.type == CHAIN_END) {
            assert(snarl_stack.back() == CHAIN_START);
            snarl_stack.pop_back();
        }
    }


    // Go through the zipcode tree and check distances and snarl tree relationships

    //Start from the end of the zip tree and walk left, checking each pair of seeds
    for (auto start_itr_left  = zip_code_tree.rbegin() ; 
         start_itr_left != zip_code_tree.rend() ; ++ start_itr_left ) {
        //Get a reverse iterator to the vector, starting from the end and going left
        if (start_itr_left->type != SEED) {
            continue;
        }

        //The seed that the iterator points to
        const Seed& start_seed = seeds->at(start_itr_left->value);

        //Do we want the distance going left in the node
        //This takes into account the position and the orientation of the tree traversal
        bool start_is_reversed = start_itr_left->is_reversed ? !is_rev(start_seed.pos) : is_rev(start_seed.pos);

        //Walk through the tree starting from the vector iterator going left, and check the distance
        for (reverse_iterator tree_itr_left (start_itr_left, zip_code_tree.rend()) ;
             tree_itr_left != reverse_iterator(zip_code_tree.rend(), zip_code_tree.rend()) ;
             ++tree_itr_left) {
            seed_result_t next_seed_result = *tree_itr_left;
            const Seed& next_seed = seeds->at(next_seed_result.seed);
            const bool next_is_reversed = next_seed_result.is_reverse ? !is_rev(next_seed.pos) : is_rev(next_seed.pos);

            size_t tree_distance = next_seed_result.distance;

            net_handle_t start_handle = distance_index.get_node_net_handle(
                                            id(start_seed.pos),
                                            is_rev(start_seed.pos) != start_is_reversed);
            net_handle_t next_handle = distance_index.get_node_net_handle(
                                            id(next_seed.pos),
                                            is_rev(next_seed.pos) != next_is_reversed);

            size_t index_distance = distance_index.minimum_distance(id(next_seed.pos), next_is_reversed,
                    next_seed_result.is_reverse ? distance_index.minimum_length(next_handle) - offset(next_seed.pos) - 1 
                                                 : offset(next_seed.pos),
                    id(start_seed.pos), start_is_reversed,
                    start_itr_left->is_reversed ? distance_index.minimum_length(start_handle) - offset(start_seed.pos) - 1 
                                                : offset(start_seed.pos)
                    );
            if (index_distance != std::numeric_limits<size_t>::max() && is_rev(next_seed.pos) != next_is_reversed) {
                //If the seed we're starting from got reversed, then subtract 1
                index_distance -= 1;
            }
            if (index_distance != std::numeric_limits<size_t>::max() && is_rev(start_seed.pos) != start_is_reversed) {
                //If the seed we ended at got reversed, then add 1
                index_distance += 1;
            }
            pos_t start_pos = is_rev(start_seed.pos) ? make_pos_t(id(start_seed.pos), false, distance_index.minimum_length(start_handle) - offset(start_seed.pos) - 1)
                                                     : start_seed.pos;
            pos_t next_pos = is_rev(next_seed.pos) ? make_pos_t(id(next_seed.pos), false, distance_index.minimum_length(next_handle) - offset(next_seed.pos) - 1)
                                                     : next_seed.pos;

            bool in_non_dag_snarl = false;
            while (!in_non_dag_snarl && !distance_index.is_root(next_handle)) {
                if ((distance_index.is_snarl(next_handle) && !distance_index.is_dag(next_handle)) 
                    || distance_index.is_root_snarl(next_handle)
                    || distance_index.is_looping_chain(next_handle)
                    || distance_index.is_multicomponent_chain(next_handle)) {
                    in_non_dag_snarl = true;
                }
                next_handle = distance_index.get_parent(next_handle);
            }
            while (!in_non_dag_snarl && !distance_index.is_root(start_handle)) {
                if ((distance_index.is_snarl(start_handle) && !distance_index.is_dag(start_handle)) 
                    || distance_index.is_root_snarl(start_handle)
                    || distance_index.is_looping_chain(start_handle)
                    || distance_index.is_multicomponent_chain(start_handle)) {
                    in_non_dag_snarl = true;
                }
                start_handle = distance_index.get_parent(start_handle);
            }

            if (!in_non_dag_snarl && index_distance < distance_limit) {
                if (start_pos == next_pos) {
                    if (tree_distance != 0 && tree_distance != index_distance) {
                        cerr << "Distance between " << next_seed.pos << (next_is_reversed ? "rev" : "") << " and " << start_seed.pos << (start_is_reversed ? "rev" : "") << endl;
                        cerr << "Tree distance: " << tree_distance << " index distance: " << index_distance << endl;
                    }
                    //This could be off by one if one of the seeds is reversed, but I'm being lazy and just checking against the index
                    assert((tree_distance == 0 || tree_distance == index_distance));
                } else {
                    if (tree_distance != index_distance) {
                        cerr << "Distance between " << next_seed.pos << (next_is_reversed ? "rev" : "") << " and " << start_seed.pos << (start_is_reversed ? "rev" : "") << endl;
                        cerr << "Tree distance: " << tree_distance << " index distance: " << index_distance << endl;
                    }
                    assert(tree_distance == index_distance);
                }
            }

        }
    }
}


ZipCodeTree::iterator::iterator(vector<tree_item_t>::const_iterator begin, vector<tree_item_t>::const_iterator end) : it(begin), end(end) {
    while (this->it != this->end && this->it->type != SEED) {
        // Immediately advance to the first seed
        ++this->it;
    }
}

auto ZipCodeTree::iterator::operator++() -> iterator& {
    ++it;
    while (it != end && it->type != SEED) {
        // Advance to the next seed, or the end.
        ++it;
    }
    return *this;
}

auto ZipCodeTree::iterator::operator==(const iterator& other) const -> bool {
    // Ends don't matter for comparison.
    return it == other.it;
}
    
auto ZipCodeTree::iterator::operator*() const -> oriented_seed_t {
    return {it->value, it->is_reversed};
}

auto ZipCodeTree::iterator::remaining_tree() const -> size_t {
    size_t to_return = end - it - 1;
#ifdef debug_parse
    std::cerr << "From " << &*it << " there are " << to_return << " slots after" << std::endl;
#endif
    return to_return;
}

auto ZipCodeTree::begin() const -> iterator {
    return iterator(zip_code_tree.begin(), zip_code_tree.end());
}

auto ZipCodeTree::end() const -> iterator {
    return iterator(zip_code_tree.end(), zip_code_tree.end());
}

ZipCodeTree::reverse_iterator::reverse_iterator(vector<tree_item_t>::const_reverse_iterator rbegin, vector<tree_item_t>::const_reverse_iterator rend, size_t distance_limit) : it(rbegin), rend(rend), distance_limit(distance_limit), stack(), current_state(S_START) {
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
        // Skip ahead to the first seed we actually want to yield, or to the end of the data.
        ++this->it;
#ifdef debug_parse
        if (this->it != rend) {
            std::cerr << "Able to do another initial tick." << std::endl;
        }
#endif
    }
    // As the end of the constructor, the iterator points to a seed that has been ticked and yielded, or is rend.
#ifdef debug_parse
    if (this->it == rend) {
        std::cerr << "Ran out of tree looking for first seed." << std::endl;
    }
#endif
}

auto ZipCodeTree::reverse_iterator::operator++() -> reverse_iterator& {
    // Invariant: the iterator points to a seed that has been ticked and yielded, or to rend.
    if (it != rend) {
#ifdef debug_parse
        std::cerr << "Skipping over a " << it->type << " which we assume was handled already." << std::endl;
#endif
        ++it;

    }
    while (it != rend && !tick()) {
        // Skip ahead to the next seed we actually want to yield, or to the end of the data.
        ++it;
    }
#ifdef debug_parse
    if (it == rend) {
        std::cerr << "Ran out of tree looking for next seed." << std::endl;
    }
#endif
    return *this;
}

auto ZipCodeTree::reverse_iterator::operator==(const reverse_iterator& other) const -> bool {
    // Ends and other state don't matter for comparison.
    return it == other.it;
}

auto ZipCodeTree::reverse_iterator::operator*() const -> seed_result_t {
    // We are always at a seed, so show that seed
    crash_unless(it != rend);
    crash_unless(it->type == SEED);
    crash_unless(!stack.empty());
    // We know the running distance to this seed will be at the top of the stack.
    seed_result_t to_return;
    to_return.seed = it->value;
    to_return.is_reverse = it->is_reversed;
    to_return.distance = stack.top();
    return to_return;
}

auto ZipCodeTree::reverse_iterator::push(size_t value) -> void {
    stack.push(value);
}

auto ZipCodeTree::reverse_iterator::pop() -> size_t {
    size_t value = stack.top();
    stack.pop();
    return value;
}

auto ZipCodeTree::reverse_iterator::top() -> size_t& {
    crash_unless(depth() > 0);
    return stack.top();
}

auto ZipCodeTree::reverse_iterator::dup() -> void {
    push(stack.top());
}

auto ZipCodeTree::reverse_iterator::depth() const -> size_t {
    return stack.size();
}

auto ZipCodeTree::reverse_iterator::swap() -> void {
    // Grab the top item
    size_t temp = stack.top();
    stack.pop();
    // Swap it with what was under it
    std::swap(temp, stack.top());
    // And put that back on top
    stack.push(temp);
}

auto ZipCodeTree::reverse_iterator::state(State new_state) -> void {
    current_state = new_state;
}

auto ZipCodeTree::reverse_iterator::halt() -> void {
#ifdef debug_parse
    std::cerr << "Halt iteration!" << std::endl;
#endif
    it = rend;
}

auto ZipCodeTree::reverse_iterator::tick() -> bool {
#ifdef debug_parse
    std::cerr << "Tick for state " << current_state << " on symbol " << it->type << " at " << &*it << std::endl;
#endif
    switch (current_state) {
    case S_START:
        // Initial state.
        //
        // Stack is empty and we must be at a seed to start at.
        switch (it->type) {
        case SEED:
#ifdef debug_parse
            std::cerr << "Skip over seed " << it->value << std::endl;
#endif
            push(0);
            state(S_SCAN_CHAIN);
            break;
        default:
            throw std::domain_error("Unimplemented symbol " + std::to_string(it->type) + " for state " + std::to_string(current_state)); 
        }
        break;
    case S_SCAN_CHAIN:
        // State where we are scanning a chain leftward up to its start.
        //
        // Stack has at the top the running distance along the chain, and under
        // that running distances to use at the other chains in the snarl, and
        // under that running distances to use for the other chains in the
        // snarl's parent snarl, etc.
        switch (it->type) {
        case SEED:
            // Emit seed here with distance at top of stack.
            crash_unless(depth() > 0);
#ifdef debug_parse
            std::cerr << "Yield seed " << it->value << ", distance " << top() << std::endl;
#endif
            return true;
            break;
        case SNARL_END:
            // Running distance along chain is on stack, and will need to be added to all the stored distances.
            state(S_STACK_SNARL); // Stack up pre-made scratch distances for all the things in the snarl.
            break;
        case CHAIN_START:
            if (depth() == 1) {
                // We never entered the parent snarl of this chain, so stack up
                // the distances left of here as options added to the
                // distance along this chain.
                //
                // Running distance along chain is on stack, and will need to
                // be added to all the stored distances.
                state(S_STACK_SNARL);
            } else {
                // We did enter the parent snarl already.
                // Discard the running distance along this chain, which no longer matters.
                pop();
                // Running distance for next chain, or running distance to cross the snarl, will be under it.
                state(S_SCAN_SNARL);
            }
            break;
        case EDGE:
            // Distance between things in a chain.
            // Add value into running distance, maxing it if value is max.
            top() = SnarlDistanceIndex::sum(top(), it->value);
            if (top() > distance_limit || top() == std::numeric_limits<size_t>::max()) {
                // Skip over the rest of this chain
                if (depth() == 1) {
                    // We never entered the parent snarl of this chain.
                    // So if the distance along the chain is too much, there
                    // are not going to be any results with a smaller distance.
                    halt();
                    // When we halt we have to return true to show the halting position.
                    return true;
                } else {
                    // We need to try the next thing in the parent snarl, so skip the rest of the chain.
                    // We're skipping in 0 nested snarls right now.
                    push(0);
                    state(S_SKIP_CHAIN);
                }
            }
            break;
        default:
            throw std::domain_error("Unimplemented symbol " + std::to_string(it->type) + " for state " + std::to_string(current_state)); 
        }
        break;
    case S_STACK_SNARL:
        // State where we are stacking up the stored edge values, the first
        // time we get to a particular snarl.
        //
        // Stack has the running distance along the parent chain, and under
        // that the stacked running distances for items in the snarl.
        switch (it->type) {
        case EDGE:
            // We need to add this actual number to parent running distance.
            // Duplicate parent running distance
            dup();
            // Add in the edge value to make a running distance for the thing this edge is for.
            // Account for if the edge is actually unreachable.
            top() = SnarlDistanceIndex::sum(top(), it->value);
            // Flip top 2 elements, so now parent running distance is on top, over edge running distance.
            swap();
            break;
        case CHAIN_END:
            // Throw out parent running distance
            pop();
            if (depth() == 0) {
                // We left a chain and immediately entered a chain without a distance.
                // This means the chains aren't actually connected.
                halt();
                // When we halt we have to return true to show the halting position.
                return true;
            } else {
                // So now we have the running distance for this next chain.
                if (top() > distance_limit || top() == std::numeric_limits<size_t>::max()) {
                    // Running distance is already too high so skip over the chain
                    push(0);
                    state(S_SKIP_CHAIN);
                } else {
                    // Do the chain
                    state(S_SCAN_CHAIN);
                }
            }
            break;
        case SNARL_START:
            // We didn't hit another chain in the snarl, we hit the start of
            // the snarl. We should have stacked exactly one distance.

            // Throw out parent running distance
            pop();

            // There should be a running distance on the stack still, and we
            // will continue with that in the parent chain.
            crash_unless(depth() > 0);
            state(S_SCAN_CHAIN);
            break;
        case NODE_COUNT:
            // We've found the node count in the snarl. We don't need it, so
            // skip it.
            // TODO: Use it if skipping the snarl.
            break;
        default:
            throw std::domain_error("Unimplemented symbol " + std::to_string(it->type) + " for state " + std::to_string(current_state)); 
        }
        break;
    case S_SCAN_SNARL:
        // State where we are going through a snarl and doing all its chains.
        //
        // Stack has at the top running distances to use for each chain still
        // to be visited in the snarl, and under those the same for the snarl
        // above that, etc.
        switch (it->type) {
        case SNARL_START:
            // Stack holds running distance along parent chain plus edge
            // distance to cross the snarl, or running distance out of chain we
            // started in plus distance to exit the snarl.
            //
            // This is the right running distance to use for the parent chain now.
            // So go back to scanning the parent chain.
            state(S_SCAN_CHAIN);
            break;
        case CHAIN_END:
            // We've encountered a chain to look at, and the running distance
            // into the chain is already on the stack.
            if (top() > distance_limit || top() == std::numeric_limits<size_t>::max()) {
                // Running distance is already too high so skip over the chain
                push(0);
                state(S_SKIP_CHAIN);
            } else {
                // Do the chain
                state(S_SCAN_CHAIN);
            }
            break;
        case EDGE:
            // We've found edge data in the snarl, but we already know the
            // running distances to everything we will encounter, so we ignore
            // it.
            break;
        case NODE_COUNT:
            // We've found the node count in the snarl. We don't need it, so
            // skip it.
            break;
        default:
            throw std::domain_error("Unimplemented symbol " + std::to_string(it->type) + " for state " + std::to_string(current_state)); 
        }
        break;
    case S_SKIP_CHAIN:
        // State where we are skipping over the rest of a chain because we hit
        // the distance limit, but we might need to do other chains in a parent
        // snarl.
        //
        // Stack has the nesting level of child snarls we are reading over
        // until we get back to the level we want to skip past the chain
        // start.
        // Under that is the running distance along the chain being skipped.
        // And under that it has the running distance for ther next thing in
        // the snarl, which had better exist or we shouldn't be trying to skip
        // the chain, we should have halted.
        switch (it->type) {
        case SEED:
            // We don't emit seeds until the chain is over
            return false;
            break;
        case SNARL_START:
            // We might now be able to match chain starts again
            top() -= 1;
            break;
        case SNARL_END:
            // We can't match chain starts until we leave the snarl
            top() += 1;
            break;
        case CHAIN_START:
            if (top() == 0) {
                // This is the start of the chain we were wanting to skip.
                pop();
                // We definitely should have entered the parent snarl of the chain, or we would have halted instead of trying to skip the rest of the chain.
                crash_unless(depth() > 1);
                // Discard the running distance along this chain, which no longer matters.
                pop();
                // Running distance for next chain, or running distance to cross the snarl, will be under it.
                state(S_SCAN_SNARL);
            }
            // Otherwise this is the start of a chain inside a child snarl we are skipping over and we ignore it.
            break;
        case CHAIN_END:
            // Ignore chain ends
            break;
        case EDGE:
            // Ignore edge values
            break;
        case NODE_COUNT:
            // Ignore node counts
            // TODO: We should read these and jump along instead!
            break;
        default:
            throw std::domain_error("Unimplemented symbol " + std::to_string(it->type) + " for state " + std::to_string(current_state)); 
        }
        break;
    default:
        throw std::domain_error("Unimplemented state " + std::to_string(current_state)); 
    }
    // Unless we yield something, we don't want to pause the scan here.
    return false;
}

auto ZipCodeTree::look_back(const iterator& from, size_t distance_limit) const -> reverse_iterator {
    return reverse_iterator(zip_code_tree.rbegin() + from.remaining_tree(), zip_code_tree.rend(), distance_limit);
}
auto ZipCodeTree::rend() const -> reverse_iterator {
    return reverse_iterator(zip_code_tree.rend(), zip_code_tree.rend(), 0);
}


std::ostream& operator<<(std::ostream& out, const ZipCodeTree::tree_item_type_t& type) {
    return out << std::to_string(type);
}

std::ostream& operator<<(std::ostream& out, const ZipCodeTree::reverse_iterator::State& state) {
    return out << std::to_string(state);
}

vector<size_t> ZipCodeForest::sort_seeds_by_zipcode(const SnarlDistanceIndex& distance_index) const {
    /*
      Sort the seeds in roughly linear/topological-ish order along the top-level chains

      This sorts the seeds top-down along the snarl tree
      Sorting is split into two different types of sort: radix sort or an n-log-n sort,
      depending on which will be more efficient
      Sorting begins at the root, with a radix sort to partition the seeds into connected component
      For each partition (a chain presumably), an n-log-n sort will be done to sort along the chain
      And so on down the snarl tree.
      The two sorters will each sort on a slice of the vector and update a new list of slices for the next 
      level in the snarl tree 
    */

    //Helper function to get the value to sort on from the zipcode
    //This doesn't take into account the orientation, except for nodes offsets in chains
    //It will actually be defined somewhere else
    //Used for sorting at the given depth, so use values at depth depth+1
    auto get_sort_value = [&] (Seed& seed, size_t depth) {
#ifdef DEBUG_ZIP_CODE_TREE
        cerr << "\tGet the sort value of seed " << seed.pos << " at depth " << depth << endl;
#endif
        ZipCode::code_type_t code_type = seed.zipcode_decoder->get_code_type(depth);
        if (code_type == ZipCode::NODE || code_type == ZipCode::ROOT_NODE || seed.zipcode_decoder->max_depth() == depth) {
#ifdef DEBUG_ZIP_CODE_TREE
            cerr << "\t\t this is a node: offset: " << ( is_rev(seed.pos) ? seed.zipcode_decoder->get_length(depth) - offset(seed.pos) - 1
                                    : offset(seed.pos)) << endl;;
#endif
            return is_rev(seed.pos) ? seed.zipcode_decoder->get_length(depth) - offset(seed.pos) - 1
                                    : offset(seed.pos);
        } else if (code_type == ZipCode::CHAIN || code_type == ZipCode::ROOT_CHAIN) {
#ifdef DEBUG_ZIP_CODE_TREE
            cerr << "\t\t this is a chain: prefix sum value x2 (and -1 if snarl): ";
#endif
            //Return the prefix sum in the chain
            //In order to accommodate nodes and snarls that may have the same prefix sum value, actually uses
            //the prefix sum value * 2, and subtracts 1 in this is a snarl, to ensure that it occurs 
            //before the node with the same prefix sum value 
            size_t prefix_sum;
            if (seed.zipcode_decoder->get_code_type(depth+1) == ZipCode::REGULAR_SNARL || seed.zipcode_decoder->get_code_type(depth+1) == ZipCode::IRREGULAR_SNARL) { 
                //If this is a snarl, then get the prefix sum value*2 - 1
                prefix_sum = (seed.zipcode_decoder->get_offset_in_chain(depth+1) * 2) - 1;
            } else {
                //If this is a node, then get the prefix sum value plus the offset in the position, and multiply by 2 
                prefix_sum = SnarlDistanceIndex::sum(seed.zipcode_decoder->get_offset_in_chain(depth+1),
                                                 seed.zipcode_decoder->get_is_reversed_in_parent(depth+1) != is_rev(seed.pos)
                                                     ? seed.zipcode_decoder->get_length(depth+1) - offset(seed.pos) - 1
                                                     : offset(seed.pos));
                prefix_sum *= 2;
            }
#ifdef DEBUG_ZIP_CODE_TREE
            cerr << prefix_sum << endl;
#endif
            return prefix_sum;
        } else {
#ifdef DEBUG_ZIP_CODE_TREE
            cerr << "\tThis is snarl, so return the rank in the snarl: " << seed.zipcode_decoder->get_rank_in_snarl(depth+1) << endl;
#endif
            // The ranks of children in irregular snarls are in a topological order, so 
            // sort on the ranks
            // The rank of children in a regular snarl is arbitrary but it doesn't matter anyway
            return seed.zipcode_decoder->get_rank_in_snarl(depth+1);
        }
    };

    //At the given depth, go through sort_order in the given interval to find the intervals for the next level 
    //and add to new_intervals
    auto find_next_intervals = [&] (const interval_and_orientation_t& interval,
                                    size_t depth, const vector<size_t>& sort_order, 
                                    vector<interval_and_orientation_t>& new_intervals,
                                    const std::function<size_t(Seed& seed, size_t depth)>& get_partitioning_value) {
        //Now that it's sorted, find runs of equivalent values for new_interval_to_sort
        //Also need to check the orientation
        size_t start_of_current_run = interval.interval_start;
        for (size_t i = interval.interval_start+1 ; i < interval.interval_end ; i++) {
            
            //If the current seed is a node and has nothing at depth+1 or is different from the previous seed at this depth
            bool is_node = seeds->at(sort_order[i]).zipcode_decoder->max_depth() == depth ||
                           seeds->at(sort_order[i]).zipcode_decoder->get_code_type(depth+1) == ZipCode::NODE;
            bool is_different_from_previous = get_partitioning_value(seeds->at(sort_order[i]), depth) 
                                              != get_partitioning_value(seeds->at(sort_order[i-1]), depth);
            bool is_last = i == interval.interval_end-1;
            if (is_different_from_previous && i-1 != start_of_current_run) {
                //If this is the end of a run of more than one thing
                //If the previous thing was a node, then start_of_current_run would have been set to i-1, so
                //it won't reach here

                bool current_is_reversed = ZipCodeTree::seed_is_reversed_at_depth(seeds->at(sort_order[i-1]), depth+1, distance_index) 
                                         ? !interval.is_reversed
                                         : interval.is_reversed;
                new_intervals.emplace_back(start_of_current_run,  i, current_is_reversed);
                
                start_of_current_run = i;
            } else if (is_last && !is_different_from_previous && !is_node) {
                //If this is the last thing in the sorted list, and the previous thing was in the same run

                bool current_is_reversed = ZipCodeTree::seed_is_reversed_at_depth(seeds->at(sort_order[i-1]), depth+1, distance_index) 
                                         ? !interval.is_reversed
                                         : interval.is_reversed;
                new_intervals.emplace_back(start_of_current_run, i+1, current_is_reversed);
            } else if (is_node || is_different_from_previous) {
                start_of_current_run = i;
            }
        }
    };

    //The sort order of the seeds. Each element is an index into seeds
    //Initialized to the current order of the seeds, and gets updated as sorting happens
    vector<size_t> zipcode_sort_order (seeds->size(), 0);
    for (size_t i = 0 ; i < zipcode_sort_order.size() ; i++) {
        zipcode_sort_order[i] = i;
    }

    //A vector of ranges in zipcode_sort_order that need to be sorted
    //This gets updated as sorting precedes through each level of the snarl tree
    vector<interval_and_orientation_t> intervals_to_sort;


    //Depth of the snarl tree
    size_t depth = 0;


    //First sort everything by connected component of the root
    // Assume that the number of connected components is small enough that radix sort is more efficient
    interval_and_orientation_t first_interval(0, zipcode_sort_order.size(), false);
    radix_sort_zipcodes(zipcode_sort_order, first_interval, 
                        false, std::numeric_limits<size_t>::max(), distance_index,
                        [&](Seed& seed, size_t depth) {
                            //Sort on the connected component number
                            return seed.zipcode_decoder->get_distance_index_address(0);
                        });

#ifdef DEBUG_ZIP_CODE_TREE
        cerr << "After root " << endl;
        for (size_t i : zipcode_sort_order) {
            cerr << i << ":" << seeds->at(i).pos << ", ";
        }
        cerr << endl;
#endif
    find_next_intervals(first_interval, std::numeric_limits<size_t>::max(), zipcode_sort_order, intervals_to_sort,
                        [&](Seed& seed, size_t depth) {
                            //Sort on the connected component number
                            return seed.zipcode_decoder->get_distance_index_address(0);
                        });

    //While there is still stuff to sort, walk down the snarl tree and sort each interval for each depth
    while (!intervals_to_sort.empty()) {
#ifdef DEBUG_ZIP_CODE_TREE
        cerr << "Sort seeds at depth " << depth << endl;
#endif

        //The intervals to sort at the next level of the snarl tree. To be filled in in this iteration
        vector<interval_and_orientation_t> new_intervals_to_sort;

        for (const interval_and_orientation_t& current_interval : intervals_to_sort) {
#ifdef DEBUG_ZIP_CODE_TREE
            cerr << "Sort seeds on interval " << current_interval.interval_start << "-" << current_interval.interval_end << endl;
#endif

            // Sorting will either be done with radix sort or with std::sort, depending on which is more efficient
            // Radix sort is linear time in the number of items it is sorting, but also linear space in the range 
            // of the values it is sorting on
            // If the range of values is greater than the n log n (in the number of things being sorted) of the default
            // sorter, then use radix

            bool use_radix;

            //One of the seeds getting sorted
            const Seed& seed_to_sort = seeds->at(zipcode_sort_order[current_interval.interval_start]);
            auto current_type = seed_to_sort.zipcode_decoder->get_code_type(depth);

            if (current_type  == ZipCode::ROOT_CHAIN) {
                //IF this is a root chain, then use the default sort, because it's probably too big for radix and we can't tell
                //anyways because we don't store the length of a root-chain
                use_radix = false;
            } else if (current_type == ZipCode::NODE || current_type == ZipCode::CHAIN) {
                //If we're sorting a node or chain, then the range of values is the minimum length of the node/chain
                // times 2 because it gets multiplied by 2 to differentiate nodes and snarls
                size_t radix_cost = seed_to_sort.zipcode_decoder->get_length(depth) * 2;
                size_t default_cost = (current_interval.interval_end - current_interval.interval_start) * std::log2(current_interval.interval_end - current_interval.interval_start);

                use_radix = radix_cost < default_cost;
            } else {
                //Otherwise, this is a snarl and the range of values is the number of children in the snarl
                //TODO: Since the zipcodes don't store this, and I'm pretty sure it will be small, for now default to radix

                use_radix = true;
            }

            bool reverse_order = (current_type == ZipCode::REGULAR_SNARL || current_type == ZipCode::IRREGULAR_SNARL) 
                                 ? false
                                 : current_interval.is_reversed; 

            if (use_radix) {
                //Sort the given interval using the value-getter and orientation
                radix_sort_zipcodes(zipcode_sort_order, current_interval, reverse_order, depth, distance_index, get_sort_value);
            } else {
                //Sort the given interval using the value-getter and orientation
                default_sort_zipcodes(zipcode_sort_order, current_interval, reverse_order, depth, distance_index, get_sort_value);
            }

            find_next_intervals(current_interval, depth, zipcode_sort_order, new_intervals_to_sort, get_sort_value);

        }

        //Update to the next depth
        intervals_to_sort = std::move(new_intervals_to_sort);
        depth++;
#ifdef DEBUG_ZIP_CODE_TREE
        cerr << "Order after depth " << depth-1 << endl;
        for (size_t i = 0 ; i < zipcode_sort_order.size() ; i++) {
            cerr << i << ":" << seeds->at(zipcode_sort_order[i]).pos << ", ";
        }
        cerr << endl;
#endif
    }

    return zipcode_sort_order;
}

void ZipCodeForest::radix_sort_zipcodes(vector<size_t>& zipcode_sort_order, const interval_and_orientation_t& interval,
                                   bool reverse_order, size_t depth, const SnarlDistanceIndex& distance_index, 
                                   const std::function<size_t(Seed& seed, size_t depth)>& get_sort_value) const {
    //Radix sort the interval of zipcode_sort_order in the given interval
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "\tradix sort" << endl;
#endif

    //Mostly copied from Jordan Eizenga

    // count up occurrences of each rank
    std::vector<size_t> counts;
    for (size_t i  = interval.interval_start ; i < interval.interval_end ; i++) {
        size_t next_rank = get_sort_value(seeds->at(zipcode_sort_order[i]), depth) + 1;

        while (counts.size() <= next_rank) {
            counts.push_back(0);
        }
        ++counts[next_rank];
    }
    
    //Make this a count of the number of things before it
    for (size_t i = 1; i < counts.size(); ++i) {
        counts[i] += counts[i - 1];
    }

    //Get the sorted order
    std::vector<size_t> sorted(interval.interval_end - interval.interval_start);
    for (size_t i = interval.interval_start ; i < interval.interval_end ; i++) {
        size_t rank = get_sort_value(seeds->at(zipcode_sort_order[i]), depth);
        sorted[counts[rank]++] = zipcode_sort_order[i];
    }
    
    //And place everything in the correct position
    for (size_t i = 0 ; i < sorted.size() ; i++) {


        //If this is reversed in the top-level chain, then the order should be backwards
        //TODO: I'm not sure how this should work for a snarl
        if (reverse_order) {
            zipcode_sort_order[interval.interval_end - i - 1] = sorted[i];
        } else {
            zipcode_sort_order[i + interval.interval_start] = sorted[i];
        }
    }

}
void ZipCodeForest::default_sort_zipcodes(vector<size_t>& zipcode_sort_order, const interval_and_orientation_t& interval,
                                   bool reverse_order, size_t depth, const SnarlDistanceIndex& distance_index, 
                                   const std::function<size_t(Seed& seed, size_t depth)>& get_sort_value) const { 
    //std::sort the interval of zipcode_sort_order between interval_start and interval_end
    
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "\tdefault sort between " << interval.interval_start  << " and " << interval.interval_end  << endl;
    cerr << "\tis rev: " << reverse_order << endl;
#endif
    //Sort using std::sort 
    std::sort(zipcode_sort_order.begin() + interval.interval_start, zipcode_sort_order.begin() + interval.interval_end, [&] (size_t a, size_t b) {
        //If this snarl tree node is reversed, then reverse the sort order
        return reverse_order ? get_sort_value(seeds->at(a), depth) > get_sort_value(seeds->at(b), depth)
                             : get_sort_value(seeds->at(a), depth) < get_sort_value(seeds->at(b), depth);
    });
}

}

namespace std {

std::string to_string(const vg::ZipCodeTree::tree_item_type_t& type) {
    switch (type) {
    case vg::ZipCodeTree::SEED:
        return "SEED";
    case vg::ZipCodeTree::SNARL_START:
        return "SNARL_START";
    case vg::ZipCodeTree::SNARL_END:
        return "SNARL_END";
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

std::string to_string(const vg::ZipCodeTree::reverse_iterator::State& state) {
    switch (state) {
    case vg::ZipCodeTree::reverse_iterator::S_START:
        return "S_START";
    case vg::ZipCodeTree::reverse_iterator::S_SCAN_CHAIN:
        return "S_SCAN_CHAIN";
    case vg::ZipCodeTree::reverse_iterator::S_STACK_SNARL:
        return "S_STACK_SNARL";
    case vg::ZipCodeTree::reverse_iterator::S_SCAN_SNARL:
        return "S_SCAN_SNARL";
    case vg::ZipCodeTree::reverse_iterator::S_SKIP_CHAIN:
        return "S_SKIP_CHAIN";
    default:
        throw std::runtime_error("Unimplemented zip code tree reverse iterator state");
    }
}



}

