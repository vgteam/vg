//#define DEBUG_ZIP_CODE_TREE
//#define PRINT_NON_DAG_SNARLS
//#define DEBUG_ZIP_CODE_SORTING

#include "zip_code_tree.hpp"

#include "crash.hpp"

//#define debug_parse

using namespace std;
namespace vg {

void ZipCodeForest::fill_in_forest(const vector<Seed>& all_seeds, const SnarlDistanceIndex& distance_index,
                               size_t distance_limit) {
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "Make a new forest with " << all_seeds.size() << " seeds with distance limit " << distance_limit << endl;
#endif
    if (all_seeds.size() == 0) {
        return;
    }
    seeds = &all_seeds;

    /*
    Make a ZipCodeForest
    Takes a vector of seeds and fills in the forest

    The zip forest is made by sorting the seeds along chains/snarls, 
    then adding each seed, snarl/chain boundary, and distance to zip_code_tree

    Sorting and tree-making is done at the same time, in a depth-first traversal of the snarl tree
    Sorting is done for node in the snarl tree, and splits the seeds up into children of that node.
    After sorting, the new children are added to a stack of children to be sorted and processed
    A child is processed by opening it in the zip tree along with any relevant distances, and uj 
    */

    //Start by initializing the state
    forest_growing_state_t forest_state;
    forest_state.active_zip_tree = std::numeric_limits<size_t>::max();

    forest_state.seed_sort_order.assign(seeds->size(), 0);
    for (size_t i = 0 ; i < forest_state.seed_sort_order.size() ; i++) {
        forest_state.seed_sort_order[i] = i;
    }

    //Start with the root
    interval_and_orientation_t first_interval (0, seeds->size(), false, ZipCode::EMPTY);
    //Get the intervals of the connected components
    vector<interval_and_orientation_t> new_intervals = sort_one_interval(forest_state.seed_sort_order, first_interval, 0, distance_index);;
    forest_state.intervals_to_process.insert(forest_state.intervals_to_process.end(),
                                             new_intervals.rbegin(),
                                             new_intervals.rend());


    while (!forest_state.intervals_to_process.empty()) {
#ifdef DEBUG_ZIP_CODE_TREE
        print_self();
#endif
        // For each unprocessed interval, process it
        // First, check if anything needs to be closed, which will happen if the interval_end in an open snarl/chains
        // gets reached or exceeded 
        // Get the intervals of this interval's children and add them in reverse order to the stack intervals_to_process
        // Then, add any extra seeds or distances between this interval and the previous child -
        //    for snarls that are children of chains, check if there are seeds that need to get added
        //    for chains that are children of snarls, add distances in snarl
        // Open the current interval's snarl/chain


        //Get the interval
        interval_and_orientation_t current_interval = std::move(forest_state.intervals_to_process.back());
        forest_state.intervals_to_process.pop_back();

        /*********
         * First, check if anything needs to be closed and close it
         ********/
#ifdef DEBUG_ZIP_CODE_TREE
        cerr << "Process interval of type " << current_interval.code_type << " with range " << current_interval.interval_start << "-" << current_interval.interval_end << endl;
        cerr << "Close anything open" << endl;
#endif
        while (!forest_state.open_intervals.empty()) {
            if (forest_state.open_intervals.back().interval_end <= current_interval.interval_start) {
                //If the range of the this interval comes after the range in the open interval,
                //close the last thing in open_intervals

#ifdef DEBUG_ZIP_CODE_TREE
cerr << "\tclose something at depth " << forest_state.open_intervals.size()-1 << endl;
#endif

                size_t depth = forest_state.open_intervals.size()-1;

                //The last seed in the thing to close 
                const interval_and_orientation_t& ancestor_interval = forest_state.open_intervals.back();
                const Seed& last_seed = seeds->at(forest_state.seed_sort_order[ancestor_interval.interval_end-1]);

                if (ancestor_interval.code_type == ZipCode::CHAIN ||
                    ancestor_interval.code_type == ZipCode::NODE ||
                    ancestor_interval.code_type == ZipCode::ROOT_CHAIN ||
                    ancestor_interval.code_type == ZipCode::ROOT_NODE) {
                    //Close a chain

                    close_chain(forest_state, distance_index, distance_limit, depth, 
                                last_seed, ancestor_interval.is_reversed); 
                } else {
#ifdef DEBUG_ZIP_CODE_TREE
                    assert(ancestor_interval.code_type == ZipCode::REGULAR_SNARL ||
                           ancestor_interval.code_type == ZipCode::IRREGULAR_SNARL ||
                           ancestor_interval.code_type == ZipCode::CYCLIC_SNARL ||
                           ancestor_interval.code_type == ZipCode::ROOT_SNARL);
#endif
                    //Close a snarl
                    close_snarl(forest_state, distance_index, depth, last_seed, 
                                ancestor_interval.is_reversed, ancestor_interval.code_type == ZipCode::CYCLIC_SNARL); 
                }

                //Clear the list of children of the thing at this level
                forest_state.sibling_indices_at_depth[depth].clear();

                //Take out this ancestor
                forest_state.open_intervals.pop_back();
            } else {
                //If the current interval is contained in this open interval, then it is also contained in all other
                // ancestors so break
                break;
            }
        }

        /************ 
         * Now start processing the current interval
         *
         *
         * Sort this interval and add the child intervals in reverse order to intervals_to_process  
         ***********/

        
        // The depth of the current interval
        size_t current_depth = forest_state.open_intervals.size();

        if (current_interval.code_type != ZipCode::NODE ) {
            //Sort the current interval and get the intervals corresponding to its children
            vector<interval_and_orientation_t> child_intervals = sort_one_interval(forest_state.seed_sort_order, current_interval, current_depth, distance_index);

            //Add the child intervals to the to_process stack, in reverse order so the first one gets popped first
            forest_state.intervals_to_process.insert(forest_state.intervals_to_process.end(),
                                                     child_intervals.rbegin(),
                                                     child_intervals.rend());
        }
        if (current_interval.code_type == ZipCode::CYCLIC_SNARL) {
            // For cyclic snarls, the orientation is set after sorting the parent chain. 
            // The orientation of a cyclic snarl is the direction that the read takes in a start-to-end traversal of
            // the snarl, but this is only necessary for sorting the snarl and finding its children. After that,
            // the snarl should have the orientation of its parent chain so that the distances will be found properly

            current_interval.is_reversed = forest_state.open_intervals.back().is_reversed;
        }
    
        
        /**********
         * Open the current interval
         * If the current interval is a snarl and a child of a chain, then add the preceding sibling seeds before the snarl
         *******/
#ifdef DEBUG_ZIP_CODE_TREE
         cerr << "Open next interval or (if the interval is for nodes), add seeds" << endl;
#endif
        if (forest_state.open_intervals.size()+1 > forest_state.sibling_indices_at_depth.size()) {
            forest_state.sibling_indices_at_depth.emplace_back();
        }
        if (forest_state.open_intervals.empty()) {
            // If there is nothing open, then this is starting a new connected component
            // Just open it
#ifdef DEBUG_ZIP_CODE_TREE
            assert(current_interval.code_type == ZipCode::ROOT_NODE ||
                   current_interval.code_type == ZipCode::NODE ||
                   current_interval.code_type == ZipCode::ROOT_CHAIN ||
                   current_interval.code_type == ZipCode::ROOT_SNARL);
#endif

            // Start a new connected component
            if (forest_state.active_zip_tree == std::numeric_limits<size_t>::max() 
                || trees[forest_state.active_zip_tree].zip_code_tree.size() != 0) {
                trees.emplace_back(seeds);
                forest_state.active_zip_tree = trees.size()-1;
            }

            if (current_interval.code_type == ZipCode::ROOT_SNARL) {
                // Open the root snarl
                open_snarl(forest_state, 0);
            } else if (current_interval.code_type == ZipCode::NODE) {
                //For a root node, just add the chain and all the seeds

                trees[forest_state.active_zip_tree].zip_code_tree.push_back({ZipCodeTree::CHAIN_START, std::numeric_limits<size_t>::max(), false});

                //Remember the start of the chain
                forest_state.sibling_indices_at_depth[0].push_back({ZipCodeTree::CHAIN_START, 0});

                //If this is a node, then the interval contains everything in it, so add the seeds and close the chain here
                for (size_t seed_i = current_interval.interval_start ; seed_i < current_interval.interval_end ; seed_i++) {

                    add_child_to_chain(forest_state, distance_index, distance_limit, current_depth, 
                                       forest_state.seed_sort_order[seed_i], current_interval.is_reversed,
                                       current_interval.is_reversed ); 
                }
                close_chain(forest_state, distance_index, distance_limit, current_depth, 
                            seeds->at(forest_state.seed_sort_order[current_interval.interval_end-1]), current_interval.is_reversed); 

                
            } else {
                // Open the root chain/node
                trees[forest_state.active_zip_tree].zip_code_tree.push_back({ZipCodeTree::CHAIN_START, std::numeric_limits<size_t>::max(), false});

                //Remember the start of the chain
                forest_state.sibling_indices_at_depth[0].push_back({ZipCodeTree::CHAIN_START, 0});
            }                       
        } else if (forest_state.open_intervals.back().code_type == ZipCode::CHAIN || 
                   forest_state.open_intervals.back().code_type == ZipCode::ROOT_CHAIN ||
                   forest_state.open_intervals.back().code_type == ZipCode::ROOT_NODE) {
            // This is the child of a chain

            if (current_interval.code_type == ZipCode::NODE) {
                // If the type of this interval is NODE, then this is a range of seeds that are on nodes on the chain,
                // not necessarily on the same node
                // Add each seed

                for (size_t seed_i = current_interval.interval_start ; seed_i < current_interval.interval_end ; seed_i++) {

                    if (current_depth-1 == seeds->at(forest_state.seed_sort_order[seed_i]).zipcode_decoder->max_depth()) {
                        //If this is getting added to a node
                        add_child_to_chain(forest_state, distance_index, distance_limit, current_depth-1, 
                                           forest_state.seed_sort_order[seed_i], current_interval.is_reversed,
                                           forest_state.open_intervals.back().is_reversed ); 
                    } else {
                        add_child_to_chain(forest_state, distance_index, distance_limit, current_depth, 
                                           forest_state.seed_sort_order[seed_i], current_interval.is_reversed,
                                           forest_state.open_intervals.back().is_reversed ); 
                    }
                }

            } else {
#ifdef DEBUG_ZIP_CODE_TREE
                assert(current_interval.code_type == ZipCode::REGULAR_SNARL || 
                       current_interval.code_type == ZipCode::IRREGULAR_SNARL ||
                       current_interval.code_type == ZipCode::CYCLIC_SNARL);
#endif

                //Add the snarl to the chain
                add_child_to_chain(forest_state, distance_index, distance_limit, current_depth,
                                   forest_state.seed_sort_order[current_interval.interval_start], 
                                   current_interval.is_reversed, forest_state.open_intervals.back().is_reversed);
            }
            

        } else {
            //If there is an open ancestor that isn't a chain, so the ancestor must be a snarl
#ifdef DEBUG_ZIP_CODE_TREE
            assert(forest_state.open_intervals.back().code_type == ZipCode::REGULAR_SNARL ||
                   forest_state.open_intervals.back().code_type == ZipCode::IRREGULAR_SNARL ||
                   forest_state.open_intervals.back().code_type == ZipCode::CYCLIC_SNARL ||
                   forest_state.open_intervals.back().code_type == ZipCode::ROOT_SNARL);
#endif

            //Open the child chain
            open_chain(forest_state, distance_index, distance_limit, forest_state.open_intervals.size(), 
                       seeds->at(forest_state.seed_sort_order[current_interval.interval_start]), current_interval.is_reversed);
            
        }

        if (current_interval.code_type != ZipCode::NODE) {
            // Add to open_intervals
            forest_state.open_intervals.emplace_back(std::move(current_interval));
        }
    }

    //Now close anything that remained open
    while(!forest_state.open_intervals.empty()) {
        interval_and_orientation_t& ancestor_interval = forest_state.open_intervals.back();
        const Seed& last_seed = seeds->at(forest_state.seed_sort_order[ancestor_interval.interval_end-1]);

        if (ancestor_interval.code_type == ZipCode::CHAIN ||
            ancestor_interval.code_type == ZipCode::ROOT_CHAIN ||
            ancestor_interval.code_type == ZipCode::ROOT_NODE) {
            //Close a chain

            close_chain(forest_state, distance_index, distance_limit, forest_state.open_intervals.size()-1, 
                        last_seed, ancestor_interval.is_reversed); 
        } else {
#ifdef DEBUG_ZIP_CODE_TREE
            assert(ancestor_interval.code_type == ZipCode::REGULAR_SNARL ||
                   ancestor_interval.code_type == ZipCode::IRREGULAR_SNARL ||
                   ancestor_interval.code_type == ZipCode::CYCLIC_SNARL ||
                   ancestor_interval.code_type == ZipCode::ROOT_SNARL);
#endif
            //Close a snarl
            close_snarl(forest_state, distance_index, forest_state.open_intervals.size()-1, 
                        last_seed, ancestor_interval.is_reversed, ancestor_interval.code_type == ZipCode::CYCLIC_SNARL); 
        }

        forest_state.open_intervals.pop_back();
    }

    if (trees[forest_state.active_zip_tree].zip_code_tree.size() == 0) {
        trees.erase(trees.begin() + forest_state.active_zip_tree);
    }
#ifdef DEBUG_ZIP_CODE_TREE
    print_self();
    validate_zip_forest(distance_index, distance_limit);
    assert(forest_state.open_chains.empty());
    assert(forest_state.open_intervals.empty());
#endif

}

void ZipCodeForest::open_chain(forest_growing_state_t& forest_state, const SnarlDistanceIndex& distance_index,
                       const size_t& distance_limit, const size_t& depth, const Seed& current_seed, bool chain_is_reversed) {
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
        if (forest_state.active_zip_tree == std::numeric_limits<size_t>::max() 
            || trees[forest_state.active_zip_tree].zip_code_tree.size() != 0) {
            //Don't add a new tree if the current one is empty
            trees.emplace_back(seeds);
            forest_state.active_zip_tree = trees.size()-1;
        }
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
                chain_is_reversed != is_rev(current_seed.pos)
                    ? current_seed.zipcode_decoder->get_length(depth) - offset(current_seed.pos)
                    : offset(current_seed.pos);
        } else {
            //Otherwise, this is really a chain, so get the prefix sum in the chain

            forest_state.sibling_indices_at_depth[depth-1].back().distances.first = chain_is_reversed 
                ? SnarlDistanceIndex::minus(current_seed.zipcode_decoder->get_length(depth) ,
                                            SnarlDistanceIndex::sum(
                                                current_seed.zipcode_decoder->get_offset_in_chain(depth+1),
                                                current_seed.zipcode_decoder->get_length(depth+1))) 
                : current_seed.zipcode_decoder->get_offset_in_chain(depth+1);

            if (depth+1 == current_max_depth) {
                //If this is a node, then add the offset of the position in the node
                bool child_is_reversed = ZipCodeTree::seed_is_reversed_at_depth(current_seed, depth+1, distance_index) 
                    ? !chain_is_reversed : chain_is_reversed;
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
                                               forest_state.sibling_indices_at_depth[depth-1].back().distances.first > distance_limit);
    }
}

void ZipCodeForest::close_chain(forest_growing_state_t& forest_state, const SnarlDistanceIndex& distance_index,
                       const size_t& distance_limit, const size_t& depth, const Seed& last_seed, bool chain_is_reversed) {

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
        while (trees[forest_state.active_zip_tree].zip_code_tree.size() > 0 &&
               trees[forest_state.active_zip_tree].zip_code_tree.back().type == ZipCodeTree::EDGE) {
            trees[forest_state.active_zip_tree].zip_code_tree.pop_back();
        }

        //Forget about the chain
        if (depth != 0) {
            forest_state.open_chains.pop_back();
        }

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
                    if (forest_state.open_chains.back().first != 0) {

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
                        while (trees[forest_state.active_zip_tree].zip_code_tree.size() > 0
                                && trees[forest_state.active_zip_tree].zip_code_tree.back().type == ZipCodeTree::EDGE) {
                            trees[forest_state.active_zip_tree].zip_code_tree.pop_back();
                        }
                    }
#ifdef DEBUG_ZIP_CODE_TREE
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
                    size_t last_length = depth == last_seed.zipcode_decoder->max_depth() ? 0 
                                                                                         : last_seed.zipcode_decoder->get_length(depth+1);

                    distance_to_chain_end = SnarlDistanceIndex::sum(distance_to_chain_end, 
                                            SnarlDistanceIndex::sum(last_edge,
                                                                    last_length));
                }
            }
            if (add_distances) {
                // If this chain (or chain slice) remains in the snarl, then add the distances
                // in the snarl

                //remember the distance to the end to be used in snarl distances
                forest_state.sibling_indices_at_depth[depth-1].back().distances.second = distance_to_chain_end;

                bool snarl_is_reversed = forest_state.open_intervals[forest_state.open_intervals.size()-2].is_reversed;
                bool is_cyclic_snarl =  forest_state.open_intervals[forest_state.open_intervals.size()-2].code_type == ZipCode::CYCLIC_SNARL;

                add_snarl_distances(forest_state, distance_index, depth-1, last_seed, chain_is_reversed, snarl_is_reversed, 
                                    false, is_cyclic_snarl);
            }
            //We've closed a chain, so take out the latest open chain
            forest_state.open_chains.pop_back();
        }
    }
}

void ZipCodeForest::add_child_to_chain(forest_growing_state_t& forest_state, const SnarlDistanceIndex& distance_index,
                       const size_t& distance_limit, const size_t& depth, const size_t& seed_index, bool child_is_reversed,
                       bool chain_is_reversed) {
    const Seed& current_seed = seeds->at(seed_index);
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
            child_is_reversed != is_rev(current_seed.pos)
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

            if (forest_state.active_zip_tree == std::numeric_limits<size_t>::max() 
                || trees[forest_state.active_zip_tree].zip_code_tree.size() != 0) {
                //Add a new tree and make sure it is the new active tree
                trees.emplace_back(seeds);
                forest_state.active_zip_tree = trees.size()-1;
            }

            //Add the start of the new chain
            trees[forest_state.active_zip_tree].zip_code_tree.push_back({ZipCodeTree::CHAIN_START, 
                                                                         std::numeric_limits<size_t>::max(), 
                                                                         false});

            //The first sibling in the chain is now the chain start, not the previous seed, so replace it
            forest_state.sibling_indices_at_depth[chain_depth].pop_back();
            forest_state.sibling_indices_at_depth[chain_depth].push_back({ZipCodeTree::CHAIN_START, 0}); 

        } else if (distance_between > distance_limit) { 
            //If this is too far from the previous thing

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
        trees[forest_state.active_zip_tree].zip_code_tree.push_back({ZipCodeTree::SEED, seed_index, child_is_reversed != is_rev(current_seed.pos)});
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
                                const size_t& depth, const Seed& last_seed, bool last_is_reversed, bool is_cyclic_snarl) {
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
        while (trees[forest_state.active_zip_tree].zip_code_tree.size() > 0
                && trees[forest_state.active_zip_tree].zip_code_tree.back().type == ZipCodeTree::EDGE) {
            trees[forest_state.active_zip_tree].zip_code_tree.pop_back();
        }
#ifdef DEBUG_ZIP_CODE_TREE
        assert(trees[forest_state.active_zip_tree].zip_code_tree.back().type == ZipCodeTree::SNARL_START);
#endif        
        //Pop the snarl start out
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

            if (depth > 0 && forest_state.open_chains.size() > 0 
                && forest_state.open_chains.back().first >= trees[forest_state.active_zip_tree].zip_code_tree.size()) {
                //If there was a chain slice that could have started at this snarl
#ifdef DEBUG_ZIP_CODE_TREE
                assert(forest_state.open_chains.back().second);
#endif
                //Find the start of the previous child
                size_t previous_index = trees[forest_state.active_zip_tree].zip_code_tree.size() - 1;
                bool found_sibling = false;
                bool opened_snarl = false;
                while (!found_sibling) {
                    if (!opened_snarl && trees[forest_state.active_zip_tree].zip_code_tree.at(previous_index).type == ZipCodeTree::SEED) {
                        found_sibling = true;
                    } else if (trees[forest_state.active_zip_tree].zip_code_tree.at(previous_index).type == ZipCodeTree::SNARL_END) {
                        opened_snarl = true;
                        previous_index--;
                    } else if ((trees[forest_state.active_zip_tree].zip_code_tree.at(previous_index).type == ZipCodeTree::SNARL_START)) {
                        found_sibling = true;
                    } else {
                        previous_index--;
                    }
                }
                if (trees[forest_state.active_zip_tree].zip_code_tree.at(previous_index-1).type == ZipCodeTree::CHAIN_START) {
                    previous_index--;
                }
#ifdef DEBUG_ZIP_CODE_TREE
                assert(( trees[forest_state.active_zip_tree].zip_code_tree.at(previous_index).type == ZipCodeTree::SEED ||
                         trees[forest_state.active_zip_tree].zip_code_tree.at(previous_index).type == ZipCodeTree::SNARL_START ||
                         trees[forest_state.active_zip_tree].zip_code_tree.at(previous_index).type == ZipCodeTree::CHAIN_START));
                cerr << "New start of previous open chain: " << previous_index << endl;;
#endif
                forest_state.open_chains.back().first = previous_index;
                
            }
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
       
        add_snarl_distances(forest_state, distance_index, depth, last_seed, last_is_reversed, last_is_reversed, true,
                            is_cyclic_snarl);

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
                                const size_t& depth, const Seed& seed, bool child_is_reversed, bool snarl_is_reversed, 
                                bool to_snarl_end, bool is_cyclic_snarl) {

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
    
    
    // If this is to the end bound, get the distance to all siblings. If it is to the last child, don't get
    // the distance to itself
    size_t sibling_count = to_snarl_end ? forest_state.sibling_indices_at_depth[depth].size() 
                                        : forest_state.sibling_indices_at_depth[depth].size()-1;
    for ( size_t sibling_i = 0 ; sibling_i < sibling_count ; sibling_i++) {
        const auto& sibling = forest_state.sibling_indices_at_depth[depth][sibling_i];

        if (sibling.type == ZipCodeTree::SNARL_START && !is_cyclic_snarl) {
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

                if (to_snarl_end && !is_cyclic_snarl) {

                    distance = SnarlDistanceIndex::sum( sibling.distances.second,
                                                            snarl_is_reversed ? sibling_seed.zipcode_decoder->get_distance_to_snarl_start(depth+1)
                                                                              : sibling_seed.zipcode_decoder->get_distance_to_snarl_end(depth+1));
                } else {

                    //If to_snarl_end is true, then we want the distance to the end (or start if snarl_is_reversed)
                    // Rank is 0 and the orientation doesn't matter
                    size_t rank2 = to_snarl_end ? (snarl_is_reversed ? 0 : 1) 
                                                : seed.zipcode_decoder->get_rank_in_snarl(depth+1);
                    bool right_side2 = child_is_reversed;

                    //If the sibling is the start, then get the distance to the appropriate bound
                    size_t rank1 = sibling.type == ZipCodeTree::SNARL_START ? (snarl_is_reversed ? 1 : 0)
                                                                            : sibling_seed.zipcode_decoder->get_rank_in_snarl(depth+1);
                    bool right_side1 = !sibling.is_reversed;

                    size_t distance_to_end_of_last_child = sibling.type == ZipCodeTree::SNARL_START ? 0
                                                     : sibling.distances.second;
                    //TODO: idk about this distance- I think the orientations need to change
                    //The bools for this are true if the distance is to/from the right side of the child
                    //We want the right side of 1 (which comes first in the dag ordering) to the left side of 2
                    //relative to the orientation of the snarl
                    net_handle_t snarl_handle = seed.zipcode_decoder->get_net_handle(depth, &distance_index);
                    distance = SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(
                        distance_index.distance_in_snarl(snarl_handle, rank1, right_side1, rank2, right_side2),
                        distance_to_chain_start),
                        distance_to_end_of_last_child);
                }
            }
            trees[forest_state.active_zip_tree].zip_code_tree.at(last_child_index - 1 - sibling_i) = {ZipCodeTree::EDGE, distance, false};
        }
    
    }

    //Remember the orientation of this child for the next time it gets used
    forest_state.sibling_indices_at_depth[depth].back().is_reversed = child_is_reversed;
}

std::pair<size_t, size_t> ZipCodeTree::dag_and_non_dag_snarl_count(const vector<Seed>& seeds, const SnarlDistanceIndex& distance_index) const {
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
                           seeds[current_item.value].zipcode_decoder->get_code_type(snarl_depth) == ZipCode::CYCLIC_SNARL ||
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
            cerr << seeds->at(item.value).pos << "/" << seeds->at(item.value).source;
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

bool ZipCodeTree::node_is_invalid(nid_t id, const SnarlDistanceIndex& distance_index, size_t distance_limit) const {
    bool is_invalid = false;
    net_handle_t net = distance_index.get_node_net_handle(id);
    while (!distance_index.is_root(net)) {
        if (distance_index.is_multicomponent_chain(net) || distance_index.is_looping_chain(net)) {
            //If this is something that we haven't handled
            is_invalid = true;
            break;
        } else if (distance_index.is_chain(distance_index.get_parent(net)) && 
                    !distance_index.is_trivial_chain(distance_index.get_parent(net))) {
            //Check if this net_handle_t could be involved in a chain loop that is smaller than the distance limit
            size_t forward_loop = distance_index.is_node(net) ? distance_index.get_forward_loop_value(net)
                                : distance_index.get_forward_loop_value(distance_index.get_node_from_sentinel(distance_index.get_bound(net, true, false)));
            size_t reverse_loop = distance_index.is_node(net) ? distance_index.get_reverse_loop_value(net)
                                : distance_index.get_reverse_loop_value(distance_index.get_node_from_sentinel(distance_index.get_bound(net, false, false)));
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

void ZipCodeTree::validate_zip_tree(const SnarlDistanceIndex& distance_index, size_t distance_limit) const {
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "Validate tree with distance limit " << distance_limit << endl;
#endif

    assert(zip_code_tree.size() != 0);

    /**********  Make sure that all snarls/chains are opened and closed in a valid order ****************/
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

    /************ Make sure that everything is in a valid order ****************/
    size_t previous_seed_index = std::numeric_limits<size_t>::max();
    bool previous_is_invalid = false;
    for (size_t i = 0 ; i < zip_code_tree.size() ; i++) {
        const tree_item_t& current_item = zip_code_tree[i];
        if (current_item.type == SEED) {
            //Check if this is worth validating
            //Use a distance limit of 0 so it will ignore looping chains
            bool current_is_invalid = node_is_invalid(id(seeds->at(current_item.value).pos), distance_index, 0);
            bool current_is_in_cyclic_snarl = node_is_in_cyclic_snarl(id(seeds->at(current_item.value).pos), distance_index);

            if (previous_seed_index != std::numeric_limits<size_t>::max() &&
                !current_is_invalid && !previous_is_invalid) {
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
                                   ? seeds->at(previous_seed_index).zipcode_decoder->get_length(depth) - offset(seeds->at(previous_seed_index).pos)
                                   : offset(seeds->at(previous_seed_index).pos);
                    size_t offset2 = is_rev(seeds->at(current_item.value).pos)
                                   ? seeds->at(current_item.value).zipcode_decoder->get_length(depth) - offset(seeds->at(current_item.value).pos)
                                   : offset(seeds->at(current_item.value).pos);
                    if (!current_is_in_cyclic_snarl) {
                        if (!a_is_reversed) {
                            //If they are in previous_seed_index snarl or they are facing forward on a chain, then order by
                            //the offset in the node
                            assert( offset1 <= offset2);
                        } else {
                            //Otherwise, the node is facing backwards in the chain, so order backwards in node
                            assert( offset2 <= offset1);
                        }
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
                    if (!current_is_in_cyclic_snarl) {

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
                    }
                } else if (seeds->at(previous_seed_index).zipcode_decoder->get_code_type(depth-1) == ZipCode::REGULAR_SNARL
                            || seeds->at(previous_seed_index).zipcode_decoder->get_code_type(depth-1) == ZipCode::IRREGULAR_SNARL) {
#ifdef DEBUG_ZIP_CODE_TREE
                    cerr << "\t they are children of a common dag snarl" << endl;
#endif
                    // Otherwise, they are children of a snarl
                    // Sort by a topological ordering from the start of the snarl
                    // The ranks of children in snarls are in a topological order, so 
                    // sort on the ranks
                    if (!current_is_in_cyclic_snarl) {
                        assert( seeds->at(previous_seed_index).zipcode_decoder->get_rank_in_snarl(depth) <=
                                seeds->at(current_item.value).zipcode_decoder->get_rank_in_snarl(depth));
                    }
                } 

            }
            previous_seed_index = current_item.value;
            previous_is_invalid = current_is_invalid;
        } else if (current_item.type == CHAIN_START) {
            //Chains can't start with edges
            assert(zip_code_tree[i+1].type != EDGE);
        } else if (current_item.type == CHAIN_END) {
            //And can't end with edges
            assert(zip_code_tree[i-1].type != EDGE);
        }
    }



    /************* Check distances and snarl tree relationships *******************/

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

        //For cyclic snarls, the tree distance isn't always guaranteed to be the same as the minimum distance
        // I think that the smallest distance between any pair of seeds will be guaranteed to be the same as the
        // actual minimum distance, so store the minimum (non infinite) distance here
        // The first pair of size_t's are indices into seeds (start then next), 
        // the second pair are the tree distance and actual distance

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

            size_t index_distance = distance_index.minimum_distance(id(next_seed.pos), is_rev(next_seed.pos), offset(next_seed.pos),
                    id(start_seed.pos), is_rev(start_seed.pos), offset(start_seed.pos), true);
            if (index_distance != std::numeric_limits<size_t>::max() && is_rev(next_seed.pos) != next_is_reversed) {
                //If the seed we're starting from got reversed, then subtract 1
                index_distance -= 1;
            }
            if (index_distance != std::numeric_limits<size_t>::max() && is_rev(start_seed.pos) != start_is_reversed) {
                //If the seed we ended at got reversed, then add 1
                index_distance += 1;
            }
            pos_t start_pos = is_rev(start_seed.pos) ? make_pos_t(id(start_seed.pos), false, distance_index.minimum_length(start_handle) - offset(start_seed.pos) )
                                                     : start_seed.pos;
            pos_t next_pos = is_rev(next_seed.pos) ? make_pos_t(id(next_seed.pos), false, distance_index.minimum_length(next_handle) - offset(next_seed.pos) )
                                                     : next_seed.pos;
            size_t start_length = distance_index.minimum_length(start_handle);
            size_t next_length = distance_index.minimum_length(next_handle);

            bool in_non_dag_snarl = node_is_in_cyclic_snarl(id(next_seed.pos), distance_index) ||
                                    node_is_in_cyclic_snarl(id(start_seed.pos), distance_index);
            bool distance_is_invalid = node_is_invalid(id(next_seed.pos), distance_index, distance_limit) ||
                               node_is_invalid(id(start_seed.pos), distance_index, distance_limit);
            if (in_non_dag_snarl) {
                //TODO: I don't actually know how to check these properly

            } else if (!distance_is_invalid && index_distance <= distance_limit) {
                if (start_pos == next_pos) {
                    if (tree_distance != 0 && tree_distance != index_distance) {
                        for (auto& seed : *seeds) {
                            cerr << seed.pos << endl;
                        }
                        cerr << "Distance between " << next_seed.pos << (next_is_reversed ? "rev" : "") << " and " << start_seed.pos << (start_is_reversed ? "rev" : "") << endl;
                        cerr << "Forward positions: " << start_pos << " " << next_pos << " and length " << start_length << endl;
                        cerr << "Tree distance: " << tree_distance << " index distance: " << index_distance << endl;
                        cerr << "With distance limit: " << distance_limit << endl;
                    }
                    //This could be off by one if one of the seeds is reversed, but I'm being lazy and just checking against the index
                    assert((tree_distance == 0 || tree_distance == index_distance));
                } else {
                    if (tree_distance != index_distance) {
                        for (auto& seed : *seeds) {
                            cerr << seed.pos << endl;
                        }
                        cerr << "Distance between " << next_seed.pos << (next_is_reversed ? "rev" : "") << " and " << start_seed.pos << (start_is_reversed ? "rev" : "") << endl;
                        cerr << "Forward positions: " << start_pos << " " << next_pos << " and lengths " << start_length << " " << next_length << endl;
                        cerr << "Tree distance: " << tree_distance << " index distance: " << index_distance << endl;
                        cerr << "With distance limit: " << distance_limit << endl;
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
                // Note that there may be 0 stored distances if we are below the top-level snarl.
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
            // the snarl. We should have stacked exactly one or zero distances.
            
            if (depth() == 1) {
                // We have hit the start of a top-level snarl
#ifdef debug_parse
                std::cerr << "Hit start of top-level snarl" << std::endl;
#endif
                halt();
                // When we halt we have to return true to show the halting position.
                return true;
            }

            // Throw out parent running distance
            pop();

            // There will be a running distance on the stack still, and we
            // will continue with that in the parent chain.
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
                // Parent snarl may be a top-level snarl.
                if (depth() == 1) {
                    // We have hit the start of a top-level snarl
#ifdef debug_parse
                    std::cerr << "Hit start of top-level snarl" << std::endl;
#endif
                    halt();
                    // When we halt we have to return true to show the halting position.
                    return true;
                }

                // This is the start of the chain we were wanting to skip.
                pop();
                crash_unless(depth() >= 1);
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

vector<ZipCodeForest::interval_and_orientation_t> ZipCodeForest::sort_one_interval(vector<size_t>& zipcode_sort_order,
    const interval_and_orientation_t& interval, size_t interval_depth, const SnarlDistanceIndex& distance_index) const {

    /*
      Sort the seeds in roughly linear/topological-ish order along the top-level chains

      Sorting is split into two different types of sort: radix sort or an n-log-n sort,
      depending on which will be more efficient
    */

    //Helper function to get the value to sort on from the zipcode
    //This doesn't take into account the orientation, except for nodes offsets in chains
    //It will actually be defined somewhere else
    //Used for sorting at the given depth, so use values at depth depth+1
    auto get_sort_value = [&] (const Seed& seed, size_t depth) {
#ifdef DEBUG_ZIP_CODE_SORTING
        cerr << "\tGet the sort value of seed " << seed.pos << " at depth " << depth << endl;
#endif
        ZipCode::code_type_t code_type = seed.zipcode_decoder->get_code_type(depth);
        if (code_type == ZipCode::NODE || code_type == ZipCode::ROOT_NODE || seed.zipcode_decoder->max_depth() == depth) {
#ifdef DEBUG_ZIP_CODE_SORTING
            cerr << "\t\t this is a node: offset: " << ( is_rev(seed.pos) ? seed.zipcode_decoder->get_length(depth) - offset(seed.pos)
                                    : offset(seed.pos)) << endl;;
#endif
            return is_rev(seed.pos) ? seed.zipcode_decoder->get_length(depth) - offset(seed.pos)
                                    : offset(seed.pos);
        } else if (code_type == ZipCode::CHAIN || code_type == ZipCode::ROOT_CHAIN) {
#ifdef DEBUG_ZIP_CODE_SORTING
            cerr << "\t\t this is a chain: prefix sum value x2 (and -1 if snarl): ";
#endif
            //Return the prefix sum in the chain
            //Since the offset stored represents the space between nucleotides, two positions on different nodes
            // could have the same offset. Similarly, a snarl could have the same prefix sum as a node.
            // For example, in this graph:
            //                2
            //               [AA]  
            //           1  /   \  3
            //          [AA] --- [AA]
            // The positions n1-0 and 3+0, and the snarl 1-3 all have the same offset of 2
            // To solve this, the prefix sum of a chain will always be multiplied by 3, and 1 will be added to snarls,
            // And 2 will be added to the node with an offset in the node of 0 (node 3 if the chain is traversed forward)

            size_t prefix_sum;
            if (seed.zipcode_decoder->get_code_type(depth+1) == ZipCode::REGULAR_SNARL 
                || seed.zipcode_decoder->get_code_type(depth+1) == ZipCode::IRREGULAR_SNARL
                || seed.zipcode_decoder->get_code_type(depth+1) == ZipCode::CYCLIC_SNARL) { 
                //If this is a snarl, then get the prefix sum value*3 + 1
                prefix_sum = SnarlDistanceIndex::sum(seed.zipcode_decoder->get_offset_in_chain(depth+1) * 3,  1);
            } else {
                //If this is a node, then get the prefix sum value plus the offset in the position, and multiply by 2 
                size_t node_offset = seed.zipcode_decoder->get_is_reversed_in_parent(depth+1) != is_rev(seed.pos)
                                                     ? seed.zipcode_decoder->get_length(depth+1) - offset(seed.pos)
                                                     : offset(seed.pos);
                prefix_sum = SnarlDistanceIndex::sum(seed.zipcode_decoder->get_offset_in_chain(depth+1), node_offset);
                prefix_sum *= 3;
                if (node_offset == 0) {
                    prefix_sum = SnarlDistanceIndex::sum(prefix_sum, 2);
                }
            }
#ifdef DEBUG_ZIP_CODE_SORTING
            cerr << prefix_sum << endl;
#endif
            return prefix_sum;
        } else {
#ifdef DEBUG_ZIP_CODE_SORTING
            cerr << "\tThis is snarl, so return the rank in the snarl: " << seed.zipcode_decoder->get_rank_in_snarl(depth+1) << endl;
#endif
            // The ranks of children in irregular snarls are in a topological order, so 
            // sort on the ranks
            // The rank of children in a regular snarl is arbitrary but it doesn't matter anyway
            return seed.zipcode_decoder->get_rank_in_snarl(depth+1);
        }
    };

    //Returns true if a cyclic snarl gets traversed end-to-start in a forward traversal of the read
    //Assumes that the parent of the cyclic snarl has been sorted
    //The snarl has seed indices in forest_state.seed_sort_order [start_of_snarl, end_of_snarl)
    auto cyclic_snarl_is_traversed_backwards = [&] (const interval_and_orientation_t& chain_interval,
                                    size_t start_of_snarl, size_t end_of_snarl, const vector<size_t>& sort_order) {
        //If this is a cyclic snarl, then check if it is being traversed forward or backward by the read
        // Take a sample of seeds before and after the snarl to get the direction

        //This contains read offsets from before the snarl (or from the snarl if there was nothing before it in its parent)
        vector<size_t> preceding_offsets;
        if (start_of_snarl == chain_interval.interval_start) {
            //If this is the first interval of the chain, then just take stuff from the snarl
            for (int check_i = start_of_snarl ; check_i < end_of_snarl && check_i - start_of_snarl < 3; check_i++) {
                preceding_offsets.emplace_back(seeds->at(sort_order[check_i]).source);
            }
        } else {
            //Otherwise, take seeds from before the snarl in the chain
            for (int check_i = start_of_snarl-1 ; check_i >= chain_interval.interval_start && start_of_snarl - check_i <= 3; check_i--) {
                preceding_offsets.emplace_back(seeds->at(sort_order[check_i]).source);
            }
        }

        //This contains read offsets from after the snarl
        vector<size_t> succeeding_offsets;
        if (end_of_snarl == chain_interval.interval_end) {
            //If there is nothing after, take from the snarl
            for (int check_i = start_of_snarl ; check_i < end_of_snarl && check_i - start_of_snarl < 3; check_i++) {
                succeeding_offsets.emplace_back(seeds->at(sort_order[check_i]).source);
            }
        } else {
            //Otherwise, take from whatever comes next in the chain
            for (int check_i = end_of_snarl ; check_i < chain_interval.interval_end && check_i < end_of_snarl+3 ; check_i++) {
                succeeding_offsets.emplace_back(seeds->at(sort_order[check_i]).source);
            }
        }
#ifdef DEBUG_ZIP_CODE_SORTING
        if (preceding_offsets.size() == 0 || succeeding_offsets.size() == 0) {
            //If there is nothing to judge by, just say it isn't reversed
            return false;
            //TODO: I don't think this will happen. If there is nothing before or after, it will fill both in with the snarl
            assert(false);
        }
#endif 
        //Take the median of each vector and see which is greater
        std::sort(preceding_offsets.begin(), preceding_offsets.end());
        size_t median_preceding = preceding_offsets[ preceding_offsets.size() / 2];

        std::sort(succeeding_offsets.begin(), succeeding_offsets.end());
        size_t median_succeeding = succeeding_offsets[ succeeding_offsets.size() / 2];
#ifdef DEBUG_ZIP_CODE_SORTING
        cerr << "Preceeding: ";
        for (auto x : preceding_offsets) { cerr << x << " ";}
        cerr << endl << "Succeeding ";
        for (auto x : succeeding_offsets) {cerr << x << " ";}
        cerr << endl;
#endif

        return median_preceding > median_succeeding;

    };

    //At the given depth, go through sort_order in the given interval to find the intervals for the next level 
    //and add to new_intervals
    auto find_next_intervals = [&] (const interval_and_orientation_t& interval,
                                    size_t depth, const vector<size_t>& sort_order, 
                                    const std::function<size_t(const Seed& seed, size_t depth)>& get_partitioning_value) {
#ifdef DEBUG_ZIP_CODE_TREE
        cerr << "Finding intervals after sorting at depth " << depth << endl;
#endif
        vector<interval_and_orientation_t> new_intervals;
        //After sorting, find runs of equivalent values for new_interval_to_sort
        //Everything gets put into a new interval, even if it is the only thing with that partitioning value
        //Since nodes are really just seeds on the same chain, runs of nodes get put together even if they are
        // actually on different nodes, as long as the nodes are facing in the same direction
        //Also need to check the orientation
        //For intervals corresponding to cyclic snarls, the orientation is based on the read, not the snarl

        if (seeds->at(sort_order[interval.interval_start]).zipcode_decoder->max_depth() == depth ) {
            //If this is a trivial chain, then just return the same interval as a node
            new_intervals.emplace_back(interval.interval_start, interval.interval_end, interval.is_reversed, ZipCode::NODE);
            return new_intervals;
        }


        //These get compared to see if the next seeds is in the same interval
        ZipCode::code_type_t first_type = seeds->at(sort_order[interval.interval_start]).zipcode_decoder->get_code_type(depth+1);

        //This is only for nodes in chains, since anything on nodes in chains are considered just children of the chain
        bool previous_is_node = first_type == ZipCode::NODE;

        //This only matters if it isn't a node
        size_t previous_sort_value = previous_is_node 
                                   ? (ZipCodeTree::seed_is_reversed_at_depth(seeds->at(sort_order[interval.interval_start]), depth+1, distance_index) ? 1 : 0)
                                   : get_partitioning_value(seeds->at(sort_order[interval.interval_start]), depth);

        //Start the first interval. The end value and is_reversed gets set when ending the interval
        new_intervals.emplace_back(interval.interval_start, interval.interval_start, interval.is_reversed, 
                                   previous_is_node ? ZipCode::NODE : first_type);
        for (size_t i = interval.interval_start+1 ; i < interval.interval_end ; i++) {
            
            //If the current seed is a node and has nothing at depth+1 or is different from the previous seed at this depth
            ZipCode::code_type_t current_type = seeds->at(sort_order[i]).zipcode_decoder->get_code_type(depth+1);
            bool is_node = current_type == ZipCode::NODE;
            size_t sort_value = is_node 
                              ? (ZipCodeTree::seed_is_reversed_at_depth(seeds->at(sort_order[i]), depth+1, distance_index) ? 1 : 0)
                              : get_partitioning_value(seeds->at(sort_order[i]), depth);
            bool is_different_from_previous = is_node != previous_is_node ? true : sort_value != previous_sort_value;
            previous_is_node = is_node;
            previous_sort_value = sort_value;

            if (is_different_from_previous) {
                //If this is the end of a run, close the previous run
                //Add its end value and orientation

                new_intervals.back().interval_end = i;

                if (new_intervals.back().code_type == ZipCode::CYCLIC_SNARL) {
                    new_intervals.back().is_reversed = cyclic_snarl_is_traversed_backwards(interval,
                                    new_intervals.back().interval_start, i, sort_order); 
                } else {
                    new_intervals.back().is_reversed = ZipCodeTree::seed_is_reversed_at_depth(seeds->at(sort_order[i-1]), depth+1, distance_index) 
                                         ? !interval.is_reversed
                                         : interval.is_reversed;
                }
 

                //Open a new run
                new_intervals.emplace_back(i, i, interval.is_reversed, is_node ? ZipCode::NODE : current_type);
            }
        }

        //Close the last run
        new_intervals.back().interval_end = interval.interval_end;

        //Get the orientation of the previous child
        if (new_intervals.back().code_type == ZipCode::CYCLIC_SNARL) {
            //For a cyclic snarl
            new_intervals.back().is_reversed = cyclic_snarl_is_traversed_backwards(interval,
                            new_intervals.back().interval_start, interval.interval_end, sort_order); 
        } else {
            new_intervals.back().is_reversed = ZipCodeTree::seed_is_reversed_at_depth(seeds->at(sort_order[interval.interval_end-1]), depth+1, distance_index) 
                                 ? !interval.is_reversed
                                 : interval.is_reversed;
        }
#ifdef DEBUG_ZIP_CODE_SORTING
        cerr << "New sort order " << endl;
        for (auto& interval : new_intervals) {
            for (size_t i = interval.interval_start ; i < interval.interval_end ; i++) {
                cerr << seeds->at(i).pos << ", ";
            }
            cerr << "|";
        }
        cerr << endl;
#endif
        return new_intervals;
    };

    if (interval.code_type == ZipCode::EMPTY) {

        // If we are sorting the root int connected components 
        // Assume that the number of connected components is small enough that radix sort is more efficient
        radix_sort_zipcodes(zipcode_sort_order, interval, 
                            false, std::numeric_limits<size_t>::max(), distance_index,
                            [&](const Seed& seed, size_t depth) {
                                //Sort on the connected component number
                                return seed.zipcode_decoder->get_distance_index_address(0);
                            });

#ifdef DEBUG_ZIP_CODE_SORTING
        cerr << "After root " << endl;
        for (size_t i : zipcode_sort_order) {
            cerr << i << ":" << seeds->at(i).pos << ", ";
        }
        cerr << endl;
#endif
        return find_next_intervals(interval, std::numeric_limits<size_t>::max(), zipcode_sort_order,
                            [&](const Seed& seed, size_t depth) {
                                //Sort on the connected component number
                                return seed.zipcode_decoder->get_distance_index_address(0);
                            });
    } else {

#ifdef DEBUG_ZIP_CODE_SORTING
        cerr << "Sort seeds on interval " << interval.interval_start << "-" << interval.interval_end << " at depth " << interval_depth << endl;
#endif


        //One of the seeds getting sorted
        const Seed& seed_to_sort = seeds->at(zipcode_sort_order[interval.interval_start]);



        if (interval.code_type == ZipCode::CYCLIC_SNARL) {
            // If this is a cyclic snarl, then the children should be sorted by both their position on the graph
            // and their offset on the read

            //First, figure out if the read flows through the snarl start-to-end or end-to-start

            //Sort the snarl and get intervals of the snarl's children
            auto new_intervals =  sort_zipcodes_on_cyclic_snarl(zipcode_sort_order, interval, interval_depth, distance_index);
            if (new_intervals.size() != 0) {
                return new_intervals;
            }
            //If finding intervals on the cyclic snarl failed, then keep going as if it wasn't cyclic
        }
        //If this either wasn't a cyclic snarl or it was a cyclic snarl that failed

        // Sorting will either be done with radix sort or with std::sort, depending on which is more efficient
        // Radix sort is linear time in the number of items it is sorting, but also linear space in the range 
        // of the values it is sorting on
        // If the range of values is greater than the n log n (in the number of things being sorted) of the default
        // sorter, then use radix

        bool use_radix;
        if (interval.code_type  == ZipCode::ROOT_CHAIN) {
            //If this is a root chain, then use the default sort, because it's probably too big for radix and we can't tell
            //anyways because we don't store the length of a root-chain
            use_radix = false;
        } else if (interval.code_type == ZipCode::NODE || interval.code_type == ZipCode::CHAIN) {
            //If we're sorting a node or chain, then the range of values is the minimum length of the node/chain
            // times 3 because it gets multiplied by 3 to differentiate nodes and snarls
            size_t radix_cost = seed_to_sort.zipcode_decoder->get_length(interval_depth) * 3;
            size_t default_cost = (interval.interval_end - interval.interval_start) * std::log2(interval.interval_end - interval.interval_start);

            use_radix = radix_cost < default_cost;
        } else {
            //Otherwise, this is a snarl and the range of values is the number of children in the snarl

            size_t radix_cost = seed_to_sort.zipcode_decoder->get_snarl_child_count(interval_depth, &distance_index);
            size_t default_cost = (interval.interval_end - interval.interval_start) * std::log2(interval.interval_end - interval.interval_start);

            use_radix = radix_cost < default_cost;
        }
        bool reverse_order = (interval.code_type == ZipCode::REGULAR_SNARL || interval.code_type == ZipCode::IRREGULAR_SNARL) 
                             ? false
                             : interval.is_reversed; 
        //For everything except a cyclic snarl, sort normally
        if (use_radix) {
            //Sort the given interval using the value-getter and orientation
            radix_sort_zipcodes(zipcode_sort_order, interval, reverse_order, interval_depth, distance_index, get_sort_value);
        } else {
            //Sort the given interval using the value-getter and orientation
            default_sort_zipcodes(zipcode_sort_order, interval, reverse_order, interval_depth, distance_index, get_sort_value);
        }
        return find_next_intervals(interval, interval_depth, zipcode_sort_order, get_sort_value);
    }

}

void ZipCodeForest::radix_sort_zipcodes(vector<size_t>& zipcode_sort_order, const interval_and_orientation_t& interval,
                                   bool reverse_order, size_t depth, const SnarlDistanceIndex& distance_index, 
                                   const std::function<size_t(const Seed& seed, size_t depth)>& get_sort_value) const {
    //Radix sort the interval of zipcode_sort_order in the given interval
#ifdef DEBUG_ZIP_CODE_SORTING
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
                                   const std::function<size_t(const Seed& seed, size_t depth)>& get_sort_value) const { 
    //std::sort the interval of zipcode_sort_order between interval_start and interval_end
    
#ifdef DEBUG_ZIP_CODE_SORTING
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

vector<ZipCodeForest::interval_and_orientation_t> ZipCodeForest::sort_zipcodes_on_cyclic_snarl(vector<size_t>& zipcode_sort_order, const interval_and_orientation_t& interval,
                             size_t depth, const SnarlDistanceIndex& distance_index) const {
#ifdef DEBUG_ZIP_CODE_SORTING
    cerr << "Sort seeds on a cyclic snarl" << endl;
#endif

    /**** First, sort by the child that the seeds are on ****/
    
    size_t radix_cost = seeds->at(zipcode_sort_order[interval.interval_start]).zipcode_decoder->get_snarl_child_count(depth, &distance_index);
    size_t default_cost = (interval.interval_end - interval.interval_start) * std::log2(interval.interval_end - interval.interval_start);

    bool use_radix = radix_cost < default_cost;
    if (use_radix) {
        radix_sort_zipcodes(zipcode_sort_order, interval, interval.is_reversed, depth, distance_index, [&] (const Seed& seed, size_t depth) {
            return  seed.zipcode_decoder->get_rank_in_snarl(depth+1);
        });
    } else {
        default_sort_zipcodes(zipcode_sort_order, interval, interval.is_reversed, depth, distance_index, [&] (const Seed& seed, size_t depth) {
            return  seed.zipcode_decoder->get_rank_in_snarl(depth+1);
        });
    }
#ifdef DEBUG_ZIP_CODE_SORTING
    cerr << "Sorted order: ";
    for (size_t i = interval.interval_start ; i < interval.interval_end ; i++) {
        cerr << seeds->at(zipcode_sort_order[i]).pos << " ";
    }
    cerr << endl;
#endif

    /****Find the intervals of the children ****/

    vector<interval_and_orientation_t> child_intervals;


    child_intervals.emplace_back(interval.interval_start, interval.interval_start, false, ZipCode::CHAIN);
    for (size_t i = interval.interval_start+1 ; i < interval.interval_end ; i++) {
           
        const Seed& current_seed = seeds->at(zipcode_sort_order[i]);
        const Seed& previous_seed = seeds->at(zipcode_sort_order[i-1]);
        //Are the seeds on different children of the snarl? 
        bool is_different_from_previous = !ZipCodeDecoder::is_equal(*current_seed.zipcode_decoder, 
                                                                    *previous_seed.zipcode_decoder, depth+1);

        if (is_different_from_previous) {

            child_intervals.back().interval_end = i;

            child_intervals.emplace_back(i,  i, false, ZipCode::CHAIN);
        }
    }
    child_intervals.back().interval_end = interval.interval_end;

#ifdef DEBUG_ZIP_CODE_SORTING
    cerr << "Intervals of children" << endl;
    for (auto& interval : child_intervals) {
        for (size_t i = interval.interval_start ; i < interval.interval_end ; i++) {
            cerr << seeds->at(zipcode_sort_order[i]).pos << ", ";
        }
        cerr << "|";
    }
    cerr << endl;
#endif

    /**** For each child interval, sort the seeds by their offset in the read ****/

    //Remember the largest and smallest read offsets, so we can determine if its faster to do radix or nlogn sort
    size_t min_read_offset = seeds->at(zipcode_sort_order[interval.interval_start]).source;
    size_t max_read_offset = min_read_offset; 

    for (const interval_and_orientation_t& child_interval : child_intervals) {

        //First, which sort should we use?
        size_t radix_cost = max_read_offset - min_read_offset;
        size_t default_cost = (child_interval.interval_end - child_interval.interval_start) * 
                              std::log2(child_interval.interval_end - child_interval.interval_start);

        bool use_radix = radix_cost < default_cost;

        if (use_radix) {
            radix_sort_zipcodes(zipcode_sort_order, child_interval, 
                        interval.is_reversed, std::numeric_limits<size_t>::max(), distance_index,
                        [&](const Seed& seed, size_t depth) {
                            //Sort on the offset in the read
                            return seed.source;
                        });
        } else {
            default_sort_zipcodes(zipcode_sort_order, child_interval, 
                        interval.is_reversed, std::numeric_limits<size_t>::max(), distance_index,
                        [&](const Seed& seed, size_t depth) {
                            //Sort on the offset in the read
                            return seed.source;
                        });
        }
    }

#ifdef DEBUG_ZIP_CODE_SORTING
    cerr << "After sorting children" << endl;
    for (auto& interval : child_intervals) {
        for (size_t i = interval.interval_start ; i < interval.interval_end ; i++) {
            cerr << seeds->at(zipcode_sort_order[i]).pos << "/" << seeds->at(zipcode_sort_order[i]).source << ", ";
        }
        cerr << "|";
    }
    cerr << endl;
#endif

    /****** Find intervals along each child where the order of the read and the order in the chain disagree  *******/

    //Helper function to get the prefix sum of the child on the chain (child of the cyclic snarl). 
    //Used for ordering the children
    auto get_prefix_sum = [&] (const Seed& seed) {
        size_t prefix_sum;
        if (seed.zipcode_decoder->max_depth() == depth+1) {
            //If this is a node pretending to be a chain

            //Just use the offset in the node
            prefix_sum = seed.zipcode_decoder->get_is_reversed_in_parent(depth+1) != is_rev(seed.pos)
                            ? seed.zipcode_decoder->get_length(depth+1) - offset(seed.pos)
                            : offset(seed.pos);
        } else if (seed.zipcode_decoder->get_code_type(depth+2) == ZipCode::REGULAR_SNARL 
            || seed.zipcode_decoder->get_code_type(depth+2) == ZipCode::IRREGULAR_SNARL
            || seed.zipcode_decoder->get_code_type(depth+2) == ZipCode::CYCLIC_SNARL) { 
            //If this is a snarl, then get the prefix sum value*3 + 1
            prefix_sum = SnarlDistanceIndex::sum(seed.zipcode_decoder->get_offset_in_chain(depth+2) * 3,  1);
        } else {
            //If this is a node, then get the prefix sum value plus the offset in the position, and multiply by 2 
            size_t node_offset = seed.zipcode_decoder->get_is_reversed_in_parent(depth+2) != is_rev(seed.pos)
                             ? seed.zipcode_decoder->get_length(depth+2) - offset(seed.pos)
                             : offset(seed.pos);
            prefix_sum = SnarlDistanceIndex::sum(seed.zipcode_decoder->get_offset_in_chain(depth+2), node_offset);
            prefix_sum *= 3;
            if (node_offset == 0) {
                prefix_sum = SnarlDistanceIndex::sum(prefix_sum, 2);
            }
        }
        return prefix_sum;
    };

    vector<interval_and_orientation_t> read_intervals;
    for (const interval_and_orientation_t& child_interval : child_intervals) {
        //For each child interval, split into new intervals if the order in the read differs from the order in the graph
        // The read may go through the child forwards, backwards, both, multiple times, etc.
        // TODO: I don't know how to deal with this properly, so for now this will find slices of seeds that are monotonically
        // increasing or decreasing along the child. Any time it switches (the next seed came before the previous in the chain),
        // then start a new segment. The new segment's orientation will be determined by the seed after it.
        // This is very easily broken but it's the best I can think of for now

        enum orientation_t {FORWARD, BACKWARD, EQUAL};
        //At first, we don't know if the current run of seeds is going forwards or backwards in the child
        orientation_t current_orientation = EQUAL;

        //Start a new read_interval, initially just the start, add the end when starting a new one
        read_intervals.emplace_back(child_interval.interval_start, child_interval.interval_start, false, ZipCode::CHAIN);
        size_t previous_prefix_sum = get_prefix_sum(seeds->at(zipcode_sort_order[child_interval.interval_start]));
        for (size_t i = child_interval.interval_start+1 ; i < child_interval.interval_end ; i++) {
            size_t current_prefix_sum = get_prefix_sum(seeds->at(zipcode_sort_order[i]));
            if (current_orientation == EQUAL) {
                // If we don't know yet what orientation this run is in, this seed will still be added to the
                // current run, and we just need to check if there is a new orientation
                if (previous_prefix_sum < current_prefix_sum) {
                    current_orientation = FORWARD;
                } else if (previous_prefix_sum > current_prefix_sum) {
                    current_orientation = BACKWARD;
                }
                                       
            } else if ((current_orientation == FORWARD && previous_prefix_sum > current_prefix_sum) ||
                       (current_orientation == BACKWARD && previous_prefix_sum < current_prefix_sum)) {
                //If we are currently traversing in a specific direction and the next seed is
                // going in the opposite direction

                //End the current run
                read_intervals.back().interval_end = i;

                //If the child chain is traversed backwards in its own local orientation
                read_intervals.back().is_reversed = current_orientation == BACKWARD;

                //Start a new run
                read_intervals.emplace_back(i, i, false, ZipCode::CHAIN);

                //We don't yet know the orientation of the next run, so leave it at EQUAL
                current_orientation = EQUAL;
            }

            previous_prefix_sum = current_prefix_sum;
        }
        //Now end the last run
        read_intervals.back().interval_end = child_interval.interval_end;
        read_intervals.back().is_reversed = current_orientation == BACKWARD;
    }

    if (read_intervals.size() > 5*child_intervals.size()) {
        //If there are more than 5 duplicates per child chain
        vector<interval_and_orientation_t> empty;
        return empty;
    }

#ifdef DEBUG_ZIP_CODE_SORTING
    cerr << "After splitting/duplicating chains " << endl;
    for (auto& interval : read_intervals) {
        for (size_t i = interval.interval_start ; i < interval.interval_end ; i++) {
            cerr << seeds->at(zipcode_sort_order[i]).pos << ", ";
        }
        cerr << interval.is_reversed << "|";
    }
    cerr << endl;
#endif

    /*****  Find the sort order of the intervals, ordered by the first seed in the read *****/
    vector<size_t> interval_sort_order(read_intervals.size(), 0);
    for (size_t i = 0 ; i < interval_sort_order.size() ; i++) {
        interval_sort_order[i] = i;
    } 
    std::sort(interval_sort_order.begin(), interval_sort_order.end(), [&] (const size_t& a, const size_t& b) {
        // Sort by the first seed in the read in the interval. Since the intervals are sorted by read position, 
        // the first seed in the read will be the first seed
        size_t seed_index_a = read_intervals[a].interval_start;

        size_t seed_index_b = read_intervals[b].interval_start;

        if (interval.is_reversed) {
            return seeds->at(zipcode_sort_order[seed_index_a]).source >= seeds->at(zipcode_sort_order[seed_index_b]).source;
        } else {
            return seeds->at(zipcode_sort_order[seed_index_a]).source < seeds->at(zipcode_sort_order[seed_index_b]).source;
        }
    });

    /****** Get the sort order of the seeds, to be copied back into the real one ********/

    //The new sort order. Values are indices into seeds, so it will be copied directly into zipcode_sort_order
    vector<size_t> new_sort_order;
    new_sort_order.reserve(interval.interval_end - interval.interval_start);

    //Get the same intervals, but this time ordered and in terms of indices into zipcode_sort_order
    //The new order might put two seeds in the same chain next to each other, when they should really be different intervals
    vector<interval_and_orientation_t> new_intervals;

    for (size_t interval_i : interval_sort_order) {
        const auto& current_interval = read_intervals[interval_i];

        //Add this interval in terms of the actual zipcode_sort_order
        new_intervals.emplace_back(new_sort_order.size() + interval.interval_start, 
                                   new_sort_order.size() + interval.interval_start + (current_interval.interval_end - current_interval.interval_start), 
                                   current_interval.is_reversed,
                                   ZipCode::CHAIN);

        //Add everything in this interval to the new sort order
        for (size_t i = current_interval.interval_start ; i < current_interval.interval_end ; i++) {
            new_sort_order.emplace_back(zipcode_sort_order[i]);
        }
    }

    //Replace everything in the interval in zipcode_sort_order with the new sorted values
    for (size_t i = 0 ; i < new_sort_order.size() ; i++) {
        zipcode_sort_order[interval.interval_start + i] = new_sort_order[i];
    }
#ifdef DEBUG_ZIP_CODE_SORTING
    cerr << "New sort order for cyclic snarl" << endl;
    for (auto& interval : new_intervals) {
        for (size_t i = interval.interval_start ; i < interval.interval_end ; i++) {
            cerr << seeds->at(zipcode_sort_order[i]).pos << ", ";
        }
        cerr << "|";
    }
    cerr << endl;
#endif

    return new_intervals;
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

