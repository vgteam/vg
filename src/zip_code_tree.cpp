//#define DEBUG_ZIP_CODE_TREE
//#define PRINT_NON_DAG_SNARLS

#include "zip_code_tree.hpp"

#include "crash.hpp"

#define debug_parse

using namespace std;
namespace vg {

void ZipCodeTree::fill_in_tree(vector<Seed>& all_seeds, const SnarlDistanceIndex& distance_index) {
    if (all_seeds.size() == 0) {
        return;
    }
    seeds = &all_seeds;

    /*
    Constructor for the ZipCodeTree
    Takes a vector of seeds and constructs the tree

    Tree construction is done by first sorting the seeds along chains/snarls
    Then, adding each seed, snarl/chain boundary, and distance to zip_code_tree
    Finally (optionally), the tree is refined to take out unnecessary edges
    */

    //////////////////// Sort the seeds


    //Helper function to get the orientation of a snarl tree node at a given depth
    //does the same thing as the zipcode decoder's get_is_reversed_in_parent, except
    //that is also considers chains that are children of irregular snarls.
    //We assume that all snarls are DAGs, so all children of snarls must only be 
    //traversable in one orientation through the snarl. In a start-to-end traversal
    //of a snarl, each node will only be traversable start-to-end or end-to-start.
    //If it is traversable end-to-start, then it is considered to be oriented
    //backwards in its parent
    auto get_is_reversed_at_depth =  [&] (const Seed& seed, size_t depth) {
        if (seed.zipcode_decoder->get_is_reversed_in_parent(depth)) {
            return true;
        } else if (depth > 0 && seed.zipcode_decoder->get_code_type(depth-1) == IRREGULAR_SNARL) {
            //If the parent is an irregular snarl, then check the orientation of the child in the snarl
            net_handle_t snarl_handle = seed.zipcode_decoder->get_net_handle(depth-1, &distance_index);
            size_t rank = seed.zipcode_decoder->get_rank_in_snarl(depth);
            if (distance_index.distance_in_snarl(snarl_handle, 0, false, rank, false) 
                        == std::numeric_limits<size_t>::max() 
                && 
                distance_index.distance_in_snarl(snarl_handle, 1, false, rank, true) 
                        == std::numeric_limits<size_t>::max()) {
                //If the distance from the start of the snarl to the start of the child is infinite
                //and the distance from the end of the snarl to the end of the child is infinite
                //then we assume that this child is "reversed" in the parent snarl
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }

    };

    //A vector of indexes into seeds
    //To be sorted along each chain/snarl the snarl tree
    vector<size_t> seed_indices (seeds->size(), 0);
    for (size_t i = 0 ; i < seed_indices.size() ; i++) {
        seed_indices[i] = i;
    }
    assert(seeds->size() == seed_indices.size());

    //Sort the indices
    std::sort(seed_indices.begin(), seed_indices.end(), [&] (const size_t& a, const size_t& b) {
        for (auto x : seed_indices) {
            assert (x < seed_indices.size());
        }
        assert(a < seeds->size());
        assert(b < seeds->size());
#ifdef DEBUG_ZIP_CODE_TREE
        cerr << "Comparing seeds " << seeds->at(a).pos << " and " << seeds->at(b).pos << endl;
#endif
        //Comparator returning a < b
        size_t depth = 0;

        //Keep track of the orientation of each seed
        //Everything should be sorted according to the orientation in the top-level structure,
        //so if things are traversed backwards, reverse the orientation
        bool a_is_reversed = false;
        bool b_is_reversed = false;
        while (depth < seeds->at(a).zipcode_decoder->max_depth() &&
               depth < seeds->at(b).zipcode_decoder->max_depth() &&
               ZipCodeDecoder::is_equal(*seeds->at(a).zipcode_decoder, *seeds->at(b).zipcode_decoder, depth)) {

            //Remember the orientation
            if (get_is_reversed_at_depth(seeds->at(a), depth)) { 
                a_is_reversed = !a_is_reversed;
            }
            if (get_is_reversed_at_depth(seeds->at(b), depth)) {
                b_is_reversed = !b_is_reversed;
            }

            depth++;
        }

        //Remember the orientation of the parent too
        size_t parent_of_a_is_reversed = a_is_reversed;

        //Check the orientations one last time
        if (get_is_reversed_at_depth(seeds->at(a), depth)) { 
            a_is_reversed = !a_is_reversed;
        }
        if (get_is_reversed_at_depth(seeds->at(b), depth)) {
            b_is_reversed = !b_is_reversed;
        }
        
#ifdef DEBUG_ZIP_CODE_TREE
        cerr << "\t different at depth " << depth << endl;
#endif
        //Either depth is the last thing in a or b, or they are different at this depth


        if ( ZipCodeDecoder::is_equal(*seeds->at(a).zipcode_decoder, *seeds->at(b).zipcode_decoder, depth)) {
#ifdef DEBUG_ZIP_CODE_TREE
            cerr << "\tthey are on the same node" << endl;
#endif
            //If they are equal, then they must be on the same node

            size_t offset1 = is_rev(seeds->at(a).pos)
                           ? seeds->at(a).zipcode_decoder->get_length(depth) - offset(seeds->at(a).pos) - 1
                           : offset(seeds->at(a).pos);
            size_t offset2 = is_rev(seeds->at(b).pos)
                           ? seeds->at(b).zipcode_decoder->get_length(depth) - offset(seeds->at(b).pos) - 1
                           : offset(seeds->at(b).pos);
            if (!a_is_reversed) {
                //If they are in a snarl or they are facing forward on a chain, then order by
                //the offset in the node
                return offset1 < offset2;
            } else {
                //Otherwise, the node is facing backwards in the chain, so order backwards in node
                return offset2 < offset1;
            }
        }  else if (depth == 0) {
#ifdef DEBUG_ZIP_CODE_TREE
            cerr << "\tThey are on different connected components" << endl;
#endif
            //If they are on different connected components, sort by connected component
            return seeds->at(a).zipcode_decoder->get_distance_index_address(0) < seeds->at(b).zipcode_decoder->get_distance_index_address(0);
            
        }  else if (seeds->at(a).zipcode_decoder->get_code_type(depth-1) == CHAIN || seeds->at(a).zipcode_decoder->get_code_type(depth-1) == ROOT_CHAIN) {
#ifdef DEBUG_ZIP_CODE_TREE
            cerr << "\t they are children of a common chain" << endl;
#endif
            //If a and b are both children of a chain
            size_t offset_a = seeds->at(a).zipcode_decoder->get_offset_in_chain(depth);
            size_t offset_b = seeds->at(b).zipcode_decoder->get_offset_in_chain(depth);

            if ( offset_a == offset_b) {
                //If they have the same prefix sum, then the snarl comes first
                //They will never be on the same child at this depth
                if (parent_of_a_is_reversed) {
                    return seeds->at(b).zipcode_decoder->get_code_type(depth) != NODE && seeds->at(a).zipcode_decoder->get_code_type(depth) == NODE;  
                } else {
                    return seeds->at(a).zipcode_decoder->get_code_type(depth) != NODE && seeds->at(b).zipcode_decoder->get_code_type(depth) == NODE;  
                }
            } else {
                //Check if the parent chain is reversed and if so, then the order should be reversed
                //The parent could be reversed if it is in an irregular snarl and the 
                if (parent_of_a_is_reversed) {
                    return offset_b < offset_a;
                } else {
                    return offset_a < offset_b;
                }
            }
        } else if (seeds->at(a).zipcode_decoder->get_code_type(depth-1) == REGULAR_SNARL) {
#ifdef DEBUG_ZIP_CODE_TREE
            cerr << "\t they are children of a common regular snarl" << endl;
#endif
            //If the parent is a regular snarl, then sort by order along the parent chain
            size_t offset1 = is_rev(seeds->at(a).pos) 
                           ? seeds->at(a).zipcode_decoder->get_length(depth) - offset(seeds->at(a).pos) - 1
                           : offset(seeds->at(a).pos); 
            size_t offset2 = is_rev(seeds->at(b).pos) 
                           ? seeds->at(b).zipcode_decoder->get_length(depth) - offset(seeds->at(b).pos) - 1
                           : offset(seeds->at(b).pos);
            if (!parent_of_a_is_reversed) {
                return offset1 < offset2;
            } else {
                return offset2 < offset1;
            }
        } else {
#ifdef DEBUG_ZIP_CODE_TREE
            cerr << "\t they are children of a common irregular snarl" << endl;
#endif
            // Otherwise, they are children of an irregular snarl
            // Sort by a topological ordering from the start of the snarl
            // The ranks of children in snarls are in a topological order, so 
            // sort on the ranks
            return seeds->at(a).zipcode_decoder->get_rank_in_snarl(depth) <
                   seeds->at(b).zipcode_decoder->get_rank_in_snarl(depth);
        } 
    });

#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "Sorted positions:" << endl;
    for (const size_t& i : seed_indices) {
        cerr << seeds->at(i).pos << endl;
    }
#endif

    //seed_indices is now sorted roughly along snarls and chains


    ///////////////////// Build the tree

    //For children of snarls, we need to remember the siblings and start bound that came before them
    //so we can record their distances
    //This holds the indices (into zip_code_tree) of each seed or start of a chain,
    // and each start and child chain start of a snarl
    //The children are stored at the depth of their parents. For example, for a root chain,
    //the vector at index 0 would have the chain start, seeds that are on the chain, and the start
    //of snarls on the chain. Similarly, for a top-level snarl at depth 1, the second vector would contain
    //the starts of chains at depth 2 
    //For the children of a chain, the value is the prefix sum in the chain (relative to the orientation 
    //of the top-level chain, not necessarily the chain itself)
    //For the children of a snarl, the value is the index of the seed
    struct child_info_t {
        tree_item_type_t type;  //the type of the item
        size_t value;  //A value associated with the item, could be offset in a chain, index of the seed

        //For the children of snarls, the distance to the left and right of the chain, that gets added to
        //edges in the snarl
        std::pair<size_t, size_t> distances; 
    };
    vector<vector<child_info_t>> sibling_indices_at_depth;

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
        //Make sure sibling_indices_at_depth has enough spaces for this zipcode
        while (sibling_indices_at_depth.size() < current_max_depth+1) {
            sibling_indices_at_depth.emplace_back();
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

        for (size_t depth = 0 ; depth <= max_depth ; depth++) {
            first_different_ancestor_depth = depth;
            
            if (get_is_reversed_at_depth(current_seed, depth)) {

                current_is_reversed = !current_is_reversed;
            }
            if (i != 0 && get_is_reversed_at_depth(previous_seed, depth)) {

                previous_is_reversed = !previous_is_reversed;
            }
            if (!ZipCodeDecoder::is_equal(*current_seed.zipcode_decoder, 
                        *previous_seed.zipcode_decoder, depth)) {
                break;
            } else if (depth == max_depth) {
                same_node = true;
            }
        }
        if (previous_max_depth > current_max_depth) {
            //We might need to update previous_is_reversed
            for (size_t depth = max_depth ; depth <= previous_max_depth ; depth++) {
                
                if (get_is_reversed_at_depth(previous_seed, depth)) {
                    previous_is_reversed = !previous_is_reversed;
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
            code_type_t previous_type = previous_seed.zipcode_decoder->get_code_type(depth);
            if (previous_type == CHAIN || previous_type == ROOT_CHAIN || previous_type == ROOT_NODE) {
#ifdef DEBUG_ZIP_CODE_TREE
                cerr << "\t\tclose a chain at depth " << depth << endl;
#endif

                //Add the end of the chain to the zip code tree
                zip_code_tree.push_back({CHAIN_END, std::numeric_limits<size_t>::max(), false});


                //The distance from the last thing in the chain to the end of the chain
                //will be added to the relevant distances in the parent snarl.
                //Remember that distance in sibling_indices_at_depth for the chain in the snarl
                //
                //If this is reversed, then the distance should be the distance to the start of 
                //the chain. Otherwise, the distance to the end
                //The value that got stored in sibling_indices_at_depth was the prefix sum
                //traversing the chain according to its orientation in the tree, so either way
                //the distance is the length of the chain - the prefix sum
                if (previous_type == CHAIN) {
#ifdef DEBUG_ZIP_CODE_TREE
                    assert(sibling_indices_at_depth[depth-1].size() > 0);
                    assert(sibling_indices_at_depth[depth-1].back().type == CHAIN_START);
#endif
                    //Only add the distance for a non-root chain
                    if ( sibling_indices_at_depth[depth].back().type == SEED) {
                        //If the last thing in the chain was a node, add 1 to include the position
                        sibling_indices_at_depth[depth-1].back().distances.second =
                        SnarlDistanceIndex::sum(1,
                            SnarlDistanceIndex::minus(previous_seed.zipcode_decoder->get_length(depth),
                                                      sibling_indices_at_depth[depth].back().value));
                    } else {
                        //If the last thing in the chain was a snarl, the distance is length-offset
                        sibling_indices_at_depth[depth-1].back().distances.second =
                            SnarlDistanceIndex::minus(previous_seed.zipcode_decoder->get_length(depth),
                                                      sibling_indices_at_depth[depth].back().value);
                    }
                }


            } else if (previous_type == REGULAR_SNARL || previous_type == IRREGULAR_SNARL) { 
#ifdef DEBUG_ZIP_CODE_TREE
                cerr << "\t\tclose a snarl at depth " << depth << endl;
#endif
                //If this is the end of the snarl, then we need to save the distances to 
                //all previous children of the snarl

                zip_code_tree.resize(zip_code_tree.size() + sibling_indices_at_depth[depth].size());

                for (size_t sibling_i = 0 ; sibling_i < sibling_indices_at_depth[depth].size() ; sibling_i++) {
                    const auto& sibling = sibling_indices_at_depth[depth][sibling_i];
                    if (sibling.type == SNARL_START) {
                        //First, the distance between ends of the snarl, which is the length
                        zip_code_tree[zip_code_tree.size() - 1 - sibling_i] = {EDGE,
                            previous_seed.zipcode_decoder->get_length(depth)+1, false};
                    } else {
                        //For the rest of the children, find the distance from the child to
                        //the end
                        //If the child is reversed relative to the top-level chain, then get the distance to start
                        //Also include the distance to the end of the child, sibling.distances.second
                        zip_code_tree[zip_code_tree.size() - 1 - sibling_i] = {EDGE,
                            SnarlDistanceIndex::sum(
                                sibling.distances.second,
                                previous_is_reversed 
                                    ? seeds->at(sibling.value).zipcode_decoder->get_distance_to_snarl_start(depth+1)
                                    : seeds->at(sibling.value).zipcode_decoder->get_distance_to_snarl_end(depth+1)),
                            false};

                    }
                }
                //Note the count of children and the end of the snarl
                zip_code_tree.push_back({NODE_COUNT, sibling_indices_at_depth[depth].size()-1, false});
                zip_code_tree.push_back({SNARL_END, std::numeric_limits<size_t>::max(), false});
            }
            //Update previous_is_reversed to the one before this
            if (depth > 0 && get_is_reversed_at_depth(previous_seed, depth-1)) {
                previous_is_reversed = !previous_is_reversed;
            }

            //Clear the list of children of the thing at this level
            sibling_indices_at_depth[depth].clear();
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
            code_type_t current_type = current_seed.zipcode_decoder->get_code_type(depth);

            if (current_type == NODE || current_type == REGULAR_SNARL || current_type == IRREGULAR_SNARL
                || current_type == ROOT_NODE) {
                //For these things, we need to remember the offset in the node/chain

                if (current_type == ROOT_NODE && sibling_indices_at_depth[depth].empty()) {
                    //If this is a root-level node and the first time we've seen it,
                    //then open the node
                    zip_code_tree.push_back({CHAIN_START, std::numeric_limits<size_t>::max(), false});
                    sibling_indices_at_depth[depth].push_back({CHAIN_START, 0});
                }

                ///////////////// Get the offset in the parent chain (or node)
                size_t current_offset;

                //If we're traversing this chain backwards, then the offset is the offset from the end
                bool current_parent_is_reversed = get_is_reversed_at_depth(current_seed, depth) 
                    ? !current_is_reversed : current_is_reversed;

                //First, get the prefix sum in the chain
                if (current_type == ROOT_NODE) {
                    //Which is 0 if this is just a node
                    current_offset = 0;
                } else {
                    //And the distance to the start or end of the chain if it's a node/snarl in a chain
                    current_offset = current_parent_is_reversed 
                            ? SnarlDistanceIndex::minus(current_seed.zipcode_decoder->get_length(depth-1) ,
                                                        SnarlDistanceIndex::sum(
                                                            current_seed.zipcode_decoder->get_offset_in_chain(depth),
                                                            current_seed.zipcode_decoder->get_length(depth))) 
                            : current_seed.zipcode_decoder->get_offset_in_chain(depth);
                }

                if (depth == current_max_depth) {
                    //If this is a node, then add the offset of the position in the node
                    current_offset = SnarlDistanceIndex::sum(current_offset, 
                        current_is_reversed != is_rev(current_seed.pos)
                            ? current_seed.zipcode_decoder->get_length(depth) - offset(current_seed.pos)
                            : offset(current_seed.pos)+1);
                }

                /////////////////////// Get the offset of the previous thing in the parent chain/node
                size_t previous_offset = depth == 0 ? sibling_indices_at_depth[depth][0].value 
                                                    : sibling_indices_at_depth[depth-1][0].value;
                tree_item_type_t previous_type = depth == 0 ? sibling_indices_at_depth[depth][0].type 
                                                    : sibling_indices_at_depth[depth-1][0].type;


#ifdef DEBUG_ZIP_CODE_TREE
                if (depth > 0) {
                    assert(sibling_indices_at_depth[depth-1].size() == 1);
                }
                //TODO: THis won't always be treu
                //assert(current_offset >= previous_offset);
#endif

                ///////////////////// Record the distance from the previous thing in the chain/node
                if (depth > 1 &&
                     sibling_indices_at_depth[depth-1][0].type == CHAIN_START){
                    //If this is the first thing in a non-root chain or node, remember the distance to the 
                    //start of the chain/node.
                    //This distance will be added to distances in the parent snarl
                    sibling_indices_at_depth[depth-2][0].distances.first = current_offset;

                } else if (!(depth == 0 && sibling_indices_at_depth[depth][0].type == CHAIN_START) &&
                    !(depth > 0 && sibling_indices_at_depth[depth-1][0].type == CHAIN_START)) {
                    //for everything except the first thing in a node/chain
                    size_t distance_between;
                    if (previous_offset > current_offset) {
                        //If the parent is a multicomponent chain, then they might be in different components
                        //TODO: This won't catch all cases of different components in the chain
                        distance_between = std::numeric_limits<size_t>::max();
                    } else {
                        //If the both are seeds or this is a snarl and the previous thing was a seed, 
                        //then add 1 to get to the positions
                        bool current_is_seed = current_type == NODE || current_type == ROOT_NODE;
                        bool previous_is_seed = previous_type == SEED;
                        distance_between = (current_is_seed && previous_is_seed) || (!current_is_seed && previous_is_seed) || (!current_is_seed && !previous_is_seed) 
                            ? current_offset - previous_offset + 1
                            : current_offset - previous_offset; 
                    }

                    zip_code_tree.push_back({EDGE, distance_between, false});
                }

                /////////////////////////////Record this thing in the chain
                if (current_type == NODE || current_type == ROOT_NODE) {
#ifdef DEBUG_ZIP_CODE_TREE
                    cerr << "\t\tContinue node/chain with seed " << seeds->at(seed_indices[i]).pos << " at depth " << depth << endl;
#endif
                    //If this was a node, just remember the seed
                    zip_code_tree.push_back({SEED, seed_indices[i], current_is_reversed != is_rev(seeds->at(seed_indices[i]).pos)});
                } else {
#ifdef DEBUG_ZIP_CODE_TREE
                    cerr << "\t\tOpen new snarl at depth " << depth << endl;
#endif
                    //If this was a snarl, record the start of the snarl
                    zip_code_tree.push_back({SNARL_START, std::numeric_limits<size_t>::max(), false});

                    //Remember the start of the snarl
                    sibling_indices_at_depth[depth].push_back({SNARL_START, std::numeric_limits<size_t>::max()});

                    //For finding the distance to the next thing in the chain, the offset
                    //stored should be the offset of the end bound of the snarl, so add the 
                    //length of the snarl
                    current_offset = SnarlDistanceIndex::sum(current_offset,
                        current_seed.zipcode_decoder->get_length(depth));

                }

                //Remember this thing for the next sibling in the chain
                if (depth == 0) {
                    sibling_indices_at_depth[depth].pop_back();
                    sibling_indices_at_depth[depth].push_back({(current_type == NODE || current_type == ROOT_NODE) ? SEED : SNARL_START, current_offset}); 
                } else {
                    sibling_indices_at_depth[depth-1].pop_back();
                    sibling_indices_at_depth[depth-1].push_back({(current_type == NODE || current_type == ROOT_NODE) ? SEED : SNARL_START, current_offset}); 
                }
                cerr << "Add sibling with type " << current_type << endl;
            } else {
                //Otherwise, this is a chain or root chain
                //If it is a chain, then it is the child of a snarl, so we need to find distances
                //to everything preceding it in the snarl
                assert(current_type == CHAIN || current_type == ROOT_CHAIN);
                if (sibling_indices_at_depth[depth].size() == 0) {
                    //If this is the start of a new chain
#ifdef DEBUG_ZIP_CODE_TREE
                    cerr << "\t\tOpen new chain at depth " << depth << endl;
#endif

                    //For each sibling in the snarl, record the distance from the sibling to this
                    if (current_type == CHAIN) {
                        //If this is the start of a non-root chain, then it is the child of a snarl and 
                        //we need to find the distances to the previous things in the snarl

                        //The distances will be added in reverse order that they were found in
                        zip_code_tree.resize(zip_code_tree.size() + sibling_indices_at_depth[depth-1].size());

                        //If the parent snarl is reversed
                        bool current_parent_is_reversed = get_is_reversed_at_depth(current_seed, depth) 
                            ? !current_is_reversed : current_is_reversed;

                        //The distances in the snarl include the distances to the ends of the child chains
                        //This is the distance to the start of this child (at depth depth+1) in the chain
                        size_t distance_to_start_of_current_child;
                        if (depth == current_max_depth) {
                            //If this is really a node, then get the distance to the start of the node
                            distance_to_start_of_current_child =
                                current_is_reversed != is_rev(current_seed.pos)
                                    ? current_seed.zipcode_decoder->get_length(depth) - offset(current_seed.pos)
                                    : offset(current_seed.pos)+1;
                        } else {
                            //Otherwise, this is really a chain
                            distance_to_start_of_current_child = current_is_reversed 
                                ? SnarlDistanceIndex::minus(current_seed.zipcode_decoder->get_length(depth) ,
                                                            SnarlDistanceIndex::sum(
                                                                current_seed.zipcode_decoder->get_offset_in_chain(depth+1),
                                                                current_seed.zipcode_decoder->get_length(depth+1))) 
                                : current_seed.zipcode_decoder->get_offset_in_chain(depth+1);
                            if (depth+1 == current_max_depth) {
                                //If this is a node, then add the offset of the position in the node
                                bool child_is_reversed = get_is_reversed_at_depth(current_seed, depth+1) 
                                    ? !current_is_reversed : current_is_reversed;
                                distance_to_start_of_current_child = SnarlDistanceIndex::sum(distance_to_start_of_current_child, 
                                    child_is_reversed != is_rev(current_seed.pos)
                                        ? current_seed.zipcode_decoder->get_length(depth+1) - offset(current_seed.pos)
                                        : offset(current_seed.pos)+1);
                            }
                        }

                        for ( size_t sibling_i = 0 ; sibling_i < sibling_indices_at_depth[depth-1].size() ; sibling_i++) {
                            const auto& sibling = sibling_indices_at_depth[depth-1][sibling_i];
                            size_t distance_to_end_of_previous_child = sibling.type == SNARL_START ? 0
                                                                     : sibling.distances.second;
                            if (sibling.type == SNARL_START) {
                                //Get the distance to the start (or end if it's reversed) of the snarl
                                zip_code_tree[zip_code_tree.size() - 1 - sibling_i] = 
                                 {EDGE, 
                                  SnarlDistanceIndex::sum(distance_to_start_of_current_child,
                                    current_parent_is_reversed
                                        ? current_seed.zipcode_decoder->get_distance_to_snarl_end(depth)
                                        : current_seed.zipcode_decoder->get_distance_to_snarl_start(depth)),
                                  false};
                            } else {
                                //Otherwise, the previous thing was another child of the snarl
                                //and we need to record the distance between these two
                                //TODO: This can be improved for simple snarls
                                size_t distance;
                                if (current_type == CHAIN && 
                                    current_seed.zipcode_decoder->get_code_type(depth-1) == REGULAR_SNARL) {
                                    //If this is the child of a regular snarl, then the distance between
                                    //any two chains is inf
                                    distance = std::numeric_limits<size_t>::max();
                                } else {
                                    net_handle_t snarl_handle = current_seed.zipcode_decoder->get_net_handle(depth-1, &distance_index);
                                    size_t rank2 = current_seed.zipcode_decoder->get_rank_in_snarl(depth);
                                    size_t rank1 = seeds->at(sibling.value).zipcode_decoder->get_rank_in_snarl(depth);
                                    //TODO: idk about this distance- I think the orientations need to change
                                    distance = SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(
                                        distance_index.distance_in_snarl(snarl_handle, rank1, false, rank2, false),
                                        distance_to_start_of_current_child),
                                        distance_to_end_of_previous_child);
                                }
                                zip_code_tree[zip_code_tree.size() - 1 - sibling_i] = {EDGE, distance, false};
                            }

                        }
                    }

                    //Now record the start of this chain
                    zip_code_tree.push_back({CHAIN_START, std::numeric_limits<size_t>::max(), false});

                    //Remember the start of the chain, with the prefix sum value
                    sibling_indices_at_depth[depth].push_back({CHAIN_START, 0});

                    //And, if it is the child of a snarl, then remember the chain as a child of the snarl
                    if (depth != 0) {
                        sibling_indices_at_depth[depth-1].push_back({CHAIN_START,
                                                                     seed_indices[i]});
                    }
                }

                if (current_type == CHAIN && depth == current_max_depth) {
                    //If this is a trivial chain, then also add the seed and the distance to the 
                    //thing before it
                    size_t current_offset = current_is_reversed
                            ? current_seed.zipcode_decoder->get_length(depth) - offset(current_seed.pos)
                            : offset(current_seed.pos)+1;

                    if (sibling_indices_at_depth[depth].back().type == CHAIN_START) {
                        //If the previous thing in the "chain" was the start, then don't add the distance,
                        //but remember it to add to snarl distances later
                        sibling_indices_at_depth[depth].back().distances.first = current_offset;
                    } else {
                        zip_code_tree.push_back({EDGE, 
                                                 current_offset - sibling_indices_at_depth[depth].back().value+1, 
                                                 false}); 
                    }
                    zip_code_tree.push_back({SEED, seed_indices[i], current_is_reversed != is_rev(seeds->at(seed_indices[i]).pos)}); 

                    //And update sibling_indices_at_depth to remember this child
                    sibling_indices_at_depth[depth].pop_back();
                    sibling_indices_at_depth[depth].push_back({SEED, current_offset});
                    
                }
            }
            
            //Finished with this depth, so update current_is_reversed to be for the next ancestor
            if (depth < current_max_depth && get_is_reversed_at_depth(current_seed, depth+1)) {
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
        if (get_is_reversed_at_depth(last_seed, depth)) {
            last_is_reversed = !last_is_reversed;
        }
    }
    for (int depth = last_max_depth ; depth >= 0 ; depth--) {
        if (sibling_indices_at_depth[depth].size() > 0) {
            code_type_t last_type = last_seed.zipcode_decoder->get_code_type(depth);
            if (last_type == CHAIN || last_type == ROOT_CHAIN || last_type == ROOT_NODE) {
#ifdef DEBUG_ZIP_CODE_TREE
                cerr << "\t\tclose a chain at depth " << depth << endl;
#endif
                //Add the end of the chain to the zip code tree
                // TODO: When we get C++20, change this to emplace_back aggregate initialization
                zip_code_tree.push_back({CHAIN_END, std::numeric_limits<size_t>::max(), false});


                //The distance from the last thing in the chain to the end of the chain
                //will be added to the relevant distances in the parent snarl.
                //Remember that distance in sibling_indices_at_depth for the chain in the snarl
                //
                //If this is reversed, then the distance should be the distance to the start of 
                //the chain. Otherwise, the distance to the end
                //The value that got stored in sibling_indices_at_depth was the prefix sum
                //traversing the chain according to its orientation in the tree, so either way
                //the distance is the length of the chain - the prefix sum
                if (last_type == CHAIN) {
#ifdef DEBUG_ZIP_CODE_TREE
                    assert(sibling_indices_at_depth[depth-1].size() > 0);
                    assert(sibling_indices_at_depth[depth-1].back().type == CHAIN_START);
#endif
                    if (sibling_indices_at_depth[depth].back().type == SEED) {
                        //If the previous child was a seed, add 1 to the distance to include the position
                        sibling_indices_at_depth[depth-1].back().distances.second =
                            SnarlDistanceIndex::sum(1,
                                SnarlDistanceIndex::minus(last_seed.zipcode_decoder->get_length(depth),
                                                      sibling_indices_at_depth[depth].back().value));
                    } else {
                        //If the previous child was a snarl, don't add 1
                        sibling_indices_at_depth[depth-1].back().distances.second =
                            SnarlDistanceIndex::minus(last_seed.zipcode_decoder->get_length(depth),
                                                      sibling_indices_at_depth[depth].back().value);
                    }
                }

            } else if (last_type == REGULAR_SNARL || last_type == IRREGULAR_SNARL) { 
#ifdef DEBUG_ZIP_CODE_TREE
               cerr << "\t\tclose a snarl at depth " << depth << endl;
#endif
                //If this is the end of the snarl, then we need to save the distances to 
                //all previous children of the snarl

                zip_code_tree.resize(zip_code_tree.size() + sibling_indices_at_depth[depth].size());

                for (size_t sibling_i = 0 ; sibling_i < sibling_indices_at_depth[depth].size() ; sibling_i++) {
                    const auto& sibling = sibling_indices_at_depth[depth][sibling_i];
                    if (sibling.type == SNARL_START) {
                        //First, the distance between ends of the snarl, which is the length
                        zip_code_tree[zip_code_tree.size() - 1 - sibling_i] = {EDGE,
                            last_seed.zipcode_decoder->get_length(depth), false};
                    } else {
                        //For the rest of the children, find the distance from the child to
                        //the end
                        //If the child is reversed relative to the top-level chain, then get the distance to start
                        //Remember to add the distance to the end of the child
                        zip_code_tree[zip_code_tree.size() - 1 - sibling_i] = {EDGE,
                            SnarlDistanceIndex::sum(
                                last_is_reversed 
                                    ? seeds->at(sibling.value).zipcode_decoder->get_distance_to_snarl_start(depth)
                                    : seeds->at(sibling.value).zipcode_decoder->get_distance_to_snarl_end(depth),
                                 sibling.distances.second),
                             false};
                    }
                }
                //Note the count of children and the end of the snarl
                zip_code_tree.push_back({NODE_COUNT, sibling_indices_at_depth[depth].size()-1, false});
                zip_code_tree.push_back({SNARL_END, std::numeric_limits<size_t>::max(), false});
            }
        }
        //Update last_is_reversed to the one before this
        if (depth > 0 && get_is_reversed_at_depth(last_seed, depth-1)) {
            last_is_reversed = !last_is_reversed;
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
        if (current_item.type == SNARL_START) {
            //For the start of a snarl, make a note of the depth to check the next seed
            snarl_depths.emplace_back(current_depth);

            //Increment the depth
            current_depth++;
        } else if (current_item.type == CHAIN_START) {
            //For the start of a chain, increment the depth
            current_depth++;
        } else if (current_item.type == CHAIN_END || current_item.type == SNARL_END) {
            //For the end of a snarl or chain, decrement the depth
            current_depth--;
        } else if (current_item.type == SEED) {
            //If this is a seed, check the snarls we've seen previously
            for (const size_t& snarl_depth : snarl_depths) {
                if (seeds[current_item.value].zipcode_decoder->get_code_type(snarl_depth) == REGULAR_SNARL) {
                    //If this is a regular snarl, then it must be a DAG too
                    dag_count++;
                } else {
                    //If this is an irregular snarl

                    //Check the snarl in the distance index
                    net_handle_t snarl_handle = seeds[current_item.value].zipcode_decoder->get_net_handle(snarl_depth, &distance_index);
#ifdef DEBUG_ZIP_CODE_TREE
                    assert(seeds[current_item.value].zipcode_decoder->get_code_type(snarl_depth) == IRREGULAR_SNARL ||
                           seeds[current_item.value].zipcode_decoder->get_code_type(snarl_depth) == ROOT_SNARL);
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

void ZipCodeTree::validate_zip_tree(const SnarlDistanceIndex& distance_index) const {
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
        bool start_is_reversed = start_itr_left->is_reversed != is_rev(start_seed.pos);

        //Walk through the tree starting from the vector iterator going left, and check the distance
        for (reverse_iterator tree_itr_left (start_itr_left, zip_code_tree.rend()) ;
             tree_itr_left != reverse_iterator(zip_code_tree.rend(), zip_code_tree.rend()) ;
             ++tree_itr_left) {
            seed_result_t next_seed_result = *tree_itr_left;
            const Seed& next_seed = seeds->at(next_seed_result.seed);
            const bool next_is_reversed = next_seed_result.is_reverse != is_rev(next_seed.pos);

            size_t tree_distance = next_seed_result.distance;

            net_handle_t start_handle = distance_index.get_node_net_handle(
                                            id(start_seed.pos),
                                            is_rev(start_seed.pos) != start_is_reversed);
            net_handle_t next_handle = distance_index.get_node_net_handle(
                                            id(next_seed.pos),
                                            is_rev(next_seed.pos) != next_is_reversed);
            cerr << "Distance between " << next_seed.pos << (next_is_reversed ? "rev" : "") << " and " << start_seed.pos << (start_is_reversed ? "rev" : "") << endl;
            cerr << "Values: " << id(next_seed.pos) << " " <<  (is_rev(next_seed.pos) != next_is_reversed ? "rev" : "fd" ) << " " << 
                    (is_rev(next_seed.pos) == next_is_reversed ? offset(next_seed.pos)
                                                                  : distance_index.minimum_length(next_handle) - offset(next_seed.pos) - 1) << " " << 
                    id(start_seed.pos) << " " <<  (is_rev(start_seed.pos) != start_is_reversed ? "rev" : "fd")<< " " << 
                    (is_rev(start_seed.pos) == start_is_reversed ? offset(start_seed.pos)
                                                                  : distance_index.minimum_length(start_handle) - offset(start_seed.pos) - 1 ) << endl;

            size_t index_distance = distance_index.minimum_distance(id(next_seed.pos), next_is_reversed,
                    is_rev(next_seed.pos) == next_is_reversed ? offset(next_seed.pos)
                                                                  : distance_index.minimum_length(next_handle) - offset(next_seed.pos) - 1 ,
                    id(start_seed.pos), start_is_reversed,
                    is_rev(start_seed.pos) == start_is_reversed ? offset(start_seed.pos)
                                                                  : distance_index.minimum_length(start_handle) - offset(start_seed.pos) - 1 
                    );
            cerr << "Tree distance: " << tree_distance << " index distance: " << index_distance << endl;
            assert(tree_distance == index_distance);
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
            // Add value into running distance.
            // Except the stored distance seems to be 1 more than the actual distance.
            // TODO: why?
            
            if(it->value == 0 || it->value == std::numeric_limits<size_t>::max()) {
                // TODO: We assume a 0 distance can't be crossed because it is really infinite.
                // TODO: Which of these are actually supposed to mean that?
                
                // Adjust top of stack to distance limit so we hit the stopping condition.
                top() = distance_limit;
            } else {
                // Add in the actual distance
                top() += (it->value - 1);
            } 
            if (top() > distance_limit) {
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
            if (it->value == std::numeric_limits<size_t>::max()) {
                // Unreachable placeholder, so push it
                push(std::numeric_limits<size_t>::max());
                // And make it be under parent running distance.
                swap();
            } else {
                // We need to add this actual number to parent running distance.
                // Duplicate parent running distance
                dup();
                // Add in the edge value to make a running distance for the thing this edge is for.
                // TODO: We subtract out 1 for snarl edge distances; should we be doing that here???
                top() += it->value;
                // Flip top 2 elements, so now parent running distance is on top, over edge running distance.
                swap();
            }
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
                if (top() > distance_limit) {
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
            if (top() > distance_limit) {
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
