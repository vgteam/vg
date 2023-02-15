#include "zip_code.hpp"

//#define DEBUG_ZIP_CODE

namespace vg{
using namespace std;

void zip_code_t::fill_in_zip_code (const SnarlDistanceIndex& distance_index, const pos_t& pos) {

    std::vector<net_handle_t> ancestors;
    net_handle_t current_handle = distance_index.get_node_net_handle(id(pos));

    //Put all ancestors of the node in a vector, starting from the node, and not including the root
    while (!distance_index.is_root(current_handle)) {
        ancestors.emplace_back(current_handle);
        current_handle = distance_index.get_parent(current_handle);
    }


    //Now add the root-level snarl or chain
    if (distance_index.is_root_snarl(current_handle)) {
        //FIrst thing is a snarl, so add the snarl's connected component number
        zip_code.add_value(0);
#ifdef DEBUG_ZIP_CODE
        cerr << "Adding code for top-level snarl " << distance_index.net_handle_as_string(current_handle) << endl;
#endif
        zip_code.add_value(distance_index.get_connected_component_number(current_handle));
    } else {
        //FIrst thing is a chain so add its connected component number and remove the chain from the stack
        zip_code.add_value(1);

        //If the root-level structure is actually a chain, then save the connected component number and take out
        //the chain from the stack
        //If the root-level structure is a trivial chain, then just store the node (as a chain, which will have the 
        //connected-component number as the rank in the snarl anyways)
        if (!distance_index.is_trivial_chain(ancestors.back())) {
#ifdef DEBUG_ZIP_CODE
        cerr << "Adding code for top-level chain" << endl;
#endif
            zip_code.add_value(distance_index.get_connected_component_number(ancestors.back()));
            ancestors.pop_back();
        }
    }

    //Go through the ancestors top (root) down and add them to the zip code
    for (int i = ancestors.size()-1 ; i >= 0 ; i--) {
        net_handle_t current_ancestor = ancestors[i];
#ifdef DEBUG_ZIP_CODE
        cerr << "Adding code for " << distance_index.net_handle_as_string(current_ancestor) << endl;
#endif
        if (distance_index.is_node(current_ancestor)) {
            vector<size_t> to_add = get_node_code(current_ancestor, distance_index);
            for (auto& x : to_add) {
                zip_code.add_value(x);
            }
        } else if (distance_index.is_chain(current_ancestor)) {
            vector<size_t> to_add = get_chain_code(current_ancestor, distance_index);
            for (auto& x : to_add) {
                zip_code.add_value(x);
            }
            if (distance_index.is_trivial_chain(current_ancestor)) {
                return;
            }
        } else if (distance_index.is_regular_snarl(current_ancestor)) {
            vector<size_t> to_add = get_regular_snarl_code(current_ancestor, ancestors[i-1], distance_index); 
            for (auto& x : to_add) {
                zip_code.add_value(x);
            }
        } else {
#ifdef DEBUG_ZIP_CODE
            assert(distance_index.is_snarl(current_ancestor));
#endif
            vector<size_t> to_add =get_irregular_snarl_code(current_ancestor, distance_index);
            for (auto& x : to_add) {
                zip_code.add_value(x);
            }
        }
    }
}

zip_code_decoder_t zip_code_t::decode() const {
    zip_code_decoder_t result;
    
    size_t zip_index, zip_value;
    std::tie(zip_value, zip_index) = zip_code.get_value_and_next_index(0);

    //Is the root a chain/node?
    bool is_chain = zip_value;
    result.emplace_back(is_chain, 0);



    //The next thing is the connected-component number
    std::tie(zip_value, zip_index) = zip_code.get_value_and_next_index(zip_index);

    //If the top-level structure is a chain, it might actually be a node, in which case
    //the only other thing that got stored is the length
    if (is_chain) {
        if (zip_code.get_value_and_next_index(zip_index).second == std::numeric_limits<size_t>::max()) {
            //If the zip code ends here, then this was a node and we're done
            return result;
        }
    }
    is_chain=!is_chain;

    //And then the codes start
    while (zip_index != std::numeric_limits<size_t>::max()) {
        //Remember this
        result.emplace_back(is_chain, zip_index);

        //And get to the next thing
        if (is_chain) {
            //If the current zip_index points to a chain (or a node)
            std::tie(zip_value, zip_index) = zip_code.get_value_and_next_index(zip_index);
            std::tie(zip_value, zip_index) = zip_code.get_value_and_next_index(zip_index);

            //This might be a node that is a child of the chain, in which case there is one
            //more thing in the zip code

            if (zip_index != std::numeric_limits<size_t>::max() &&
                zip_code.get_value_and_next_index(zip_index).second == std::numeric_limits<size_t>::max()) {
                //If the zip code ends here, then this was a node and we're done
                return result;
            }

        } else {
            //If the last zip_index pointed to a chain, then this should point to a snarl, unless it is
            //the last thing in the code, in which case it is a node in a chain
            //So if there are only 3 things left in the zip code, then this is a node

            //The regular/irregular snarl tag
            std::tie(zip_value, zip_index) = zip_code.get_value_and_next_index(zip_index);

            zip_index = zip_code.get_value_and_next_index(zip_index).second;

            if (zip_value) {
                //Regular snarl, so 2 remaining things in the code
                std::tie(zip_value, zip_index) = zip_code.get_value_and_next_index(zip_index);
                if (zip_index == std::numeric_limits<size_t>::max()) {
                    //If the zip code ends here, then this was a node, not a snarl, so
                    //take out the last snarl and replace it with a node
                    size_t last_index = result.back().second;
                    result.pop_back();
                    result.emplace_back(true, last_index);
                    return result;
                }
                std::tie(zip_value, zip_index) = zip_code.get_value_and_next_index(zip_index);
            } else {
                //If it was an irregular snarl, then we're already at the end but check to see if this was 
                //actually a node at the end of the zip code
                if (zip_code.get_value_and_next_index(zip_index).second == std::numeric_limits<size_t>::max()) {
                    //If the zip code ends here, then this was a node, not a snarl, so
                    //take out the last snarl and replace it with a node
                    size_t last_index = result.back().second;
                    result.pop_back();
                    result.emplace_back(true, last_index);
                    return result;
                }
            }
        }
        is_chain = !is_chain;
    }
    return result;
}

decoded_code_t zip_code_t::decode_one_code(size_t index, const code_type_t& code_type, const SnarlDistanceIndex& distance_index) const {
    if (code_type == ROOT_CHAIN || code_type == ROOT_SNARL || 
        ((code_type == CHAIN || code_type == REGULAR_SNARL || code_type == IRREGULAR_SNARL) && index == 0)) {
        //Only need the rank
        size_t rank = zip_code.get_value_and_next_index(zip_code.get_value_and_next_index(index).second).first;
        if (code_type == ROOT_CHAIN || code_type == CHAIN ) {
            return decoded_code_t {distance_index.get_root(),
                                   std::numeric_limits<size_t>::max(),
                                   rank,
                                   ROOT_CHAIN, 
                                   false};
        } else {
            return decoded_code_t {distance_index.get_handle_from_connected_component(rank),
                                   std::numeric_limits<size_t>::max(),
                                   rank,
                                   ROOT_SNARL, 
                                   false};
        }
    } else if (code_type == ROOT_NODE || (code_type == NODE && index == 0)) {
        size_t rank;
        //Get the second thing (rank) and the index of the next thing (length)
        std::tie(rank, index) = zip_code.get_value_and_next_index(zip_code.get_value_and_next_index(index).second);
        size_t length = zip_code.get_value_and_next_index(index).first;
        return decoded_code_t { distance_index.get_root(),
                                (length == 0 ? std::numeric_limits<size_t>::max() : length-1),
                                rank,
                                ROOT_NODE, false};
    } else if (code_type == NODE) {
        size_t prefix_sum;
        std::tie(prefix_sum, index) = zip_code.get_value_and_next_index(index);
        size_t length; 
        std::tie(length, index) = zip_code.get_value_and_next_index(index);
        bool is_rev = zip_code.get_value_and_next_index(index).first;
        return decoded_code_t {distance_index.get_root(),
                               (length == 0 ? std::numeric_limits<size_t>::max() : length-1), 
                               (prefix_sum == 0 ? std::numeric_limits<size_t>::max() : prefix_sum-1), 
                               code_type, is_rev}; 
    } else if (code_type == CHAIN) {
        size_t rank;
        std::tie(rank, index) = zip_code.get_value_and_next_index(index);
        size_t length = zip_code.get_value_and_next_index(index).first;
        return decoded_code_t {distance_index.get_root(),
                               (length == 0 ? std::numeric_limits<size_t>::max() : length-1), 
                               rank, 
                               code_type, false}; 
    } else if (code_type == REGULAR_SNARL || code_type == IRREGULAR_SNARL) {
        net_handle_t handle = distance_index.get_root();
        bool is_regular; 
        size_t rank;
        size_t length; 
        bool is_rev;
        std::tie(is_regular, index) = zip_code.get_value_and_next_index(index);
        std::tie(rank, index) = zip_code.get_value_and_next_index(index);
        if (is_regular) {
            //If this is a regular snarl, then the values are found from the zip code
            std::tie(length, index) = zip_code.get_value_and_next_index(index);
            if (length == 0) {
                length = std::numeric_limits<size_t>::max();
            } else {
                length -= 1;
            }
            is_rev = zip_code.get_value_and_next_index(index).first;
            if (rank == 0) {
                rank = std::numeric_limits<size_t>::max();
            } else {
                rank -= 1;
            }
        } else {
            //If it's irregular, then they are found from the distance index
            //The rank stored was actually the location in the distance index
            handle = distance_index.get_net_handle_from_values(
                rank, SnarlDistanceIndex::START_END, SnarlDistanceIndex::SNARL_HANDLE);

            net_handle_t start_node = distance_index.get_node_from_sentinel(distance_index.get_bound(handle, false, false));
            rank = distance_index.get_prefix_sum_value(start_node) + distance_index.minimum_length(start_node);

            length = distance_index.minimum_length(handle);
            is_rev = false;
        }
        return decoded_code_t {handle, 
                               length, 
                               rank, 
                               is_regular ? REGULAR_SNARL : IRREGULAR_SNARL, 
                               is_rev}; 
    } else {
        throw std::runtime_error("zipcode: invalid code type");
    }
}

vector<size_t> zip_code_t::get_node_code(const net_handle_t& node, const SnarlDistanceIndex& distance_index) {
#ifdef DEBUG_ZIP_CODE
    assert(!distance_index.is_trivial_chain(node));
    assert((distance_index.is_chain(distance_index.get_parent(node)) || distance_index.is_root(distance_index.get_parent(node))));
#endif
    //Node code is: offset in chain, length, is reversed
    vector<size_t> node_code;
    //Assume this node is in a regular chain
    size_t prefix_sum = distance_index.get_prefix_sum_value(node); 
    node_code.emplace_back(prefix_sum == std::numeric_limits<size_t>::max() ? 0 : prefix_sum+1);
    node_code.emplace_back(distance_index.minimum_length(node)+1);
    node_code.emplace_back(distance_index.is_reversed_in_parent(node));
    return node_code;

}
vector<size_t> zip_code_t::get_chain_code(const net_handle_t& chain, const SnarlDistanceIndex& distance_index) {
    //Chain code is: rank in snarl, length
    vector<size_t> chain_code;
    chain_code.emplace_back(distance_index.get_rank_in_parent(chain));
    size_t len = distance_index.minimum_length(chain);
    chain_code.emplace_back(len == std::numeric_limits<size_t>::max() ? 0 : len+1);
    return chain_code;

}
vector<size_t> zip_code_t::get_regular_snarl_code(const net_handle_t& snarl, const net_handle_t& snarl_child, const SnarlDistanceIndex& distance_index) {
    //Regular snarl code is 1, offset in chain, length, is reversed
    vector<size_t> snarl_code;

    //Tag to say that it's a regular snarl
    snarl_code.emplace_back(1);

    //Chain prefix sum value for the start of the snarl, which is the prefix sum of the start node + length of the start node
    net_handle_t start_node = distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, false, false));
    size_t prefix_sum = SnarlDistanceIndex::sum(distance_index.get_prefix_sum_value(start_node), distance_index.minimum_length(start_node));
    snarl_code.emplace_back(prefix_sum == std::numeric_limits<size_t>::max() ? 0 : prefix_sum+1);

    //Length of the snarl
    size_t len = distance_index.minimum_length(snarl);
    snarl_code.emplace_back(len == std::numeric_limits<size_t>::max() ? 0 : len+1);

    //Is the child of the snarl reversed in the snarl
#ifdef DEBUG_ZIP_CODE
    assert(distance_index.is_chain(snarl_child));
#endif
    snarl_code.emplace_back(distance_index.distance_in_parent(snarl, distance_index.get_bound(snarl, false, true),
                                                              distance_index.flip(distance_index.canonical(snarl_child))) != 0);

    return snarl_code;

}
vector<size_t> zip_code_t::get_irregular_snarl_code(const net_handle_t& snarl, const SnarlDistanceIndex& distance_index) {
    //Regular snarl code is 0, snarl record offset
    vector<size_t> snarl_code;

    //Tag to say that it's an irregular snarl
    snarl_code.emplace_back(0);

    //Record offset to look up distances in the index later
    snarl_code.emplace_back(distance_index.get_record_offset(snarl));

    return snarl_code;

}

size_t zip_code_t::minimum_distance_between(const zip_code_t& zip1, const pos_t& pos1,   
    const zip_code_t& zip2, const pos_t& pos2, const SnarlDistanceIndex& distance_index){
#ifdef DEBUG_ZIP_CODE
    zip_code_t check_zip1;
    check_zip1.fill_in_zip_code(distance_index, pos1);
    assert(zip1 == check_zip1);

    zip_code_t check_zip2;
    check_zip2.fill_in_zip_code(distance_index, pos2);
    assert(zip2 == check_zip2);
#endif

    //Helper function to update the distances to the ends of the parent
    //distance_start and distance_end get updated
    auto update_distances_to_ends_of_parent = [&] (const decoded_code_t& child_code, const decoded_code_t& parent_code, 
                                            size_t& distance_to_start, size_t& distance_to_end) {
        //The distances from the start/end of current child to the start/end(left/right) of the parent
        size_t distance_start_left, distance_start_right, distance_end_left, distance_end_right;
        if (parent_code.code_type == IRREGULAR_SNARL) {
            distance_start_left = distance_index.distance_in_snarl(parent_code.net_handle, 
                                    child_code.rank_or_offset, false, 0, false);
            distance_start_right = distance_index.distance_in_snarl(parent_code.net_handle, 
                                    child_code.rank_or_offset, false, 1, false);
            distance_end_right = distance_index.distance_in_snarl(parent_code.net_handle, 
                                    child_code.rank_or_offset, true, 1, false);
            distance_end_left = distance_index.distance_in_snarl(parent_code.net_handle, 
                                    child_code.rank_or_offset, true, 0, false);
#ifdef DEBUG_ZIP_CODE
            cerr << "Distances to parent irregular snarl: " << distance_start_left << " " << distance_start_right << " " << distance_end_left << " " << distance_end_right << endl;
#endif
        } else if (parent_code.code_type == REGULAR_SNARL) {
            //If its a regular snarl, then the distances to the ends are either 0 or inf
            //For a regular snarl, the snarl stores if the child was reversed, rather than the child
            if (parent_code.is_reversed) {
                distance_start_left = std::numeric_limits<size_t>::max();
                distance_start_right = 0;
                distance_end_right = std::numeric_limits<size_t>::max();
                distance_end_left = 0;
            } else {
                distance_start_left = 0;
                distance_start_right = std::numeric_limits<size_t>::max();
                distance_end_right = 0;
                distance_end_left = std::numeric_limits<size_t>::max();
            }
#ifdef DEBUG_ZIP_CODE
            cerr << "Distances to parent regular snarl: " << distance_start_left << " " << distance_start_right << " " << distance_end_left << " " << distance_end_right << endl;
#endif
        } else if (parent_code.code_type == CHAIN) {
            if (child_code.code_type == NODE && child_code.is_reversed){ 
                distance_start_left = std::numeric_limits<size_t>::max();
                distance_end_right = std::numeric_limits<size_t>::max();
                //Prefix sum of the child
                distance_end_left = child_code.rank_or_offset;
                //Length of the chain - prefix sum of the child - length of the child
                distance_start_right = SnarlDistanceIndex::minus(SnarlDistanceIndex::minus(
                        parent_code.length, child_code.rank_or_offset), child_code.length);
            } else {
                distance_end_left = std::numeric_limits<size_t>::max();
                distance_start_right = std::numeric_limits<size_t>::max();
                //Prefix sum of the child
                distance_start_left = child_code.rank_or_offset;
                //Length of the chain - prefix sum of the child - length of the child
                distance_end_right = SnarlDistanceIndex::minus(SnarlDistanceIndex::minus(
                        parent_code.length, child_code.rank_or_offset), child_code.length);
            }
#ifdef DEBUG_ZIP_CODE
            cerr << "Distances to parent chain: " << distance_start_left << " " << distance_start_right << " " << distance_end_left << " " << distance_end_right << endl;
#endif
        }

        size_t new_distance_to_start = std::min(SnarlDistanceIndex::sum(distance_start_left, distance_to_start),
                                      SnarlDistanceIndex::sum(distance_end_left, distance_to_end));
        size_t new_distance_to_end = std::min(SnarlDistanceIndex::sum(distance_start_right, distance_to_start),
                                      SnarlDistanceIndex::sum(distance_end_right, distance_to_end));
        distance_to_start = new_distance_to_start;
        distance_to_end = new_distance_to_end;

    };
    size_t zip_index1 = 0; size_t zip_index2 = 0;
    size_t zip_value1 = std::numeric_limits<size_t>::max();
    size_t zip_value2 = std::numeric_limits<size_t>::max();

    //If the two positions aren't on the same connected component, then we're done
    std::tie(zip_value1, zip_index1) = zip1.zip_code.get_value_and_next_index(0);
    std::tie(zip_value2, zip_index2) = zip2.zip_code.get_value_and_next_index(0);
    if (zip_value1 != zip_value2) {
#ifdef DEBUG_ZIP_CODE
        cerr << "Zip codes are on different connected components" << endl;
#endif
        return std::numeric_limits<size_t>::max();
    }
    std::tie(zip_value1, zip_index1) = zip1.zip_code.get_value_and_next_index(zip_index1);
    std::tie(zip_value2, zip_index2) = zip2.zip_code.get_value_and_next_index(zip_index2);
    if (zip_value1 != zip_value2) {
#ifdef DEBUG_ZIP_CODE
        cerr << "Zip codes are on different connected components" << endl;
#endif
        return std::numeric_limits<size_t>::max();
    }

    //The two positions are in the same connected component so now try to find the distance
    zip_code_decoder_t zip1_decoder = zip1.decode();
    zip_code_decoder_t zip2_decoder = zip2.decode();

    //Now find the lowest common ancestor of the two zipcodes
    size_t lowest_common_ancestor_index;
    for (size_t i = 0 ; i < zip1_decoder.size() ; i++) {
        if (i >= zip2_decoder.size()) {
            //Don't go beyond the end of the second zip code
            break;
        } else if (i == zip1_decoder.size()-1 && i == zip2_decoder.size()-1) {
            //If this is the node for both zip codes, then they are the same if the node ids are the same
            if (id(pos1) == id(pos2)) {
                lowest_common_ancestor_index = i;
            } else {
                break;
            }
        } else if (zip1_decoder[i] == zip2_decoder[i]){
            decoded_code_t decoded1 = zip1.decode_one_code(zip1_decoder[i].second, zip1_decoder[i].first ? (zip1_decoder.size() == 1 || (i > 0 && zip1_decoder[i-1].first) ? NODE : CHAIN) 
                                                                                    : REGULAR_SNARL, distance_index); 
            decoded_code_t decoded2 = zip2.decode_one_code(zip2_decoder[i].second, zip2_decoder[i].first ? (zip2_decoder.size() == 1 || (i > 0 && zip2_decoder[i-1].first) ? NODE : CHAIN) 
                                                                                    : REGULAR_SNARL, distance_index); 
            if ( decoded1 == decoded2) {
                lowest_common_ancestor_index = i;
            } else {
                break;
            }
        } else {
            //If they are different, stop looking
            break;
        }
    }
#ifdef DEBUG_ZIP_CODE
    vector<net_handle_t> ancestors;
    net_handle_t ancestor = distance_index.get_node_net_handle(id(pos1));
    while (!distance_index.is_root(ancestor)) {
        ancestors.push_back(ancestor);
        ancestor = distance_index.get_parent(ancestor);
    }
    ancestors.push_back(ancestor);
    cerr << "The lowest common ancestor is the " << lowest_common_ancestor_index << "th thing from the root" << endl;
    cerr << "That should be " << distance_index.net_handle_as_string(ancestors[ancestors.size() - lowest_common_ancestor_index - 1]) << endl; 
#endif

    //Get the decoded node (or technically chain if it's a trivial chain in a snarl)
    decoded_code_t current_code1 = zip1.decode_one_code(zip1_decoder.back().second, 
                zip1_decoder.size() == 1 ? ROOT_NODE : (
                    zip1_decoder[zip1_decoder.size()-2].first ? NODE : CHAIN), distance_index);
    decoded_code_t current_code2 = zip2.decode_one_code(zip2_decoder.back().second, 
                zip2_decoder.size() == 1 ? ROOT_NODE : (
                    zip2_decoder[zip2_decoder.size()-2].first ? NODE : CHAIN), distance_index); 

    size_t distance_to_start1 = is_rev(pos1) ? current_code1.length - offset(pos1) : offset(pos1) + 1;
    size_t distance_to_end1 = is_rev(pos1) ? offset(pos1) + 1 : current_code1.length - offset(pos1);
    size_t distance_to_start2 = is_rev(pos2) ? current_code2.length - offset(pos2) : offset(pos2) + 1;
    size_t distance_to_end2 = is_rev(pos2) ? offset(pos2) + 1 : current_code2.length - offset(pos2);

    //These are directed distances so set backwards distances to inf
    if (is_rev(pos1)) {
        distance_to_end1 = std::numeric_limits<size_t>::max();
    } else {
        distance_to_start1 = std::numeric_limits<size_t>::max();
    }
    if (is_rev(pos2)) {
        distance_to_start2 = std::numeric_limits<size_t>::max();
    } else {
        distance_to_end2 = std::numeric_limits<size_t>::max();
    }

#ifdef DEBUG_ZIP_CODE
cerr << "Distances in nodes: " << distance_to_start1 << " " << distance_to_end1 << " " << distance_to_start2 << " " << distance_to_end2 << endl;
cerr << "Finding distances to ancestors of first position" << endl;
#endif


    //Now walk up the snarl tree from each position to one level below the lowest common ancestor
    for (int i = zip1_decoder.size()-2 ; i > 0 && i > lowest_common_ancestor_index ; i--) {
        //current_code1 is the child of parent_code1, which is at index i
        //The distances are currently to the ends of current_code1
        //FInd the distances to the ends of parent_code1

        decoded_code_t parent_code1 = zip1.decode_one_code(zip1_decoder[i].second,
            zip1_decoder[i].first ? CHAIN : REGULAR_SNARL, distance_index);
#ifdef DEBUG_ZIP_CODE
        assert(parent_code1.code_type != NODE);
        assert(parent_code1.code_type != ROOT_NODE);
        assert(parent_code1.code_type != ROOT_SNARL);
        assert(parent_code1.code_type != ROOT_CHAIN);
#endif
        update_distances_to_ends_of_parent(current_code1, parent_code1, distance_to_start1, distance_to_end1);
        current_code1 = std::move(parent_code1);
    }
#ifdef DEBUG_ZIP_CODE
cerr << "Finding distances to ancestors of second position" << endl;
#endif
    //The same thing for the second position
    for (int i = zip2_decoder.size()-2 ; i > 0 && i > lowest_common_ancestor_index ; i--) {
        //current_code2 is the child of parent_code2, which is at index i
        //The distances are currently to the ends of current_code2
        //FInd the distances to the ends of parent_code2

        decoded_code_t parent_code2 = zip2.decode_one_code(zip2_decoder[i].second,
            zip2_decoder[i].first ? CHAIN : REGULAR_SNARL, distance_index);
#ifdef DEBUG_ZIP_CODE
        assert(parent_code2.code_type != NODE);
        assert(parent_code2.code_type != ROOT_NODE);
        assert(parent_code2.code_type != ROOT_SNARL);
        assert(parent_code2.code_type != ROOT_CHAIN);
#endif
        update_distances_to_ends_of_parent(current_code2, parent_code2, distance_to_start2, distance_to_end2);
        current_code2 = std::move(parent_code2);
    }


    //Distances are now the distances to the ends of a child of the common ancestor

#ifdef DEBUG_ZIP_CODE
    cerr << "Distances in children of common ancestor: " << distance_to_start1 << " " << distance_to_end1 << " " << distance_to_start2 << " " << distance_to_end2 << endl;
    //Check that the current nodes are actually children of the lca
    if (lowest_common_ancestor_index != zip1_decoder.size() - 1) { 
        pair<bool, size_t> zip1_index = zip1_decoder[lowest_common_ancestor_index+1];
        assert(current_code1 == zip1.decode_one_code(zip1_index.second,
                 zip1_index.first ? (zip1_decoder[lowest_common_ancestor_index].first ? NODE : CHAIN) : REGULAR_SNARL, distance_index));
                                
    }
    if (lowest_common_ancestor_index != zip2_decoder.size() - 1) { 
        pair<bool, size_t> zip2_index = zip2_decoder[lowest_common_ancestor_index+1];
        assert(current_code2 == zip2.decode_one_code(zip2_index.second,
                 zip2_index.first ? (zip2_decoder[lowest_common_ancestor_index].first ? NODE : CHAIN) : REGULAR_SNARL, distance_index));
                                
    }
#endif

    //Find the distance between them in the lowest common ancestor

    size_t distance_between = std::numeric_limits<size_t>::max();

    //Walk up the snarl tree from the lca and find the distance between the common ancestor
    for (int i = lowest_common_ancestor_index ; i >= 0 ; i--) {
#ifdef DEBUG_ZIP_CODE
        cerr << "At " << i << "st/th ancestor" << endl;
        cerr << "\tdistances are " << distance_to_start1 << " " << distance_to_end1 << " " << distance_to_start2 << " " << distance_to_end2 << endl;
#endif
        decoded_code_t parent_code;
        if (i == zip1_decoder.size()-1) {
            //If the lca is a node that both positions are on

#ifdef DEBUG_ZIP_CODE
            //If the lca is a node, then both the current_codex's should be the same node
            assert(current_code1 == current_code2);
            assert(i == zip2_decoder.size()-1);
            cerr << "\tAncestor should be a node" << endl;
#endif
            size_t d1 = SnarlDistanceIndex::sum(distance_to_end1, distance_to_start2);
            size_t d2 = SnarlDistanceIndex::sum(distance_to_end2, distance_to_start1);
            if (d1 > current_code1.length) {
                distance_between = std::min(distance_between,
                                            SnarlDistanceIndex::minus(SnarlDistanceIndex::minus(d1, current_code1.length),1));
            } 
            if (d2 > current_code1.length) {
                distance_between = std::min(distance_between,
                                            SnarlDistanceIndex::minus(SnarlDistanceIndex::minus(d2, current_code1.length),1));
            }
            parent_code = std::move(current_code1); 
        } else if ( zip1_decoder[i].first) {
#ifdef DEBUG_ZIP_CODE
            cerr << "\tancestor should be a chain" << endl;
#endif
            //If this ancestor is a chain
            parent_code = zip1.decode_one_code(zip1_decoder[i].second, CHAIN, distance_index);

            //If the children are reversed in the chain, then flip their distances
            if (current_code1.code_type == NODE && current_code1.is_reversed) {
#ifdef DEBUG_ZIP_CODE
                cerr << "Reverse child1 distances" << endl;
#endif
                size_t temp = distance_to_start1;
                distance_to_start1 = distance_to_end1;
                distance_to_end1 = temp;
            }
            if (current_code2.code_type == NODE && current_code2.is_reversed) {
#ifdef DEBUG_ZIP_CODE
                cerr << "Reverse child2 distances" << endl;
#endif
                size_t temp = distance_to_start2;
                distance_to_start2 = distance_to_end2;
                distance_to_end2 = temp;
            }

            //If they are the same child, then there is no path between them in the chain because we don't allow loops
            if (!(current_code1 == current_code2 || (current_code1.code_type == NODE && id(pos1) == id(pos2)))) {
                if (current_code1.rank_or_offset < current_code2.rank_or_offset ||
                    (current_code1.rank_or_offset == current_code2.rank_or_offset &&
                     (current_code1.code_type == REGULAR_SNARL || current_code1.code_type == IRREGULAR_SNARL)
                     && current_code2.code_type == NODE)) {
                    //First child comes first in the chain
                    
                    if (current_code1.code_type == REGULAR_SNARL || current_code1.code_type == IRREGULAR_SNARL) {
                        //If the first thing is a snarl, then we need to take into account the length of the snarl
                        //(prefix sum 2 + distance left 2) - (prefix sum 1 + length 1) + distance right 1

#ifdef DEBUG_ZIP_CODE
                        cerr << "First child comes first in the chain and it is a snarl" << endl;
                        cerr << "Find distances from : " << current_code2.rank_or_offset << " " << distance_to_start2 << " " << current_code1.rank_or_offset << " " << current_code1.length  << " " << distance_to_end1 << endl;
#endif
                        if (distance_to_start2 != std::numeric_limits<size_t>::max()
                            && distance_to_end1 != std::numeric_limits<size_t>::max()) {
                            distance_between = std::min(distance_between,
                                                        SnarlDistanceIndex::minus(SnarlDistanceIndex::sum(
                                                            SnarlDistanceIndex::minus(
                                                                SnarlDistanceIndex::sum(current_code2.rank_or_offset, 
                                                                                        distance_to_start2), 
                                                                SnarlDistanceIndex::sum(current_code1.rank_or_offset,
                                                                                        current_code1.length)),
                                                             distance_to_end1),1));
                        }
                    } else {
                        //Otherwise, all that matters is the prefix sums
                        //(Prefix sum 2  + distance left 2) - (prefix sum1+ length 1) + distance right 1
#ifdef DEBUG_ZIP_CODE
                        cerr << "First child comes first in the chain and it isn't a snarl" << endl;
                        cerr << "Find distances from : " << current_code2.rank_or_offset << " " << distance_to_start2 << " " << current_code1.rank_or_offset << " " << distance_to_end1 << endl;
#endif
                        if (distance_to_start2 != std::numeric_limits<size_t>::max()
                            && distance_to_end1 != std::numeric_limits<size_t>::max()) {
                            distance_between = std::min(distance_between,
                                                        SnarlDistanceIndex::minus(
                                                        SnarlDistanceIndex::sum(
                                                        SnarlDistanceIndex::minus(
                                                            SnarlDistanceIndex::sum(current_code2.rank_or_offset, 
                                                                                    distance_to_start2),
                                                            SnarlDistanceIndex::sum(current_code1.rank_or_offset,
                                                                                    current_code1.length)), 

                                                            distance_to_end1),1) );
                        }
                    }
                } else {
                    //Second child comes first in the chain, or they are the same (doesn't matter)
                    if (current_code2.code_type == REGULAR_SNARL || current_code2.code_type == IRREGULAR_SNARL) {
                        //If the first thing is a snarl, then we need to take into account the length of the snarl
                        //(prefix sum 1 + distance left 1) - (prefix sum 2 + length 2) + distance right 2
#ifdef DEBUG_ZIP_CODE
                        cerr << "Second child comes first in the chain and it is a snarl" << endl;
                        cerr << "Find distances from : " << current_code1.rank_or_offset << " " << distance_to_start1 << " " << current_code2.rank_or_offset << " " << current_code2.length  << " " << distance_to_end2 << endl;
#endif
                        if (distance_to_start1 != std::numeric_limits<size_t>::max() 
                             && distance_to_end2 != std::numeric_limits<size_t>::max() ){
                            distance_between = std::min(distance_between,
                                                        SnarlDistanceIndex::minus(SnarlDistanceIndex::sum(
                                                            SnarlDistanceIndex::minus(
                                                                SnarlDistanceIndex::sum(current_code1.rank_or_offset, 
                                                                                        distance_to_start1), 
                                                                SnarlDistanceIndex::sum(current_code2.rank_or_offset,
                                                                                        current_code2.length)),
                                                             distance_to_end2), 1));
                        }
                    } else {
                        //Otherwise, all that matters is the prefix sums
                        //(Prefix sum 1  + distance left 1) - (prefix sum2 + length 2) + distance right 2
#ifdef DEBUG_ZIP_CODE
                        cerr << "Second child comes first in the chain and it isn't a snarl" << endl;
                        cerr << "Find distances from : " << current_code1.rank_or_offset << " " << distance_to_start1 << " " << current_code2.rank_or_offset << " " << distance_to_end2 << endl;
#endif
                        if (distance_to_start1 != std::numeric_limits<size_t>::max() 
                             && distance_to_end2 != std::numeric_limits<size_t>::max() ){
                            distance_between = std::min(distance_between,
                                                        SnarlDistanceIndex::minus(
                                                        SnarlDistanceIndex::sum(
                                                        SnarlDistanceIndex::minus(
                                                            SnarlDistanceIndex::sum(current_code1.rank_or_offset, 
                                                                                    distance_to_start1),
                                                            SnarlDistanceIndex::sum(current_code2.rank_or_offset,
                                                                                    current_code2.length)), 

                                                            distance_to_end2),1) );
                        }
                    }
                }
            }
        } else {

#ifdef DEBUG_ZIP_CODE
            cerr << "\tancestor is a snarl" << endl;
#endif
            //If the ancestor is a snarl
            parent_code = zip1.decode_one_code(zip1_decoder[i].second, REGULAR_SNARL, distance_index);
            
            //If the parent is a regular snarl, then there is no path between them so
            //just update the distances to the ends of the parent 
            if (parent_code.code_type != REGULAR_SNARL) {
                //Parent may be an irregular snarl or a root snarl (which is also irregular)
#ifdef DEBUG_ZIP_CODE
                cerr << "irregular snarl so find distances in the distance index: " << distance_index.net_handle_as_string(parent_code.net_handle) << endl;
                cerr << "\t at offset " << distance_index.get_record_offset(parent_code.net_handle) << endl;
                cerr << "ranks: " << current_code1.rank_or_offset << " and " << current_code2.rank_or_offset << endl;
#endif

                size_t distance_start_start = distance_index.distance_in_snarl(parent_code.net_handle, 
                                    current_code1.rank_or_offset, false, current_code2.rank_or_offset, false);
                size_t distance_start_end = distance_index.distance_in_snarl(parent_code.net_handle, 
                                    current_code1.rank_or_offset, false, current_code2.rank_or_offset, true);
                size_t distance_end_start = distance_index.distance_in_snarl(parent_code.net_handle, 
                                    current_code1.rank_or_offset, true, current_code2.rank_or_offset, false);
                size_t distance_end_end = distance_index.distance_in_snarl(parent_code.net_handle, 
                                    current_code1.rank_or_offset, true, current_code2.rank_or_offset, true);
                size_t distance_between_snarl = std::min( SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(
                                                distance_to_start1, distance_to_start2), distance_start_start),
                                   std::min( SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(
                                                distance_to_start1, distance_to_end2), distance_start_end),
                                   std::min( SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(
                                                distance_to_end1, distance_to_start2), distance_end_start),
                                             SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(
                                                distance_to_end1, distance_to_end2), distance_end_end))));

                distance_between = std::min(distance_between,
                                            SnarlDistanceIndex::minus(distance_between_snarl, 1));
            }
#ifdef DEBUG_ZIP_CODE
            else {
                cerr << "\tAncestor is a regular snarl so there is no path between the children" << endl;
            }
#endif
            update_distances_to_ends_of_parent(current_code1, parent_code, distance_to_start1, distance_to_end1);
            update_distances_to_ends_of_parent(current_code2, parent_code, distance_to_start2, distance_to_end2);
        }
        current_code1 = parent_code;
        current_code2 = std::move(parent_code);
#ifdef DEBUG_ZIP_CODE
        cerr << "distance in ancestor: " << distance_between << endl;
#endif
    }

    
                                        


    return distance_between;
}

bool zip_code_t::is_farther_than(const zip_code_t& zip1, const zip_code_t& zip2, const size_t& limit){
#ifdef DEBUG_ZIP_CODE
    cerr << "Checking if two zip codes are farther than " << limit << endl;
#endif

    size_t zip_index1 = 0; size_t zip_index2 = 0;
    size_t zip_value1 = std::numeric_limits<size_t>::max();
    size_t zip_value2 = std::numeric_limits<size_t>::max();

    //If the two positions aren't on the same connected component, then we're done
    std::tie(zip_value1, zip_index1) = zip1.zip_code.get_value_and_next_index(0);
    std::tie(zip_value2, zip_index2) = zip2.zip_code.get_value_and_next_index(0);
    if (zip_value1 != zip_value2) {
#ifdef DEBUG_ZIP_CODE
        cerr << "Zip codes are on different connected components" << endl;
#endif
        return true;
    }

    bool is_top_level_chain = zip_value1;
    std::tie(zip_value1, zip_index1) = zip1.zip_code.get_value_and_next_index(zip_index1);
    std::tie(zip_value2, zip_index2) = zip2.zip_code.get_value_and_next_index(zip_index2);
    if (zip_value1 != zip_value2) {
#ifdef DEBUG_ZIP_CODE
        cerr << "Zip codes are on different connected components" << endl;
#endif
        return true;
    }

    if (!is_top_level_chain) {
        //If the top-level thing is a snarl, then check if the zips are in the same chain. 
        //If they are, then proceed from the shared chain

        //The next thing will be the identifier for the chain
        std::tie(zip_value1, zip_index1) = zip1.zip_code.get_value_and_next_index(zip_index1);
        std::tie(zip_value2, zip_index2) = zip2.zip_code.get_value_and_next_index(zip_index2);
        if (zip_value1 != zip_value2) {
            //We can't tell
            return false;
        }
        //Next is the length of the chain
        std::tie(zip_value1, zip_index1) = zip1.zip_code.get_value_and_next_index(zip_index1);
        std::tie(zip_value2, zip_index2) = zip2.zip_code.get_value_and_next_index(zip_index2);
        if (zip_value1 < limit) {
            return true;
        }

        //The zips now point to the children of the shared chain, so we can proceed as if the top-level
        //structure was a chain

    }

    //Both zips now point to a thing in a shared chain
    //Get the minimum possible distance between the structures on the chain
    //For a lower bound, this assumes that the positions are as close as they can be on the structure in the chain
    size_t prefix_sum1, prefix_sum2, length1, length2;

    //The next thing could either be a snarl or a node. If it is a node, 
    vector<size_t> next_values;
    for (size_t i = 0 ; i < 3 ; i++ ) {
#ifdef DEBUG_ZIP_CODE
        assert(zip_index1 != std::numeric_limits<size_t>::max());
#endif
        std::tie(zip_value1, zip_index1) = zip1.zip_code.get_value_and_next_index(zip_index1);
        next_values.emplace_back(zip_value1);
    }
    if (zip_index1 == std::numeric_limits<size_t>::max()) {
#ifdef DEBUG_ZIP_CODE
        cerr << "zip1 is a node in a chain" << endl;
#endif
        //If the last thing was a node
        prefix_sum1 = next_values[0];
        length1 = next_values[1];
        prefix_sum1 = prefix_sum1 == 0 ? std::numeric_limits<size_t>::max() : prefix_sum1-1;
        length1 = length1 == 0 ? std::numeric_limits<size_t>::max() : length1-1;
    } else {
#ifdef DEBUG_ZIP_CODE
        cerr << "zip1 is in a snarl in a chain" << endl;
#endif
        //If the last thing was a snarl
        if (next_values[0]) {
            //If the next thing was a regular snarl
            prefix_sum1 = next_values[1];
            length1 = next_values[2];
            prefix_sum1 = prefix_sum1 == 0 ? std::numeric_limits<size_t>::max() : prefix_sum1-1;
            length1 = length1 == 0 ? std::numeric_limits<size_t>::max() : length1-1;
        } else {
            //If the next thing was an irregular snarl
            //TODO: If it's an irregular snarl, then we don't actually store the relevant values so we can't tell. Could look it up in the distance index or store it
            return false;
        }
    }
    
    //Do the same for the other zip
    next_values.clear();
    for (size_t i = 0 ; i < 3 ; i++ ) {
#ifdef DEBUG_ZIP_CODE
        assert(zip_index2 != std::numeric_limits<size_t>::max());
#endif
        std::tie(zip_value2, zip_index2) = zip2.zip_code.get_value_and_next_index(zip_index2);
        next_values.emplace_back(zip_value2);
    }
    if (zip_index2 == std::numeric_limits<size_t>::max()) {
#ifdef DEBUG_ZIP_CODE
        cerr << "zip2 is a node in a chain" << endl;
#endif
        //If the last thing was a node
        prefix_sum2 = next_values[0];
        length2 = next_values[1];
        prefix_sum2 = prefix_sum2 == 0 ? std::numeric_limits<size_t>::max() : prefix_sum2-1;
        length2 = length2 == 0 ? std::numeric_limits<size_t>::max() : length2-1;
    } else {
#ifdef DEBUG_ZIP_CODE
        cerr << "zip2 is in a snarl in a chain" << endl;
#endif
        //If the last thing was a snarl
        if (next_values[0]) {
            //If the next thing was a regular snarl
            prefix_sum2 = next_values[1];
            length2 = next_values[2];
            prefix_sum2 = prefix_sum2 == 0 ? std::numeric_limits<size_t>::max() : prefix_sum2-1;
            length2 = length2 == 0 ? std::numeric_limits<size_t>::max() : length2-1;
        } else {
            //If the next thing was an irregular snarl
            //TODO: If it's an irregular snarl, then we don't actually store the relevant values so we can't tell. Could look it up in the distance index or store it
            return false;
        }
    }
#ifdef DEBUG_ZIP_CODE
    cerr << "Finding distance in chain between " << prefix_sum1 << " " << length1 << " and " << prefix_sum2 << " and " << length2 << endl;
#endif

    if (prefix_sum1 == std::numeric_limits<size_t>::max() ||
        prefix_sum2 == std::numeric_limits<size_t>::max() ||
        length1 == std::numeric_limits<size_t>::max() ||
        length2 == std::numeric_limits<size_t>::max()) {
        //If anything is infinite, then we can't tell
        return false;
    }


    if (prefix_sum1 < prefix_sum2) {
        //If 1 comes first

        if (prefix_sum1 + length1 > prefix_sum2) {
            //They might be close
            return false;
        } else {
            //Return true if the distance between is greater than the limit
            return prefix_sum2 - (prefix_sum1 + length1) > limit; 
        }
    } else {
        //If 2 comes first

        if (prefix_sum2 + length2 > prefix_sum1) {
            //They might be close
            return false;
        } else {
            //Return true if the distance between is greater than the limit
            return prefix_sum1 - (prefix_sum2 + length2) > limit; 
        }
    }
}


}
