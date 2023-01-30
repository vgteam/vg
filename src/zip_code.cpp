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

vector<size_t> zip_code_t::get_node_code(const net_handle_t& node, const SnarlDistanceIndex& distance_index) {
#ifdef DEBUG_ZIP_CODE
    assert(!distance_index.is_trivial_chain(node));
    assert((distance_index.is_chain(distance_index.get_parent(node)) || distance_index.is_root(distance_index.get_parent(node))));
#endif
    //Node code is: offset in chain, length, is reversed
    vector<size_t> node_code;
    //Assume this node is in a regular chain
    node_code.emplace_back(distance_index.get_prefix_sum_value(node));
    node_code.emplace_back(distance_index.minimum_length(node));
    node_code.emplace_back(distance_index.is_reversed_in_parent(node));
    return node_code;

}
vector<size_t> zip_code_t::get_chain_code(const net_handle_t& chain, const SnarlDistanceIndex& distance_index) {
    //Chain code is: rank in snarl, length
    vector<size_t> chain_code;
    chain_code.emplace_back(distance_index.get_rank_in_parent(chain));
    chain_code.emplace_back(distance_index.minimum_length(chain));
    return chain_code;

}
vector<size_t> zip_code_t::get_regular_snarl_code(const net_handle_t& snarl, const net_handle_t& snarl_child, const SnarlDistanceIndex& distance_index) {
    //Regular snarl code is 1, offset in chain, length, is reversed
    vector<size_t> snarl_code;

    //Tag to say that it's a regular snarl
    snarl_code.emplace_back(1);

    //Chain prefix sum value for the start of the snarl, which is the prefix sum of the start node + length of the start node
    net_handle_t start_node = distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, false, false));
    snarl_code.emplace_back(distance_index.get_prefix_sum_value(start_node) + distance_index.minimum_length(start_node));

    //Length of the snarl
    snarl_code.emplace_back(distance_index.minimum_length(snarl));

    //Is the child of the snarl reversed in the snarl
#ifdef DEBUG_ZIP_CODE
    assert(distance_index.is_chain(snarl_child));
#endif
    snarl_code.emplace_back(distance_index.is_reversed_in_parent(snarl_child));

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
    size_t zip_index1 = 0; size_t zip_index2 = 0;
    size_t zip_value1 = std::numeric_limits<size_t>::max();
    size_t zip_value2 = std::numeric_limits<size_t>::max();

    //If the two positions aren't on the same connected component, then we're done
    std::tie(zip_value1, zip_index1) = zip1.zip_code.get_value_and_next_index(0);
    std::tie(zip_value2, zip_index2) = zip2.zip_code.get_value_and_next_index(0);
    if (zip_value1 != zip_value2) {
        return std::numeric_limits<size_t>::max();
    }
    std::tie(zip_value1, zip_index1) = zip1.zip_code.get_value_and_next_index(zip_index1);
    std::tie(zip_value2, zip_index2) = zip2.zip_code.get_value_and_next_index(zip_index2);
    if (zip_value1 != zip_value2) {
        return std::numeric_limits<size_t>::max();
    }

    //The two positions are in the same connected component so now try to find the distance
    zip_code_decoder_t decoded_zip1 = zip1.decode();
    zip_code_decoder_t decoded_zip2 = zip2.decode();

    return std::numeric_limits<size_t>::max();


}

}
