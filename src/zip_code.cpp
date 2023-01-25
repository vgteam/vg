#include "zip_code.hpp"

#define DEBUG_ZIP_CODE

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
        cerr << "Adding code for top-level snarl" << endl;
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
    cerr << "ADDING NODE CODE " << node_code[0] << " " << node_code[1] << " " << node_code[2] << endl;
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
}
