#ifndef VG_SNARL_DISTANCE_HPP_INCLUDED
#define VG_SNARL_DISTANCE_HPP_INCLUDED

#include <bdsg/snarl_distance_index.hpp>
#include "snarls.hpp"
#include <structures/union_find.hpp>
#include "hash_map.hpp"


namespace vg { 

using namespace sdsl;
using namespace handlegraph;
using namespace bdsg;
//Fill in the index
void fill_in_distance_index(SnarlDistanceIndex* distance_index, const HandleGraph* graph, const HandleGraphSnarlFinder* snarl_finder, size_t size_limit = 500);

//Fill in the temporary snarl record with distances
void populate_snarl_index(SnarlDistanceIndex::TemporaryDistanceIndex& temp_index, 
    pair<SnarlDistanceIndex::temp_record_t, size_t> snarl_index, size_t size_limit, const HandleGraph* graph) ;

SnarlDistanceIndex::TemporaryDistanceIndex make_temporary_distance_index(const HandleGraph* graph, const HandleGraphSnarlFinder* snarl_finder, size_t size_limit);

//Define wang_hash for net_handle_t's so that we can use a hash_map
template<> struct wang_hash<handlegraph::net_handle_t> {
public:
    inline size_t operator()(const net_handle_t& net_handle) const {
        return wang_hash_64(as_integer(net_handle));
    }
};

//TODO: This is outdated
//Given a position, return distances that can be stored by a minimizer
//
//If the position is on a boundary node of a top level chain, then return true, and 
//a unique identifier for the connected component that the node is on and
//the offset of the position in the root chain - the minimum distance from the beginning of the chain to 
//the position
//The second bool will be false and the remaining size_t's will be 0
//
//If the position is on a child node of a top-level simple bubble (bubble has no children and nodes connect only to boundaries)
//return false, 0, 0, true, and the rank of the bubble in its chain, the length of the start
//node of the snarl, the length of the end node (relative to a fd traversal of the chain), and
//the length of the node
//
//If the position is not on a root node (that is, a boundary node of a snarl in a root chain), returns
//false and MIPayload::NO_VALUE for all values
//


//Given a position, return the distances that can be stored by a minimizer
//node length, connected component, prefix sum, chain component, is-reversed 
tuple<size_t, size_t, size_t, size_t, bool>  get_minimizer_distances (const SnarlDistanceIndex& distance_index, pos_t pos);



//Given an alignment to a graph and a range, find the set of nodes in the
//graph for which the minimum distance from the position to any position
//in the node is within the given distance range
//If look_forward is true, then start from the start of the path forward,
//otherwise start from the end going backward
void subgraph_in_distance_range(const SnarlDistanceIndex& distance_index, const Path& path, const HandleGraph* super_graph, size_t min_distance,
                                size_t max_distance, std::unordered_set<nid_t>& subgraph, bool look_forward);
///Helper for subgraph_in_distance_range
///Given starting handles in the super graph and the distances to each handle (including the start position and
//the first position in the handle), add all nodes within the distance range, excluding nodes in seen_nodes
void add_nodes_in_distance_range(const HandleGraph* super_graph, size_t min_distance, size_t max_distance,
                        std::unordered_set<nid_t>& subgraph, vector<pair<handle_t, size_t>>& start_nodes,
                        hash_set<pair<nid_t, bool>>& seen_nodes);


//Add nodes to the subgraph if they are near the path in the snarl tree
//Walks up the snarl tree from either end of the path, then takes everything
//in between them
void subgraph_containing_path_snarls(const SnarlDistanceIndex& distance_index, const HandleGraph* graph, const Path& path, std::unordered_set<nid_t>& subgraph);


//Helper function for subgraph_containing_path_snarls
//Add all the nodes in the parent to the subgraph
void add_descendants_to_subgraph(const SnarlDistanceIndex& distance_index, const net_handle_t& parent, std::unordered_set<nid_t>& subgraph);



///**
//TODO: I'm not doing this anymore
// * The encoding of distances for positions in top-level chains or top-level simple bubbles.
// * Either stores (chain id, chain offset) for a position on a top-level chain, or
// * (snarl rank, node length, start length, end length) for a position on a simple bubble
// * We store this information in the minimizer index.
// */
///*
//Simple bubble: 
//    
// 8 bit  |     1    |        24           |    10     |     10   |    10     |    1
//  ---   |  is rev  | snarl rank in chain | start len | end len  | node len  |  is_node
//   
//Top level chain 
//     
//    31 bit   |    32    |     1
//component id |  offset  |  is_node
//
//
//is_node is true if it is a top-level chain node, false if it is a simple bubble
//*/

/**
 * The encoding of distances for positions
 * Stores 
 *     11 bit   |            13         |     31       |         8         |      1
 *  node length |  connected component  |  prefix sum  |  chain_component  | is_reversed
 *
 * We store this information in the minimizer index.
 */
struct MIPayload {
    typedef std::uint64_t code_type; // We assume that this fits into gbwtgraph::payload_type.

    constexpr static code_type NO_CODE = std::numeric_limits<code_type>::max();
    constexpr static std::size_t NO_VALUE = std::numeric_limits<size_t>::max(); 

    const static size_t CHAIN_COMPONENT_OFFSET = 1;
    const static size_t CHAIN_COMPONENT_WIDTH = 8;
    const static code_type CHAIN_COMPONENT_MASK = (static_cast<code_type>(1) << CHAIN_COMPONENT_WIDTH) - 1;

    const static size_t PREFIX_SUM_OFFSET = 9;
    const static size_t PREFIX_SUM_WIDTH = 31;
    const static code_type PREFIX_SUM_MASK = (static_cast<code_type>(1) << PREFIX_SUM_WIDTH) - 1;

    const static size_t COMPONENT_OFFSET = 40;
    const static size_t COMPONENT_WIDTH = 13;
    const static code_type COMPONENT_MASK = (static_cast<code_type>(1) << COMPONENT_WIDTH) - 1;

    const static size_t NODE_LENGTH_OFFSET = 53;
    const static size_t NODE_LENGTH_WIDTH = 11;
    const static code_type NODE_LENGTH_MASK = (static_cast<code_type>(1) << NODE_LENGTH_WIDTH) - 1;



    static code_type encode(tuple<size_t, size_t, size_t, size_t, bool> info) {

        size_t node_length = std::get<0>(info);
        size_t component = std::get<1>(info);
        size_t prefix_sum = std::get<2>(info);
        size_t chain_component = std::get<3>(info);
        bool is_reversed = std::get<4>(info);

        if ( node_length > NODE_LENGTH_MASK 
             || component > COMPONENT_MASK
             || prefix_sum > PREFIX_SUM_MASK
             || chain_component > CHAIN_COMPONENT_MASK) {
            //If there aren't enough bits to represent one of the values
            return NO_CODE;
        }

        return (static_cast<code_type>(node_length) << NODE_LENGTH_OFFSET) 
             | (static_cast<code_type>(component) << COMPONENT_OFFSET)
             | (static_cast<code_type>(prefix_sum) << PREFIX_SUM_OFFSET)
             | (static_cast<code_type>(chain_component) << CHAIN_COMPONENT_OFFSET)
             | is_reversed;

    }

    static std::tuple<size_t, size_t, size_t, size_t, bool> decode(code_type code) {
        if (code == NO_CODE) {
            return make_tuple(NO_VALUE, NO_VALUE, NO_VALUE, NO_VALUE, false);
        } else {

            return std::tuple<size_t, size_t, size_t, size_t, bool> (
                code >> NODE_LENGTH_OFFSET & NODE_LENGTH_MASK,
                code >> COMPONENT_OFFSET & COMPONENT_MASK,
                code >> PREFIX_SUM_OFFSET & PREFIX_SUM_MASK,
                code >> CHAIN_COMPONENT_OFFSET & CHAIN_COMPONENT_MASK,
                code & 1);
        }
    }
};

}

#endif
