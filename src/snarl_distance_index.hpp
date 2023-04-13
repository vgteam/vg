#ifndef VG_SNARL_DISTANCE_HPP_INCLUDED
#define VG_SNARL_DISTANCE_HPP_INCLUDED

#include <bdsg/snarl_distance_index.hpp>
#include "snarls.hpp"
#include <structures/union_find.hpp>
#include "hash_map.hpp"
#include <gbwtgraph/minimizer.h>


namespace vg { 

using namespace sdsl;
using namespace handlegraph;
using namespace bdsg;

//Minimum distance taking a pos instead of id/orientation/offset
size_t minimum_distance(const SnarlDistanceIndex& distance_index, pos_t pos1, pos_t pos2,
                        bool unoriented_distance = false, const HandleGraph* graph=nullptr); 
//Maximum distance taking a pos instead of id/orientation/offset
size_t maximum_distance(const SnarlDistanceIndex& distance_index, pos_t pos1, pos_t pos2); 

//Fill in the index
//size_limit is a limit on the number of nodes in a snarl, after which the index won't store pairwise distances
void fill_in_distance_index(SnarlDistanceIndex* distance_index, const HandleGraph* graph, const HandleGraphSnarlFinder* snarl_finder, size_t size_limit = 10000);

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
//traversal_start is the node that we started the search from, because we can traverse it a second time but 
//we don't want to include a loop distance to any node after it
void subgraph_in_distance_range_walk_graph(const HandleGraph* super_graph, size_t min_distance, size_t max_distance,
                        std::unordered_set<nid_t>& subgraph, vector<pair<handle_t, size_t>>& start_nodes,
                        hash_set<pair<nid_t, bool>>& seen_nodes, const pair<nid_t, bool>& traversal_start);

//Helper function for subgraph_in_distance_range
//Given a node that is a child of a chain and the distance to it, walk forward from the
//current child and add nodes to search_start_nodes
void subgraph_in_distance_range_walk_across_chain (const SnarlDistanceIndex& distance_index, 
        const HandleGraph* super_graph,std::unordered_set<nid_t>& subgraph,
        net_handle_t current_node, size_t current_distance, 
        vector<pair<handle_t, size_t>>& search_start_nodes,
        hash_set<pair<nid_t, bool>>& seen_nodes, 
        const size_t& min_distane, const size_t& max_distance, bool checked_loop=false);


//Add nodes to the subgraph if they are near the path in the snarl tree
//Walks up the snarl tree from either end of the path, then takes everything
//in between them
void subgraph_containing_path_snarls(const SnarlDistanceIndex& distance_index, const HandleGraph* graph, const Path& path, std::unordered_set<nid_t>& subgraph);


//Helper function for subgraph_containing_path_snarls
//Add all the nodes in the parent to the subgraph
void add_descendants_to_subgraph(const SnarlDistanceIndex& distance_index, const net_handle_t& parent, std::unordered_set<nid_t>& subgraph);



//The distance values that get stored in an MIPayload
struct MIPayloadValues{

    //The record offset of the node
    size_t record_offset; 

    //The record offset of the parent
    size_t parent_record_offset; 

    //The node record offset of the node (eg, which node in a trivial snarl)
    size_t node_record_offset;

    size_t node_length;

    //Is the node reversed in its parent
    bool is_reversed;

    bool is_trivial_chain;

    bool parent_is_chain;

    bool parent_is_root; 

    size_t prefix_sum;

    size_t chain_component;
};

///
// The encoding of distances for positions in top-level chains
// We store this information in the minimizer index.
// 
// This gets stored in two separate uint64_t's
//
//           32 bits       |           32
//   record offset of node |  record offset of parent  
//     
//    8 bits          |    12 bit   |      1      |          1       |     1           |        1       |   32       |         8        
// node record offset | node length | is_reversed | is trivial chain | parent is chain | parent is root | prefix sum |  chain_component 
//
//
// These values are en/de-coded from the raw values in the order above
//
// If no values are stored, then the two uint64_t's will both be inf
// bools are always stored, everything else is all 1's if it is not stored
// 

struct MIPayload {
    typedef std::uint64_t code_type; // We assume that this fits into gbwtgraph::payload_type.
    //typedef std::pair<code_type, code_type> payload_type;

    
    constexpr static gbwtgraph::payload_type NO_CODE = {0, 0};
    constexpr static std::size_t NO_VALUE = std::numeric_limits<size_t>::max(); 


    //Static values for the offset from the right side of the uint64_t storing the values, the width of each value, and a bit mask for the value
    const static size_t PARENT_RECORD_OFFSET = 0;
    const static size_t PARENT_RECORD_WIDTH = 32;
    const static code_type PARENT_RECORD_MASK = (static_cast<code_type>(1) << PARENT_RECORD_WIDTH) - 1;

    const static size_t NODE_RECORD_OFFSET = 32;
    const static size_t NODE_RECORD_WIDTH = 32;
    const static code_type NODE_RECORD_MASK = (static_cast<code_type>(1) << NODE_RECORD_WIDTH) - 1;


    const static size_t CHAIN_COMPONENT_OFFSET = 0;
    const static size_t CHAIN_COMPONENT_WIDTH = 8;
    const static code_type CHAIN_COMPONENT_MASK = (static_cast<code_type>(1) << CHAIN_COMPONENT_WIDTH) - 1;
    
    const static size_t PREFIX_SUM_OFFSET = 8;
    const static size_t PREFIX_SUM_WIDTH = 32;
    const static code_type PREFIX_SUM_MASK = (static_cast<code_type>(1) << PREFIX_SUM_WIDTH) - 1;

    const static size_t PARENT_IS_ROOT_OFFSET = 40;
    const static size_t PARENT_IS_CHAIN_OFFSET = 41;
    const static size_t IS_TRIVIAL_CHAIN_OFFSET = 42;
    const static size_t IS_REVERSED_OFFSET = 43;
    
    const static size_t NODE_LENGTH_OFFSET = 44;
    const static size_t NODE_LENGTH_WIDTH = 12;
    const static code_type NODE_LENGTH_MASK = (static_cast<code_type>(1) << NODE_LENGTH_WIDTH) - 1;
    
    const static size_t NODE_RECORD_OFFSET_OFFSET = 56;
    const static size_t NODE_RECORD_OFFSET_WIDTH = 8;
    const static code_type NODE_RECORD_OFFSET_MASK = (static_cast<code_type>(1) << NODE_RECORD_OFFSET_WIDTH) - 1;

    //Encode and decode from the following values:
    //record offset of node, record offset of parent, node record offset, node length, is_reversed, parent is chain, prefix sum, chain_component 
    static gbwtgraph::payload_type encode(MIPayloadValues info) {

        if ( info.record_offset > NODE_RECORD_MASK 
             || info.parent_record_offset > PARENT_RECORD_MASK
             || info.node_record_offset > NODE_RECORD_OFFSET_MASK
             || info.node_length > NODE_LENGTH_MASK
             || info.prefix_sum > PREFIX_SUM_MASK
             || info.chain_component > CHAIN_COMPONENT_MASK) {
            //If there aren't enough bits to represent one of the values
            return NO_CODE;
        }

        code_type encoded1 = (static_cast<code_type>(info.record_offset)           << NODE_RECORD_OFFSET)
                           | (static_cast<code_type>(info.parent_record_offset)  << PARENT_RECORD_OFFSET);

        code_type encoded2 = (static_cast<code_type>(info.node_record_offset) << NODE_RECORD_OFFSET_OFFSET)
                           | (static_cast<code_type>(info.node_length)        << NODE_LENGTH_OFFSET)
                           | (static_cast<code_type>(info.is_reversed)        << IS_REVERSED_OFFSET)
                           | (static_cast<code_type>(info.is_trivial_chain)   << IS_TRIVIAL_CHAIN_OFFSET)
                           | (static_cast<code_type>(info.parent_is_chain)    << PARENT_IS_CHAIN_OFFSET)
                           | (static_cast<code_type>(info.parent_is_root)    << PARENT_IS_ROOT_OFFSET)
                           | (static_cast<code_type>(info.prefix_sum)         << PREFIX_SUM_OFFSET)
                           | (static_cast<code_type>(info.chain_component)    << CHAIN_COMPONENT_OFFSET);

        return {encoded1, encoded2};

    }

    //Set the values of a code. Mutate the given code 
    static void set_record_offset(gbwtgraph::payload_type& code, size_t record_offset) {
        //Set everything in node_record slot to 0's
        code.first = code.first & ~(NODE_RECORD_MASK << NODE_RECORD_OFFSET); 
        //And | with the value to set it
        code.first = code.first | (static_cast<code_type>(record_offset) << NODE_RECORD_OFFSET); 
    }
    static void set_parent_record_offset(gbwtgraph::payload_type& code, size_t parent_record_offset) {
        code.first = code.first & ~(PARENT_RECORD_MASK << PARENT_RECORD_OFFSET); 
        code.first = code.first | (static_cast<code_type>(parent_record_offset) << PARENT_RECORD_OFFSET); 
    }
    static void set_node_record_offset(gbwtgraph::payload_type& code, size_t node_record_offset) {
        code.second = code.second & ~(NODE_RECORD_OFFSET_MASK << NODE_RECORD_OFFSET_OFFSET);
        code.second = code.second | (static_cast<code_type>(node_record_offset) << NODE_RECORD_OFFSET_OFFSET);
    }
    static void set_node_length(gbwtgraph::payload_type& code, size_t node_length) {
        code.second = code.second & ~(NODE_LENGTH_MASK << NODE_LENGTH_OFFSET);
        code.second = code.second | (static_cast<code_type>(node_length) << NODE_LENGTH_OFFSET);
    }
    static void set_is_reversed(gbwtgraph::payload_type& code, bool is_reversed) {
        code.second = code.second & ~(static_cast<code_type>(1) << IS_REVERSED_OFFSET);
        code.second = code.second | (static_cast<code_type>(is_reversed) << IS_REVERSED_OFFSET);
    }
    static void set_is_trivial_chain(gbwtgraph::payload_type& code, bool is_trivial_chain) {
        code.second = code.second & ~(static_cast<code_type>(1) << IS_TRIVIAL_CHAIN_OFFSET);
        code.second = code.second | (static_cast<code_type>(is_trivial_chain)   << IS_TRIVIAL_CHAIN_OFFSET);
    }
    static void set_parent_is_chain(gbwtgraph::payload_type& code, bool parent_is_chain) {
        code.second = code.second & ~(static_cast<code_type>(1) << PARENT_IS_CHAIN_OFFSET);
        code.second = code.second | (static_cast<code_type>(parent_is_chain) << PARENT_IS_CHAIN_OFFSET);
    }
    static void set_parent_is_root(gbwtgraph::payload_type& code, bool parent_is_root) {
        code.second = code.second & ~(static_cast<code_type>(1) << PARENT_IS_ROOT_OFFSET);
        code.second = code.second | (static_cast<code_type>(parent_is_root) << PARENT_IS_ROOT_OFFSET);
    }
    static void set_prefix_sum(gbwtgraph::payload_type& code, size_t prefix_sum) {
        code.second = code.second & ~(PREFIX_SUM_MASK << PREFIX_SUM_OFFSET);
        code.second = code.second | (static_cast<code_type>(prefix_sum) << PREFIX_SUM_OFFSET);
    }
    static void set_chain_component(gbwtgraph::payload_type& code, size_t chain_component) {
        code.second = code.second & ~(CHAIN_COMPONENT_MASK << CHAIN_COMPONENT_OFFSET);
        code.second = code.second | (static_cast<code_type>(chain_component) << CHAIN_COMPONENT_OFFSET);
    }


    //How do decode the code
    static size_t record_offset(const gbwtgraph::payload_type code) { 
        if (code == NO_CODE) {
            return NO_VALUE;
        }
        return (size_t) (code.first  >> NODE_RECORD_OFFSET & NODE_RECORD_MASK);
    }
    static size_t parent_record_offset(const gbwtgraph::payload_type code) {
        if (code == NO_CODE) {
            return NO_VALUE;
        } 
        return (size_t) (code.first  >> PARENT_RECORD_OFFSET & PARENT_RECORD_MASK);
    }

    static size_t node_record_offset(const gbwtgraph::payload_type code) { 
        if (code == NO_CODE) {
            return NO_VALUE;
        }
        return (size_t) (code.second >> NODE_RECORD_OFFSET_OFFSET & NODE_RECORD_OFFSET_MASK);
    }
    static size_t node_length(const gbwtgraph::payload_type code) { 
        if (code == NO_CODE) {
            return NO_VALUE;
        }
        return (size_t) (code.second >> NODE_LENGTH_OFFSET & NODE_LENGTH_MASK);
    }
    static bool is_reversed(const gbwtgraph::payload_type code) {
        if (code == NO_CODE) {
            return false;
        }
        return (bool) (code.second >> IS_REVERSED_OFFSET & 1);
    }
    static bool is_trivial_chain (const gbwtgraph::payload_type code) { 
        if (code == NO_CODE) {
            return false;
        }
        return (bool) (code.second >> IS_TRIVIAL_CHAIN_OFFSET   & 1);
    }
    static bool parent_is_chain(const gbwtgraph::payload_type code) { 
        if (code == NO_CODE) {
            return false;
        }
        return (bool) (code.second >> PARENT_IS_CHAIN_OFFSET    & 1);
    }
    static bool parent_is_root (const gbwtgraph::payload_type code) { 
        if (code == NO_CODE) {
            return false;
        }
        return (bool) (code.second >> PARENT_IS_ROOT_OFFSET     & 1);
    }
    static size_t prefix_sum (const gbwtgraph::payload_type code) { 
        if (code == NO_CODE) {
            return NO_VALUE;
        }
        return (size_t) (code.second >> PREFIX_SUM_OFFSET & PREFIX_SUM_MASK);
    }
    static size_t chain_component (const gbwtgraph::payload_type code) { 
        if (code == NO_CODE) {
            return NO_VALUE;
        }
        return (size_t) (code.second >> CHAIN_COMPONENT_OFFSET    & CHAIN_COMPONENT_MASK);
    }

    

    static MIPayloadValues decode(gbwtgraph::payload_type code) {
        if (code == NO_CODE) {
            return {NO_VALUE, NO_VALUE, NO_VALUE, NO_VALUE, false, false, false, false, NO_VALUE, NO_VALUE};
        } else {
            return {
                record_offset(code),
                parent_record_offset(code),
                node_record_offset(code),
                node_length(code),
                is_reversed(code),
                is_trivial_chain(code),
                parent_is_chain(code),
                parent_is_root(code),
                prefix_sum(code),
                chain_component(code)};


        }
    }

};

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
//record offset of node, record offset of parent, node record offset, node length, is_reversed, is_trivial_chain, parent is chain, prefix sum, chain_component 
MIPayloadValues get_minimizer_distances (const SnarlDistanceIndex& distance_index, pos_t pos);



}

#endif
