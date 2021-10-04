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
        return wang_hash_64(reinterpret_cast<const size_t>(as_integer(net_handle)));
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
tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool>  get_minimizer_distances (const SnarlDistanceIndex& distance_index, pos_t pos);


//Empty distances so that the minimizers won't cache anything
tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool>  get_empty_minimizer_distances ();



/**
 * The encoding of distances for positions in top-level chains or top-level simple bubbles.
 * Either stores (chain id, chain offset) for a position on a top-level chain, or
 * (snarl rank, node length, start length, end length) for a position on a simple bubble
 * We store this information in the minimizer index.
 */
/*
Simple bubble: 
    
 8 bit  |     1    |        24           |    10     |     10   |    10     |    1
  ---   |  is rev  | snarl rank in chain | start len | end len  | node len  |  is_node
   
Top level chain 
     
    31 bit   |    32    |     1
component id |  offset  |  is_node


is_node is true if it is a top-level chain node, false if it is a simple bubble
*/

struct MIPayload {
    typedef std::uint64_t code_type; // We assume that this fits into gbwtgraph::payload_type.

    constexpr static code_type NO_CODE = std::numeric_limits<code_type>::max();
    constexpr static size_t NO_VALUE = std::numeric_limits<size_t>::max(); // From offset_in_root_chain().

    constexpr static size_t NODE_LEN_OFFSET = 1;
    constexpr static size_t END_LEN_OFFSET = 11;
    constexpr static size_t START_LEN_OFFSET = 21;
    constexpr static size_t RANK_OFFSET = 31;
    constexpr static size_t REV_OFFSET = 55;

    
    constexpr static size_t LENGTH_WIDTH = 10;
    constexpr static size_t RANK_WIDTH = 24;
    constexpr static code_type LENGTH_MASK = (static_cast<code_type>(1) << LENGTH_WIDTH) - 1;
    constexpr static code_type RANK_MASK = (static_cast<code_type>(1) << RANK_WIDTH) - 1;
    

    
    constexpr static size_t ID_OFFSET = 33;
    constexpr static size_t ID_WIDTH = 31;
    constexpr static size_t OFFSET_WIDTH = 32;
    constexpr static code_type OFFSET_MASK = (static_cast<code_type>(1) << OFFSET_WIDTH) - 1;

    static code_type encode(std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_pos) {
        bool is_top_level_node; size_t parent_offset; size_t offset_in_parent; //Values for a top level chain
        bool is_top_level_snarl; size_t snarl_rank; size_t node_length; size_t start_length; size_t end_length; bool is_rev;//values for a bubble
        std::tie(is_top_level_node, parent_offset, offset_in_parent, is_top_level_snarl, snarl_rank, start_length, end_length, node_length, is_rev) = chain_pos;

        if (!is_top_level_node && ! is_top_level_snarl) {

            return NO_CODE;

        } else if (is_top_level_node) {
            //Top level node in chain

            if (parent_offset >= (static_cast<code_type>(1) << 31) - 1 
                || offset_in_parent >= static_cast<size_t>(OFFSET_MASK) ) {
                //If the values are too large to be stored
                return NO_CODE;
            }

            return (parent_offset << ID_OFFSET) | (offset_in_parent << 1) | static_cast<code_type>(true);

        } else {
            //Top level simple bubble

            if (snarl_rank >= static_cast<size_t>(RANK_MASK)
                || start_length >= static_cast<size_t>(LENGTH_MASK)
                || end_length >=  static_cast<size_t>(LENGTH_MASK)
                || node_length >= static_cast<size_t>(LENGTH_MASK) ){
                //If the values are too large to be stored
                return NO_CODE;
            }

            return (static_cast<code_type>(is_rev) << REV_OFFSET) | (snarl_rank << RANK_OFFSET) | (start_length << START_LEN_OFFSET) | (end_length << END_LEN_OFFSET) | (node_length << NODE_LEN_OFFSET) ;
        }
    }

    static std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> decode(code_type code) {
        if (code == NO_CODE) {
            return std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool>(false, NO_VALUE, NO_VALUE, false, NO_VALUE, NO_VALUE, NO_VALUE, NO_VALUE, false);
        } else if ((code & (static_cast<code_type>(1))) == (static_cast<code_type>(1))) {
            //This is a top-level chain
            return std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool>
                    (true, code >> ID_OFFSET, code >> 1 & OFFSET_MASK, false, NO_VALUE, NO_VALUE, NO_VALUE, NO_VALUE, false);
        } else {
            //This is a top-level bubble
            return std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool>
                    (false, NO_VALUE, NO_VALUE, true,
                      code >> RANK_OFFSET & RANK_MASK, 
                      code >> START_LEN_OFFSET & LENGTH_MASK, 
                      code >> END_LEN_OFFSET & LENGTH_MASK, 
                      code >> NODE_LEN_OFFSET & LENGTH_MASK,
                      code >> REV_OFFSET & static_cast<code_type>(1));
        }
    }
};

}

#endif
