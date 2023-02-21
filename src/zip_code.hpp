#ifndef VG_ZIP_CODE_HPP_INCLUDED
#define VG_ZIP_CODE_HPP_INCLUDED

#include "varint.hpp"
#include "snarl_distance_index.hpp"
#include <gbwtgraph/minimizer.h>

namespace vg{
using namespace std;

//A decoded zip code as a vector of pair<is_chain, index>
//where is_chain indicates whether it's a chain/node, and index
//is the index of the node/snarl/chain code in the varint_vector_t
typedef std::vector<pair<bool, size_t>> zip_code_decoder_t;

enum code_type_t { NODE = 1, CHAIN, REGULAR_SNARL, IRREGULAR_SNARL, ROOT_SNARL, ROOT_CHAIN, ROOT_NODE};

/// A struct that represents a decoded node/snarl/chain code
/// Not all fields may be filled, code_type is always filled
/// node: length, rank_or_offset = prefix sum, is_reversed
/// chain: length, rank_or_offset = rank in snarl
/// regular_snarl: length, rank_or_offset = prefix sum, is_reversed (of the child)
/// irregular snarl: net_handle, length, rank_or_offset = prefix sum
/// root snarl: net_handle, rank or offset = connected component number
/// root chain: rank or offset = connected component number
/// root node: length, rank_or_offset = connected component number
struct decoded_code_t {
    net_handle_t net_handle;
    size_t length;
    size_t rank_or_offset;
    code_type_t code_type;
    bool is_reversed;

    /// Equality operator
    /// Do the two decoded_code_t's represent the same snarl tree node, assuming that all ancestors were the same
    /// All values must be the same, except for is_reversed in regular snarls, since this value refers to the 
    /// child of the regular snarl
    inline bool operator== (const decoded_code_t& other) const {
        return net_handle == net_handle &&
               length == other.length &&
               rank_or_offset == other.rank_or_offset &&
               code_type == other.code_type &&
               (code_type == REGULAR_SNARL || is_reversed == other.is_reversed);
    }

};

struct MIPayload;

/* Zip codes store the snarl decomposition location and distance information for a position on a graph
 * A zip code will contain all the information necessary to compute the minimum distance between two 
 * positions, with minimal queries to the distance index
 */

struct zip_code_t {

    public:
        typedef std::uint64_t code_type; // We assume that this fits into gbwtgraph::payload_type.

        //Constructor for a position and a distance index
        void fill_in_zip_code (const SnarlDistanceIndex& distance_index, const vg::pos_t& pos);

        //Get a decoder for interpreting the zip code
        zip_code_decoder_t decode() const;

        //Decode just one node/chain/snarl code given the index of its start in the varint_vector_t
        //And the code type of the actual code (ie, a chain if it is a trivial chain thats really a node)
        //It should be able to figure out what it is from NODE, SNARL, or CHAIN- if the index is
        //0 then it is assumed to be a root
        //It doesn't matter if it's given REGULAR or IRREGULAR_SNARL, the correct code type will be inferred from the actual code
        decoded_code_t decode_one_code(size_t index, const code_type_t& code_type, const SnarlDistanceIndex& distance_index) const;

        //Get the exact minimum distance between two positions and their zip codes
        static size_t minimum_distance_between(const zip_code_t& zip1, const pos_t& pos1, 
                                       const zip_code_t& zip2, const pos_t& pos2,
                                       const SnarlDistanceIndex& distance_index);

        //Return true if the minimum distance between the zip codes is definitely greater than limit
        //A false result is inconclusive
        static bool is_farther_than(const zip_code_t& zip1, const zip_code_t& zip2, const size_t& limit);


        //////////////////Functions to work with minimizer payloads for clustering 
        // Since we're sill using the old implementation, we need to be able to 
        // switch from zipcodes to payloads and back

        constexpr static gbwtgraph::payload_type NO_PAYLOAD = {0,0};

        //Encode zip code so it can be stored in the payload
        gbwtgraph::payload_type get_payload_from_zip() const;

        //Decode the zip code that got stored in the payload
        void fill_in_zip_code_from_payload(const gbwtgraph::payload_type& payload); 

        //This re-formats the new payload into the old payload format so it can be used 
        //for clustering
        gbwtgraph::payload_type get_old_payload_from_zipcode(const SnarlDistanceIndex& distance_index,
                                                                    const nid_t& id);


        size_t byte_count() const {
            return zip_code.byte_count();
        }

    //TODO: Make this private:
        varint_vector_t zip_code;


        /// Equality operator
        inline bool operator== (const zip_code_t& other) const {
            return zip_code == other.zip_code;
    }

    private:

        /* Functions for getting the zip code for each snarl/chain/node
         * Distances will be stored as distance+1, 0 will be reserved for inf
         */
        //Return a vector of size_ts that will represent the node in the zip code
        inline vector<size_t> get_node_code(const net_handle_t& node, const SnarlDistanceIndex& distance_index);
        //Return a vector of size_ts that will represent the chain in the zip code
        inline vector<size_t> get_chain_code(const net_handle_t& chain, const SnarlDistanceIndex& distance_index);
        //Return a vector of size_ts that will represent the snarl in the zip code
        inline vector<size_t> get_regular_snarl_code(const net_handle_t& snarl, const net_handle_t& snarl_child, 
                                                            const SnarlDistanceIndex& distance_index);
        //Return a vector of size_ts that will represent the snarl in the zip code
        inline vector<size_t> get_irregular_snarl_code(const net_handle_t& snarl, const SnarlDistanceIndex& distance_index);
};



/**
    The payload for the minimizer index. This stores distance information that gets used in clustering
    The payload now uses zip codes, so this gets used to go from a zip code to distance information
    usable by the clusterer, which expects the old payload format
*/
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

    
}; 
}

#endif
