#ifndef VG_ZIP_CODE_HPP_INCLUDED
#define VG_ZIP_CODE_HPP_INCLUDED

#include "varint.hpp"
#include "snarl_distance_index.hpp"
#include <gbwtgraph/minimizer.h>

namespace vg{
using namespace std;

//A decoder for interpreting a zipcode
//Can interpret the values for a snarl tree node given the depth (index into the vector)
struct zipcode_decoder_t;

enum code_type_t { NODE = 1, CHAIN, REGULAR_SNARL, IRREGULAR_SNARL, ROOT_SNARL, ROOT_CHAIN, ROOT_NODE};

///A struct to interpret the minimizer payload
///I want to use zipcodes as the payload but at the moment clustering still expects the old payload
///This can interpret zipcodes to format them as the old payload
struct MIPayload;

/* Zip codes store the snarl decomposition location and distance information for a position on a graph
 * A zip code will contain all the information necessary to compute the minimum distance between two 
 * positions, with minimal queries to the distance index
 */

struct zipcode_t {

    public:
        typedef std::uint64_t code_type; // We assume that this fits into gbwtgraph::payload_type.

        //Constructor for a position and a distance index
        void fill_in_zipcode (const SnarlDistanceIndex& distance_index, const vg::pos_t& pos);

        //Get the exact minimum distance between two positions and their zip codes
        static size_t minimum_distance_between(const zipcode_t& zip1, const pos_t& pos1, 
                                       const zipcode_t& zip2, const pos_t& pos2,
                                       const SnarlDistanceIndex& distance_index, bool directed_distance=true);

        //Return true if the minimum distance between the zip codes is definitely greater than limit
        //A false result is inconclusive
        static bool is_farther_than(const zipcode_t& zip1, const zipcode_t& zip2, const size_t& limit);

        //Get a tuple of the top-level structure id, prefix sum of the child of the top-level chain, and 
        //the length of the child of the top-level chain
        //This gets used to quickly compare the two zip codes for is_farther_than
        static tuple<size_t, size_t, size_t> get_top_level_chain_offset();


        //////////////////Functions to work with minimizer payloads for clustering 
        // Since we're sill using the old implementation, we need to be able to 
        // switch from zipcodes to payloads and back

        //Encode zip code so it can be stored in the payload
        gbwtgraph::payload_type get_payload_from_zip() const;

        //Decode the zip code that got stored in the payload
        void fill_in_zipcode_from_payload(const gbwtgraph::payload_type& payload); 


        size_t byte_count() const {
            return zipcode.byte_count();
        }

    //TODO: Make this private:
        varint_vector_t zipcode;


        /// Equality operator
        inline bool operator== (const zipcode_t& other) const {
            return zipcode == other.zipcode;
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
    friend class zipcode_decoder_t;
};


///A struct for decoding a zipcode
struct zipcode_decoder_t {

    ///The decoder as a vector of pair<is_chain, index>, one for each snarl tree node in the zip
    ///where is_chain indicates whether it's a chain/node, and index
    ///is the index of the node/snarl/chain code in the varint_vector_t
    std::vector<pair<bool, size_t>> decoder;

    ///The zipcode that this is decoding
    const zipcode_t* zipcode;


    ///Constructor that goes through the zipcode and decodes it to fill in decoder
    ///If a depth is given, then only fill in up to depth snarl tree nodes
    ///Otherwise, fill in the whole zipcode
    zipcode_decoder_t(const zipcode_t* zipcode, const size_t& depth=std::numeric_limits<size_t>::max());

    //Go through the entire zipcode and fill in the decoder
    void fill_in_full_decoder();

    //Fill in one more item in the decoder
    //Returns true if this is the last thing in the zipcode and false if there is more to decode
    bool fill_in_next_decoder();
    ///What type of snarl tree node is at the given depth (index into the zipcode)
    code_type_t get_code_type(const size_t& depth) ;

    ///Get the length of a snarl tree node given the depth in the snarl tree
    ///This requires the distance index for irregular snarls (except for a top-level snarl)
    ///Throws an exception if the distance index is not given when it is needed
    ///Doesn't use a given distance index if it isn't needed
    size_t get_length(const size_t& depth, const SnarlDistanceIndex* distance_index=nullptr) ;

    ///Get the rank of a node/snarl in a snarl. Throw an exception if it isn't the child of a snarl
    size_t get_rank_in_snarl(const size_t& depth) ;

    ///Get the prefix sum of a child of a chain
    ///This requires the distance index for irregular snarls (except for a top-level snarl)
    ///Throws an exception if the distance index is not given when it is needed
    ///Doesn't use a given distance index if it isn't needed
    size_t get_offset_in_chain(const size_t& depth, const SnarlDistanceIndex* distance_index=nullptr) ;

    ///Get the handle of the thing at the given depth. This can only be used for
    ///Root-level structures or irregular snarls
    bool get_is_reversed_in_parent(const size_t& depth);

    ///Get the handle of the thing at the given depth. This can only be used for
    ///Root-level structures or irregular snarls
    net_handle_t get_net_handle(const size_t& depth, const SnarlDistanceIndex* distance_index) ;
    ///Get the information that was stored to get the address in the distance index
    ///This is the connected component number for a root structure, or the address of
    ///an irregular snarl. Throws an error for anything else
    ///This is used for checking equality without looking at the distance index.
    ///Use get_net_handle for getting the actual handle
    size_t get_distance_index_address(const size_t& depth) ;


    ///Are the two decoders pointing to the same snarl tree node at the given depth
    ///This only checks if the values in the zipcode are the same at the given depth, 
    ///so if the preceeding snarl tree nodes are different, 
    ///then this might actually refer to different things
    static inline bool is_equal(zipcode_decoder_t& decoder1, zipcode_decoder_t& decoder2,
                                const size_t& depth);

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


    //How do decode the zipcode to get the old payload values
    static size_t record_offset(const zipcode_t& zip, const SnarlDistanceIndex& distance_index, const nid_t& id);

    static size_t parent_record_offset(const zipcode_t& zip, const SnarlDistanceIndex& distance_index, const nid_t& id);

    static size_t node_record_offset(const zipcode_t& zip, const SnarlDistanceIndex& distance_index, const nid_t& id);

    static size_t node_length(const zipcode_t& zip);

    static bool is_reversed(const zipcode_t& zip, const SnarlDistanceIndex& distance_index, const nid_t& id);

    static bool is_trivial_chain (const zipcode_t& zip);

    static bool parent_is_chain(const zipcode_t& zip, const SnarlDistanceIndex& distance_index, const nid_t& id);

    static bool parent_is_root (const zipcode_t& zip);

    static size_t prefix_sum (const zipcode_t& zip, const SnarlDistanceIndex& distance_index, const nid_t& id);

    static size_t chain_component (const zipcode_t& zip, const SnarlDistanceIndex& distance_index, const nid_t& id);

    
}; 
}

#endif
