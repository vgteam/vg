#ifndef VG_ZIP_CODE_HPP_INCLUDED
#define VG_ZIP_CODE_HPP_INCLUDED

#include "varint.hpp"
#include "snarl_distance_index.hpp"
#include <gbwtgraph/minimizer.h>

namespace vg{
using namespace std;

/**
 * Zipcodes are structures that store distance index information for a node in a graph.
 * Their basic structure is a vector of "codes", with one code for each snarl tree node 
 * (node/snarl/chain) that is the ancestor of the node, starting with the root-level
 * structure and going down to the node.
 * Each code has an identifier and information used to calculate distances.
 *
 * A ZipCode stores the information and can be used to create a zipcode. It can be used
 * to calculate the distance between zipcodes
 *
 * A ZipCodeDecoder is used for interpreting zipcodes to find specific values that were
 * stored in the ZipCode. A ZipCodeDecoder must be constructed from a specific zipcode.
 * Construction of a decoder occurs one code at a time, starting from the root snarl or chain, 
 * so it is possible to have a partially constructed ZipCodeDecoder, to avoid having to 
 * walk through the entire ZipCode to get the values for things higher in the snarl tree.
 * The full decoder must be constructed to get values for the node.
 */

///A decoder for interpreting a zipcode
///Can interpret the values for a snarl tree node given the depth 
///(depth in the snarl tree, also the index into the zipcode vector)
class ZipCodeDecoder;


///The type of codes that can be stored in the zipcode
///TOP_LEVEL_IRREGULAR_SNARL is kind of a special case of an irregular snarl that is the child of a top-level chain
enum code_type_t { NODE = 1, CHAIN, REGULAR_SNARL, IRREGULAR_SNARL, ROOT_SNARL, ROOT_CHAIN, ROOT_NODE, TOP_LEVEL_IRREGULAR_SNARL};

///A struct to interpret the minimizer payload
///I want to use zipcodes as the payload but at the moment clustering still expects the old payload
///This can interpret zipcodes to format them as the old payload
struct MIPayload;

/* Zip codes store the snarl decomposition location and distance information for a position on a graph
 * A zip code will contain all the information necessary to compute the minimum distance between two 
 * positions, with minimal queries to the distance index
 */
class ZipCode {

    public:

        //Fill in an empty zipcode given a position
        void fill_in_zipcode (const SnarlDistanceIndex& distance_index, const vg::pos_t& pos);

        //Fill in an empty zipcode using the information that was stored in a payload
        void fill_in_zipcode_from_payload(const gbwtgraph::Payload& payload); 

        //Get the exact minimum distance between two positions and their zip codes
        //If distance_limit is set, return std::numeric_limits<size_t>::max() if the distance
        //will be greater than the distance limit
        //static size_t minimum_distance_between(const ZipCode& zip1, const pos_t& pos1, 
        //                               const ZipCode& zip2, const pos_t& pos2,
        //                               const SnarlDistanceIndex& distance_index, 
        //                               size_t distance_limit = std::numeric_limits<size_t>::max(),
        //                               bool directed_distance=true, 
        //                               const HandleGraph* graph = nullptr);

        //The same thing but using a zipcode decoder (which also has a pointer to the zipcode)
        //This is faster because otherwise the zipcode would need to be decoded
        //The decoders may or may not be filled in, and may be filled in when this is run
        //If distance_limit is set, return std::numeric_limits<size_t>::max() if the distance
        //will be greater than the distance limit
        static size_t minimum_distance_between(ZipCodeDecoder& zip_decoder1, const pos_t& pos1, 
                                       ZipCodeDecoder& zip_decoder2, const pos_t& pos2,
                                       const SnarlDistanceIndex& distance_index, 
                                       size_t distance_limit = std::numeric_limits<size_t>::max(),
                                       bool undirected_distance=false, 
                                       const HandleGraph* graph = nullptr);

        //Return true if the minimum distance between the zip codes is definitely greater than limit
        //A false result is inconclusive
        static bool is_farther_than(const ZipCode& zip1, const ZipCode& zip2, const size_t& limit);

        //Get a tuple of the top-level structure id, prefix sum of the child of the top-level chain, and 
        //the length of the child of the top-level chain
        //This gets used to quickly compare the two zip codes for is_farther_than
        static tuple<size_t, size_t, size_t> get_top_level_chain_offset();


        //////////////////Functions to work with minimizer payloads for clustering 
        // Since we're sill using the old implementation, we need to be able to 
        // switch from zipcodes to payloads and back

        //Encode zip code so it can be stored in the payload
        gbwtgraph::Payload get_payload_from_zip() const;
        typedef std::uint64_t code_type; // We assume that this fits into gbwtgraph::Payload.


        ///How many bytes were used to store this zipcode?
        size_t byte_count() const {
            return zipcode.byte_count();
        }

    //TODO: Make this private:
        //The actual data for a zipcode is a vector of ints
        varint_vector_t zipcode;


        /// Equality operator
        inline bool operator== (const ZipCode& other) const {
            return zipcode == other.zipcode;
        }

        /// Dump to a normal vector
        std::vector<size_t> to_vector() const;

        /// Load from a normal vector
        void from_vector(const std::vector<size_t>& values);

    private:

        /* These offsets are used to define each type of "code"
        */
        //TODO: I still access these in order so the order can't change

        ///Offsets of values in a root chain or snarl code
        ///Roots have a bool for is_chain and an identifier, which is the
        ///connected component number from the distance index
        const static size_t ROOT_CHAIN_OR_SNARL_SIZE = 2;
        const static size_t ROOT_IS_CHAIN_OFFSET = 0;
        const static size_t ROOT_IDENTIFIER_OFFSET = 1;

        //If the zipcode is for a root-level node, then there are only three things
        //in the zipcode, and the last is the length of the node
        const static size_t ROOT_NODE_SIZE = 3;
        const static size_t ROOT_NODE_LENGTH_OFFSET = 2;

        ///Offsets for chain codes
        const static size_t CHAIN_SIZE = 2;
        const static size_t CHAIN_RANK_IN_SNARL_OFFSET = 0;
        const static size_t CHAIN_LENGTH_OFFSET = 1;

        ///Offsets for snarl codes
        const static size_t REGULAR_SNARL_SIZE = 4;
        const static size_t IRREGULAR_SNARL_SIZE = 2;
        const static size_t TOP_LEVEL_IRREGULAR_SNARL_SIZE = 4;
        const static size_t SNARL_IS_REGULAR_OFFSET = 0;

        const static size_t REGULAR_SNARL_OFFSET_IN_CHAIN_OFFSET = 1;
        const static size_t REGULAR_SNARL_LENGTH_OFFSET = 2;
        const static size_t REGULAR_SNARL_IS_REVERSED_OFFSET = 3;

        const static size_t IRREGULAR_SNARL_RECORD_OFFSET = 1;

        //These are only for top-level irregular snarls
        const static size_t IRREGULAR_SNARL_OFFSET_IN_CHAIN_OFFSET = 2;
        const static size_t IRREGULAR_SNARL_LENGTH_OFFSET = 3;

        ///Offsets for nodes
        const static size_t NODE_SIZE = 3;
        const static size_t NODE_OFFSET_OR_RANK_OFFSET = 0;
        const static size_t NODE_LENGTH_OFFSET = 1;
        const static size_t NODE_IS_REVERSED_OFFSET = 2;


        /* Functions for getting the code for each snarl/chain/node
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
        inline vector<size_t> get_top_level_irregular_snarl_code(const net_handle_t& snarl, const SnarlDistanceIndex& distance_index);
        //Return a vector of size_ts that will represent the snarl in the zip code
        inline vector<size_t> get_irregular_snarl_code(const net_handle_t& snarl, const SnarlDistanceIndex& distance_index);
    friend class ZipCodeDecoder;
};

//A struct for holding a vector of zipcodes
//This is really just used for serializing
struct zipcode_vector_t {
    vector<ZipCode>* zipcodes;
    zipcode_vector_t (vector<ZipCode>* z) {
        zipcodes = z;
    }

    void serialize(std::ostream& out) const;
    void deserialize(std::istream& in);
};


/*
 * Struct for interpreting a ZipCode
 */
class ZipCodeDecoder {

    public:
    //TODO: Make the decoder and zipcode private, still need it for unit testing
    ///The decoder as a vector of pair<is_chain, index>, one for each snarl tree node in the zip
    ///where is_chain indicates whether it's a chain/node, and index
    ///is the index of the node/snarl/chain code in the varint_vector_t
    std::vector<pair<bool, size_t>> decoder;

    ///The zipcode that this is decoding
    const ZipCode* zipcode;

    public:

    ///Constructor that goes through the zipcode and decodes it to fill in decoder
    ///If a depth is given, then only fill in up to depth snarl tree nodes
    ///Otherwise, fill in the whole zipcode
    ZipCodeDecoder(const ZipCode* zipcode, const size_t& depth=std::numeric_limits<size_t>::max());

    ///Go through the entire zipcode and fill in the decoder
    void fill_in_full_decoder();

    ///Fill in one more item in the decoder
    ///Returns true if this is the last thing in the zipcode and false if there is more to decode
    bool fill_in_next_decoder();

    ///How many codes in the zipcode have been decoded?
    size_t decoder_length() {return decoder.size();}

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
    static inline bool is_equal(ZipCodeDecoder& decoder1, ZipCodeDecoder& decoder2,
                                const size_t& depth);

    /// Dump a ZipCodeDecoder to a stream so that it can be reconstructed for a
    /// unit test from the resulting information.
    void dump(std::ostream& out) const;

};

std::ostream& operator<<(std::ostream& out, const ZipCodeDecoder& decoder); 


/**
    The payload for the minimizer index. This stores distance information that gets used in clustering
    The payload now uses zip codes, so this gets used to go from a zip code to distance information
    usable by the clusterer, which expects the old payload format
*/
struct MIPayload {
    typedef std::uint64_t code_type; // We assume that this fits into gbwtgraph::Payload.
    //typedef std::pair<code_type, code_type> payload_type;

    
    constexpr static gbwtgraph::Payload NO_CODE = {0, 0};
    constexpr static std::size_t NO_VALUE = std::numeric_limits<size_t>::max(); 


    //How do decode the zipcode to get the old payload values
    static size_t record_offset(const ZipCode& zip, const SnarlDistanceIndex& distance_index, const nid_t& id);

    static size_t parent_record_offset(const ZipCode& zip, const SnarlDistanceIndex& distance_index, const nid_t& id);

    static size_t node_record_offset(const ZipCode& zip, const SnarlDistanceIndex& distance_index, const nid_t& id);

    static size_t node_length(const ZipCode& zip);

    static bool is_reversed(const ZipCode& zip, const SnarlDistanceIndex& distance_index, const nid_t& id);

    static bool is_trivial_chain (const ZipCode& zip);

    static bool parent_is_chain(const ZipCode& zip, const SnarlDistanceIndex& distance_index, const nid_t& id);

    static bool parent_is_root (const ZipCode& zip);

    static size_t prefix_sum (const ZipCode& zip, const SnarlDistanceIndex& distance_index, const nid_t& id);

    static size_t chain_component (const ZipCode& zip, const SnarlDistanceIndex& distance_index, const nid_t& id);

    
}; 
}

#endif