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
 * A decoder is used for interpreting zipcodes to find specific values that were
 * stored in the ZipCode. 
 * Construction of a decoder occurs one code at a time, starting from the root snarl or chain, 
 * so it is possible to have a partially constructed decoder, to avoid having to 
 * walk through the entire ZipCode to get the values for things higher in the snarl tree.
 * The full decoder must be constructed to get values for the node.
 */



///A struct to interpret the minimizer payload
///I want to use zipcodes as the payload but at the moment clustering still expects the old payload
///This can interpret zipcodes to format them as the old payload
struct MIPayload;


/// A struct to be used as a unique identifier for a snarl tree node (node/snarl/chain)
/// using information from the zipcodes.
/// It should be unique and hashable
typedef std::string net_identifier_t;


/* Zip codes store the snarl decomposition location and distance information for a position on a graph
 * A zip code will contain all the information necessary to compute the minimum distance between two 
 * positions, with minimal queries to the distance index
 */
class ZipCode {

    ///structs to store an unpacked version of one node/snarl/chain code
    public:
        struct node_code_t;
        struct chain_code_t;
        struct snarl_code_t;


    /// The type of codes that can be stored in the zipcode
    /// Trivial chains that are children of snarls get saved as a chain with no child node
    /// EMPTY doesn't actually mean anything, it's used to catch errors
    /// Snarls can be regular, irregular, or cyclic. 
    /// Regular snarls are bubbles. Irregular snarls are snarls that aren't bubbles but are dags
    /// Cyclic snarls are non-dags. They are stored the same as irregular snarls. Only the type is different
    public:
        enum code_type_t { NODE = 1, CHAIN, REGULAR_SNARL, IRREGULAR_SNARL, CYCLIC_SNARL, ROOT_SNARL, ROOT_CHAIN, ROOT_NODE, EMPTY };

    public:

        //Fill in an empty zipcode given a position
        void fill_in_zipcode (const SnarlDistanceIndex& distance_index, const vg::pos_t& pos, bool fill_in_decoder = true);

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
        static size_t minimum_distance_between(ZipCode& zip1, const pos_t& pos1, 
                                       ZipCode& zip2, const pos_t& pos2,
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
        const static size_t ROOT_SNARL_SIZE = 2;
        const static size_t ROOT_IS_CHAIN_OFFSET = 0;
        const static size_t ROOT_IDENTIFIER_OFFSET = 1;

        //FOr a chain, also include the component count
        const static size_t ROOT_CHAIN_SIZE = 4;
        const static size_t ROOT_CHAIN_COMPONENT_COUNT_OFFSET = 2;
        //This is a bitvector storing if there is connectivity between the bounds of the node/chain
        const static size_t ROOT_NODE_OR_CHAIN_CONNECTIVITY_OFFSET = 3;

        //If the zipcode is for a root-level node, then there are only three things
        //in the zipcode, and the last is the length of the node
        const static size_t ROOT_NODE_SIZE = 4;
        const static size_t ROOT_NODE_LENGTH_OFFSET = 2;

        ///Offsets for chain codes
        const static size_t CHAIN_SIZE = 3;
        const static size_t CHAIN_RANK_IN_SNARL_OFFSET = 0;
        //This is the distance index's chain_minimum_length, meaning that if it's a multicomponent chain,
        //then it is the length of the last component.
        const static size_t CHAIN_LENGTH_OFFSET = 1;

        //This tells us if the chain is a multicomponent chain, how many components it has, and if the chain loops
        //The value is the component of the last node in the chain * 2, +1 if the chain loops
        //So 0 means normal chain, 1 means one component but the chain loops, 2 means 2 components, 3 means 2 components with a loop...
        //This is maybe not the most efficient way of storing it but since it is pretty rare for the chains to 
        //be multicomponent chains and rarer for them to loop, and the multicomponent chains probably won't have
        //a lot of components anyway, this is more efficient for the majority of cases when the value will be 0
        const static size_t CHAIN_COMPONENT_COUNT_OFFSET = 2;

        ///Offsets for snarl codes
        const static size_t REGULAR_SNARL_SIZE = 6;
        const static size_t IRREGULAR_SNARL_SIZE = 10;

        //Both regular and irregular snarls have these

        // This will be 0 for irregular snarl, 1 for regular, and 2 for non-dag irregular snarls
        // cyclic snarls will be identical to irregular snarls except for SNARL_IS_REGULAR
        const static size_t SNARL_IS_REGULAR_OFFSET = 0; 
        const static size_t SNARL_OFFSET_IN_CHAIN_OFFSET = 1;
        const static size_t SNARL_LENGTH_OFFSET = 2;
        const static size_t SNARL_CHILD_COUNT_OFFSET = 3;
        //THis will be the lower of the two component numbers if the snarl spans two
        //This only happens if the snarl is not start-end connected, which we'll know from the length
        const static size_t SNARL_CHAIN_COMPONENT_OFFSET = 4;

        //Only for regular snarls
        const static size_t REGULAR_SNARL_IS_REVERSED_OFFSET = 5;

        //Only for irregular snarls
        const static size_t IRREGULAR_SNARL_RECORD_OFFSET = 5;
        //Distance from the left side of the child to the start of the snarl
        const static size_t IRREGULAR_SNARL_DISTANCE_LEFT_START_OFFSET = 6;
        const static size_t IRREGULAR_SNARL_DISTANCE_LEFT_END_OFFSET = 7;
        const static size_t IRREGULAR_SNARL_DISTANCE_RIGHT_START_OFFSET = 8;
        const static size_t IRREGULAR_SNARL_DISTANCE_RIGHT_END_OFFSET = 9;

        ///Offsets for nodes
        const static size_t NODE_SIZE = 4;
        const static size_t NODE_OFFSET_OFFSET = 0;
        const static size_t NODE_LENGTH_OFFSET = 1;
        const static size_t NODE_IS_REVERSED_OFFSET = 2;
        const static size_t NODE_CHAIN_COMPONENT_OFFSET = 3;


        /* Functions for getting the code for each snarl/chain/node
         * Distances will be stored as distance+1, 0 will be reserved for inf
         */
        //Return a node_code_t that will represent the node in the zip code
        inline node_code_t get_node_code(const net_handle_t& node, const SnarlDistanceIndex& distance_index);
        //Return a chain_code_t that will represent the chain in the zip code
        //The actual values being stored, not the raw values
        inline chain_code_t get_chain_code(const net_handle_t& chain, const SnarlDistanceIndex& distance_index);
        //Return a vector of size_ts that will represent the snarl in the zip code
        inline snarl_code_t get_regular_snarl_code(const net_handle_t& snarl, const net_handle_t& snarl_child, 
                                                            const SnarlDistanceIndex& distance_index);
        //Return a vector of size_ts that will represent the snarl in the zip code
        inline snarl_code_t get_irregular_snarl_code(const net_handle_t& snarl, const net_handle_t& snarl_child, const SnarlDistanceIndex& distance_index);

    public:

        /* Functions to get the values out of the zipcode for one code
           The decoded code might not have all the values set*/

        // Get a node_code_t for the given level
        //For a root node, use a chain
        node_code_t unpack_node_code(size_t zipcode_level);
        //Return a chain_code_t that will represent the chain in the zip code
        //The actual values being stored, not the raw values
        chain_code_t unpack_chain_code(size_t zipcode_level);
        //Return a vector of size_ts that will represent the snarl in the zip code
        snarl_code_t unpack_snarl_code(size_t zipcode_level);


    //////////////////////////////// Stuff for decoding the zipcode

    public:
        //TODO: Make the decoder and zipcode private, still need it for unit testing
        ///The decoder as a vector of pair<is_chain, index>, one for each snarl tree node in the zip
        ///where is_chain indicates whether it's a chain/node, and index
        ///is the index of the node/snarl/chain code in the varint_vector_t
        struct decoder_t {
            bool is_chain : 1;
            size_t offset : 15;
            decoder_t(bool is_chain, size_t offset) : is_chain(is_chain), offset(offset) {}
            inline bool operator==(const decoder_t& other) const {
                return is_chain == other.is_chain && offset == other.offset;
            }
        };
        std::vector<decoder_t> decoder;

        ///Did we fill in the entire decoder
        ///TODO: I'm making it fill in the decoder automatically because it seems to be faster that way, instead of 
        /// waiting to see which parts are actually needed
        bool finished_decoding = false;

        public:

        ///Go through the entire zipcode and fill in the decoder
        void fill_in_full_decoder();

        ///Fill in one more item in the decoder
        ///Returns true if this is the last thing in the zipcode and false if there is more to decode
        bool fill_in_next_decoder();

        ///What is the maximum depth of this zipcode?
        size_t max_depth() const;

        ///How many codes in the zipcode have been decoded?
        size_t decoder_length() const {return decoder.size();}

        ///What type of snarl tree node is at the given depth (index into the zipcode)
        ZipCode::code_type_t get_code_type(const size_t& depth) const ;

        ///Get the length of a snarl tree node given the depth in the snarl tree
        ///This requires the distance index for irregular snarls (except for a top-level snarl)
        ///Throws an exception if the distance index is not given when it is needed
        ///Doesn't use a given distance index if it isn't needed
        size_t get_length(const size_t& depth, const SnarlDistanceIndex* distance_index=nullptr, bool get_chain_component_length = false) const ;

        ///Get the rank of a node/snarl in a snarl. Throw an exception if it isn't the child of a snarl
        size_t get_rank_in_snarl(const size_t& depth) const ;

        ///Get the number of children in a snarl. Throw an exception if it isn't a snarl
        size_t get_snarl_child_count(const size_t& depth, const SnarlDistanceIndex* distance_index=nullptr) const ;

        ///Get the prefix sum of a child of a chain
        ///This requires the distance index for irregular snarls (except for a top-level snarl)
        ///Throws an exception if the distance index is not given when it is needed
        ///Doesn't use a given distance index if it isn't needed
        size_t get_offset_in_chain(const size_t& depth, const SnarlDistanceIndex* distance_index=nullptr) const ;

        ///Get the chain component of a chain child.
        ///For snarls, this will be the component of the start node
        size_t get_chain_component(const size_t& depth) const ;

        ///Get the chain component of the last node in the chain
        /// This behaves like the distance index get_chain_component- 
        /// for looping chains it returns the last component if get_end is true,
        /// and 0 if it is false
        size_t get_last_chain_component(const size_t& depth, bool get_end = false) const ;
        bool get_is_looping_chain(const size_t& depth) const ;

        ///Is the snarl tree node backwards relative to its parent
        bool get_is_reversed_in_parent(const size_t& depth) const;

        ///Get the handle of the thing at the given depth. This can only be used for
        ///Root-level structures or irregular snarls
        net_handle_t get_net_handle(const size_t& depth, const SnarlDistanceIndex* distance_index)  const;

        ///Get the handle of the thing at the given depth. This can be used for anything but is slow,
        /// even for roots and irregular/cyclic snarls. It's a separate function to make sure I
        /// remember that it's slow
        net_handle_t get_net_handle_slow(nid_t id, const size_t& depth, const SnarlDistanceIndex* distance_index) const;

        ///Get the information that was stored to get the address in the distance index
        ///This is the connected component number for a root structure, or the address of
        ///an irregular snarl. Throws an error for anything else
        ///This is used for checking equality without looking at the distance index.
        ///Use get_net_handle for getting the actual handle
        size_t get_distance_index_address(const size_t& depth)  const;

        /// The minimum distance from start or end of the snarl to the left or right side of the child
        size_t get_distance_to_snarl_bound(const size_t& depth, bool snarl_start, bool left_side) const;

        bool is_externally_start_end_connected(const size_t& depth) const;
        bool is_externally_start_start_connected(const size_t& depth) const;
        bool is_externally_end_end_connected(const size_t& depth) const;


        ///Are the two decoders pointing to the same snarl tree node at the given depth
        ///This only checks if the values in the zipcode are the same at the given depth, 
        ///so if the preceeding snarl tree nodes are different, 
        ///then this might actually refer to different things
        const static bool is_equal(const ZipCode& zip1, const ZipCode& zip2,
                                    const size_t& depth);

        /// Dump a ZipCode to a stream so that it can be reconstructed for a
        /// unit test from the resulting information.
        void dump(std::ostream& out) const;

        //TODO: I want to make a struct for holding all values of a code as real values

        ///Fill in a payload with values from the zipcode
        MIPayload get_payload_from_zipcode(nid_t id, const SnarlDistanceIndex& distance_index) const;

        /// Get an identifier for the snarl tree node at this depth. If the snarl tree node at this depth
        /// would be the node, also include the node id
        net_identifier_t get_identifier(size_t depth) const;
        const static net_identifier_t get_parent_identifier(const net_identifier_t& child);

};

/// Print a code type to a stream
std::ostream& operator<<(std::ostream& out, const ZipCode::code_type_t& type);


//A structure for holding a vector of zipcodes
//This is really just used for serializing
class ZipCodeCollection {
    private:
    vector<ZipCode> zipcodes;

    public:
    ZipCodeCollection () {}

    void serialize(std::ostream& out) const;
    void deserialize(std::istream& in);
    bool empty() const {return zipcodes.empty();}
    ZipCode at(size_t i) const {return zipcodes.at(i);}
    void emplace_back(ZipCode zip) {zipcodes.emplace_back(zip);}
    size_t size() const { return zipcodes.size();}

    private:

    //magic number to identify the file
    const static uint32_t magic_number = 0x5a495053; //ZIPS
    const static uint32_t version = 3;

    public:
    const static std::uint32_t get_magic_number() {return magic_number;}
    const static std::string get_magic_number_as_string() {
        std::uint32_t num = get_magic_number();
        return std::string(reinterpret_cast<char*>(&num), sizeof(num));
    }


};



/**
    An unpacked version of one node code
    The values actually stored are the same ones that get stored in the zipcode
    This has getters and setters for getting the actual value, 
    and getters and setters for getting the raw values 
*/
struct ZipCode::node_code_t {
    private:
    //Prefix sum for a nested node, address for a root node
    size_t prefix_sum_or_identifier ;
    size_t chain_component : 32;
    size_t length : 31;
    bool is_reversed;

    public:

    ////// Raw getters
    size_t get_raw_prefix_sum_or_identifier() {return prefix_sum_or_identifier;}
    size_t get_raw_chain_component() {return chain_component;}
    size_t get_raw_length() {return length;}
    bool get_raw_is_reversed() {return is_reversed;}

    ///// Raw setters
    void set_raw_prefix_sum_or_identifier(size_t val) {prefix_sum_or_identifier = val;}
    void set_raw_chain_component(size_t val) {chain_component = val;}
    void set_raw_length(size_t val) {length = val;}
    void set_raw_is_reversed(bool val) {is_reversed = val;}

    //// Real value setters
    size_t get_prefix_sum_or_identifier() {return prefix_sum_or_identifier == 0 ? numeric_limits<size_t>::max() : prefix_sum_or_identifier-1;}
    size_t get_chain_component() {return chain_component;}
    size_t get_length() {return length-1;}
    bool get_is_reversed() {return is_reversed;}

    ////Real value getters
    void set_prefix_sum_or_identifier(size_t val) {prefix_sum_or_identifier = val == std::numeric_limits<size_t>::max() ? 0 : val+1;}
    void set_chain_component(size_t val) {chain_component = val == std::numeric_limits<size_t>::max() ? 0 : val;}
    void set_length(size_t val) {length = val+1;}
    void set_is_reversed(bool val) {is_reversed = val;}
}; 

/**
    An unpacked version of one chain code
    The values actually stored are the same ones that get stored in the zipcode
    This has getters and setters for getting the actual value, 
    and getters and setters for getting the raw values 
*/
struct ZipCode::chain_code_t {


    private:
    //The length of the last component of the chain (which may be the whole chain)
    size_t length;
    //The rank in the parent snarl or, if it is a root chain, the identifier
    size_t snarl_rank_or_identifier : 32;

    //This stores the component and is_looping_chain
    size_t last_component : 16;

    //For root chain/nodes, a bitvector representing the connectivity
    size_t connectivity : 4;


    public:
    size_t get_raw_length() {return length;}
    size_t get_raw_snarl_rank_or_identifier() {return snarl_rank_or_identifier;}
    size_t get_raw_last_component() {return last_component;}
    size_t get_raw_connectivity() {return connectivity;}
    void set_raw_length(size_t val) {length = val;}
    void set_raw_snarl_rank_or_identifier(size_t val) {snarl_rank_or_identifier = val;}
    void set_raw_last_component(size_t val) {last_component = val;}
    void set_raw_connectivity (size_t val){connectivity = val;}

    size_t get_length() {
        return length == 0 ? std::numeric_limits<size_t>::max() : length-1;
    }
    size_t get_snarl_rank_or_identifier() {return snarl_rank_or_identifier;}
    size_t get_last_component() {
        if (last_component % 2 ) {
            return (last_component-1) / 2;
        } else {
            return last_component / 2;
        }
    }

    size_t get_connectivity() {return connectivity;}
    bool get_is_looping_chain() {return last_component % 2;}

    void set_length(size_t val) {
        length = val == std::numeric_limits<size_t>::max() ? 0 : val+1;
    }
    void set_snarl_rank_or_identifier(size_t val) {
        snarl_rank_or_identifier = val;
    }
    void set_last_component(size_t comp, bool loops) {
        comp = comp == std::numeric_limits<size_t>::max() ? 0 : comp*2;
        if (loops) { comp ++;}
        last_component = comp;
    }
    void set_connectivity(size_t val) {connectivity = val;}
}; 

/**
    An unpacked version of one snarl code
    The values actually stored are the same ones that get stored in the zipcode
    This has getters and setters for getting the actual value, 
    and getters and setters for getting the raw values 
*/
struct ZipCode::snarl_code_t {

    private:
        size_t length;
        //Prefix sum for a nested snarl, identifier for a root snarl
        size_t prefix_sum_or_identifier;

        size_t distance_start_left;
        size_t distance_start_right;
        size_t distance_end_left;
        size_t distance_end_right;

        size_t record_offset : 32;

        size_t child_count : 16; 
        size_t chain_component : 16;

        size_t code_type : 4;

        bool is_reversed;

    public:
        //We use getters and setters to deal with things that are max() but stored as 0
        //and getters and setters for the raw values. These are sometimes redundant

        size_t get_raw_length() {return length;}
        size_t get_raw_prefix_sum_or_identifier () {return prefix_sum_or_identifier;}
        size_t get_raw_distance_start_left () {return distance_start_left;}
        size_t get_raw_distance_start_right () {return distance_start_right;}
        size_t get_raw_distance_end_left () {return distance_end_left;}
        size_t get_raw_distance_end_right () {return distance_end_right;}
        size_t get_raw_record_offset () { return record_offset;}
        size_t get_raw_child_count() {return child_count;}
        size_t get_raw_chain_component() {return chain_component;}
        size_t get_raw_code_type() {return code_type;}
        bool get_raw_is_reversed() {return is_reversed;}

        void set_raw_length(size_t val) {length = val;}
        void set_raw_prefix_sum_or_identifier (size_t val) {prefix_sum_or_identifier = val;}
        void set_raw_distance_start_left (size_t val) {distance_start_left = val;}
        void set_raw_distance_start_right (size_t val) {distance_start_right = val;}
        void set_raw_distance_end_left (size_t val) {distance_end_left = val;}
        void set_raw_distance_end_right (size_t val) {distance_end_right = val;}
        void set_raw_record_offset (size_t val) { record_offset = val;}
        void set_raw_child_count(size_t val) {child_count = val;}
        void set_raw_chain_component(size_t val) {chain_component = val;}
        void set_raw_code_type(size_t val) {code_type = val;}
        void set_raw_is_reversed(bool val) {is_reversed = val;}



        //// Getters
        size_t get_length() {
            return length == 0 ? std::numeric_limits<size_t>::max() : length-1;
        }
        size_t get_prefix_sum_or_identifier() {
            return prefix_sum_or_identifier == 0 ? std::numeric_limits<size_t>::max() : prefix_sum_or_identifier-1;
        }

        //distance from the left side of the child to the start of the snarl
        //or, for root nodes/chains, start-start connected
        //start-right and end-left are the same for root nodes/chains
        size_t get_distance_start_left() {
            return distance_start_left == 0 ? std::numeric_limits<size_t>::max() : distance_start_left-1;
        }
        size_t get_distance_start_right() {
            return distance_start_right == 0 ? std::numeric_limits<size_t>::max() : distance_start_right-1;
        }
        size_t get_distance_end_left() {
            return distance_end_left == 0 ? std::numeric_limits<size_t>::max() : distance_end_left-1;
        }
        size_t get_distance_end_right() {
            return distance_end_right == 0 ? std::numeric_limits<size_t>::max() : distance_end_right-1;
        }

        size_t get_record_offset() {return record_offset;}

        size_t get_child_count() {return child_count;} 
        size_t get_chain_component() {return chain_component;}

        size_t get_code_type() {return code_type;}

        bool get_is_reversed() {return is_reversed;}

        //////// Setters
        void set_length(size_t val) {
            length = val == std::numeric_limits<size_t>::max() ? 0 : val+1;
        }
        void set_prefix_sum_or_identifier(size_t val) {
            prefix_sum_or_identifier = val == std::numeric_limits<size_t>::max() ? 0 : val+1;
        }

        void set_distance_start_left(size_t val) {
            distance_start_left = val == std::numeric_limits<size_t>::max() ? 0 : val+1;
        }
        void set_distance_start_right(size_t val) {
            distance_start_right = val == std::numeric_limits<size_t>::max() ? 0 : val+1;
        }
        void set_distance_end_left(size_t val) {
            distance_end_left = val == std::numeric_limits<size_t>::max() ? 0 : val+1;
        }
        void set_distance_end_right(size_t val) {
            distance_end_right = val == std::numeric_limits<size_t>::max() ? 0 : val+1;
        }

        void set_record_offset(size_t val) {
            record_offset = val;
        }

        void set_child_count(size_t val) {
            child_count = val;
        }

        void set_chain_component(size_t val) {
            chain_component = val == std::numeric_limits<size_t>::max() ? 0 : val;
        }

        void set_code_type(size_t val) {
            code_type = val;
        }

        void set_is_reversed(bool val) {
            is_reversed = val;
        }

};


template<>
struct wang_hash<net_identifier_t> {
    size_t operator()(const net_identifier_t& id) const {
        return wang_hash<std::string>()(id);
    }
};

std::ostream& operator<<(std::ostream& out, const ZipCode& decoder); 


/**
    The payload for the minimizer index. This stores distance information that gets used in clustering
    The payload now uses zip codes, so this gets used to go from a zip code to distance information
    usable by the clusterer
*/
struct MIPayload {    
    typedef std::uint64_t code_type; // We assume that this fits into gbwtgraph::Payload.
    //typedef std::pair<code_type, code_type> payload_type;


    constexpr static gbwtgraph::Payload NO_CODE = {0, 0};
    constexpr static std::size_t NO_VALUE = std::numeric_limits<size_t>::max();


    net_handle_t node_handle;
    net_handle_t parent_handle;

    size_t node_length = std::numeric_limits<size_t>::max();
    size_t prefix_sum = 0;
    size_t chain_component = 0;
    //Depth according to the distance index
    size_t parent_depth = 0;
    size_t parent_record_offset = 0;

    ZipCode::code_type_t parent_type = ZipCode::EMPTY;
    bool is_reversed = false;
    bool is_trivial_chain = false;
    bool parent_is_chain = false;
    bool parent_is_root = false;
}; 
}

#endif
