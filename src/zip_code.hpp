#ifndef VG_ZIP_CODE_HPP_INCLUDED
#define VG_ZIP_CODE_HPP_INCLUDED

#include "varint.hpp"
#include "snarl_distance_index.hpp"

namespace vg{
using namespace std;

//A decoded zip code as a vector of pair<is_chain, index>
//where is_chain indicates whether it's a chain/node, and index
//is the index of the node/snarl/chain code in the varint_vector_t
typedef std::vector<pair<bool, size_t>> zip_code_decoder_t;


/* Zip codes store the snarl decomposition location and distance information for a position on a graph
 * A zip code will contain all the information necessary to compute the minimum distance between two 
 * positions, with minimal queries to the distance index
 */

struct zip_code_t {

    public:

        //Constructor for a position and a distance index
        void fill_in_zip_code (const SnarlDistanceIndex& distance_index, const vg::pos_t& pos);

        zip_code_decoder_t decode() const;

        //Get the exact minimum distance between two positions and their zip codes
        static inline size_t minimum_distance_between(const zip_code_t& zip1, const pos_t& pos1, 
                                       const zip_code_t& zip2, const pos_t& pos2,
                                       const SnarlDistanceIndex& distance_index);

        //Return true if the minimum distance between the zip codes is definitely greater than limit
        //A false result is inconclusive
        static inline bool is_farther_than(const zip_code_t& zip1, const zip_code_t& zip2, const size_t& limit);

    //TODO: Make this private:
        varint_vector_t zip_code;

    private:

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

}

#endif
