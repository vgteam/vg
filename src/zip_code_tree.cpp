#define DEBUG_ZIP_CODE_TREE

#include "zip_code_tree.hpp"

using namespace std;
namespace vg {

ZipCodeTree::ZipCodeTree(vector<SnarlDistanceIndexClusterer::Seed>& seeds) :
    seeds(seeds) {

    /*
    Constructor for the ZipCodeTree
    Takes a vector of seeds and constructs the tree

    Tree construction is done by first sorting the seeds along chains/snarls
    Then, adding each seed, snarl/chain boundary, and distance to zip_code_tree
    Finally (optionally), the tree is refined to take out unnecessary edges
    */

    //////////////////// Sort the seeds

    //A vector of indexes into seeds
    //To be sorted along each chain/snarl the snarl tree
    vector<size_t> seed_indices (seeds.size(), 0);
    for (size_t i = 0 ; i < seed_indices.size() ; i++) {
        seed_indices[i] = i;
    }

    //Sort the indices
    std::sort(seed_indices.begin(), seed_indices.end(), [&] (const size_t& a, const size_t& b) {
        //Comparator returning a < b
#ifdef DEBUG_ZIP_CODE_TREE
        cerr << "Comparing seeds " << seeds[a].pos << " and " << seeds[b].pos << endl;
#endif
        size_t depth = 0;
        while (depth < seeds[a].zipcode_decoder->max_depth() &&
               depth < seeds[b].zipcode_decoder->max_depth() &&
               ZipCodeDecoder::is_equal(*seeds[a].zipcode_decoder, *seeds[b].zipcode_decoder, depth)) {
            depth++;
        }
#ifdef DEBUG_ZIP_CODE_TREE
        cerr << "\t different at depth " << depth << endl;
#endif
        //Either depth is the last thing in a or b, or they are different at this depth


        if ( ZipCodeDecoder::is_equal(*seeds[a].zipcode_decoder, *seeds[b].zipcode_decoder, depth)) {
#ifdef DEBUG_ZIPCODE_CLUSTERING
            cerr << "\tthey are on the same node" << endl;
#endif
            //If they are equal, then they must be on the same node

            size_t offset1 = is_rev(seeds[a].pos)
                           ? seeds[a].zipcode_decoder->get_length(depth) - offset(seeds[a].pos) - 1
                           : offset(seeds[a].pos);
            size_t offset2 = is_rev(seeds[b].pos)
                           ? seeds[b].zipcode_decoder->get_length(depth) - offset(seeds[b].pos) - 1
                           : offset(seeds[b].pos);
            if (!seeds[a].zipcode_decoder->get_is_reversed_in_parent(depth)) {
                //If they are in a snarl or they are facing forward on a chain, then order by
                //the offset in the node
                return offset1 < offset2;
            } else {
                //Otherwise, the node is facing backwards in the chain, so order backwards in node
                return offset2 < offset1;
            }
        }  else if (depth == 0) {
#ifdef DEBUG_ZIPCODE_CLUSTERING
            cerr << "\tThey are on different connected components" << endl;
#endif
            //If they are on different connected components, sort by connected component
            return seeds[a].zipcode_decoder->get_distance_index_address(0) < seeds[b].zipcode_decoder->get_distance_index_address(0);
            
        }  else if (seeds[a].zipcode_decoder->get_code_type(depth-1) == CHAIN || seeds[a].zipcode_decoder->get_code_type(depth-1) == ROOT_CHAIN) {
#ifdef DEBUG_ZIPCODE_CLUSTERING
            cerr << "\t they are children of a common chain" << endl;
#endif
            //If a and b are both children of a chain
            size_t offset_a = seeds[a].zipcode_decoder->get_offset_in_chain(depth);
            size_t offset_b = seeds[b].zipcode_decoder->get_offset_in_chain(depth);
            if ( offset_a == offset_b) {
                //If they have the same prefix sum, then the snarl comes first
                //They will never be on the same child at this depth
                return seeds[a].zipcode_decoder->get_code_type(depth) != NODE && seeds[b].zipcode_decoder->get_code_type(depth) == NODE;  
            } else {
                return offset_a < offset_b;
            }
        } else if (seeds[a].zipcode_decoder->get_code_type(depth-1) == REGULAR_SNARL) {
#ifdef DEBUG_ZIPCODE_CLUSTERING
            cerr << "\t they are children of a common regular snarl" << endl;
#endif
            //If the parent is a regular snarl, then sort by order along the parent chain
            size_t offset1 = is_rev(seeds[a].pos) 
                           ? seeds[a].zipcode_decoder->get_length(depth) - offset(seeds[a].pos) - 1
                           : offset(seeds[a].pos); 
            size_t offset2 = is_rev(seeds[b].pos) 
                           ? seeds[b].zipcode_decoder->get_length(depth) - offset(seeds[b].pos) - 1
                           : offset(seeds[b].pos);
            if (seeds[a].zipcode_decoder->get_is_reversed_in_parent(depth)) {
                return offset1 < offset2;
            } else {
                return offset2 < offset1;
            }
        } else {
#ifdef DEBUG_ZIPCODE_CLUSTERING
            cerr << "\t they are children of a common irregular snarl" << endl;
#endif
            //Otherwise, they are children of an irregular snarl
            //Sort by the distance to the start of the irregular snarl
            size_t distance_to_start_a = seeds[a].zipcode_decoder->get_distance_to_snarl_start(depth);
            size_t distance_to_start_b = seeds[b].zipcode_decoder->get_distance_to_snarl_start(depth);
            if (distance_to_start_a == distance_to_start_b) {
                //If they are equi-distant to the start of the snarl, then put the one that is
                //farther from the end first

                return seeds[a].zipcode_decoder->get_distance_to_snarl_end(depth) >
                         seeds[b].zipcode_decoder->get_distance_to_snarl_end(depth);
            } else {
                return distance_to_start_a < distance_to_start_b;
            }
        } 
    });


}

}
