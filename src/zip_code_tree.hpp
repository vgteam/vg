#ifndef VG_ZIP_CODE_TREE_HPP_INCLUDED

#define VG_ZIP_CODE_TREE_HPP_INCLUDED

#include "zip_code.hpp"
#include "snarl_seed_clusterer.hpp"

namespace vg{
using namespace std;

/**

A ZipCodeTree takes a set of SnarlDistanceIndexCluserer::Seed's (seed alignments between a read and reference) 
and provides an iterator that, given a seed and a distance limit, iterates through seeds that are
reachable within the distance limit

Generally, this will take a collection of seeds and build a tree structure representing the connectivity
of the seeds, based on the snarl decomposition
Edges are labelled with distance values.
The tree can be traversed to find distances between seeds
*/
class ZipCodeTree {

    public:

    /**
     * Constructor
     * The constructor creates a tree of the input seeds that is used for calculating distances
     */
    ZipCodeTree(vector<SnarlDistanceIndexClusterer::Seed>& seeds);

    /** Given a seed, run iteratee on all seeds that are reachable with a minimum distance
     * less than or equal to the distance_limit
     *
     * iteratee returns false to stop ant true to continue
     * returns false if the iteration stopped early, and true if it completed
     */
    bool for_each_seed_within_distance_range(SnarlDistanceIndexClusterer::Seed& seed, 
            const size_t distance_limit,
            const std::function<bool(SnarlDistanceIndexClusterer::Seed&, size_t)>& iteratee) const;

    private:

    //The seeds to that are taken as input
    //The order of the seeds will never change, but the vector is not const because the zipcodes
    //decoders may change
    vector<SnarlDistanceIndexClusterer::Seed>& seeds;


    /*
      The tree will represent the seeds' placement in the snarl tree 
      Each node in the tree is either a seed (position on the graph) or the boundary of a snarl
      Edges are labelled with the distance between the two nodes

      This graph is actually represented as a vector of the nodes and edges
      Each item in the vector represents either a node (seed or boundary) or an edge (distance)
      TODO: Fill in a description once it's finalized more
     */

    enum tree_item_t {SEED, SNARL_START, SNARL_END, CHAIN_START, CHAIN_END, EDGE, NODE_COUNT};
    struct tree_item {

        //Is this a seed, boundary, or an edge
        tree_item_t type;

        //For a seed, the index into seeds
        //For an edge, the distance value
        //Empty for a bound
        size_t value;
    };

    //The actual tree structure
    vector<tree_item> zip_code_tree;


    
};
}
#endif
