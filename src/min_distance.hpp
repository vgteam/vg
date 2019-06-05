#ifndef VG_MIN_DISTANCE_HPP_INCLUDED
#define VG_MIN_DISTANCE_HPP_INCLUDED

#include "snarls.hpp"
#include "hash_map.hpp"
#include "hash_graph.hpp"
#include "split_strand_graph.hpp"
#include "algorithms/dagify.hpp"
#include "algorithms/split_strands.hpp"
#include "algorithms/topological_sort.hpp"
#include "algorithms/is_acyclic.hpp"
#include "algorithms/is_single_stranded.hpp"

using namespace sdsl;
namespace vg { 

class MinimumDistanceIndex {

    /**
     * The distance index. Stores minimum distances among nodes in each 
     * netgraph and chain. Used for calculation of the minimum distance between
     * two positions and for a maximum distance estimation. The maximum distance
     * estimation is at least as large as the maximum distance between two 
     * positions up to a specified cap
     */

    public: 
    //Constructor 
    //Cap is the distance up to which the maximum distance will give a reliable bound - if there is a path with length greater than cap, then cap may be returned
    //If the cap is set to 0 (default), then the maximum distance index is not
    //included
    MinimumDistanceIndex (const HandleGraph* graph, const SnarlManager* snarl_manager,
                            int64_t cap = 0);

    
    //Constructor to load index from serialization
    MinimumDistanceIndex (istream& in);
    
    //Default constructor; load() must be called next.
    MinimumDistanceIndex ();

    //Serialize object into out
    void serialize(ostream& out) const;

    //Load serialized object from in. Does not rely on the internal graph or 
    //snarl manager pointers.
    void load(istream& in);
    
    //Get the minimum distance between two positions
    // Distance includes only one of the positions. The distance from a 
    // position to itself would be 1
    //If there is no path between the two positions then the distance is -1
    
    int64_t minDistance( pos_t pos1, pos_t pos2);

    //Get a maximum distance bound between the positions, ignoring direction
    //Returns a positive value even if the two nodes are unreachable
    int64_t maxDistance(pos_t pos1, pos_t pos2);

    //Helper function to find the minimum value that is not -1
    static int64_t minPos(vector<int64_t> vals);

    /*print the distance index for debugging*/
    void printSelf();
    /*Prints the number of nodes in each snarl netgraph and number of snarls in each chain*/
    void printSnarlStats();

    protected:
    class SnarlIndex {
        
        /* Stores distance information for nodes in a snarl.
           visitToIndex maps each visit (node_id, negative if reverse) to an int
           distances stores all the distance bewteen each pair of visits in a 
           snarl
        */

        public:
        
            //Constructor
            //if inChain is true, parentIndex and revInParent are for a chain
            //otherwise, for the parent snarl
            SnarlIndex(id_t parent_id, bool rev_in_parent, 
                       id_t id_in_parent, bool is_unary_snarl, size_t depth,
                       size_t num_nodes, bool in_chain);
           
            //Construct and empty SnarlIndex. Must call load after construction to populate it 
            SnarlIndex();

            //Load data from serialization
            void load(istream& in);

            /*Serialize the snarl
              Stored as start node id + end node id + 
                        [visit to index as list of node ids in order of index] +
                        [distances]
            */
            void serialize(ostream& out) const;
            
            ///Distance between start and end, not including the lengths of
            ///the two nodes
            ///start and end are the ranks of the node+direction, given by
            ///primary_snarls and secondary_snarl
            ///Only works for nodes heading their chains (which represent the 
            ///chains), or snarl boundaries.
            ///Rank 0 is the start node and rank num_nodes*2-1 is the end node.
            int64_t snarlDistance(size_t start, size_t end);


            //Add the distance from start to end to the index
            void insertDistance(size_t start, size_t end, int64_t dist);
             
            //Length of a node in the netgraph of the snarl
            //If it is a node, then the length of the node. If it is a snarl or
            //chain, then the shortest distance between the boundaries
            // i is the index of the node in the snarl
            int64_t nodeLength(size_t i);
        
            //Total length of the snarl-shortest distance from start to end
            //including the lengths of boundary nodes
            int64_t snarlLength();

            /*Given distances from a position to either end of a node, find the
              shortest distance from that position to the start and end nodes of
              the snarl
              node is the index of the node in the forward direction
            */
            pair<int64_t, int64_t> distToEnds(size_t rank, 
                                              int64_t distL, int64_t distR);

            void printSelf();

        protected:
 
             /*Store the distance between every pair nodes, not including the 
             lengths of the nodes. 
             The lengths of each of the nodes are stored as the first n entries
             Distances stored are 1 greater than actual distances 
             -1 (stored as 0) indicates no path
             For child snarls that are unary or only connected to one node
             in the snarl, distances between that node leaving the snarl
             and any other node is -1
             
             */
            int_vector<> distances;

            bool in_chain; //True if this snarl is in a chain

            /// id of the parent snarl or chain of this snarl 
            //0 if this is a top level snarl with no chain
            id_t parent_id;
            bool rev_in_parent;
            //Id of this snarl in the parent. If the parent is a chain, then 
            //the id of the boundary node that occurs first in the order of 
            //the chain
            id_t id_in_parent;

            //Number of nodes in the snarl
            size_t num_nodes;

            //Depth in the snarl tree - 0 for root
            size_t depth;
            
            //True if this snarl is a unary snarl
            //Since the start and end node are the same, the last ranking
            //node is no longer the end node
            bool is_unary_snarl;

            //The index into distances for distance start->end
            size_t index(size_t start, size_t end);


        friend class MinimumDistanceIndex;
        friend class SnarlSeedClusterer;
        friend class TestMinDistanceIndex;
    };

    class ChainIndex {
        /*Stores distances between snarls in a chain*/

        public:
        
            //Constructor
            //Takes the index into snarlIndexes of the parent snarl, 
            //whether it is reversed in the snarl, and the chain length
            ChainIndex(size_t parent_id, size_t id_in_parent, bool rev_in_parent,
                       bool loops, size_t length);

            //Constructor from vector of ints after serialization
            ChainIndex();
            //Load data from serialization
            void load(istream& in);

            /*Serialize the chain
             * stored as chainStartID +  chainEndID,
             * node_id1, prefixsum1 start,
                           prefixsum1 end, loopfd1, loopfd2, node_id2, ...]
            */
            void serialize(ostream& out) const;
       
            /** 
             * Distance between two node sides in a chain. size_t values specify
             * the nodes - rank of the node in the chain, 
             * bools specify the side and orientation.
             *
             * Bool specifies the direction the node is traversed i.e. if both
             * bools are false, then the distance is between the left side of
             * the start node traversing it forward to the left side of the 
             * end node traversing forward. 
             *
             * If the chain loops, then chainDistance must be called twice to 
             * include the distance by taking the loop 
             */
            int64_t chainDistance(pair<size_t, bool> start, 
                                  pair<size_t, bool> end, int64_t startLen, 
                                  int64_t endLen, bool check_loop=false);

            int64_t loopDistance(pair<size_t, bool> start, 
                                  pair<size_t, bool> end, int64_t startLen, 
                                  int64_t endLen);

            //Length of entire chain
            int64_t chainLength();

            void printSelf();

        protected:


            /*Dist from start of chain to start and end of each boundary node of
              all snarls in the chain
              The first value should be 0 according to this scheme but it is
              the length of the first node in the chain. Similarly, an extra
              value is stored at the end of the vector that is the length of the
              entire chain*/
            int_vector<> prefix_sum;

            /*For each boundary node of snarls in the chain, the distance
               from the start of the node traversing forward to the end of 
               the same node traversing backwards -directions relative to the 
               direction the node is traversed in the chain*/
            int_vector<> loop_fd;
    
            /*For each boundary node of snarls in the chain, the distance
               from the end of the node traversing backward to the start of 
               the same node traversing forward*/
            int_vector<> loop_rev;

            /// id of parent snarl of the chain 
            //0 if top level chain
            id_t parent_id;
            bool rev_in_parent;
            id_t id_in_parent; //Id of the start node of this chain 
            bool is_looping_chain; //True if the chain loops



        friend class MinimumDistanceIndex;   
        friend class SnarlSeedClusterer;
        friend class TestMinDistanceIndex;
    }; 


    ///Class used to find the maximum distance between two positions in the 
    //graph
    //Disregards direction
    class MaxDistanceIndex{

        public:

            MaxDistanceIndex();

            void serialize(ostream& out) const;
            void load(istream& in);

            //Get the maximum distance between two nodes
            int64_t maxDistance(id_t id1, id_t id2);
            

        private:
            //For each node in the graph, store the minimum and maximum
            //distances from a tip to the node
            sdsl::int_vector<> min_distances;
            sdsl::int_vector<> max_distances;

            void calculateMaxIndex(const HandleGraph* graph, int64_t cap); 


            


        friend class MinimumDistanceIndex;
    };

    ///////// Data members of overall distance index


    //vector of all snarl indexes
    vector<SnarlIndex> snarl_indexes;

    //vector of all chain indexes
    vector< ChainIndex> chain_indexes;


    //Each of the ints in these vectors are offset by 1: 0 is stored as 1, etc.
    //This is so that we can store -1 as 0 instead of int max

    //Vector of length max node id - min node id
    //For each node, stores the index into snarlIndexes for the snarl
    //containing the node and the rank of the node in the snarlIndex 
    //used for looking up distances in the matrix. The index for
    //the distance matrix is always for the fd direction, rev direction
    //is the index + 1. For the start and end nodes of the snarl, the inward
    //pointing node is stored as the first and last elements in the distance
    //matrix
    sdsl::int_vector<> primary_snarl_assignments;
    sdsl::int_vector<> primary_snarl_ranks;

    //Similar to primary snarls, stores snarl and index of secondary snarl
    //each node belongs to, if any.
    //Secondary snarl can be a node that represents a snarl/chain in the
    //netgraph of the parent snarl or a node that participates in multiple
    //snarls in a chain. The primary snarl will always
    //be the snarl that occurs first in the chain
    sdsl::int_vector<> secondary_snarl_assignments;
    sdsl::int_vector<> secondary_snarl_ranks;
    
    //For each node, stores 1 if the node is in a secondary snarl and 0
    //otherwise. Use rank to find which index into secondary_snarls
    //a node's secondary snarl is at
    sdsl::bit_vector has_secondary_snarl_bv;
    sdsl::rank_support_v<1> has_secondary_snarl;

    //For each node, store the index and rank for the chain that the node
    //belongs to, if any
    sdsl::int_vector<> chain_assignments;
    sdsl::int_vector<> chain_ranks;
    sdsl::bit_vector has_chain_bv;
    sdsl::rank_support_v<1> has_chain;

    id_t min_node_id; //minimum node id of the graph
    id_t max_node_id; //maximum node id of the graph


    //The depth of the snarl tree, starting from 0
    size_t tree_depth;


    //True if we are including the maximum distance index
    bool include_maximum;
    //Index for maximum distances
    MaxDistanceIndex max_index;



    ////// Private helper functions
 




    //Helper function for constructor - populate the minimum distance index
    //Given the top level snarls
    int64_t calculateMinIndex(const HandleGraph* graph, 
                      const SnarlManager* snarl_manager, const Chain* chain, 
                       size_t parent_id, bool rev_in_parent, 
                       bool trivial_chain, size_t depth); 


    /*Helper function for distance calculation
      Returns the distance to the start of and end of the child snarl of
      common ancestor containing snarl, commonAncestor if snarl is
      the common ancestor and the node id of the node in the common ancestor.
      The node in the common ancestor is a pair of <node id, rev>
      rev is false if the pos is the start pos ( if it must be traversed 
      forward) and false if it is the end pos (if it must be reached in 
      the direction of pos)
    */
    tuple<int64_t, int64_t, pair<id_t, bool>> distToCommonAncestor(
                pair<size_t, bool> common_ancestor, pos_t& pos, bool rev); 


    //Get the index into chain_indexes/rank in chain of node i
    size_t getPrimaryAssignment(id_t i);
    size_t getPrimaryRank(id_t i);

    size_t getChainAssignment(id_t i);
    size_t getChainRank(id_t i);

    size_t getSecondaryAssignment(id_t i);
    size_t getSecondaryRank(id_t i);


    friend class SnarlIndex;
    friend class ChainIndex;
    friend class SnarlSeedClusterer;
    friend class TestMinDistanceIndex;


};
 
}

#endif
