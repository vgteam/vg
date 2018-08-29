#include "snarls.hpp"
#include "hash_map.hpp"
using namespace sdsl;
namespace vg { 

class DistanceIndex {

    /*The distance index. Used for calculation of the minimum distance between
       two positions and for a maximum distance estimation. The maximum distance
       estimation is at least as large as the maximum distance between two 
       positions up to a specified cap*/

    public: 
    //Constructor 
    DistanceIndex (VG* vg, SnarlManager* snarlManager);

    //Constructor to load index from serialization 
    DistanceIndex (VG* vg, SnarlManager* snarlManager, istream& in);
  
    //Serialize object into out
    void serialize(ostream& out);

    //Load serialized object from in
    void load(istream& in);

    int64_t sizeOf();

    /*Get the distance between two positions
      pos1 must be on a node contained in snarl1 and not on any children of
      snarl1. The same for pos2 and snarl2
    */
    int64_t distance( pos_t& pos1, pos_t& pos2);
    int64_t distance( 
         const Snarl* snarl1, const Snarl* snarl2, pos_t& pos1, pos_t& pos2);

  
    //Helper function to find the minimum value that is not -1
    static int64_t minPos(vector<int64_t> vals);


    //Given a node, find the snarl containing it
    const Snarl* snarlOf(id_t nodeID);

    protected:
    void printSelf();
    class SnarlDistances {
        
        /* Stores distance information for nodes in a snarl.
           visitToIndex maps each visit (node_id, negative if reverse) to an int
           distances stores all the distance bewteen each pair of visits in a 
           snarl
        */

        public:
        
            //Constructor
            SnarlDistances(DistanceIndex* di,
                          unordered_set<pair<id_t, bool>>& allNodes, 
                          pair<id_t, bool> start, pair<id_t, bool> end);
           
            //Construct from vector - inverse of toVector
            SnarlDistances(DistanceIndex* di, vector<int64_t> v);

            /*Store contents of object as a vector of ints for serialization
              Stored as [# nodes, start node id, end node id, snarl length] + 
                        [visit to index as list of node ids in order of index] +
                        [distances]
            */
            vector<int64_t>  toVector();
            
            //Distance between beginning of node start and beginning of node end
            int64_t snarlDistance(VG* graph,NetGraph* ng,pair<id_t, bool> start,
                                                         pair<id_t, bool> end);

 
            //Distance between end of node start and beginning of node end
            int64_t snarlDistanceShort(pair<id_t, bool> start, 
                                       pair<id_t, bool> end); 

            //Add the distance from start to end to the index
            void insertDistance(pair<id_t, bool> start, pair<id_t, bool> curr, 
                     int64_t dist);
             
            //Length of a node
            int64_t nodeLength(VG*graph, NetGraph* ng,  pair<id_t, bool> node);
        
            //Total length of the snarl
            int64_t snarlLength(VG* graph, NetGraph* ng);

            /*Given distances from a position to either end of a node, find the
              shortest distance from that position to the start and end nodes of
              the snarl
            */
            pair<int64_t, int64_t> distToEnds(VG* graph, NetGraph* ng, 
                             id_t node, bool rev, int64_t distL, int64_t distR);

            void printSelf();

        protected:

            //Maps node to index to get its distance
            hash_map< pair<id_t, bool>, size_t> visitToIndex;
 
             /*Store the distance between every pair nodes, -1 indicates no path
             For child snarls that are unary or only connected to one node
             in the snarl, distances between that node leaving the snarl
             and any other node is -1
             Distance from a node to itself is -1 unless there is a path leaving
             that node and reaching it again
             Distances are stored with an offset of 1 (-1 dist stored as 0, 0
               distance stored as 1)
             */
            int_vector<> distances;

            //ID of the first node in the snarl, also key for distance index 
            pair<id_t, bool> snarlStart;
 
           //End facing out of snarl
            pair<id_t, bool> snarlEnd;           

            //The index into distances for distance start->end
            size_t index(pair<id_t, bool> start, pair<id_t, bool> end);

        private: 
            DistanceIndex* distIndex; 



        friend class DistanceIndex;
        friend class TestDistanceIndex;
    };

    class ChainDistances {
        /*Stores distances between snarls in a chain*/

        public:
        
            //Constructor
            ChainDistances(hash_map<id_t, size_t> s, vector<int64_t> p,
                        vector<int64_t> fd, vector<int64_t> rev );

            //Constructor from vector of ints after serialization
            ChainDistances(vector<int64_t> v);

            /*Convert contents into vector of ints for serialization
               stored as [node_id1, prefixsum1 start, prefixsum1 end,
                                            loopfd1, loopfd2, node_id2, ...]
            */
            vector<int64_t> toVector();
       
            /*Distance between two snarls starting from the beginning of the 
              start node to the beginning of the end node.
              bool is true if traversing reverse relative to the start of 
              the chain */
            int64_t chainDistance(pair<id_t, bool> start, pair<id_t, bool> end);

            /*Distance between two snarls starting from the beginning of 
              the node after start to the beginning of end */ 
            int64_t chainDistanceShort(VG* graph, pair<id_t, bool> start, 
                                                        pair<id_t, bool> end);
            //Length of entire chain
            int64_t chainLength();

            /* Given the distance from a position to either end of a snarl in 
               the chain, find the shortest distance from the position to 
               either end of the chain
            */ 
            pair<int64_t, int64_t> distToEnds(pair<id_t, bool> start,
                                              int64_t distL, int64_t distR);
            void printSelf();

        protected:

            hash_map<id_t, size_t> snarlToIndex; 

            /*Dist from start of chain to start and end of each boundary node of
              all snarls in the chain*/
            int_vector<> prefixSum;

            /*For each boundary node of snarls in the chain, the distance
               from the start of the node traversing forward to the end of 
               the same node traversing backwards*/
            int_vector<> loopFd;
    
            /*For each boundary node of snarls in the chain, the distance
               from the end of the node traversing backward to the start of 
               the same node traversing forward*/

            int_vector<> loopRev;


        
            /*Helper function for finding distances*/
            int64_t chainDistanceHelper(pair<size_t, bool> start, 
                                       pair<size_t, bool> end     );

            /*Returns true if the snarl is reversed in the chain*/
            bool isReverse(const Snarl* snarl, SnarlManager* sm);

        friend class DistanceIndex;   
        friend class TestDistanceIndex;
    }; 

    ///////// Data members of overall index

    //map each node to connected component for max distance estimation 
    hash_map<id_t, size_t> nodeToCycles;

    //map from start node of a snarl to its index
    unordered_map<pair<id_t, bool>, SnarlDistances> snarlIndex;

    //map from node id of first node in snarl to that chain's index
    hash_map<id_t, ChainDistances> chainIndex;

    //Graph and snarl manager for this index
    VG* graph;

    SnarlManager* sm;


    /*Index to find the snarl containing a node
      The start node id of the snarl containing each node - negative if 
       the start node is reverse
    TODO: Maybe put this somewhere else*/
    dac_vector<> nodeToSnarl;
    id_t minNodeID; //minimum node id of the graph




    ////// Private helper functions
 




    //Helper function for constructor - populate the minimum distance index
    int64_t calculateMinIndex(const Chain* chain); 

    //Helper function for constructor - populate the minimum distance index
    void calculateMaxIndex(const Chain* chain, int64_t cap); 

    //Helper function for constructor - populate node to snarl
    int_vector<> calculateNodeToSnarl(VG* vg, SnarlManager* sm);

    //Flag each node with true if it is in a cycle that has minimum length
    //smaller than cap
    void flagCycles(const Snarl* snarl, bit_vector& inCycle,
                    int64_t cap);
    //Assign each node to a connected component of cycles
    pair<int_vector<>, uint64_t> findCycleComponents();

    //Minimum distance of a loop that involves node
    int64_t loopDistance(pair<id_t, bool> node1, pair<id_t, bool> node2); 
    int64_t loopDistance( const Snarl* snarl1, const Snarl* snarl2,
                          pair<id_t, bool> node1, pair<id_t, bool> node2); 

    /*Helper function for distance calculation
      Returns the distance to the start of and end of the child snarl of
      common ancestor containing snarl, commonAncestor if snarl is
      the common ancestor
    */
    pair<pair<int64_t, int64_t>, const Snarl*> distToCommonAncestor(
                const Snarl* snarl, const Snarl* commonAncestor, pos_t& pos); 



    // Methods for testing
    int64_t checkChainDist(id_t snarl, size_t index);
    int64_t checkChainLoopFd(id_t snarl, size_t index);
    int64_t checkChainLoopRev(id_t snarl, size_t index);
    friend class SnarlDistances;
    friend class ChainDistances;


};
 
}
