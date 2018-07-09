#include "snarls.hpp"


namespace vg { 

class DistanceIndex {

    /*The distance index*/

    public: 
    //Constructor 
    DistanceIndex (VG* vg, SnarlManager* snarlManager);


    /*Get the distance between two positions
      pos1 must be on a node contained in snarl1 and not on any children of
      snarl1. The same for pos2 and snarl2
    */
    int64_t distance( 
         const Snarl* snarl1, const Snarl* snarl2, pos_t& pos1, pos_t& pos2);
  
    //Helper function to find the minimum value that is not -1
    static int64_t minPos(vector<int64_t> vals);


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
            SnarlDistances(unordered_set<pair<id_t, bool>>& allNodes, 
                          pair<id_t, bool> start, pair<id_t, bool> end);
           

            //Distance between beginning of node start and beginning of node end
            int64_t snarlDistance(pair<id_t, bool> start, pair<id_t, bool> end);

 
            //Distance between end of node start and beginning of node end
            int64_t snarlDistanceShort(VG* graph,NetGraph* ng,
                     pair<id_t, bool> start, pair<id_t, bool> end); 

            //Add the distance from start to end to the index
            void insertDistance(pair<id_t, bool> start, pair<id_t, bool> curr, 
                     int64_t dist);
             
            //Length of a node
            int64_t nodeLength(VG*graph, NetGraph* ng,  pair<id_t, bool> node);
        
            //Total length of the snarl
            int64_t snarlLength();

            /*Given distances from a position to either end of a node, find the
              shortest distance from that position to the start and end nodes of
              the snarl
            */
            pair<int64_t, int64_t> distToEnds(id_t node, bool rev, 
                                                 int64_t distL, int64_t distR);

            void printSelf();

        protected:

            //Maps node to index to get its distance
            unordered_map< pair<id_t, bool>, size_t> visitToIndex;
 
             //Store the distance between every pair nodes, -1 indicates no path
             //For child snarls that are unary or only connected to one node
             //in the snarl, distances between that node leaving the snarl
             //and any other node is -1
            vector<int64_t> distances;

            //ID of the first node in the snarl, also key for distance index 
            pair<id_t, bool> snarlStart;
 
           //End facing out of snarl
            pair<id_t, bool> snarlEnd;           

            //Total length of the snarl- start to end plus length of end
            int64_t length; 

            //The index into distances for distance start->end
            size_t index(pair<id_t, bool> start, pair<id_t, bool> end);

        private: 

            int64_t snarlDistanceShortHelp(VG* graph,NetGraph* ng,
                     pair<id_t, bool> start, pair<id_t, bool> end); 


        friend class DistanceIndex;
        friend class TestDistanceIndex;
    };

    class ChainDistances {
        /*Stores distances between snarls in a chain*/

        public:
        
            //Constructor
            ChainDistances(unordered_map<id_t, size_t> s, vector<int64_t> p,
                        vector<int64_t> fd, vector<int64_t> rev );
       
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
            unordered_map<id_t, size_t> snarlToIndex; 

            /*Dist from start of chain to start and end of each boundary node of
              all snarls in the chain*/
            vector<int64_t> prefixSum;

            /*For each boundary node of snarls in the chain, the distance
               from the start of the node traversing forward to the end of 
               the same node traversing backwards*/
            vector<int64_t> loopFd;
    
            /*For each boundary node of snarls in the chain, the distance
               from the end of the node traversing backward to the start of 
               the same node traversing forward*/
            vector<int64_t> loopRev;
        
            /*Helper function for finding distances*/
            int64_t chainDistanceHelper(pair<size_t, bool> start, 
                                       pair<size_t, bool> end     );

            /*Returns true if the snarl is reversed in the chain*/
            bool isReverse(const Snarl* snarl, SnarlManager* sm);

        friend class DistanceIndex;   
        friend class TestDistanceIndex;
    }; 


    //map from start node of a snarl to its index
    unordered_map<pair<id_t, bool>, SnarlDistances> snarlIndex;

    //map from node id of first node in snarl to that chain's index
    unordered_map<id_t, ChainDistances> chainIndex;

    //Graph and snarl manager for this index
    VG* graph;

    SnarlManager* sm;
 
    //Helper function for constructor
    int64_t calculateIndex(const Chain* chain); 

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
