#include "snarls.hpp"

namespace vg {

    class ChainDistances {
    /*Stores distances between the boundary nodes of two snarls in a chain*/
    public:
        
        //Constructor
        ChainDistances(unordered_map<id_t, size_t>* s, vector<int64_t>* p);
       
        /*Distance between two snarls starting from the beginning of the start
          node to the beginning of the end node */
        int64_t chainDistance(pair<id_t, bool> start, pair<id_t, bool> end);

        /*Distance between two snarls starting from the beginning of the node
          after start to the beginning of end */ 
        int64_t chainDistance_short(VG* graph, pair<id_t, bool> start, 
                                                        pair<id_t, bool> end);
        //Length of entire chain
        int64_t chainLength();

    private:
         unordered_map<id_t, size_t>* snarl_to_index; 
         vector<int64_t>* prefix_sum;
         /*Dist from start of chain to start of each snarl, last element is 
            length of entire chain*/
    };    

    class SnarlDistances {
        
    /* Stores distance information for nodes in a snarl.
       visit_to_index maps each visit (node_id, negative if reverse) to an int
       distances stores all the distance bewteen each pair of visits in a snarl
    */

    public:
        
        //Constructor
        SnarlDistances(unordered_map<pair<id_t, bool>, size_t> * vti,
                       vector<int64_t>* dists);

        //Distance between beginning of node start and beginning of node end
        int64_t snarlDistance(pair<id_t, bool> start, pair<id_t, bool> end);

 
        //Distance between end of node start and beginning of node end
        int64_t snarlDistance_short(VG* graph, pair<id_t, bool> start, 
                                                 pair<id_t, bool> end); 

        //Add the distance from start to end
        void insertDistance(pair<id_t, bool> start, pair<id_t, bool> curr, 
                 int64_t dist);

    private:
        unordered_map< pair<id_t, bool>, size_t>* visit_to_index;
        vector<int64_t>* distances;

        //The index into distances at which the distance from start to end is
        size_t index(pair<id_t, bool> start, pair<id_t, bool> end);
    };


    struct DistanceIndex {

        /*The distance index containing snarl and chain information*/
        unordered_map<id_t, SnarlDistances*> sd;
             //map from node id of first nodein snarl
        unordered_map<id_t, ChainDistances*> cd;
    };
      
    
    DistanceIndex makeDistanceIndex(VG* graph, SnarlManager* sm);
    
}
