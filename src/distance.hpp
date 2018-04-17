#include "snarls.hpp"

namespace vg {
    
    struct ChainDistances {
        /*Stores distance information for snarls in a chain*/
        unordered_map<id_t, size_t>* snarl_to_index;
        vector<int64_t>* prefix_sum;//Dist from start of chain to each snarl
    };


    struct SnarlDistances {
        
        /* Stores distance information for nodes in a snarl.
        visit_to_index maps each visit (node_id, negative if reverse) to an int
        distances stores all the distance bewteen each pair of visits in a snarl
        */
        unordered_map<pair<id_t, bool>, size_t>* visit_to_index;
        vector<int64_t>* distances;
    };


    struct DistanceIndex {

        /*The distance index containing snarl and chain information*/
        unordered_map<id_t, SnarlDistances*> sd;
             //map from node id of first nodein snarl
        unordered_map<id_t, ChainDistances*> cd;
    };
      
    
    DistanceIndex makeDistanceIndex(VG* graph, SnarlManager& sm);
}
