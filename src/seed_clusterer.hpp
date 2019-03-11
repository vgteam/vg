#include "snarls.hpp"
#include "distance.hpp"
#include "hash_map.hpp"

namespace vg{
class SnarlSeedClusterer {

    public:

        SnarlSeedClusterer();
        vector<hash_set<size_t>> cluster_seeds ( vector<pos_t> seeds,
               size_t distance_limit, SnarlManager& snarl_manager,
               DistanceIndex& dist_index);
    private:

        struct seed_index_list {
        /* linked list item that contains the index of a seed and a pointer to
         * the next item in the list
        */
            size_t index;
            seed_index_list* next = NULL;
 
        };
        struct cluster_t {
            //Minimum distance from any seed to the end of the snarl or chain
            //containing this cluster
            int64_t dist_left = -1;
            int64_t dist_right = -1;

            //Start and end of list containing seeds
            seed_index_list* seeds_start = NULL;
            seed_index_list* seeds_end = NULL;
        };

        vector<cluster_t> get_clusters_node(vector<pos_t>& seeds, 
                             vector<seed_index_list>& seed_list,
                             hash_map< id_t, vector<size_t> >& node_to_seed,
                             size_t distance_limit, SnarlManager& snarl_manager,
                             DistanceIndex& dist_index, id_t root, bool rev,
                             int64_t node_length);

        vector<cluster_t> get_clusters_chain(vector<pos_t>& seeds,
                             vector<seed_index_list>& seed_list,
                             hash_set< const Chain* >& chains_with_seeds,
                             hash_set< const Snarl* >& snarls_with_seeds,
                             hash_map< id_t, vector<size_t> >& node_to_seed,
                             size_t distance_limit, SnarlManager& snarl_manager,
                             DistanceIndex& dist_index, const Chain* root);

        vector<cluster_t> get_clusters_snarl(vector<pos_t>& seeds,
                             vector<seed_index_list>& seed_list,
                             hash_set< const Chain* >& chains_with_seeds,
                             hash_set< const Snarl* >& snarls_with_seeds,
                             hash_map< id_t, vector<size_t> >& node_to_seed,
                             size_t distance_limit, SnarlManager& snarl_manager,
                             DistanceIndex& dist_index, const Snarl* root,
                             bool include_boundaries);
};
}
