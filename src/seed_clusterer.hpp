#include "snarls.hpp"
#include "distance.hpp"
#include "hash_map.hpp"
#include <structures/union_find.hpp>

namespace vg{
class SnarlSeedClusterer {

    public:

        SnarlSeedClusterer();
        vector<hash_set<size_t>> cluster_seeds ( vector<pos_t> seeds,
               size_t distance_limit, SnarlManager& snarl_manager,
               DistanceIndex& dist_index);
    private:

        //Maps each chain to an ordered vector of its component snarls that contain
        //seeds 
        typedef hash_map< const Chain*, vector<const Snarl*>> chains_to_snarl_t;

        //Maps each snarl to its child nodes that contain seeds
        typedef hash_map< const Snarl*, hash_set<pair<id_t, bool>>> snarls_to_node_t;

        //Maps each node to the indices of seeds
        typedef hash_map< id_t, vector<size_t> > node_to_seed_t;

        struct seed_index_list {
        /* linked list item that contains the index of a seed and a pointer to
         * the next item in the list
        */
            size_t index;
            seed_index_list* next = nullptr;
 
        };
        struct cluster_t {
            //Minimum distance from any seed to the end of the snarl or chain
            //containing this cluster
            int64_t dist_left = -1;
            int64_t dist_right = -1;

            //Start and end of list containing seeds
            seed_index_list* seeds_start = nullptr;
            seed_index_list* seeds_end = nullptr;
        };

        void seed2subtree( const vector<pos_t>& seeds, 
                           const SnarlManager& snarl_manager, 
                           DistanceIndex& dist_index,
                           chains_to_snarl_t& chains_to_snarl, 
                           snarls_to_node_t& snarls_to_node, 
                           node_to_seed_t& node_to_seed,
                           vector<seed_index_list>& seed_list);

        vector<cluster_t> get_clusters_node(const vector<pos_t>& seeds, 
                             vector<seed_index_list>& seed_list,
                             const node_to_seed_t& node_to_seed,
                             size_t distance_limit, const SnarlManager& snarl_manager,
                             DistanceIndex& dist_index, id_t root, 
                             int64_t node_length);

        vector<cluster_t> get_clusters_chain(const vector<pos_t>& seeds,
                             vector<seed_index_list>& seed_list,
                             const chains_to_snarl_t& chains_to_snarl,
                             const snarls_to_node_t& snarls_to_node,
                             const node_to_seed_t& node_to_seed,
                             size_t distance_limit, const SnarlManager& snarl_manager,
                             DistanceIndex& dist_index, const Chain* root);

        vector<cluster_t> get_clusters_snarl(const vector<pos_t>& seeds,
                             vector<seed_index_list>& seed_list,
                             const chains_to_snarl_t& chains_to_snarl,
                             const snarls_to_node_t& snarls_to_node,
                             const node_to_seed_t& node_to_seed,
                             size_t distance_limit, const SnarlManager& snarl_manager,
                             DistanceIndex& dist_index, const Snarl* root,
                             bool rev);
};
}
