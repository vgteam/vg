#ifndef VG_SEED_CLUSTERER_HPP_INCLUDED
#define VG_SEED_CLUSTERER_HPP_INCLUDED

#include "snarls.hpp"
#include "distance.hpp"
#include "hash_map.hpp"
#include <structures/union_find.hpp>

namespace vg{
class SnarlSeedClusterer {

    public:

        SnarlSeedClusterer();
        vector<vector<size_t>> cluster_seeds ( vector<pos_t> seeds,
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



        void seed2subtree( const vector<pos_t>& seeds, 
                           const SnarlManager& snarl_manager, 
                           DistanceIndex& dist_index,
                           chains_to_snarl_t& chains_to_snarl, 
                           snarls_to_node_t& snarls_to_node, 
                           node_to_seed_t& node_to_seed);

        hash_set<size_t> get_clusters_node(const vector<pos_t>& seeds, 
                             structures::UnionFind& union_find_clusters,
                             hash_map<size_t, pair<int64_t, int64_t>>& cluster_dists,
                             const node_to_seed_t& node_to_seed,
                             size_t distance_limit, const SnarlManager& snarl_manager,
                             DistanceIndex& dist_index, id_t root, 
                             int64_t node_length);

        hash_set<size_t> get_clusters_chain(const vector<pos_t>& seeds,
                             structures::UnionFind& union_find_clusters,
                             hash_map<size_t, pair<int64_t, int64_t>>& cluster_dists,
                             const chains_to_snarl_t& chains_to_snarl,
                             const snarls_to_node_t& snarls_to_node,
                             const node_to_seed_t& node_to_seed,
                             size_t distance_limit, const SnarlManager& snarl_manager,
                             DistanceIndex& dist_index, const Chain* root);

        hash_set<size_t> get_clusters_snarl(const vector<pos_t>& seeds,
                             structures::UnionFind& union_find_clusters,
                             hash_map<size_t, pair<int64_t, int64_t>>& cluster_dists,
                             const chains_to_snarl_t& chains_to_snarl,
                             const snarls_to_node_t& snarls_to_node,
                             const node_to_seed_t& node_to_seed,
                             size_t distance_limit, const SnarlManager& snarl_manager,
                             DistanceIndex& dist_index, const Snarl* root,
                             bool rev);
};
}

#endif
