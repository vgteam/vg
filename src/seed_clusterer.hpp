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

        //Maps each chain to an ordered map of its component snarls that 
        //contain seeds 
        typedef hash_map< const Chain*, map<size_t, const Snarl*>> 
                                                          chains_to_snarl_t;

        //Maps each snarl to its child nodes that contain seeds
        //nodes are a pair of id and orientation and an indicator for the 
        //node type: 0 = chain, 1 = snarl, 2 = node
        typedef hash_map< const Snarl*, 
                          hash_set<pair<pair<id_t, bool>, int64_t>>> 
                snarls_to_node_t;

        //Maps each node to the indices of seeds
        typedef hash_map< id_t, vector<size_t> > node_to_seed_t;



        vector<const Snarl*> seed2subtree( const vector<pos_t>& seeds, 
                           const SnarlManager& snarl_manager, 
                           DistanceIndex& dist_index,
                           chains_to_snarl_t& chains_to_snarl, 
                           snarls_to_node_t& snarls_to_node, 
                           node_to_seed_t& node_to_seed);

        tuple<hash_set<size_t>, int64_t, int64_t> 
             get_clusters_node(const vector<pos_t>& seeds, 
                             structures::UnionFind& union_find_clusters,
                             vector<pair<int64_t, int64_t>>& cluster_dists,
                             const node_to_seed_t& node_to_seed,
                             size_t distance_limit, const SnarlManager& snarl_manager,
                             DistanceIndex& dist_index, id_t root, 
                             int64_t node_length);

        tuple<hash_set<size_t>, int64_t, int64_t> get_clusters_chain(
                             const vector<pos_t>& seeds,
                             structures::UnionFind& union_find_clusters,
                             vector<pair<int64_t, int64_t>>& cluster_dists,
                             const chains_to_snarl_t& chains_to_snarl,
                             const snarls_to_node_t& snarls_to_node,
                             const node_to_seed_t& node_to_seed,
                             size_t distance_limit, const SnarlManager& snarl_manager,
                             DistanceIndex& dist_index, const Chain* root);

        tuple<hash_set<size_t>, int64_t, int64_t> get_clusters_snarl(
                             const vector<pos_t>& seeds,
                             structures::UnionFind& union_find_clusters,
                             vector<pair<int64_t, int64_t>>& cluster_dists,
                             const chains_to_snarl_t& chains_to_snarl,
                             const snarls_to_node_t& snarls_to_node,
                             const node_to_seed_t& node_to_seed,
                             size_t distance_limit, const SnarlManager& snarl_manager,
                             DistanceIndex& dist_index, const Snarl* root,
                             bool rev);
};
}

#endif
