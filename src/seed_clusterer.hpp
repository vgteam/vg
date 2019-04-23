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

        //list of children of a snarl. 
        //pair<id_t, bool> is the node id and orientation
        // int64_t is 0 if its a chain, 1 for snarl, 2 for node
        typedef vector<pair<pair<id_t, bool>, int64_t>> child_node_list;



        void get_nodes( const vector<pos_t>& seeds, 
                        const SnarlManager& snarl_manager, 
                        DistanceIndex& dist_index,
                        hash_map<id_t, vector<size_t>>& node_to_seeds,
                        vector<hash_map<const Snarl*, child_node_list>>& 
                                                            snarl_to_nodes);

        tuple<hash_set<size_t>, int64_t, int64_t> 
             get_clusters_node(const vector<pos_t>& seeds, 
                             structures::UnionFind& union_find_clusters,
                             vector<pair<int64_t, int64_t>>& cluster_dists,
                             vector<size_t>& seed_indices,
                             size_t distance_limit, 
                             const SnarlManager& snarl_manager, id_t root,
                             int64_t node_length); 

        tuple<hash_set<size_t>, int64_t, int64_t> get_clusters_chain(
                             const vector<pos_t>& seeds,
                             structures::UnionFind& union_find_clusters,
                             vector<pair<int64_t, int64_t>>& cluster_dists,
                             vector<pair<size_t, pair<const Snarl*, 
                              DistanceIndex::SnarlIndex*>>>& snarls_in_chain,
                             hash_map<const Snarl*, child_node_list>& 
                                                        curr_snarl_children,
                             hash_map<pair<id_t, bool>,
                                 tuple<hash_set<size_t>, int64_t, int64_t>>&
                                  child_clusters,
                             hash_map<id_t, vector<size_t>>& node_to_seeds,
                             size_t distance_limit, 
                             const SnarlManager& snarl_manager,
                             DistanceIndex& dist_index, const Chain* root);

        tuple<hash_set<size_t>, int64_t, int64_t> get_clusters_snarl(
                             const vector<pos_t>& seeds,
                             structures::UnionFind& union_find_clusters,
                             vector<pair<int64_t, int64_t>>& cluster_dists,
                             child_node_list&  child_nodes,
                             hash_map<pair<id_t, bool>,
                                 tuple<hash_set<size_t>, int64_t, int64_t>>&
                                                         child_clusters_set,
                             hash_map<id_t, vector<size_t>>& node_to_seeds,
                             size_t distance_limit, 
                             const SnarlManager& snarl_manager,
                             DistanceIndex::SnarlIndex& snarl_index, const Snarl* root,
                             bool rev) ;

};
}

#endif
