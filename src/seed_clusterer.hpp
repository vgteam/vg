#ifndef VG_SEED_CLUSTERER_HPP_INCLUDED
#define VG_SEED_CLUSTERER_HPP_INCLUDED

#include "snarls.hpp"
#include "distance.hpp"
#include "hash_map.hpp"
#include <structures/union_find.hpp>

namespace vg{

class SnarlSeedClusterer {

    public:

        SnarlSeedClusterer(SnarlManager& snarl_manager, 
                           DistanceIndex& dist_index);

        //Given a vector of seeds (pos_t) and a distance limit, 
        //cluster the seeds such that two seeds whose minimum distance
        //between them is less than the distance limit are in the same
        //cluster
        //Returns a vector of clusters. The cluster is a vector of
        //indices into seeds
        vector<vector<size_t>> cluster_seeds ( vector<pos_t> seeds,
               size_t distance_limit);
    private:
        SnarlManager& snarl_manager;
        DistanceIndex& dist_index;


        enum ChildNodeType {CHAIN, SNARL, NODE};

        //children of a snarl. 
        //pair<id_t, bool> is the node id and orientation
        // int64_t is 0 if its a chain, 1 for snarl, 2 for node
        typedef pair<pair<id_t, bool>, ChildNodeType> child_node_t;

        //A cluster in the context of a snarl tree node
        //set of the seed indices in the cluster and left and right distance
        //from a seed in the cluster to the ends of the snarl tree node
        //containing the cluster
        typedef tuple<hash_set<size_t>, int64_t, int64_t> child_cluster_t;



        //Find which nodes contain seeds and assign those nodes to the 
        //snarls that contain them. Also find the depth of each snarl
        void get_nodes( const vector<pos_t>& seeds,
                        hash_map<id_t, vector<size_t>>& node_to_seeds,
                        vector<hash_map<const Snarl*, 
                                  vector<pair<child_node_t, child_cluster_t>>>>&
                                                            snarl_to_nodes);

        //Given a node and the seeds on that node, cluster the seeds
        child_cluster_t get_clusters_node(const vector<pos_t>& seeds, 
                             structures::UnionFind& union_find_clusters,
                             vector<pair<int64_t, int64_t>>& cluster_dists,
                             vector<size_t>& seed_indices,
                             size_t distance_limit, id_t root,
                             int64_t node_length); 

        //Cluster the seeds in a chain
        child_cluster_t get_clusters_chain(
                             const vector<pos_t>& seeds,
                             structures::UnionFind& union_find_clusters,
                             vector<pair<int64_t, int64_t>>& cluster_dists,
                             vector<pair<size_t, pair<const Snarl*, 
                              DistanceIndex::SnarlIndex*>>>& snarls_in_chain,
                             hash_map<const Snarl*, 
                                  vector<pair<child_node_t, child_cluster_t>>>&
                                                        curr_snarl_children,
                             hash_map<id_t, vector<size_t>>& node_to_seeds,
                             size_t distance_limit,  const Chain* root);

        //Cluster the seeds in a snarl 
        child_cluster_t get_clusters_snarl(
                             const vector<pos_t>& seeds,
                             structures::UnionFind& union_find_clusters,
                             vector<pair<int64_t, int64_t>>& cluster_dists,
                             vector<pair<child_node_t, child_cluster_t>>&  
                                                              child_nodes,
                             hash_map<id_t, vector<size_t>>& node_to_seeds,
                             size_t distance_limit, 
                             DistanceIndex::SnarlIndex& snarl_index, 
                             const Snarl* root,
                             bool rev) ;

};
}

#endif
