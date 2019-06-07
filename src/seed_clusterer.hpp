#ifndef VG_SEED_CLUSTERER_HPP_INCLUDED
#define VG_SEED_CLUSTERER_HPP_INCLUDED

#include "snarls.hpp"
#include "min_distance.hpp"
#include "hash_map.hpp"
#include <structures/union_find.hpp>

namespace vg{

class SnarlSeedClusterer {

    public:

        SnarlSeedClusterer(MinimumDistanceIndex& dist_index);

        //Given a vector of seeds (pos_t) and a distance limit, 
        //cluster the seeds such that two seeds whose minimum distance
        //between them (including both of the positions) is less than
        // the distance limit are in the same
        //cluster
        //Returns a vector of clusters. The cluster is a vector of
        //indices into seeds
        vector<vector<size_t>> cluster_seeds ( vector<pos_t> seeds,
               int64_t distance_limit);
    private:
        MinimumDistanceIndex& dist_index;


        enum ChildNodeType {CHAIN, SNARL, NODE};
        
        static inline string typeToString(ChildNodeType t) {
            switch (t) {
            case CHAIN:
                return "CHAIN";
            case SNARL:
                return "SNARL";
            case NODE:
                return "NODE";
            default:
                return "OUT_OF_BOUNDS";
            }
        }

        //child nodes of a snarl's netgraph 
        //size_t is the node id if the node is just a node, index into
        //dist_index's snarl_indexes/chain_index if it is a snarl/chain
        typedef pair<size_t, ChildNodeType> child_node_t;

        //A cluster in the context of a snarl tree node
        //set of the seed indices in the cluster and left and right distance
        //from a seed in the cluster to the ends of the snarl tree node
        //containing the cluster
        typedef tuple<hash_set<size_t>, int64_t, int64_t> child_cluster_t;



        //Find which nodes contain seeds and assign those nodes to the 
        //snarls that contain them. Also find the depth of each snarl
        void get_nodes( const vector<pos_t>& seeds,
                        hash_map<id_t, vector<size_t>>& node_to_seeds,
                        vector<hash_map<size_t, 
                                  vector<pair<child_node_t, child_cluster_t>>>>&
                                                            snarl_to_nodes);

        //Given a node and the indices of seeds on that node, root, 
        //cluster the seeds
        child_cluster_t get_clusters_node(const vector<pos_t>& seeds, 
                             structures::UnionFind& union_find_clusters,
                             vector<pair<int64_t, int64_t>>& cluster_dists,
                             vector<size_t>& seed_indices,
                             int64_t distance_limit, id_t root,
                             int64_t node_length); 

        //Cluster the seeds in a chain
        //snarls_in_chain is an unordered vector of snarls where the snarl
        //  is represented as its index into dist_index.snarl_indexes
        //curr_snarl_children maps each snarl to a vector of its children and 
        //clusters on the children
        child_cluster_t get_clusters_chain(
                             const vector<pos_t>& seeds,
                             structures::UnionFind& union_find_clusters,
                             vector<pair<int64_t, int64_t>>& cluster_dists,
                             vector<size_t>& snarls_in_chain,
                             hash_map<size_t, 
                                  vector<pair<child_node_t, child_cluster_t>>>&
                                                        curr_snarl_children,
                             hash_map<id_t, vector<size_t>>& node_to_seeds,
                             int64_t distance_limit,  size_t chain_index_i);

        //Cluster the seeds in a snarl 
        //child_nodes is a vector of the children of root and their clusters
        child_cluster_t get_clusters_snarl(
                             const vector<pos_t>& seeds,
                             structures::UnionFind& union_find_clusters,
                             vector<pair<int64_t, int64_t>>& cluster_dists,
                             vector<pair<child_node_t, child_cluster_t>>&  
                                                              child_nodes,
                             hash_map<id_t, vector<size_t>>& node_to_seeds,
                             int64_t distance_limit, 
                             size_t snarl_index_i, bool rev) ;

};
}

#endif
