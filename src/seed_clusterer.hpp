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
        //node_id is the node id if the node is just a node, index into
        //dist_index's snarl_indexes/chain_index if it is a snarl/chain
        struct child_node_t {
            size_t node_id; 
            ChildNodeType node_type;

            id_t id_in_parent(MinimumDistanceIndex& dist_index) { 
                //Get the id of this node in the parent netgraph
                switch (node_type) {
                case NODE:
                    return node_id;
                case SNARL:
                    return dist_index.snarl_indexes[node_id].id_in_parent;
                case CHAIN:
                    return dist_index.chain_indexes[node_id].id_in_parent;
                default:
                    return 0;
                }
            }

            size_t rank_in_parent(MinimumDistanceIndex& dist_index, id_t id) {
                size_t rank = node_type == NODE ? 
                                    dist_index.getPrimaryRank(id) :
                                    dist_index.getSecondaryRank(id);
                if ( (node_type == SNARL && 
                      dist_index.snarl_indexes[dist_index.getPrimaryAssignment(id)].rev_in_parent) ||
                     (node_type == CHAIN && 
                      dist_index.chain_indexes[dist_index.getChainAssignment(id)].rev_in_parent)) {
                    rank = rank % 2 == 0 ? rank + 1 : rank - 1;
                }
                return rank;
            }
        };

//TODO: Put all these arguments into a struct
//One loop going up tree
//One method for each node/snarl/chain level

        //A cluster in the context of a snarl tree node
        //set of the seed indices in the cluster and left and right distance
        //from a seed in the cluster to the ends of the snarl tree node
        //containing the cluster
        typedef tuple<hash_set<size_t>, int64_t, int64_t> child_cluster_t;


        struct tree_state_t {
            //Hold all the tree relationships, seed locations, and cluster info
            //for the current level of the snarl tree and the parent level

            //Vector of all the seeds
            vector<pos_t> seeds; 

            //The minimum distance between nodes for them to be put in the
            //same cluster
            int64_t distance_limit;


            //////////Data structures to hold clustering information

            //Structure to hold the clustering of the seeds
            structures::UnionFind* union_find_clusters;

            //For each seed, store the distances to the left and right ends
            //of the netgraph node of the cluster it belongs to
            //These values are only relevant for seeds that represent a cluster
            //in union_find_clusters
            vector<pair<int64_t, int64_t>> cluster_dists;



            //////////Data structures to hold snarl tree relationships

            //Maps each node to a vector of the indices into seeds of the
            //seeds that are contained in it
            hash_map<id_t, vector<size_t>> node_to_seeds;

            //Map from snarl (index into dist_index.snarl_indexes) i
            //to the netgraph nodes contained in the snarl as well as the 
            //clusters at the node
            hash_map<size_t,vector<pair<child_node_t,child_cluster_t>>>
                                                            snarl_to_nodes;
            
            //Map each chain to the snarls (only ones that contain seeds) that
            //comprise it. 
            //Snarls are order by their order in the chain
            //Snarls and chains represented as their indexes into 
            //dist_index.chain/snarl_indexes
            hash_map<size_t, vector<size_t>> chain_to_snarls;


            //Same structure as snarl_to_nodes but for the level of the snarl
            //tree above the current one
            hash_map<size_t,vector<pair<child_node_t,child_cluster_t>>>
                                                          parent_snarl_to_nodes;


        };

        //Find which nodes contain seeds and assign those nodes to the 
        //snarls that contain them. Also find the depth of each snarl
        void get_nodes( const vector<pos_t>& seeds,
                        hash_map<id_t, vector<size_t>>& node_to_seeds,
                        vector<hash_map<size_t, 
                                  vector<pair<child_node_t, child_cluster_t>>>>&
                                                            snarl_to_nodes);

        //Given a node and the indices of seeds on that node, root, 
        //cluster the seeds
        child_cluster_t get_clusters_node(tree_state_t& tree_state, 
                                          id_t node_id, int64_t node_length); 

        //Cluster the seeds in a chain
        //snarls_in_chain is an unordered vector of snarls where the snarl
        //  is represented as its index into dist_index.snarl_indexes
        //curr_snarl_children maps each snarl to a vector of its children and 
        //clusters on the children
        child_cluster_t get_clusters_chain(tree_state_t& tree_state,
                                           size_t chain_index_i);

        //Cluster the seeds in a snarl 
        //child_nodes is a vector of the children of root and their clusters
        child_cluster_t get_clusters_snarl(tree_state_t& tree_state,
                             size_t snarl_index_i, bool rev) ;

};
}

#endif
