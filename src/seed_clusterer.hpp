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

        ///Given a vector of seeds (pos_t) and a distance limit, 
        //cluster the seeds such that two seeds whose minimum distance
        //between them (including both of the positions) is less than
        // the distance limit are in the same cluster
        //
        //Returns a vector of clusters. Each cluster is a vector of
        //indices into seeds
        vector<vector<size_t>> cluster_seeds ( 
                vector<pos_t> seeds, int64_t read_distance_limit) const;
        
        ///The same thing, but for paired end reads.
        //Given seeds from multiple reads of a fragment, cluster each set of seeds
        //by the read distance and all seeds by the fragment distance limit
        //fragment_distance_limit must be greater than read_distance_limit
        //Returns clusters for each read and clusters of all the seeds in all reads
        //The read clusters refer to seeds by their indexes in the input vectors of seeds
        //The fragment clusters give seeds the index they would get if the vectors of
        // seeds were appended to each other in the order given
        tuple<vector<vector<vector<size_t>>>,vector<vector<size_t>>> cluster_seeds ( 
                vector<vector<pos_t>> all_seeds,
                int64_t read_distance_limit, int64_t fragment_distance_limit=0) const;

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

        struct NetgraphNode {
            //child nodes of a snarl's netgraph 
            //node_id is the node id if the node is just a node, index into
            //dist_index's snarl_indexes/chain_index if it is a snarl/chain
            size_t node_id; 

            //Is this a node, snarl, or chain
            ChildNodeType node_type;

            NetgraphNode(size_t node_id, ChildNodeType node_type) :
                node_id(node_id), node_type(node_type) {
            }

            id_t id_in_parent(MinimumDistanceIndex& dist_index) const { 
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

            size_t rank_in_parent(MinimumDistanceIndex& dist_index, id_t id) const {
                //Get the forward rank of this node in the parent's netgraph
                //to look up distances

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

        struct NodeClusters {
            //Clusters in the context of a snarl tree node
            //The node containing this struct may be an actual node,
            // snarl/chain that is a node the parent snarl's netgraph,
            // or a snarl in a chain

            //set of the indices of heads of clusters (group ids in the 
            //union find)
            hash_set<size_t> cluster_heads;

            //The shortest distance from any seed in any cluster to the 
            //left/right end of the snarl tree node that contains these
            //clusters
            int64_t best_left;
            int64_t best_right;

            NodeClusters() :
                best_left(-1), best_right(-1) {}
        };


        struct TreeState {
            //Hold all the tree relationships, seed locations, and cluster info
            //for the current level of the snarl tree and the parent level
            //As clustering occurs at the current level, the parent level
            //is updated to know about its children

            //Vector of all the seeds for each read
            vector<vector<pos_t>>* all_seeds; 

            //Vector of the offset of indices for each seed
            vector<size_t> seed_index_offsets;

            //The minimum distance between nodes for them to be put in the
            //same cluster
            int64_t read_distance_limit;
            int64_t fragment_distance_limit;


            //////////Data structures to hold clustering information

            //Structure to hold the clustering of the seeds
            vector<structures::UnionFind> read_union_find;
            structures::UnionFind fragment_union_find;

            //For each seed, store the distances to the left and right ends
            //of the netgraph node of the cluster it belongs to
            //These values are only relevant for seeds that represent a cluster
            //in union_find_reads
            vector<pair<int64_t, int64_t>> read_cluster_dists;



            //////////Data structures to hold snarl tree relationships
            //The snarls and chains get updated as we move up the snarl tree

            //Maps each node to a vector of the seeds that are contained in it
            //seeds are represented by indexes into the seeds vector
            //The array is sorted.
            vector<tuple<id_t, size_t, size_t>> node_to_seeds;

            //Map from snarl (index into dist_index.snarl_indexes) i
            //to the netgraph nodes contained in the snarl as well as the 
            //clusters at the node
            hash_map<size_t,vector<pair<NetgraphNode, NodeClusters>>>
                                                            snarl_to_nodes;
            
            //Map each chain to the snarls (only ones that contain seeds) that
            //comprise it. 
            //Snarls and chains represented as their indexes into 
            //dist_index.chain/snarl_indexes
            //Map maps the rank of the snarl to the snarl and snarl's clusters
            //  Since maps are ordered, it will be in the order of traversal
            //  of the snarls in the chain
            hash_map<size_t, std::map<size_t, pair<size_t, NodeClusters>>>
                                                         chain_to_snarls;


            //Same structure as snarl_to_nodes but for the level of the snarl
            //tree above the current one
            //This gets updated as the current level is processed
            hash_map<size_t,vector<pair<NetgraphNode,NodeClusters>>>
                                                          parent_snarl_to_nodes;

            //Constructor takes in a pointer to the seeds and the distance limit 
            TreeState (vector<vector<pos_t>>* all_seeds, int64_t read_distance_limit, 
                       int64_t fragment_distance_limit) :
                all_seeds(all_seeds),
                read_cluster_dists(seeds->size(), make_pair(-1, -1)),
                read_union_find (seeds->size(), false),
                fragment_union_find (seeds->size(), false),
                read_distance_limit(read_distance_limit),
                fragment_distance_limit(fragment_distance_limit){
                    seed_index_offsets.push_back(0);
                    for (auto& v : all_seeds) {
                        size_t offset = seed_index_offsets.back() + v.size();
                        seed_index_offsets.push_back(offset);
                    }
            }
        };

        //Find which nodes contain seeds and assign those nodes to the 
        //snarls that contain them
        //Update the tree state's node_to_seed
        //and snarl_to_nodes_by_level, which assigns each node that contains
        //seeds to a snarl, organized by the level of the snarl in the snarl 
        //tree. snarl_to_nodes_by_level will be used to populate snarl_to_nodes
        //in the tree state as each level is processed
        void get_nodes( TreeState& tree_state,
                        vector<hash_map<size_t, 
                                  vector<pair<NetgraphNode, NodeClusters>>>>&
                                                       snarl_to_nodes_by_level) const;

        //Cluster all the snarls at the current level and update the tree_state
        //to add each of the snarls to the parent level
        void cluster_snarls(TreeState& tree_state, size_t depth) const;

        //Cluster all the chains at the current level
        void cluster_chains(TreeState& tree_state, size_t depth) const;

        //Given a node and the indices of seeds on that node, root, 
        //cluster the seeds
        NodeClusters cluster_one_node(TreeState& tree_state, 
                                          id_t node_id, int64_t node_length) const; 

        //Cluster the seeds in a chain given by chain_index_i, an index into
        //dist_index.chain_indexes
        NodeClusters cluster_one_chain(TreeState& tree_state,
                                           size_t chain_index_i) const;

        //Cluster the seeds in a snarl given by snarl_index_i, an index into
        //dist_index.snarl_indexes
        //rev is true if this snarl is reversed in its parent
        NodeClusters cluster_one_snarl(TreeState& tree_state,
                                       size_t snarl_index_i) const;

};
}

#endif
