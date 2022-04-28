#ifndef VG_SEED_CLUSTERER_HPP_INCLUDED
#define VG_SEED_CLUSTERER_HPP_INCLUDED

#include "snarls.hpp"
#include "snarl_distance_index.hpp"
#include "hash_map.hpp"
#include "small_bitset.hpp"
#include <structures/union_find.hpp>


namespace vg{

class NewSnarlSeedClusterer {



    public:

        /// Seed information used in Giraffe.
        struct Seed {
            pos_t  pos;
            size_t source; // Source minimizer.


            //Cached values from the minimizer
            //node length, root component, prefix sum, chain component, is_reversed
            tuple<size_t, size_t, size_t, size_t, bool> minimizer_cache  = 
                make_tuple(MIPayload::NO_VALUE, MIPayload::NO_VALUE, MIPayload::NO_VALUE, MIPayload::NO_VALUE, false);

            //The distances to the left and right of whichever cluster this seed represents
            //This gets updated as clustering proceeds
            size_t distance_left = std::numeric_limits<size_t>::max();
            size_t distance_right = std::numeric_limits<size_t>::max();

            net_handle_t node_handle;

        };

        /// Cluster information used in Giraffe.
        struct Cluster {
            std::vector<size_t> seeds; // Seed ids.
            size_t fragment; // Fragment id.
            double score; // Sum of scores of distinct source minimizers of the seeds.
            double coverage; // Fraction of read covered by the seeds.
            SmallBitset present; // Minimizers that are present in the cluster.
        };

        NewSnarlSeedClusterer(const SnarlDistanceIndex& distance_index, const HandleGraph* graph);

        //TODO: I don't want to be too tied to the minimizer_mapper implementation with seed structs

        ///Given a vector of seeds and a distance limit, 
        //cluster the seeds such that two seeds whose minimum distance
        //between them (including both of the positions) is less than
        // the distance limit are in the same cluster

        vector<Cluster> cluster_seeds ( vector<Seed>& seeds, size_t read_distance_limit) const;
        
        ///The same thing, but for paired end reads.
        //Given seeds from multiple reads of a fragment, cluster each read
        //by the read distance and all seeds by the fragment distance limit.
        //fragment_distance_limit must be greater than read_distance_limit
        //Returns clusters for each read and clusters of all the seeds in all reads
        //The read clusters refer to seeds by their indexes in the input vectors of seeds
        //The fragment clusters give seeds the index they would get if the vectors of
        // seeds were appended to each other in the order given
        // TODO: Fix documentation
        // Returns: For each read, a vector of clusters.

        vector<vector<Cluster>> cluster_seeds ( 
                vector<vector<Seed>>& all_seeds, size_t read_distance_limit, size_t fragment_distance_limit=0) const;

    private:


        //Actual clustering function that takes a vector of pointers to seeds
        tuple<vector<structures::UnionFind>, structures::UnionFind> cluster_seeds_internal ( 
                vector<vector<Seed>*>& all_seeds,
                size_t read_distance_limit, size_t fragment_distance_limit=0) const;

        const SnarlDistanceIndex& distance_index;
        const HandleGraph* graph;

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


        struct NodeClusters {
            //All clusters of a snarl tree node
            //The node containing this struct may be an actual node,
            // snarl/chain that is a node the parent snarl's netgraph,
            // or a snarl in a chain

            //The snarl tree node that the clusters are on
            net_handle_t containing_net_handle; 

            //Only set these for nodes or snarls in chains
            nid_t node_id = 0;
            bool is_reversed_in_parent = false;

            //Minimum length of a node or snarl
            //If it is a chain, then it is distance_index.chain_minimum_length(), which is
            //the expected length for a normal chain, and the length of the 
            //last component for a multicomponent chain 
            size_t node_length = std::numeric_limits<size_t>::max();             size_t prefix_sum_value = std::numeric_limits<size_t>::max(); //of node or first node in snarl
            size_t chain_component_start = 0; //of node or start of snarl
            size_t chain_component_end = 0; //of node or end of snarl

            //Only set this one for a chain
            bool is_looping_chain = false;
            size_t chain_last_component = std::numeric_limits<size_t>::max();
            size_t chain_last_child_offset = std::numeric_limits<size_t>::max();

            //This one gets set for a (nontrivial) chain or snarl
            net_handle_t start_in;
            net_handle_t end_in;
            

            //set of the indices of heads of clusters (group ids in the 
            //union find)
            //TODO: Add cluster distances here
            //maps pair of <read index, seed index> to pair of <left distance, right distance>
            hash_set<pair<size_t, size_t>> read_cluster_heads;

            //The shortest distance from any seed in any cluster to the 
            //left/right end of the snarl tree node that contains these
            //clusters
            size_t fragment_best_left = std::numeric_limits<size_t>::max();
            size_t fragment_best_right = std::numeric_limits<size_t>::max();
            vector<size_t> read_best_left;
            vector<size_t> read_best_right;

            //Distance from the start of the parent to the left of this node, etc
            size_t distance_start_left = std::numeric_limits<size_t>::max();
            size_t distance_start_right = std::numeric_limits<size_t>::max();
            size_t distance_end_left = std::numeric_limits<size_t>::max();
            size_t distance_end_right = std::numeric_limits<size_t>::max();

            //Constructor
            //read_count is the number of reads in a fragment (2 for paired end)
            NodeClusters( net_handle_t net, size_t read_count, size_t seed_count, const SnarlDistanceIndex& distance_index) :
                containing_net_handle(std::move(net)),
                fragment_best_left(std::numeric_limits<size_t>::max()), fragment_best_right(std::numeric_limits<size_t>::max()),
                read_best_left(read_count, std::numeric_limits<size_t>::max()), 
                read_best_right(read_count, std::numeric_limits<size_t>::max()){
                read_cluster_heads.reserve(seed_count);
                if (distance_index.is_chain(containing_net_handle) && !distance_index.is_trivial_chain(containing_net_handle)) {
                    is_looping_chain = distance_index.is_looping_chain(containing_net_handle);
                    node_length = distance_index.chain_minimum_length(containing_net_handle);
                    start_in = distance_index.get_bound(containing_net_handle, false, true);
                    end_in = distance_index.get_bound(containing_net_handle, true, true);
                    chain_last_component = distance_index.get_chain_component(end_in, true);
                } else if (distance_index.is_snarl(containing_net_handle)) {
                    node_length = distance_index.minimum_length(containing_net_handle);
                    start_in = distance_index.get_node_from_sentinel(distance_index.get_bound(containing_net_handle, false, true));
                    end_in =   distance_index.get_node_from_sentinel(distance_index.get_bound(containing_net_handle, true, true));
                }
            }
            NodeClusters( net_handle_t net, size_t read_count, size_t seed_count, bool is_reversed_in_parent, nid_t node_id, size_t node_length, size_t prefix_sum, size_t component) :
                containing_net_handle(net),
                is_reversed_in_parent(is_reversed_in_parent),
                node_length(node_length),
                prefix_sum_value(prefix_sum),
                chain_component_start(component),
                chain_component_end(component),
                node_id(node_id),
                fragment_best_left(std::numeric_limits<size_t>::max()), fragment_best_right(std::numeric_limits<size_t>::max()),
                read_best_left(read_count, std::numeric_limits<size_t>::max()), 
                read_best_right(read_count, std::numeric_limits<size_t>::max()){
                    read_cluster_heads.reserve(seed_count);
                }

        };

        struct ParentToChildMap {
            //Struct for storing a map from a parent net_handle_t to a list of it's children
            //The children are represented as an index into all_clusters

            //The actual data that gets stored
            //The first size_t is the parent, as an index into all_clusters
            //The bool is true if the child is a snarl and false if it is a single seed on a node
            //The second size_t is for the child; If it is a snarl, then an index into all_seeds
            // if it is a seed, then the two indies into all_seeds
            //This stores every child as a separate pair
            vector<tuple<size_t, bool, size_t, size_t>> parent_to_children;

            //is parent_to_children sorted?
            //Each time we look up the children of a parent, sort and look it up
            //If it's already sorted, skip sorting
            //Also set this to false any time something gets added
            bool is_sorted = false;
            bool is_sorted_children = false;

            void add_child(size_t parent_index, bool is_snarl, size_t child_index, size_t child_index2) {
                parent_to_children.emplace_back(parent_index, is_snarl, child_index, child_index2);
                is_sorted=false;
            }
            void reserve(size_t size) {
                parent_to_children.reserve(size);
            }

            //Sort the parent_to_children vector first by parent, and second by the order
            //of the children determined by comparator
            void sort(const std::function<bool(const tuple<size_t, bool, size_t, size_t>&,
                                               const tuple<size_t, bool, size_t, size_t>&)>& comparator) {
                if (!is_sorted) {
                    std::sort(parent_to_children.begin(), parent_to_children.end(),
                    [&] (const tuple<size_t, bool, size_t, size_t>& a,
                         const tuple<size_t, bool, size_t, size_t>& b)->bool {
                        if (std::get<0>(a) == std::get<0>(b)) {
                            return comparator(a, b);
                        } else {
                            return std::get<0>(a) < std::get<0>(b);
                        }
                    });
                    is_sorted = true;
                }
            }

            //Get a list of the children of this parent
            //Equivalent of map[parent]
            //Does this by sorting (if necessary) the vector parent_to_children of the parent
            //and then finding the first occurrence of the parent using std::lower_bound and walking
            //through the vector
            //The vector of children will not be sorted
            vector<tuple<bool, size_t, size_t>> get_children(const size_t& parent,
                    const std::function<bool(const tuple<size_t, bool, size_t, size_t>&,
                                             const tuple<size_t, bool, size_t, size_t>&)>& comparator) {
                //We need to sort the vector first to find everything with the right parent
                if (!is_sorted) {
                    sort(comparator);
                }
                vector<tuple<bool, size_t, size_t>> children;
                auto iter_start = std::lower_bound(parent_to_children.begin(), parent_to_children.end(),
                        std::tuple<size_t, bool, size_t, size_t>(parent, (size_t)0, (size_t)0, (size_t)0));
                for (auto iter = iter_start ; iter != parent_to_children.end() && std::get<0>(*iter) == parent ; ++iter) {
                    children.emplace_back(std::get<1>(*iter), std::get<2>(*iter), std::get<3>(*iter));
                }
                return children;
            }
        };

        struct TreeState {
            //Hold all the tree relationships, seed locations, and cluster info

            //for the current level of the snarl tree and the parent level
            //As clustering occurs at the current level, the parent level
            //is updated to know about its children

            //Vector of all the seeds for each read
            vector<vector<Seed>*>* all_seeds; 

            //prefix sum vector of the number of seeds per read
            //To get the index of a seed for the fragment clusters
            //Also se this so that data structures that store information per seed can be single
            //vectors, instead of a vector of vectors following the structure of all_seeds 
            //since it uses less memory allocation to use a single vector
            vector<size_t> seed_count_prefix_sum;

            //The minimum distance between nodes for them to be put in the
            //same cluster
            size_t read_distance_limit;
            size_t fragment_distance_limit;


            //////////Data structures to hold clustering information

            //Structure to hold the clustering of the seeds
            vector<structures::UnionFind> read_union_find;
            structures::UnionFind fragment_union_find;



            //////////Data structures to hold snarl tree relationships
            //The snarls and chains get updated as we move up the snarl tree

            //Maps each node to a vector of the seeds that are contained in it
            //seeds are represented by indexes into the seeds vector (read_num, seed_num)
            //The array is sorted.
            vector<tuple<id_t, size_t, size_t>> node_to_seeds;

            //This stores all the node clusters so we stop spending all our time allocating lots of vectors of NodeClusters
            vector<NodeClusters> all_node_clusters;

            //Map each net_handle to its index in all_node_clusters
            hash_map<net_handle_t, size_t> net_handle_to_index;

            //Map from each snarl's net_handle_t to it's child net_handle_ts 
            //clusters at the node
            //size_t is the index into all_node_clusters
            std::unordered_multimap<size_t, size_t> snarl_to_children;
            
            //Map each chain to the snarls (only ones that contain seeds) that
            //comprise it. 
            //Snarls and chains represented as their indexes into 
            //distance_index.chain/snarl_indexes
            //Map maps the rank of the snarl to the snarl and snarl's clusters
            //  Since maps are ordered, it will be in the order of traversal
            //  of the snarls in the chain
            //  size_t is the index into all_node_clusters
            ParentToChildMap* chain_to_children;


            //Same structure as chain_to_children but for the level of the snarl
            //tree above the current one
            //This gets updated as the current level is processed
            //size_t is the index into all_node_clusters
            ParentToChildMap* parent_chain_to_children;

            //This holds all the child clusters of the root
            //each size_t is the index into all_node_clusters
            //Each pair is the parent and the child. This will be sorted by parent before
            //clustering so it
            vector<pair<size_t, size_t>> root_children;


            /////////////////////////////////////////////////////////

            //Constructor takes in a pointer to the seeds, the distance limits, and 
            //the total number of seeds in all_seeds
            TreeState (vector<vector<Seed>*>* all_seeds, size_t read_distance_limit, 
                       size_t fragment_distance_limit, size_t seed_count) :
                all_seeds(all_seeds),
                read_distance_limit(read_distance_limit),
                fragment_distance_limit(fragment_distance_limit),
                fragment_union_find (seed_count, false),
                seed_count_prefix_sum(1,0){

                for (size_t i = 0 ; i < all_seeds->size() ; i++) {
                    size_t size = all_seeds->at(i)->size();
                    size_t offset = seed_count_prefix_sum.back() + size;
                    seed_count_prefix_sum.push_back(offset);
                    read_union_find.emplace_back(size, false);

                }

                all_node_clusters.reserve(5*seed_count);
                net_handle_to_index.reserve(5*seed_count);
                root_children.reserve(seed_count);
            }
        };

        //Find which nodes contain seeds and assign those nodes to the 
        //snarls that contain them
        //Update the tree state's node_to_seed
        //and snarl_to_nodes_by_level, which assigns each node that contains
        //seeds to a snarl, organized by the level of the snarl in the snarl 
        //tree. snarl_to_nodes_by_level will be used to populate snarl_to_nodes
        //in the tree state as each level is processed
        //size_t is the index into all_node_clusters
        void get_nodes( TreeState& tree_state,
                        vector<ParentToChildMap>& chain_to_children_by_level) const;


        //Cluster all the snarls at the current level and update the tree_state
        //to add each of the snarls to the parent level
        //void cluster_snarl_level(TreeState& tree_state) const;

        //Cluster all the chains at the current level
        void cluster_chain_level(TreeState& tree_state) const;

        //Cluster the seeds on the specified node
        void cluster_one_node(TreeState& tree_state, NodeClusters& node_clusters) const; 

        //Cluster the seeds in a snarl given by its net handle
        //Snarl_cluster_index is the index into tree_state.all_node_clusters
        void cluster_one_snarl(TreeState& tree_state, size_t snarl_clusters_index, size_t parent_clusters_index) const;

        //Cluster the seeds in a chain given by chain_index_i, an index into
        //distance_index.chain_indexes
        //If the depth is 0, also incorporate the top-level seeds from tree_state.top_level_seed_clusters
        //Chain children are tuples<is_snarl, (child index, inf) or (seed read num, seed index)>
        //If the children of the chain are only seeds on nodes, then cluster as if it is a node
        void cluster_one_chain(TreeState& tree_state, size_t chain_clusters_index, vector<tuple<bool, size_t, size_t>>& children_in_chain, bool only_seeds) const;

        //Cluster in the root 
        void cluster_root(TreeState& tree_state) const;

        //Cluster a list of seeds (SeedIndexes) that are on a single linear structure (node or chain)
        //Requires that the list of seeds are sorted relative to their position on the structure
        //The list of seeds is everything in the list between range_start and range_end
        //This can be called on a chain if there are no nested seeds on the chain
        //get_offset_from_seed_index returns a tuple of <read_num, seed_num, left offset> indices into all_seeds from whatever
        //SeedIndex is used to store the seeds
        //left offset is the distance from the left side of the structure
        template <typename SeedIndex>
        void cluster_seeds_on_linear_structure(TreeState& tree_state, NodeClusters& node_clusters, vector<SeedIndex>& seed_indices, 
                size_t structure_length, std::function<std::tuple<size_t, size_t, size_t>(const SeedIndex&)>& get_offset_from_seed_index) const;

        //Compare two children of the parent and combine their clusters, to create clusters in the parent
        //This assumes that the first node hasn't been seen before but the second one has, so all of the
        //first node's clusters get added to the parent but assume that all of the second ones are already
        //part of the parent
        //old_distances contains the distances for cluster heads in the children, 
        //since the distances in tree_state.read_cluster_heads_to_distances will get updated
        void compare_and_combine_cluster_on_child_structures(TreeState& tree_state, NodeClusters& child_clusters1,
                NodeClusters& child_clusters2, NodeClusters& parent_clusters, 
                const vector<pair<size_t, size_t>>& child_distances, bool is_root, bool get_distances_to_parent) const;

        //The same as above, but compare clusters on a single child
        //This assumes that the child is the child of the root and not a root snarl
        //so we just look at external distances 
        void compare_and_combine_cluster_on_one_child(TreeState& tree_state, NodeClusters& child_clusters) const;

};
}

#endif
