#ifndef VG_SNARL_SEED_CLUSTERER_HPP_INCLUDED
#define VG_SNARL_SEED_CLUSTERER_HPP_INCLUDED

#include "snarls.hpp"
#include "snarl_distance_index.hpp"
#include "zip_code.hpp"
#include "hash_map.hpp"
#include "small_bitset.hpp"
#include <structures/union_find.hpp>


namespace vg{


/**
 * SnarlDistanceIndexClusterer is used for clustering seeds (positions on the graph)
 * A "cluster" is a partition of seeds that is based on the minimum distance between them in the graph
 * Consider a graph where each seed is a node and two seeds are connected if the minimum distance
 * between them is smaller than a given distance limit. Each connected component of this graph is a cluster 
 *
 * The clustering algorithm is based on the snarl tree
 * Clusters are formed on nodes of the snarl tree, which represent nodes/snarls/chains
 * Each node/snarl/chain represents a subgraph of the variation graph
 * A clustered snarl tree node contains all seeds that occur on its subgraph, and the seeds have been partitioned into clusters
 * Each cluster knows the shortest distance from any seed it contains to both ends of the snarl tree node containing it
 * Clustering is done progressively by walking up the snarl tree and forming clusters on each snarl tree node (only visiting nodes that have seeds on them)
 * At each snarl tree node, assume that its children have already been clustered. 
 * The clusters of the children are compared to each other, and any pair that are close enough 
 * are combined to produce clusters on the parent 
 * The distances from each cluster to the ends of the parent are updated
 *
 * The algorithm starts by assigning each seed to its node on the snarl tree
 * Since nodes are all on chains, this fills in all the children of chains that are nodes
 * It then walks up the snarl tree, level by level, and clusters each snarl tree node that contains seeds
 * At a given level, first cluster each chain in the level. After clustering a chain, assign it 
 * to its parent snarl. Then, go through each of the snarls that have just been given children, and
 * cluster the snarls. Each snarl then gets assigned to its parent chain
 * This completes one level of the snarl tree. Each chain in the next level has just been populated by the snarls
 * from this level, and already knew about its nodes from the first step, so it is ready to be clustered 
 *
 * Every time the clusterer is run, a ClusteringProblem is made to store information about the state of the clusterer
 * The ClusteringProblem keeps track of which level of the snarl tree is currently being clustered, and
 * keeps track of the children of the current and next level of the snarl tree. 
 * Each snarl tree node that contains seeds is represented by a SnarlTreeNodeProblem.
 * The SnarlTreeNodeProblem represents the problem of clustering one snarl tree node. 
 * It knows the identities of its children and keeps track of its cluster heads
 * 
 * 
 *
 */
class SnarlDistanceIndexClusterer {



    public:

        /// Seed information used in Giraffe.
        struct Seed {
            /// Position of the seed.
            ///
            /// If the minimizer is from the read sequence's forward strand,
            /// this corresponds to the first base in the read that is part of
            /// the minimizer occurrence, and points in the read's forward
            /// direction.
            ///
            /// If the minimizer is from the read sequence's reverse strand,
            /// this corresponds to the *last* base in the read that is part of
            /// the minimizer occurrence, but *still* points in the read's
            /// *forward* direction.
            pos_t  pos;
            size_t source; // Source minimizer.

            //zipcode for distance information, optionally stored in the minimizer payload
            // Clustering requires that the zipcode and its decoder are filled in
            ZipCode zipcode; 

            Seed() = default;
            Seed(pos_t pos, size_t source, ZipCode zipcode) : pos(pos), source(source), zipcode(zipcode) {
                zipcode.fill_in_full_decoder();
            }

            //Move constructor
            Seed (Seed&& other) :
                pos(std::move(other.pos)),
                source(std::move(other.source)),
                zipcode(std::move(other.zipcode)){}

            //Move assignment operator
            Seed& operator=(Seed&& other) {
                pos = std::move(other.pos);
                source = std::move(other.source);
                zipcode = std::move(other.zipcode);
                return *this;
            }
        };

        /// Seed information used for clustering
        // Corresponds to one seed and stores the minimizer payload and distance information 
        // that gets updated during clustering
        // TODO: This will copy information from the seed, since we need per-seed information anyways
        // and some of it needs to be mutable, it's simpler than keeping around two collections of Seeds
        struct SeedCache{
            const Seed* seed;

            //TODO: I think I can skip the zipcode now since I have the payload
            MIPayload payload;

            //The distances to the left and right of whichever cluster this seed represents
            //This gets updated as clustering proceeds
            //For a seed in a chain, distance_left is the left of the chain, right is the distance
            //to the right side of the node, relative to the chain
            size_t distance_left = std::numeric_limits<size_t>::max();
            size_t distance_right = std::numeric_limits<size_t>::max();
            //Values from the payload that we're saving
            size_t payload_prefix_sum = std::numeric_limits<size_t>::max();
            size_t payload_node_length = std::numeric_limits<size_t>::max();

        };

        /// Cluster information used in Giraffe.
        struct Cluster {
            std::vector<size_t> seeds; // Seed ids.
            size_t fragment; // Fragment id.
            double score; // Sum of scores of distinct source minimizers of the seeds.
            double coverage; // Fraction of read covered by the seeds.
            SmallBitset present; // Minimizers that are present in the cluster.
        };

        SnarlDistanceIndexClusterer(const SnarlDistanceIndex& distance_index, const HandleGraph* graph);
        SnarlDistanceIndexClusterer(const SnarlDistanceIndex* distance_index, const HandleGraph* graph);
        SnarlDistanceIndexClusterer(const SnarlDistanceIndex& distance_index);
        SnarlDistanceIndexClusterer(const SnarlDistanceIndex* distance_index);


        /*Given a vector of seeds and a distance limit, 
         *cluster the seeds such that two seeds whose minimum distance
         *between them (including both of the positions) is less than
         *the distance limit are in the same cluster
         *This produces a vector of clusters
         *
         * Requires that the zipcodes and decoders in the seeds are filled in
         */
        vector<Cluster> cluster_seeds ( const vector<Seed>& seeds, size_t read_distance_limit) const;
        
        /* The same thing, but for paired end reads.
         * Given seeds from multiple reads of a fragment, cluster each read
         * by the read distance and all seeds by the fragment distance limit.
         * fragment_distance_limit must be greater than read_distance_limit
         * Returns a vector clusters for each read, where each cluster also has an assignment
         * to a fragment cluster
         *
         * Requires that there are only two reads per fragment (all_seeds.size() == 2, meaning paired end reads)
         *    this requirement is just because I used std::pairs to represent two reads, but could be changed to a vector if we every have to map more than two reads per fragment
         *
         * Requires that the zipcodes and decoders in the seeds are filled in
         */

        vector<vector<Cluster>> cluster_seeds ( 
                const vector<vector<Seed>>& all_seeds, 
                size_t read_distance_limit, size_t fragment_distance_limit=0) const;


    private:


        //Actual clustering function that takes a vector of pointers to seeds
        //fragment_distance_limit defaults to 0, meaning that we don't cluster by fragment
        tuple<vector<structures::UnionFind>, structures::UnionFind> cluster_seeds_internal ( 
                vector<vector<SeedCache>*>& all_seeds,
                size_t read_distance_limit, size_t fragment_distance_limit=0) const;

        const SnarlDistanceIndex& distance_index;
        const HandleGraph* graph;


        /*
         * This struct is used to store the clustering information about one 
         * snarl tree node (node/snarl/chain)
         *
         * It knows the cluster heads of the clusters on the node 
         * and the minimum distance from any seed in each cluster to the ends of the node
         * If the node is a snarl, then the distances stored are to the boundary nodes but
         * don't include the lengths of the boundary nodes; if the node is a node or chain,
         * then the distances include the boundary nodes
         *
         * Relevant children of the snarl tree node are stored as SnarlTreeChild's, which 
         * may represent a seed (if the parent is a chain) or another snarl tree node
         * The list of children is unsorted and must be sorted before clustering a chain
         *
         * This also stores additional information about the snarl tree node from the distance index
         * including the distance from the ends of the node to the ends of the parent
         */
        struct SnarlTreeNodeProblem {

            //set of the indices of heads of clusters (group ids in the 
            //union find)
            //pair of <read index, seed index>
            hash_set<pair<size_t, size_t>> read_cluster_heads;

            //Struct to store one child, which may be a seed, node, snarl, or chain
            struct SnarlTreeChild {
                //If the net_handle is a node, then the child is a seed, otherwise the handle 
                //is used to find the problem
                net_handle_t net_handle;
                pair<size_t, size_t> seed_indices;

                //The values used to sort the children of a chain
                //Storing it here is faster than looking it up each time
                size_t chain_component;
                size_t prefix_sum;
                //Is this child a seed
                //This is redundant with net_handle because any net_handle_t that is a node will really be a seed,
                //but it's faster than looking it up in the distance index
                bool is_seed;
                //Have chain_component and prefix_sum been set?
                //For a seed, it gets set when the child is made, otherwise the first time this 
                //child is seen when sorting
                bool has_chain_values;
            };
            //The children of this snarl tree node
            //Initially unsorted, sort before clustering for chains
            vector<SnarlTreeChild> children;

            //The shortest distance from any seed in any cluster to the 
            //left/right end of the snarl tree node that contains these
            //clusters
            pair<size_t, size_t> read_best_left = make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
            pair<size_t, size_t> read_best_right = make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
            size_t fragment_best_left = std::numeric_limits<size_t>::max();
            size_t fragment_best_right = std::numeric_limits<size_t>::max();

            //Distance from the start of the parent to the left of this node, etc
            size_t distance_start_left = std::numeric_limits<size_t>::max();
            size_t distance_start_right = std::numeric_limits<size_t>::max();
            size_t distance_end_left = std::numeric_limits<size_t>::max();
            size_t distance_end_right = std::numeric_limits<size_t>::max();

            //The snarl tree node that the clusters are on
            net_handle_t containing_net_handle; 




            //The parent and grandparent of containing_net_handle, which might or might not be set
            //This is just to store information from the minimizer cache
            net_handle_t parent_net_handle;
            net_handle_t grandparent_net_handle;

            //One representative seed so we can get the zipcode and stuff
            const SeedCache* seed;
            size_t zipcode_depth;

            //Minimum length of a node or snarl
            //If it is a chain, then it is distance_index.chain_minimum_length(), which is
            //the expected length for a normal chain, and the length of the 
            //last component for a multicomponent chain 
            size_t node_length = std::numeric_limits<size_t>::max();             
            size_t prefix_sum_value = std::numeric_limits<size_t>::max(); //of node or first node in snarl
            size_t chain_component_start = 0; //of node or start of snarl
            size_t chain_component_end = 0; //of node or end of snarl

            size_t loop_left = std::numeric_limits<size_t>::max();
            size_t loop_right = std::numeric_limits<size_t>::max();

            //These are sometimes set if the value was in the cache
            bool has_parent_handle = false;
            bool has_grandparent_handle = false;

            //Only set this for nodes or snarls in chains
            bool is_reversed_in_parent = false;

            bool is_trivial_chain = false;
            bool is_looping_chain = false;
            



            //Constructor
            //read_count is the number of reads in a fragment (2 for paired end)
            SnarlTreeNodeProblem( net_handle_t net, size_t read_count, size_t seed_count, const SnarlDistanceIndex& distance_index, 
                                  const SeedCache* seed, size_t zipcode_depth) :
                containing_net_handle(std::move(net)),
                fragment_best_left(std::numeric_limits<size_t>::max()), fragment_best_right(std::numeric_limits<size_t>::max()),
                seed(seed),
                zipcode_depth(zipcode_depth) {
                read_cluster_heads.reserve(seed_count);
            }
            //Constructor for a node or trivial chain, used to remember information from the cache
            SnarlTreeNodeProblem( net_handle_t net, size_t read_count, size_t seed_count, bool is_reversed_in_parent, 
                                 size_t node_length, size_t prefix_sum, size_t component, const SeedCache* seed, size_t zipcode_depth) :
                containing_net_handle(net),
                is_reversed_in_parent(is_reversed_in_parent),
                node_length(node_length),
                prefix_sum_value(prefix_sum),
                chain_component_start(component),
                chain_component_end(component),
                fragment_best_left(std::numeric_limits<size_t>::max()), fragment_best_right(std::numeric_limits<size_t>::max()),
                seed(seed),
                zipcode_depth(zipcode_depth) {
                    read_cluster_heads.reserve(seed_count);
            }

            //Set the values needed to cluster a chain
            void set_chain_values(const SnarlDistanceIndex& distance_index) {
                ZipCode::chain_code_t chain_code = seed->seed->zipcode.unpack_chain_code(zipcode_depth);
                is_looping_chain = chain_code.get_is_looping_chain();
                node_length = zipcode_depth == 0 ? distance_index.chain_minimum_length(containing_net_handle)
                                                 : chain_code.get_length();
                chain_component_end = chain_code.get_last_component();
                is_reversed_in_parent = seed->seed->zipcode.get_is_reversed_in_parent(zipcode_depth);
            }

            //Set the values needed to cluster a snarl
            void set_snarl_values(const SnarlDistanceIndex& distance_index) {
                ZipCode::snarl_code_t snarl_code = seed->seed->zipcode.unpack_snarl_code(zipcode_depth);
                node_length = snarl_code.get_length();
                chain_component_start = snarl_code.get_chain_component();
                chain_component_end = node_length == std::numeric_limits<size_t>::max() ? chain_component_start+1
                                                                                      : chain_component_start;
                prefix_sum_value = snarl_code.get_prefix_sum_or_identifier();

                net_handle_t start_in = distance_index.get_node_from_sentinel(distance_index.get_bound(containing_net_handle, false, true));
                net_handle_t end_in = distance_index.get_node_from_sentinel(distance_index.get_bound(containing_net_handle, true, true));
                loop_right = SnarlDistanceIndex::sum(distance_index.get_forward_loop_value(end_in),
                                                             2*distance_index.minimum_length(end_in));
                //Distance to go backward in the chain and back
                loop_left = SnarlDistanceIndex::sum(distance_index.get_reverse_loop_value(start_in),
                                                            2*distance_index.minimum_length(start_in));


            }

        };

        //These will be the cluster heads and distances for a cluster
        struct ClusterHead {
            size_t read_num = std::numeric_limits<size_t>::max();
            size_t cluster_num = 0;
            size_t distance_left = 0;
            size_t distance_right = 0;

            inline ClusterHead() {}
            inline ClusterHead(const size_t& read_num, const size_t& cluster_num, 
                           const size_t& distance_left, const size_t& distance_right) :
                read_num(read_num), cluster_num(cluster_num), 
                distance_left(distance_left), distance_right(distance_right) {} 
        };


        /* Hold all the tree relationships, seed locations, and cluster info
         * for the current level of the snarl tree and the parent level
         * As clustering occurs at the current level, the parent level
         * is updated to know about its children
         *
         * One "level" is the chains at that level, and their parent snarls.
         * Clustering one level means clustering the chains and then clustering the 
         * parent snarls. The parent snarls then get assigned to their parent chains,
         * and ClusteringProblem gets reset for the next level (parent chains)
         */
        struct ClusteringProblem {

            //Vector of all the seeds for each read
            vector<vector<SeedCache>*>* all_seeds; 

            //prefix sum vector of the number of seeds per read
            //Used to get the index of a seed for the fragment clusters
            //Also use this so that data structures that store information per seed can be single
            //vectors, instead of a vector of vectors following the structure of all_seeds 
            //since it uses less memory allocation to use a single vector
            vector<size_t> seed_count_prefix_sum;

            //The distance limits.
            //If the minimum distance between two seeds is less than this, 
            //they get put in the same cluster
            size_t read_distance_limit;
            size_t fragment_distance_limit;


            //////////Data structures to hold clustering information

            //Structure to hold the clustering of the seeds
            vector<structures::UnionFind> read_union_find;
            //The indices of seeds in the union find are the indices if you appended each of
            //the vectors of seeds for the fragment (i.e. if a seed in the second read is
            //at index x in the second vector of seeds, then its index in fragment_union_find
            //is x + the length of the first vector of seeds)
            structures::UnionFind fragment_union_find;



            //////////Data structures to hold snarl tree relationships
            //The snarls and chains get updated as we move up the snarl tree

            //Maps each net_handle_t to an index to its node problem, in all_node_problems
            hash_map<net_handle_t, size_t> net_handle_to_node_problem_index;
            //This stores all the snarl tree nodes and their clustering scratch work 
            vector<SnarlTreeNodeProblem> all_node_problems;
           
            //All chains for the current level of the snarl tree and gets updated as the algorithm
            //moves up the snarl tree. At one iteration, the algorithm will go through each chain
            //in chain to children and cluster the chain using clusters on the children
            vector<net_handle_t>* current_chains;


            //Same as current_chains but for the level of the snarl
            //tree above the current one
            //This gets updated as the current level is processed - the snarls from this level
            //are added as children to parent_chain_to_children.
            //After processing one level, this becomes the next chain_to_children
            vector<net_handle_t>* parent_chains;

            //All snarls for the current level of the snarl tree 
            //(chains from chain_to_children get added to their parent snarls, snarls get added to parent_snarls
            //then all snarls in snarl_to_children are clustered and added to parent_chain_to_children)
            vector<net_handle_t> parent_snarls;


            //This holds all the child problems of the root
            //Each pair is the parent and the child. This will be sorted by parent before
            //clustering
            vector<pair<net_handle_t, net_handle_t>> root_children;


            /////////////////////////////////////////////////////////

            //Constructor takes in a pointer to the seeds, the distance limits, and 
            //the total number of seeds in all_seeds
            ClusteringProblem (vector<vector<SeedCache>*>* all_seeds, 
                       size_t read_distance_limit, size_t fragment_distance_limit, size_t seed_count) :
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

                net_handle_to_node_problem_index.reserve(5*seed_count);
                all_node_problems.reserve(5*seed_count);
                parent_snarls.reserve(seed_count);
                root_children.reserve(seed_count);
            }
        };

        //Go through all the seeds and assign them to their parent chains or roots
        //If a node is in a chain, then assign it to its parent chain and add the parent
        //chain to chain_to_children_by_level
        //If a node is a child of the root or of a root snarl, then add cluster it and
        //remember to cluster the root snarl 
        void get_nodes( ClusteringProblem& clustering_problem,
                        vector<vector<net_handle_t>>& chains_by_level) const;


        //Cluster all the snarls at the current level
        void cluster_snarl_level(ClusteringProblem& clustering_problem) const;

        //Cluster all the chains at the current level
        //also assigns each chain to its parent and saves the distances to the ends of the parent
        //for each chain
        void cluster_chain_level(ClusteringProblem& clustering_problem, size_t depth) const;

        //Cluster the seeds on the specified node
        //The seeds are unsorted
        void cluster_one_node(ClusteringProblem& clustering_problem, SnarlTreeNodeProblem* node_problem) const; 

        //Cluster the seeds in a snarl
        void cluster_one_snarl(ClusteringProblem& clustering_problem, SnarlTreeNodeProblem* snarl_problem) const;

        //Cluster the seeds in a chain
        //The children are unsorted, they get sorted before clustering
        //If the children of the chain are only seeds on nodes, then cluster as if it is a node
        void cluster_one_chain(ClusteringProblem& clustering_problem, SnarlTreeNodeProblem* chain_problem,
            bool is_top_level_chain) const;

        //Helper function for adding the next seed to the chain clusters
        void add_seed_to_chain_problem(ClusteringProblem& clustering_problem, SnarlTreeNodeProblem* chain_problem,
                                        SnarlTreeNodeProblem::SnarlTreeChild& last_child,
                                        size_t& last_prefix_sum, size_t& last_length, size_t& last_chain_component_end, 
                                        vector<ClusterHead>& cluster_heads_to_add_again,
                                        bool& found_first_node, pair<bool, bool>& found_first_node_by_read,
                                        const SnarlTreeNodeProblem::SnarlTreeChild& current_child, bool is_first_child, bool is_last_child,
                                        bool skip_distances_to_ends) const;

        //Helper function for adding the next snarl to the chain clusters
        void add_snarl_to_chain_problem(ClusteringProblem& clustering_problem, SnarlTreeNodeProblem* chain_problem,
                                        SnarlTreeNodeProblem::SnarlTreeChild& last_child, 
                                        size_t& last_prefix_sum, size_t& last_length, size_t& last_chain_component_end, 
                                        vector<ClusterHead>& cluster_heads_to_add_again,
                                        bool& found_first_node, pair<bool, bool>& found_first_node_by_read,
                                        const SnarlTreeNodeProblem::SnarlTreeChild& current_child, bool is_first_child, bool is_last_child, 
                                        bool skip_distances_to_ends) const;

        //Cluster in the root - everything in clustering_problem.root_children 
        void cluster_root(ClusteringProblem& clustering_problem) const;

        //Cluster a list of seeds (SeedIndexes) that are on a single linear structure (node or chain)
        //Requires that the list of seeds are sorted relative to their position on the structure
        //This can be called on a chain if there are no nested seeds on the chain
        //left offset is the distance from the left side of the structure
        //if include_prefix_sum, then this is being called on a chain and the prefix sum must be added to the
        //distance_left of the seeds
        void cluster_seeds_on_linear_structure(ClusteringProblem& clustering_problem, 
                SnarlTreeNodeProblem* problem, size_t structure_length, bool include_prefix_sum, 
                bool skip_distances_to_ends) const;

        //Compare two children of the parent and combine their clusters, to create clusters in the parent
        //child_distances contains the distances for cluster heads in the children, 
        //since the distances in the seeds will get updated to be the distances in the parent
        //First child is true if this is the first time we see child_problem1. If first_child is true and this is 
        //a snarl, then we need to update the snarl's distances to its parents
        void compare_and_combine_cluster_on_child_structures(ClusteringProblem& clustering_problem, SnarlTreeNodeProblem* child_problem1,
                SnarlTreeNodeProblem* child_problem2, SnarlTreeNodeProblem* parent_problem, 
                const vector<pair<size_t, size_t>>& child_distances, bool is_root, bool first_child) const;

        //The same as above, but compare clusters on a single child
        //This assumes that the child is the child of the root and not a root snarl
        //so we just look at external distances 
        void compare_and_combine_cluster_on_one_child(ClusteringProblem& clustering_problem, SnarlTreeNodeProblem* child_problem) const;

};
}

#endif
