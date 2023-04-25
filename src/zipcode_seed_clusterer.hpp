#ifndef VG_ZIPCODE_SEED_CLUSTERER_HPP_INCLUDED
#define VG_ZIPCODE_SEED_CLUSTERER_HPP_INCLUDED

#include "snarl_seed_clusterer.hpp"

namespace vg {

    class ZipcodeClusterer{
        public:

        typedef SnarlDistanceIndexClusterer::Seed Seed;
        typedef SnarlDistanceIndexClusterer::Cluster Cluster;

        //Given a vector of seeds, coarsely cluster the seeds based on the distance in the graph
        //This is guaranteed to put seeds that are closer than the distance limit into the same
        //bucket, but may also put seeds that are far away in the same bucket
        vector<Cluster> coarse_cluster_seeds(const vector<Seed>& seeds, size_t distance_limit);

        private:
        const SnarlDistanceIndex* distance_index;
        const HandleGraph* graph;

        public:

        ZipcodeClusterer (const SnarlDistanceIndex* distance_index, const HandleGraph* graph) :
            distance_index(distance_index),
            graph(graph) {};

        ZipcodeClusterer (const SnarlDistanceIndex& distance_index, const HandleGraph& graph) :
            distance_index(&distance_index),
            graph(&graph) {};

        private:

        /*
         * Coarse clustering is done by partitioning the zipcodes
         * The zipcodes can be partially ordered along chains and snarls, so partitioning will be
         * done by walking along ordered lists of seeds and splitting the lists into different partitions 
         * Partitioning is done in a bfs traversal of the snarl tree
         */
 
        /////////////////////////////////// DATA STRUCTURES ////////////////////////////////////////////////

        /*
         * The partitions are stored using doubly linked lists. Each item in the list represents one seed,
         * which is represented as an index into the vector of seeds
         * Because partitioning is done top-down, the list will only change the current snarl tree node,
         * but the descendants will remain the same 
         */


        /// A node in a doubly linked list representing one seed
        struct partition_item_t {
            size_t seed;  //The index of the seed in a vector of seeds
            size_t prev;  //The index of the previous item in the list, as an index in the backing vector
            size_t next;  //The index of the next item in the linked list, std::numeric_limits<size_t>::max if it is the last

            //We need to be able to jump from the first seed in a snarl tree node to the last seed in the same node,
            //so that we don't traverse the whole list when partitioning its parent
            //If this is the first seed in a child with multiple seeds, then last_item_at_depth stores the index of the
            //last item in the child, as a pair of <depth, index>
            //Because the snarl tree is processed top-down, and all sorts on the vector will be stable sorts,
            // the index of the last thing will always be the same by the time we get to it
            //It gets stored in reverse order of depth (bottom up) so the back of the vector can be popped
            //TODO: This should maybe be the length, in case sorting gets messed up
            vector<pair<size_t, size_t>> last_item_at_depth; 
        };


        /// A partition_set_t stores a set of partitions of some data
        /// Each partition is a doubly linked list, and gets stored as the first thing in the list
        /// The ends of the lists aren't stored, but can be identified because their next pointers will
        ///  be std::numeric_limits<size_t>::max()
        /// The actual data is stored in a vector of partition_item_t's
        ///
        /// It is intended to be used for putting all data in at once, sorting all the data, then partitioning
        class partition_set_t {

            public:

            partition_set_t();

            //Add a new item to its own partition
            void add_new_item(size_t value);

            //Reserve space for the list
            void reserve(const size_t& size);

            ///Get the index of the next thing in a linked list, skipping to the next child at the same depth
            /// Returns std::numeric_limits<size_t>::max() if it's the end
            size_t get_last_index_at_depth( const size_t& current_index, const size_t& depth);

            /// Sorts everything in the range [range_start, range_end) using the comparator
            /// The range is specified by the index into data, not the index in a linked list
            /// Assumes that everything in the range is in the same partition, and keeps connections
            /// to whatever was attached outside of the range
            /// This changes the order of the vector between range_start and range_end. 
            /// Nothing else will be affected
            /// Uses std::stable_sort
            void sort (size_t range_start, size_t range_end,
                       std::function<bool(const partition_item_t& a, const partition_item_t& b)> cmp,
                       bool sort_everything=false);

            ///Split the partition containing range_start, to create a new partition
            ///starting at range_start
            ///Splitting changes the linked list, but not the order of the vector
            void split_partition (size_t range_start);

            ///Split the partition containing range_start and range_end,
            ///creating a new partition containing range_start and range_end 
            void split_partition (size_t range_start, size_t range_end);


            /////////////////////// DATA //////////////////////////////

            ///The actual data
            ///The order of nodes in the vector doesn't matter except when sorting
            vector<partition_item_t> data;

            /// The partitions of the data
            /// This stores the first node in the linked list of each partition
            /// as an index into data
            vector<size_t> partition_heads;
        };

        ///This holds the information of a new snarl/chain that needs to be partitioned
        ///range_start and range_end are indices into the data field of a partition_set_t 
        ///that specify a range of seeds that all belong to the same snarl/chain at the given depth 
        ///These get put in a queue of things that need to be partitioned, which is updated as the 
        ///algorithm proceeds
        struct partitioning_problem_t {
            size_t range_start;
            size_t range_end;
            size_t depth;
        };


        private:

        /*
         * The helper functions for doing the work of partitioning
         * coarse_cluster_seeds() will call these to coordinate partitioning 
         * Partitioning is split up by snarl/chain
         * These functions will pass around references to a partitioning_set_t of all partitions, 
         * and a queue of partitioning problems that need to be solved
         * Each will partition the given snarl or chain, and added partitioning problems for each child
         */

        /// Partition the seeds on a chain, specified by the current_problem 
        /// Each new partition that is made must be added to all_partitions, and
        ///  any children of the chain that need to be partitioned further must
        ///  be added to to_partition
        /// Assumes that the seeds in the range are sorted along the chain
        /// Doesn't alter the order of anything in all_partitions.data
        /// This should also handle nodes
        void partition_by_chain(const vector<Seed>& seeds,  
            const partitioning_problem_t& current_problem, 
            partition_set_t& all_partitions,
            std::list<partitioning_problem_t>& to_partition,
            const size_t& distance_limit);

        /// Partition the seeds on a snarl, specified by the current_problem 
        /// Each new partition that is made must be added to all_partitions, and
        ///  any children of the snarl that need to be partitioned further must
        ///  be added to to_partition
        /// Assumes that the seeds in the snarl are sorted by the distance to
        ///  the start of the snarl
        /// This may change the order of the snarl's children in the vector all_partitions.data,
        /// but the order of seeds within the children will remain the same
        void partition_by_snarl(const vector<Seed>& seeds,  
            const partitioning_problem_t& current_problem, 
            partition_set_t& all_partitions,
            std::list<partitioning_problem_t>& to_partition,
            const size_t& distance_limit);
    };
}
#endif
