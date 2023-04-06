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

    };
}
#endif

