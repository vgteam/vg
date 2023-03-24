#include "zipcode_seed_clusterer.hpp"

//#define DEBUG_ZIPCODE_CLUSTERING

namespace vg {

vector<ZipcodeClusterer::Cluster> ZipcodeClusterer::cluster_seeds(const vector<Seed>& seeds, size_t distance_limit ) {
#ifdef DEBUG_ZIPCODE_CLUSTERING
    cerr << endl << endl << "New zipcode clustering of " << seeds.size() << " seeds with distance limit" <<  distance_limit << endl;
#endif
    //Bucket the seeds roughly by their distance along the top-level chain

    vector<Cluster> clusters;

    /*First, sort the seeds by their connected component, and by the distance along the top-level chain (or other long chain)
    */

    //This will hold information from a seed for sorting and partitioning
    struct seed_values_t {
        size_t index;               //Index into seeds
        size_t connected_component; //Connected component identifier
        size_t prefix_sum;          //Prefix sum of the thing on the top-level chain
        size_t length;              //length of the thing on the top-level chain
    };

    //Make a vector of seed_value_t's and fill in the index of the seed and distance values 
    vector<seed_values_t> sorted_indices (seeds.size());
    for (size_t i = 0 ; i < sorted_indices.size() ; i++) {
        if (seeds[i].zipcode.byte_count() == 0) {
            //If the zipcode is empty, then fill it in
            cerr << "warning: Can't cluster empty zipcodes" << endl;
            return clusters;
        }
        sorted_indices[i].index = i;
        sorted_indices[i].connected_component = seeds[i].zipcode_decoder->get_distance_index_address(0);

        if (seeds[i].zipcode_decoder->get_code_type(0) == ROOT_CHAIN) { 
            //If this is in a top-level chain, then store the offset and length
            sorted_indices[i].prefix_sum = seeds[i].zipcode_decoder->get_offset_in_chain(1);
            sorted_indices[i].length = seeds[i].zipcode_decoder->get_length(1);
        } else {
            //If this is in a top-level snarl, then it all goes into the same cluster so these don't matter
            sorted_indices[i].prefix_sum = std::numeric_limits<size_t>::max();
            sorted_indices[i].length = std::numeric_limits<size_t>::max();
        }
    }

    //Sort
    std::sort(sorted_indices.begin(), sorted_indices.end(), [&] (const seed_values_t& a, const seed_values_t& b) {
        //Comparator for sorting. Returns a < b
        if (a.connected_component == b.connected_component){ 
            //If they are on the same connected component, then check the offset in the top-level chain
            //If this is a top-level snarl, then both prefix sum values are max(), because the order 
            //doesn't matter
            return a.prefix_sum < b.prefix_sum;
        } else if (a.connected_component < b.connected_component) {
            return true;
        } else {
            return false;
        }
    });
#ifdef DEBUG_ZIPCODE_CLUSTERING
    cerr << "Sorted seeds:" << endl; 
    for (const seed_values_t& this_seed : sorted_indices) {
        cerr << seeds[this_seed.index].pos << " " << this_seed.prefix_sum << " " << this_seed.length << endl;
    }
    cerr << endl;
#endif

    /*Next, walk through the sorted list of seeds and partition
    */
    seed_values_t empty_seed = {std::numeric_limits<size_t>::max(),
                              std::numeric_limits<size_t>::max(),
                              std::numeric_limits<size_t>::max(),
                              std::numeric_limits<size_t>::max()}; 
    seed_values_t& last_seed = empty_seed;
    for (const seed_values_t& this_seed : sorted_indices) {
#ifdef DEBUG_ZIPCODE_CLUSTERING
        cerr << "At seed " << seeds[this_seed.index].pos << endl;
#endif
        if (last_seed.index == std::numeric_limits<size_t>::max()) {
            //If this is the first seed in the sorted list, then make a new cluster
#ifdef DEBUG_ZIPCODE_CLUSTERING
            cerr << "\tthis is the first seed so make a new cluster" << endl;
#endif
            clusters.emplace_back();
            clusters.back().seeds.emplace_back(this_seed.index);
        } else if (last_seed.connected_component != this_seed.connected_component) {
            //If this is on a new connected component, make a new cluster
#ifdef DEBUG_ZIPCODE_CLUSTERING
            cerr << "\tthis is on a new connected component so make a new cluster" << endl;
#endif
            clusters.emplace_back();
            clusters.back().seeds.emplace_back(this_seed.index);
        } else if (last_seed.prefix_sum == std::numeric_limits<size_t>::max() ||
                   this_seed.prefix_sum == std::numeric_limits<size_t>::max()) {

#ifdef DEBUG_ZIPCODE_CLUSTERING
            cerr << "\tone or both prefix sums are max() so put them in the same cluster" << endl;
#endif
            clusters.back().seeds.emplace_back(this_seed.index);
        } else if (SnarlDistanceIndex::minus(this_seed.prefix_sum,
                                             SnarlDistanceIndex::sum(last_seed.prefix_sum, last_seed.length)) 
                   > distance_limit) {
#ifdef DEBUG_ZIPCODE_CLUSTERING
            cerr << "\tthis is too far from the last seed so make a new cluster" << endl;
            cerr << "Last prefix sum: " << last_seed.prefix_sum << " last length " << last_seed.length << " this prefix sum: " << this_seed.prefix_sum << endl;
#endif
            //If too far from the last seed, then put it in a new cluster
            clusters.emplace_back();
            clusters.back().seeds.emplace_back(this_seed.index);
        } else {
            //If they are on the same component and close enough, add this seed to the last cluster
#ifdef DEBUG_ZIPCODE_CLUSTERING
            cerr << "\tthis was close to the last seed so add it to the previous cluster" << endl;
            cerr << "Last prefix sum: " << last_seed.prefix_sum << " last length " << last_seed.length << " this prefix sum: " << this_seed.prefix_sum << endl;
#endif
            clusters.back().seeds.emplace_back(this_seed.index);
        }
        last_seed = this_seed;
    }

    return clusters;
}

}
