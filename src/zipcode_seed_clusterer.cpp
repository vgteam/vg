#include "zipcode_seed_clusterer.hpp"

//#define DEBUG_ZIPCODE_CLUSTERING

namespace vg {


/*
 * Coarsely cluster the seeds using their zipcodes
 * All seeds start out in the same partition and are split into different partitions according to their position on the snarl tree
 * Seeds are first ordered recursively along the snarl tree - along chains and according to the distance to the start of a snarl.
 * Snarls/chains are found by walking along the ordered list of seeds and processed in a bfs traversal of the snarl tree
 * This is accomplished using a queue of partitioning_problem_t's, which represent the next snarl tree node to partition.
 * All partitions are maintained in a partition_set_t, which is processed into clusters at the end
 */
vector<ZipcodeClusterer::Cluster> ZipcodeClusterer::coarse_cluster_seeds(const vector<Seed>& seeds, size_t distance_limit ) {
#ifdef DEBUG_ZIPCODE_CLUSTERING
    cerr << endl << endl << "New zipcode clustering of " << seeds.size() << " seeds with distance limit" <<  distance_limit << endl;
#endif

    //This holds all the partitions found. It gets processed into clusters at the end
    partition_set_t all_partitions;

    //A queue of everything that needs to be partitioned. Each item represents the seeds in a single snarl tree node
    //The snarl tree gets processed in a bfs traversal
    std::list<partitioning_problem_t> to_partition;

    /* First, initialize the problem with one partition for each connected component
     *
     * Sort the seeds by their position in the snarl tree
     * The seeds are sorted first by connected component, by position along a chain, by the distance to the start of a snarl,
     * and by the rank in the snarl. 
     * Then walk through the ordered list of seeds and add last_item_at_depth for skipping to the ends of snarl tree nodes, 
     * and split by connected component and create a new partitioning_problem_t in to_partition for each connected component
     */

    //This is the first partition containing all the seeds
    all_partitions.reserve(seeds.size());
    for (size_t i = 0 ; i < seeds.size() ; i++) {
        all_partitions.add_new_item(i);
    }

    //Sort
    all_partitions.sort(0, seeds.size(), [&] (const partition_item_t& a, const partition_item_t& b) {
        //Comparator for sorting. Returns a < b
        size_t depth = 0;
        while (depth < seeds[a.seed].zipcode_decoder->decoder_length()-1 &&
               depth < seeds[b.seed].zipcode_decoder->decoder_length()-1 &&
               ZipCodeDecoder::is_equal(*seeds[a.seed].zipcode_decoder, *seeds[b.seed].zipcode_decoder, depth)) {
            depth++;
        }
        //Either depth is the last thing in a or b, or they are different at this depth 
        if ( ZipCodeDecoder::is_equal(*seeds[a.seed].zipcode_decoder, *seeds[b.seed].zipcode_decoder, depth)) {
            //If they are equal
            return false;
        } else if (depth == 0) {
            //If they are on different connected components, sort by connected component
            return seeds[a.seed].zipcode_decoder->get_distance_index_address(0) < seeds[b.seed].zipcode_decoder->get_distance_index_address(0);
            
        } else if (seeds[a.seed].zipcode_decoder->get_code_type(depth-1) == CHAIN || seeds[a.seed].zipcode_decoder->get_code_type(depth-1) == ROOT_CHAIN) {
            //If a and b are both children of a chain
            size_t offset_a = seeds[a.seed].zipcode_decoder->get_offset_in_chain(depth);
            size_t offset_b = seeds[b.seed].zipcode_decoder->get_offset_in_chain(depth);
            if ( offset_a == offset_b) {
                //If they have the same prefix sum, then the snarl comes first
                return seeds[a.seed].zipcode_decoder->get_code_type(depth) != NODE && seeds[b.seed].zipcode_decoder->get_code_type(depth) == NODE;  
            } else {
                return offset_a < offset_b;
            }
        } else if (seeds[a.seed].zipcode_decoder->get_code_type(depth-1) == REGULAR_SNARL) {
            //If the parent is a regular snarl, then sort by child number
            return seeds[a.seed].zipcode_decoder->get_rank_in_snarl(depth) < seeds[b.seed].zipcode_decoder->get_rank_in_snarl(depth);
        } else {
            //Otherwise, they are children of an irregular snarl
            return seeds[a.seed].zipcode_decoder->get_distance_to_snarl_start(depth) < seeds[b.seed].zipcode_decoder->get_distance_to_snarl_start(depth);
        }
    });

#ifdef DEBUG_ZIPCODE_CLUSTERING
    cerr << "Sorted seeds:" << endl; 
    for (auto& index : all_partitions.data) {
        size_t this_seed = all_partitions.data[index].seed;
        cerr << seeds[this_seed.index].pos << " " << this_seed.prefix_sum << " " << this_seed.length << endl;
    }
    cerr << endl;
#endif

    //Partition by connected_component and create a new partitioning_problem_t for each
    //Also update last_item_at_depth for each item. For each seed that is the first seed for a particular child, 
    //store the length of that child and its depth

    //A list of the index of the first seed in a snarl tree node at each depth. This is used to fill in last_item_at_depth
    //Initialized to be 0 for all snarl tree nodes of the first seed
    std::vector<size_t> first_zipcode_at_depth (seeds[all_partitions.data[0].seed].zipcode_decoder->decoder_length(), 0);

    //The beginning of the connected component we're currently on 
    size_t last_connected_component_start = 0;

    for (size_t i = 1 ; i < all_partitions.data.size() ; i++ ) {

        auto& current_decoder = *seeds[all_partitions.data[i].seed].zipcode_decoder;
        size_t current_depth = current_decoder.decoder_length();

        //For any snarl tree node that ends here, add it's last_item_at_depth
        for (int depth = first_zipcode_at_depth.size() ; depth >= 0 ; depth--) {
            if (current_depth > depth ||
                !ZipCodeDecoder::is_equal(current_decoder, *seeds[all_partitions.data[i-1].seed].zipcode_decoder, depth)) {
                //If the previous thing was in a different snarl tree node at this depth

                if (first_zipcode_at_depth[depth] != i-1 ) {
                    //If the first seed in this child wasn't the seed right before this one
                    //Add the number of things that were in that snarl tree node
                    all_partitions.data[first_zipcode_at_depth[depth]].last_item_at_depth.emplace_back(depth, i - first_zipcode_at_depth[depth]);
                }
                first_zipcode_at_depth[depth] = i;

            }
        }
        if (current_depth > first_zipcode_at_depth.size()) {
            //We need to add things
            while (first_zipcode_at_depth.size() <= current_depth) {
                first_zipcode_at_depth.emplace_back(i);
            }
        } else if (current_depth > first_zipcode_at_depth.size()) {
            //We need to remove things
            while (first_zipcode_at_depth.size() > current_depth+1) {
                first_zipcode_at_depth.pop_back();
            }
        }

        //Now check if this is the start of a new connected component
        if (!ZipCodeDecoder::is_equal(*seeds[all_partitions.data[i-1].seed].zipcode_decoder, 
                                      current_decoder, 0)) {
            //If these are on different connected components

            //Make a new partition at i
            all_partitions.split_partition(i);

            //Remember to partition everything from the start to i-1
            to_partition.push_back({last_connected_component_start, i-1, 0});

            //i is the new start of the current partition
            last_connected_component_start = i;
             

            //Update the first zipcode at each depth
            first_zipcode_at_depth.assign (current_decoder.decoder_length(), i);
        }
    }

    /*
     * Now go through all the partitioning_problem_t's and solve them
     * partition_by_chain/snarl will add to to_partition as they go
     */
     
    while (!to_partition.empty()) {

        //Get the next problem from the front of the queue 
        const auto& current_problem = to_partition.front();
        //Remove it from the queue
        to_partition.pop_front();

        code_type_t code_type = seeds[all_partitions.data[current_problem.range_start].seed].zipcode_decoder->get_code_type(current_problem.depth);

        if (code_type == CHAIN || code_type == NODE) {
            partition_by_chain(seeds, current_problem, all_partitions, to_partition, distance_limit);
        } else {
            partition_by_snarl(seeds, current_problem, all_partitions, to_partition, distance_limit);
        }

    }

    
    /* When there is nothing left in to_partition, partitioning is done.
     * Go through all partitions and create clusters
     */
     vector<Cluster> all_clusters;
     all_clusters.reserve(all_partitions.partition_heads.size());
     for (const size_t& cluster_head : all_partitions.partition_heads) {
         all_clusters.emplace_back();

         partition_item_t& current_item = all_partitions.data[cluster_head];
         while (current_item.next != std::numeric_limits<size_t>::max()){
            all_clusters.back().seeds.emplace_back(current_item.seed);
            current_item = all_partitions.data[current_item.next]; 
         }
         all_clusters.back().seeds.emplace_back(current_item.seed);
     }

    return all_clusters;
}

/* Partition the given problem along a chain
 * The seeds in the current_problem must be sorted along the chain
 * Chains are split when the distance between subsequent seeds is definitely larger than the distance_limit
 */

void ZipcodeClusterer::partition_by_chain(const vector<Seed>& seeds, const partitioning_problem_t& current_problem,
    partition_set_t& all_partitions, std::list<partitioning_problem_t>& to_partition,
    const size_t& distance_limit){

    const size_t& depth = current_problem.depth;

    //We're going to walk through the seeds on children of the chain, starting from the second one
    size_t previous_index = current_problem.range_start;
    partition_item_t& previous_item = all_partitions.data[previous_index];

    //First, check if we actually have to do any work
    if (previous_item.next == std::numeric_limits<size_t>::max() ||
        seeds[previous_item.seed].zipcode_decoder->get_length(depth) <= distance_limit) {
        //If there was only one seed, or the chain is too short, then don't do anything
        return;
    }

    //Get the index of the next partition_item_t in the chain
    size_t current_index = all_partitions.get_last_index_at_depth(previous_index, depth);

    //If the first seed was in a snarl with other seeds, then remember to partition the snarl
    if (all_partitions.data[current_index].prev != previous_index) {
        to_partition.push_back({previous_index, all_partitions.data[current_index].prev, depth+1});
    }

    /*Walk through the sorted list of seeds and partition
    */
    while (current_index != std::numeric_limits<size_t>::max()) {

#ifdef DEBUG_ZIPCODE_CLUSTERING
        cerr << "At seed " << seeds[all_partitions.data[current_index].seed].pos << endl;
#endif

        //Get the values we need to calculate distance
        size_t current_prefix_sum = seeds[all_partitions.data[current_index].seed].zipcode_decoder->get_offset_in_chain(depth);
        size_t previous_prefix_sum = seeds[all_partitions.data[previous_index].seed].zipcode_decoder->get_offset_in_chain(depth);
        size_t previous_length = seeds[all_partitions.data[previous_index].seed].zipcode_decoder->get_length(depth);

        if (previous_prefix_sum != std::numeric_limits<size_t>::max() &&
            current_prefix_sum != std::numeric_limits<size_t>::max() &&
            SnarlDistanceIndex::minus(current_prefix_sum,
                                      SnarlDistanceIndex::sum(previous_prefix_sum, previous_length)) 
                   > distance_limit) {

#ifdef DEBUG_ZIPCODE_CLUSTERING
            cerr << "\tthis is too far from the last seed so make a new cluster" << endl;
            cerr << "\tLast prefix sum: " << previous_prefix_sum << " last length " << previous_length << " this prefix sum: " << current_prefix_sum << endl;
#endif
            //If too far from the last seed, then split off a new cluster
            all_partitions.split_partition(current_index);
        }
#ifdef DEBUG_ZIPCODE_CLUSTERING
        else {
            cerr << "\tthis is close enough to the last thing, so it is in the same cluster" << endl;
            cerr << "\tLast prefix sum: " << previous_prefix_sum << " last length " << previous_length << " this prefix sum: " << current_prefix_sum << endl;
        }
#endif

        //Update to the next thing in the list
        previous_index = current_index;

        //Check if this was the last thing in the range
        if (current_index == current_problem.range_end) {
            //If this is the last thing we wanted to process
            current_index = std::numeric_limits<size_t>::max();
        } else {
            //Otherwise, get the next thing, skipping other things in the same child at this depth
            current_index = all_partitions.get_last_index_at_depth(previous_index, depth+1);

            //If this skipped a snarl in the chain, then remember to cluster it later
            if (all_partitions.data[current_index].prev != previous_index) {
                to_partition.push_back({previous_index, all_partitions.data[current_index].prev, depth+1});
            }

#ifdef DEBUG_ZIPCODE_CLUSTERING
            if (current_index == std::numeric_limits<size_t>::max()) {
                assert(previous_index == current_problem.range_end);
            }
#endif
        }
    }

    return;
}

/*
 * Snarls are split in two passes over the seeds. First, they are sorted by the distance to the start of the snarl and
 * split if the difference between the distances to the start is greater than the distance limit
 * For each child, x, in a snarl, we know the minimum distance to the start and end boundary nodes of the snarl (x_start and x_end)
 * For two children of the snarl, x and y, assume that x_start <= y_start.
 * Then there can be no path from x to y that is less than (y_start - x_start), otherwise y_start would be smaller. So y_start-x_start is a lower bound of the distance from x to y
 */

void ZipcodeClusterer::partition_by_snarl(const vector<Seed>& seeds, const partitioning_problem_t& current_problem,
    partition_set_t& all_partitions, std::list<partitioning_problem_t>& to_partition,
    const size_t& distance_limit){

    const size_t& depth = current_problem.depth;

    //We're going to walk through the seeds on children of the chain, starting from the second one
    size_t previous_index = current_problem.range_start;
    partition_item_t& previous_item = all_partitions.data[previous_index];

    //First, check if we actually have to do any work
    if (previous_item.next == std::numeric_limits<size_t>::max() ||
        seeds[previous_item.seed].zipcode_decoder->get_length(depth) <= distance_limit) {
        //If there was only one seed, or the chain is too short, then don't do anything
        return;
    }

    //Get the index of the next partition_item_t in the chain
    size_t current_index = all_partitions.get_last_index_at_depth(previous_index, depth);

    //If the first seed was in a snarl with other seeds, then remember to partition the snarl
    if (all_partitions.data[current_index].prev != previous_index) {
        to_partition.push_back({previous_index, all_partitions.data[current_index].prev, depth+1});
    }
}

ZipcodeClusterer::partition_set_t::partition_set_t() {
}

//Move constructor
//ZipcodeClusterer::partition_set_t::partition_set_t(partition_set_t&& other) :
//    data(std::move(other.data)), head(other.head), tail(other.tail) {
//    other.data = std::vector<partition_item_t>(0);
//    other.head = nullptr;
//    other.tail = nullptr; 
//}

void ZipcodeClusterer::partition_set_t::add_new_item(size_t value) {
    data.push_back({value, 
                      std::numeric_limits<size_t>::max(), 
                      std::numeric_limits<size_t>::max()});
}
void ZipcodeClusterer::partition_set_t::reserve(const size_t& size) {
    data.reserve(size);
}


size_t ZipcodeClusterer::partition_set_t::get_last_index_at_depth(const size_t& current_index, const size_t& depth) {
    partition_item_t& current_item = data[current_index];
    if (current_item.next == std::numeric_limits<size_t>::max()) {
        return std::numeric_limits<size_t>::max();
    } else if (current_item.last_item_at_depth.empty() || 
            current_item.last_item_at_depth.back().first < depth) {
        //If there are no other children at this depth
        return current_item.next;
    } else {
        while (current_item.last_item_at_depth.back().first > depth) {
            current_item.last_item_at_depth.pop_back();
        }
        const pair<size_t, size_t>& last = current_item.last_item_at_depth.back();
        current_item.last_item_at_depth.pop_back();
        return data[current_index + last.second - 1].next;
    }
}


void ZipcodeClusterer::partition_set_t::sort(size_t range_start, size_t range_end, std::function<bool(const partition_item_t& a, const partition_item_t& b)> cmp, bool sort_everything) {


    //Sort the vector
    std::stable_sort(data.begin()+range_start, data.begin()+range_end, cmp);

    //Connections to outside of the range. May be max() if the start or end of a list was in the range
    size_t prev, next;

    //If the start of the range was in the range, then we need to replace it as the start of a list in partitions
    size_t old_start = std::numeric_limits<size_t>::max();

    
    //Make sure that everything points to the proper thing
    for (size_t i = 0 ; i < data.size() ; i++) {

        if (!sort_everything) {
            //Remember if anything pointed to outside the range
            if (data[i].prev == std::numeric_limits<size_t>::max()) {
                old_start = i;
                prev = std::numeric_limits<size_t>::max();
            } else if (data[i].prev < range_start) {
                prev = data[i].prev;
            }
            if (data[i].next > range_end || data[i].next == std::numeric_limits<size_t>::max()) {
                next = data[i].next;
            }
        }
    
        data[i].prev = i == 0 ? std::numeric_limits<size_t>::max() : i-1;
        data[i].next = i == data.size()-1 ? std::numeric_limits<size_t>::max() : i+1;
    }

    if (sort_everything) {
        //If we sorted the whole list, then everything is in the same partition
        partition_heads.clear();
        partition_heads.emplace_back(0);
    } else {
        if (prev != std::numeric_limits<size_t>::max()) {
            //If the start of the list was outside the range

            //Make sure the list is connected from the start
            data[prev].next = range_start;
            data[range_start].prev = prev;
        } else {
            //If the start of the list was in the range, then we need to replace the start of the linked list in partition_heads 
            for (size_t i = 0 ; i < partition_heads.size() ; i++) {
                if (partition_heads[i] == old_start) {
                    partition_heads[i] = range_start;
                    break;
                }
            }
        }

        if (next != std::numeric_limits<size_t>::max()) {
            // If the end of the list was outside the range, update the end
            data[next].prev = range_end;
            data[range_end].next = next;
        }
    }

    return;
}

void ZipcodeClusterer::partition_set_t::split_partition(size_t range_start) {
    if (data[range_start].prev == std::numeric_limits<size_t>::max()) {
        //If this is the first thing in a list
        return;
    } else {
        //Otherwise, tell the previous thing that it's now the end of a linked list, and add this one as a new partition

        //Update previous to be the last thing in it's list
        data[data[range_start].prev].next = std::numeric_limits<size_t>::max();

        //Tell range_start that it's the start of a list
        data[range_start].prev = std::numeric_limits<size_t>::max();

        //Add range_start as a new partition
        partition_heads.emplace_back(range_start);
        
    }
}

void ZipcodeClusterer::partition_set_t::split_partition(size_t range_start, size_t range_end) {
    if (data[range_start].prev == std::numeric_limits<size_t>::max() && data[range_end].next == std::numeric_limits<size_t>::max()) {
        //If this is the whole list
        return;
    } else if (data[range_start].prev == std::numeric_limits<size_t>::max()) {
        //If this is the start of an existing list, then start a new one after range_end

        //Update the next head to know it's a head
        data[ data[range_end].next ].prev = std::numeric_limits<size_t>::max();

        //Tell range_end that it's now the end
        data[range_end].next = std::numeric_limits<size_t>::max();

        //Add the next thing as a new partition
        partition_heads.emplace_back(range_end+1);
    } else if (data[range_end].next == std::numeric_limits<size_t>::max()) {
        //This is the end of a partition
        split_partition(range_start);
    } else {
        //Otherwise, this is in the middle of a partition and we need to update the previous and next things to point to each other

        //Update previous and next things to point to each other
        size_t previous = data[range_start].prev;
        size_t next = data[range_end].next;

        data[previous].next = next;
        data[next].prev = previous;

        //Tell range_start and range end that they're the start/end of a list
        data[range_start].prev = std::numeric_limits<size_t>::max();
        data[range_end].next = std::numeric_limits<size_t>::max();

        //Add range_start as a new partition
        partition_heads.emplace_back(range_start);
        
    }
}


}
