#include "snarl_seed_clusterer.hpp"

#include <algorithm>

//#define DEBUG_CLUSTER
//#define debug_distances
namespace vg {

SnarlDistanceIndexClusterer::SnarlDistanceIndexClusterer( const SnarlDistanceIndex& distance_index, const HandleGraph* graph) :
                                        distance_index(distance_index),
                                        graph(graph){
};
SnarlDistanceIndexClusterer::SnarlDistanceIndexClusterer( const SnarlDistanceIndex* distance_index, const HandleGraph* graph) :
                                        distance_index(*distance_index),
                                        graph(graph){
};
SnarlDistanceIndexClusterer::SnarlDistanceIndexClusterer( const SnarlDistanceIndex& distance_index) :
                                        distance_index(distance_index),
                                        graph(nullptr){
};
SnarlDistanceIndexClusterer::SnarlDistanceIndexClusterer( const SnarlDistanceIndex* distance_index) :
                                        distance_index(*distance_index),
                                        graph(nullptr){
};

vector<SnarlDistanceIndexClusterer::Cluster> SnarlDistanceIndexClusterer::cluster_seeds (vector<Seed>& seeds, size_t read_distance_limit) const {
    //Wrapper for single ended

    vector<vector<Seed>*> all_seeds;
    all_seeds.push_back(&seeds);
    std::vector<std::vector<size_t>> all_clusters =
        std::get<0>(cluster_seeds_internal(all_seeds, read_distance_limit, 0))[0].all_groups();

    std::vector<Cluster> result;
    result.reserve(all_clusters.size());
    for (auto& cluster : all_clusters) {
        result.emplace_back();
        result.back().seeds = std::move(cluster);
    }
    //TODO: Sorting fixes determinism issues but seems unecessary
    std::sort(result.begin(), result.end(), [&] (Cluster& cluster1, Cluster& cluster2) {
        return cluster1.seeds.front() < cluster2.seeds.front();
    });

    return result;
};

vector<vector<SnarlDistanceIndexClusterer::Cluster>> SnarlDistanceIndexClusterer::cluster_seeds (
              vector<vector<Seed>>& all_seeds, size_t read_distance_limit,
              size_t fragment_distance_limit) const {
    //Wrapper for paired end

    if (all_seeds.size() > 2) {
        throw std::runtime_error("Clusterer: We can't handle more than paired end mapping");
    }
    vector<vector<Seed>*> seed_pointers;
    seed_pointers.reserve(all_seeds.size());
    for (vector<Seed>& v : all_seeds) seed_pointers.push_back(&v);

    //Actually cluster the seeds
    auto union_finds = cluster_seeds_internal(seed_pointers, read_distance_limit, fragment_distance_limit);

    vector<structures::UnionFind>* read_union_finds = &std::get<0>(union_finds);
    structures::UnionFind* fragment_union_find = &std::get<1>(union_finds);

    std::vector<std::vector<Cluster>> result (all_seeds.size());
    //Map the old group heads to new indices
    size_t curr_index = 0;
    size_t read_num_offset = 0;
    hash_map<size_t, size_t> old_to_new_cluster_index;

    for (size_t read_num = 0 ; read_num < read_union_finds->size() ; read_num++) {
        vector<vector<size_t>> read_clusters = read_union_finds->at(read_num).all_groups();
        result[read_num].reserve(read_clusters.size());
        for (vector<size_t>& cluster : read_clusters) {
            result[read_num].emplace_back();
            Cluster& curr = result[read_num].back();
            curr.seeds = std::move(cluster);
        }
    //TODO: Sorting fixes determinism issues but seems unecessary
        std::sort(result[read_num].begin(), result[read_num].end(), [&] (Cluster& cluster1, Cluster& cluster2) {
            return cluster1.seeds.front() < cluster2.seeds.front();
        });
    }
    for (size_t read_num = 0 ; read_num < result.size() ; read_num++) {
        for (Cluster& cluster : result[read_num]) {
            size_t fragment_index = read_num_offset + cluster.seeds[0];
            size_t fragment_cluster_head = fragment_union_find->find_group(fragment_index);
            if (old_to_new_cluster_index.count(fragment_cluster_head) == 0) {
                old_to_new_cluster_index.emplace(fragment_cluster_head, curr_index);
                fragment_cluster_head = curr_index;
                curr_index++;
            } else {
                fragment_cluster_head = old_to_new_cluster_index[fragment_cluster_head];
            }
            cluster.fragment = fragment_cluster_head;
        }
        read_num_offset += all_seeds[read_num].size();
    }


    return result;
}


tuple<vector<structures::UnionFind>, structures::UnionFind> SnarlDistanceIndexClusterer::cluster_seeds_internal (
              vector<vector<Seed>*>& all_seeds, size_t read_distance_limit,
              size_t fragment_distance_limit) const {
    /* Given a vector of seeds and a limit, find a clustering of seeds where
     * seeds that are closer than the limit cluster together.
     * Returns a vector of clusters
     */
#ifdef DEBUG_CLUSTER
cerr << endl << endl << endl << endl << "New cluster calculation:" << endl;
cerr << "\tread distance limit: " << read_distance_limit << " and fragment distance limit: " << fragment_distance_limit << endl;
#endif
    if (fragment_distance_limit != 0 &&
        fragment_distance_limit < read_distance_limit) {
        throw std::runtime_error("Fragment distance limit must be greater than read distance limit");
    }

    //For each level of the snarl tree, maps chains at that level to children belonging to the chain.
    //Children are represented as a size_t index into all_node_clusters, or two size_t indexes into all_seeds
    //Starts out with nodes from get_nodes, 
    //then adds snarls as we walk up the snarl tree in 
    vector<ParentToChildMap> chain_to_children_by_level;
    chain_to_children_by_level.reserve(distance_index.get_max_tree_depth()+1);



    //This stores all the tree relationships and cluster information
    //for a single level of the snarl tree as it is being processed
    //It also keeps track of the parents of the current level
    size_t seed_count = 0;
    for (auto v : all_seeds) seed_count+= v->size();
    TreeState tree_state (&all_seeds, read_distance_limit, fragment_distance_limit, seed_count);


    //Initialize chain_to_children_by_level with all the seeds on chains
    //Also clusters seeds on nodes in the root or root snarls and adds them to the root snarls
    get_nodes(tree_state, chain_to_children_by_level);

    //Initialize the tree state to the bottom level
    tree_state.chain_to_children = &chain_to_children_by_level[chain_to_children_by_level.size() - 1];

    for (int depth = chain_to_children_by_level.size() - 1 ; depth >= 0 ; depth --) {
        // Go through each level of the tree, bottom up, and cluster that level. 
        // When we reach a level, we know all the children of the chains at that level
        // Cluster each chain, assign chains to parent snarls
        // Cluster each snarl, assign snarls to parent chains
        // Reset the current level to point to the parent chains

#ifdef DEBUG_CLUSTER
assert(tree_state.seed_count_prefix_sum[0] == 0);
for (size_t i = 1 ; i < tree_state.all_seeds->size() ; i++) {
    assert (tree_state.seed_count_prefix_sum[i] + tree_state.all_seeds->at(i)->size() == tree_state.seed_count_prefix_sum[i+1]);
}
#endif
        if (depth != 0) {
            tree_state.parent_chain_to_children = &chain_to_children_by_level[depth-1];
        }


        //Cluster all the chains at this depth
        //Also records which chains are in snarls and the parents of these
        //chains in tree_state.parent_chain_to_children
        cluster_chain_level(tree_state, depth);

        //And cluster all the snarls, record the parents of these snarls
        cluster_snarl_level(tree_state);


        // Swap buffer over for the next level
        tree_state.chain_to_children = tree_state.parent_chain_to_children;
    }
    //There may be some connectivity in the root, so also try to cluster in the root
    cluster_root(tree_state);
   


#ifdef DEBUG_CLUSTER

    cerr << "Found read clusters : " << endl;
    for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
        cerr << "\t read num " << read_num << ": " ;
        for (auto group : tree_state.read_union_find[read_num].all_groups()){
            cerr << "\t\t";
            for (size_t c : group) {
               cerr << tree_state.all_seeds->at(read_num)->at(c).pos << " ";
            }
            cerr << endl;
        }
        cerr << endl;
    }
    vector<Seed> ordered_seeds;
    for (size_t i = 0 ; i < tree_state.all_seeds->size() ; i++) {
        const auto v = tree_state.all_seeds->at(i);
        for ( auto x : *v) {
            ordered_seeds.push_back(x);
        }
    }
    cerr << "Found fragment clusters : " << endl;
    for (auto group : tree_state.fragment_union_find.all_groups()){
        cerr << "\t";
        for (size_t c : group) {
           cerr << ordered_seeds[c].pos << " ";
        }
        cerr << endl;
    }

/*
    //CHeck read clusters
    for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
        auto all_groups =  tree_state.read_union_find[read_num].all_groups();
        for (size_t g1 = 0 ; g1 < all_groups.size() ; g1 ++ ){
            auto group = all_groups[g1];
            structures::UnionFind uf(group.size(), false);
            for (size_t i1 = 0 ; i1 < group.size() ; i1++) {
                size_t c = group[i1];
                pos_t pos1 = tree_state.all_seeds->at(read_num)->at(c).pos;
                pos_t rev1 = make_pos_t(get_id(pos1), !is_rev(pos1), distance_index.node_length(get_id(pos1)) - get_offset(pos1) - 1);

                for (size_t i2 = 0 ; i2 < i1 ; i2++) {
                
                    size_t d = group[i2];
                
                    pos_t pos2 = tree_state.all_seeds->at(read_num)->at(d).pos;
                    pos_t rev2 = make_pos_t(get_id(pos2), !is_rev(pos2), distance_index.node_length(get_id(pos2))- get_offset(pos2) - 1);
                    size_t d1 = distance_index.min_distance(pos1, pos2);
                    size_t d2 = std::min(d1, distance_index.min_distance(pos1, rev2));
                    size_t d3 = std::min(d2, distance_index.min_distance(rev1, rev2));
                    size_t d4 = std::min(d3, distance_index.min_distance(rev1, pos2));
                    if (d4 != -1 && d4 <= tree_state.read_distance_limit) {
                    
                         uf.union_groups(i1, i2);
                    }
                }
                for (size_t g2 = 0 ; g2 < all_groups.size() ; g2 ++) {
                    if (g2 != g1) {
                        auto group2 = all_groups[g2];
                        for (size_t d : group2) {
                           pos_t pos2 = tree_state.all_seeds->at(read_num)->at(d).pos;
                           pos_t rev2 = make_pos_t(get_id(pos2), !is_rev(pos2), distance_index.node_length(get_id(pos2)) - get_offset(pos2) - 1);
                           
                           size_t d1 = distance_index.min_distance(pos1, pos2);
                           size_t d2 = std::min(d1, distance_index.min_distance(pos1, rev2));
                           size_t d3 = std::min(d2, distance_index.min_distance(rev1, rev2));
                           size_t d4 = std::min(d3, distance_index.min_distance(rev1, pos2));
                           
                           assert (d4 == -1 || d4 > tree_state.read_distance_limit);
                        }  
                    
                    }
                }
            }
            if (uf.all_groups().size() != 1) {
                cerr << "These should be separate clusters: " << endl;
                for (auto uf_group : uf.all_groups()) {
                    for (size_t i : uf_group) {
                        size_t c = group[i];
                        cerr << tree_state.all_seeds->at(read_num)->at(c).pos << ":" << tree_state.all_seeds->at(read_num)->at(c).component << ":"
                             << tree_state.all_seeds->at(read_num)->at(c).offset << ", ";
                    }
                    cerr << endl;
                }

            }
            assert (uf.all_groups().size() == 1);
        }
    }
    */


#endif
    return make_tuple(std::move(tree_state.read_union_find), std::move(tree_state.fragment_union_find));

};


//Go through all the seeds and assign them to their parent chains or roots
//If a node is in a chain, then assign it to its parent chain and add the parent
//chain to chain_to_children_by_level
//If a node is a child of the root or of a root snarl, then add cluster it and
//remember to cluster the root snarl 
void SnarlDistanceIndexClusterer::get_nodes( TreeState& tree_state, vector<ParentToChildMap>& chain_to_children_by_level) const {
#ifdef DEBUG_CLUSTER
cerr << "Add all seeds to nodes: " << endl;
#endif

    //This is to remember the nodes that we are going to cluster at the end of get_nodes
    //these will be the nodes that are children of the root or root snarl. 
    //All other seeds are added directly to their parent chains as children
    //Used as a map from each node (as an index into all_node_clusters) to all of its seeds
    //seeds are represented by indexes into the seeds vector (read_num, seed_num)
    vector<std::tuple<size_t,size_t, size_t>> node_to_seeds;


    //Map the parent NodeClusters to its depth so we don't use get_depth() as much
    vector<size_t> parent_node_cluster_offset_to_depth;
    parent_node_cluster_offset_to_depth.reserve(tree_state.seed_count_prefix_sum.back());


    //All nodes we've already assigned
    hash_set<id_t> seen_nodes;
    seen_nodes.reserve(tree_state.seed_count_prefix_sum.back());

    for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++){ 
        vector<Seed>* seeds = tree_state.all_seeds->at(read_num);
        for (size_t i = 0; i < seeds->size(); i++) {
            Seed& seed = seeds->at(i);
            pos_t pos = seed.pos;
            id_t id = get_id(pos);
            

#ifdef DEBUG_CLUSTER
            cerr << "\t" << read_num << ":" << pos << ", ";
#endif


            //We are going to add the seed to its parent.
            //If the node is in the root, then cluster it after going through all seeds and forget about it
            //If the parent is a proper chain, then add the seed directly to the parent chain

            
            //Remember the new cache. We will replace it with the actual values for the seed as necessary
            //cached values are:
            //(0)record offset of node, (1)record offset of parent, (2)node record offset, (3)node length, (4)is_reversed, 
            // (5)is_trivial_chain, (6)parent is chain, (7)parent is root, (8)prefix sum, (9)chain_component
            MIPayloadValues old_cache = seed.minimizer_cache;

            //TODO: For now, we're either storing all values or none
            bool has_cached_values = old_cache.record_offset != MIPayload::NO_VALUE;


            //Get the net_handle for the node the seed is on
            net_handle_t node_net_handle = !has_cached_values ? distance_index.get_node_net_handle(id) 
                                            : distance_index.get_net_handle_from_values(old_cache.record_offset, 
                                                             SnarlDistanceIndex::START_END, 
                                                             SnarlDistanceIndex::NODE_HANDLE, 
                                                             old_cache.node_record_offset); 


            //Get the parent of the node
            net_handle_t parent;
            //If the grandparent is a root/root snarl, then make it the parent and the node a trivial chain 
            //because they will be clustered here and added to the root instead of being added to the 
            //snarl tree to be clustered
            if (has_cached_values) {
                if (old_cache.is_trivial_chain) {
                    //If the node is a trivial chain, then the parent is just the node but recorded as a chain in the net handle
                    parent = distance_index.get_net_handle_from_values (distance_index.get_record_offset(node_net_handle),
                                                            SnarlDistanceIndex::START_END,
                                                            SnarlDistanceIndex::CHAIN_HANDLE,
                                                            old_cache.node_record_offset);
                    if (old_cache.parent_record_offset == 0) {
                        //If the parent offset stored in the cache is the root, then this is a trivial chain
                        //child of the root not in a root snarl, so remember the root as the parent and the 
                        //trivial chain as th enode
                        node_net_handle = parent;
                        parent = distance_index.get_root();
                    } else if (old_cache.parent_is_root && !old_cache.parent_is_chain) {
                        //If the parent is a root snarl, then the node becomes the trivial chain 
                        //and we get the parent root snarl from the cache
                        node_net_handle = parent;
                        parent = distance_index.get_net_handle_from_values(old_cache.parent_record_offset,
                                                                           SnarlDistanceIndex::START_END,
                                                                           SnarlDistanceIndex::ROOT_HANDLE);
                    }
                } else if (old_cache.parent_record_offset == 0) {
                    //The parent is just the root
                    parent = distance_index.get_root();
                } else if (old_cache.parent_is_root && !old_cache.parent_is_chain) {
                    //If the parent is a root snarl
                    parent = distance_index.get_net_handle_from_values(old_cache.parent_record_offset,
                                                           SnarlDistanceIndex::START_END,
                                                           SnarlDistanceIndex::ROOT_HANDLE);
                } else {
                    //Otherwise the parent is an actual chain and we use the value from the cache
                    parent = distance_index.get_net_handle_from_values(old_cache.parent_record_offset,
                                                           SnarlDistanceIndex::START_END,
                                                           SnarlDistanceIndex::CHAIN_HANDLE);
                }
            } else {
                parent = distance_index.start_end_traversal_of(distance_index.get_parent(node_net_handle));
                if (distance_index.is_trivial_chain(parent)){
                    net_handle_t grandparent = distance_index.get_parent(parent);
                    if (distance_index.is_root(grandparent)){
                        node_net_handle = parent;
                        parent = distance_index.start_end_traversal_of(grandparent);
                    }
                }
            }

#ifdef DEBUG_CLUSTER
cerr << old_cache.is_trivial_chain << " " << old_cache.parent_is_chain << " " << old_cache.parent_is_root << endl;
cerr << distance_index.net_handle_as_string(node_net_handle) << " parent: " << distance_index.net_handle_as_string(parent) << endl;
            if (!distance_index.is_root(parent)) {
                cerr << "Parent should be " << distance_index.net_handle_as_string(distance_index.start_end_traversal_of(distance_index.get_parent(node_net_handle))) << endl; 
                assert( distance_index.start_end_traversal_of(parent) == distance_index.start_end_traversal_of(distance_index.get_parent(node_net_handle)));
            }
#endif
            if (!distance_index.is_root(parent)) {
                //If the parent is not the root and not a root snarl (it is a chain or trivial chain)

#ifdef DEBUG_CLUSTER
                cerr << "\tchild of a chain " << distance_index.net_handle_as_string(parent) << endl;
#endif

                //Add the seed to its parent
                //Also update the minimizer_cache on the seed 



                //Seed payload is: 
                //record offset of node, record offset of parent, node record offset, node length, is_reversed, is_trivial_chain, parent is chain, parent is root, prefix sum, chain_component

                bool is_trivial_chain = has_cached_values ? old_cache.is_trivial_chain 
                                                      : distance_index.is_trivial_chain(parent);
                size_t prefix_sum = old_cache.prefix_sum;
                size_t node_length = old_cache.node_length;
                bool is_reversed_in_parent = old_cache.is_reversed;

                if (!has_cached_values) {
                    //If we didn't store information in the seed, then get it from the distance index
                    //and remember it in the seed's cache

                    //prefix sum
                    prefix_sum = is_trivial_chain ? std::numeric_limits<size_t>::max() 
                                                  : distance_index.get_prefix_sum_value(node_net_handle);
                    seed.minimizer_cache.prefix_sum = prefix_sum;

                    //component
                    seed.minimizer_cache.chain_component = distance_index.is_multicomponent_chain(parent) 
                            ? distance_index.get_chain_component(node_net_handle)
                            : 0;

                    //node length
                    node_length = distance_index.minimum_length(node_net_handle);
                    seed.minimizer_cache.node_length = node_length;

                    //is_reversed_in_parent
                    is_reversed_in_parent = is_trivial_chain ? distance_index.is_reversed_in_parent(parent)
                                                             : distance_index.is_reversed_in_parent(node_net_handle);
                    seed.minimizer_cache.is_reversed = is_reversed_in_parent;

                }
#ifdef DEBUG_CLUSTER
                //assert(prefix_sum == (is_trivial_chain ? std::numeric_limits<size_t>::max() 
                //                                  : distance_index.get_prefix_sum_value(node_net_handle)));
                assert(node_length == distance_index.minimum_length(node_net_handle));

                assert(is_reversed_in_parent == (is_trivial_chain ? distance_index.is_reversed_in_parent(parent)
                                                             : distance_index.is_reversed_in_parent(node_net_handle)));
#endif

                //Add the parent chain or trivial chain
                size_t parent_index;
                bool new_parent = false;
                size_t depth;
                if (old_cache.is_trivial_chain && old_cache.parent_is_chain && old_cache.parent_is_root) {
                    //If the node is a trivial chain, and the parent we stored is a chain and root,
                    //then the node is in a simple snarl on the root-level chain
                    depth = 2;
                } else if (old_cache.parent_is_root) {
                    //If the parent is a root (or root-level chain)
                    depth = 1;
                } else {
                    //Otherwise get it later from parent_node_cluster_offset_to_depth
                    depth = std::numeric_limits<size_t>::max();
                }
                if (tree_state.net_handle_to_index.count(parent) == 0) {
                    //If we haven't seen the parent chain before, make a new NodeClusters for it
                    parent_index = tree_state.all_node_clusters.size();
                    tree_state.net_handle_to_index[parent] = parent_index;
                    if (is_trivial_chain ) {
                        tree_state.all_node_clusters.emplace_back(parent, tree_state.all_seeds->size(),
                                                     tree_state.seed_count_prefix_sum.back(),
                                                     false, id, node_length, std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()); 
                        tree_state.all_node_clusters.back().is_trivial_chain = true;
                    } else {
                        //The parent is an actual chain
                        tree_state.all_node_clusters.emplace_back(parent, tree_state.all_seeds->size(),
                                                              tree_state.seed_count_prefix_sum.back(), distance_index);
                    }

                    //Get the depth from the parent if we didn't cache it
                    if (depth == std::numeric_limits<size_t>::max()) {
                        depth = distance_index.get_depth(parent);
                    }
                    parent_node_cluster_offset_to_depth.emplace_back(depth);
                    new_parent = true;
                } else {
                    //If we've seen the parent before, just find its index into all_node_clusters and its depth
                    parent_index = tree_state.net_handle_to_index[parent];
                    if (depth == std::numeric_limits<size_t>::max()) {
                        depth = parent_node_cluster_offset_to_depth[parent_index];
                    }
                }
#ifdef DEBUG_CLUSTER
                assert(depth == distance_index.get_depth(parent));
#endif


                //If chain_to_children_by_level isn't big enough for this depth, resize it and reserve space at each level
                if (depth+1 > chain_to_children_by_level.size()) {
                    size_t to_add = (depth+1) - chain_to_children_by_level.size(); 
                    for (size_t i = 0 ; i < to_add ; i++) {
                        chain_to_children_by_level.emplace_back();
                        chain_to_children_by_level.back().reserve(tree_state.seed_count_prefix_sum.back());
                    }
                }

                //Make sure the seed's distances are relative to the orientation in the parent
                seed.distance_left = is_reversed_in_parent != is_rev(pos) ? node_length- get_offset(pos) 
                                                                          : get_offset(pos) + 1;
                seed.distance_right = is_reversed_in_parent != is_rev(pos) ? get_offset(pos) + 1 
                                                                           : node_length- get_offset(pos);
                //Add this seed to its parent cluster
                chain_to_children_by_level[depth].add_child(parent_index, node_net_handle, read_num, i, 
                              seed.minimizer_cache.chain_component,  //component
                              is_trivial_chain ? seed.distance_left : SnarlDistanceIndex::sum(prefix_sum, seed.distance_left));


                //If the parent is a trivial chain and not in the root, then we also stored the identity of the snarl, so add it here too
                if (new_parent && has_cached_values) {
                    if (is_trivial_chain && !old_cache.parent_is_root) {
                        bool grandparent_is_simple_snarl = old_cache.parent_is_chain;
                        tree_state.all_node_clusters[parent_index].has_parent_handle = true;
                        tree_state.all_node_clusters[parent_index].parent_net_handle = grandparent_is_simple_snarl 
                                                            ? distance_index.get_net_handle_from_values(distance_index.get_record_offset(node_net_handle),
                                                                                            SnarlDistanceIndex::START_END,
                                                                                            SnarlDistanceIndex::SNARL_HANDLE,
                                                                                            1)
                                                            : distance_index.get_net_handle_from_values(old_cache.parent_record_offset,
                                                                                            SnarlDistanceIndex::START_END,
                                                                                            SnarlDistanceIndex::SNARL_HANDLE);

                        if (grandparent_is_simple_snarl) {
                            //If the grandparent is a simple snarl, then we also stored the identity of its parent chain, so add it here too
                            tree_state.all_node_clusters[parent_index].has_grandparent_handle = true;
                            tree_state.all_node_clusters[parent_index].grandparent_net_handle = distance_index.get_net_handle_from_values(old_cache.parent_record_offset,
                                                                                                SnarlDistanceIndex::START_END,
                                                                                                SnarlDistanceIndex::CHAIN_HANDLE);
                        }
                    } else if (old_cache.parent_is_root && old_cache.parent_is_chain && !is_trivial_chain) {
                        //The parent chain is a child of the root
                        tree_state.all_node_clusters[parent_index].has_parent_handle = true;
                        tree_state.all_node_clusters[parent_index].parent_net_handle = 
                                distance_index.get_net_handle_from_values(0, SnarlDistanceIndex::START_END, SnarlDistanceIndex::ROOT_HANDLE);
                    }
                }


            } else {
                //Otherwise, the parent is the root or a root snarl, and the node_net_handle is a node


                //Get the values from the seed. Some may be infinite and need to be re-set
                size_t node_length = has_cached_values ? old_cache.node_length
                                                       : distance_index.minimum_length(node_net_handle);
                bool is_reversed_in_parent = has_cached_values ? old_cache.is_reversed
                                                         : distance_index.is_reversed_in_parent(node_net_handle);


                size_t node_index;
                //Create a new NodeClusters for this node
                if (seen_nodes.count(id) == 0) {
                    node_index = tree_state.all_node_clusters.size();
                    tree_state.all_node_clusters.emplace_back(
                                NodeClusters(std::move(node_net_handle), tree_state.all_seeds->size(),
                                             tree_state.seed_count_prefix_sum.back(),
                                             false, id, node_length, std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
                    //Remember the parent of this node, since it will be needed to remember the root snarl later
                    tree_state.all_node_clusters.back().parent_net_handle = parent;

                    //The depth is never used but add it just to match the length of all_node_clusters
                    parent_node_cluster_offset_to_depth.emplace_back(0);

                    seen_nodes.insert(id);

                    //Remember the index of the NodeClusters
                     tree_state.net_handle_to_index[node_net_handle] = node_index; 
                } else {
                    node_index = tree_state.net_handle_to_index[node_net_handle]; 
                }

                //Remember this seed as a child of the node
                node_to_seeds.emplace_back(node_index, read_num, i);
                seed.distance_left = is_reversed_in_parent != is_rev(pos) ? node_length- get_offset(pos) : get_offset(pos) + 1;
                seed.distance_right = is_reversed_in_parent != is_rev(pos) ? get_offset(pos) + 1 : node_length- get_offset(pos);
            }
        }
    }
    //Sort node_to_seeds by the node and the seeds by their offset
    std::sort(node_to_seeds.begin(), node_to_seeds.end(), [&](const auto& a, const auto b) -> bool {
        if (std::get<0>(a) == std::get<0>(b) ) { 
            return  tree_state.all_seeds->at(std::get<1>(a))->at(std::get<2>(a)).distance_left <
                    tree_state.all_seeds->at(std::get<1>(b))->at(std::get<2>(b)).distance_left;
        } else {
            return std::get<0>(a) < std::get<0>(b);
       } 
    });

#ifdef DEBUG_CLUSTER
    cerr << endl;
#endif

    //Go through and cluster nodes that are children of the root or root snarls
    std::vector<std::tuple<size_t, size_t, size_t>>::iterator seed_range_start = node_to_seeds.begin();
    for(auto current_iterator = node_to_seeds.begin() ; current_iterator != node_to_seeds.end() ; ++current_iterator) {
        if (std::next(current_iterator) == node_to_seeds.end() 
             || std::get<0>(*current_iterator) != std::get<0>(*seed_range_start)) {
            //If the current iterator is pointing to a different node than seed_range_start, then we just passed the last
            //seed for that node and we can cluster it
            NodeClusters& node_clusters = tree_state.all_node_clusters[std::get<0>(*seed_range_start)];

            //Cluster the node. Give it the range in node_to_seeds, which is from seed_range_start
            //to either current_iterator (if current_iterator is a different node), or the end of node_to_seeds
            //if current_iterator is the last thing in the list and the same node
            cluster_one_node(tree_state, node_clusters, seed_range_start, 
                              (std::get<0>(*current_iterator) != std::get<0>(*seed_range_start) ? current_iterator
                                                                               : node_to_seeds.end()));

            net_handle_t parent = node_clusters.parent_net_handle;

            if (distance_index.is_root_snarl(parent)) {
                //If this is a root snarl, then remember it to cluster in the root
                size_t parent_index;
                if (tree_state.net_handle_to_index.count(parent) == 0) {
                    parent_index = tree_state.all_node_clusters.size();
                    tree_state.net_handle_to_index[parent] = parent_index;
                    tree_state.all_node_clusters.emplace_back(parent, tree_state.all_seeds->size(),
                                                              tree_state.seed_count_prefix_sum.back(), distance_index);
                } else {
                    parent_index = tree_state.net_handle_to_index[parent];
                }
                tree_state.root_children.emplace_back(parent_index, std::get<0>(*seed_range_start));
            } else {
                //Otherwise, just compare the single child's external connectivity
                compare_and_combine_cluster_on_one_child(tree_state, tree_state.all_node_clusters.back());
            }
            if (std::next(current_iterator) == node_to_seeds.end() 
                 && std::get<0>(*current_iterator) != std::get<0>(*seed_range_start)) {
                //If this is the last thing in the list and is a different node than seed_range_start,
                //then deal cluster the last node. Even though it will just be one cluster,
                //we still need to cluster it to get the distances to the ends of the node and 
                //check connectivity in the root

                NodeClusters& node_clusters = tree_state.all_node_clusters[std::get<0>(*current_iterator)];

                //Cluster the node
                cluster_one_node(tree_state, node_clusters, current_iterator, node_to_seeds.end());

                net_handle_t parent = node_clusters.parent_net_handle;

                if (distance_index.is_root_snarl(parent)) {
                    //If this is a root snarl, then remember it to cluster in the root
                    size_t parent_index;
                    if (tree_state.net_handle_to_index.count(parent) == 0) {
                        parent_index = tree_state.all_node_clusters.size();
                        tree_state.net_handle_to_index[parent] = parent_index;
                        tree_state.all_node_clusters.emplace_back(parent, tree_state.all_seeds->size(),
                                                                  tree_state.seed_count_prefix_sum.back(), distance_index);
                    } else {
                        parent_index = tree_state.net_handle_to_index[parent];
                    }
                    tree_state.root_children.emplace_back(parent_index, std::get<0>(*current_iterator));
                } else {
                    //Otherwise, just compare the single child's external connectivity
                    compare_and_combine_cluster_on_one_child(tree_state, tree_state.all_node_clusters.back());
                }
            }

            seed_range_start = current_iterator;
        }
    }

    if (chain_to_children_by_level.empty()) {
        chain_to_children_by_level.resize(1);
    }
}



//Cluster all of the snarls in tree_state from the same depth
//Assumes that all the children of the snarls have been clustered already and are present in tree_state.snarls_to_children
void SnarlDistanceIndexClusterer::cluster_snarl_level(TreeState& tree_state) const {

    
    //An iterator into snarl_to_children (which is a multimap from snarl to children, each as indices into all_node_clusters)
    //This gets updated to the next snarl in each iteration until it reaches the end
    std::multimap<size_t, size_t>::iterator child_list_iterator = tree_state.snarl_to_children.begin();
    std::multimap<size_t, size_t>::iterator child_list_iterator_end = tree_state.snarl_to_children.end();

    while (child_list_iterator != child_list_iterator_end){
        //Go through each of the snarls at this level, cluster them,
        //and find which chains they belong to, if any

        size_t snarl_index = (*child_list_iterator).first;

        //An iterator to the last entry in the multimap with this snarl index (exclusive)
        std::multimap<size_t, size_t>::iterator next_child_list_iterator = tree_state.snarl_to_children.upper_bound(snarl_index);

#ifdef DEBUG_CLUSTER
        cerr << "Cluster one snarl " << distance_index.net_handle_as_string(tree_state.all_node_clusters[snarl_index].containing_net_handle) << endl;
#endif

        //Cluster the snarl
        NodeClusters& snarl_clusters = tree_state.all_node_clusters[snarl_index];
        cluster_one_snarl(tree_state, snarl_index, child_list_iterator, next_child_list_iterator);

        /*Now add the snarl to its parent. Only do so if the clusters are close enough to the boundaries that it can be clustered*/

        //Check the best distance of any seed to the ends of the snarl
        //Is the distance small enough that we can cluster it with something else?
        bool reachable_right = snarl_clusters.fragment_best_right <= 
            (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit);
        bool reachable_left = snarl_clusters.fragment_best_left <= 
            (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit);


        if (reachable_left || reachable_right) {

            //Make a new NodeClusters for the parent

            net_handle_t snarl_parent = snarl_clusters.has_parent_handle
                                      ? snarl_clusters.parent_net_handle
                                      : distance_index.start_end_traversal_of(distance_index.get_parent(snarl_clusters.containing_net_handle));
            size_t parent_index;
            if (tree_state.net_handle_to_index.count(snarl_parent) == 0) {
                parent_index = tree_state.all_node_clusters.size();
                tree_state.net_handle_to_index[snarl_parent] = parent_index;
                tree_state.all_node_clusters.emplace_back(snarl_parent, tree_state.all_seeds->size(),
                                                      tree_state.seed_count_prefix_sum.back(), distance_index);
                if (snarl_clusters.has_grandparent_handle) {
                    tree_state.all_node_clusters.back().has_parent_handle = true;
                    tree_state.all_node_clusters.back().parent_net_handle = snarl_clusters.grandparent_net_handle;
                }
            } else {
                parent_index = tree_state.net_handle_to_index[snarl_parent];
            }

            //Add the snarl to its parent
            if (distance_index.is_root(snarl_parent)) {
                 if(distance_index.is_root_snarl(snarl_parent)) {
                    //If the parent is a root snarl, then remember it to be compared in the root
                    tree_state.root_children.emplace_back(parent_index, snarl_index);
                 } else {
                     //Otherwise, compare it to itself using external connectivity and forget about it since we're done
                     compare_and_combine_cluster_on_one_child(tree_state, snarl_clusters);
                 }
            } else {
                //Add the snarl to its parent chain
                tree_state.parent_chain_to_children->add_child(parent_index, snarl_clusters.containing_net_handle, 
                                                               snarl_index, std::numeric_limits<size_t>::max(), 
                                                               snarl_clusters.chain_component_start, snarl_clusters.prefix_sum_value);
            }
        }
        child_list_iterator = std::move(next_child_list_iterator);

#ifdef DEBUG_CLUSTER
        cerr << "\tRecording snarl " << distance_index.net_handle_as_string(snarl_clusters.containing_net_handle)  << " as a child of "
              << distance_index.net_handle_as_string(distance_index.get_parent(snarl_clusters.containing_net_handle)) << endl;
#endif

    }
    tree_state.snarl_to_children.clear();
}


void SnarlDistanceIndexClusterer::cluster_chain_level(TreeState& tree_state, size_t depth) const {

    //Go through chain_to_children, which is a vector of chain, child pairs. Start by sorting by parent chain
    tree_state.chain_to_children->sort(distance_index);
    vector<ParentToChildMap::ParentChildValues>& chain_to_children = tree_state.chain_to_children->parent_to_children;

    if (chain_to_children.empty()) {
        return;
    }


    //Go through the list of children, where the children are represented as a tuple of parent, child
    //Keep an iterator to the  first item of the current chain, and cluster a chain when we find its last child
    std::vector<ParentToChildMap::ParentChildValues>::iterator chain_range_start = chain_to_children.begin();
    bool only_seeds=true;

    for (std::vector<ParentToChildMap::ParentChildValues>::iterator current_iterator = chain_to_children.begin() ; current_iterator != chain_to_children.end() ; ++ current_iterator) {


        if (current_iterator+1 == chain_to_children.end() || chain_range_start->parent_index != current_iterator->parent_index) {
            //If the current iterator is for a different chain than chain_range_start, 
            //or if it's the last thing in the list, then cluster everything from chain_range_start to
            //the current iterator (inclusive if it's the last thing in the list and the same chain)

            net_handle_t chain_handle = tree_state.all_node_clusters[chain_range_start->parent_index].containing_net_handle;
#ifdef DEBUG_CLUSTER
            cerr << "Cluster one chain " <<  distance_index.net_handle_as_string(chain_handle) << endl;
#endif


            net_handle_t parent = tree_state.all_node_clusters[chain_range_start->parent_index].has_parent_handle
                                ? tree_state.all_node_clusters[chain_range_start->parent_index].parent_net_handle
                                : distance_index.start_end_traversal_of(distance_index.get_parent(chain_handle));
            bool is_root = distance_index.is_root(parent);
            bool is_root_snarl = is_root ? distance_index.is_root_snarl(parent) : false;

            //This is used to determine if we need to remember the distances to the ends of the chain, since
            //for a top level chain it doesn't matter
            bool is_top_level_chain = (depth == 1) && !is_root_snarl &&
                             !distance_index.is_externally_start_start_connected(chain_handle) &&
                             !distance_index.is_externally_start_end_connected(chain_handle) &&
                             !distance_index.is_externally_end_end_connected(chain_handle) &&
                             !distance_index.is_looping_chain(chain_handle);

            // Compute the clusters for the chain
            if (current_iterator+1 == chain_to_children.end() 
                && chain_range_start->parent_index == current_iterator->parent_index) {
                //If this is the last thing in chain_to_children and the chain is the same, 
                //then the last iterator should be the end since it has to be exclusive
                only_seeds = only_seeds && (current_iterator->child_index2 != std::numeric_limits<size_t>::max());
                cluster_one_chain(tree_state, chain_range_start->parent_index, chain_range_start, 
                                  chain_to_children.end(), only_seeds, is_top_level_chain);
            } else {
                cluster_one_chain(tree_state, chain_range_start->parent_index, chain_range_start, 
                                  current_iterator, only_seeds, is_top_level_chain);
            }

            //Add the chain to its parent
            if (is_root) {
                //If the parent is the root, remember the index of this chain in all_node_clusters
                if (is_root_snarl) {
                    //If the parent is a root snarl, then remember it to cluster in the root
                    size_t parent_index;
                    if (tree_state.net_handle_to_index.count(parent) == 0) {
                        parent_index = tree_state.all_node_clusters.size();
                        tree_state.net_handle_to_index[parent] = parent_index;
                        tree_state.all_node_clusters.emplace_back(parent, tree_state.all_seeds->size(),
                                                              tree_state.seed_count_prefix_sum.back(), distance_index);
                    } else {
                        parent_index = tree_state.net_handle_to_index[parent];
                    }
                    tree_state.root_children.emplace_back(parent_index,chain_range_start->parent_index);
                } else if (!is_top_level_chain) {
                    //Otherwise, cluster it with itself using external connectivity only
                    //is_top_level_chain also includes external connectivity, so if it's true we don't need to check this
                     compare_and_combine_cluster_on_one_child(tree_state, tree_state.all_node_clusters[chain_range_start->parent_index]);
                }
            } else if (!is_top_level_chain) {
                //If the parent is just a snarl

                //Remember the distances to the ends of the parent 

                NodeClusters& chain_clusters = tree_state.all_node_clusters[chain_range_start->parent_index];
                chain_clusters.distance_start_left = 
                        distance_index.distance_to_parent_bound(parent, true, distance_index.flip(chain_handle),
                            std::make_tuple(SnarlDistanceIndex::SNARL_HANDLE, 
                                            SnarlDistanceIndex::SNARL_HANDLE,
                                            (chain_clusters.is_trivial_chain ? SnarlDistanceIndex::NODE_HANDLE 
                                                                             : SnarlDistanceIndex::CHAIN_HANDLE), 
                                            SnarlDistanceIndex::CHAIN_HANDLE));

                chain_clusters.distance_start_right = 
                        distance_index.distance_to_parent_bound(parent, true, chain_handle, 
                            std::make_tuple(SnarlDistanceIndex::SNARL_HANDLE, 
                                            SnarlDistanceIndex::SNARL_HANDLE,
                                            (chain_clusters.is_trivial_chain ? SnarlDistanceIndex::NODE_HANDLE 
                                                                             : SnarlDistanceIndex::CHAIN_HANDLE), 
                                            SnarlDistanceIndex::CHAIN_HANDLE));

                chain_clusters.distance_end_left =
                        distance_index.distance_to_parent_bound(parent, false, distance_index.flip(chain_handle), 
                            std::make_tuple(SnarlDistanceIndex::SNARL_HANDLE, 
                                            SnarlDistanceIndex::SNARL_HANDLE,
                                            (chain_clusters.is_trivial_chain ? SnarlDistanceIndex::NODE_HANDLE 
                                                                             : SnarlDistanceIndex::CHAIN_HANDLE), 
                                            SnarlDistanceIndex::CHAIN_HANDLE));

                chain_clusters.distance_end_right = 
                        distance_index.distance_to_parent_bound(parent, false, chain_handle,
                            std::make_tuple(SnarlDistanceIndex::SNARL_HANDLE, 
                                            SnarlDistanceIndex::SNARL_HANDLE,
                                            (chain_clusters.is_trivial_chain ? SnarlDistanceIndex::NODE_HANDLE 
                                                                             : SnarlDistanceIndex::CHAIN_HANDLE), 
                                            SnarlDistanceIndex::CHAIN_HANDLE));
#ifdef DEBUG_CLUSTER
                cerr << "This child has distances to end : " << chain_clusters.distance_start_left << " " << chain_clusters.distance_start_right 
                     << " " << chain_clusters.distance_end_left << " " << chain_clusters.distance_end_right << endl;
#endif
                //And add it to its parent snarl
                size_t parent_index;
                if (tree_state.net_handle_to_index.count(parent) == 0) {
                    parent_index = tree_state.all_node_clusters.size();
                    tree_state.net_handle_to_index[parent] = parent_index;
                    tree_state.all_node_clusters.emplace_back(parent, tree_state.all_seeds->size(),
                                                              tree_state.seed_count_prefix_sum.back(), distance_index);
                    if (chain_clusters.has_grandparent_handle) {
                        tree_state.all_node_clusters.back().has_parent_handle = true;
                        tree_state.all_node_clusters.back().parent_net_handle = chain_clusters.grandparent_net_handle;
                    }
                } else {
                    parent_index = tree_state.net_handle_to_index[parent];
                }
                tree_state.snarl_to_children.emplace(parent_index, chain_range_start->parent_index);


            }

            if (current_iterator+1 == chain_to_children.end() && chain_range_start->parent_index != current_iterator->parent_index) {
                //If this is the last thing and it isn't the same chain as the previous thing, then we cluster just the last thing 
                //TODO: This is very repetative but I don't know a better way of doing it 

                net_handle_t chain_handle = tree_state.all_node_clusters[current_iterator->parent_index].containing_net_handle;
#ifdef DEBUG_CLUSTER
                cerr << "Cluster one chain " <<  distance_index.net_handle_as_string(chain_handle) << endl;
#endif


                net_handle_t parent = tree_state.all_node_clusters[current_iterator->parent_index].has_parent_handle
                                    ? tree_state.all_node_clusters[current_iterator->parent_index].parent_net_handle
                                    : distance_index.start_end_traversal_of(distance_index.get_parent(chain_handle));
                bool is_root = distance_index.is_root(parent);
                bool is_root_snarl = is_root ? distance_index.is_root_snarl(parent) : false;
                bool is_top_level_chain = (depth == 1) && !is_root_snarl &&
                                 !distance_index.is_externally_start_start_connected(chain_handle) &&
                                 !distance_index.is_externally_start_end_connected(chain_handle) &&
                                 !distance_index.is_externally_end_end_connected(chain_handle) &&
                                 !distance_index.is_looping_chain(chain_handle);

                only_seeds =current_iterator->child_index2 != std::numeric_limits<size_t>::max();
                cluster_one_chain(tree_state, current_iterator->parent_index, current_iterator, 
                                  chain_to_children.end(), only_seeds, is_top_level_chain);

                //Add the chain to its parent
                if (is_root) {
                    //If the parent is the root, remember the index of this chain in all_node_clusters
                    if (is_root_snarl) {
                        //If the parent is a root snarl, then remember it to cluster in the root
                        size_t parent_index;
                        if (tree_state.net_handle_to_index.count(parent) == 0) {
                            parent_index = tree_state.all_node_clusters.size();
                            tree_state.net_handle_to_index[parent] = parent_index;
                            tree_state.all_node_clusters.emplace_back(parent, tree_state.all_seeds->size(),
                                                                  tree_state.seed_count_prefix_sum.back(), distance_index);
                        } else {
                            parent_index = tree_state.net_handle_to_index[parent];
                        }
                        tree_state.root_children.emplace_back(parent_index,current_iterator->parent_index);
                    } else if (!is_top_level_chain) {
                        //Otherwise, cluster it with itself using external connectivity only
                        //is_top_level_chain also includes external connectivity, so if it's true we don't need to check this
                         compare_and_combine_cluster_on_one_child(tree_state, tree_state.all_node_clusters[current_iterator->parent_index]);
                    }
                } else if (!is_top_level_chain) {
                    //If the parent is just a snarl

                    //Remember the distances to the ends of the parent 

                    NodeClusters& chain_clusters = tree_state.all_node_clusters[current_iterator->parent_index];
                    chain_clusters.distance_start_left = 
                            distance_index.distance_to_parent_bound(parent, true, distance_index.flip(chain_handle), 
                            std::make_tuple(SnarlDistanceIndex::SNARL_HANDLE, 
                                            SnarlDistanceIndex::SNARL_HANDLE,
                                            (chain_clusters.is_trivial_chain ? SnarlDistanceIndex::NODE_HANDLE 
                                                                             : SnarlDistanceIndex::CHAIN_HANDLE), 
                                            SnarlDistanceIndex::CHAIN_HANDLE));

                    chain_clusters.distance_start_right = 
                            distance_index.distance_to_parent_bound(parent, true, chain_handle, 
                            std::make_tuple(SnarlDistanceIndex::SNARL_HANDLE, 
                                            SnarlDistanceIndex::SNARL_HANDLE,
                                            (chain_clusters.is_trivial_chain ? SnarlDistanceIndex::NODE_HANDLE 
                                                                             : SnarlDistanceIndex::CHAIN_HANDLE), 
                                            SnarlDistanceIndex::CHAIN_HANDLE));

                    chain_clusters.distance_end_left =
                            distance_index.distance_to_parent_bound(parent, false, distance_index.flip(chain_handle), 
                            std::make_tuple(SnarlDistanceIndex::SNARL_HANDLE, 
                                            SnarlDistanceIndex::SNARL_HANDLE,
                                            (chain_clusters.is_trivial_chain ? SnarlDistanceIndex::NODE_HANDLE 
                                                                             : SnarlDistanceIndex::CHAIN_HANDLE), 
                                            SnarlDistanceIndex::CHAIN_HANDLE));

                    chain_clusters.distance_end_right = 
                            distance_index.distance_to_parent_bound(parent, false, chain_handle, 
                            std::make_tuple(SnarlDistanceIndex::SNARL_HANDLE, 
                                            SnarlDistanceIndex::SNARL_HANDLE,
                                            (chain_clusters.is_trivial_chain ? SnarlDistanceIndex::NODE_HANDLE 
                                                                             : SnarlDistanceIndex::CHAIN_HANDLE), 
                                            SnarlDistanceIndex::CHAIN_HANDLE));
#ifdef DEBUG_CLUSTER
                    cerr << "This child has distances to end : " << chain_clusters.distance_start_left << " " << chain_clusters.distance_start_right 
                         << " " << chain_clusters.distance_end_left << " " << chain_clusters.distance_end_right << endl;
#endif
                    //And add it to its parent snarl
                    size_t parent_index;
                    if (tree_state.net_handle_to_index.count(parent) == 0) {
                        parent_index = tree_state.all_node_clusters.size();
                        tree_state.net_handle_to_index[parent] = parent_index;
                        tree_state.all_node_clusters.emplace_back(parent, tree_state.all_seeds->size(),
                                                                  tree_state.seed_count_prefix_sum.back(), distance_index);
                        if (chain_clusters.has_grandparent_handle) {
                            tree_state.all_node_clusters.back().has_parent_handle = true;
                            tree_state.all_node_clusters.back().parent_net_handle = chain_clusters.grandparent_net_handle;
                        }
                    } else {
                        parent_index = tree_state.net_handle_to_index[parent];
                    }
                    tree_state.snarl_to_children.emplace(parent_index, current_iterator->parent_index);


                }
            }
            chain_range_start = current_iterator;
            only_seeds = true;
        }

        if (current_iterator->child_index2 == std::numeric_limits<size_t>::max()) {
            //If this is a snarl
            only_seeds = false;
        }
    }
}


void SnarlDistanceIndexClusterer::cluster_one_node(
                   TreeState& tree_state, NodeClusters& node_clusters, 
                   const std::vector<std::tuple<size_t, size_t, size_t>>::iterator& seed_range_start,
                   const std::vector<std::tuple<size_t, size_t, size_t>>::iterator& seed_range_end) const {
#ifdef DEBUG_CLUSTER
    cerr << "Finding clusters on node " << distance_index.net_handle_as_string(node_clusters.containing_net_handle) << endl;
    for (auto x = seed_range_start ; x != seed_range_end ; ++x) {
        NodeClusters& seed_clusters = tree_state.all_node_clusters[std::get<0>(*x)];
        assert(&node_clusters == &seed_clusters);
        assert(node_clusters.containing_net_handle == seed_clusters.containing_net_handle);
    }
#endif

    size_t node_length = node_clusters.node_length;
    nid_t node_id = node_clusters.node_id;


    std::function<std::tuple<size_t, size_t, size_t>(const std::tuple<size_t, size_t, size_t>&)> get_offset_from_indices = 
        [&](const std::tuple<size_t, size_t, size_t>& seed_index){
            //This function returns a tuple of <read num, seed num, left offset>
            return std::make_tuple(std::get<1>(seed_index), std::get<2>(seed_index),
                    tree_state.all_seeds->at(std::get<1>(seed_index))->at(std::get<2>(seed_index)).distance_left); 
    };
    cluster_seeds_on_linear_structure(tree_state, node_clusters, seed_range_start, seed_range_end, 
                                      node_length, get_offset_from_indices, false);

#ifdef DEBUG_CLUSTER

    cerr << "\tFound read clusters on node " << node_id << endl;

    bool got_left = false;
    bool got_right = false;
    for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
        cerr << "\t for read num " << read_num << " best left: " << (read_num == 0 ? node_clusters.read_best_left.first : node_clusters.read_best_left.second) << " best right: " << (read_num == 0 ? node_clusters.read_best_right.first : node_clusters.read_best_right.second) << endl;
        bool got_read_left=false;
        bool got_read_right = false;
        for (pair<size_t,size_t> c : node_clusters.read_cluster_heads) {
            if (c.first == read_num) {
                pair<size_t, size_t> dists (tree_state.all_seeds->at(c.first)->at(c.second).distance_left, tree_state.all_seeds->at(c.first)->at(c.second).distance_right);
                cerr << "\t\t" << c.first << ":"<<c.second << ": left: " << dists.first << " right : " << dists.second << ": ";
                bool has_seeds = false;
                for (size_t x = 0 ; x < tree_state.all_seeds->at(c.first)->size() ; x++) {
                    if (tree_state.read_union_find[c.first].find_group(x) == c.second) {
                        cerr << tree_state.all_seeds->at(c.first)->at(x).pos << " ";
                        has_seeds = true;
                    }
                }
                assert(dists.first == std::numeric_limits<size_t>::max() || dists.first >= (read_num == 0 ? node_clusters.read_best_left.first : node_clusters.read_best_left.second));
 ;               assert(dists.second == std::numeric_limits<size_t>::max() || dists.second >= (read_num == 0 ? node_clusters.read_best_right.first : node_clusters.read_best_right.second));
                assert(dists.first == std::numeric_limits<size_t>::max() || dists.first >= node_clusters.fragment_best_left);
                assert(dists.second == std::numeric_limits<size_t>::max() || dists.second >= node_clusters.fragment_best_right);
                if (dists.first == node_clusters.fragment_best_left) {got_left = true;}
                if (dists.second == node_clusters.fragment_best_right) {got_right = true;}
                if (dists.first == (read_num == 0 ? node_clusters.read_best_left.first : node_clusters.read_best_left.second)) {got_read_left = true;}
                if (dists.second == (read_num == 0 ? node_clusters.read_best_right.first : node_clusters.read_best_right.second)) {got_read_right = true;}
                cerr << endl;
                assert(has_seeds);
            }
        }
        //assert(got_read_left || (read_num == 0 ? node_clusters.read_best_left.first : node_clusters.read_best_left.second) == std::numeric_limits<size_t>::max());
        //assert(got_read_right || (read_num == 0 ? node_clusters.read_best_right.first : node_clusters.read_best_right.second) == std::numeric_limits<size_t>::max());
    }
    //assert(got_left);
    //assert(got_right);
    for (pair<size_t, size_t> group_id : node_clusters.read_cluster_heads) {
        assert (group_id.second == tree_state.read_union_find[group_id.first].find_group(group_id.second));
    }
#endif
    return;
        
};


//Go through pairs of clusters of the two children and see which ones can be combined
//The first child may not have been seen before, so all of it's clusters may be added to the parent, then
//anything that was combined gets removed and only the cluster heads get added.
//For the second child, everything is already in the parent so remove ones that were combined then
//add the head of the combined clusters
//
//If this is the first time we see the first child, then also update the best distances to the ends of the 
//parent for the parent clusters
void SnarlDistanceIndexClusterer::compare_and_combine_cluster_on_child_structures(TreeState& tree_state, NodeClusters& child_clusters1, 
    NodeClusters& child_clusters2, NodeClusters& parent_clusters, 
    const vector<pair<size_t, size_t>> & child_distances, bool is_root, bool first_child) const {
#ifdef DEBUG_CLUSTER
    cerr << "\tCompare " << distance_index.net_handle_as_string(child_clusters1.containing_net_handle) 
         << " and " << distance_index.net_handle_as_string(child_clusters2.containing_net_handle)
         << " which are children of " << distance_index.net_handle_as_string(parent_clusters.containing_net_handle) << endl;
#endif

    net_handle_t& parent_handle = parent_clusters.containing_net_handle;
    net_handle_t& child_handle1 = child_clusters1.containing_net_handle;
    net_handle_t& child_handle2 = child_clusters2.containing_net_handle;



    //Get the distances between the two sides of the children in the parent
    size_t distance_left_left = distance_index.distance_in_parent(parent_handle, distance_index.flip(child_handle1), 
                                            distance_index.flip(child_handle2), graph,
                                            (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit 
                                                                                     : tree_state.fragment_distance_limit));
    size_t distance_left_right = distance_index.distance_in_parent(parent_handle, distance_index.flip(child_handle1), 
                                            child_handle2, graph,
                                            (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit 
                                                                                     : tree_state.fragment_distance_limit));
    size_t distance_right_right = distance_index.distance_in_parent(parent_handle, child_handle1, child_handle2, graph,
                                            (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit 
                                                                                     : tree_state.fragment_distance_limit));
    size_t distance_right_left = distance_index.distance_in_parent(parent_handle, child_handle1, 
                                            distance_index.flip(child_handle2), graph,
                                            (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit 
                                                                                     : tree_state.fragment_distance_limit));

    //If it's the root and nothing can be clustered, return here
    //If it's not the root, then we still need add the new clusters
    if (is_root){
        if (distance_left_left == std::numeric_limits<size_t>::max() &&
            distance_left_right == std::numeric_limits<size_t>::max() &&
            distance_right_right == std::numeric_limits<size_t>::max() &&
            distance_right_left == std::numeric_limits<size_t>::max()) {
            return;
        }
    }
 
   
#ifdef DEBUG_CLUSTER
    cerr << "\t\tFound distances between the two children: " << distance_left_left << " " << distance_left_right << " " << distance_right_right << " " << distance_right_left << endl;
    cerr << "\t\tBest left and right distances for the two children: " << child_clusters1.fragment_best_left << " " << child_clusters1.fragment_best_right << " and " << child_clusters2.fragment_best_left << " " << child_clusters2.fragment_best_right << endl; 
    cerr << "\t\tAnd distances from the ends of child1 to ends of parent: " << child_clusters1.distance_start_left << " " 
         << child_clusters1.distance_start_right << " " << child_clusters1.distance_end_left << " " << child_clusters1.distance_end_right << endl;
#endif
    /*
     * We're going to go through all clusters to see which can get combined. There will be up to four combined clusters (per read), 
     * one for each path between the two sides of the two nodes
     *
     * If a cluster from the first child can be combined with a cluster of the second by taking the left-right path, 
     * then any cluster of the second child with a right distance that is less than the distance limit - the best left distance 
     * of the first will be combined with the cluster with the best left distances. So in the end there will only be four combined
     * clusters, one for each path
     * This strategy ends up unioning a cluster with itself but it only goes through the clusters once so I think
     * it's efficient
     */

    //The cluster heads that will be removed from the parent's read_cluster_heads
    vector<pair<size_t, size_t>> to_erase;

    //Helper function that will compare two clusters
    //Given the read num and seed_num of the cluster head, the distance to the other node side we're looking at, 
    //the distances to the ends of the parent for the cluster head, a reference
    //to the current cluster head and distances of the potential combined cluster (pair<pair<>pair<>> which will be updated if it gets combined),
    //the relevant combined cluster head for the fragment
    //Returns true if this cluster got combined
    auto compare_and_combine_clusters = [&] (const size_t& read_num, const size_t& cluster_num, const size_t& distance_between_reads, 
            const size_t& distance_between_fragments, pair<size_t, size_t>& old_distances, 
            ClusterIndices& new_cluster_head_and_distances, size_t& new_cluster_head_fragment){
        if ((read_num == new_cluster_head_and_distances.read_num 
                && cluster_num ==  new_cluster_head_and_distances.cluster_num) || 
           ( distance_between_fragments == std::numeric_limits<size_t>::max())) {
            //If this is the same as the old cluster head, or the distances are infinite,
            //then don't bother trying to compare
            return false;
        }
        size_t distance_reads = SnarlDistanceIndex::minus(distance_between_reads, 1);
        size_t distance_fragments = SnarlDistanceIndex::minus(distance_between_fragments, 1);
        bool combined = false;

        if (distance_reads <= tree_state.read_distance_limit) {
            //If this can be combined with the given combined cluster
            if (new_cluster_head_and_distances.read_num == std::numeric_limits<size_t>::max()){
                //new cluster head
                new_cluster_head_and_distances.read_num =read_num;
                new_cluster_head_and_distances.cluster_num = cluster_num;
                new_cluster_head_and_distances.distance_left = old_distances.first;
                new_cluster_head_and_distances.distance_right = old_distances.second;
            } else {
                //Combine with old cluster head
                size_t new_cluster_head = tree_state.read_union_find.at(read_num).union_groups(cluster_num, new_cluster_head_and_distances.cluster_num);

                //Update distances
                size_t new_best_left = std::min(old_distances.first, new_cluster_head_and_distances.distance_left);
                size_t new_best_right = std::min(old_distances.second, new_cluster_head_and_distances.distance_right);

                //And remember new head and distances
                new_cluster_head_and_distances.read_num =read_num;
                new_cluster_head_and_distances.cluster_num = new_cluster_head;
                new_cluster_head_and_distances.distance_left = new_best_left; 
                new_cluster_head_and_distances.distance_right = new_best_right;

                //Remember these distances because we might need to check them later
                old_distances = make_pair(new_best_left, new_best_right);
                tree_state.all_seeds->at(read_num)->at(cluster_num).distance_left = new_best_left;
                tree_state.all_seeds->at(read_num)->at(cluster_num).distance_right = new_best_right;
            }
            //Remember to erase this cluster head
            to_erase.emplace_back(read_num, cluster_num);
            combined = true;

#ifdef DEBUG_CLUSTER
            cerr << "\t\t\tCombining read/cluster " << read_num << "/" << cluster_num << "... new cluster head:" << tree_state.all_seeds->at(read_num)->at(new_cluster_head_and_distances.cluster_num).pos << endl; 
            cerr << "\t\t\t\t Best distances for this cluster: " << old_distances.first << " and " << old_distances.second << endl;
            cerr << "\t\t\t\t New best distances for combined cluster: " << new_cluster_head_and_distances.distance_left << " and " << new_cluster_head_and_distances.distance_right << endl;
#endif
        }
        if (tree_state.fragment_distance_limit != 0 && 
                    distance_fragments <= tree_state.fragment_distance_limit ) {
            //Just union the fragment
            if (new_cluster_head_fragment == std::numeric_limits<size_t>::max()) {
                new_cluster_head_fragment =cluster_num+tree_state.seed_count_prefix_sum[read_num];
            } else {
                new_cluster_head_fragment = tree_state.fragment_union_find.union_groups(cluster_num+tree_state.seed_count_prefix_sum[read_num], 
                                                                    new_cluster_head_fragment);
            }
#ifdef DEBUG_CLUSTER
            cerr << "\t\t\tCombining fragment" << endl;
#endif
        }
        return combined;
    };
    /*
     * Go through all clusters on the first child and see if they can be combined with clusters on the second child
     */

    //Did any cluster get combined? If it didn't, we don't need to go through the second node's clusters
    bool combined_anything = false;


    if (first_child || distance_left_left != std::numeric_limits<size_t>::max()
                    || distance_left_right != std::numeric_limits<size_t>::max()
                    || distance_right_left != std::numeric_limits<size_t>::max()
                    || distance_right_right != std::numeric_limits<size_t>::max()){
        //These will be the cluster heads and distances of everything combined by taking the indicated path
        //one cluster head per read
        //The default value will be ((inf, 0), (0,0)). Only the inf gets checked to see if it's a real value so I filled it in with 0 so I wouldn't have to type out inf 
        pair<ClusterIndices, ClusterIndices> new_cluster_left_left_by_read;
        pair<ClusterIndices, ClusterIndices> new_cluster_left_right_by_read;
        pair<ClusterIndices, ClusterIndices> new_cluster_right_right_by_read;
        pair<ClusterIndices, ClusterIndices> new_cluster_right_left_by_read;

        //And the new cluster heads for the fragment
        //These are the values of the cluster heads in the union finds, which include the values from read_index_offset
        size_t new_cluster_left_left_fragment = std::numeric_limits<size_t>::max();
        size_t new_cluster_left_right_fragment = std::numeric_limits<size_t>::max();
        size_t new_cluster_right_right_fragment = std::numeric_limits<size_t>::max();
        size_t new_cluster_right_left_fragment = std::numeric_limits<size_t>::max();


        for (auto& child_cluster_head : child_clusters1.read_cluster_heads) {

            bool combined = false;
            size_t read_num = child_cluster_head.first;
            ClusterIndices& new_cluster_left_left = (read_num == 0 ? new_cluster_left_left_by_read.first : new_cluster_left_left_by_read.second);
            ClusterIndices& new_cluster_left_right = (read_num == 0 ? new_cluster_left_right_by_read.first : new_cluster_left_right_by_read.second);
            ClusterIndices& new_cluster_right_right = (read_num == 0 ? new_cluster_right_right_by_read.first : new_cluster_right_right_by_read.second);
            ClusterIndices& new_cluster_right_left = (read_num == 0 ? new_cluster_right_left_by_read.first : new_cluster_right_left_by_read.second);
            size_t cluster_num = tree_state.read_union_find[read_num].find_group(child_cluster_head.second);

            //Distances to the ends of the child
            pair<size_t, size_t> distances = child_distances[child_cluster_head.second + tree_state.seed_count_prefix_sum[read_num]];

            //Distances to the parent
            size_t new_dist_left = std::min(SnarlDistanceIndex::sum(distances.first, child_clusters1.distance_start_left),
                                            SnarlDistanceIndex::sum(distances.second, child_clusters1.distance_start_right));
            size_t new_dist_right= std::min(SnarlDistanceIndex::sum(distances.first, child_clusters1.distance_end_left),
                                            SnarlDistanceIndex::sum(distances.second, child_clusters1.distance_end_right));

            pair<size_t, size_t> distances_to_parent = make_pair(new_dist_left, new_dist_right);
            //If this is already in the parent, take the minimum of the parent distances
            if (parent_clusters.read_cluster_heads.count(make_pair(read_num, cluster_num)) > 0) {
                distances_to_parent = make_pair(
                    std::min(new_dist_left, tree_state.all_seeds->at(read_num)->at(cluster_num).distance_left),
                    std::min(new_dist_right, tree_state.all_seeds->at(read_num)->at(cluster_num).distance_right));
            }
                

            //Check if the left of 1 can connect with the left of 2
            combined = combined | compare_and_combine_clusters (read_num, cluster_num, 
                SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(distances.first, distance_left_left), (read_num == 0 ? child_clusters2.read_best_left.first : child_clusters2.read_best_left.second)), 
                SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(distances.first, distance_left_left), child_clusters2.fragment_best_left), 
                distances_to_parent, new_cluster_left_left, 
                new_cluster_left_left_fragment);

            //Check if the left of 1 can connect with the right of 2
            combined = combined | compare_and_combine_clusters (read_num, cluster_num, 
                SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(distances.first,distance_left_right), (read_num == 0 ? child_clusters2.read_best_right.first : child_clusters2.read_best_right.second)), 
                SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(distances.first,distance_left_right), child_clusters2.fragment_best_right), 
                distances_to_parent, new_cluster_left_right, 
                new_cluster_left_right_fragment);

            //Check if the right of 1 can connect with the right of 2
            combined = combined | compare_and_combine_clusters (read_num, cluster_num, 
                SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(distances.second,distance_right_right), (read_num == 0 ? child_clusters2.read_best_right.first : child_clusters2.read_best_right.second)), 
                SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(distances.second,distance_right_right), child_clusters2.fragment_best_right), 
                distances_to_parent, new_cluster_right_right, 
                new_cluster_right_right_fragment);

            //Check if the right of 1 can connect with the left of 2
            combined = combined | compare_and_combine_clusters (read_num, cluster_num, 
                SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(distances.second,distance_right_left), (read_num == 0 ? child_clusters2.read_best_left.first : child_clusters2.read_best_left.second)), 
                SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(distances.second,distance_right_left), child_clusters2.fragment_best_left), 
                distances_to_parent, new_cluster_right_left, 
                new_cluster_right_left_fragment);

            //Is the distance small enough that we can cluster it with something else?
            bool reachable_left = distances_to_parent.first <= 
                (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit);
            bool reachable_right = distances_to_parent.second <= 
                (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit);
            //If this cluster wasn't combined and hasn't been seen before and its reachable from other clusters, add it to the parent
            if (first_child && !combined && (reachable_left || reachable_right)) {
                parent_clusters.read_cluster_heads.emplace(read_num, cluster_num);
                tree_state.all_seeds->at(read_num)->at(cluster_num).distance_left = distances_to_parent.first;
                tree_state.all_seeds->at(read_num)->at(cluster_num).distance_right = distances_to_parent.second;
            }
            combined_anything |= combined;
        }

        if (combined_anything || new_cluster_left_left_fragment != std::numeric_limits<size_t>::max()
                              || new_cluster_left_right_fragment != std::numeric_limits<size_t>::max()
                              || new_cluster_right_left_fragment != std::numeric_limits<size_t>::max()
                              || new_cluster_right_right_fragment != std::numeric_limits<size_t>::max()) {
            /* If anything got combined, then we have to go through the second child
             * and see if any cluster can be combined with clusters on the first child
             */
            for (auto& child_cluster_head : child_clusters2.read_cluster_heads) {

                size_t read_num = child_cluster_head.first;
                size_t cluster_num = tree_state.read_union_find[read_num].find_group(child_cluster_head.second);

                ClusterIndices& new_cluster_left_left = (read_num == 0 ? new_cluster_left_left_by_read.first 
                                                                       : new_cluster_left_left_by_read.second);
                ClusterIndices& new_cluster_left_right = (read_num == 0 ? new_cluster_left_right_by_read.first 
                                                                        : new_cluster_left_right_by_read.second);
                ClusterIndices& new_cluster_right_right = (read_num == 0 ? new_cluster_right_right_by_read.first 
                                                                         : new_cluster_right_right_by_read.second);
                ClusterIndices& new_cluster_right_left = (read_num == 0 ? new_cluster_right_left_by_read.first 
                                                                        : new_cluster_right_left_by_read.second);

                pair<size_t, size_t> distances = child_distances[child_cluster_head.second 
                                                    + tree_state.seed_count_prefix_sum[read_num]];
                size_t new_dist_left = std::min(SnarlDistanceIndex::sum(distances.first,child_clusters2.distance_start_left), 
                                                SnarlDistanceIndex::sum(distances.second,child_clusters2.distance_start_right));
                size_t new_dist_right = std::min(SnarlDistanceIndex::sum(distances.first,child_clusters2.distance_end_left), 
                                                SnarlDistanceIndex::sum(distances.second,child_clusters2.distance_end_right));
                pair<size_t, size_t> distances_to_parent = make_pair(new_dist_left, new_dist_right);

                if (parent_clusters.read_cluster_heads.count(make_pair(read_num, cluster_num)) > 0) {
                    distances_to_parent = make_pair(
                        std::min(new_dist_left, tree_state.all_seeds->at(read_num)->at(cluster_num).distance_left),
                        std::min(new_dist_right, tree_state.all_seeds->at(read_num)->at(cluster_num).distance_right));
                }

                //Check if the left of 1 can connect with the left of 2
                compare_and_combine_clusters (read_num, cluster_num, 
                    SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(distances.first,distance_left_left), (read_num == 0 ? child_clusters1.read_best_left.first : child_clusters1.read_best_left.second)), 
                    SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(distances.first,distance_left_left),child_clusters1.fragment_best_left), 
                    distances_to_parent, new_cluster_left_left, new_cluster_left_left_fragment);

                //Check if the left of 1 can connect with the right of 2
                compare_and_combine_clusters (read_num, cluster_num, 
                    SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(distances.second,distance_left_right), (read_num == 0 ? child_clusters1.read_best_left.first : child_clusters1.read_best_left.second)),
                    SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(distances.second,distance_left_right),child_clusters1.fragment_best_left),
                    distances_to_parent, new_cluster_left_right, new_cluster_left_right_fragment);

                //Check if the right of 1 can connect with the right of 2
                compare_and_combine_clusters (read_num, cluster_num, 
                    SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(distances.second,distance_right_right), (read_num == 0 ? child_clusters1.read_best_right.first : child_clusters1.read_best_right.second)),
                    SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(distances.second,distance_right_right),child_clusters1.fragment_best_right),
                    distances_to_parent, new_cluster_right_right, new_cluster_right_right_fragment);

                //Check if the right of 1 can connect with the left of 2
                compare_and_combine_clusters (read_num, cluster_num, 
                    SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(distances.first,distance_right_left), (read_num == 0 ? child_clusters1.read_best_right.first : child_clusters1.read_best_right.second)),
                    SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(distances.first,distance_right_left),child_clusters1.fragment_best_right),
                    distances_to_parent, new_cluster_right_left, new_cluster_right_left_fragment);
            }

            /*then remove all clusters that got erase, then add back in the cluster heads
             */

            //remove cluster heads that got combined with new ones
            for (pair<size_t, size_t>& cluster_head : to_erase) {
                parent_clusters.read_cluster_heads.erase(cluster_head);
            }

            //And add back in the new cluster heads
            for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
                ClusterIndices& new_cluster_left_left = (read_num == 0 ? new_cluster_left_left_by_read.first 
                                                                       : new_cluster_left_left_by_read.second);
                ClusterIndices& new_cluster_left_right = (read_num == 0 ? new_cluster_left_right_by_read.first 
                                                                        : new_cluster_left_right_by_read.second);
                ClusterIndices& new_cluster_right_right = (read_num == 0 ? new_cluster_right_right_by_read.first 
                                                                         : new_cluster_right_right_by_read.second);
                ClusterIndices& new_cluster_right_left = (read_num == 0 ? new_cluster_right_left_by_read.first 
                                                                        : new_cluster_right_left_by_read.second);

                //If the new cluster is clusterable, then add the new cluster_left_left
                if (new_cluster_left_left.read_num != std::numeric_limits<size_t>::max()){

                    //Check the old distances for this cluster head, because it may include the distance from another combined cluster
                    //head
                    pair<size_t, size_t> old_distances = 
                                                parent_clusters.read_cluster_heads.count(make_pair(new_cluster_left_left.read_num, new_cluster_left_left.cluster_num)) == 0
                                                ? make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max())
                                                : make_pair(tree_state.all_seeds->at(read_num)->at(new_cluster_left_left.cluster_num).distance_left,
                                                            tree_state.all_seeds->at(read_num)->at(new_cluster_left_left.cluster_num).distance_right);
                     //Is the distance small enough that we can cluster it with something else?
                     size_t best_left = std::min(new_cluster_left_left.distance_left, old_distances.first);
                     size_t best_right = std::min(new_cluster_left_left.distance_right, old_distances.second);
                     bool reachable_left = best_left <= 
                                           (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit);
                     bool reachable_right = best_right <= 
                                            (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit);
                    if ((reachable_left || reachable_right)) {
                        parent_clusters.read_cluster_heads.emplace(new_cluster_left_left.read_num, new_cluster_left_left.cluster_num);
                        tree_state.all_seeds->at(read_num)->at(new_cluster_left_left.cluster_num).distance_left = best_left;
                        tree_state.all_seeds->at(read_num)->at(new_cluster_left_left.cluster_num).distance_right = best_right;
                    } else {
                        parent_clusters.read_cluster_heads.erase(std::make_pair( new_cluster_left_left.read_num, new_cluster_left_left.cluster_num ));
                    }
                }
                //Add the new cluster_right_right
                if (new_cluster_right_right.read_num != std::numeric_limits<size_t>::max()){


                    pair<size_t, size_t> old_distances = parent_clusters.read_cluster_heads.count(make_pair(new_cluster_right_right.read_num, new_cluster_right_right.cluster_num)) == 0
                        ? make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max())
                        : make_pair(tree_state.all_seeds->at(read_num)->at( new_cluster_right_right.cluster_num).distance_left,
                                     tree_state.all_seeds->at(read_num)->at( new_cluster_right_right.cluster_num).distance_right);
                    size_t best_left =  std::min(new_cluster_right_right.distance_left, 
                                                 old_distances.first);
                    size_t best_right = std::min( new_cluster_right_right.distance_right, 
                                                  old_distances.second);

                    bool reachable_left = best_left <= 
                        (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit);
                    bool reachable_right = best_right <= 
                        (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit);

                    if ((reachable_left || reachable_right)){
                        parent_clusters.read_cluster_heads.emplace( new_cluster_right_right.read_num, new_cluster_right_right.cluster_num);
                        tree_state.all_seeds->at(read_num)->at(new_cluster_right_right.cluster_num).distance_left = best_left;
                        tree_state.all_seeds->at(read_num)->at(new_cluster_right_right.cluster_num).distance_right = best_right;
                    } else {
                        parent_clusters.read_cluster_heads.erase(make_pair(new_cluster_right_right.read_num, new_cluster_right_right.cluster_num));
                    }
                }

                //Add the new cluster_left_right
                if (new_cluster_left_right.read_num != std::numeric_limits<size_t>::max()) {
                    pair<size_t, size_t> old_distances = parent_clusters.read_cluster_heads.count(
                                make_pair(new_cluster_left_right.read_num, new_cluster_left_right.cluster_num)) == 0
                           ? make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max())
                           : make_pair(tree_state.all_seeds->at(read_num)->at(new_cluster_left_right.cluster_num).distance_left,
                    tree_state.all_seeds->at(read_num)->at(new_cluster_left_right.cluster_num).distance_right);

                    size_t best_left =  std::min(new_cluster_left_right.distance_left, old_distances.first);
                    size_t best_right = std::min( new_cluster_left_right.distance_right, old_distances.second);

                    bool reachable_left = best_left <= 
                        (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit);
                    bool reachable_right = best_right <= 
                        (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit);

                    if ((reachable_left || reachable_right) ){

                        parent_clusters.read_cluster_heads.emplace(  new_cluster_left_right.read_num, new_cluster_left_right.cluster_num);
                        tree_state.all_seeds->at(read_num)->at( new_cluster_left_right.cluster_num).distance_left = best_left;
                        tree_state.all_seeds->at(read_num)->at(new_cluster_left_right.cluster_num).distance_right = best_right;
                    } else {
                        parent_clusters.read_cluster_heads.erase(make_pair(new_cluster_left_right.read_num, new_cluster_left_right.cluster_num));
                    }
                }
                //Add the new cluster_right_left
                if ( new_cluster_right_left.read_num != std::numeric_limits<size_t>::max()) {

                    pair<size_t, size_t> old_distances = parent_clusters.read_cluster_heads.count(make_pair(new_cluster_right_left.read_num, new_cluster_right_left.cluster_num)) == 0
                        ? make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max())
                        : make_pair(tree_state.all_seeds->at(read_num)->at(new_cluster_right_left.cluster_num).distance_left,
                                    tree_state.all_seeds->at(read_num)->at(new_cluster_right_left.cluster_num).distance_right);
                    size_t best_left = std::min(new_cluster_right_left.distance_left, old_distances.first);
                    size_t best_right = std::min(new_cluster_right_left.distance_right, old_distances.second);

                    bool reachable_left = best_left <= 
                        (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit);
                    bool reachable_right = best_right <= 
                        (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit);

                    if ((reachable_left || reachable_right)) { 
                        parent_clusters.read_cluster_heads.emplace(new_cluster_right_left.read_num, new_cluster_right_left.cluster_num);
                        tree_state.all_seeds->at(read_num)->at(new_cluster_right_left.cluster_num).distance_left = best_left;
                        tree_state.all_seeds->at(read_num)->at(new_cluster_right_left.cluster_num).distance_right = best_right;
                    } else {
                        parent_clusters.read_cluster_heads.erase(make_pair(new_cluster_right_left.read_num, new_cluster_right_left.cluster_num));
                    }
                }
            }
        }
    }


    /*Update the parent's best left and right distances, only looking at the first child since we've already seen the second one
     */
    if (first_child) {
        //Update the parent's fragment best distances
        parent_clusters.fragment_best_left = std::min(parent_clusters.fragment_best_left,
                                             std::min(SnarlDistanceIndex::sum(child_clusters1.distance_start_left, child_clusters1.fragment_best_left),
                                                      SnarlDistanceIndex::sum(child_clusters1.distance_start_right , child_clusters1.fragment_best_right)));
        parent_clusters.fragment_best_right = std::min(parent_clusters.fragment_best_right,
                                              std::min(SnarlDistanceIndex::sum(child_clusters1.distance_end_left , child_clusters1.fragment_best_left),
                                                       SnarlDistanceIndex::sum(child_clusters1.distance_end_right , child_clusters1.fragment_best_right)));


        //Update the best distances in the parent for each read num
        for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num ++) {
            //Find the best distances to the ends of the parent from child1
            size_t best_start = std::min(SnarlDistanceIndex::sum(child_clusters1.distance_start_left, (read_num == 0 ? child_clusters1.read_best_left.first : child_clusters1.read_best_left.second)), 
                                        SnarlDistanceIndex::sum(child_clusters1.distance_start_right, (read_num == 0 ? child_clusters1.read_best_right.first : child_clusters1.read_best_right.second)));
            size_t best_end = std::min(SnarlDistanceIndex::sum(child_clusters1.distance_end_left, (read_num == 0 ? child_clusters1.read_best_left.first : child_clusters1.read_best_left.second)), 
                                        SnarlDistanceIndex::sum(child_clusters1.distance_end_right, (read_num == 0 ? child_clusters1.read_best_right.first : child_clusters1.read_best_right.second)));
            //And update the distances in the parent
            if (read_num == 0) {
                parent_clusters.read_best_left.first = std::min(best_start, parent_clusters.read_best_left.first);
                parent_clusters.read_best_right.first = std::min(best_end, parent_clusters.read_best_right.first);
            } else {
                parent_clusters.read_best_left.second = std::min(best_start, parent_clusters.read_best_left.second);
                parent_clusters.read_best_right.second = std::min(best_end, parent_clusters.read_best_right.second);
            }
        }
    }
    
//#ifdef DEBUG_CLUSTER
//        cerr << "\tIntermediate clusters on " << distance_index.net_handle_as_string(parent_clusters.containing_net_handle);
//        cerr << "   with best left and right values: " << parent_clusters.fragment_best_left << " "
//             << parent_clusters.fragment_best_right << endl;
//        bool got_left = false;
//        bool got_right = false;
//        for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
//            cerr << "\t\t\tfor read num " << read_num << " best left: " << (read_num == 0 ? parent_clusters.read_best_left.first : parent_clusters.read_best_lef.second) << " best right: " << (read_num == 0 ? parent_clusters.read_best_right.first : parent_clusters.read_best_right.second) << endl;
//            for (pair<pair<size_t, size_t>, pair<size_t,size_t>> c : parent_clusters.read_cluster_heads) {
//                if (c.first.first == read_num) {
//                    pair<size_t, size_t> dists = c.second;
//                    cerr << "\t\t\t" << tree_state.all_seeds->at(c.first.first)->at(c.first.second).pos << " (" << c.first.first << ":"<<c.first.second << ")  left: " << dists.first << " right : " << dists.second << ": ";
//                    for (size_t x = 0 ; x < tree_state.all_seeds->at(c.first.first)->size() ; x++) {
//                        if (tree_state.read_union_find[c.first.first].find_group(x) == c.first.second) {
//                            cerr << tree_state.all_seeds->at(c.first.first)->at(x).pos << " ";
//                        }
//                    }
//                }
//                cerr << endl;
//            }
//        }
//#endif
}

void SnarlDistanceIndexClusterer::compare_and_combine_cluster_on_one_child(TreeState& tree_state, NodeClusters& child_clusters) const {
#ifdef DEBUG_CLUSTER
    cerr << "\tCompare " << distance_index.net_handle_as_string(child_clusters.containing_net_handle) 
         << " to itself in the root" << endl;
#endif

    net_handle_t& handle = child_clusters.containing_net_handle;


    //Get the distances between the two sides of the child
    size_t distance_left_left = distance_index.is_externally_start_start_connected(handle) ? 0 : std::numeric_limits<size_t>::max();
    size_t distance_left_right = distance_index.is_externally_start_end_connected(handle) ? 0 : std::numeric_limits<size_t>::max();
    size_t distance_right_right = distance_index.is_externally_end_end_connected(handle) ? 0 : std::numeric_limits<size_t>::max();
    if (distance_left_left == std::numeric_limits<size_t>::max() &&
        distance_left_right == std::numeric_limits<size_t>::max() &&
        distance_right_right == std::numeric_limits<size_t>::max()) {
        //If there is no external connectivity
        return;
    }

#ifdef DEBUG_CLUSTER
    cerr << "\t\tFound distances between the two children: " << distance_left_left << " " << distance_left_right << " " << distance_right_right <<  endl;
#endif
    /*
     * We're going to go through all clusters to see which can get combined. There are at most three new clusters,
     * because there are three possible paths between the two ends of the single node
     */

    //These will be the cluster heads of everything combined by taking the indicated path
    //one cluster head per read
    //pair< pair<read_num, seed_num>>
    //The default value will be (inf, 0). Only the inf gets checked to see if it's a real value so I filled it in with 0 so I wouldn't have to type out inf 
    pair<pair<size_t, size_t>, pair<size_t, size_t>> new_cluster_left_left_by_read = std::make_pair(
            std::make_pair(std::numeric_limits<size_t>::max(), 0), std::make_pair(std::numeric_limits<size_t>::max(), 0));
    pair<pair<size_t, size_t>, pair<size_t, size_t>> new_cluster_left_right_by_read = std::make_pair(
            std::make_pair(std::numeric_limits<size_t>::max(), 0), std::make_pair(std::numeric_limits<size_t>::max(), 0));
    pair<pair<size_t, size_t>, pair<size_t, size_t>> new_cluster_right_right_by_read = std::make_pair(
            std::make_pair(std::numeric_limits<size_t>::max(), 0), std::make_pair(std::numeric_limits<size_t>::max(), 0));

    //And the new cluster heads for the fragment
    //These are the values of the cluster heads in the union finds, which include the values from read_index_offset
    size_t new_cluster_left_left_fragment = std::numeric_limits<size_t>::max();
    size_t new_cluster_left_right_fragment = std::numeric_limits<size_t>::max();
    size_t new_cluster_right_right_fragment = std::numeric_limits<size_t>::max();

    //Helper function that will compare two clusters
    //Given the read num and seed_num of the cluster head, the distance to the other node side we're looking at, 
    //the distances to the ends of the parent for the cluster head, a reference
    //to the current cluster head and distances of the potential combined cluster (pair<pair<>pair<>> which will be updated if it gets combined),
    //the relevant combined cluster head for the fragment
    //Returns true if this cluster got combined
    auto compare_and_combine_clusters = [&] (size_t read_num, size_t cluster_num, size_t distance_between_reads, 
            size_t distance_between_fragments, 
            pair<size_t, size_t>& new_cluster_head, size_t& new_cluster_head_fragment){

        if (read_num == new_cluster_head.first && cluster_num ==  new_cluster_head.second) {
            //If this is the same as the old cluster head, then don't bother trying to compare
            return;
        }
        distance_between_reads = SnarlDistanceIndex::minus(distance_between_reads, 1);
        distance_between_fragments = SnarlDistanceIndex::minus(distance_between_fragments, 1);

        if (distance_between_reads <= tree_state.read_distance_limit) {
            //If this can be combined with the given combined cluster
            if (new_cluster_head.first == std::numeric_limits<size_t>::max()){
                //new cluster head
                new_cluster_head = make_pair(read_num, cluster_num);
            } else {
                //Combine with old cluster head
                new_cluster_head = make_pair(read_num, 
                                             tree_state.read_union_find.at(read_num).union_groups(cluster_num, new_cluster_head.second));
            }
            //Remember to erase this cluster head

#ifdef DEBUG_CLUSTER
            cerr << "\t\t\tCombining read/cluster " << read_num << "/" << cluster_num << "... new cluster head:" << new_cluster_head.second << endl; 
#endif
        }
        if (tree_state.fragment_distance_limit != 0 && 
                    distance_between_fragments <= tree_state.fragment_distance_limit ) {
            //Just union the fragment
            if (new_cluster_head_fragment == std::numeric_limits<size_t>::max()) {
                new_cluster_head_fragment =cluster_num+tree_state.seed_count_prefix_sum[read_num];
            } else {
                new_cluster_head_fragment = tree_state.fragment_union_find.union_groups(cluster_num+tree_state.seed_count_prefix_sum[read_num], 
                                                                    new_cluster_head_fragment);
            }
#ifdef DEBUG_CLUSTER
            cerr << "\t\t\tCombining fragment" << endl;
#endif
        }
        return;
    };

    /*
     * Go through all clusters and see if they can be combined with anything
     */
    for (auto& child_cluster_head : child_clusters.read_cluster_heads) {

        size_t read_num = child_cluster_head.first;
        size_t cluster_num = tree_state.read_union_find[read_num].find_group(child_cluster_head.second);

        pair<size_t, size_t>& new_cluster_left_left = read_num == 0 ? new_cluster_left_left_by_read.first 
                                                              : new_cluster_left_left_by_read.second;
        pair<size_t, size_t>& new_cluster_left_right = read_num == 0 ? new_cluster_left_right_by_read.first 
                                                               : new_cluster_left_right_by_read.second;
        pair<size_t, size_t>& new_cluster_right_right = read_num == 0 ? new_cluster_right_right_by_read.first 
                                                                : new_cluster_right_right_by_read.second;

        //Distances to the ends of the child
        pair<size_t, size_t> distances (tree_state.all_seeds->at(read_num)->at(child_cluster_head.second).distance_left,
                                        tree_state.all_seeds->at(read_num)->at(child_cluster_head.second).distance_right);
            

        //Check if this can be combined left-left
        compare_and_combine_clusters (read_num, cluster_num, 
            SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(distances.first,distance_left_left), (read_num == 0 ? child_clusters.read_best_left.first : child_clusters.read_best_left.second)), 
            SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(distances.first,distance_left_left), child_clusters.fragment_best_left), 
                new_cluster_left_left, new_cluster_left_left_fragment);

        //Check if the left of this can be combined with the right of anything else
        compare_and_combine_clusters (read_num, cluster_num, 
            SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(distances.first,distance_left_right), (read_num == 0 ? child_clusters.read_best_right.first : child_clusters.read_best_right.second)), 
            SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(distances.first,distance_left_right), child_clusters.fragment_best_right), 
              new_cluster_left_right, new_cluster_left_right_fragment);

        //Check if the right of this can be combined with the left of anything else
        compare_and_combine_clusters (read_num, cluster_num, 
            SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(distances.second,distance_left_right), (read_num == 0 ? child_clusters.read_best_left.first : child_clusters.read_best_left.second)), 
            SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(distances.second,distance_left_right), child_clusters.fragment_best_left), 
              new_cluster_left_right, new_cluster_left_right_fragment);

        //Check if this can be combined right-right
        compare_and_combine_clusters (read_num, cluster_num, 
            SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(distances.second,distance_right_right), (read_num == 0 ? child_clusters.read_best_right.first : child_clusters.read_best_right.second)), 
            SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(distances.second,distance_right_right), child_clusters.fragment_best_right), 
            new_cluster_right_right, new_cluster_right_right_fragment);
    }

}


void SnarlDistanceIndexClusterer::cluster_one_snarl(TreeState& tree_state, size_t snarl_clusters_index, 
        std::multimap<size_t, size_t>::iterator child_range_start, std::multimap<size_t, size_t>::iterator child_range_end) const { 
    //Get the clusters on this snarl, assumes that all of the snarls children have been clustered already.
    

    NodeClusters& snarl_clusters = tree_state.all_node_clusters[snarl_clusters_index]; 
    snarl_clusters.set_snarl_values(distance_index);
    net_handle_t& snarl_handle = snarl_clusters.containing_net_handle;

#ifdef DEBUG_CLUSTER
        cerr << "Finding clusters on snarl " << distance_index.net_handle_as_string(snarl_handle) << endl;
#endif


    //If the snarl is a simple snarl, then there is no clustering to do because there is no path between
    //the nodes. Otherwise, compare the children of the snarl
    if (!distance_index.is_simple_snarl(snarl_handle)) {
        //Get the children of this snarl and their clusters

        //The old distances from clusters to the bounds of the children, since we will be updating the distances
        //to represent distances to the parent
        vector<pair<size_t, size_t>> child_distances (tree_state.seed_count_prefix_sum.back(), 
                                                      make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));


        for (auto iter = child_range_start; iter != child_range_end; ++iter) {
            //Go through each child node of the netgraph

            NodeClusters& child_clusters = tree_state.all_node_clusters[iter->second];

            if (child_clusters.fragment_best_left > (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit) &&  
                child_clusters.fragment_best_right > (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit)) {
                //If everything is too far away to cluster, then skip it
                continue;
            }

            //This is true if this is the first time we see the outer loop's child. Used so we know if we need to calculate the distance to parents
            bool first_child = true;

            //Remember the distances for this child since they will get overwritten
            for (const pair<size_t, size_t>& head : child_clusters.read_cluster_heads) {
                child_distances[head.second + tree_state.seed_count_prefix_sum[head.first]] = 
                        make_pair(tree_state.all_seeds->at(head.first)->at(head.second).distance_left,
                                 tree_state.all_seeds->at(head.first)->at(head.second).distance_right);
            }

            for (auto inner_iter = child_range_start ; inner_iter != std::next(iter) ; ++inner_iter){
                //Go through other child net graph nodes up to and including i

                //Get the other node and its clusters
                NodeClusters& child_clusters_j = tree_state.all_node_clusters[inner_iter->second];

                if (child_clusters_j.fragment_best_left > (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit) &&  
                    child_clusters_j.fragment_best_right > (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit)) {
                    continue;
                }

#ifdef DEBUG_CLUSTER
                cerr << "\tComparing two children of " << distance_index.net_handle_as_string(snarl_handle) << ": " 
                     << distance_index.net_handle_as_string(child_clusters.containing_net_handle) << " and " 
                     << distance_index.net_handle_as_string(child_clusters_j.containing_net_handle) << endl;
                     


#endif

                compare_and_combine_cluster_on_child_structures(tree_state, child_clusters, 
                        child_clusters_j, snarl_clusters, child_distances, false, first_child);
                first_child = false;
            }
        }
    } else {
        //IF this is a simple snarl

        for (auto iter = child_range_start; iter != child_range_end; ++iter) {
            //Go through each child node of the netgraph and add its clusters to the snarl

            NodeClusters& node_clusters = tree_state.all_node_clusters[iter->second]; 

            //Add the cluster heads
            for (auto& cluster_head : node_clusters.read_cluster_heads) {
                snarl_clusters.read_cluster_heads.emplace(cluster_head);
            }

            //Update the distances
            //Because the orientation of the nodes was determined by the orientation of the chain,
            //the orientation relative to the snarl is correct
            for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
                if (read_num == 0) {
                    snarl_clusters.read_best_left.first = std::min(snarl_clusters.read_best_left.first,
                                                                   node_clusters.read_best_left.first);
                    snarl_clusters.read_best_right.first = std::min(snarl_clusters.read_best_right.first,
                                                                    node_clusters.read_best_right.first);
                } else {
                    snarl_clusters.read_best_left.second = std::min(snarl_clusters.read_best_left.second,
                                                                   node_clusters.read_best_left.second);
                    snarl_clusters.read_best_right.second = std::min(snarl_clusters.read_best_right.second,
                                                                    node_clusters.read_best_right.second);
                }
            }
            snarl_clusters.fragment_best_left = std::min(snarl_clusters.fragment_best_left,
                                                          node_clusters.fragment_best_left);
            snarl_clusters.fragment_best_right = std::min(snarl_clusters.fragment_best_right,
                                                           node_clusters.fragment_best_right);


        }
    }

#ifdef DEBUG_CLUSTER
    cerr << "\tFound clusters on " << distance_index.net_handle_as_string(snarl_handle) << endl;
    cerr << "\t   with best left and right values: " << snarl_clusters.fragment_best_left << " "
         << snarl_clusters.fragment_best_right << endl;
    bool got_left = false;
    bool got_right = false;
    for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
        cerr << "\t\tfor read num " << read_num << " best left: " << (read_num == 0 ? snarl_clusters.read_best_left.first : snarl_clusters.read_best_left.second) 
             << " best right: " << (read_num == 0 ? snarl_clusters.read_best_right.first : snarl_clusters.read_best_right.second) << endl;
        bool got_read_left=false;
        bool got_read_right = false;
        bool any_clusters = false;
        for (pair<size_t,size_t> c : snarl_clusters.read_cluster_heads) {
            if (c.first == read_num) {
                any_clusters = true;
                pair<size_t, size_t> dists (tree_state.all_seeds->at(c.first)->at(c.second).distance_left, 
                                            tree_state.all_seeds->at(c.first)->at(c.second).distance_right);
                cerr << "\t\t" << c.first << ":"<<c.second << ": left: " << dists.first << " right : " << dists.second << ": ";
                bool has_seeds = false;
                for (size_t x = 0 ; x < tree_state.all_seeds->at(c.first)->size() ; x++) {
                    if (tree_state.read_union_find[c.first].find_group(x) == c.second) {
                        cerr << tree_state.all_seeds->at(c.first)->at(x).pos << " ";
                        has_seeds = true;
                    }
                }
                assert(dists.first == std::numeric_limits<size_t>::max() || dists.first >= (read_num == 0 ? snarl_clusters.read_best_left.first : snarl_clusters.read_best_left.second));
                assert(dists.second == std::numeric_limits<size_t>::max() || dists.second >= (read_num == 0 ? snarl_clusters.read_best_right.first : snarl_clusters.read_best_right.second));
                assert(dists.first == std::numeric_limits<size_t>::max() || dists.first >= snarl_clusters.fragment_best_left);
                assert(dists.second == std::numeric_limits<size_t>::max() || dists.second >= snarl_clusters.fragment_best_right);
                if (dists.first == snarl_clusters.fragment_best_left) {got_left = true;}
                if (dists.second == snarl_clusters.fragment_best_right) {got_right = true;}
                if (dists.first == (read_num == 0 ? snarl_clusters.read_best_left.first : snarl_clusters.read_best_left.second)) {got_read_left = true;}
                if (dists.second == (read_num == 0 ? snarl_clusters.read_best_right.first : snarl_clusters.read_best_right.second)) {got_read_right = true;}
                cerr << endl;
                //assert(has_seeds);
            }
        }
        //assert(!any_clusters ||got_read_left ||  (read_num == 0 ? snarl_clusters.read_best_left.first : snarl_clusters.read_best_left.second) == std::numeric_limits<size_t>::max());
        //assert(!any_clusters ||got_read_right ||  (read_num == 0 ? snarl_clusters.read_best_right.first : snarl_clusters.read_best_right.second) == std::numeric_limits<size_t>::max());
    }
    //assert(got_left);
    //assert(got_right);

    //for (pair<pair<size_t, size_t>, pair<size_t, size_t>> group_id : snarl_clusters.read_cluster_heads) {
    //    assert (group_id.first.second == tree_state.read_union_find[group_id.first.first].find_group(group_id.first.second));
    //}
#endif
};



void SnarlDistanceIndexClusterer::cluster_one_chain(TreeState& tree_state, size_t chain_clusters_index, 
        const std::vector<ParentToChildMap::ParentChildValues>::iterator& chain_range_start,
        const std::vector<ParentToChildMap::ParentChildValues>::iterator& chain_range_end,
        bool only_seeds, bool is_top_level_chain) const {
#ifdef DEBUG_CLUSTERS
    assert(distance_index.is_chain(chain_clusters.containing_net_handle));
    //if (only_seeds) {
    //    for (auto child : children_in_chain) {
    //        assert(!std::get<3>(child));
    //    }
    //} else {
    //    bool is_only_seeds = true;
    //    for (auto child : children_in_chain) {
    //        if (std::get<3>(child)) {
    //            is_only_seeds=false;
    //        }
    //    }
    //    assert(!is_only_seeds);
    //}
#endif

    NodeClusters& chain_clusters = tree_state.all_node_clusters[chain_clusters_index];
    net_handle_t& chain_handle = chain_clusters.containing_net_handle;


    if (!chain_clusters.is_trivial_chain && ! is_top_level_chain) {
        //If we need it, get the values from the distance index:
        //is_looping_chain, node_length, the end boundary node, and the end component
        //THese only get used if we need the distances to the ends of the chain
        chain_clusters.set_chain_values(distance_index);
    }


    if (only_seeds && !chain_clusters.is_looping_chain && 
        (chain_clusters.chain_component_end == 0 
           || chain_clusters.chain_component_end == std::numeric_limits<size_t>::max())) {
        //If there are only seeds in the chain (and the chain doesn't loop and isn't a multicomponent chain), 
        //then cluster by walking through the seeds
        //This also does the work of clustering a trivial chain (which is just a node), which should be the same amount of work as using cluster_one_node

        std::function<std::tuple<size_t, size_t, size_t>(const ParentToChildMap::ParentChildValues&)> get_offset_from_indices = 
            [&](const ParentToChildMap::ParentChildValues& seed_index){
                //This function returns a tuple of <read num, seed num, left offset>
                return std::make_tuple(seed_index.child_index1, seed_index.child_index2, seed_index.child_offset); 
        };
        cluster_seeds_on_linear_structure(tree_state, chain_clusters, chain_range_start, chain_range_end, chain_clusters.node_length, get_offset_from_indices, is_top_level_chain);
#ifdef DEBUG_CLUSTER
    cerr << "\tFound clusters on " << distance_index.net_handle_as_string(chain_handle) << endl;
    cerr << "\t   with best left and right values: " << chain_clusters.fragment_best_left << " "
         << chain_clusters.fragment_best_right << endl;
    bool got_left = false;
    bool got_right = false;
    for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
        cerr << "\t for read num " << read_num << " best left: " << (read_num == 0 ? chain_clusters.read_best_left.first : chain_clusters.read_best_left.second) << " best right: " << (read_num == 0 ? chain_clusters.read_best_right.first : chain_clusters.read_best_right.second) << endl;
        bool got_read_left=false;
        bool got_read_right = false;
        bool any_clusters = false;
        for (pair<size_t,size_t> c : chain_clusters.read_cluster_heads) {
            if (c.first == read_num) {
                any_clusters = true;
                pair<size_t, size_t> dists (tree_state.all_seeds->at(c.first)->at(c.second).distance_left,
                                            tree_state.all_seeds->at(c.first)->at(c.second).distance_right);
                cerr << "\t\t" << c.first << ":"<<c.second << ": left: " << dists.first << " right : " << dists.second << ": ";
                bool has_seeds = false;
                for (size_t x = 0 ; x < tree_state.all_seeds->at(c.first)->size() ; x++) {
                    if (tree_state.read_union_find[c.first].find_group(x) == c.second) {
                        cerr << tree_state.all_seeds->at(c.first)->at(x).pos << " ";
                        has_seeds = true;
                    }
                }
                //assert(dists.first == std::numeric_limits<size_t>::max() || dists.first >= (read_num == 0 ? chain_clusters.read_best_left.first : chain_clusters.read_best_left.second));
                //assert(dists.second == std::numeric_limits<size_t>::max() || dists.second >= (read_num == 0 ? chain_clusters.read_best_right.first : chain_clusters.read_best_right.second));
                //assert(dists.first == std::numeric_limits<size_t>::max() || dists.first >= chain_clusters.fragment_best_left);
                //assert(dists.second == std::numeric_limits<size_t>::max() || dists.second >= chain_clusters.fragment_best_right);
                if (dists.first == chain_clusters.fragment_best_left) {got_left = true;}
                if (dists.second == chain_clusters.fragment_best_right) {got_right = true;}
                if (dists.first == (read_num == 0 ? chain_clusters.read_best_left.first : chain_clusters.read_best_left.second)) {got_read_left = true;}
                if (dists.second == (read_num == 0 ? chain_clusters.read_best_right.first : chain_clusters.read_best_right.second)) {got_read_right = true;}
                cerr << endl;
                //assert(has_seeds);
            }
        }
        //assert(!any_clusters ||got_read_left ||  (read_num == 0 ? chain_clusters.read_best_left.first : chain_clusters.read_best_left.second) == std::numeric_limits<size_t>::max());
        //assert(!any_clusters ||got_read_right ||  (read_num == 0 ? chain_clusters.read_best_right.first : chain_clusters.read_best_right.second) == std::numeric_limits<size_t>::max());
    }
    //assert(got_left);
    //assert(got_right);

    for (pair<size_t, size_t> group_id : chain_clusters.read_cluster_heads) {
        //assert (group_id.first.second == tree_state.read_union_find[group_id.first.first].find_group(group_id.first.second));
    }
#endif
        return;

    }


#ifdef DEBUG_CLUSTER
    cerr << "Cluster chain " << distance_index.net_handle_as_string(chain_handle) << endl;
    cerr << "\t chain has " << chain_range_end - chain_range_start << " children" << endl;
#endif

    /*Go through the chain child by child 
     *
     * As we walk through the chain, keep track of all clusters found up to the current child. 
     * So after we saw a child, the chain knows all clusters with distances up to the right side of the child
    * 
    * For each child, 
    *  - check if the clusters in the child can be combined with each other by walking out then back in through the chain
    *  - For snarls, update distances to the ends of the child (taking into account the distances to loop around in the chain)
    *  - compare and combine with clusters of the chain that get build as we walk along it
    * 
    *  - after combining clusters of the current child, remove redundant cluster heads from the chain clusters
    *
     * As we walk along the chain, we maintain clusters of the chain up to the last node we saw (the later
     * boundary node of a snarl relative to the chain)
     *
     * The clusters in the chain have left distances to the beginning of the chain and right distances
     * to the last thing we saw
     */


    //These are the things that we update as we walk through the chain: Information about the last child we saw
    //Initialized to the first child in the chain  

    //The last child we saw
    ParentToChildMap::ParentChildValues last_child = *chain_range_start;

    //And values we need to save from the last child
    //If the last child is a snarl, get it from the NodeClusters otherwise from the seed's cache
    size_t last_prefix_sum = last_child.child_offset;
    size_t last_length = last_child.child_index2 == std::numeric_limits<size_t>::max() 
                                          ? tree_state.all_node_clusters[last_child.child_index1].node_length
                                          : tree_state.all_seeds->at(last_child.child_index1)->at(last_child.child_index2).minimizer_cache.node_length;
    size_t last_chain_component_end = last_child.child_index2 == std::numeric_limits<size_t>::max() 
                                          ? tree_state.all_node_clusters[last_child.child_index1].chain_component_end
                                          : last_child.child_chain_component;

    //These are clusters that we don't want to consider as we walk through the chain but that 
    //we want to remember after we're done with the chain because the left distance is small
    vector<ClusterIndices> cluster_heads_to_add_again;

    //For remembering the best left distances of the chain, we only need to check for the smallest chain distance left
    //for the children up to the first node
    bool found_first_node = false;
    pair<bool, bool> found_first_node_by_read = std::make_pair(false, false);


    for (auto current_iter = chain_range_start ; current_iter != chain_range_end ; ++current_iter) {
        /*
         * Snarls and nodes are in the order that they are traversed in the chain
         * For each child, compare all of the clusters on the child to clusters of the chain so far
         * The clusters of the chain have right distances up to the end of the last child seen
         */

        ParentToChildMap::ParentChildValues& child_clusters = *current_iter;

        if (child_clusters.child_index2 == std::numeric_limits<size_t>::max()){

            //If this is a snarl, then cluster the children here
            add_snarl_to_chain_clusters(tree_state, chain_clusters, last_child, last_prefix_sum, last_length,
                                       last_chain_component_end, cluster_heads_to_add_again, found_first_node, found_first_node_by_read,
                                       child_clusters, current_iter == chain_range_start, std::next(current_iter) == chain_range_end, is_top_level_chain);
        } else {

            add_seed_to_chain_clusters(tree_state, chain_clusters, last_child,  last_prefix_sum, last_length,
                                       last_chain_component_end, cluster_heads_to_add_again, found_first_node, found_first_node_by_read,
                                       child_clusters, current_iter == chain_range_start, std::next(current_iter) == chain_range_end, is_top_level_chain);
        }

#ifdef DEBUG_CLUSTER
    cerr << "\tintermediate clusters on " << distance_index.net_handle_as_string(chain_handle) << " after child " << distance_index.net_handle_as_string(child_clusters.child_handle) << endl;
    cerr << "\t   with best left and right values: " << chain_clusters.fragment_best_left << " "
         << chain_clusters.fragment_best_right << endl;
    bool got_left = false;
    bool got_right = false;
    for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
        cerr << "\t for read num " << read_num << " best left: " << (read_num == 0 ? chain_clusters.read_best_left.first : chain_clusters.read_best_left.second) << " best right: " << (read_num == 0 ? chain_clusters.read_best_right.first : chain_clusters.read_best_right.second) << endl;
        bool got_read_left=false;
        bool got_read_right = false;
        bool any_clusters = false;
        for (pair<size_t,size_t> c : chain_clusters.read_cluster_heads) {
            if (c.first == read_num) {
                any_clusters = true;
                pair<size_t, size_t> dists (tree_state.all_seeds->at(c.first)->at(c.second).distance_left,
                                           tree_state.all_seeds->at(c.first)->at(c.second).distance_right);
                cerr << "\t\t" << c.first << ":"<<c.second << ": left: " << dists.first << " right : " << dists.second << ": ";
                bool has_seeds = false;
                for (size_t x = 0 ; x < tree_state.all_seeds->at(c.first)->size() ; x++) {
                    if (tree_state.read_union_find[c.first].find_group(x) == c.second) {
                        cerr << tree_state.all_seeds->at(c.first)->at(x).pos << " ";
                        has_seeds = true;
                    }
                }
            }
        }
    }
    vector<Seed> ordered_seeds;
    for (size_t i = 0 ; i < tree_state.all_seeds->size() ; i++) {
        const auto v = tree_state.all_seeds->at(i);
        for ( auto x : *v) {
            ordered_seeds.push_back(x);
        }
    }
    cerr << "Found intermediate fragment clusters : " << endl;
    for (auto group : tree_state.fragment_union_find.all_groups()){
        cerr << "\t";
        for (size_t c : group) {
           cerr << ordered_seeds[c].pos << " ";
        }
        cerr << endl;
    }
#endif
    }
    //Add back clusters we skipped
    for (auto& cluster_head : cluster_heads_to_add_again) {
        chain_clusters.read_cluster_heads.emplace(cluster_head.read_num, cluster_head.cluster_num);
        tree_state.all_seeds->at(cluster_head.read_num)->at(cluster_head.cluster_num).distance_left = cluster_head.distance_left;
        tree_state.all_seeds->at(cluster_head.read_num)->at(cluster_head.cluster_num).distance_right = cluster_head.distance_right;
        chain_clusters.fragment_best_left = std::min(chain_clusters.fragment_best_left, cluster_head.distance_left);
        if (cluster_head.read_num == 0) {
            chain_clusters.read_best_left.first = std::min(chain_clusters.read_best_left.first, cluster_head.distance_left);
        } else {
            chain_clusters.read_best_left.second = std::min(chain_clusters.read_best_left.second, cluster_head.distance_left);
        }

    }


    //If the chain loops, then we also have to compare the first thing we saw to the last things

    if (chain_clusters.is_looping_chain){
#ifdef DEBUG_CLUSTER
        cerr << "Check connectivity around a looping chain" << endl;
    cerr << "\tFound clusters on " << distance_index.net_handle_as_string(chain_handle) << endl;
    cerr << "\t   with best left and right values: " << chain_clusters.fragment_best_left << " "
         << chain_clusters.fragment_best_right << endl;
    bool got_left = false;
    bool got_right = false;
    for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
        cerr << "\t for read num " << read_num << " best left: " << (read_num == 0 ? chain_clusters.read_best_left.first : chain_clusters.read_best_left.second) << " best right: " << (read_num == 0 ? chain_clusters.read_best_right.first : chain_clusters.read_best_right.second) << endl;
        bool got_read_left=false;
        bool got_read_right = false;
        bool any_clusters = false;
        for (pair<size_t,size_t> c : chain_clusters.read_cluster_heads) {
            if (c.first == read_num) {
                any_clusters = true;
                pair<size_t, size_t> dists (tree_state.all_seeds->at(c.first)->at(c.second).distance_left,
                                            tree_state.all_seeds->at(c.first)->at(c.second).distance_right);
                cerr << "\t\t" << c.first << ":"<<c.second << ": left: " << dists.first << " right : " << dists.second << ": ";
                bool has_seeds = false;
                for (size_t x = 0 ; x < tree_state.all_seeds->at(c.first)->size() ; x++) {
                    if (tree_state.read_union_find[c.first].find_group(x) == c.second) {
                        cerr << tree_state.all_seeds->at(c.first)->at(x).pos << " ";
                        has_seeds = true;
                    }
                }
                //assert(dists.first == std::numeric_limits<size_t>::max() || dists.first >= (read_num == 0 ? chain_clusters.read_best_left.first : chain_clusters.read_best_left.second));
                //assert(dists.second == std::numeric_limits<size_t>::max() || dists.second >= (read_num == 0 ? chain_clusters.read_best_right.first : chain_clusters.read_best_right.second));
                //assert(dists.first == std::numeric_limits<size_t>::max() || dists.first >= chain_clusters.fragment_best_left);
                //assert(dists.second == std::numeric_limits<size_t>::max() || dists.second >= chain_clusters.fragment_best_right);
                if (dists.first == chain_clusters.fragment_best_left) {got_left = true;}
                if (dists.second == chain_clusters.fragment_best_right) {got_right = true;}
                if (dists.first == (read_num == 0 ? chain_clusters.read_best_left.first : chain_clusters.read_best_left.second)) {got_read_left = true;}
                if (dists.second == (read_num == 0 ? chain_clusters.read_best_right.first : chain_clusters.read_best_right.second)) {got_read_right = true;}
                cerr << endl;
                //assert(has_seeds);
            }
        }
        //assert(!any_clusters ||got_read_left ||  (read_num == 0 ? chain_clusters.read_best_left.first : chain_clusters.read_best_left.second) == std::numeric_limits<size_t>::max());
        //assert(!any_clusters ||got_read_right ||  (read_num == 0 ? chain_clusters.read_best_right.first : chain_clusters.read_best_right.second) == std::numeric_limits<size_t>::max());
    }
    //assert(got_left);
    //assert(got_right);

    for (pair<size_t, size_t> group_id : chain_clusters.read_cluster_heads) {
        //assert (group_id.first.second == tree_state.read_union_find[group_id.first.first].find_group(group_id.first.second));
    }
#endif
        pair<size_t, size_t> combined_cluster_by_read = std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()); 
        size_t combined_cluster_fragment = std::numeric_limits<size_t>::max();
        for (auto& cluster_head : chain_clusters.read_cluster_heads) {
            size_t read_num = cluster_head.first;
            size_t cluster_num=cluster_head.second;
            size_t dist_left = tree_state.all_seeds->at(read_num)->at(cluster_num).distance_left;
            size_t dist_right = tree_state.all_seeds->at(read_num)->at(cluster_num).distance_right;

            size_t distance_between_left_right =SnarlDistanceIndex::minus(SnarlDistanceIndex::sum(dist_left, (read_num == 0 ? chain_clusters.read_best_right.first : chain_clusters.read_best_right.second)), 1); 
            size_t distance_between_right_left = SnarlDistanceIndex::minus( SnarlDistanceIndex::sum(dist_right, (read_num == 0 ? chain_clusters.read_best_left.first : chain_clusters.read_best_left.second)), 1);
            if (distance_between_left_right <= tree_state.read_distance_limit
             || distance_between_right_left <= tree_state.read_distance_limit) {
                //If we can combine the read
                size_t& combined_cluster_num = (read_num == 0 ? combined_cluster_by_read.first : combined_cluster_by_read.second);
                if (combined_cluster_num == std::numeric_limits<size_t>::max()) {
                    combined_cluster_num = cluster_num;
                } else {
                    combined_cluster_num = tree_state.read_union_find[read_num].union_groups(
                                                combined_cluster_num, cluster_num);
                }
            }
            size_t distance_between_left_right_fragment = SnarlDistanceIndex::minus(SnarlDistanceIndex::sum(dist_left, chain_clusters.fragment_best_right) ,1);
            size_t distance_between_right_left_fragment = SnarlDistanceIndex::minus(SnarlDistanceIndex::sum(dist_right, chain_clusters.fragment_best_left) ,1);
            if (tree_state.fragment_distance_limit != 0 &&
                (distance_between_left_right_fragment <= tree_state.fragment_distance_limit
              || distance_between_right_left_fragment <= tree_state.fragment_distance_limit)) {

                if (combined_cluster_fragment != std::numeric_limits<size_t>::max()) {
                    combined_cluster_fragment = tree_state.fragment_union_find.union_groups(combined_cluster_fragment, 
                            cluster_num + tree_state.seed_count_prefix_sum[read_num]);
                } else {
                    combined_cluster_fragment = cluster_num + tree_state.seed_count_prefix_sum[read_num];
                }
            }
        }

    }

#ifdef DEBUG_CLUSTER
    cerr << "\tFound clusters on " << distance_index.net_handle_as_string(chain_handle) << endl;
    cerr << "\t   with best left and right values: " << chain_clusters.fragment_best_left << " "
         << chain_clusters.fragment_best_right << endl;
    bool got_left = false;
    bool got_right = false;
    for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
        cerr << "\t for read num " << read_num << " best left: " << (read_num == 0 ? chain_clusters.read_best_left.first : chain_clusters.read_best_left.second) << " best right: " << (read_num == 0 ? chain_clusters.read_best_right.first : chain_clusters.read_best_right.second) << endl;
        bool got_read_left=false;
        bool got_read_right = false;
        bool any_clusters = false;
        for (pair<size_t,size_t> c : chain_clusters.read_cluster_heads) {
            if (c.first == read_num) {
                any_clusters = true;
                pair<size_t, size_t> dists (tree_state.all_seeds->at(c.first)->at(c.second).distance_left,
                                            tree_state.all_seeds->at(c.first)->at(c.second).distance_right);
                cerr << "\t\t" << c.first << ":"<<c.second << ": left: " << dists.first << " right : " << dists.second << ": ";
                bool has_seeds = false;
                for (size_t x = 0 ; x < tree_state.all_seeds->at(c.first)->size() ; x++) {
                    if (tree_state.read_union_find[c.first].find_group(x) == c.second) {
                        cerr << tree_state.all_seeds->at(c.first)->at(x).pos << " ";
                        has_seeds = true;
                    }
                }
                //assert(dists.first == std::numeric_limits<size_t>::max() || dists.first >= (read_num == 0 ? chain_clusters.read_best_left.first : chain_clusters.read_best_left.second));
                //assert(dists.second == std::numeric_limits<size_t>::max() || dists.second >= (read_num == 0 ? chain_clusters.read_best_right.first : chain_clusters.read_best_right.second));
                //assert(dists.first == std::numeric_limits<size_t>::max() || dists.first >= chain_clusters.fragment_best_left);
                //assert(dists.second == std::numeric_limits<size_t>::max() || dists.second >= chain_clusters.fragment_best_right);
                if (dists.first == chain_clusters.fragment_best_left) {got_left = true;}
                if (dists.second == chain_clusters.fragment_best_right) {got_right = true;}
                if (dists.first == (read_num == 0 ? chain_clusters.read_best_left.first : chain_clusters.read_best_left.second)) {got_read_left = true;}
                if (dists.second == (read_num == 0 ? chain_clusters.read_best_right.first : chain_clusters.read_best_right.second)) {got_read_right = true;}
                cerr << endl;
                //assert(has_seeds);
            }
        }
        //assert(!any_clusters ||got_read_left ||  (read_num == 0 ? chain_clusters.read_best_left.first : chain_clusters.read_best_left.second) == std::numeric_limits<size_t>::max());
        //assert(!any_clusters ||got_read_right ||  (read_num == 0 ? chain_clusters.read_best_right.first : chain_clusters.read_best_right.second) == std::numeric_limits<size_t>::max());
    }
    //assert(got_left);
    //assert(got_right);

    for (pair<size_t, size_t> group_id : chain_clusters.read_cluster_heads) {
        //assert (group_id.first.second == tree_state.read_union_find[group_id.first.first].find_group(group_id.first.second));
    }
#endif
}
void SnarlDistanceIndexClusterer::add_seed_to_chain_clusters(TreeState& tree_state, NodeClusters& chain_clusters,
                                ParentToChildMap::ParentChildValues& last_child,
                                size_t& last_prefix_sum, size_t& last_length, size_t& last_chain_component_end,
                                vector<ClusterIndices>& cluster_heads_to_add_again,
                                bool& found_first_node, pair<bool, bool>& found_first_node_by_read,
                                const ParentToChildMap::ParentChildValues& current_child, bool is_first_child, 
                                bool is_last_child, bool skip_distances_to_ends) const {
    size_t read_num = current_child.child_index1;
    size_t cluster_num = current_child.child_index2;
    net_handle_t& chain_handle = chain_clusters.containing_net_handle;
    Seed& current_child_seed = tree_state.all_seeds->at(read_num)->at(cluster_num);
    /*
    Get a bunch of distances from the current child that will be used to calculate distance
    from the last child
    */
    
    //The distance from the right side of the last child to the left side of this child 
    //(relative to the orientation of the chain
    size_t distance_from_last_child_to_current_child = std::numeric_limits<size_t>::max();
    if (!is_first_child) {
        //If this isn't the first child we're looking at
        if (last_child.child_handle == current_child.child_handle) {
            //This can happen if the last thing was also a seed on the same node
            distance_from_last_child_to_current_child = 0; 
        } else if ( last_chain_component_end == current_child_seed.minimizer_cache.chain_component) {
            //If this child is in the same component as the last one
            if (last_length == std::numeric_limits<size_t>::max()) {
                //If the last length is infinite, then is must be a snarl that is not start-end reachable, so the distance
                //from the last child is the same as the distance from the start of the chain (the start of this compnent)
                distance_from_last_child_to_current_child = current_child_seed.minimizer_cache.prefix_sum;
            } else {
                size_t distance_from_chain_start_to_last_node = SnarlDistanceIndex::sum(last_prefix_sum,last_length);
    
                //Distance is the current node's prefix sum minus the distance from the start of the chain to the last node
                distance_from_last_child_to_current_child = SnarlDistanceIndex::minus(current_child_seed.minimizer_cache.prefix_sum, 
                                                distance_from_chain_start_to_last_node); 
            }
        }
    }
    
    
    
    //The distance to add to get to the end of the chain. Only matters if this is the last thing in the chain
    //The distances will include the distance to the end of a trivial chain,
    //so we can't rely on distance_in_parent to know when the distance should be 0
    
    size_t distance_from_current_end_to_end_of_chain;
    if (!is_last_child || skip_distances_to_ends) {
        //If this isn't the last child in the chain, then we only want the distance to the end of the current child
    
        distance_from_current_end_to_end_of_chain = 0;
    } else if (SnarlDistanceIndex::get_record_offset(current_child.child_handle) == SnarlDistanceIndex::get_record_offset(chain_clusters.end_in)) {
        //If this is the last node in the chain
        if (chain_clusters.chain_component_end != current_child_seed.minimizer_cache.chain_component) { 
            //If they aren't in the same component
            distance_from_current_end_to_end_of_chain = std::numeric_limits<size_t>::max();
        } else {
            distance_from_current_end_to_end_of_chain = 0;
        }
    } else if (chain_clusters.chain_component_end != current_child_seed.minimizer_cache.chain_component) { 
        //If they aren't in the same component
        distance_from_current_end_to_end_of_chain = std::numeric_limits<size_t>::max();
    } else {
    
        //Length of the chain - (prefix sum + node length of the current node)
        distance_from_current_end_to_end_of_chain = SnarlDistanceIndex::minus(chain_clusters.node_length, 
                    SnarlDistanceIndex::sum(current_child_seed.minimizer_cache.prefix_sum, 
                                            current_child_seed.minimizer_cache.node_length));
    
    }

#ifdef DEBUG_CLUSTER
    cerr << "\tDistance from last child to this one: " << distance_from_last_child_to_current_child << endl;
    cerr << "\tDistance from start of chain to the left side of this one: " << (current_child_seed.minimizer_cache.chain_component != 0 ? std::numeric_limits<size_t>::max() : current_child_seed.minimizer_cache.prefix_sum) << endl;
    cerr << "\tDistance to get to the end of the chain: " << distance_from_current_end_to_end_of_chain << endl;
#endif


    if (last_child.child_handle != current_child.child_handle &&
        SnarlDistanceIndex::sum(distance_from_last_child_to_current_child, chain_clusters.fragment_best_right) 
          > (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit)) {
#ifdef DEBUG_CLUSTER
            cerr << "This child is too far away from the last one to cluster anything" << endl;
#endif
        //If the distance from the last cluster is too far to cluster anything
        if (!skip_distances_to_ends) {
            for (auto& cluster_head : chain_clusters.read_cluster_heads) {
                //For each of the chain clusters, remember the ones that are still reachable from the left side of the chain
                pair<size_t, size_t> dists (tree_state.all_seeds->at(cluster_head.first)->at(cluster_head.second).distance_left,
                                            tree_state.all_seeds->at(cluster_head.first)->at(cluster_head.second).distance_right);
                if (dists.first <= (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit)) {
                    //If this cluster can be clustered outside of the chain, remember to add it back
                    cluster_heads_to_add_again.emplace_back();
                    cluster_heads_to_add_again.back().read_num = cluster_head.first;
                    cluster_heads_to_add_again.back().cluster_num = cluster_head.second; 
                    cluster_heads_to_add_again.back().distance_left = dists.first; 
                    cluster_heads_to_add_again.back().distance_right = std::numeric_limits<size_t>::max();
                }
            }
        }
    
        //Now clear the chain's list of clusters
        chain_clusters.read_cluster_heads.clear();
        
    
        //Update the distances stored in the seed to reach the ends of the chain
        //The distance left and right of the seed are currently oriented relative to the chain
    
        //The current left distance is infinite if it is not in the first component of a multicomponent chain
        if (current_child_seed.minimizer_cache.chain_component != 0) {
            //If this node isn't in the first component of the chain
            current_child_seed.distance_left = std::numeric_limits<size_t>::max();
        } else {
            //Prefix sum + offset of the seed in the node
            current_child_seed.distance_left = SnarlDistanceIndex::sum(current_child_seed.distance_left, 
                                                                       current_child_seed.minimizer_cache.prefix_sum);
        }
        current_child_seed.distance_right = SnarlDistanceIndex::sum(current_child_seed.distance_right, 
                                                       distance_from_current_end_to_end_of_chain);
    
        //Add the cluster to the chain
        chain_clusters.read_cluster_heads.emplace(read_num, cluster_num);
    
        //Update the best distances on the chain
        if (!found_first_node) {
            chain_clusters.fragment_best_left = std::min(chain_clusters.fragment_best_left, current_child_seed.distance_left);
        }
        if (! (read_num == 0 ? found_first_node_by_read.first : found_first_node_by_read.second)) {
            (read_num == 0 ? chain_clusters.read_best_left.first : chain_clusters.read_best_left.second) = std::min((read_num == 0 ? chain_clusters.read_best_left.first : chain_clusters.read_best_left.second), current_child_seed.distance_left);
        }
        //Since this child is a seed on a node, it's right distance will be the best one for the chain so far
        chain_clusters.fragment_best_right = current_child_seed.distance_right;
        (read_num == 0 ? chain_clusters.read_best_right.first : chain_clusters.read_best_right.second) = current_child_seed.distance_right;
    
        //Also update the best right distances to the end of this node for the other read
        //Since it was too far away from this node, it will be too far away from anything else and it can just be infinite
        for (size_t chain_read_num = 0 ; chain_read_num < tree_state.all_seeds->size() ; chain_read_num++) {
            if (chain_read_num != read_num) {
                if (chain_read_num == 0) {
                    chain_clusters.read_best_right.first = std::numeric_limits<size_t>::max();
                } else {
                    chain_clusters.read_best_right.second = std::numeric_limits<size_t>::max();
                }
            }
        }
    
    } else {
        //Otherwise, check to see if anything on the current child can be combined with 
        //anything in the chain thus far

        //The distance from the right side of the last child to the right side of this child, which is
        //the distance we need to update the chain clusters to the end of this child
        //This isn't quite right for the first thing in the chain but it doesn't matter because it only
        //gets added to chain clusters
        //IF it gets calculated, then it's the distance from the last child to this node + the length
        //of this node (the first value in the cache)
        size_t distance_from_last_child_to_current_end = 
                distance_from_last_child_to_current_child == std::numeric_limits<size_t>::max() 
                        ? std::numeric_limits<size_t>::max() : 
                (last_child.child_handle == current_child.child_handle ? 0 
                    : SnarlDistanceIndex::sum(distance_from_last_child_to_current_child, current_child_seed.minimizer_cache.node_length));
    
        //The new distances from this child to the start of the chain and the end of this child (or the end of the chain if it's the last child)
        //Left distance is the prefix sum (or inf if the node isn't in the first component of the chain) + offset of seed in node
        //Right distance is the right offst of the seed in the node + the distance from the end of the node to the end of the chain 
        // (or 0 if it isn't the last thing in the chain)
        pair<size_t, size_t> new_distances = make_pair(
                SnarlDistanceIndex::sum(current_child_seed.distance_left, 
                        current_child_seed.minimizer_cache.chain_component != 0 ? std::numeric_limits<size_t>::max() 
                                                                             : current_child_seed.minimizer_cache.prefix_sum),
                SnarlDistanceIndex::sum(current_child_seed.distance_right, distance_from_current_end_to_end_of_chain)); 
    
    
        //Cluster heads to remove because they got combined with the current seed
        vector<pair<size_t, size_t>> to_remove;
        //And the new cluster containing the current seed, and possibly anything that gets combined with it
        ClusterIndices new_cluster (read_num, cluster_num, new_distances.first, new_distances.second);
    
        /**Go through the clusters on the chain up to this point and see if anything can
        be combined with the clusters on the child
        Also update the distances of the chain clusters to reach the end of this node
        */
        for (auto& chain_cluster_head : chain_clusters.read_cluster_heads) {
            //Each has distances up to the previous node
    
            const size_t chain_cluster_read_num = chain_cluster_head.first;
            const size_t chain_cluster_cluster_num = chain_cluster_head.second;
    
            //The distances of the chain cluster
            pair<size_t, size_t> chain_cluster_distances (tree_state.all_seeds->at(chain_cluster_read_num)->at(chain_cluster_cluster_num).distance_left,
                                        tree_state.all_seeds->at(chain_cluster_read_num)->at(chain_cluster_cluster_num).distance_right);
    
    
            //The distance between the current seed and the current chain cluster
            size_t distance_between = SnarlDistanceIndex::minus(
                    SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(chain_cluster_distances.second, 
                                             distance_from_last_child_to_current_child), 
                                             current_child_seed.distance_left),
                    1);
            if (!is_first_child && last_child.child_handle == current_child.child_handle) {
                //If the last child was the same as this child (seeds on the same node), 
                //then the distances right are including the current node, so subtract
                //the length of this node
                distance_between -= current_child_seed.minimizer_cache.node_length;
            }

#ifdef DEBUG_CLUSTER
            cerr << "\t\tCombine this seed " << read_num << ":" << cluster_num << endl;
#endif

                             
            if (chain_cluster_read_num == read_num && distance_between <= tree_state.read_distance_limit) {
#ifdef DEBUG_CLUSTER
                        cerr << "\t\tCombine chain cluster " << read_num << ":" << chain_cluster_cluster_num << endl;
#endif
                //Union the two clusters and remember the new cluster head
                new_cluster.cluster_num = tree_state.read_union_find.at(read_num).union_groups(chain_cluster_cluster_num, cluster_num);

                //Find the best distances of the two. The best right distance will always be the current seed's distance
                //And remember the new combined cluster head
                new_cluster.distance_left = std::min(new_cluster.distance_left, chain_cluster_distances.first);

                to_remove.emplace_back(chain_cluster_read_num, chain_cluster_cluster_num);

                //Try to union the fragment
                if (tree_state.fragment_distance_limit != 0 && distance_between <= tree_state.fragment_distance_limit) {
                    tree_state.fragment_union_find.union_groups(cluster_num + tree_state.seed_count_prefix_sum[read_num], 
                                chain_cluster_cluster_num + tree_state.seed_count_prefix_sum[chain_cluster_read_num]);
                }
    
            } else if (tree_state.fragment_distance_limit != 0 && distance_between <= tree_state.fragment_distance_limit) {
                //If we can union the fragments, then union them and keep the cluster around, updating the right distance

                tree_state.fragment_union_find.union_groups(cluster_num + tree_state.seed_count_prefix_sum[read_num], 
                            chain_cluster_cluster_num + tree_state.seed_count_prefix_sum[chain_cluster_read_num]);

                tree_state.all_seeds->at(chain_cluster_read_num)->at(chain_cluster_cluster_num).distance_right = 
                        SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(chain_cluster_distances.second, 
                                                 distance_from_last_child_to_current_end),
                                                 distance_from_current_end_to_end_of_chain);
            } else {
                //If this chain cluster doesn't get combined, then it is too far away to combine with anything later in the chain, 
                //so we remove it but remember to add it again if the left distance is small enough
                
                if (chain_cluster_distances.first <=
                        (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit)
                    && !skip_distances_to_ends) { 
                    //If the current chain cluster can still be reached from the left
                    tree_state.all_seeds->at(chain_cluster_read_num)->at(chain_cluster_cluster_num).distance_right = std::numeric_limits<size_t>::max();
                    cluster_heads_to_add_again.emplace_back(chain_cluster_read_num, chain_cluster_cluster_num,
                                                       chain_cluster_distances.first, std::numeric_limits<size_t>::max());
                }
                to_remove.emplace_back(chain_cluster_read_num, chain_cluster_cluster_num);
            }
              
        }
    
        //Remove all chain clusters that got combined with the current seed
        for (pair<size_t, size_t>& cluster_head : to_remove) {
            chain_clusters.read_cluster_heads.erase(cluster_head);
        }
    
        //Add the cluster of the current seed which may or may not have been combined
        chain_clusters.read_cluster_heads.emplace(new_cluster.read_num, new_cluster.cluster_num);
        tree_state.all_seeds->at(new_cluster.read_num)->at(new_cluster.cluster_num).distance_left = new_cluster.distance_left;
        tree_state.all_seeds->at(new_cluster.read_num)->at(new_cluster.cluster_num).distance_right = new_cluster.distance_right;
    
    
        //Update the best distances
        //Only update the left distances if we haven't seen a node in the chain yet
        if (!found_first_node) {
            chain_clusters.fragment_best_left = std::min(chain_clusters.fragment_best_left, new_distances.first);
        }

        //If we haven't found the first node for this read
        if (!(read_num == 0 ? found_first_node_by_read.first : found_first_node_by_read.second)){
            //Update the best left distance
            (read_num == 0 ? chain_clusters.read_best_left.first : chain_clusters.read_best_left.second) = std::min((read_num == 0 ? chain_clusters.read_best_left.first : chain_clusters.read_best_left.second), new_distances.first);
        }

        //Since this is a node, the best right distance will be this distance
        chain_clusters.fragment_best_right = new_distances.second; 
        (read_num == 0 ? chain_clusters.read_best_right.first : chain_clusters.read_best_right.second) = new_distances.second;

        //Also update the best right distances to the end of this node for clusters of the other read
        for (size_t chain_read_num = 0 ; chain_read_num < tree_state.all_seeds->size() ; chain_read_num++) {
            if (chain_read_num != read_num) {
                if (chain_read_num == 0) {
                    chain_clusters.read_best_right.first = SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(
                                                                    chain_clusters.read_best_right.first,
                                                                    distance_from_last_child_to_current_end),
                                                                    distance_from_current_end_to_end_of_chain);
                } else {
                    chain_clusters.read_best_right.second = SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(
                                                                    chain_clusters.read_best_right.second,
                                                                    distance_from_last_child_to_current_end),
                                                                    distance_from_current_end_to_end_of_chain);
                }
            }
        }
    }
    
    found_first_node = true;
    if (read_num == 0) {
        found_first_node_by_read.first = true;
    } else {
        found_first_node_by_read.second = true;
    }
    
    
    //Update the last node we saw to this one
    last_child = current_child;
    last_prefix_sum = current_child_seed.minimizer_cache.prefix_sum;
    last_length = current_child_seed.minimizer_cache.node_length;
    last_chain_component_end = current_child_seed.minimizer_cache.chain_component;

}

void SnarlDistanceIndexClusterer::add_snarl_to_chain_clusters(TreeState& tree_state, NodeClusters& chain_clusters,
                                ParentToChildMap::ParentChildValues& last_child, 
                                size_t& last_prefix_sum, size_t& last_length, size_t& last_chain_component_end,
                                vector<ClusterIndices>& cluster_heads_to_add_again,
                                bool& found_first_node, pair<bool, bool>& found_first_node_by_read,
                                const ParentToChildMap::ParentChildValues& current_child, bool is_first_child, 
                                bool is_last_child, bool skip_distances_to_ends) const {

    /*Define a helper function to update the distances in a child using the loop distances
     * in the chain
     */
     auto update_distances_on_same_child = [&] (NodeClusters& child_clusters) {
         //Distance to go forward (relative to the child) in the chain and back
         if (child_clusters.loop_left == std::numeric_limits<size_t>::max() && child_clusters.loop_right == std::numeric_limits<size_t>::max()) {
             return;
         }


         //Combined clusters in case we can combine anything
         vector<pair<pair<size_t, size_t>, pair<size_t, size_t>>> combined_left (tree_state.all_seeds->size(),
                make_pair(make_pair(std::numeric_limits<size_t>::max(), 0), make_pair(0,0)));
         vector<pair<pair<size_t, size_t>, pair<size_t, size_t>>> combined_right (tree_state.all_seeds->size(),
               make_pair(make_pair(std::numeric_limits<size_t>::max(), 0), make_pair(0,0))); 
         size_t combined_fragment_left = std::numeric_limits<size_t>::max();
         size_t combined_fragment_right = std::numeric_limits<size_t>::max();
         vector<pair<size_t, size_t>> to_erase;

         for (auto& child_cluster_head : child_clusters.read_cluster_heads) {
            //Go through each of the clusters on this child
            size_t read_num = child_cluster_head.first;
            size_t cluster_num = child_cluster_head.second;
            size_t old_left =  tree_state.all_seeds->at(read_num)->at(cluster_num).distance_left;
            size_t old_right = tree_state.all_seeds->at(read_num)->at(cluster_num).distance_right;
            //Get the new best distances for the cluster considering chain loops
            size_t updated_left = std::min(old_left, SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(old_right, child_clusters.loop_right), child_clusters.node_length));
            size_t updated_right = std::min(old_right, SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(old_left, child_clusters.loop_left), child_clusters.node_length));



            if (updated_left < old_left || updated_right < old_right ) {
                //Update the distances
                tree_state.all_seeds->at(read_num)->at(cluster_num).distance_left = updated_left;
                tree_state.all_seeds->at(read_num)->at(cluster_num).distance_right = updated_right;
                
                child_clusters.fragment_best_left = std::min(child_clusters.fragment_best_left,
                                                             updated_left);
                child_clusters.fragment_best_right = std::min(child_clusters.fragment_best_right,
                                                             updated_right);
                (read_num == 0 ? child_clusters.read_best_left.first : child_clusters.read_best_left.second) = std::min((read_num == 0 ? child_clusters.read_best_left.first : child_clusters.read_best_left.second),updated_left);
                (read_num == 0 ? child_clusters.read_best_right.first : child_clusters.read_best_right.second) = std::min((read_num == 0 ? child_clusters.read_best_right.first : child_clusters.read_best_right.second),
                                                                       updated_right);
            }

            //Now see if we can combine this cluster with anything else 
            //The distance between this cluster and anything else taking the left loop
            size_t distance_between_left = SnarlDistanceIndex::minus(
                    SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(updated_left, 
                                             child_clusters.loop_left),
                                             (read_num == 0 ? child_clusters.read_best_left.first : child_clusters.read_best_left.second)),
                    1); 
            size_t distance_between_right = SnarlDistanceIndex::minus(
                     SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(updated_right, 
                                              child_clusters.loop_right),
                                              (read_num == 0 ? child_clusters.read_best_right.first : child_clusters.read_best_right.second)),
                     1); 
            size_t distance_between_left_fragment = SnarlDistanceIndex::minus(
                      SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(updated_left, 
                                                child_clusters.loop_left),
                                                child_clusters.fragment_best_left),
                      1); 
            size_t distance_between_right_fragment = SnarlDistanceIndex::minus(
                       SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(updated_right, 
                                                child_clusters.loop_right),
                                                child_clusters.fragment_best_right),
                       1); 
            pair<size_t, size_t> cluster_head = make_pair(read_num, cluster_num);
            if (distance_between_left <= tree_state.read_distance_limit) {
                //Combine it left
                to_erase.emplace_back(cluster_head);
                if (combined_left[read_num].first.first == std::numeric_limits<size_t>::max()){
                    combined_left[read_num] = make_pair(cluster_head, 
                         make_pair(updated_left, updated_right));
                } else {
                    to_erase.emplace_back(combined_left[read_num].first);
                    combined_left[read_num] = make_pair(
                         make_pair(read_num, 
                                   tree_state.read_union_find.at(read_num).union_groups(cluster_num, 
                                                                                        combined_left[read_num].first.second)), 
                         make_pair(std::min(updated_left, combined_left[read_num].second.first), 
                                   std::min(updated_right, combined_left[read_num].second.second))); 
                    
                }
            }
            if (distance_between_right <= tree_state.read_distance_limit) {
                //Combine it right
                to_erase.emplace_back(cluster_head);
                if (combined_right[read_num].first.first == std::numeric_limits<size_t>::max()){
                    combined_right[read_num] = make_pair(cluster_head, make_pair(updated_left, updated_right));
                } else {
                    to_erase.emplace_back(combined_right[read_num].first);
                    combined_right[read_num] =make_pair(
                         make_pair(read_num, 
                                   tree_state.read_union_find.at(read_num).union_groups(cluster_num, 
                                                                                        combined_right[read_num].first.second)), 
                         make_pair(std::min(updated_left, combined_right[read_num].second.first), 
                                   std::min(updated_right, combined_right[read_num].second.second))); 
                    
                }
            }
            if (tree_state.fragment_distance_limit != 0 &&
                  distance_between_left_fragment <= tree_state.fragment_distance_limit) {
                //Combine the fragment
                if (combined_fragment_left != std::numeric_limits<size_t>::max()) {
                    combined_fragment_left = tree_state.fragment_union_find.union_groups(combined_fragment_left, 
                            cluster_num + tree_state.seed_count_prefix_sum[read_num]);
                } else {
                    combined_fragment_left = cluster_num + tree_state.seed_count_prefix_sum[read_num];
                }
            }
            if (tree_state.fragment_distance_limit != 0 &&
                  distance_between_right_fragment <= tree_state.fragment_distance_limit) {
                //Combine the fragment
                if (combined_fragment_right != std::numeric_limits<size_t>::max()) {
                    combined_fragment_right = tree_state.fragment_union_find.union_groups(combined_fragment_right, 
                            cluster_num + tree_state.seed_count_prefix_sum[read_num]);
                } else {
                    combined_fragment_right = cluster_num + tree_state.seed_count_prefix_sum[read_num];
                }
            }
         }
        for (pair<size_t, size_t>& cluster_head : to_erase) {
            child_clusters.read_cluster_heads.erase(cluster_head);
        }
        //Add new clusters that were combined
        for (pair<pair<size_t, size_t>,pair<size_t, size_t>>& cluster : combined_left) {
            if (cluster.first.first != std::numeric_limits<size_t>::max()){
                child_clusters.read_cluster_heads.emplace(cluster.first);
                tree_state.all_seeds->at(cluster.first.first)->at(cluster.first.second).distance_left = cluster.second.first;
                tree_state.all_seeds->at(cluster.first.first)->at(cluster.first.second).distance_right = cluster.second.second;
            }
        }
        for (pair<pair<size_t, size_t>, pair<size_t, size_t>>& cluster : combined_right) {
            if (cluster.first.first != std::numeric_limits<size_t>::max()){
                child_clusters.read_cluster_heads.emplace(cluster.first);
                tree_state.all_seeds->at(cluster.first.first)->at(cluster.first.second).distance_left = cluster.second.first;
                tree_state.all_seeds->at(cluster.first.first)->at(cluster.first.second).distance_right = cluster.second.second;
            }
        }
    };


    net_handle_t& chain_handle = chain_clusters.containing_net_handle;
    NodeClusters& child_clusters = tree_state.all_node_clusters[current_child.child_index1];
    
    //Skip this child if its seeds are all too far away
    bool skip_snarl = false;
    if (child_clusters.fragment_best_left > (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit) &&  
        child_clusters.fragment_best_right > (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit)) {
        skip_snarl = true;
    } else {
        //See if this clusters of this child can be combined with each other
        //Also updates the minimum distances to include loops in the chain
        //This only matters for snarls, since any path would have to pass through a node anyway
        update_distances_on_same_child(child_clusters);
    }
#ifdef DEBUG_CLUSTER
            cerr << "At child " << distance_index.net_handle_as_string(current_child.child_handle) << endl;
#endif

    /*
    Get a bunch of distances from the current child that will be used to calculate distance
    from the last child
    */
    
    
    //The distance from the right side of the last child to the left side of this child 
    //(relative to the orientation of the chain
    //If this is a looping chain, then find the distance normally. Otherwise use the prefix sums
    size_t distance_from_last_child_to_current_child = std::numeric_limits<size_t>::max();
    if (!is_first_child) {
        //If this isn't the first child we're looking at
        if ( last_chain_component_end == child_clusters.chain_component_start) {
            //If this child is in the same component as the last one
            if (last_length == std::numeric_limits<size_t>::max() && last_chain_component_end ) {
                //If the last length is infinite, then is must be a snarl that is not start-end reachable, so the distance
                //from the last child is the same as the distance from the start of the chain (the start of this compnent)
                distance_from_last_child_to_current_child = child_clusters.prefix_sum_value;
            } else {
                size_t distance_from_chain_start_to_last_node = SnarlDistanceIndex::sum(last_prefix_sum,last_length);
                distance_from_last_child_to_current_child = SnarlDistanceIndex::minus(child_clusters.prefix_sum_value, 
                                                distance_from_chain_start_to_last_node); 
            }
        }
    }
    
    
    //The distance from the right side of the last child to the right side of this child, which is
    //the distance we need to update the chain clusters to the end of this child
    //This isn't quite right for the first thing in the chain but it doesn't matter because it only
    //gets added to chain clusters
    //If it gets calculated, it is the distance from the last child to the start of this child snarl + the length of the child snarl
    size_t distance_from_last_child_to_current_end = 
            distance_from_last_child_to_current_child == std::numeric_limits<size_t>::max() 
                    ? std::numeric_limits<size_t>::max() : 
            (last_child.child_handle == current_child.child_handle ? 0 
                : SnarlDistanceIndex::sum(distance_from_last_child_to_current_child, 
                                           child_clusters.node_length));
    
    //The distance to add to get to the end of the chain. Only matters if this is the last thing in the chain
    //The distances will include the distance to the end of a trivial chain,
    //so we can't rely on distance_in_parent to know when the distance should be 0
    
    size_t distance_from_current_end_to_end_of_chain;
    if (!is_last_child || skip_distances_to_ends) {
        //If this isn't the last child in the chain, then we only want the distance to the end of the current child
    
        distance_from_current_end_to_end_of_chain = 0;
    } else if (SnarlDistanceIndex::get_record_offset(current_child.child_handle) == SnarlDistanceIndex::get_record_offset(chain_clusters.end_in)) {
        //If this is the last node in the chain
        if (chain_clusters.chain_component_end != child_clusters.chain_component_end) { 
            //If they aren't in the same component
            distance_from_current_end_to_end_of_chain = std::numeric_limits<size_t>::max();
        } else {
            distance_from_current_end_to_end_of_chain = 0;
        }
    } else if (chain_clusters.is_looping_chain) {
        //TODO: I think I should be able to do this without the distance index but none of our graphs so far have loops 
        //      so I'm not going to bother
        //If it's a looping chain then use the distance index
        distance_from_current_end_to_end_of_chain = distance_index.distance_in_parent(chain_handle, chain_clusters.end_in, 
                 current_child.child_handle);
    } else if (child_clusters.node_length == std::numeric_limits<size_t>::max() ) {
            //If the node length is infinite, then it is a snarl that isn't start-end connected, so the start
            //and end of the snarl are in different components of the chain. Since it reached here, the end
            //node of the snarl is in the same component as the end of the chain, so the distance to the
            //end of the chain is just the length of the last component of the chain, which is
            //chain_clusters.node_length
            distance_from_current_end_to_end_of_chain = chain_clusters.node_length;
    
    } else {
        distance_from_current_end_to_end_of_chain = SnarlDistanceIndex::minus(chain_clusters.node_length, 
                    SnarlDistanceIndex::sum(child_clusters.prefix_sum_value, child_clusters.node_length));
    
    }

#ifdef DEBUG_CLUSTER
cerr << "\tDistance from last child to this one: " << distance_from_last_child_to_current_child << endl;
cerr << "\tDistance from start of chain to the left side of this one: " << (child_clusters.chain_component_start != 0
                                                             ? std::numeric_limits<size_t>::max() : child_clusters.prefix_sum_value) << endl;
cerr << "\tDistance from the last child to the right side of this one: " << distance_from_last_child_to_current_end << endl;
cerr << "\tDistance to get to the end of the chain: " << distance_from_current_end_to_end_of_chain << endl;
#endif

    //Clusters to remove from the chain because they got combined
    vector<pair<size_t, size_t>> to_erase;

    //And new clusters to add that didn't get combined
    vector<pair<pair<size_t, size_t>, pair<size_t, size_t>>> to_add;
    
    //There is at most one new cluster per read
    pair<ClusterIndices, ClusterIndices> new_cluster_by_read;
    //And one new fragment cluster
    size_t new_cluster_head_fragment = std::numeric_limits<size_t>::max();
    
    bool child_is_reversed = child_clusters.is_reversed_in_parent;
    
    //Remember the current best chain distances, and reset them to inf since we need to update them
    size_t old_best_right = std::move(chain_clusters.fragment_best_right);
    chain_clusters.fragment_best_right = std::numeric_limits<size_t>::max(); 
    pair<size_t, size_t> old_best_right_by_read = std::move(chain_clusters.read_best_right);
    chain_clusters.read_best_right = std::make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
    
    
    if (last_child.child_handle != current_child.child_handle &&
        SnarlDistanceIndex::sum(distance_from_last_child_to_current_child, old_best_right) 
          > (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit)) {
#ifdef DEBUG_CLUSTER
        cerr << "This child is too far away from the last one to cluster anything" << endl;
#endif
        if (!skip_distances_to_ends) {
            //If we care about the distances to the ends, then remember which clusters might be reachable from the
            //left side of the chain
            for (auto& cluster_head : chain_clusters.read_cluster_heads) {
                //For each of the chain clusters
                pair<size_t, size_t> dists (tree_state.all_seeds->at(cluster_head.first)->at(cluster_head.second).distance_left,
                                            tree_state.all_seeds->at(cluster_head.first)->at(cluster_head.second).distance_right);
                if (dists.first <= (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit)) {
                    //If this cluster can be clustered outside of the chain, remember to add it back
                    cluster_heads_to_add_again.emplace_back(cluster_head.first, cluster_head.second, 
                                                            dists.first, std::numeric_limits<size_t>::max());
                }
            }
        }
    
        //Now clear the chain's list of clusters
        chain_clusters.read_cluster_heads.clear();
        
        //If the current child snarl has combinable clusters, add them to the chain
        if (!skip_snarl) {
            for (auto& cluster_head : child_clusters.read_cluster_heads) {
                //Add the clusters from this child to the chain
                size_t read_num = cluster_head.first;
                pair<size_t, size_t> dists (tree_state.all_seeds->at(read_num)->at(cluster_head.second).distance_left,
                                           tree_state.all_seeds->at(read_num)->at(cluster_head.second).distance_right);
                size_t dist_left = child_clusters.is_reversed_in_parent ? dists.second : dists.first;
                size_t dist_right = child_clusters.is_reversed_in_parent ? dists.first : dists.second;
    
                //Distances to the start of the chain, and the end of this node
                //If this is the last thing in the chain, then the distance to the end of the chain
                //If the snarl is isn't in the first component of the chain, then the left distance is infinite
                pair<size_t, size_t> new_distances = make_pair(
                         child_clusters.chain_component_start != 0 ? std::numeric_limits<size_t>::max() 
                                                                   : SnarlDistanceIndex::sum(dist_left, child_clusters.prefix_sum_value),
                         SnarlDistanceIndex::sum(dist_right, distance_from_current_end_to_end_of_chain));
    
                //Add this to the chain
                chain_clusters.read_cluster_heads.emplace(cluster_head); 
                tree_state.all_seeds->at(cluster_head.first)->at(cluster_head.second).distance_left = new_distances.first;
                tree_state.all_seeds->at(cluster_head.first)->at(cluster_head.second).distance_right = new_distances.second;
                //And update the best distances
                if (!found_first_node) {
                    chain_clusters.fragment_best_left = std::min(chain_clusters.fragment_best_left, new_distances.first);
                }
                if (! (read_num == 0 ? found_first_node_by_read.first : found_first_node_by_read.second)) {
                    (read_num == 0 ? chain_clusters.read_best_left.first : chain_clusters.read_best_left.second) = 
                        std::min((read_num == 0 ? chain_clusters.read_best_left.first : chain_clusters.read_best_left.second), new_distances.first);
                }
                chain_clusters.fragment_best_right = std::min(chain_clusters.fragment_best_right, new_distances.second); 
                (read_num == 0 ? chain_clusters.read_best_right.first : chain_clusters.read_best_right.second) = 
                        std::min((read_num == 0 ? chain_clusters.read_best_right.first : chain_clusters.read_best_right.second), new_distances.second);
            }
        }
    
    
    } else if (!skip_snarl) {
        //Otherwise, check to see if anything on the current child can be combined with 
        //anything in the chain thus far
    
    
        /**First, go through the clusters of the current child and see what can be combined
        */
    
        for (auto& child_cluster_head : child_clusters.read_cluster_heads) {
            //Go through all clusters of the current child and see if they can be combined with anything on the chain
            const size_t read_num = child_cluster_head.first;
            const size_t cluster_num = child_cluster_head.second;
            ClusterIndices& new_cluster = read_num == 0 ?  new_cluster_by_read.first : new_cluster_by_read.second;
            pair<size_t, size_t> dists (tree_state.all_seeds->at(read_num)->at(cluster_num).distance_left,
                                       tree_state.all_seeds->at(read_num)->at(cluster_num).distance_right);
            const size_t distance_left = child_is_reversed ? dists.second : dists.first;
            const size_t distance_right = child_is_reversed ? dists.first : dists.second;
            //Distance between this cluster and a cluster on the same read from the previous child
            size_t distance_between = SnarlDistanceIndex::minus(
                     SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(distance_left, 
                                              distance_from_last_child_to_current_child),
                                              (read_num == 0 ? old_best_right_by_read.first : old_best_right_by_read.second)),
                      1);
            //Distance between this cluster and any cluster on the previous child
            size_t fragment_distance_between = SnarlDistanceIndex::minus(
                      SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(distance_left, 
                                               distance_from_last_child_to_current_child),  
                                               old_best_right),
                      1);
    
            //The new distances from this child to the start of the chain and the end of this child
            pair<size_t, size_t> new_distances = make_pair(
                     child_clusters.chain_component_start != 0 ? std::numeric_limits<size_t>::max() 
                                                               : SnarlDistanceIndex::sum(distance_left, child_clusters.prefix_sum_value),
                    SnarlDistanceIndex::sum(distance_right, distance_from_current_end_to_end_of_chain)); 
    
            if (distance_between <= tree_state.read_distance_limit) {
#ifdef DEBUG_CLUSTER
                cerr << "\t\tCombine child cluster " << read_num << ":" << cluster_num << endl;
#endif
                //If this cluster can be merged with anything on the chain
                if (new_cluster.read_num == std::numeric_limits<size_t>::max()){
                    //If nothing is in the combined cluster yet, this is the new combined cluster
                    new_cluster.read_num =child_cluster_head.first;
                    new_cluster.cluster_num =child_cluster_head.second;
                    new_cluster.distance_left = new_distances.first;
                    new_cluster.distance_right = new_distances.second;
                } else {
                    //Otherwise, remember to forget about the old cluster
                    //Union the two clusters
                    new_cluster.cluster_num = tree_state.read_union_find.at(read_num).union_groups(cluster_num, new_cluster.cluster_num);
                    //And find the best distances of the two
                    new_cluster.distance_left = std::min(new_cluster.distance_left, new_distances.first); 
                    new_cluster.distance_right = std::min(new_cluster.distance_right, new_distances.second); 
                }
            } else {
                //If it didn't get combined, remember to add it at the end
                size_t distance_limit = tree_state.fragment_distance_limit == 0
                                         ? tree_state.read_distance_limit 
                                         : tree_state.fragment_distance_limit; 
                if (new_distances.first <= distance_limit || new_distances.second <= distance_limit){
                    //But only if the distances are small enough
                    to_add.emplace_back(make_pair(read_num, cluster_num), new_distances);
                }
            }
            //If we can combine the fragments
            if (tree_state.fragment_distance_limit != 0 && fragment_distance_between <= tree_state.fragment_distance_limit){
                if (new_cluster_head_fragment != std::numeric_limits<size_t>::max()) {
                    new_cluster_head_fragment = tree_state.fragment_union_find.union_groups(new_cluster_head_fragment, 
                            cluster_num + tree_state.seed_count_prefix_sum[read_num]);
                } else {
                    new_cluster_head_fragment = cluster_num + tree_state.seed_count_prefix_sum[read_num];
                }
            }
        
        
            //Update the best distances
            if (!found_first_node) {
                chain_clusters.fragment_best_left = std::min(chain_clusters.fragment_best_left, new_distances.first);
            }
            if (!(read_num == 0 ? found_first_node_by_read.first : found_first_node_by_read.second)) {
                (read_num == 0 ? chain_clusters.read_best_left.first : chain_clusters.read_best_left.second) = std::min((read_num == 0 ? chain_clusters.read_best_left.first : chain_clusters.read_best_left.second), new_distances.first);
            }
            chain_clusters.fragment_best_right = std::min(chain_clusters.fragment_best_right, new_distances.second); 
            (read_num == 0 ? chain_clusters.read_best_right.first : chain_clusters.read_best_right.second) = std::min((read_num == 0 ? chain_clusters.read_best_right.first : chain_clusters.read_best_right.second), new_distances.second);
        }
        
        
        
        /**Next, go through the clusters on the chain up to this point and see if anything can
           be combined with the clusters on the child
           */
        for (auto& chain_cluster_head : chain_clusters.read_cluster_heads) {
            //Each has distances up to the previous node
        
            const size_t read_num = chain_cluster_head.first;
            const size_t cluster_num = chain_cluster_head.second;

            ClusterIndices& new_cluster = read_num == 0 ?  new_cluster_by_read.first : new_cluster_by_read.second;

            //The distances for the chain cluster
            pair<size_t, size_t> chain_cluster_distances (tree_state.all_seeds->at(read_num)->at(cluster_num).distance_left,
                                        tree_state.all_seeds->at(read_num)->at(cluster_num).distance_right);
        
            //Best distance to the left side (relative to the chain) of the current child
            const size_t current_distance_left =  (read_num == 0 ? child_clusters.read_best_left.first : child_clusters.read_best_left.second);
            const size_t current_fragment_distance_left =  child_clusters.fragment_best_left;

            size_t distance_between = SnarlDistanceIndex::minus(
                    SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(chain_cluster_distances.second, 
                                             distance_from_last_child_to_current_child), 
                                             current_distance_left),
                    1);
            size_t distance_between_fragment = SnarlDistanceIndex::minus(
                    SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(chain_cluster_distances.second, 
                                             distance_from_last_child_to_current_child), 
                                             current_fragment_distance_left),
                    1);
            pair<size_t, size_t> new_distances = make_pair(
                chain_cluster_distances.first,
                SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(chain_cluster_distances.second, 
                                         distance_from_last_child_to_current_end),
                                         distance_from_current_end_to_end_of_chain));
        
                 
            if (distance_between <= tree_state.read_distance_limit) {
#ifdef DEBUG_CLUSTER
                cerr << "\t\tCombine chain cluster " << read_num << ":" << cluster_num << endl;
#endif
                //If we can union the reads
                if ( new_cluster.read_num == std::numeric_limits<size_t>::max()) {
                    new_cluster.read_num =chain_cluster_head.first;
                    new_cluster.cluster_num =chain_cluster_head.second;
                    new_cluster.distance_left = new_distances.first;
                    new_cluster.distance_right = new_distances.second;
                } else {
                    //Union the two clusters
                    new_cluster.cluster_num = tree_state.read_union_find.at(read_num).union_groups(cluster_num, 
                                                                      new_cluster.cluster_num);
                    //Find the best distances of the two
                    new_cluster.distance_left = std::min( new_cluster.distance_left, new_distances.first);
                    new_cluster.distance_right = std::min(new_cluster.distance_right, new_distances.second); 
                   
                }
                //Remember to erase the combined cluster. The new cluster head will be added later
                to_erase.emplace_back(read_num, cluster_num);
            } else {
                to_add.emplace_back(make_pair(read_num, cluster_num), new_distances);
            }
            if (tree_state.fragment_distance_limit != 0 && distance_between_fragment <= tree_state.fragment_distance_limit) {
                //If we can union the fragments
                if (new_cluster_head_fragment != std::numeric_limits<size_t>::max()) {
                    new_cluster_head_fragment = tree_state.fragment_union_find.union_groups(new_cluster_head_fragment, 
                            cluster_num + tree_state.seed_count_prefix_sum[read_num]);
                } else {
                    new_cluster_head_fragment = cluster_num + tree_state.seed_count_prefix_sum[read_num];
                }
            }
            chain_clusters.fragment_best_right = std::min(chain_clusters.fragment_best_right, new_distances.second); 
            (read_num == 0 ? chain_clusters.read_best_right.first : chain_clusters.read_best_right.second) = std::min((read_num == 0 ? chain_clusters.read_best_right.first : chain_clusters.read_best_right.second), new_distances.second);
        }
    
        //Remove clusters that got combined
        for (pair<size_t, size_t>& cluster_head : to_erase) {
            chain_clusters.read_cluster_heads.erase(cluster_head);
        }
        //Add new clusters that weren't combined
        for (pair<pair<size_t, size_t>, pair<size_t, size_t>>& cluster : to_add) {
            chain_clusters.read_cluster_heads.emplace(cluster.first);
            tree_state.all_seeds->at(cluster.first.first)->at(cluster.first.second).distance_left = cluster.second.first;
            tree_state.all_seeds->at(cluster.first.first)->at(cluster.first.second).distance_right = cluster.second.second;
        }
        //Add new clusters that were combined
        if (new_cluster_by_read.first.read_num != std::numeric_limits<size_t>::max()){
            chain_clusters.read_cluster_heads.emplace(new_cluster_by_read.first.read_num, new_cluster_by_read.first.cluster_num);
            tree_state.all_seeds->at(new_cluster_by_read.first.read_num)->at(new_cluster_by_read.first.cluster_num).distance_left = new_cluster_by_read.first.distance_left;
            tree_state.all_seeds->at(new_cluster_by_read.first.read_num)->at(new_cluster_by_read.first.cluster_num).distance_right = new_cluster_by_read.first.distance_right;
        }
        if (new_cluster_by_read.second.read_num != std::numeric_limits<size_t>::max()){
            chain_clusters.read_cluster_heads.emplace(new_cluster_by_read.second.read_num, new_cluster_by_read.second.cluster_num);
            tree_state.all_seeds->at(new_cluster_by_read.second.read_num)->at(new_cluster_by_read.second.cluster_num).distance_left = new_cluster_by_read.second.distance_left;
            tree_state.all_seeds->at(new_cluster_by_read.second.read_num)->at(new_cluster_by_read.second.cluster_num).distance_right = new_cluster_by_read.second.distance_right;
        }
    } else {
#ifdef DEBUG_CLUSTER
       cerr << "The snarl was too big to combine with anything, so go through the current clusters of the chain and add distances" << endl;
       cerr << distance_from_last_child_to_current_end << " and " << distance_from_current_end_to_end_of_chain <<  " to the distances right" << endl;
#endif
        //If this was a snarl that we skipped because the distances were too big, then we need to update the distances
        chain_clusters.fragment_best_right = SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(old_best_right,
                                                                      distance_from_last_child_to_current_end),
                                                                      distance_from_current_end_to_end_of_chain);
        chain_clusters.read_best_right.first = SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(old_best_right_by_read.first,
                                                                  distance_from_last_child_to_current_end),
                                                                  distance_from_current_end_to_end_of_chain);
        chain_clusters.read_best_right.second = SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(old_best_right_by_read.second,
                                                                  distance_from_last_child_to_current_end),
                                                                  distance_from_current_end_to_end_of_chain);
        for (pair<size_t, size_t> cluster_head : chain_clusters.read_cluster_heads) {
            tree_state.all_seeds->at(cluster_head.first)->at(cluster_head.second).distance_right = SnarlDistanceIndex::sum(SnarlDistanceIndex::sum(
                                                        tree_state.all_seeds->at(cluster_head.first)->at(cluster_head.second).distance_right,
                                                        distance_from_last_child_to_current_end),
                                                        distance_from_current_end_to_end_of_chain);
        }
    }
    
    
    //Update the last node we saw to this one
    last_child = current_child;
    last_prefix_sum = child_clusters.prefix_sum_value;
    last_length = child_clusters.node_length; //The length of this snarl
    last_chain_component_end = child_clusters.chain_component_end;//The component of the end node of this snarl
}

//Cluster the root
//all children of the root will be in tree_state.root_children
//This is basically cluster_one_snarl except the snarl is the root, which has no boundary nodes
void SnarlDistanceIndexClusterer::cluster_root(TreeState& tree_state) const { 
#ifdef DEBUG_CLUSTER
    cerr << "Finding clusters on the root with " << tree_state.root_children.size() << " children" << endl;
#endif
    if (tree_state.root_children.size() == 0) {
        return;
    }

    //Keep track of all clusters on the root
    NodeClusters root_clusters(distance_index.get_root(), tree_state.all_seeds->size(),
                               tree_state.seed_count_prefix_sum.back(), distance_index);

    //Remember old distances
    vector<pair<size_t, size_t>> child_distances (tree_state.seed_count_prefix_sum.back(), 
                                                  make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));

 
    //Sort the root children by parent, the order of the children doesn't matter
    std::sort(tree_state.root_children.begin(), tree_state.root_children.end());   

    //Go through the list of parent child pairs. Once we reach a new parent, cluster all children found up to this point
    net_handle_t current_parent = tree_state.all_node_clusters[tree_state.root_children.front().first].containing_net_handle;
    vector<size_t> children;
    for (size_t root_child_i = 0 ; root_child_i < tree_state.root_children.size() ; root_child_i++) {
        const pair<size_t, size_t>& parent_to_child = tree_state.root_children[root_child_i];
        net_handle_t& parent = tree_state.all_node_clusters[parent_to_child.first].containing_net_handle;

        if (current_parent == parent || root_child_i == 0) {
            children.emplace_back(parent_to_child.second);
        }
        if (current_parent != parent || root_child_i == tree_state.root_children.size()-1) {
#ifdef DEBUG_CLUSTER
            cerr << "Clustering root snarl " << distance_index.net_handle_as_string(parent) << " with " << children.size() << " chidlren" << endl;
#endif

            if (children.size() > 0) {

                for (size_t i = 0; i < children.size() ; i++) {
                    //Go through each child node of the netgraph

                    NodeClusters& child_clusters = tree_state.all_node_clusters[children[i]];
                    for (const pair<size_t, size_t>& head : child_clusters.read_cluster_heads) {
                        child_distances[head.second + tree_state.seed_count_prefix_sum[head.first]] = 
                            make_pair(tree_state.all_seeds->at(head.first)->at(head.second).distance_left,
                                     tree_state.all_seeds->at(head.first)->at(head.second).distance_right);
                    }

                    for (size_t j = 0 ; j <= i ; j++){
                        //Go through other child net graph nodes up to and including i

                        //Get the other node and its clusters
                        NodeClusters& child_clusters_j = tree_state.all_node_clusters[children[j]];



                        compare_and_combine_cluster_on_child_structures(tree_state, child_clusters,
                                    child_clusters_j, root_clusters, child_distances, true, false);

                    }
                }
            }
            current_parent = parent;
            children.clear();
            children.emplace_back(parent_to_child.second);
        }

    }
#ifdef DEBUG_CLUSTER
    cerr << "\tFound clusters on the root" << endl;
    for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
        cerr << "\t for read num " << read_num << endl;
        for (pair<size_t,size_t> c : root_clusters.read_cluster_heads) {
            if (c.first == read_num) {
                cerr << "\t\t" << c.first << ":"<<c.second << ":  ";
                for (size_t x = 0 ; x < tree_state.all_seeds->at(c.first)->size() ; x++) {
                    if (tree_state.read_union_find[c.first].find_group(x) == c.second) {
                        cerr << tree_state.all_seeds->at(c.first)->at(x).pos << " ";
                    }
                }
                cerr << endl;
            }
        }
    }

    for (pair<size_t, size_t> group_id : root_clusters.read_cluster_heads) {
        assert (group_id.second == tree_state.read_union_find[group_id.first].find_group(group_id.second));
    }
#endif
}


//Cluster all the seeds on a node or chain of only seeds
//Seeds are assumed to be sorted
//Since the seeds can be linearly arranged, they can be clustered just by walking along an ordered list of seeds
template<typename SeedIndex>
void SnarlDistanceIndexClusterer::cluster_seeds_on_linear_structure(TreeState& tree_state, NodeClusters& node_clusters, 
    const typename vector<SeedIndex>::iterator& range_start, const typename vector<SeedIndex>::iterator& range_end,
    size_t structure_length, std::function<tuple<size_t, size_t, size_t>(const SeedIndex&)>& get_offset_from_seed_index, 
    bool skip_distances_to_ends) const{
    if (range_start == range_end) {
        return;
    }
#ifdef DEBUG_CLUSTER
        cerr << "Cluster " << (range_end - range_start) << " seeds on a single structure " << distance_index.net_handle_as_string(node_clusters.containing_net_handle) << endl;
        cerr << "\t with node length " << structure_length << endl;
#endif

    if (tree_state.read_distance_limit >= structure_length) {
        //If the limit is greater than the node length, then all the
        //seeds on this node must be in the same cluster

        //The cluster heads of the new cluster
        size_t fragment_group_id = std::numeric_limits<size_t>::max();
        pair<size_t, size_t> group_ids (std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());

        for (typename vector<SeedIndex>::iterator iter = range_start ; iter != range_end ; ++iter) {
            //Go through all seeds in the range
            const SeedIndex& seed_index_template = *iter;

            //For each seed on this node, add it to the cluster

            //Get the index of the seed 
            tuple<size_t, size_t, size_t> seed_index = get_offset_from_seed_index(seed_index_template);
            size_t read_num = std::get<0>(seed_index);
            size_t seed_i = std::get<1>(seed_index);

            //And the distances for this seed
            size_t dist_left = std::get<2>(seed_index);
            size_t dist_right = structure_length + 1 - dist_left;

            //Find the new best distances for anything in this cluster
            //Since we're traversing the seeds in order, the best left and right will be the first and last 
            //ones we see for each read
            if ((read_num == 0 ? node_clusters.read_best_left.first : node_clusters.read_best_left.second) == std::numeric_limits<size_t>::max()) {
                //Only update the best left if it hasn't been set yet
                (read_num == 0 ? node_clusters.read_best_left.first : node_clusters.read_best_left.second) = dist_left;
            }
            //Best right will always be the most recent thing seen
            (read_num == 0 ? node_clusters.read_best_right.first : node_clusters.read_best_right.second) = dist_right;
            node_clusters.fragment_best_right = dist_right;


            //Put this seed in the cluster for the node
            if ((read_num == 0 ? group_ids.first : group_ids.second) == std::numeric_limits<size_t>::max()) {
                (read_num == 0 ? group_ids.first : group_ids.second) = seed_i;
            } else {
                (read_num == 0 ? group_ids.first : group_ids.second) = tree_state.read_union_find[read_num].union_groups((read_num == 0 ? group_ids.first : group_ids.second), seed_i);
            }
            if (tree_state.fragment_distance_limit != 0 ) {
                if (fragment_group_id == std::numeric_limits<size_t>::max() ) {
                    fragment_group_id = fragment_group_id = seed_i + tree_state.seed_count_prefix_sum[read_num];
                } else {
                    fragment_group_id = tree_state.fragment_union_find.union_groups(
                        fragment_group_id, seed_i + tree_state.seed_count_prefix_sum[read_num]);
                }
            }
        }
        if (!skip_distances_to_ends) {
            node_clusters.fragment_best_left = std::get<2>(get_offset_from_seed_index(*range_start));

            //Record the new cluster
            for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++ ) {
                if ((read_num == 0 ? group_ids.first : group_ids.second) != std::numeric_limits<size_t>::max()) {
                    size_t group_id = (read_num == 0 ? group_ids.first : group_ids.second);
                    node_clusters.read_cluster_heads.emplace(read_num, group_id);
                    tree_state.all_seeds->at(read_num)->at(group_id).distance_left = (read_num == 0 ? node_clusters.read_best_left.first : node_clusters.read_best_left.second);
                    tree_state.all_seeds->at(read_num)->at(group_id).distance_right = (read_num == 0 ? node_clusters.read_best_right.first : node_clusters.read_best_right.second);
                }

            }
        }
        
#ifdef DEBUG_CLUSTER
        cerr << "\t" << distance_index.net_handle_as_string(node_clusters.containing_net_handle) <<  " is shorter than the distance limit so just one cluster" << endl;

#endif
        return;
    }


    //The seeds may form multiple clusters on the node
    //Walk through a sorted list of seeds and split into clusters
        
    //Offset of the first seed for the cluster we're building
    //One for each read
    pair<size_t, size_t> read_first_offset (std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
    //Offset of the latest seed for the cluster we're building
    pair<size_t, size_t> read_last_offset (std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
    //And the same for the fragment 
    size_t fragment_last_offset = std::numeric_limits<size_t>::max();
    size_t fragment_last_cluster = std::numeric_limits<size_t>::max();

    pair<size_t, size_t> read_last_cluster (std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());


    for (typename vector<SeedIndex>::iterator iter = range_start ; iter != range_end ; ++iter) {
        const SeedIndex& seed_index_template = *iter;
        //<index of read, index of seed, offset of seed>
        tuple<size_t, size_t, size_t> s = get_offset_from_seed_index(seed_index_template);

        //For each seed, in order of position in the node,
        //see if it belongs to a new read/fragment cluster - if it is
        //close enough to the previous seed
        size_t read_num = std::get<0>(s);
        size_t seed_num = std::get<1>(s);
        size_t offset = std::get<2>(s);

        //First and last offset and last cluster head for this read
        size_t& first_offset = read_num == 0 ? read_first_offset.first : read_first_offset.second;
        size_t& last_offset = read_num == 0 ? read_last_offset.first : read_last_offset.second;
        size_t& last_cluster = read_num == 0 ? read_last_cluster.first : read_last_cluster.second;


        if (first_offset == std::numeric_limits<size_t>::max()) {
            //If this is the first seed we've seen of this read
            first_offset = offset;
            (read_num == 0 ? node_clusters.read_best_left.first : node_clusters.read_best_left.second) = offset;
        }

        if (last_offset != std::numeric_limits<size_t>::max() &&
            offset - last_offset <= tree_state.read_distance_limit) {
            //If this seed is in the same read cluster as the previous one,
            //union them

            last_cluster = tree_state.read_union_find[read_num].union_groups(seed_num, last_cluster);
            last_offset = offset;

            if (tree_state.fragment_distance_limit != 0) {
                //If we are also clustering paired end reads by fragment distance,
                //cluster these together
                fragment_last_cluster = tree_state.fragment_union_find.union_groups(seed_num+tree_state.seed_count_prefix_sum[read_num], fragment_last_cluster);
            }
        } else {
            //This becomes a new read cluster
            if (!skip_distances_to_ends && last_cluster != std::numeric_limits<size_t>::max()) {
                //Record the previous cluster
                node_clusters.read_cluster_heads.emplace(read_num, last_cluster);
                tree_state.all_seeds->at(read_num)->at(last_cluster).distance_left = first_offset;
                tree_state.all_seeds->at(read_num)->at(last_cluster).distance_right =  structure_length - last_offset + 1;
            }
            last_cluster = seed_num;
            first_offset = offset;
            last_offset = offset;
            if (tree_state.fragment_distance_limit != 0) {
                if (fragment_last_offset != std::numeric_limits<size_t>::max() &&
                    offset - fragment_last_offset <= tree_state.fragment_distance_limit) {
                    //If this is a new read cluster but the same fragment cluster
                    fragment_last_cluster = tree_state.fragment_union_find.union_groups(seed_num+tree_state.seed_count_prefix_sum[read_num], fragment_last_cluster);

                } else {
                    //If this is a new fragment cluster as well
                    fragment_last_cluster = seed_num+tree_state.seed_count_prefix_sum[read_num];
                }
            }
        }
        fragment_last_offset = offset;
    }
    if (!skip_distances_to_ends) {
        //If we want to remember, record the best distances to the ends of this structure
        for (size_t i = 0 ; i < tree_state.all_seeds->size() ; i++) {
            if ((i == 0 ? read_last_cluster.first : read_last_cluster.second) != std::numeric_limits<size_t>::max()) {
                node_clusters.read_cluster_heads.emplace(i, (i == 0 ? read_last_cluster.first : read_last_cluster.second));
                tree_state.all_seeds->at(i)->at((i == 0 ? read_last_cluster.first : read_last_cluster.second)).distance_left = (i == 0 ? read_first_offset.first : read_first_offset.second);
                tree_state.all_seeds->at(i)->at((i == 0 ? read_last_cluster.first : read_last_cluster.second)).distance_right = structure_length-(i == 0 ? read_last_offset.first : read_last_offset.second)+1;

                if (i == 0) {
                    node_clusters.read_best_right.first = structure_length-(i == 0 ? read_last_offset.first : read_last_offset.second) +1;
                } else {
                    node_clusters.read_best_right.second = structure_length-(i == 0 ? read_last_offset.first : read_last_offset.second)+1;
                }
            }
        }

        //Get the best left and right values of the node from the first and last seeds
        node_clusters.fragment_best_left = std::get<2>(get_offset_from_seed_index(*range_start));
        node_clusters.fragment_best_right = structure_length-fragment_last_offset+1;
    }
    return;
}

size_t SnarlDistanceIndexClusterer::distance_between_seeds(Seed& seed1, Seed& seed2, bool stop_at_lowest_common_ancestor) const {

    /*Helper function to walk up the snarl tree
     * Given a net handle, its parent,  and the distances to the start and end of the handle, 
     * update the distances to reach the ends of the parent and update the handle and its parent
     * If the parent is a chain, then the new distances include the boundary nodes of the chain.
     * If it is a snarl, it does not*/
    auto update_distances = [&](net_handle_t& net, net_handle_t& parent, size_t& dist_start, size_t& dist_end) {
#ifdef debug_distances
        cerr << "     Updating distance from node " << distance_index.net_handle_as_string(net) << " at parent " << distance_index.net_handle_as_string(parent) << " from " << dist_start << " " << dist_end << endl;
#endif

        if (distance_index.is_trivial_chain(parent)) {
            //Don't update distances for the trivial chain
            return;
        } else if (distance_index.is_simple_snarl(parent)) {
            //If it's a simple snarl just check if they should be reversed
            if (distance_index.is_reversed_in_parent (net)) {
                size_t tmp = dist_start;
                dist_start = dist_end;
                dist_end = tmp;
            }
            return;
        }

        net_handle_t start_bound = distance_index.get_bound(parent, false, true);
        net_handle_t end_bound = distance_index.get_bound(parent, true, true);

        //The lengths of the start and end nodes of net
        //This is only needed if net is a snarl, since the boundary nodes are not technically part of the snarl
        size_t start_length = distance_index.is_chain(parent) ? distance_index.node_length(start_bound) : 0;
        size_t end_length = distance_index.is_chain(parent) ? distance_index.node_length(end_bound) : 0;

        //Get the distances from the bounds of the parent to the node we're looking at
        size_t distance_start_start = start_bound == net ? 0
                : SnarlDistanceIndex::sum(start_length, distance_index.distance_in_parent(parent, start_bound, distance_index.flip(net), graph));
        size_t distance_start_end = start_bound == distance_index.flip(net) ? 0
                : SnarlDistanceIndex::sum(start_length, distance_index.distance_in_parent(parent, start_bound, net, graph));
        size_t distance_end_start = end_bound == net ? 0
                : SnarlDistanceIndex::sum(end_length, distance_index.distance_in_parent(parent, end_bound, distance_index.flip(net), graph));
        size_t distance_end_end = end_bound == distance_index.flip(net) ? 0
                : SnarlDistanceIndex::sum(end_length, distance_index.distance_in_parent(parent, end_bound, net, graph));

        size_t distance_start = dist_start;
        size_t distance_end = dist_end;

        dist_start = std::min(SnarlDistanceIndex::sum(distance_start_start, distance_start),
                              SnarlDistanceIndex::sum(distance_start_end , distance_end));
        dist_end = std::min(SnarlDistanceIndex::sum(distance_end_start , distance_start),
                            SnarlDistanceIndex::sum(distance_end_end , distance_end));
#ifdef debug_distances
        cerr << "        ...new distances to start and end: " << dist_start << " " << dist_end << endl;
#endif
        return;
    };

    /*
     * Get net handles for the two nodes and the distances from each position to the ends of the handles
     */

    bool has_cached_values1 = seed1.minimizer_cache.record_offset != MIPayload::NO_VALUE;
    bool has_cached_values2 = seed2.minimizer_cache.record_offset != MIPayload::NO_VALUE;
    net_handle_t net1 = has_cached_values1 ? distance_index.get_net_handle_from_values(seed1.minimizer_cache.record_offset, 
                                                             SnarlDistanceIndex::START_END, 
                                                             SnarlDistanceIndex::NODE_HANDLE, 
                                                             seed1.minimizer_cache.node_record_offset)
                                           : distance_index.get_node_net_handle(get_id(seed1.pos));
    net_handle_t net2 = has_cached_values2 ? distance_index.get_net_handle_from_values(seed2.minimizer_cache.record_offset, 
                                                             SnarlDistanceIndex::START_END, 
                                                             SnarlDistanceIndex::NODE_HANDLE, 
                                                             seed2.minimizer_cache.node_record_offset)
                                           : distance_index.get_node_net_handle(get_id(seed2.pos));

    size_t minimum_distance = std::numeric_limits<size_t>::max();
    if (net1 == net2) {
        //If the two positions are on the same node, get the distance between them
        size_t node_length = has_cached_values1 ? seed1.minimizer_cache.node_length 
                                                 : distance_index.node_length(net1);
        size_t distance_to_start1 = is_rev(seed1.pos) ? node_length - get_offset(seed1.pos) : get_offset(seed1.pos) + 1;
        size_t distance_to_end1 =   is_rev(seed1.pos) ? get_offset(seed1.pos) + 1 : node_length - get_offset(seed1.pos);
        size_t distance_to_start2 = is_rev(seed2.pos) ? node_length - get_offset(seed2.pos) : get_offset(seed2.pos) + 1;
        size_t distance_to_end2 =   is_rev(seed2.pos) ? get_offset(seed2.pos) + 1 : node_length - get_offset(seed2.pos);

        if (distance_to_start1 < distance_to_start2) {
            //IF 1 comes before 2
            minimum_distance = SnarlDistanceIndex::minus(SnarlDistanceIndex::sum(distance_to_end1 , distance_to_start2), node_length);
        } else {
            minimum_distance = SnarlDistanceIndex::minus(SnarlDistanceIndex::sum(distance_to_end2 , distance_to_start1), node_length);
        }
        if (stop_at_lowest_common_ancestor) {
            //If we only care about the lowest common ancestor, then return
            return SnarlDistanceIndex::minus(minimum_distance, 1);
        }
        
    }

    /*
     * Since we want to use the minimizer payload, go up one level of the snarl tree here, before using the
     * distance index.
     * Find the parent and the distances to the ends of the parent using the payload
     */

    //Get the parents of the nodes
    net_handle_t parent1;
    //If the grandparent is a root/root snarl, then make it the parent and the node a trivial chain 
    //because they will be clustered here and added to the root instead of being added to the 
    //snarl tree to be clustered
    if (has_cached_values1) {
        if (seed1.minimizer_cache.is_trivial_chain) {
            //If the node is a trivial chain, then the parent is just the node but recorded as a chain in the net handle
            parent1 = distance_index.get_net_handle_from_values (distance_index.get_record_offset(net1),
                                                    SnarlDistanceIndex::START_END,
                                                    SnarlDistanceIndex::CHAIN_HANDLE,
                                                    seed1.minimizer_cache.node_record_offset);
            if (seed1.minimizer_cache.parent_record_offset == 0) {
                //If the parent offset stored in the cache is the root, then this is a trivial chain
                //child of the root not in a root snarl, so remember the root as the parent and the 
                //trivial chain as th enode
                net1 = parent1;
                parent1 = distance_index.get_root();
            } else if (seed1.minimizer_cache.parent_is_root && !seed1.minimizer_cache.parent_is_chain) {
                //If the parent is a root snarl, then the node becomes the trivial chain 
                //and we get the parent root snarl from the cache
                net1 = parent1;
                parent1 = distance_index.get_net_handle_from_values(seed1.minimizer_cache.parent_record_offset,
                                                                   SnarlDistanceIndex::START_END,
                                                                   SnarlDistanceIndex::ROOT_HANDLE);
            }
        } else if (seed1.minimizer_cache.parent_record_offset == 0) {
            //The parent is just the root
            parent1 = distance_index.get_root();
        } else if (seed1.minimizer_cache.parent_is_root && !seed1.minimizer_cache.parent_is_chain) {
            //If the parent is a root snarl
            parent1 = distance_index.get_net_handle_from_values(seed1.minimizer_cache.parent_record_offset,
                                                   SnarlDistanceIndex::START_END,
                                                   SnarlDistanceIndex::ROOT_HANDLE);
        } else {
            //Otherwise the parent is an actual chain and we use the value from the cache
            parent1 = distance_index.get_net_handle_from_values(seed1.minimizer_cache.parent_record_offset,
                                                   SnarlDistanceIndex::START_END,
                                                   SnarlDistanceIndex::CHAIN_HANDLE);
        }
    } else {
        parent1 = distance_index.start_end_traversal_of(distance_index.get_parent(net1));
        if (distance_index.is_trivial_chain(parent1)){
            net_handle_t grandparent = distance_index.get_parent(parent1);
            if (distance_index.is_root(grandparent)){
                net1 = parent1;
                parent1 = distance_index.start_end_traversal_of(grandparent);
            }
        }
    }

    net_handle_t parent2;
    //If the grandparent is a root/root snarl, then make it the parent and the node a trivial chain 
    //because they will be clustered here and added to the root instead of being added to the 
    //snarl tree to be clustered
    if (has_cached_values2) {
        if (seed2.minimizer_cache.is_trivial_chain) {
            //If the node is a trivial chain, then the parent is just the node but recorded as a chain in the net handle
            parent2 = distance_index.get_net_handle_from_values (distance_index.get_record_offset(net2),
                                                    SnarlDistanceIndex::START_END,
                                                    SnarlDistanceIndex::CHAIN_HANDLE,
                                                    seed2.minimizer_cache.node_record_offset);
            if (seed2.minimizer_cache.parent_record_offset == 0) {
                //If the parent offset stored in the cache is the root, then this is a trivial chain
                //child of the root not in a root snarl, so remember the root as the parent and the 
                //trivial chain as th enode
                net2 = parent2;
                parent2 = distance_index.get_root();
            } else if (seed2.minimizer_cache.parent_is_root && !seed2.minimizer_cache.parent_is_chain) {
                //If the parent is a root snarl, then the node becomes the trivial chain 
                //and we get the parent root snarl from the cache
                net2 = parent2;
                parent2 = distance_index.get_net_handle_from_values(seed2.minimizer_cache.parent_record_offset,
                                                                   SnarlDistanceIndex::START_END,
                                                                   SnarlDistanceIndex::ROOT_HANDLE);
            }
        } else if (seed2.minimizer_cache.parent_record_offset == 0) {
            //The parent is just the root
            parent2 = distance_index.get_root();
        } else if (seed2.minimizer_cache.parent_is_root && !seed2.minimizer_cache.parent_is_chain) {
            //If the parent is a root snarl
            parent2 = distance_index.get_net_handle_from_values(seed2.minimizer_cache.parent_record_offset,
                                                   SnarlDistanceIndex::START_END,
                                                   SnarlDistanceIndex::ROOT_HANDLE);
        } else {
            //Otherwise the parent is an actual chain and we use the value from the cache
            parent2 = distance_index.get_net_handle_from_values(seed2.minimizer_cache.parent_record_offset,
                                                   SnarlDistanceIndex::START_END,
                                                   SnarlDistanceIndex::CHAIN_HANDLE);
        }
    } else {
        parent2 = distance_index.start_end_traversal_of(distance_index.get_parent(net2));
        if (distance_index.is_trivial_chain(parent2)){
            net_handle_t grandparent = distance_index.get_parent(parent2);
            if (distance_index.is_root(grandparent)){
                net2 = parent2;
                parent2 = distance_index.start_end_traversal_of(grandparent);
            }
        }
    }



#ifdef debug_distances
        cerr << "Found parents "  << distance_index.net_handle_as_string(parent1) << " and " << distance_index.net_handle_as_string(parent2) << endl;
#endif

    pair<net_handle_t, bool> lowest_ancestor = distance_index.lowest_common_ancestor(parent1, parent2);
    //The lowest common ancestor of the two positions
    net_handle_t common_ancestor = distance_index.start_end_traversal_of(lowest_ancestor.first);

#ifdef debug_distances
        cerr << "Found the lowest common ancestor " << distance_index.net_handle_as_string(common_ancestor) << endl;
#endif

    //These are the distances to the ends of the node, including the position
    size_t node_length1 = has_cached_values1 ? seed1.minimizer_cache.node_length 
                                             : distance_index.minimum_length(net1);
    size_t node_length2 = has_cached_values2 ? seed2.minimizer_cache.node_length 
                                             : distance_index.minimum_length(net2);
    size_t distance_to_start1 = is_rev(seed1.pos) ? node_length1 - get_offset(seed1.pos) : get_offset(seed1.pos) + 1;
    size_t distance_to_end1 =   is_rev(seed1.pos) ? get_offset(seed1.pos) + 1 : node_length1 - get_offset(seed1.pos);
    size_t distance_to_start2 = is_rev(seed2.pos) ? node_length2 - get_offset(seed2.pos) : get_offset(seed2.pos) + 1;
    size_t distance_to_end2 =   is_rev(seed2.pos) ? get_offset(seed2.pos) + 1 : node_length2 - get_offset(seed2.pos);

#ifdef debug_distances
        cerr << "Reached node " << distance_index.net_handle_as_string(net1) << " for position 1" << endl;
        cerr << "   with distances to ends " << distance_to_start1 << " and " << distance_to_end1 << endl;
        cerr << "Reached node " << distance_index.net_handle_as_string(net2) << " for position 2" << endl;
        cerr << "   with distances to ends " << distance_to_start2 << " and " << distance_to_end2 << endl;
#endif
    /* get the distance from the ends of the nodes to the ends of the parent, and update the nodes to their parent*/

    if (distance_index.start_end_traversal_of(parent1) == distance_index.start_end_traversal_of(parent2)) {
        //If the parents are the same, then just find the distance between the nodes and return
        //Find the minimum distance between the two children (net1 and net2)
        if ( has_cached_values1 && seed1.minimizer_cache.parent_is_chain) {
            if (seed1.minimizer_cache.prefix_sum < seed2.minimizer_cache.prefix_sum) {
                //If seed1 comes before seed2
                size_t distance_between = SnarlDistanceIndex::minus( SnarlDistanceIndex::minus(seed2.minimizer_cache.prefix_sum,
                                                                                               seed1.minimizer_cache.prefix_sum),
                                                                    seed1.minimizer_cache.node_length);
                minimum_distance = SnarlDistanceIndex::sum(distance_between, 
                            SnarlDistanceIndex::sum(seed1.minimizer_cache.is_reversed ? distance_to_start1 : distance_to_end1,
                                               seed2.minimizer_cache.is_reversed ? distance_to_end2 : distance_to_start2));
            } else {
                size_t distance_between = SnarlDistanceIndex::minus( SnarlDistanceIndex::minus(seed1.minimizer_cache.prefix_sum,
                                                                                               seed2.minimizer_cache.prefix_sum),
                                                                    seed2.minimizer_cache.node_length);
                minimum_distance = SnarlDistanceIndex::sum(distance_between, 
                            SnarlDistanceIndex::sum(seed2.minimizer_cache.is_reversed ? distance_to_start2 : distance_to_end2,
                                               seed1.minimizer_cache.is_reversed ? distance_to_end1 : distance_to_start1));
            }
        } else { 
            //Otherwise, the parent is a snarl and the distances are found with the index
            size_t distance_start_start = distance_index.distance_in_parent(parent1, distance_index.flip(net1), distance_index.flip(net2), graph);
            size_t distance_start_end = distance_index.distance_in_parent(parent1, distance_index.flip(net1), net2, graph);
            size_t distance_end_start = distance_index.distance_in_parent(parent1, net1, distance_index.flip(net2), graph);
            size_t distance_end_end = distance_index.distance_in_parent(parent1, net1, net2, graph);

            //And add those to the distances we've found to get the minimum distance between the positions
            minimum_distance = std::min(SnarlDistanceIndex::sum({distance_start_start , distance_to_start1 , distance_to_start2}),
                   std::min(SnarlDistanceIndex::sum({distance_start_end , distance_to_start1 , distance_to_end2}),
                   std::min(SnarlDistanceIndex::sum({distance_end_start , distance_to_end1 , distance_to_start2}),
                            SnarlDistanceIndex::sum({distance_end_end , distance_to_end1 , distance_to_end2})))); 
        }
        if (stop_at_lowest_common_ancestor) {
            return minimum_distance == std::numeric_limits<size_t>::max() ? std::numeric_limits<size_t>::max() 
                                                                        : minimum_distance - 1;
        }
    }

    //Otherwise, find the distances to the ends of the parents, update them, and continue
    //only if the parent isn't the common ancestor
    if (parent1 != common_ancestor && !distance_index.is_root(parent1)) {
        if (has_cached_values1 && seed1.minimizer_cache.parent_is_chain && !seed1.minimizer_cache.is_trivial_chain) {
            size_t distance_to_chain_start = seed1.minimizer_cache.prefix_sum;
            size_t distance_to_chain_end = SnarlDistanceIndex::minus(SnarlDistanceIndex::minus(distance_index.minimum_length(parent1),
                                                                seed1.minimizer_cache.prefix_sum), seed1.minimizer_cache.node_length);
            size_t old_distance_to_start = distance_to_start1;
            size_t old_distance_to_end = distance_to_end1;
#ifdef debug_distances
            cerr << "\tUsing cache to update to ends of chain1 using distances " << distance_to_chain_start << " and " << distance_to_chain_end << endl;
#endif

            distance_to_start1 = SnarlDistanceIndex::sum(distance_to_chain_start, 
                                        seed1.minimizer_cache.is_reversed ? old_distance_to_end : old_distance_to_start);
            distance_to_end1 = SnarlDistanceIndex::sum(distance_to_chain_end, 
                                        seed1.minimizer_cache.is_reversed ? old_distance_to_start : old_distance_to_end);
        } else {
            update_distances(net1, parent1, distance_to_start1, distance_to_end1);
        }
        net1 = std::move(parent1);
    }
    if (parent2 != common_ancestor && !distance_index.is_root(parent2)) {
        if (has_cached_values2 && seed2.minimizer_cache.parent_is_chain && !seed2.minimizer_cache.is_trivial_chain) {
            size_t distance_to_chain_start = seed2.minimizer_cache.prefix_sum;
            size_t distance_to_chain_end = SnarlDistanceIndex::minus(SnarlDistanceIndex::minus(distance_index.minimum_length(parent2),
                                                                seed2.minimizer_cache.prefix_sum), seed2.minimizer_cache.node_length);
            size_t old_distance_to_start = distance_to_start2;
            size_t old_distance_to_end = distance_to_end2;
#ifdef debug_distances
            cerr << "\tUsing cache to update to ends of chain2 using distances " << distance_to_chain_start << " and " << distance_to_chain_end << endl;
#endif

            distance_to_start2 = SnarlDistanceIndex::sum(distance_to_chain_start, 
                                        seed2.minimizer_cache.is_reversed ? old_distance_to_end : old_distance_to_start);
            distance_to_end2 = SnarlDistanceIndex::sum(distance_to_chain_end, 
                                        seed2.minimizer_cache.is_reversed ? old_distance_to_start : old_distance_to_end);

        } else {
            update_distances(net2, parent2, distance_to_start2, distance_to_end2);
        }
        net2 = std::move(parent2);
    }



#ifdef debug_distances
        cerr << "Updated to parents" << endl;
        cerr << "Reached node " << distance_index.net_handle_as_string(net1) << " for position 1" << endl;
        cerr << "   with distances to ends " << distance_to_start1 << " and " << distance_to_end1 << endl;
        cerr << "Reached node " << distance_index.net_handle_as_string(net2) << " for position 2" << endl;
        cerr << "   with distances to ends " << distance_to_start2 << " and " << distance_to_end2 << endl;
#endif



    if (!lowest_ancestor.second) {
        //If these are not in the same connected component
#ifdef debug_distances
        cerr << "These are in different connected components" << endl;
#endif
        return std::numeric_limits<size_t>::max();
    }

    /*
     * Walk up the snarl tree until net1 and net2 are children of the lowest common ancestor
     * Keep track of the distances to the ends of the net handles as we go
     */

    if (distance_index.start_end_traversal_of(net1) == distance_index.start_end_traversal_of(net2)){
        if (SnarlDistanceIndex::sum(distance_to_end1 , distance_to_start2) > distance_index.minimum_length(net1) &&
            SnarlDistanceIndex::sum(distance_to_end1 , distance_to_start2) != std::numeric_limits<size_t>::max()) {
            //If the positions are on the same node and are pointing towards each other, then
            //check the distance between them in the node
            minimum_distance = SnarlDistanceIndex::minus(SnarlDistanceIndex::sum(distance_to_end1 , distance_to_start2), 
                                                         distance_index.minimum_length(net1));
        }
        if (SnarlDistanceIndex::sum({distance_to_start1 , distance_to_end2}) > distance_index.minimum_length(net1) &&
            SnarlDistanceIndex::sum({distance_to_start1 , distance_to_end2}) != std::numeric_limits<size_t>::max()) {
            minimum_distance = std::min(SnarlDistanceIndex::minus(SnarlDistanceIndex::sum(distance_to_start1 , distance_to_end2), 
                                                                  distance_index.minimum_length(net1)), 
                                        minimum_distance);
        }
        if (!stop_at_lowest_common_ancestor) {
            common_ancestor = distance_index.start_end_traversal_of(distance_index.get_parent(net1));
        }


    } else {

        //Get the distance from position 1 up to the ends of a child of the common ancestor
#ifdef debug_distances
        cerr << "Reaching the children of the lowest common ancestor for first position..." << endl;
#endif
        while (distance_index.start_end_traversal_of(distance_index.get_parent(net1)) != common_ancestor && !distance_index.is_root(distance_index.get_parent(net1))) {
            net_handle_t parent = distance_index.start_end_traversal_of(distance_index.get_parent(net1));
            update_distances(net1, parent, distance_to_start1, distance_to_end1);
            net1 = parent;
        }
#ifdef debug_distances
        cerr << "Reached node " << distance_index.net_handle_as_string(net1) << " for position 1" << endl;
        cerr << "   with distances to ends " << distance_to_start1 << " and " << distance_to_end1 << endl;
        cerr << "Reaching the children of the lowest common ancestor for position 2..." << endl;
#endif
        //And the same for position 2
        while (distance_index.start_end_traversal_of(distance_index.get_parent(net2)) != distance_index.start_end_traversal_of(common_ancestor) && !distance_index.is_root(distance_index.get_parent(net2))) {
            net_handle_t parent = distance_index.start_end_traversal_of(distance_index.get_parent(net2));
            update_distances(net2, parent, distance_to_start2, distance_to_end2);
            net2 = parent;
        }
#ifdef debug_distances
        cerr << "Reached node " << distance_index.net_handle_as_string(net2) << " for position 2" << endl;
        cerr << "   with distances to ends " << distance_to_start2 << " and " << distance_to_end2 << endl;
#endif
    }
    if (stop_at_lowest_common_ancestor) {
        
        return minimum_distance == std::numeric_limits<size_t>::max() ? std::numeric_limits<size_t>::max() : minimum_distance-1;    
    }

    /* 
     * common_ancestor is now the lowest common ancestor of both net handles, and 
     * net1 and net2 are both children of common_ancestor
     * Walk up to the root and check for distances between the positions within each
     * ancestor
     */

    while (!distance_index.is_root(net1)){
#ifdef debug_distances
            cerr << "At common ancestor " << distance_index.net_handle_as_string(common_ancestor) <<  endl;
            cerr << "  with distances for child 1 (" << distance_index.net_handle_as_string(net1) << "): " << distance_to_start1 << " "  << distance_to_end1 << endl;
            cerr << "                     child 2 (" << distance_index.net_handle_as_string(net2) << "): " << distance_to_start2 << " " <<  distance_to_end2 << endl;
#endif

        //Find the minimum distance between the two children (net1 and net2)
        size_t distance_start_start = distance_index.distance_in_parent(common_ancestor, distance_index.flip(net1), distance_index.flip(net2), graph);
        size_t distance_start_end = distance_index.distance_in_parent(common_ancestor, distance_index.flip(net1), net2, graph);
        size_t distance_end_start = distance_index.distance_in_parent(common_ancestor, net1, distance_index.flip(net2), graph);
        size_t distance_end_end = distance_index.distance_in_parent(common_ancestor, net1, net2, graph);

        //And add those to the distances we've found to get the minimum distance between the positions
        minimum_distance = std::min(minimum_distance,
                           std::min(SnarlDistanceIndex::sum({distance_start_start , distance_to_start1 , distance_to_start2}),
                           std::min(SnarlDistanceIndex::sum({distance_start_end , distance_to_start1 , distance_to_end2}),
                           std::min(SnarlDistanceIndex::sum({distance_end_start , distance_to_end1 , distance_to_start2}),
                                    SnarlDistanceIndex::sum({distance_end_end , distance_to_end1 , distance_to_end2})))));

#ifdef debug_distances
            cerr << "    Found distances between nodes: " << distance_start_start << " " << distance_start_end << " " << distance_end_start << " "      << distance_end_end << endl;
            cerr << "  best distance is " << minimum_distance << endl;
#endif
        if (!distance_index.is_root(common_ancestor)) {
            //Update the distances to reach the ends of the common ancestor
            update_distances(net1, common_ancestor, distance_to_start1, distance_to_end1);
            update_distances(net2, common_ancestor, distance_to_start2, distance_to_end2);

            //Update which net handles we're looking at
            net1 = common_ancestor;
            net2 = common_ancestor;
            common_ancestor = distance_index.start_end_traversal_of(distance_index.get_parent(common_ancestor));
        } else {
            //Just update this one to break out of the loop
            net1 = common_ancestor;
        }
    }

    //minimum distance currently includes both positions
    return minimum_distance == std::numeric_limits<size_t>::max() ? std::numeric_limits<size_t>::max() : minimum_distance-1;
}


}

