#include "snarl_seed_clusterer.hpp"

#include <algorithm>

//#define DEBUG_CLUSTER
namespace vg {

NewSnarlSeedClusterer::NewSnarlSeedClusterer( const SnarlDistanceIndex& distance_index, const HandleGraph* graph) :
                                        distance_index(distance_index),
                                        graph(graph){
};

vector<NewSnarlSeedClusterer::Cluster> NewSnarlSeedClusterer::cluster_seeds (const vector<Seed>& seeds, size_t read_distance_limit) const {
    //Wrapper for single ended

    vector<const vector<Seed>*> all_seeds;
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

vector<vector<NewSnarlSeedClusterer::Cluster>> NewSnarlSeedClusterer::cluster_seeds (
              const vector<vector<Seed>>& all_seeds, size_t read_distance_limit,
              size_t fragment_distance_limit) const {

    //Wrapper for paired end
    vector<const vector<Seed>*> seed_pointers;
    seed_pointers.reserve(all_seeds.size());
    for (const vector<Seed>& v : all_seeds) seed_pointers.push_back(&v);

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


tuple<vector<structures::UnionFind>, structures::UnionFind> NewSnarlSeedClusterer::cluster_seeds_internal (
              const vector<const vector<Seed>*>& all_seeds, size_t read_distance_limit,
              size_t fragment_distance_limit) const {
    /* Given a vector of seeds and a limit, find a clustering of seeds where
     * seeds that are closer than the limit cluster together.
     * Returns a vector of cluster assignments
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
    //Also maps to the chain's NodeClusters in all_node_clusters
    //Children are represented as a size_t index into all_node_clusters
    //Starts out with nodes from get_nodes, 
    //then adds snarls as we walk up the snarl tree in 
    //This is later used to populate snarl_to_node in the tree state
    vector<hash_map<net_handle_t, pair<size_t, vector<size_t>>>> chain_to_children_by_level;

    chain_to_children_by_level.reserve(distance_index.get_max_tree_depth()+1);



    //This stores all the tree relationships and cluster information
    //for a single level of the snarl tree as it is being processed
    //It also keeps track of the parents of the current level
    size_t seed_count = 0;
    for (auto v : all_seeds) seed_count+= v->size();
    TreeState tree_state (&all_seeds, read_distance_limit, fragment_distance_limit, seed_count);


    //Populate tree_state.node_to_seeds (mapping each node to the seeds it
    //contains) and chain_to_children_by_level
    get_nodes(tree_state, chain_to_children_by_level);

    //Initialize the tree state to the bottom level
    tree_state.chain_to_children = &chain_to_children_by_level[chain_to_children_by_level.size() - 1];

    for (int depth = chain_to_children_by_level.size() - 1 ; depth >= 0 ; depth --) {
        //Go through each level of the tree, bottom up, and cluster that
        // level. Each level includes the snarl at that level, the nodes
        // belonging to those snarls, and the chains comprised of them
        //
        // tree_state knows all children of the snarls at this level

#ifdef DEBUG_CLUSTER
assert(tree_state.read_index_offsets[0] == 0);
for (size_t i = 1 ; i < tree_state.all_seeds->size() ; i++) {
    assert (tree_state.read_index_offsets[i] + tree_state.all_seeds->at(i)->size() == tree_state.read_index_offsets[i+1]);
}
#endif
        if (depth != 0) {
            tree_state.parent_chain_to_children = &chain_to_children_by_level[depth-1];
        }

        //Cluster all the chains at this depth
        //Also records which chains are in snarls and the parents of these
        //chains in tree_state.parent_chain_to_children
        cluster_chain_level(tree_state);

        //And cluster all the snarls, record the parents of these snarls
        //cluster_snarl_level(tree_state);


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


//Find which nodes contain seeds and assign those nodes to the
//snarls that contain them
//Update the tree state's node_to_seed
//and snarl_to_nodes_by_level, which assigns each node that contains
//seeds to a snarl, organized by the level of the snarl in the snarl
//tree. chain_to_children_by_level will be used to populate snarl_to_nodes
//in the tree state as each level is processed
void NewSnarlSeedClusterer::get_nodes( TreeState& tree_state,
          vector<hash_map<net_handle_t,pair<size_t, vector<size_t>>>>& chain_to_children_by_level) const {
#ifdef DEBUG_CLUSTER
cerr << "Add all seeds to nodes: " << endl << "\t";
#endif


    // Assign each seed to a node.
    hash_set<id_t> seen_nodes;
    for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++){ 
        const vector<Seed>* seeds = tree_state.all_seeds->at(read_num);
        vector<pair<id_t, size_t>>& node_to_seeds = tree_state.node_to_seeds.at(read_num);
        for (size_t i = 0; i < seeds->size(); i++) {
            const Seed& seed = seeds->at(i);
            pos_t pos = seed.pos;
            id_t id = get_id(pos);
            

            node_to_seeds.emplace_back(id, i);
#ifdef DEBUG_CLUSTER
            cerr << read_num << ":" << pos << ", ";
#endif

             //And the node to a chain
            if (seen_nodes.count(id) < 1) {
                 seen_nodes.insert(id);
                 net_handle_t node_net_handle = distance_index.get_node_net_handle(id) ; 
                 net_handle_t parent = std::get<1>(seed.minimizer_cache) == MIPayload::NO_VALUE 
                                      ? distance_index.get_parent(node_net_handle)
                                      : distance_index.get_handle_from_connected_component(std::get<1>(seed.minimizer_cache));

                 //Get the values from the seed. Some may be infinite and need to be re-set
                 size_t depth = std::get<1>(seed.minimizer_cache) == MIPayload::NO_VALUE ? distance_index.get_depth(parent) : 1;
                 size_t node_length = std::get<0>(seed.minimizer_cache);
                 size_t prefix_sum = std::get<2>(seed.minimizer_cache);
                 size_t component = std::get<3>(seed.minimizer_cache);
                 bool is_reversed_in_parent = std::get<4>(seed.minimizer_cache);

                 //Seed payload is: node length, root component, prefix sum, chain component, is_reversed 
                 //Node length and is_reversed are always set
                 if (prefix_sum == MIPayload::NO_VALUE) {
                     //If we didn't store information in the seed, then get it from the distance index

                    prefix_sum = distance_index.is_trivial_chain(parent) 
                            ? std::numeric_limits<size_t>::max() 
                            : distance_index.get_prefix_sum_value(node_net_handle);
                    component = distance_index.is_multicomponent_chain(parent) 
                            ? distance_index.get_chain_component(node_net_handle)
                            : 0;
                    if (node_length == MIPayload::NO_VALUE) {
                        node_length = distance_index.minimum_length(node_net_handle);
                        is_reversed_in_parent = distance_index.is_reversed_in_parent(node_net_handle);
                    }
                 }
                 if (depth+1 > chain_to_children_by_level.size()) {
                     chain_to_children_by_level.resize(depth+1);
                     chain_to_children_by_level.back().reserve(seeds->size()*2);
                 }



                 tree_state.all_node_clusters.emplace_back(
                                 NodeClusters(std::move(node_net_handle), tree_state.all_seeds->size(),
                                              is_reversed_in_parent, id, node_length, prefix_sum, component));
                 if (distance_index.is_root(parent)) {
                     //If this is a child of the root, then cluster it now
                     cluster_one_node(tree_state, tree_state.all_node_clusters.back());
                     if (distance_index.is_root_snarl(parent)) {
                         //If this is a root snarl, then remember it to cluster in the root
                        if (tree_state.root_children.count(parent) == 0) {
                           tree_state.root_children.emplace(parent, vector<size_t>(0));
                        }
                        tree_state.root_children[parent].emplace_back(tree_state.all_node_clusters.size()-1);
                     } else {
                         //Otherwise, just compare the single child's external connectivity
                         compare_and_combine_cluster_on_one_child(tree_state, tree_state.all_node_clusters.back());

                     }
                 } else {
                     //Otherwise, add this node to its parent
                    add_child_to_vector(tree_state, &chain_to_children_by_level[depth], parent, tree_state.all_node_clusters.size() - 1);
                 }
            }
        }
        std::sort(node_to_seeds.begin(), node_to_seeds.end());
    }
#ifdef DEBUG_CLUSTER
    cerr << endl;
#endif

    if (chain_to_children_by_level.empty()) {
        chain_to_children_by_level.resize(1);
    }
}



//Cluster all of the snarls in tree_state from the same depth (all present in snarl_to_children)
//Assumes that all the children of the snarls have been clustered already and are present in tree_state.snarls_to_children
//void NewSnarlSeedClusterer::cluster_snarl_level(TreeState& tree_state) const {
//
//    for (auto& kv : tree_state.snarl_to_children){
//        //Go through each of the snarls at this level, cluster them,
//        //and find which chains they belong to, if any
//        //key is the index of the snarl and value is a vector of pair of
//        // net_handle_t, NodeClusters
//
//        net_handle_t snarl_handle = kv.first;
//
//#ifdef DEBUG_CLUSTER
//        cerr << "Cluster one snarl " << distance_index.net_handle_as_string(snarl_handle) << endl;
//#endif
//        net_handle_t snarl_parent = distance_index.get_parent(snarl_handle);
//        //If this snarl is in a chain, cluster and add let the
//        //tree state know which chain it belongs to (in order)
//
//        //Cluster the snarl
//        cluster_one_snarl(tree_state, tree_state.all_node_clusters[kv.second.first]);
//        NodeClusters& snarl_clusters = tree_state.all_node_clusters[kv.second.first];
//
//        //Now check the best distance of any seed to the ends of the snarl
//        //Is the distance small enough that we can cluster it with something else?
//        bool reachable_right = snarl_clusters.fragment_best_right <= 
//            (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit);
//        bool reachable_left = snarl_clusters.fragment_best_left <= 
//            (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit);
//
//
//        if (reachable_left || reachable_right) {
//            //If we can reach the ends of the snarl, add it to it's parent 
//            if (distance_index.is_root(snarl_parent)) {
//                 if(distance_index.is_root_snarl(snarl_parent)){
//                     //If the parent is a root snarl, then remember it to be compared in the root
//                    if (tree_state.root_children.count(snarl_parent) == 0) {
//                       tree_state.root_children.emplace(snarl_parent, vector<size_t>(0));
//                    }
//                    tree_state.root_children[snarl_parent].emplace_back(kv.second.first);
//                 } else {
//                     //Otherwise, compare it to itself using external connectivity
//                         compare_and_combine_cluster_on_one_child(tree_state, snarl_clusters);
//                 }
//            } else {
//                add_child_to_vector(tree_state, tree_state.parent_chain_to_children, snarl_parent, kv.second.first);
//            }
//        }
//
//#ifdef DEBUG_CLUSTER
//        cerr << "\tRecording snarl " << distance_index.net_handle_as_string(snarl_handle)  << " as a child of "
//              << distance_index.net_handle_as_string(distance_index.get_parent(snarl_handle)) << endl;
//#endif
//
//    }
//    tree_state.snarl_to_children.clear();
//}


void NewSnarlSeedClusterer::cluster_chain_level(TreeState& tree_state) const {
    for (auto& kv : (*tree_state.chain_to_children)) {
        //For each chain at this level that has relevant child snarls in it,
        //find the clusters.


        // Get the chain's number
        net_handle_t chain_handle = kv.first;

#ifdef DEBUG_CLUSTER
        cerr << "Cluster one chain " <<  distance_index.net_handle_as_string(chain_handle) << endl;
#endif

        // Compute the clusters for the chain
        cluster_one_chain(tree_state, tree_state.all_node_clusters[kv.second.first]);

        //Add the chain to its parent
        net_handle_t parent = distance_index.get_parent(chain_handle);
        if (distance_index.is_root(parent)) {
            //If the parent is the root, remember the index of this chain in all_node_clusters
            if(distance_index.is_root_snarl(parent)) {
                //If the parent is a root snarl, then remember it to cluster in the root
                if (tree_state.root_children.count(parent) == 0) {
                   tree_state.root_children.emplace(parent, vector<size_t>(0));
                }
                tree_state.root_children[parent].emplace_back(kv.second.first);
            } else {
                //Otherwise, cluster it with itself using external connectivity only
                 compare_and_combine_cluster_on_one_child(tree_state, tree_state.all_node_clusters[kv.second.first]);
            }
        } else {
            //If the parent is just a snarl, add it to its parent snarl
            size_t parent_snarl_index = add_child_to_vector(tree_state, &tree_state.snarl_to_children, parent, kv.second.first);

            //And add the snarl to its parent chain
            net_handle_t grandparent_chain = distance_index.get_parent(parent);
            if (distance_index.is_root(grandparent_chain)) {
                assert(false);
                 //if(distance_index.is_root_snarl(grandparent_chain)){
                 //    //If the parent is a root snarl, then remember it to be compared in the root
                 //   if (tree_state.root_children.count(grandparent_chain) == 0) {
                 //      tree_state.root_children.emplace(grandparent_chain, vector<size_t>(0));
                 //   }
                 //   tree_state.root_children[grandparent_chain].emplace_back(parent_snarl_index);
                 //} else {
                 //    //Otherwise, compare it to itself using external connectivity
                 //    cluster_one_snarl(tree_state,
                 //    compare_and_combine_cluster_on_one_child(tree_state, snarl_clusters);
                 //}
            } else {
                if (parent_snarl_index != std::numeric_limits<size_t>::max()) {
                    //If this is the first time we see this parent snarl, then add it to the grandparent chain
                    add_child_to_vector(tree_state, tree_state.parent_chain_to_children, grandparent_chain, parent_snarl_index);
                }
            }
        }


    }
}


void NewSnarlSeedClusterer::cluster_one_node(
                   TreeState& tree_state, NodeClusters& node_clusters) const {
#ifdef DEBUG_CLUSTER
    cerr << "Finding clusters on node " << distance_index.net_handle_as_string(node_clusters.containing_net_handle) << endl;
#endif

    size_t node_length = node_clusters.node_length;//TODO:distance_index.minimum_length(node_clusters.containing_net_handle);
    nid_t node_id = node_clusters.node_id;


    if (tree_state.read_distance_limit >= node_length) {
        //If the limit is greater than the node length, then all the
        //seeds on this node must be in the same cluster

        size_t fragment_group_id = std::numeric_limits<size_t>::max();
        for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
            auto seed_range_start = std::lower_bound(
                tree_state.node_to_seeds[read_num].begin(), tree_state.node_to_seeds[read_num].end(),
                std::pair<id_t, size_t>(node_id, 0));
            if (seed_range_start != tree_state.node_to_seeds[read_num].end() 
                    && seed_range_start->first == node_id) {

                size_t group_id = seed_range_start->second;

                for (auto iter = seed_range_start; iter != tree_state.node_to_seeds.at(read_num).end() 
                                                  && iter->first == node_id; ++iter) {
                    //For each seed on this node, add it to the cluster
                    //And find the shortest distance from any seed to both
                    //ends of the node


                    //Get the seed we're looking at
                    size_t& seed_i = iter->second;
                    pos_t seed = tree_state.all_seeds->at(read_num)->at(seed_i).pos;
                    //And the distances for this seed
                    size_t dist_left = is_rev(seed) ? node_length- get_offset(seed) : get_offset(seed) + 1;
                    size_t dist_right = is_rev(seed) ? get_offset(seed) + 1 : node_length - get_offset(seed);

                    //Find the new best distances for anything in this cluster
                    node_clusters.read_best_left[read_num] = std::min(dist_left,
                                                          node_clusters.read_best_left[read_num]);
                    node_clusters.read_best_right[read_num] = std::min(dist_right,
                                                          node_clusters.read_best_right[read_num]);
                    node_clusters.fragment_best_left = std::min(dist_left,
                                                          node_clusters.fragment_best_left);
                    node_clusters.fragment_best_right = std::min(dist_right,
                                                          node_clusters.fragment_best_right);

                    //Put this seed in the cluster for the node
                    tree_state.read_union_find[read_num].union_groups(group_id, seed_i);
                    if (tree_state.fragment_distance_limit != 0 ) {
                        if (fragment_group_id == std::numeric_limits<size_t>::max() ) {
                            fragment_group_id = seed_range_start->second + tree_state.read_index_offsets[read_num];
                        }
                        tree_state.fragment_union_find.union_groups(
                                fragment_group_id, seed_i + tree_state.read_index_offsets[read_num]);
                    }

                }

                //Record the new cluster
                group_id = tree_state.read_union_find[read_num].find_group(group_id);
                node_clusters.read_cluster_heads.emplace(read_num, group_id);
                tree_state.read_cluster_heads_to_distances[read_num][group_id] =  
                    make_pair(node_clusters.read_best_left[read_num], node_clusters.read_best_right[read_num]);

                if (tree_state.fragment_distance_limit != 0) {
                    fragment_group_id = tree_state.fragment_union_find.find_group(fragment_group_id);
                }
            }
        }
#ifdef DEBUG_CLUSTER
        cerr << "\tFound single cluster on node " << node_id << " with fragment dists " 
                << node_clusters.fragment_best_left << " " << node_clusters.fragment_best_right << endl;

        bool got_left = false;
        bool got_right = false;
        for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
            cerr << "\t for read num " << read_num << " best left: " << node_clusters.read_best_left[read_num] << " best right: " << node_clusters.read_best_right[read_num] << endl;
            bool got_read_left=false;
            bool got_read_right = false;
            for (pair<size_t,size_t> c : node_clusters.read_cluster_heads) {
                if (c.first == read_num) {
                    pair<size_t, size_t> dists = tree_state.read_cluster_heads_to_distances[c.first][c.second];
                    cerr << "\t\t" << c.first << ":"<<c.second << ": left: " << dists.first << " right : " << dists.second << ": ";
                    bool has_seeds = false;
                    for (size_t x = 0 ; x < tree_state.all_seeds->at(c.first)->size() ; x++) {
                        if (tree_state.read_union_find[c.first].find_group(x) == c.second) {
                            cerr << tree_state.all_seeds->at(c.first)->at(x).pos << " ";
                            has_seeds = true;
                        }
                    }
                    //assert(dists.first == std::numeric_limits<size_t>::max() || dists.first >= node_clusters.read_best_left[read_num]);
                    //assert(dists.second == std::numeric_limits<size_t>::max() || dists.second >= node_clusters.read_best_right[read_num]);
                    //assert(dists.first == std::numeric_limits<size_t>::max() || dists.first >= node_clusters.fragment_best_left);
                    //assert(dists.second == std::numeric_limits<size_t>::max() || dists.second >= node_clusters.fragment_best_right);
                    if (dists.first == node_clusters.fragment_best_left) {got_left = true;}
                    if (dists.second == node_clusters.fragment_best_right) {got_right = true;}
                    if (dists.first == node_clusters.read_best_left[read_num]) {got_read_left = true;}
                    if (dists.second == node_clusters.read_best_right[read_num]) {got_read_right = true;}
                    cerr << endl;
                    //assert(has_seeds);
                }
            }
            //assert(got_read_left || node_clusters.read_best_left[read_num] == std::numeric_limits<size_t>::max());
            //assert(got_read_right || node_clusters.read_best_right[read_num] == std::numeric_limits<size_t>::max());
        }
        //assert(got_left);
        //assert(got_right);

#endif
        return;
    }


    //The seeds may for multiple clusters on the node
    //Sort the seeds by their offset in the node and split into clusters
        
    //<index of read, index of seed, offset of seed> for all seeds
    vector<tuple<size_t,size_t, size_t>> seed_offsets;
    for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
        auto seed_range_start = std::lower_bound(
            tree_state.node_to_seeds[read_num].begin(),tree_state.node_to_seeds[read_num].end(),
            std::pair<id_t, size_t>(node_id, 0));
        if (seed_range_start != tree_state.node_to_seeds[read_num].end() && seed_range_start->first == node_id) {
            for (auto iter = seed_range_start; iter != tree_state.node_to_seeds[read_num].end() 
                                                && iter->first == node_id; ++iter) {
                //For each seed, find its offset
                pos_t seed = tree_state.all_seeds->at(read_num)->at(iter->second).pos;
                size_t offset = is_rev(seed) ? node_length - get_offset(seed) : get_offset(seed) + 1;

                node_clusters.fragment_best_left = std::min(offset, node_clusters.fragment_best_left);
                node_clusters.fragment_best_right = std::min(node_length-offset+1, node_clusters.fragment_best_right);
                node_clusters.read_best_left[read_num] = std::min(offset, node_clusters.read_best_left[read_num]);
                node_clusters.read_best_right[read_num] = std::min(node_length-offset+1, node_clusters.read_best_right[read_num]);

                seed_offsets.emplace_back(read_num, iter->second, offset);

            }
        }
    }
    //Sort seeds by their position in the node
    std::sort(seed_offsets.begin(), seed_offsets.end(),
                 [&](const auto a, const auto b) -> bool {
                      return  std::get<2>(a) < std::get<2>(b);
                  } );

    //Offset of the first seed for the cluster we're building
    //One for each read
    vector<size_t> read_first_offset (tree_state.all_seeds->size(), std::numeric_limits<size_t>::max());
    //Offset of the latest seed for the cluster we're building
    vector<size_t> read_last_offset (tree_state.all_seeds->size(), std::numeric_limits<size_t>::max());
    //And the same for the fragment 
    size_t fragment_last_offset = std::numeric_limits<size_t>::max();
    size_t fragment_last_cluster = std::numeric_limits<size_t>::max();
    vector<size_t> read_last_cluster (tree_state.all_seeds->size(), std::numeric_limits<size_t>::max());

    //Each s is <index of read, index of seed, offset of seed>
    for ( tuple<size_t, size_t, size_t> s : seed_offsets) {
        //For each seed, in order of position in the node,
        //see if it belongs to a new read/fragment cluster - if it is
        //close enough to the previous seed
        size_t read_num = std::get<0>(s);
        size_t seed_num = std::get<1>(s);
        size_t offset = std::get<2>(s);
        if (read_first_offset[read_num] == std::numeric_limits<size_t>::max()) {
            read_first_offset[read_num] = offset;
        }

        if (read_last_offset[read_num] != std::numeric_limits<size_t>::max() &&
            offset - read_last_offset[read_num] <= tree_state.read_distance_limit) {
            //If this seed is in the same read cluster as the previous one,
            //union them

            tree_state.read_union_find[read_num].union_groups(seed_num, read_last_cluster[read_num]);
            read_last_cluster[read_num] = tree_state.read_union_find[read_num].find_group(seed_num);
            read_last_offset[read_num] = offset;

            if (tree_state.fragment_distance_limit != 0) {
                //If we are also clustering paired end reads by fragment distance,
                //cluster these together
                tree_state.fragment_union_find.union_groups(seed_num+tree_state.read_index_offsets[read_num], fragment_last_cluster);
                fragment_last_cluster = tree_state.fragment_union_find.find_group(seed_num+tree_state.read_index_offsets[read_num]);
            }
        } else {
            //This becomes a new read cluster
            if (read_last_cluster[read_num] != std::numeric_limits<size_t>::max()) {
                //Record the previous cluster
                node_clusters.read_cluster_heads.emplace(read_num, read_last_cluster[read_num]);
                tree_state.read_cluster_heads_to_distances[read_num][read_last_cluster[read_num]] = 
                    make_pair(read_first_offset[read_num], node_length - read_last_offset[read_num] + 1);
            }
            read_last_cluster[read_num] = seed_num;
            read_first_offset[read_num] = offset;
            read_last_offset[read_num] = offset;
            if (tree_state.fragment_distance_limit != 0) {
                if (fragment_last_offset != std::numeric_limits<size_t>::max() &&
                    offset - fragment_last_offset <= tree_state.fragment_distance_limit) {
                    //If this is a new read cluster but the same fragment cluster
                    tree_state.fragment_union_find.union_groups(seed_num+tree_state.read_index_offsets[read_num], fragment_last_cluster);
                    fragment_last_cluster = tree_state.fragment_union_find.find_group(fragment_last_cluster);

                } else {
                    //If this is a new fragment cluster as well
                    fragment_last_cluster = seed_num+tree_state.read_index_offsets[read_num];
                }
            }
        }
        fragment_last_offset = offset;
    }
    for (size_t i = 0 ; i < read_last_cluster.size() ; i++) {
        if (read_last_cluster[i] != std::numeric_limits<size_t>::max()) {
            node_clusters.read_cluster_heads.emplace(i, read_last_cluster[i]);
            tree_state.read_cluster_heads_to_distances[i][read_last_cluster[i]] =
                make_pair(read_first_offset[i], node_length-read_last_offset[i]+1);
        }
    }

#ifdef DEBUG_CLUSTER

    cerr << "\tFound read clusters on node " << node_id << endl;

    bool got_left = false;
    bool got_right = false;
    for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
        cerr << "\t for read num " << read_num << " best left: " << node_clusters.read_best_left[read_num] << " best right: " << node_clusters.read_best_right[read_num] << endl;
        bool got_read_left=false;
        bool got_read_right = false;
        for (pair<size_t,size_t> c : node_clusters.read_cluster_heads) {
            if (c.first == read_num) {
                pair<size_t, size_t> dists = tree_state.read_cluster_heads_to_distances[c.first][c.second];
                cerr << "\t\t" << c.first << ":"<<c.second << ": left: " << dists.first << " right : " << dists.second << ": ";
                bool has_seeds = false;
                for (size_t x = 0 ; x < tree_state.all_seeds->at(c.first)->size() ; x++) {
                    if (tree_state.read_union_find[c.first].find_group(x) == c.second) {
                        cerr << tree_state.all_seeds->at(c.first)->at(x).pos << " ";
                        has_seeds = true;
                    }
                }
                assert(dists.first == std::numeric_limits<size_t>::max() || dists.first >= node_clusters.read_best_left[read_num]);
                assert(dists.second == std::numeric_limits<size_t>::max() || dists.second >= node_clusters.read_best_right[read_num]);
                assert(dists.first == std::numeric_limits<size_t>::max() || dists.first >= node_clusters.fragment_best_left);
                assert(dists.second == std::numeric_limits<size_t>::max() || dists.second >= node_clusters.fragment_best_right);
                if (dists.first == node_clusters.fragment_best_left) {got_left = true;}
                if (dists.second == node_clusters.fragment_best_right) {got_right = true;}
                if (dists.first == node_clusters.read_best_left[read_num]) {got_read_left = true;}
                if (dists.second == node_clusters.read_best_right[read_num]) {got_read_right = true;}
                cerr << endl;
                assert(has_seeds);
            }
        }
        //assert(got_read_left || node_clusters.read_best_left[read_num] == std::numeric_limits<size_t>::max());
        //assert(got_read_right || node_clusters.read_best_right[read_num] == std::numeric_limits<size_t>::max());
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
//Also updates the distances to the ends of the parent for the second child's clusters
//TODO: Make sure to add the first child's clusters to the parent before looking at pairs and calling this
void NewSnarlSeedClusterer::compare_and_combine_cluster_on_child_structures(TreeState& tree_state, NodeClusters& child_clusters1, 
    NodeClusters& child_clusters2, NodeClusters& parent_clusters, 
    vector<vector<pair<size_t, size_t>>> & child_distances, bool is_root) const {
#ifdef DEBUG_CLUSTER
    cerr << "\tCompare " << distance_index.net_handle_as_string(child_clusters1.containing_net_handle) 
         << " and " << distance_index.net_handle_as_string(child_clusters2.containing_net_handle)
         << " which are children of " << distance_index.net_handle_as_string(parent_clusters.containing_net_handle) << endl;
#endif

    net_handle_t& parent_handle = parent_clusters.containing_net_handle;
    net_handle_t& child_handle1 = child_clusters1.containing_net_handle;
    net_handle_t& child_handle2 = child_clusters2.containing_net_handle;


    //Get the distances between the two sides of the children in the parent
    size_t distance_left_left = distance_index.distance_in_parent(parent_handle, distance_index.flip(child_handle1), distance_index.flip(child_handle2), graph,
            (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit));
    size_t distance_left_right = distance_index.distance_in_parent(parent_handle, distance_index.flip(child_handle1), child_handle2, graph,
            (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit));
    size_t distance_right_right = distance_index.distance_in_parent(parent_handle, child_handle1, child_handle2, graph,
            (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit));
    size_t distance_right_left = distance_index.distance_in_parent(parent_handle, child_handle1, distance_index.flip(child_handle2), graph,
            (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit));

    /* Find the distances from the start/end of the parent to the left/right of the first child
     * If this is the root, then the distances are infinite since it has no start/end
     */

    if (is_root){
        if (distance_left_left == std::numeric_limits<size_t>::max() &&
            distance_left_right == std::numeric_limits<size_t>::max() &&
            distance_right_right == std::numeric_limits<size_t>::max() &&
            distance_right_left == std::numeric_limits<size_t>::max()) {
            return;
        }

        child_clusters1.distance_start_left = std::numeric_limits<size_t>::max();
        child_clusters1.distance_start_right = std::numeric_limits<size_t>::max();
        child_clusters1.distance_end_left = std::numeric_limits<size_t>::max();
        child_clusters1.distance_end_right = std::numeric_limits<size_t>::max();

    } else {
        child_clusters1.distance_start_left = 
            distance_index.distance_in_parent(parent_handle, parent_clusters.start_in, distance_index.flip(child_handle1));

        child_clusters1.distance_start_right = 
            distance_index.distance_in_parent(parent_handle, parent_clusters.start_in, child_handle1);

        child_clusters1.distance_end_left =
            distance_index.distance_in_parent(parent_handle, parent_clusters.end_in, distance_index.flip(child_handle1));

        child_clusters1.distance_end_right = 
            distance_index.distance_in_parent(parent_handle, parent_clusters.end_in, child_handle1);

    }
   
#ifdef DEBUG_CLUSTER
    cerr << "\t\tFound distances between the two children: " << distance_left_left << " " << distance_left_right << " " << distance_right_right << " " << distance_right_left << endl;
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

    //These will be the cluster heads and distances of everything combined by taking the indicated path
    //one cluster head per read
    //pair< pair<read_num, seed_num>, pair<left dist, right_dist>>
    //The default value will be ((inf, 0), (0,0)). Only the inf gets checked to see if it's a real value so I filled it in with 0 so I wouldn't have to type out inf 
    vector<pair<pair<size_t, size_t>, pair<size_t, size_t>>> new_cluster_left_left_by_read(tree_state.all_seeds->size(), 
            make_pair(make_pair(std::numeric_limits<size_t>::max(), 0), make_pair(0, 0)));
    vector<pair<pair<size_t, size_t>, pair<size_t, size_t>>> new_cluster_left_right_by_read(tree_state.all_seeds->size(),
        make_pair(make_pair(std::numeric_limits<size_t>::max(), 0), make_pair(0, 0)));
    vector<pair<pair<size_t, size_t>, pair<size_t, size_t>>> new_cluster_right_right_by_read(tree_state.all_seeds->size(),
        make_pair(make_pair(std::numeric_limits<size_t>::max(), 0), make_pair(0, 0)));
    vector<pair<pair<size_t, size_t>, pair<size_t, size_t>>> new_cluster_right_left_by_read(tree_state.all_seeds->size(),
        make_pair(make_pair(std::numeric_limits<size_t>::max(), 0), make_pair(0, 0)));

    //And the new cluster heads for the fragment
    //These are the values of the cluster heads in the union finds, which include the values from read_index_offset
    size_t new_cluster_left_left_fragment = std::numeric_limits<size_t>::max();
    size_t new_cluster_left_right_fragment = std::numeric_limits<size_t>::max();
    size_t new_cluster_right_right_fragment = std::numeric_limits<size_t>::max();
    size_t new_cluster_right_left_fragment = std::numeric_limits<size_t>::max();

    //Helper function that will compare two clusters
    //Given the read num and seed_num of the cluster head, the distance to the other node side we're looking at, 
    //the distances to the ends of the parent for the cluster head, a reference
    //to the current cluster head and distances of the potential combined cluster (pair<pair<>pair<>> which will be updated if it gets combined),
    //the relevant combined cluster head for the fragment
    //Returns true if this cluster got combined
    auto compare_and_combine_clusters = [&] (size_t read_num, size_t cluster_num, size_t distance_between_reads, 
            size_t distance_between_fragments, pair<size_t, size_t> old_distances, 
            pair<pair<size_t, size_t>, pair<size_t, size_t>>& new_cluster_head_and_distances, size_t& new_cluster_head_fragment){
        if (read_num == new_cluster_head_and_distances.first.first && cluster_num ==  new_cluster_head_and_distances.first.second) {
            //If this is the same as the old cluster head, then don't bother trying to compare
            return false;
        }
        distance_between_reads = SnarlDistanceIndex::minus(distance_between_reads, 1);
        distance_between_fragments = SnarlDistanceIndex::minus(distance_between_fragments, 1);
        bool combined = false;

        if (distance_between_reads <= tree_state.read_distance_limit) {
            //If this can be combined with the given combined cluster
            if (new_cluster_head_and_distances.first.first == std::numeric_limits<size_t>::max()){
                //new cluster head
                new_cluster_head_and_distances = make_pair(make_pair(read_num, cluster_num),
                                                           old_distances);
            } else {
                //Combine with old cluster head
                tree_state.read_union_find.at(read_num).union_groups(cluster_num, new_cluster_head_and_distances.first.second);

                //Update distances
                size_t new_best_left = std::min(old_distances.first, new_cluster_head_and_distances.second.first);
                size_t new_best_right = std::min(old_distances.second, new_cluster_head_and_distances.second.second);

                //And remember new head and distances
                new_cluster_head_and_distances = make_pair(make_pair(read_num, tree_state.read_union_find.at(read_num).find_group(cluster_num)),
                                                           make_pair(new_best_left, new_best_right));
            }
            //Remember to erase this cluster head
            to_erase.emplace_back(read_num, cluster_num);
            combined = true;

#ifdef DEBUG_CLUSTER
            cerr << "\t\t\tCombining read/cluster " << read_num << "/" << cluster_num << "... new cluster head:" << new_cluster_head_and_distances.first.second << endl; 
            cerr << "\t\t\t\t Best distances for this cluster: " << old_distances.first << " and " << old_distances.second << endl;
            cerr << "\t\t\t\t New best distances for combined cluster: " << new_cluster_head_and_distances.second.first << " and " << new_cluster_head_and_distances.second.second << endl;
#endif
        }
        if (tree_state.fragment_distance_limit != 0 && 
                    distance_between_fragments <= tree_state.fragment_distance_limit ) {
            //Just union the fragment
            if (new_cluster_head_fragment == std::numeric_limits<size_t>::max()) {
                new_cluster_head_fragment =cluster_num+tree_state.read_index_offsets[read_num];
            } else {
                tree_state.fragment_union_find.union_groups(cluster_num+tree_state.read_index_offsets[read_num], 
                                                                    new_cluster_head_fragment);
                new_cluster_head_fragment = tree_state.fragment_union_find.find_group(new_cluster_head_fragment);
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
    for (auto& child_cluster_head : child_clusters1.read_cluster_heads) {

        bool combined = false;
        size_t read_num = child_cluster_head.first;
        size_t cluster_num = tree_state.read_union_find[read_num].find_group(child_cluster_head.second);

        //Distances to the ends of the child
        pair<size_t, size_t> distances = child_distances[read_num][child_cluster_head.second];

        //Distances to the parent
        size_t new_dist_left = std::min(SnarlDistanceIndex::sum({distances.first, child_clusters1.distance_start_left}),
                                        SnarlDistanceIndex::sum({distances.second, child_clusters1.distance_start_right}));
        size_t new_dist_right= std::min(SnarlDistanceIndex::sum({distances.first, child_clusters1.distance_end_left}),
                                        SnarlDistanceIndex::sum({distances.second, child_clusters1.distance_end_right}));
        pair<size_t, size_t> distances_to_parent = make_pair(new_dist_left, new_dist_right);
        //If this is already in the parent, take the minimum of the parent distances
        if (parent_clusters.read_cluster_heads.count(make_pair(read_num, cluster_num)) > 0) {
            distances_to_parent = make_pair(
                std::min(new_dist_left, tree_state.read_cluster_heads_to_distances[read_num][cluster_num].first),
                std::min(new_dist_right, tree_state.read_cluster_heads_to_distances[read_num][cluster_num].second));
        }
            

        //Check if the left of 1 can connect with the left of 2
        combined = combined | compare_and_combine_clusters (read_num, cluster_num, 
            SnarlDistanceIndex::sum({distances.first,distance_left_left, child_clusters2.read_best_left[read_num]}), 
            SnarlDistanceIndex::sum({distances.first,distance_left_left, child_clusters2.fragment_best_left}), 
            distances_to_parent,
            new_cluster_left_left_by_read[read_num], new_cluster_left_left_fragment);

        //Check if the left of 1 can connect with the right of 2
        combined = combined | compare_and_combine_clusters (read_num, cluster_num, 
            SnarlDistanceIndex::sum({distances.first,distance_left_right, child_clusters2.read_best_right[read_num]}), 
            SnarlDistanceIndex::sum({distances.first,distance_left_right, child_clusters2.fragment_best_right}), 
            distances_to_parent, new_cluster_left_right_by_read[read_num], new_cluster_left_right_fragment);

        //Check if the right of 1 can connect with the right of 2
        combined = combined | compare_and_combine_clusters (read_num, cluster_num, 
            SnarlDistanceIndex::sum({distances.second,distance_right_right, child_clusters2.read_best_right[read_num]}), 
            SnarlDistanceIndex::sum({distances.second,distance_right_right, child_clusters2.fragment_best_right}), 
            distances_to_parent,new_cluster_right_right_by_read[read_num], new_cluster_right_right_fragment);

        //Check if the right of 1 can connect with the left of 2
        combined = combined | compare_and_combine_clusters (read_num, cluster_num, 
            SnarlDistanceIndex::sum({distances.second,distance_right_left, child_clusters2.read_best_left[read_num]}), 
            SnarlDistanceIndex::sum({distances.second,distance_right_left, child_clusters2.fragment_best_left}), 
            distances_to_parent, new_cluster_right_left_by_read[read_num], new_cluster_right_left_fragment);

        size_t cluster_head = tree_state.read_union_find[read_num].find_group(cluster_num); 

        //Is the distance small enough that we can cluster it with something else?
        bool reachable_left = distances_to_parent.first <= 
            (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit);
        bool reachable_right = distances_to_parent.second <= 
            (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit);
        //If this cluster wasn't combined and hasn't been seen before and its reachable from other clusters, add it to the parent
        if (!combined && parent_clusters.read_cluster_heads.count(make_pair(read_num, cluster_head)) == 0 &&
            (reachable_left || reachable_right)) {
            parent_clusters.read_cluster_heads.emplace(read_num, cluster_head);
            tree_state.read_cluster_heads_to_distances[read_num][cluster_head] = distances_to_parent;
        }
    }

    /*Now go through clusters on the second child, and see if they can be combined with clusters on the first child
     */
    for (auto& child_cluster_head : child_clusters2.read_cluster_heads) {

        size_t read_num = child_cluster_head.first;
        size_t cluster_num = tree_state.read_union_find[read_num].find_group(child_cluster_head.second);
        pair<size_t, size_t> distances = child_distances[read_num][child_cluster_head.second];
        size_t new_dist_left = std::min(SnarlDistanceIndex::sum({distances.first,child_clusters2.distance_start_left}), 
                                        SnarlDistanceIndex::sum({distances.second,child_clusters2.distance_start_right}));
        size_t new_dist_right = std::min(SnarlDistanceIndex::sum({distances.first,child_clusters2.distance_end_left}), 
                                        SnarlDistanceIndex::sum({distances.second,child_clusters2.distance_end_right}));
        pair<size_t, size_t> distances_to_parent = make_pair(new_dist_left, new_dist_right);

        if (parent_clusters.read_cluster_heads.count(make_pair(read_num, cluster_num)) > 0) {
            distances_to_parent = make_pair(
                std::min(new_dist_left, tree_state.read_cluster_heads_to_distances[read_num][cluster_num].first),
                std::min(new_dist_right,tree_state.read_cluster_heads_to_distances[read_num][cluster_num].second));
        }

        //Check if the left of 1 can connect with the left of 2
        compare_and_combine_clusters (read_num, cluster_num, 
            SnarlDistanceIndex::sum({distances.first,distance_left_left,child_clusters1.read_best_left[read_num]}), 
            SnarlDistanceIndex::sum({distances.first,distance_left_left,child_clusters1.fragment_best_left}), 
            distances_to_parent, new_cluster_left_left_by_read[read_num], new_cluster_left_left_fragment);

        //Check if the left of 1 can connect with the right of 2
        compare_and_combine_clusters (read_num, cluster_num, 
            SnarlDistanceIndex::sum({distances.second,distance_left_right,child_clusters1.read_best_left[read_num]}),
            SnarlDistanceIndex::sum({distances.second,distance_left_right,child_clusters1.fragment_best_left}),
            distances_to_parent, new_cluster_left_right_by_read[read_num], new_cluster_left_right_fragment);

        //Check if the right of 1 can connect with the right of 2
        compare_and_combine_clusters (read_num, cluster_num, 
            SnarlDistanceIndex::sum({distances.second,distance_right_right,child_clusters1.read_best_right[read_num]}),
            SnarlDistanceIndex::sum({distances.second,distance_right_right,child_clusters1.fragment_best_right}),
            distances_to_parent, new_cluster_right_right_by_read[read_num], new_cluster_right_right_fragment);

        //Check if the right of 1 can connect with the left of 2
        compare_and_combine_clusters (read_num, cluster_num, 
            SnarlDistanceIndex::sum({distances.first,distance_right_left,child_clusters1.read_best_right[read_num]}),
            SnarlDistanceIndex::sum({distances.first,distance_right_left,child_clusters1.fragment_best_right}),
            distances_to_parent, new_cluster_right_left_by_read[read_num], new_cluster_right_left_fragment);
    }

    /*then remove all clusters that got erase, then add back in the cluster heads
     */

    //remove cluster heads that got combined with new ones
    for (pair<size_t, size_t>& cluster_head : to_erase) {
        parent_clusters.read_cluster_heads.erase(cluster_head);
        //tree_state.read_cluster_heads_to_distances[cluster_head.first][cluster_head.scond] =std::numeric_limits<size_t>::max();
    }

    //And add back in the new cluster heads
    for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
        if (new_cluster_left_left_by_read.at(read_num).first.first != std::numeric_limits<size_t>::max()) {
            pair<size_t, size_t> old_distances = parent_clusters.read_cluster_heads.count(new_cluster_left_left_by_read.at(read_num).first) == 0
                ? make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max())
                : tree_state.read_cluster_heads_to_distances[read_num][new_cluster_left_left_by_read.at(read_num).first.second];

            //Is the distance small enough that we can cluster it with something else?
            size_t best_left = std::min(new_cluster_left_left_by_read.at(read_num).second.first, old_distances.first);
            size_t best_right = std::min(new_cluster_left_left_by_read.at(read_num).second.second, old_distances.second);
            bool reachable_left = best_left <= 
                (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit);
            bool reachable_right = best_right <= 
                (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit);

            if (reachable_left || reachable_right) {
                parent_clusters.read_cluster_heads.emplace(new_cluster_left_left_by_read.at(read_num).first);
                tree_state.read_cluster_heads_to_distances[read_num][new_cluster_left_left_by_read.at(read_num).first.second] = 
                    make_pair(best_left,best_right);
            } else {
                parent_clusters.read_cluster_heads.erase(new_cluster_left_left_by_read.at(read_num).first);
                //tree_state.read_cluster_heads_to_distances[read_num][new_cluster_left_left_by_read.at(read_num).first] = std::numeric_limits<size_t>::max();
            }
        }
        if (new_cluster_right_right_by_read.at(read_num).first.first != std::numeric_limits<size_t>::max()) {
            pair<size_t, size_t> old_distances = parent_clusters.read_cluster_heads.count(new_cluster_right_right_by_read.at(read_num).first) == 0
                ? make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max())
                : tree_state.read_cluster_heads_to_distances[read_num][new_cluster_right_right_by_read.at(read_num).first.second];

            size_t best_left = std::min(new_cluster_right_right_by_read.at(read_num).second.first, old_distances.first);
            size_t best_right = std::min(new_cluster_right_right_by_read.at(read_num).second.second, old_distances.second);
            bool reachable_left = best_left <= 
                (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit);
            bool reachable_right = best_right <= 
                (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit);

            if (reachable_left || reachable_right) {
                parent_clusters.read_cluster_heads.emplace(new_cluster_right_right_by_read.at(read_num).first);
                tree_state.read_cluster_heads_to_distances[read_num][new_cluster_right_right_by_read.at(read_num).first.second] = 
                    make_pair(best_left,best_right);
            } else {
                parent_clusters.read_cluster_heads.erase(new_cluster_right_right_by_read.at(read_num).first);
                //tree_state.read_cluster_heads_to_distances.erase(new_cluster_right_right_by_read.at(read_num).first);
            }
        }
        if (new_cluster_left_right_by_read.at(read_num).first.first != std::numeric_limits<size_t>::max()) {
            pair<size_t, size_t> old_distances = parent_clusters.read_cluster_heads.count(new_cluster_left_right_by_read.at(read_num).first) == 0
                ? make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max())
                : tree_state.read_cluster_heads_to_distances[read_num][new_cluster_left_right_by_read.at(read_num).first.second];

            size_t best_left = std::min(new_cluster_left_right_by_read.at(read_num).second.first, old_distances.first);
            size_t best_right = std::min(new_cluster_left_right_by_read.at(read_num).second.second, old_distances.second);
            bool reachable_left = best_left <= 
                (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit);
            bool reachable_right = best_right <= 
                (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit);

            if (reachable_left || reachable_right) {
                parent_clusters.read_cluster_heads.emplace(new_cluster_left_right_by_read.at(read_num).first);
                tree_state.read_cluster_heads_to_distances[read_num][new_cluster_left_right_by_read.at(read_num).first.second] = 
                    make_pair(best_left, best_right);
            } else {
                parent_clusters.read_cluster_heads.erase(new_cluster_left_right_by_read.at(read_num).first);
                //tree_state.read_cluster_heads_to_distances.erase(new_cluster_left_right_by_read.at(read_num).first);
            }
        }
        if (new_cluster_right_left_by_read.at(read_num).first.first != std::numeric_limits<size_t>::max()) {
            pair<size_t, size_t> old_distances = parent_clusters.read_cluster_heads.count(new_cluster_right_left_by_read.at(read_num).first) == 0
                ? make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max())
                : tree_state.read_cluster_heads_to_distances[read_num][new_cluster_right_left_by_read.at(read_num).first.second];
            size_t best_left = std::min(new_cluster_right_left_by_read.at(read_num).second.first, old_distances.first);
            size_t best_right = std::min(new_cluster_right_left_by_read.at(read_num).second.second, old_distances.second);
            bool reachable_left = best_left <= 
                (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit);
            bool reachable_right = best_right <= 
                (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit);

            if (reachable_left || reachable_right) {
                parent_clusters.read_cluster_heads.emplace(new_cluster_right_left_by_read.at(read_num).first);
                tree_state.read_cluster_heads_to_distances[read_num][new_cluster_right_left_by_read.at(read_num).first.second] = 
                    make_pair(best_left, best_right);
            } else {
                parent_clusters.read_cluster_heads.erase(new_cluster_right_left_by_read.at(read_num).first);
                //tree_state.read_cluster_heads_to_distances.erase(new_cluster_right_left_by_read.at(read_num).first);
            }
        }
    }


    /*Update the parent's best left and right distances, only looking at the first child since we've already seen the second one
     */
    //Update the parent's fragment best distances
    parent_clusters.fragment_best_left = std::min(parent_clusters.fragment_best_left,
                                         std::min(SnarlDistanceIndex::sum({child_clusters1.distance_start_left, child_clusters1.fragment_best_left}),
                                                  SnarlDistanceIndex::sum({child_clusters1.distance_start_right , child_clusters1.fragment_best_right})));
    parent_clusters.fragment_best_right = std::min(parent_clusters.fragment_best_right,
                                          std::min(SnarlDistanceIndex::sum({child_clusters1.distance_end_left , child_clusters1.fragment_best_left}),
                                                   SnarlDistanceIndex::sum({child_clusters1.distance_end_right , child_clusters1.fragment_best_right})));


    //Update the best distances in the parent for each read num
    for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num ++) {
        //Find the best distances to the ends of the parent from child1
        size_t best_start = std::min(SnarlDistanceIndex::sum({child_clusters1.distance_start_left,child_clusters1.read_best_left.at(read_num)}), 
                                    SnarlDistanceIndex::sum({child_clusters1.distance_start_right,child_clusters1.read_best_right.at(read_num)}));
        size_t best_end = std::min(SnarlDistanceIndex::sum({child_clusters1.distance_end_left,child_clusters1.read_best_left.at(read_num)}), 
                                    SnarlDistanceIndex::sum({child_clusters1.distance_end_right,child_clusters1.read_best_right.at(read_num)}));
        //And update the distances in the parent
        parent_clusters.read_best_left.at(read_num) = std::min(best_start, parent_clusters.read_best_left.at(read_num));
        parent_clusters.read_best_right.at(read_num) = std::min(best_end, parent_clusters.read_best_right.at(read_num));
    }
    
//#ifdef DEBUG_CLUSTER
//        cerr << "\tIntermediate clusters on " << distance_index.net_handle_as_string(parent_clusters.containing_net_handle);
//        cerr << "   with best left and right values: " << parent_clusters.fragment_best_left << " "
//             << parent_clusters.fragment_best_right << endl;
//        bool got_left = false;
//        bool got_right = false;
//        for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
//            cerr << "\t\t\tfor read num " << read_num << " best left: " << parent_clusters.read_best_left[read_num] << " best right: " << parent_clusters.read_best_right[read_num] << endl;
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

void NewSnarlSeedClusterer::compare_and_combine_cluster_on_one_child(TreeState& tree_state, NodeClusters& child_clusters) const {
#ifdef DEBUG_CLUSTER
    cerr << "\tCompare " << distance_index.net_handle_as_string(child_clusters.containing_net_handle) 
         << " to itself in the root" << endl;
#endif

    net_handle_t& handle = child_clusters.containing_net_handle;


    //Get the distances between the two sides of the child
    size_t distance_left_left = distance_index.is_externally_start_start_connected(handle) ? 0 : std::numeric_limits<size_t>::max();
    size_t distance_left_right = distance_index.is_externally_start_end_connected(handle) ? 0 : std::numeric_limits<size_t>::max();
    size_t distance_right_right = distance_index.is_externally_end_end_connected(handle) ? 0 : std::numeric_limits<size_t>::max();

#ifdef DEBUG_CLUSTER
    cerr << "\t\tFound distances between the two children: " << distance_left_left << " " << distance_left_right << " " << distance_right_right <<  endl;
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

    //These will be the cluster heads and distances of everything combined by taking the indicated path
    //one cluster head per read
    //pair< pair<read_num, seed_num>, pair<left dist, right_dist>>
    //The default value will be ((inf, 0), (0,0)). Only the inf gets checked to see if it's a real value so I filled it in with 0 so I wouldn't have to type out inf 
    vector<pair<size_t, size_t>> new_cluster_left_left_by_read(tree_state.all_seeds->size(), 
            make_pair(std::numeric_limits<size_t>::max(), 0));
    vector<pair<size_t, size_t>> new_cluster_left_right_by_read(tree_state.all_seeds->size(),
        make_pair(std::numeric_limits<size_t>::max(), 0));
    vector<pair<size_t, size_t>> new_cluster_right_right_by_read(tree_state.all_seeds->size(),
        make_pair(std::numeric_limits<size_t>::max(), 0));

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
            return false;
        }
        distance_between_reads = SnarlDistanceIndex::minus(distance_between_reads, 1);
        distance_between_fragments = SnarlDistanceIndex::minus(distance_between_fragments, 1);
        bool combined = false;

        if (distance_between_reads <= tree_state.read_distance_limit) {
            //If this can be combined with the given combined cluster
            if (new_cluster_head.first == std::numeric_limits<size_t>::max()){
                //new cluster head
                new_cluster_head = make_pair(read_num, cluster_num);
            } else {
                //Combine with old cluster head
                tree_state.read_union_find.at(read_num).union_groups(cluster_num, new_cluster_head.second);

                //Update distances

                //And remember new head and distances
                new_cluster_head = make_pair(read_num, tree_state.read_union_find.at(read_num).find_group(cluster_num));
            }
            //Remember to erase this cluster head
            combined = true;

#ifdef DEBUG_CLUSTER
            cerr << "\t\t\tCombining read/cluster " << read_num << "/" << cluster_num << "... new cluster head:" << new_cluster_head.second << endl; 
#endif
        }
        if (tree_state.fragment_distance_limit != 0 && 
                    distance_between_fragments <= tree_state.fragment_distance_limit ) {
            //Just union the fragment
            if (new_cluster_head_fragment == std::numeric_limits<size_t>::max()) {
                new_cluster_head_fragment =cluster_num+tree_state.read_index_offsets[read_num];
            } else {
                tree_state.fragment_union_find.union_groups(cluster_num+tree_state.read_index_offsets[read_num], 
                                                                    new_cluster_head_fragment);
                new_cluster_head_fragment = tree_state.fragment_union_find.find_group(new_cluster_head_fragment);
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
    for (auto& child_cluster_head : child_clusters.read_cluster_heads) {

        bool combined = false;
        size_t read_num = child_cluster_head.first;
        size_t cluster_num = tree_state.read_union_find[read_num].find_group(child_cluster_head.second);

        //Distances to the ends of the child
        pair<size_t, size_t> distances = tree_state.read_cluster_heads_to_distances[read_num][child_cluster_head.second];
            

        //Check if the left of 1 can connect with the left of 2
        combined = combined | compare_and_combine_clusters (read_num, cluster_num, 
            SnarlDistanceIndex::sum({distances.first,distance_left_left, child_clusters.read_best_left[read_num]}), 
            SnarlDistanceIndex::sum({distances.first,distance_left_left, child_clusters.fragment_best_left}), 
            new_cluster_left_left_by_read[read_num], new_cluster_left_left_fragment);

        //Check if the left of 1 can connect with the right of 2
        combined = combined | compare_and_combine_clusters (read_num, cluster_num, 
            SnarlDistanceIndex::sum({distances.first,distance_left_right, child_clusters.read_best_right[read_num]}), 
            SnarlDistanceIndex::sum({distances.first,distance_left_right, child_clusters.fragment_best_right}), 
             new_cluster_left_right_by_read[read_num], new_cluster_left_right_fragment);

        //Check if the right of 1 can connect with the right of 2
        combined = combined | compare_and_combine_clusters (read_num, cluster_num, 
            SnarlDistanceIndex::sum({distances.second,distance_right_right, child_clusters.read_best_right[read_num]}), 
            SnarlDistanceIndex::sum({distances.second,distance_right_right, child_clusters.fragment_best_right}), 
            new_cluster_right_right_by_read[read_num], new_cluster_right_right_fragment);

        size_t cluster_head = tree_state.read_union_find[read_num].find_group(cluster_num); 

    }
//
//    /*Now go through clusters on the second child, and see if they can be combined with clusters on the first child
//     */
//    for (auto& child_cluster_head : child_clusters2.read_cluster_heads) {
//
//        size_t read_num = child_cluster_head.first;
//        size_t cluster_num = tree_state.read_union_find[read_num].find_group(child_cluster_head.second);
//        pair<size_t, size_t> distances = child_distances[read_num][child_cluster_head.second];
//        size_t new_dist_left = std::min(SnarlDistanceIndex::sum({distances.first,child_clusters2.distance_start_left}), 
//                                        SnarlDistanceIndex::sum({distances.second,child_clusters2.distance_start_right}));
//        size_t new_dist_right = std::min(SnarlDistanceIndex::sum({distances.first,child_clusters2.distance_end_left}), 
//                                        SnarlDistanceIndex::sum({distances.second,child_clusters2.distance_end_right}));
//        pair<size_t, size_t> distances_to_parent = make_pair(new_dist_left, new_dist_right);
//
//        if (parent_clusters.read_cluster_heads.count(make_pair(read_num, cluster_num)) > 0) {
//            distances_to_parent = make_pair(
//                std::min(new_dist_left, tree_state.read_cluster_heads_to_distances[read_num][cluster_num].first),
//                std::min(new_dist_right,tree_state.read_cluster_heads_to_distances[read_num][cluster_num].second));
//        }
//
//        //Check if the left of 1 can connect with the left of 2
//        compare_and_combine_clusters (read_num, cluster_num, 
//            SnarlDistanceIndex::sum({distances.first,distance_left_left,child_clusters1.read_best_left[read_num]}), 
//            SnarlDistanceIndex::sum({distances.first,distance_left_left,child_clusters1.fragment_best_left}), 
//            distances_to_parent, new_cluster_left_left_by_read[read_num], new_cluster_left_left_fragment);
//
//        //Check if the left of 1 can connect with the right of 2
//        compare_and_combine_clusters (read_num, cluster_num, 
//            SnarlDistanceIndex::sum({distances.second,distance_left_right,child_clusters1.read_best_left[read_num]}),
//            SnarlDistanceIndex::sum({distances.second,distance_left_right,child_clusters1.fragment_best_left}),
//            distances_to_parent, new_cluster_left_right_by_read[read_num], new_cluster_left_right_fragment);
//
//        //Check if the right of 1 can connect with the right of 2
//        compare_and_combine_clusters (read_num, cluster_num, 
//            SnarlDistanceIndex::sum({distances.second,distance_right_right,child_clusters1.read_best_right[read_num]}),
//            SnarlDistanceIndex::sum({distances.second,distance_right_right,child_clusters1.fragment_best_right}),
//            distances_to_parent, new_cluster_right_right_by_read[read_num], new_cluster_right_right_fragment);
//
//        //Check if the right of 1 can connect with the left of 2
//        compare_and_combine_clusters (read_num, cluster_num, 
//            SnarlDistanceIndex::sum({distances.first,distance_right_left,child_clusters1.read_best_right[read_num]}),
//            SnarlDistanceIndex::sum({distances.first,distance_right_left,child_clusters1.fragment_best_right}),
//            distances_to_parent, new_cluster_right_left_by_read[read_num], new_cluster_right_left_fragment);
//    }

}


void NewSnarlSeedClusterer::cluster_one_snarl(TreeState& tree_state, NodeClusters& snarl_clusters, const NodeClusters& parent_clusters) const { 
    //Get the clusters on this snarl, assumes that all of the snarls children have been clustered already.
    
#ifdef DEBUG_CLUSTER
        cerr << "Finding clusters on snarl " << distance_index.net_handle_as_string(snarl_clusters.containing_net_handle) << endl;
#endif

    //Keep track of all clusters on this snarl
    net_handle_t& snarl_handle = snarl_clusters.containing_net_handle;

    //Get the children of this snarl and their clusters
    vector<size_t> children = std::move(tree_state.snarl_to_children[snarl_handle].second);
    //The distances within the children, since we will be updating read_cluster_head_to_distances
    //to represent distances in the parent
    vector<vector<pair<size_t, size_t>>> child_distances (tree_state.all_seeds->size());
    for (size_t i = 0 ; i < tree_state.all_seeds->size() ; i++) {
        child_distances[i] = std::vector<pair<size_t, size_t>>(tree_state.all_seeds->at(i)->size(),
            make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
    }

    for (size_t i = 0; i < children.size() ; i++) {
        //Go through each child node of the netgraph

        //Remember the distances for this child since they will get overwritten
        NodeClusters& child_clusters_i = tree_state.all_node_clusters[children[i]];
        if (child_clusters_i.fragment_best_left > (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit) &&  
            child_clusters_i.fragment_best_right > (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit)) {
            continue;
        }

        for (const pair<size_t, size_t>& head : child_clusters_i.read_cluster_heads) {
            child_distances[head.first][head.second] = tree_state.read_cluster_heads_to_distances[head.first][head.second];
        }

        for (size_t j = 0 ; j <= i ; j++){
            //Go through other child net graph nodes up to and including i

            //Get the other node and its clusters
            NodeClusters& child_clusters_j = tree_state.all_node_clusters[children[j]];

            if (child_clusters_j.fragment_best_left > (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit) &&  
                child_clusters_j.fragment_best_right > (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit)) {
                continue;
            }

#ifdef DEBUG_CLUSTER
            cerr << "\tComparing two children of " << distance_index.net_handle_as_string(snarl_handle) << ": " 
                 << distance_index.net_handle_as_string(child_clusters_i.containing_net_handle) << " and " 
                 << distance_index.net_handle_as_string(child_clusters_j.containing_net_handle) << endl;
                 


#endif

            compare_and_combine_cluster_on_child_structures(tree_state, child_clusters_i, 
                    child_clusters_j, snarl_clusters, child_distances, false);

        }
    }
    //Get the values needed for clustering a chain
    //Prefix sum up to the start of the snarl
    snarl_clusters.prefix_sum_value = SnarlDistanceIndex::sum({
            distance_index.get_prefix_sum_value(distance_index.get_node_from_sentinel(snarl_clusters.start_in)),
            distance_index.minimum_length(      distance_index.get_node_from_sentinel(snarl_clusters.start_in))});

    //Length of the snarl
    snarl_clusters.node_length = distance_index.minimum_length(snarl_clusters.containing_net_handle);

    //Chain component of the start and end nodes of the snarl
    if (distance_index.is_multicomponent_chain(distance_index.get_parent(snarl_clusters.containing_net_handle))) {
        snarl_clusters.chain_component_start = distance_index.get_chain_component(distance_index.get_node_from_sentinel(distance_index.get_bound(snarl_clusters.containing_net_handle, false, false)));
        snarl_clusters.chain_component_end = distance_index.get_chain_component(distance_index.get_node_from_sentinel(distance_index.get_bound(snarl_clusters.containing_net_handle, true, false)));
    }
#ifdef DEBUG_CLUSTER
//    cerr << "\tFound clusters on " << distance_index.net_handle_as_string(snarl_handle) << endl;
//    cerr << "\t   with best left and right values: " << snarl_clusters.fragment_best_left << " "
//         << snarl_clusters.fragment_best_right << endl;
//    bool got_left = false;
//    bool got_right = false;
//    for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
//        cerr << "\t\tfor read num " << read_num << " best left: " << snarl_clusters.read_best_left[read_num] << " best right: " << snarl_clusters.read_best_right[read_num] << endl;
//        bool got_read_left=false;
//        bool got_read_right = false;
//        bool any_clusters = false;
//        for (pair<size_t,size_t> c : snarl_clusters.read_cluster_heads) {
//            if (c.first == read_num) {
//                any_clusters = true;
//                pair<size_t, size_t> dists = tree_state.read_cluster_heads_to_distances[c.first][c.second];
//                cerr << "\t\t" << c.first << ":"<<c.second << ": left: " << dists.first << " right : " << dists.second << ": ";
//                bool has_seeds = false;
//                for (size_t x = 0 ; x < tree_state.all_seeds->at(c.first)->size() ; x++) {
//                    if (tree_state.read_union_find[c.first].find_group(x) == c.second) {
//                        cerr << tree_state.all_seeds->at(c.first)->at(x).pos << " ";
//                        has_seeds = true;
//                    }
//                }
//                assert(dists.first == std::numeric_limits<size_t>::max() || dists.first >= snarl_clusters.read_best_left[read_num]);
//                assert(dists.second == std::numeric_limits<size_t>::max() || dists.second >= snarl_clusters.read_best_right[read_num]);
//                assert(dists.first == std::numeric_limits<size_t>::max() || dists.first >= snarl_clusters.fragment_best_left);
//                assert(dists.second == std::numeric_limits<size_t>::max() || dists.second >= snarl_clusters.fragment_best_right);
//                if (dists.first == snarl_clusters.fragment_best_left) {got_left = true;}
//                if (dists.second == snarl_clusters.fragment_best_right) {got_right = true;}
//                if (dists.first == snarl_clusters.read_best_left[read_num]) {got_read_left = true;}
//                if (dists.second == snarl_clusters.read_best_right[read_num]) {got_read_right = true;}
//                cerr << endl;
//                assert(has_seeds);
//            }
//        }
//        //assert(!any_clusters ||got_read_left ||  snarl_clusters.read_best_left[read_num] == std::numeric_limits<size_t>::max());
//        //assert(!any_clusters ||got_read_right ||  snarl_clusters.read_best_right[read_num] == std::numeric_limits<size_t>::max());
//    }
//    //assert(got_left);
//    //assert(got_right);
//
//    //for (pair<pair<size_t, size_t>, pair<size_t, size_t>> group_id : snarl_clusters.read_cluster_heads) {
//    //    assert (group_id.first.second == tree_state.read_union_find[group_id.first.first].find_group(group_id.first.second));
//    //}
#endif
};



void NewSnarlSeedClusterer::cluster_one_chain(TreeState& tree_state, NodeClusters& chain_clusters) const {
#ifdef DBUG_CLUSTERS
    assert(distance_index.is_chain(chain_clusters.containing_net_handle));
#endif
    
    net_handle_t& chain_handle = chain_clusters.containing_net_handle;
    //Note: This is the first time we deal with the net_handle_t for this chain, so it will be empty
    if (tree_state.chain_to_children->at(chain_handle).second.size() == 1 &&
        SnarlDistanceIndex::get_record_offset(chain_handle) == 
        SnarlDistanceIndex::get_record_offset(tree_state.all_node_clusters[tree_state.chain_to_children->at(chain_handle).second.back()].containing_net_handle)){
        //If this is just a node pretending to be a chain, cluster the node and claim it's a chain
#ifdef DEBUG_CLUSTER
        cerr << "\tClustering a chain that is really just a node, so just cluster the node" << endl;
#endif

        
        //This will hold the clusters of the chain
        assert(tree_state.chain_to_children->at(chain_handle).second.size() == 1); 


        NodeClusters& node_clusters = tree_state.all_node_clusters[tree_state.chain_to_children->at(chain_handle).second.back()];
        cluster_one_node(tree_state, node_clusters);

        //Copy everything from the node
        chain_clusters.read_cluster_heads = std::move(node_clusters.read_cluster_heads);
        chain_clusters.fragment_best_left = std::move(node_clusters.fragment_best_left);
        chain_clusters.fragment_best_right = std::move(node_clusters.fragment_best_right);
        chain_clusters.read_best_left = std::move(node_clusters.read_best_left);
        chain_clusters.read_best_right = std::move(node_clusters.read_best_right);
        return;

    }

    /*Define a helper function to update the distances in a child using the loop distances
     * in the chain
     */

    //These are used in update_distance_on_same_child but I"m declaring them here because I seem to be spending
    //a lot of time allocating and deallocating memory and I want to do it just once
    size_t loop_right, loop_left, combined_fragment;
    vector<pair<pair<size_t, size_t>, pair<size_t, size_t>>> combined_left,combined_right; 
    vector<pair<size_t, size_t>> to_erase;
    size_t read_num, cluster_num, old_left, old_right, updated_left, updated_right, distance_between_left, 
           distance_between_right, distance_between_left_fragment, distance_between_right_fragment;
    pair<size_t, size_t> cluster_head;

     auto update_distances_on_same_child = [&] (NodeClusters& child_clusters) {
         //Distance to go forward (relative to the child) in the chain and back
         loop_right = distance_index.distance_in_parent(chain_handle, child_clusters.containing_net_handle, 
                                  child_clusters.containing_net_handle);
         //Distance to go backward in the chain and back
         loop_left = distance_index.distance_in_parent(chain_handle, distance_index.flip(child_clusters.containing_net_handle), 
                                  distance_index.flip(child_clusters.containing_net_handle));


         //Combined clusters in case we can combine anything
         combined_left.assign(tree_state.all_seeds->size(),
                make_pair(make_pair(std::numeric_limits<size_t>::max(), 0), make_pair(0,0)));
         combined_right.assign(tree_state.all_seeds->size(),
                make_pair(make_pair(std::numeric_limits<size_t>::max(), 0), make_pair(0,0))); 
         combined_fragment = std::numeric_limits<size_t>::max();
         to_erase.clear();

         pair<size_t, size_t> cluster_head;
         for (auto& child_cluster_head : child_clusters.read_cluster_heads) {
            //Go through each of the clusters on this child
            read_num = child_cluster_head.first;
            cluster_num = tree_state.read_union_find[read_num].find_group(child_cluster_head.second);
            old_left = tree_state.read_cluster_heads_to_distances[read_num][cluster_num].first;
            old_right = tree_state.read_cluster_heads_to_distances[read_num][cluster_num].second;
            //Get the new best distances for the cluster considering chain loops
            updated_left = std::min(old_left, SnarlDistanceIndex::sum({old_right, loop_right, distance_index.minimum_length(child_clusters.containing_net_handle)}));
            updated_right = std::min(old_right, SnarlDistanceIndex::sum({old_left, loop_left, distance_index.minimum_length(child_clusters.containing_net_handle)}));



            if (updated_left < old_left || updated_right < old_right ) {
                //Update the distances
                tree_state.read_cluster_heads_to_distances[read_num][cluster_num] = 
                        std::make_pair(updated_left, updated_right);
                
                child_clusters.fragment_best_left = std::min(child_clusters.fragment_best_left,
                                                             updated_left);
                child_clusters.fragment_best_right = std::min(child_clusters.fragment_best_right,
                                                             updated_right);
                child_clusters.read_best_left[read_num] = std::min(child_clusters.read_best_left[read_num],
                                                                   updated_left);
                child_clusters.read_best_right[read_num] = std::min(child_clusters.read_best_right[read_num],
                                                                       updated_right);
            }

            //Now see if we can combine this cluster with anything else 
            //The distance between this cluster and anything else taking the left loop
            distance_between_left = SnarlDistanceIndex::minus(
                    SnarlDistanceIndex::sum({updated_left, 
                                             loop_left,
                                             child_clusters.read_best_left[read_num]}),
                    1); 
            distance_between_right = SnarlDistanceIndex::minus(
                     SnarlDistanceIndex::sum({updated_right, 
                                              loop_right,
                                              child_clusters.read_best_right[read_num]}),
                     1); 
            distance_between_left_fragment = SnarlDistanceIndex::minus(
                      SnarlDistanceIndex::sum({updated_left, 
                                                loop_left,
                                                child_clusters.fragment_best_left}),
                      1); 
            distance_between_right_fragment = SnarlDistanceIndex::minus(
                       SnarlDistanceIndex::sum({updated_right, 
                                                loop_right,
                                                child_clusters.fragment_best_right}),
                       1); 
            cluster_head = make_pair(read_num, cluster_num);
            if (distance_between_left <= tree_state.read_distance_limit) {
                //Combine it left
                to_erase.emplace_back(cluster_head);
                if (combined_left[read_num].first.first == std::numeric_limits<size_t>::max()){
                    combined_left[read_num] = make_pair(cluster_head, 
                         make_pair(updated_left, updated_right));
                } else {
                    to_erase.emplace_back(combined_left[read_num].first);
                    tree_state.read_union_find.at(read_num).union_groups(cluster_num, combined_left[read_num].first.second);
                    combined_left[read_num] = make_pair(
                         make_pair(read_num, 
                            tree_state.read_union_find.at(read_num).find_group(cluster_num)), 
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
                    tree_state.read_union_find.at(read_num).union_groups(cluster_num, combined_right[read_num].first.second);
                    combined_right[read_num] =make_pair(
                         make_pair(read_num, 
                            tree_state.read_union_find.at(read_num).find_group(cluster_num)), 
                         make_pair(std::min(updated_left, combined_right[read_num].second.first), 
                                   std::min(updated_right, combined_right[read_num].second.second))); 
                    
                }
            }
            if (tree_state.fragment_distance_limit != 0 &&
                (distance_between_right_fragment <= tree_state.fragment_distance_limit
                 || distance_between_left_fragment <= tree_state.fragment_distance_limit)) {
                //Combine the fragment
                if (combined_fragment != std::numeric_limits<size_t>::max()) {
                    tree_state.fragment_union_find.union_groups(combined_fragment, 
                            cluster_num + tree_state.read_index_offsets[read_num]);
                }
                combined_fragment = cluster_num + tree_state.read_index_offsets[read_num];
            }
         }
        for (pair<size_t, size_t>& cluster_head : to_erase) {
            child_clusters.read_cluster_heads.erase(cluster_head);
            //tree_state.read_cluster_heads_to_distances.erase(cluster_head);
        }
        //Add new clusters that were combined
        for (pair<pair<size_t, size_t>,pair<size_t, size_t>>& cluster : combined_left) {
            if (cluster.first.first != std::numeric_limits<size_t>::max()){
                child_clusters.read_cluster_heads.emplace(cluster.first);
                tree_state.read_cluster_heads_to_distances[cluster.first.first][cluster.first.second] = cluster.second;
            }
        }
        for (pair<pair<size_t, size_t>, pair<size_t, size_t>>& cluster : combined_right) {
            if (cluster.first.first != std::numeric_limits<size_t>::max()){
                child_clusters.read_cluster_heads.emplace(cluster.first);
                tree_state.read_cluster_heads_to_distances[cluster.first.first][cluster.first.second] = cluster.second;
            }
        }
    };


#ifdef DEBUG_CLUSTER
    cerr << "Cluster chain " << distance_index.net_handle_as_string(chain_handle) << endl;
    cerr << "\t chain has " << tree_state.chain_to_children->at(chain_handle).second.size() << " children" << endl;
#endif

    /*Go through the chain child by child 
     *
     * As we walk through the chain, keep track of all clusters found up to the current child. 
     * So after we saw a child, the chain knows all clusters with distances up to the right side of the child
    * 
    * For each child, 
    *  - check if the clusters in the child can be combined with each other by walking out then back in through the chain
    *       -Update distances to the ends of the child (taking into account the distances to loop around in the chain)
    *  - compare and combine with clusters of the chain that get build as we walk along it
    * 
    *  - after combining clusters of the current child, remove redundant cluster heads from the chain clusters
    * 
    *  Note that every snarl cluster must be compared to the current child, even if there was a snarl cluster that came
    *  after it in the chain (but not if there was a node cluster after it). If there were clusters on snarls 1, 2, and 4
    *  in the chain, and 1 got combined with 2, the minimum distance from 1 to 4 might still be smaller than the minimum 
    *  distance from 2 to 4, so comparing only 2 to 4 might miss the distance from 1 to 4.
    *  If multiple clusters in the same snarl get combined, then they the redundant cluster head can be removed
    * 
    */

    //Get the children of this chain from the tree state. They will be ordered by their order in the chain
    vector<size_t>& children_in_chain = tree_state.chain_to_children->at(chain_handle).second;

    /*
     * As we walk along the chain, we maintain clusters of the chain up to the last node we saw (the later
     * boundary node of a snarl relative to the chain)
     *
     * The clusters in the chain have left distances to the beginning of the chain and right distances
     * to the last thing we saw
     */

    //The last node we saw is initialized to the first node in the chain
    NodeClusters& last_child = tree_state.all_node_clusters[children_in_chain.front()];
    //These are clusters that we don't want to consider as we walk through the chain but that 
    //we want to remember after we're done with the chain because the left distance is small
    vector<pair<pair<size_t, size_t>, pair<size_t, size_t>>> cluster_heads_to_add_again;

   
    //These get cleared every loop but I"m declaring them here so I don't keep reallocating memory
    vector<pair<pair<size_t, size_t>, pair<size_t, size_t>>> to_add;
    vector<pair<pair<size_t, size_t>, pair<size_t, size_t>>> new_cluster_by_read;
    vector<size_t> old_best_left_by_read;
    vector<size_t> old_best_right_by_read;
    size_t new_cluster_head_fragment = std::numeric_limits<size_t>::max();
    bool child_is_reversed;
    size_t old_best_left;
    size_t old_best_right;

    for (size_t i = 0 ; i < children_in_chain.size() ; i++) {
        /*
         * Snarls and nodes are in the order that they are traversed in the chain
         */
        size_t& child_clusters_i = children_in_chain[i];

        //The NodeClusters we're currently looking at
        NodeClusters& child_clusters = tree_state.all_node_clusters[child_clusters_i];
        net_handle_t& child_handle = child_clusters.containing_net_handle;

        //If cluster the children here
        if (distance_index.is_node(child_handle)){
            cluster_one_node(tree_state, child_clusters);
        } else {
            cluster_one_snarl(tree_state, child_clusters, chain_clusters);
        }

        //Skip this child if its seeds are all too far away
        if (child_clusters.fragment_best_left > (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit) &&  
            child_clusters.fragment_best_right > (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit)) {
            continue;
        }
#ifdef DEBUG_CLUSTER
        cerr << "At child " << distance_index.net_handle_as_string(child_handle) << endl;
#endif

        //See if this clusters of this child can be combined with each other
        //Also updates the minimum distances to include loops in the chain
        //This only matters for snarls, since any path would have to pass through a node anyway
        if (distance_index.is_snarl(child_clusters.containing_net_handle)){
            update_distances_on_same_child(child_clusters);
        }

        //The distance from the start of the chain to the left side of the current node, which is the
        //distance needed to update the clusters in the child to the start of the chain
        size_t distance_from_chain_start_to_current_node = child_clusters.prefix_sum_value;

        size_t distance_from_chain_start_to_last_node = SnarlDistanceIndex::sum({
            last_child.prefix_sum_value,last_child.node_length});

        //The distance from the right side of the last child to the left side of this child 
        //(relative to the orientation of the chain
        //If this is a looping chain, then find the distance normally. Otherwise use the prefix sums

        size_t distance_from_last_child_to_current_child = i==0 ? std::numeric_limits<size_t>::max() :
            ( chain_clusters.is_looping_chain 
                ? distance_index.distance_in_parent(chain_handle, 
                        distance_index.is_reversed_in_parent(last_child.containing_net_handle) ? distance_index.flip(last_child.containing_net_handle) : last_child.containing_net_handle,
                        distance_index.is_reversed_in_parent(child_handle) ? child_handle : distance_index.flip(child_handle))
                : (last_child.chain_component_end == child_clusters.chain_component_start) ?SnarlDistanceIndex::minus(distance_from_chain_start_to_current_node, distance_from_chain_start_to_last_node) : std::numeric_limits<size_t>::max());
            
        if (child_clusters.chain_component_start != 0) {
            distance_from_chain_start_to_current_node = std::numeric_limits<size_t>::max();
        }


        //The distance from the right side of the last child to the right side of this child, which is
        //the distance we need to update the chain clusters to the end of this child
        //This isn't quite right for the first thing in the chain but it doesn't matter because it only
        //gets added to chain clusters
        size_t distance_from_last_child_to_current_end = SnarlDistanceIndex::sum({distance_from_last_child_to_current_child ,
                 child_clusters.node_length});

        //The distance to add to get to the end of the chain. Only matters if this is the last thing in the chain
        //The distances will include the distance to the end of a trivial chain,
        //so we can't rely on distance_in_parent to know when the distance should be 0

        size_t distance_from_current_end_to_end_of_chain;
        if (i != children_in_chain.size() - 1 || 
                     SnarlDistanceIndex::get_record_offset(child_handle) == SnarlDistanceIndex::get_record_offset(chain_clusters.end_in)) {
            //If this is the last node in the chain
            distance_from_current_end_to_end_of_chain = 0;
        } else if (chain_clusters.is_looping_chain) {
            //If it's a looping chain then use the distance index
            distance_from_current_end_to_end_of_chain = distance_index.distance_in_parent(chain_handle, chain_clusters.end_in, 
                    (distance_index.is_reversed_in_parent(child_handle) ? distance_index.flip(child_handle) : child_handle));
        } else if (chain_clusters.chain_last_component != child_clusters.chain_component_end) { 
            //If they aren't in the same component
            distance_from_current_end_to_end_of_chain = std::numeric_limits<size_t>::max();
        } else {
            distance_from_current_end_to_end_of_chain = SnarlDistanceIndex::minus(chain_clusters.chain_length, 
                            SnarlDistanceIndex::sum({distance_from_chain_start_to_current_node, distance_index.minimum_length(child_handle)}));
        }

#ifdef DEBUG_CLUSTER
cerr << "\tDistance from last child to this one: " << distance_from_last_child_to_current_child << endl;
cerr << "\tDistance from start of chain to the left side of this one: " << distance_from_chain_start_to_current_node << endl;
cerr << "\tDistance from the last child to the right side of this one: " << distance_from_last_child_to_current_end << endl;
cerr << "\tDistance to get to the end of the chain: " << distance_from_current_end_to_end_of_chain << endl;
#endif

        //Clusters to remove from the chain because they got combined
        //And new clusters to add
        to_add.clear();

        //There is at most one new cluster per read
        new_cluster_by_read.assign(tree_state.all_seeds->size(), 
            make_pair(make_pair(std::numeric_limits<size_t>::max(), 0), make_pair(0,0)) );
        //And one new fragment cluster
        new_cluster_head_fragment = std::numeric_limits<size_t>::max();

        bool child_is_reversed = child_clusters.is_reversed_in_parent;

        //Remember the current best chain distances, and reset them to inf since we need to update them
        old_best_left = std::move(chain_clusters.fragment_best_left);
        chain_clusters.fragment_best_left = std::numeric_limits<size_t>::max();
        old_best_right = std::move(chain_clusters.fragment_best_right);
        chain_clusters.fragment_best_right = std::numeric_limits<size_t>::max(); 
        old_best_left_by_read = std::move(chain_clusters.read_best_left);
        chain_clusters.read_best_left.assign(old_best_left_by_read.size(), std::numeric_limits<size_t>::max());
        old_best_right_by_read = std::move(chain_clusters.read_best_right);
        chain_clusters.read_best_right.assign(old_best_right_by_read.size(), std::numeric_limits<size_t>::max());


        if (SnarlDistanceIndex::sum({distance_from_last_child_to_current_child, old_best_right}) 
              > (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit)) {
#ifdef DEBUG_CLUSTER
            cerr << "This child is too far away from the last one to cluster anything" << endl;
#endif
            //If the distance from the last cluster is too far to cluster anything
            for (auto& cluster_head : chain_clusters.read_cluster_heads) {
                //For each of the chain clusters, remember the ones that are still reachable from the left side of the chain
                pair<size_t, size_t> dists = tree_state.read_cluster_heads_to_distances[cluster_head.first][cluster_head.second];
                if (dists.first <= (tree_state.fragment_distance_limit == 0 ? tree_state.read_distance_limit : tree_state.fragment_distance_limit)) {
                    //If this cluster can be clustered outside of the chain, remember to add it back
                    cluster_heads_to_add_again.emplace_back(cluster_head, 
                                make_pair(dists.first, std::numeric_limits<size_t>::max()));
                }
            }

            //Now clear the chain's list of clusters
            chain_clusters.read_cluster_heads.clear();
            for (auto& cluster_head : child_clusters.read_cluster_heads) {
                //Add the clusters from this child to the chain
                size_t read_num = cluster_head.first;
                pair<size_t, size_t> dists = tree_state.read_cluster_heads_to_distances[read_num][cluster_head.second];
                size_t dist_left = child_clusters.is_reversed_in_parent ? dists.second : dists.first;
                size_t dist_right = child_clusters.is_reversed_in_parent ? dists.first : dists.second;

                //Distances to the start of the chain, and the end of this node
                //If this is the last thing in the chain, then the distance to the end of the chain
                pair<size_t, size_t> new_distances = make_pair(
                         SnarlDistanceIndex::sum({dist_left, distance_from_chain_start_to_current_node}),
                         SnarlDistanceIndex::sum({dist_right, distance_from_current_end_to_end_of_chain}));

                //Add this to the chain
                chain_clusters.read_cluster_heads.emplace(cluster_head); 
                tree_state.read_cluster_heads_to_distances[cluster_head.first][cluster_head.second]= new_distances; 
                //And update the best distances
                chain_clusters.fragment_best_left = std::min(chain_clusters.fragment_best_left, new_distances.first);
                chain_clusters.fragment_best_right = std::min(chain_clusters.fragment_best_right, new_distances.second); 
                chain_clusters.read_best_left[read_num] = std::min(chain_clusters.read_best_left[read_num], new_distances.first);
                chain_clusters.read_best_right[read_num] = std::min(chain_clusters.read_best_right[read_num], new_distances.second);
            }


        } else {
            //Otherwise, check to see if anything on the current child can be combined with 
            //anything in the chain thus far


            for (auto& child_cluster_head : child_clusters.read_cluster_heads) {
                //Go through all clusters of the current child and see if they can be combined with anything on the chain
                const size_t read_num = child_cluster_head.first;
                const size_t cluster_num = tree_state.read_union_find[read_num].find_group(child_cluster_head.second);
                pair<size_t, size_t> dists = tree_state.read_cluster_heads_to_distances[read_num][cluster_num];
                const size_t distance_left = child_is_reversed ? dists.second : dists.first;
                const size_t distance_right = child_is_reversed ? dists.first : dists.second;
                //Distance between this cluster and a cluster on the same read from the previous child
                size_t distance_between = SnarlDistanceIndex::minus(
                         SnarlDistanceIndex::sum({distance_left, 
                                                  distance_from_last_child_to_current_child,
                                                  old_best_right_by_read[read_num]}),
                          1);
                //Distance between this cluster and any cluster on the previous child
                size_t fragment_distance_between = SnarlDistanceIndex::minus(
                          SnarlDistanceIndex::sum({distance_left, 
                                                   distance_from_last_child_to_current_child,  
                                                   old_best_right}),
                          1);
                //The new distances from this child to the start of the chain and the end of this child
                pair<size_t, size_t> new_distances = make_pair(
                        SnarlDistanceIndex::sum({distance_left, distance_from_chain_start_to_current_node}),
                        SnarlDistanceIndex::sum({distance_right, distance_from_current_end_to_end_of_chain})); 

                if (distance_between <= tree_state.read_distance_limit) {
#ifdef DEBUG_CLUSTER
                    cerr << "\t\tCombine child cluster " << read_num << ":" << cluster_num << endl;
#endif
                    //If this cluster can be merged with anything on the chain
                    if (new_cluster_by_read[read_num].first.first == std::numeric_limits<size_t>::max()){
                        //If nothing is in the combined cluster yet, this is the new combined cluster
                        new_cluster_by_read[read_num] = make_pair(child_cluster_head, new_distances);
                    } else {
                        //Otherwise, remember to forget about the old cluster
                        //Union the two clusters
                        tree_state.read_union_find.at(read_num).union_groups(cluster_num, new_cluster_by_read[read_num].first.second);
                        //Find the best distances of the two
                        size_t new_best_left= std::min(new_cluster_by_read[read_num].second.first, new_distances.first);
                        size_t new_best_right= std::min(new_cluster_by_read[read_num].second.second, new_distances.second);
                        //And remember the new combined cluster head
                        new_cluster_by_read[read_num] = make_pair(
                                make_pair(read_num, tree_state.read_union_find.at(read_num).find_group(cluster_num)),
                               make_pair(new_best_left, new_best_right)); 
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
                        tree_state.fragment_union_find.union_groups(new_cluster_head_fragment, 
                                cluster_num + tree_state.read_index_offsets[read_num]);
                    }
                    new_cluster_head_fragment = cluster_num + tree_state.read_index_offsets[read_num];
                }


                //Update the best distances
                chain_clusters.fragment_best_left = std::min(chain_clusters.fragment_best_left, new_distances.first);
                chain_clusters.fragment_best_right = std::min(chain_clusters.fragment_best_right, new_distances.second); 
                chain_clusters.read_best_left[read_num] = std::min(chain_clusters.read_best_left[read_num], new_distances.first);
                chain_clusters.read_best_right[read_num] = std::min(chain_clusters.read_best_right[read_num], new_distances.second);
            }

            for (auto& chain_cluster_head : chain_clusters.read_cluster_heads) {
                //Each has distances up to the previous node
                const size_t read_num = chain_cluster_head.first;
                const size_t cluster_num = tree_state.read_union_find[read_num].find_group(chain_cluster_head.second);
                pair<size_t, size_t> dists = tree_state.read_cluster_heads_to_distances[read_num][cluster_num];

                //Best distance to the left side (relative to the chain) of the current child
                const size_t current_distance_left = child_is_reversed ? child_clusters.read_best_right[read_num]
                                                                              : child_clusters.read_best_left[read_num];
                const size_t current_fragment_distance_left = child_is_reversed ? child_clusters.fragment_best_right 
                                                                                       : child_clusters.fragment_best_left;
                const size_t distance_right = dists.second;
                const size_t distance_between = SnarlDistanceIndex::minus(
                        SnarlDistanceIndex::sum({distance_right, 
                                                 distance_from_last_child_to_current_child, 
                                                 current_distance_left}),
                        1);
                const size_t distance_between_fragment = SnarlDistanceIndex::minus(
                        SnarlDistanceIndex::sum({distance_right, 
                                                 distance_from_last_child_to_current_child, 
                                                 current_fragment_distance_left}),
                        1);
                const pair<size_t, size_t> new_distances = make_pair(
                    dists.first,
                    SnarlDistanceIndex::sum({dists.second, 
                                             distance_from_last_child_to_current_end,
                                             distance_from_current_end_to_end_of_chain}));
                     
                if (distance_between <= tree_state.read_distance_limit) {
#ifdef DEBUG_CLUSTER
                    cerr << "\t\tCombine chain cluster " << read_num << ":" << cluster_num << endl;
#endif
                    //If we can union the reads
                    if (new_cluster_by_read[read_num].first.first == std::numeric_limits<size_t>::max()) {
                        new_cluster_by_read[read_num] = make_pair(chain_cluster_head, new_distances);
                    } else {
                        //Otherwise, remember to forget about the old cluster
                        //Union the two clusters
                        tree_state.read_union_find.at(read_num).union_groups(cluster_num, new_cluster_by_read[read_num].first.second);
                        //Find the best distances of the two
                        size_t new_best_left= std::min(new_cluster_by_read[read_num].second.first, new_distances.first);
                        size_t new_best_right= std::min(new_cluster_by_read[read_num].second.second, new_distances.second);
                        //And remember the new combined cluster head
                        new_cluster_by_read[read_num] = make_pair(
                                make_pair(read_num, tree_state.read_union_find.at(read_num).find_group(cluster_num)),
                               make_pair(new_best_left, new_best_right)); 
                    }
                } else {
                    size_t distance_limit = tree_state.fragment_distance_limit == 0
                                             ? tree_state.read_distance_limit 
                                             : tree_state.fragment_distance_limit; 
                    if (new_distances.first <= distance_limit || new_distances.second <= distance_limit){
                        to_add.emplace_back(make_pair(read_num, cluster_num), new_distances);
                    }
                }
                if (tree_state.fragment_distance_limit != 0 && distance_between_fragment <= tree_state.fragment_distance_limit) {
                    //If we can union the fragments
                    if (new_cluster_head_fragment != std::numeric_limits<size_t>::max()) {
                        tree_state.fragment_union_find.union_groups(new_cluster_head_fragment, 
                                cluster_num + tree_state.read_index_offsets[read_num]);
                    }
                    new_cluster_head_fragment = cluster_num + tree_state.read_index_offsets[read_num];
                }
                chain_clusters.fragment_best_left = std::min(chain_clusters.fragment_best_left, new_distances.first);
                chain_clusters.fragment_best_right = std::min(chain_clusters.fragment_best_right, new_distances.second); 
                chain_clusters.read_best_left[read_num] = std::min(chain_clusters.read_best_left[read_num], new_distances.first);
                chain_clusters.read_best_right[read_num] = std::min(chain_clusters.read_best_right[read_num], new_distances.second);
            }

            //Remove clusters that got combined
            chain_clusters.read_cluster_heads.clear();
            //Add new clusters that weren't combined
            for (pair<pair<size_t, size_t>, pair<size_t, size_t>>& cluster : to_add) {
                chain_clusters.read_cluster_heads.emplace(cluster.first);
                tree_state.read_cluster_heads_to_distances[cluster.first.first][cluster.first.second]=cluster.second;
            }
            //Add new clusters that were combined
            for (pair<pair<size_t, size_t>, pair<size_t, size_t>>& cluster : new_cluster_by_read) {
                if (cluster.first.first != std::numeric_limits<size_t>::max()){
                    chain_clusters.read_cluster_heads.emplace(cluster.first);
                    tree_state.read_cluster_heads_to_distances[cluster.first.first][cluster.first.second]=cluster.second;
                }
            }
        }
        

        //Update the last node we saw to this one
        last_child = child_clusters;

#ifdef DEBUG_CLUSTER
    cerr << "\tintermediate clusters on " << distance_index.net_handle_as_string(chain_handle) << " after child " << distance_index.net_handle_as_string(child_handle) << endl;
    cerr << "\t   with best left and right values: " << chain_clusters.fragment_best_left << " "
         << chain_clusters.fragment_best_right << endl;
    bool got_left = false;
    bool got_right = false;
    for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
        cerr << "\t for read num " << read_num << " best left: " << chain_clusters.read_best_left[read_num] << " best right: " << chain_clusters.read_best_right[read_num] << endl;
        bool got_read_left=false;
        bool got_read_right = false;
        bool any_clusters = false;
        for (pair<size_t,size_t> c : chain_clusters.read_cluster_heads) {
            if (c.first == read_num) {
                any_clusters = true;
                pair<size_t, size_t> dists = tree_state.read_cluster_heads_to_distances[c.first][c.second];
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
        chain_clusters.read_cluster_heads.emplace(cluster_head.first);
        tree_state.read_cluster_heads_to_distances[cluster_head.first.first][cluster_head.first.second] = cluster_head.second;
        chain_clusters.fragment_best_left = std::min(chain_clusters.fragment_best_left, cluster_head.second.first);
        chain_clusters.read_best_left[cluster_head.first.first] = std::min(chain_clusters.read_best_left[cluster_head.first.first], cluster_head.second.first);
    }
    //If the chain loops, then we also have to compare the first thing we saw to the last things
    bool is_looping_chain = distance_index.get_record_offset(distance_index.get_bound(chain_handle, false, true)) == distance_index.get_record_offset(distance_index.get_bound(chain_handle, true, true)); 
    if (is_looping_chain){
        vector<size_t> combined_cluster_by_read (tree_state.all_seeds->size(), std::numeric_limits<size_t>::max()); 
        size_t combined_cluster_fragment = std::numeric_limits<size_t>::max();
        for (auto& cluster_head : chain_clusters.read_cluster_heads) {
            size_t read_num = cluster_head.first;
            size_t cluster_num=cluster_head.second;
            size_t dist_left = tree_state.read_cluster_heads_to_distances[read_num][cluster_num].first;
            size_t dist_right = tree_state.read_cluster_heads_to_distances[read_num][cluster_num].second;

            size_t distance_between_left_right =SnarlDistanceIndex::minus(SnarlDistanceIndex::sum({dist_left, chain_clusters.read_best_right[read_num]}), 1); 
            size_t distance_between_right_left = SnarlDistanceIndex::minus( SnarlDistanceIndex::sum({dist_right, chain_clusters.read_best_left[read_num]}), 1);
            if (distance_between_left_right <= tree_state.read_distance_limit
             || distance_between_right_left <= tree_state.read_distance_limit) {
                //If we can combine the read
                if (combined_cluster_by_read[read_num] == std::numeric_limits<size_t>::max()) {
                    combined_cluster_by_read[read_num] = cluster_num;
                } else {
                    tree_state.read_union_find[read_num].union_groups(combined_cluster_by_read[read_num], cluster_num);
                }
            }
            size_t distance_between_left_right_fragment = SnarlDistanceIndex::minus(SnarlDistanceIndex::sum({dist_left, chain_clusters.fragment_best_right}) ,1);
            size_t distance_between_right_left_fragment = SnarlDistanceIndex::minus(SnarlDistanceIndex::sum({dist_right, chain_clusters.fragment_best_left}) ,1);
            if (tree_state.fragment_distance_limit != 0 &&
                (distance_between_left_right_fragment <= tree_state.fragment_distance_limit
              || distance_between_right_left_fragment <= tree_state.fragment_distance_limit)) {

                if (combined_cluster_fragment != std::numeric_limits<size_t>::max()) {
                    tree_state.fragment_union_find.union_groups(combined_cluster_fragment, 
                            cluster_num + tree_state.read_index_offsets[read_num]);
                }
                combined_cluster_fragment = cluster_num + tree_state.read_index_offsets[read_num];
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
        cerr << "\t for read num " << read_num << " best left: " << chain_clusters.read_best_left[read_num] << " best right: " << chain_clusters.read_best_right[read_num] << endl;
        bool got_read_left=false;
        bool got_read_right = false;
        bool any_clusters = false;
        for (pair<size_t,size_t> c : chain_clusters.read_cluster_heads) {
            if (c.first == read_num) {
                any_clusters = true;
                pair<size_t, size_t> dists = tree_state.read_cluster_heads_to_distances[c.first][c.second];
                cerr << "\t\t" << c.first << ":"<<c.second << ": left: " << dists.first << " right : " << dists.second << ": ";
                bool has_seeds = false;
                for (size_t x = 0 ; x < tree_state.all_seeds->at(c.first)->size() ; x++) {
                    if (tree_state.read_union_find[c.first].find_group(x) == c.second) {
                        cerr << tree_state.all_seeds->at(c.first)->at(x).pos << " ";
                        has_seeds = true;
                    }
                }
                //assert(dists.first == std::numeric_limits<size_t>::max() || dists.first >= chain_clusters.read_best_left[read_num]);
                //assert(dists.second == std::numeric_limits<size_t>::max() || dists.second >= chain_clusters.read_best_right[read_num]);
                //assert(dists.first == std::numeric_limits<size_t>::max() || dists.first >= chain_clusters.fragment_best_left);
                //assert(dists.second == std::numeric_limits<size_t>::max() || dists.second >= chain_clusters.fragment_best_right);
                if (dists.first == chain_clusters.fragment_best_left) {got_left = true;}
                if (dists.second == chain_clusters.fragment_best_right) {got_right = true;}
                if (dists.first == chain_clusters.read_best_left[read_num]) {got_read_left = true;}
                if (dists.second == chain_clusters.read_best_right[read_num]) {got_read_right = true;}
                cerr << endl;
                //assert(has_seeds);
            }
        }
        //assert(!any_clusters ||got_read_left ||  chain_clusters.read_best_left[read_num] == std::numeric_limits<size_t>::max());
        //assert(!any_clusters ||got_read_right ||  chain_clusters.read_best_right[read_num] == std::numeric_limits<size_t>::max());
    }
    //assert(got_left);
    //assert(got_right);

    for (pair<size_t, size_t> group_id : chain_clusters.read_cluster_heads) {
        //assert (group_id.first.second == tree_state.read_union_find[group_id.first.first].find_group(group_id.first.second));
    }
#endif
}

//Cluster the root
//all children of the root will be in tree_state.root_children
//This is basically cluster_one_snarl except the snarl is the root, which has no boundary nodes
void NewSnarlSeedClusterer::cluster_root(TreeState& tree_state) const { 
#ifdef DEBUG_CLUSTER
    cerr << "Finding clusters on the root with " << tree_state.root_children.size() << " children" << endl;
#endif

    //Keep track of all clusters on the root
    NodeClusters root_clusters(distance_index.get_root(), tree_state.all_seeds->size(), distance_index);

    vector<vector<pair<size_t, size_t>>> child_distances (tree_state.all_seeds->size());
    for (size_t i = 0 ; i < tree_state.all_seeds->size() ; i++) {
        child_distances[i] = std::vector<pair<size_t, size_t>>(tree_state.all_seeds->at(i)->size(),
            make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
    }

    for (auto& kv : tree_state.root_children) {
        vector<size_t>& children = kv.second;

        for (size_t i = 0; i < children.size() ; i++) {
            //Go through each child node of the netgraph

            NodeClusters& child_clusters_i = tree_state.all_node_clusters[children[i]];
            for (const pair<size_t, size_t>& head : child_clusters_i.read_cluster_heads) {
                child_distances[head.first][head.second] = tree_state.read_cluster_heads_to_distances[head.first][head.second];
            }

            for (size_t j = 0 ; j <= i ; j++){
                //Go through other child net graph nodes up to and including i

                //Get the other node and its clusters
                NodeClusters& child_clusters_j = tree_state.all_node_clusters[children[j]];



                compare_and_combine_cluster_on_child_structures(tree_state, child_clusters_i,
                            child_clusters_j, root_clusters, child_distances, true);

            }
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

size_t NewSnarlSeedClusterer::add_child_to_vector(TreeState& tree_state, hash_map<net_handle_t, pair<size_t, vector<size_t>>>* parent_to_child_map, const net_handle_t& parent,size_t child_index) const {

    //THe indexes of the parent in all_node_clusters, to be returned
    size_t parent_index = std::numeric_limits<size_t>::max();
    if (parent_to_child_map->count(parent) == 0) {
        //If we haven't seen the parent before

        //Make a new NodeClusters for the parent
        tree_state.all_node_clusters.emplace_back(parent, tree_state.all_seeds->size(), distance_index);

        //Add the parent to the map
        parent_to_child_map->emplace(parent, 
               make_pair(tree_state.all_node_clusters.size()-1, vector<size_t>(0)));
        parent_index = tree_state.all_node_clusters.size()-1;
    } 

    if (distance_index.is_snarl(parent)){
        //If the parent is a snarl, then the order doesn't matter
        parent_to_child_map->at(parent).second.emplace_back(child_index);
    } else {
        assert(distance_index.is_chain(parent));
        //Helper function to insert a NodeClusters into a sorted vector<NodeClusters> where the vector contains only 
        //children of the same chain. The items will be sorted by their position in the chain
        auto insert_in_order = [&](vector<size_t>& input_vector, size_t item) {
        
            //The rank of the item we're inserting is its rank in its parent
            //The rank is really the offset of the record in it's parent, so it doesn't matter what 
            //the value is except that it will be increasing as we walk along the chain
            //Get an iterator to where the thing we're inserting should go in the vector 
            std::vector<size_t>::iterator insert_itr = std::upper_bound(input_vector.begin(), input_vector.end(),
               0, [&](size_t val, size_t vector_item) {
                   //This should return true if the vector_item goes after the item with rank val
                   //based on the "rank" in the cached net handle, which is th e  offset in the chain

                   return distance_index.is_ordered_in_chain(tree_state.all_node_clusters[item].containing_net_handle,
                        tree_state.all_node_clusters[vector_item].containing_net_handle);
               });
            //Add the item to the vector in sorted order
            input_vector.insert(insert_itr, item);
        };
        insert_in_order(parent_to_child_map->at(parent).second, child_index);
    }
    return parent_index;
}

}
