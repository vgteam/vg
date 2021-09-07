#include "snarl_seed_clusterer.hpp"

#include <algorithm>

#define DEBUG_CLUSTER
namespace vg {

    NewSnarlSeedClusterer::NewSnarlSeedClusterer( SnarlDistanceIndex& distance_index) :
                                            distance_index(distance_index){
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

        //For each level of the snarl tree, maps chains (net_handle_t of the snarl) at that level to children 
        //(net_handle_t of the child) belonging to the chain. Starts out with nodes from get_nodes, 
        //then adds snarls as we walk up the snarl tree
        //This is later used to populate snarl_to_node in the tree state
        //The outer vector has one entry for each read
        //The outer hash_map maps a net_handle_t of a chain to a vector of pairs, where each pair is the child of the 
        vector<hash_map<net_handle_t, vector<NodeClusters>>> chain_to_children_by_level;
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

        //Initialize the tree state to be the bottom level
        tree_state.chain_to_children = std::move(chain_to_children_by_level[chain_to_children_by_level.size() - 1]);

        for (int depth = chain_to_children_by_level.size() - 1 ; depth >= 0 ; depth --) {
            //Go through each level of the tree, bottom up, and cluster that
            // level. Each level includes the snarl at that level, the nodes
            // belonging to those snarls, and the chains comprised of them
            //
            // tree_state knows all children of the snarls at this level

            // Bring in the direct child nodes that come in at this level in the snarl tree.
            // They only ever occur below the root.
            if (depth != 0) {
                tree_state.parent_chain_to_children = std::move(chain_to_children_by_level[depth - 1]);
            }

#ifdef DEBUG_CLUSTER
assert(tree_state.read_index_offsets[0] == 0);
for (size_t i = 1 ; i < tree_state.all_seeds->size() ; i++) {
    assert (tree_state.read_index_offsets[i] + tree_state.all_seeds->at(i)->size() == tree_state.read_index_offsets[i+1]);
}
#endif
            //Cluster all the chains at this depth
            //Also records which chains are in snarls and the parents of these
            //chains in tree_state.parent_chain_to_children
            cluster_chain_level(tree_state);

            //And cluster all the snarls, record the parents of these snarls
            cluster_snarl_level(tree_state);


            // Swap buffer over for the next level
            tree_state.chain_to_children = std::move(tree_state.parent_chain_to_children);
            tree_state.snarl_to_children.clear();
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
              vector<hash_map<net_handle_t,vector<NodeClusters>>>& chain_to_children_by_level) const {
#ifdef DEBUG_CLUSTER
cerr << "Add all seeds to nodes: " << endl << "\t";
#endif

        //Helper function to insert a NodeClusters into a sorted vector<NodeClusters> where the vector contains only 
        //children of the same chain. The items will be sorted by their position in the chain
        auto insert_in_order = [&](vector<NodeClusters>& input_vector, NodeClusters item) {
        
            //The rank of the item we're inserting is its rank in its parent
            //The rank is really the offset of the record in it's parent, so it doesn't matter what 
            //the value is except that it will be increasing as we walk along the chain
            //Get an iterator to where the thing we're inserting should go in the vector 
            std::vector<NodeClusters>::iterator insert_itr = std::upper_bound(input_vector.begin(), input_vector.end(),
               0, [&](size_t val, NodeClusters vector_item) {
                   //This should return true if the vector_item goes after the item with rank val
                   return distance_index.is_ordered_in_chain(item.containing_net_handle, vector_item.containing_net_handle);
               });
            //Add the item to the vector in sorted order
            input_vector.insert(insert_itr, std::move(item));
        };

        // Assign each seed to a node.
        hash_set<id_t> seen_nodes;
        for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++){ 
            const vector<Seed>* seeds = tree_state.all_seeds->at(read_num);
            vector<pair<id_t, size_t>>& node_to_seeds = tree_state.node_to_seeds.at(read_num);
            for (size_t i = 0; i < seeds->size(); i++) {
                pos_t pos = seeds->at(i).pos;
                id_t id = get_id(pos);
                

                //Assign the seed to a node
//                if (!seeds->at(i).is_top_level_node && !seeds->at(i).is_top_level_snarl) {
                    //If this seed is not on a top-level chain or top-level simple bubble
                    //A seed can still be added here if it is on a top-level chain
                node_to_seeds.emplace_back(id, i);
#ifdef DEBUG_CLUSTER
                cerr << read_num << ":" << pos << ", ";
#endif

                 //And the node to a chain
                if (seen_nodes.count(id) < 1) {
                     seen_nodes.insert(id);
                     net_handle_t parent = distance_index.get_parent(distance_index.get_node_net_handle(id));
                     size_t depth = distance_index.get_depth(parent);
                     if (depth+1 > chain_to_children_by_level.size()) {
                         chain_to_children_by_level.resize(depth+1);
                     }
                     //Create a new NodeClusters for this node and add it to it's parent chain
                     if (chain_to_children_by_level[depth].count(parent) == 0) {
                         vector<NodeClusters> empty_vector;
                         chain_to_children_by_level[depth].emplace(parent, std::move(empty_vector));
                     }
                     insert_in_order(chain_to_children_by_level[depth][parent],
                                     NodeClusters(distance_index.get_node_net_handle(id), tree_state.all_seeds->size()));
                 } else {
                 }
//                }
                /* TODO: This uses cached distance index information, which I might put back but hopefully it won't be necessary anymore
                 * else if (seeds->at(i).is_top_level_node) {
                    //If this is a top-level seed, defer clustering until we reach the top-level
                    size_t component = seeds->at(i).component;
                    if (tree_state.component_to_index.count(component) == 0) {
                        tree_state.component_to_index.emplace(component, tree_state.top_level_seed_clusters.size());
                        tree_state.top_level_seed_clusters.emplace_back();
                        tree_state.simple_snarl_to_nodes_by_component.emplace_back();
                    }
                    tree_state.top_level_seed_clusters[tree_state.component_to_index.at(component)].emplace_back(read_num, i);

                } else if (seeds->at(i).is_top_level_snarl) {
                    //This is a top-level simple snarl
#ifdef DEBUG_CLUSTER
                    cerr << "(bubble)" << read_num << ":" << pos << " " << seeds->at(i).start_length << " " << seeds->at(i).end_length << " " << seeds->at(i).snarl_rank << " " << seeds->at(i).rev_in_chain << ", ";
#endif

                    node_to_seeds.emplace_back(id, i);

                    size_t component = distance_index.get_connected_component(id);
                    if (tree_state.component_to_index.count(component) == 0) {
                        tree_state.component_to_index.emplace(component, tree_state.top_level_seed_clusters.size());
                        tree_state.top_level_seed_clusters.emplace_back();
                        tree_state.simple_snarl_to_nodes_by_component.emplace_back();
                    }
                    hash_map<size_t, vector<tuple<id_t, bool, size_t, size_t, size_t>>>& simple_snarl_to_nodes = 
                        tree_state.simple_snarl_to_nodes_by_component[tree_state.component_to_index.at(component)];
                    if (seen_nodes.count(id) < 1) {
                        seen_nodes.insert(id);
                        simple_snarl_to_nodes[seeds->at(i).snarl_rank].emplace_back(
                            id, seeds->at(i).rev_in_chain, seeds->at(i).start_length, seeds->at(i).end_length, seeds->at(i).node_length);
                    }
                }
                */
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
    void NewSnarlSeedClusterer::cluster_snarl_level(TreeState& tree_state) const {

        for (auto& kv : tree_state.snarl_to_children){
            //Go through each of the snarls at this level, cluster them,
            //and find which chains they belong to, if any
            //key is the index of the snarl and value is a vector of pair of
            // net_handle_t, NodeClusters

            net_handle_t snarl_handle = kv.first;

#ifdef DEBUG_CLUSTER
            cerr << "Cluster one snarl " << distance_index.net_handle_as_string(snarl_handle) << endl;
#endif
            if (!distance_index.is_root(distance_index.get_parent(snarl_handle))){
                //If this snarl is in a chain, cluster and add let the
                //tree state know which chain it belongs to (in order)

                add_child_to_vector(tree_state.parent_chain_to_children, distance_index.get_parent(snarl_handle), cluster_one_snarl(tree_state, snarl_handle));

#ifdef DEBUG_CLUSTER
                cerr << "\tRecording snarl " << distance_index.net_handle_as_string(snarl_handle)  << " as a child of "
                      << distance_index.net_handle_as_string(distance_index.get_parent(snarl_handle)) << endl;
#endif

            } else {
                //TODO: I think I don't need to check if the parent is the root because it should always be a chain
                assert(false);
            }
        }
    }

    void NewSnarlSeedClusterer::cluster_chain_level(TreeState& tree_state) const {
        for (auto& kv : tree_state.chain_to_children) {
            //For each chain at this level that has relevant child snarls in it,
            //find the clusters.


            // Get the chain's number
            net_handle_t chain_handle = kv.first;

#ifdef DEBUG_CLUSTER
            cerr << "Cluster one chain " <<  distance_index.net_handle_as_string(chain_handle) << endl;
#endif

            // Compute the clusters for the chain
            if (distance_index.is_root(distance_index.get_parent(chain_handle))) {
                tree_state.root_children.emplace_back(cluster_one_chain(tree_state, chain_handle));
            } else {
                add_child_to_vector(tree_state.snarl_to_children, distance_index.get_parent(chain_handle), cluster_one_chain(tree_state, chain_handle));
#ifdef DEBUG_CLUSTER
                cerr << "\tRecording " << distance_index.net_handle_as_string(chain_handle)
                    << " as a child of " << distance_index.net_handle_as_string(distance_index.get_parent(chain_handle)) << endl;
#endif
            }
        }
    }


    void NewSnarlSeedClusterer::cluster_one_node(
                       TreeState& tree_state, NodeClusters& node_clusters) const {
#ifdef DEBUG_CLUSTER
        cerr << "Finding clusters on node " << distance_index.net_handle_as_string(node_clusters.containing_net_handle) << endl;
#endif

        size_t node_length = distance_index.node_length(node_clusters.containing_net_handle);
        nid_t node_id = distance_index.node_id(node_clusters.containing_net_handle);


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
                    node_clusters.read_cluster_heads[make_pair(read_num, group_id)] =  
                        make_pair(node_clusters.read_best_left[read_num],node_clusters.read_best_right[read_num]);

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
                for (pair<pair<size_t, size_t>, pair<size_t,size_t>> c : node_clusters.read_cluster_heads) {
                    if (c.first.first == read_num) {
                        pair<size_t, size_t> dists = c.second;
                        cerr << "\t\t" << c.first.first << ":"<<c.first.second << ": left: " << c.second.first << " right : " << c.second.second << ": ";
                        bool has_seeds = false;
                        for (size_t x = 0 ; x < tree_state.all_seeds->at(c.first.first)->size() ; x++) {
                            if (tree_state.read_union_find[c.first.first].find_group(x) == c.first.second) {
                                cerr << tree_state.all_seeds->at(c.first.first)->at(x).pos << " ";
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
                assert(got_read_left || node_clusters.read_best_left[read_num] == std::numeric_limits<size_t>::max());
                assert(got_read_right || node_clusters.read_best_right[read_num] == std::numeric_limits<size_t>::max());
            }
            assert(got_left);
            assert(got_right);

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
                    fragment_last_offset = offset;
                }
            } else {
                //This becomes a new read cluster
                if (read_last_cluster[read_num] != std::numeric_limits<size_t>::max()) {
                    //Record the previous cluster
                    node_clusters.read_cluster_heads[make_pair(read_num, read_last_cluster[read_num])] = 
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
                    fragment_last_offset = offset;
                }
            }
        }
        for (size_t i = 0 ; i < read_last_cluster.size() ; i++) {
            if (read_last_cluster[i] != std::numeric_limits<size_t>::max()) {
                node_clusters.read_cluster_heads[make_pair(i, read_last_cluster[i])] =
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
            for (pair<pair<size_t, size_t>, pair<size_t,size_t>> c : node_clusters.read_cluster_heads) {
                if (c.first.first == read_num) {
                    pair<size_t, size_t> dists = c.second;
                    cerr << "\t\t" << c.first.first << ":"<<c.first.second << ": left: " << dists.first << " right : " << dists.second << ": ";
                    bool has_seeds = false;
                    for (size_t x = 0 ; x < tree_state.all_seeds->at(c.first.first)->size() ; x++) {
                        if (tree_state.read_union_find[c.first.first].find_group(x) == c.first.second) {
                            cerr << tree_state.all_seeds->at(c.first.first)->at(x).pos << " ";
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
            assert(got_read_left || node_clusters.read_best_left[read_num] == std::numeric_limits<size_t>::max());
            assert(got_read_right || node_clusters.read_best_right[read_num] == std::numeric_limits<size_t>::max());
        }
        assert(got_left);
        assert(got_right);
        for (pair<pair<size_t, size_t>, pair<size_t, size_t>> group_id : node_clusters.read_cluster_heads) {
            assert (group_id.first.second == tree_state.read_union_find[group_id.first.first].find_group(group_id.first.second));
        }
#endif
        return;
        
    };


/*
    //Given a structure and its parent, try to combine the clusters on the children
    //TODO: I think this is correct:
    //This doesn't update the distances to the ends of the child structure, since any updates would be found using paths in
    //the parent structure, so these distances would be taken into account later anyway
    //TODO: I think I might not need a separate function for things on the same structure, I could just do the one on different child structures but give it the same structure twice
    void NewSnarlSeedClusterer::compare_and_combine_cluster_on_same_structure(TreeState& tree_state, NodeClusters& child_clusters, 
        NodeClusters& parent_clusters) {
        net_handle_t& child_handle = child_clusters.containin_net_handle;
        net_handle_t& parent_handle = parent_clusters.containing_net_handle;

        //Get the distances between the two sides of the child in the parent
        size_t distance_left_left = distance_index.distance_in_parent(parent_handle, distance_index.flip(child_handle), distance_index.flip(child_handle));
        size_t distance_right_right = distance_index.distance_in_parent(parent_handle, child_handle, child_handle);
        size_t distance_left_right = distance_index.distance_in_parent(parent_handle, distance_index.flip(child_handle), child_handle);

        //The cluster heads that will be removed from read_cluster_heads
        vector<pair<size_t, size_t>> to_erase;

        //These will be the cluster heads and distances of everything combined by taking the indicated path
        //one cluster head per read
        //pair< pair<read_num, seed_num>, pair<left dist, right_dist>>
        //The default value will be ((inf, 0), (0,0)). Only the inf gets checked to see if it's a real value so I filled it in with 0 so I wouldn't have to type out inf 
        vector<pair<pair<size_t, size_t>, pair<size_t, size_t>> new_cluster_head_left_left_by_read(tree_state.all-seeds->size(), 
                make_pair(make_pair(std::numeric_limits<size_t>::max(), 0), make_pair(0, 0)));
        vector<pair<pair<size_t, size_t>, pair<size_t, size_t>> new_cluster_head_right_right_by_read(tree_state.all-seeds->size(),
            make_pair(make_pair(std::numeric_limits<size_t>::max(), 0), make_pair(0, 0)));
        vector<pair<pair<size_t, size_t>, pair<size_t, size_t>> new_cluster_head_left_right_by_read(tree_state.all-seeds->size(),
            make_pair(make_pair(std::numeric_limits<size_t>::max(), 0), make_pair(0, 0)));

        //And the new cluster heads for the fragment
        //These are the values of the cluster heads in the union finds, which include the values from read_index_offset
        size_t new_cluster_head_left_left_fragment;
        size_t new_cluster_head_right_right_fragment;
        size_t new_cluster_head_left_right_fragment;

        //Go through the clusters once, and see which ones can be put in each of the four possible combined clusters
        //If a cluster can be combined with another cluster by taking the right-right loop, then any cluster 
        //with a right distance that is less than the distance limit - the best right distance for the child will
        //be combined with the cluster with the best right distances. So in the end there will only be three combined
        //clusters, one for each path
        //This strategy ends up unioning a cluster with itself but it only goes through the clusters once so I think
        //it's efficient
        for (auto& cluster_head_and_distance : child_clusters.read_cluster_heads) {
            size_t read_num = cluster_head_and_distance.first.first;
            size_t cluster_num = cluster_head_and_distance.first.second;
            pair<size_t, size_t>& distances = cluster_head_and_distance.second;

            //Try to combine going left-left

            //TODO: Use sum to sum infinities better
            if (distances.first + distance_left_left + child_clusters.read_best_left.at(read_num) <= tree_state.read_distance_limit) {
                //If this can be combined with the new left-left read cluster
                if (new_cluster_head_left_left_by_read.at(read_num) == std::numeric_limits<size_t>::max()){
                    //new cluster head
                    new_cluster_head_left_left_by_read.at(read_num) = make_pair(read_num, cluster_num);
                } else {
                    //Combine with old cluster head
                    tree_state.read_union_find.at(read_num).union_groups(cluster_num, new_cluster_head_left_left.at(read_num));
                    new_cluster_head_left_left_by_read.at(read_num) = make_pair(read_num, tree_state.read_union_find.at(read_num).find_group(cluster_num));
                }
                to_erase.emplace_back(read_num, cluster_num);
                tree_state.fragment_union_find.union_groups(cluster_num+tree_state.read_index_offsets[read_num], 
                                                                    new_cluster_head_left_left_fragment);
                new_cluster_head_left_left_fragment = tree_state.fragment_union_find.find_group(new_cluster_head_left_left_fragment);

            } else if (tree_state.fragment_distance_limit != 0 && 
                        distances.first + distance_left_left + child_clusters.fragment_best_left <= tree_state.fragment_distance_limit ) {
                //Just union the fragment
                tree_state.fragment_union_find.union_groups(cluster_num+tree_state.read_index_offsets[read_num], 
                                                                    new_cluster_head_left_left_fragment);
                new_cluster_head_left_left_fragment = tree_state.fragment_union_find.find_group(new_cluster_head_left_left_fragment);
            }

            //Try to combine going right-right

            //TODO: Use sum to sum infinities better
            if (distances.second + distance_right_right + child_clusters.read_best_right.at(read_num) <= tree_state.read_distance_limit) {
                //If this can be combined with the new right-right read cluster
                if (new_cluster_head_right_right_by_read.at(read_num) == std::numeric_limits<size_t>::max()){
                    //new cluster head
                    new_cluster_head_right_right_by_read.at(read_num) = make_pair(read_num, cluster_num);
                } else {
                    //Combine with old cluster head
                    tree_state.read_union_find.at(read_num).union_groups(cluster_num, new_cluster_head_right_right.at(read_num));
                    new_cluster_head_right_right_by_read.at(read_num) = make_pair(read_num, tree_state.read_union_find.at(read_num).find_group(cluster_num));
                }
                to_erase.emplace_back(read_num, cluster_num);
                tree_state.fragment_union_find.union_groups(cluster_num+tree_state.read_index_offsets[read_num], 
                                                                    new_cluster_head_right_right_fragment);
                new_cluster_head_right_right_fragment = tree_state.fragment_union_find.find_group(new_cluster_head_right_right_fragment);

            } else if (tree_state.fragment_distance_limit != 0 && 
                        distances.second + distance_right_right + child_clusters.fragment_best_right <= tree_state.fragment_distance_limit ) {
                //Just union the fragment
                tree_state.fragment_union_find.union_groups(cluster_num+tree_state.read_index_offsets[read_num], 
                                                                    new_cluster_head_right_right_fragment);
                new_cluster_head_right_right_fragment = tree_state.fragment_union_find.find_group(new_cluster_head_right_right_fragment);
            }
            //Try to combine going left-right

            //TODO: Use sum to sum infinities better
            if (distances.first + distance_left_right + child_clusters.read_best_right.at(read_num) <= tree_state.read_distance_limit ||
                distances.second + distance_left_right + child_clusters.read_best_left.at(read_num) <= tree_state.read_distance_limit ) {
                //If this can be combined with the new left-right read cluster
                if (new_cluster_head_left_right_by_read.at(read_num) == std::numeric_limits<size_t>::max()){
                    //new cluster head
                    new_cluster_head_left_right_by_read.at(read_num) = make_pair(read_num, cluster_num);
                } else {
                    //Combine with old cluster head
                    tree_state.read_union_find.at(read_num).union_groups(cluster_num, new_cluster_head_left_right.at(read_num));
                    new_cluster_head_left_right_by_read.at(read_num) = make_pair(read_num, tree_state.read_union_find.at(read_num).find_group(cluster_num));
                }
                to_erase.emplace_back(read_num, cluster_num);
                tree_state.fragment_union_find.union_groups(cluster_num+tree_state.read_index_offsets[read_num], 
                                                                    new_cluster_head_left_right_fragment);
                new_cluster_head_left_right_fragment = tree_state.fragment_union_find.find_group(new_cluster_head_left_right_fragment);

            } else if (tree_state.fragment_distance_limit != 0 && 
                        distances.first + distance_left_right + child_clusters.fragment_best_right <= tree_state.fragment_distance_limit ||
                        distances.second + distance_left_right + child_clusters.fragment_best_left <= tree_state.fragment_distance_limit) {
                //Just union the fragment
                tree_state.fragment_union_find.union_groups(cluster_num+tree_state.read_index_offsets[read_num], 
                                                                    new_cluster_head_left_right_fragment);
                new_cluster_head_left_right_fragment = tree_state.fragment_union_find.find_group(new_cluster_head_left_right_fragment);
            }
        }

        //Now remove cluster heads that got combined with new ones
        for (pair<size_t, size_t>& cluster_head : to_erase) {
            child_clusters.read_cluster_heads.erase(cluster_head);
        }
        //And add back in the new cluster heads
        for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
            if (new_cluster_head_left_left_by_read.at(read_num).first != std::numeric_limits<size_t>::max()) {
                read_cluster_heads.insert(new_cluster_head_left_left_by_read.at(read_num));
            }
            if (new_cluster_head_right_right_by_read.at(read_num).first != std::numeric_limits<size_t>::max()) {
                read_cluster_heads.insert(new_cluster_head_right_right_by_read.at(read_num));
            }
            if (new_cluster_head_left_right_by_read.at(read_num).first != std::numeric_limits<size_t>::max()) {
                read_cluster_heads.insert(new_cluster_head_left_right_by_read.at(read_num));
            }
        }

        


    }
    */


    //Go through pairs of clusters of the two children and see which ones can be combined
    //The first child may not have been seen before, so all of it's clusters may be added to the parent, then
    //anything that was combined gets removed and only the cluster heads get added.
    //For the second child, everything is already in the parent so remove ones that were combined then
    //add the head of the combined clusters
    //
    //
    //This does basically the same thing as compare_and_combine_cluster_on_same_structure, except that
    //it updates the parent clusters
    //
    //Be sure to update the distances to the ends of the parent for the second child's clusters
    //TODO: Make sure to add the first child's clusters to the parent before looking at pairs and calling this
    void NewSnarlSeedClusterer::compare_and_combine_cluster_on_child_structures(TreeState& tree_state, NodeClusters& child_clusters1, 
        NodeClusters& child_clusters2, NodeClusters& parent_clusters, bool is_root) const {
#ifdef DEBUG_CLUSTER
        cerr << "\tCompare " << distance_index.net_handle_as_string(child_clusters1.containing_net_handle) 
             << " and " << distance_index.net_handle_as_string(child_clusters2.containing_net_handle)
             << " which are children of " << distance_index.net_handle_as_string(parent_clusters.containing_net_handle) << endl;
#endif

        net_handle_t& parent_handle = parent_clusters.containing_net_handle;
        net_handle_t& child_handle1 = child_clusters1.containing_net_handle;
        net_handle_t& child_handle2 = child_clusters2.containing_net_handle;


        //Get the distances between the two sides of the children in the parent
        size_t distance_left_left = distance_index.distance_in_parent(parent_handle, distance_index.flip(child_handle1), distance_index.flip(child_handle2));
        size_t distance_left_right = distance_index.distance_in_parent(parent_handle, distance_index.flip(child_handle1), child_handle2);
        size_t distance_right_right = distance_index.distance_in_parent(parent_handle, child_handle1, child_handle2);
        size_t distance_right_left = distance_index.distance_in_parent(parent_handle, child_handle1, distance_index.flip(child_handle2));

        /* Find the distances from the start/end of the parent to the left/right of the first child
         * If this is the root, then the distances are infinite since it has no start/end
         */
        if (is_root){

            child_clusters1.distance_start_left = std::numeric_limits<size_t>::max();
            child_clusters1.distance_start_right = std::numeric_limits<size_t>::max();
            child_clusters1.distance_end_left = std::numeric_limits<size_t>::max();
            child_clusters1.distance_end_right = std::numeric_limits<size_t>::max();

        } else {
            net_handle_t parent_start_in = distance_index.get_bound(parent_handle, false, true);
            net_handle_t parent_end_in = distance_index.get_bound(parent_handle, true, true);
            size_t start_length = !distance_index.is_chain(parent_handle) ? 0 
                        : distance_index.node_length(distance_index.get_bound(parent_handle, false, false));
                size_t end_length =  !distance_index.is_chain(parent_handle) ? 0 
                        : distance_index.node_length(distance_index.get_bound(parent_handle, true, false));

            child_clusters1.distance_start_left = distance_index.flip(child_handle1) == parent_start_in ? 0 :
                SnarlDistanceIndex::sum({start_length, distance_index.distance_in_parent(parent_handle, parent_start_in, distance_index.flip(child_handle1))});
            child_clusters1.distance_start_right = child_handle1 == parent_start_in ? 0 :
                 SnarlDistanceIndex::sum({start_length, distance_index.distance_in_parent(parent_handle, parent_start_in, child_handle1)});
            child_clusters1.distance_end_left = distance_index.flip(child_handle1) == parent_end_in ? 0 :
                 SnarlDistanceIndex::sum({end_length, distance_index.distance_in_parent(parent_handle, parent_end_in, distance_index.flip(child_handle1))});
            child_clusters1.distance_end_right = child_handle1 == parent_end_in ? 0  :
                SnarlDistanceIndex::sum({end_length, distance_index.distance_in_parent(parent_handle, parent_end_in, child_handle1)});

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
                size_t distance_between_fragments, pair<size_t, size_t> old_distances, pair<pair<size_t, size_t>, 
                pair<size_t, size_t>>& new_cluster_head_and_distances, size_t& new_cluster_head_fragment){
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
        for (auto& cluster_head_and_distance : child_clusters1.read_cluster_heads) {

            bool combined = false;
            size_t read_num = cluster_head_and_distance.first.first;
            size_t cluster_num = tree_state.read_union_find[read_num].find_group(cluster_head_and_distance.first.second);
            pair<size_t, size_t> distances = cluster_head_and_distance.second;

            //TODO: This gets done many times, maybe we could save it somewhere - in the NodeClusters as an optional field maybe
            size_t new_dist_left = std::min(SnarlDistanceIndex::sum({distances.first,child_clusters1.distance_start_left}),
                                            SnarlDistanceIndex::sum({distances.second,child_clusters1.distance_start_right}));
            size_t new_dist_right= std::min(SnarlDistanceIndex::sum({distances.first,child_clusters1.distance_end_left}),
                                            SnarlDistanceIndex::sum({distances.second,child_clusters1.distance_end_right}));
            pair<size_t, size_t> distances_to_parent = make_pair(new_dist_left, new_dist_right);
            if (parent_clusters.read_cluster_heads.count(make_pair(read_num, cluster_num)) > 0) {
                distances_to_parent = make_pair(
                    std::min(new_dist_left, parent_clusters.read_cluster_heads[make_pair(read_num, cluster_num)].first),
                    std::min(new_dist_right, parent_clusters.read_cluster_heads[make_pair(read_num, cluster_num)].second));
            }
                

            //Check if the left of 1 can connect with the left of 2
//#ifdef DEBUG_CLUSTER
//            cerr << "\t\tCheck any clusters in 1 going left to left of 2" << endl;
//#endif
            combined = combined | compare_and_combine_clusters (read_num, cluster_num, 
                SnarlDistanceIndex::sum({distances.first,distance_left_left,child_clusters2.read_best_left[read_num]}), 
                SnarlDistanceIndex::sum({distances.first,distance_left_left,child_clusters2.fragment_best_left}), 
                distances_to_parent,
                new_cluster_left_left_by_read[read_num], new_cluster_left_left_fragment);
            //Check if the left of 1 can connect with the right of 2
//#ifdef DEBUG_CLUSTER
//            cerr << "\t\tCheck any clusters in 1 going left to right of 2" << endl;
//#endif
            combined = combined | compare_and_combine_clusters (read_num, cluster_num, 
                SnarlDistanceIndex::sum({distances.first,distance_left_right,child_clusters2.read_best_right[read_num]}), 
                SnarlDistanceIndex::sum({distances.first,distance_left_right,child_clusters2.fragment_best_right}), 
                distances_to_parent, new_cluster_left_right_by_read[read_num], new_cluster_left_right_fragment);
            //Check if the right of 1 can connect with the right of 2
//#ifdef DEBUG_CLUSTER
//            cerr << "\t\tCheck any clusters in 1 going right to right of 2" << endl;
//#endif
            combined = combined | compare_and_combine_clusters (read_num, cluster_num, 
                SnarlDistanceIndex::sum({distances.second,distance_right_right,child_clusters2.read_best_right[read_num]}), 
                SnarlDistanceIndex::sum({distances.second,distance_right_right,child_clusters2.fragment_best_right}), 
                distances_to_parent,new_cluster_right_right_by_read[read_num], new_cluster_right_right_fragment);
            //Check if the right of 1 can connect with the left of 2
//#ifdef DEBUG_CLUSTER
//            cerr << "\t\tCheck any clusters in 1 going right to left of 2" << endl;
//#endif
            combined = combined | compare_and_combine_clusters (read_num, cluster_num, 
                SnarlDistanceIndex::sum({distances.second,distance_right_left,child_clusters2.read_best_left[read_num]}), 
                SnarlDistanceIndex::sum({distances.second,distance_right_left,child_clusters2.fragment_best_left}), 
                distances_to_parent, new_cluster_right_left_by_read[read_num], new_cluster_right_left_fragment);

            //If this cluster wasn't combined and hasn't been seen before, add it to the parent
            size_t cluster_head = tree_state.read_union_find[read_num].find_group(cluster_num); 
            if (!combined && parent_clusters.read_cluster_heads.count(make_pair(read_num, cluster_head)) == 0) {
                parent_clusters.read_cluster_heads[make_pair(read_num, cluster_head)] = distances_to_parent;
            }
        }

        /*Now go through clusters on the second child, and see if they can be combined with clusters on the first child
         */
        for (auto& cluster_head_and_distance : child_clusters2.read_cluster_heads) {

            size_t read_num = cluster_head_and_distance.first.first;
            size_t cluster_num = tree_state.read_union_find[read_num].find_group(cluster_head_and_distance.first.second);
            pair<size_t, size_t> distances = cluster_head_and_distance.second;
            size_t new_dist_left = std::min(SnarlDistanceIndex::sum({distances.first,child_clusters2.distance_start_left}), 
                                            SnarlDistanceIndex::sum({distances.second,child_clusters2.distance_start_right}));
            size_t new_dist_right = std::min(SnarlDistanceIndex::sum({distances.first,child_clusters2.distance_end_left}), 
                                            SnarlDistanceIndex::sum({distances.second,child_clusters2.distance_end_right}));
            pair<size_t, size_t> distances_to_parent = make_pair(new_dist_left, new_dist_right);

            if (parent_clusters.read_cluster_heads.count(make_pair(read_num, cluster_num)) > 0) {
                distances_to_parent = make_pair(
                    std::min(new_dist_left, parent_clusters.read_cluster_heads[make_pair(read_num, cluster_num)].first),
                    std::min(new_dist_right, parent_clusters.read_cluster_heads[make_pair(read_num, cluster_num)].second));
            }

            //Check if the left of 1 can connect with the left of 2
//#ifdef DEBUG_CLUSTER
//            cerr << "\t\tCheck any clusters in 2 going left to left of 1" << endl;
//#endif
            compare_and_combine_clusters (read_num, cluster_num, 
                SnarlDistanceIndex::sum({distances.first,distance_left_left,child_clusters1.read_best_left[read_num]}), 
                SnarlDistanceIndex::sum({distances.first,distance_left_left,child_clusters1.fragment_best_left}), 
                distances_to_parent, new_cluster_left_left_by_read[read_num], new_cluster_left_left_fragment);
            //Check if the left of 1 can connect with the right of 2
//#ifdef DEBUG_CLUSTER
//            cerr << "\t\tCheck any clusters in 2 going right to left of 1" << endl;
//#endif
            compare_and_combine_clusters (read_num, cluster_num, 
                SnarlDistanceIndex::sum({distances.second,distance_left_right,child_clusters1.read_best_left[read_num]}),
                SnarlDistanceIndex::sum({distances.second,distance_left_right,child_clusters1.fragment_best_left}),
                distances_to_parent, new_cluster_left_right_by_read[read_num], new_cluster_left_right_fragment);
            //Check if the right of 1 can connect with the right of 2
//#ifdef DEBUG_CLUSTER
//            cerr << "\t\tCheck any clusters in 2 going right to right of 1" << endl;
//#endif
            compare_and_combine_clusters (read_num, cluster_num, 
                SnarlDistanceIndex::sum({distances.second,distance_right_right,child_clusters1.read_best_right[read_num]}),
                SnarlDistanceIndex::sum({distances.second,distance_right_right,child_clusters1.fragment_best_right}),
                distances_to_parent, new_cluster_right_right_by_read[read_num], new_cluster_right_right_fragment);
            //Check if the right of 1 can connect with the left of 2
//#ifdef DEBUG_CLUSTER
//            cerr << "\t\tCheck any clusters in 2 going left to right of 1" << endl;
//#endif
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
        }

        //And add back in the new cluster heads
        for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
            if (new_cluster_left_left_by_read.at(read_num).first.first != std::numeric_limits<size_t>::max()) {
                pair<size_t, size_t> old_distances = parent_clusters.read_cluster_heads.count(new_cluster_left_left_by_read.at(read_num).first) == 0
                    ? make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max())
                    : parent_clusters.read_cluster_heads.at(new_cluster_left_left_by_read.at(read_num).first);
                parent_clusters.read_cluster_heads[new_cluster_left_left_by_read.at(read_num).first] = 
                    make_pair(std::min(new_cluster_left_left_by_read.at(read_num).second.first, old_distances.first),
                              std::min(new_cluster_left_left_by_read.at(read_num).second.second, old_distances.second));
            }
            if (new_cluster_right_right_by_read.at(read_num).first.first != std::numeric_limits<size_t>::max()) {
                pair<size_t, size_t> old_distances = parent_clusters.read_cluster_heads.count(new_cluster_right_right_by_read.at(read_num).first) == 0
                    ? make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max())
                    : parent_clusters.read_cluster_heads.at(new_cluster_right_right_by_read.at(read_num).first);
                parent_clusters.read_cluster_heads[new_cluster_right_right_by_read.at(read_num).first] = 
                    make_pair(std::min(new_cluster_right_right_by_read.at(read_num).second.first, old_distances.first),
                              std::min(new_cluster_right_right_by_read.at(read_num).second.second, old_distances.second));
            }
            if (new_cluster_left_right_by_read.at(read_num).first.first != std::numeric_limits<size_t>::max()) {
                pair<size_t, size_t> old_distances = parent_clusters.read_cluster_heads.count(new_cluster_left_right_by_read.at(read_num).first) == 0
                    ? make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max())
                    : parent_clusters.read_cluster_heads.at(new_cluster_left_right_by_read.at(read_num).first);
                parent_clusters.read_cluster_heads[new_cluster_left_right_by_read.at(read_num).first] = 
                    make_pair(std::min(new_cluster_left_right_by_read.at(read_num).second.first, old_distances.first),
                              std::min(new_cluster_left_right_by_read.at(read_num).second.second, old_distances.second));
            }
            if (new_cluster_right_left_by_read.at(read_num).first.first != std::numeric_limits<size_t>::max()) {
                pair<size_t, size_t> old_distances = parent_clusters.read_cluster_heads.count(new_cluster_right_left_by_read.at(read_num).first) == 0
                    ? make_pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max())
                    : parent_clusters.read_cluster_heads.at(new_cluster_right_left_by_read.at(read_num).first);
                parent_clusters.read_cluster_heads[new_cluster_right_left_by_read.at(read_num).first] = 
                    make_pair(std::min(new_cluster_right_left_by_read.at(read_num).second.first, old_distances.first),
                              std::min(new_cluster_right_left_by_read.at(read_num).second.second, old_distances.second));
            }
        }

/*
        for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
            if (new_cluster_left_left_by_read.at(read_num).first.first != std::numeric_limits<size_t>::max()) {
                parent_clusters.read_cluster_heads[new_cluster_left_left_by_read.at(read_num).first] = new_cluster_left_left_by_read.at(read_num).second;
            }
            if (new_cluster_right_right_by_read.at(read_num).first.first != std::numeric_limits<size_t>::max()) {
                parent_clusters.read_cluster_heads[new_cluster_right_right_by_read.at(read_num).first] = new_cluster_right_right_by_read.at(read_num).second;
            }
            if (new_cluster_left_right_by_read.at(read_num).first.first != std::numeric_limits<size_t>::max()) {
                parent_clusters.read_cluster_heads[new_cluster_left_right_by_read.at(read_num).first] = new_cluster_left_right_by_read.at(read_num).second;
            }
            if (new_cluster_right_left_by_read.at(read_num).first.first != std::numeric_limits<size_t>::max()) {
                parent_clusters.read_cluster_heads[new_cluster_right_left_by_read.at(read_num).first] = new_cluster_right_left_by_read.at(read_num).second;
            }
        }
        */


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


    NewSnarlSeedClusterer::NodeClusters NewSnarlSeedClusterer::cluster_one_snarl(
                    TreeState& tree_state, net_handle_t snarl_handle) const {
        //Get the clusters on this snarl, assumes that all of the snarls children have been clustered already.
        
#ifdef DEBUG_CLUSTER
        cerr << "Finding clusters on snarl " << distance_index.net_handle_as_string(snarl_handle) << endl;
#endif

        //Keep track of all clusters on this snarl
        NodeClusters snarl_clusters(snarl_handle, tree_state.all_seeds->size());



        //Get the children of this snarl and their clusters
        vector< NodeClusters>& children = tree_state.snarl_to_children[snarl_handle];

        for (size_t i = 0; i < children.size() ; i++) {
            //Go through each child node of the netgraph

            NodeClusters& child_clusters_i = children[i];

            for (size_t j = 0 ; j <= i ; j++){
                //Go through other child net graph nodes up to and including i

                //Get the other node and its clusters
                NodeClusters& child_clusters_j = children[j];

#ifdef DEBUG_CLUSTER
                cerr << "\tComparing two children of " << distance_index.net_handle_as_string(snarl_handle) << ": " 
                     << distance_index.net_handle_as_string(child_clusters_i.containing_net_handle) << " and " 
                     << distance_index.net_handle_as_string(child_clusters_j.containing_net_handle) << endl;
                     


#endif

                compare_and_combine_cluster_on_child_structures(tree_state, child_clusters_i, 
                        child_clusters_j, snarl_clusters);

            }
        }
#ifdef DEBUG_CLUSTER
        cerr << "\tFound clusters on " << distance_index.net_handle_as_string(snarl_handle) << endl;
        cerr << "\t   with best left and right values: " << snarl_clusters.fragment_best_left << " "
             << snarl_clusters.fragment_best_right << endl;
        bool got_left = false;
        bool got_right = false;
        for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
            cerr << "\t\tfor read num " << read_num << " best left: " << snarl_clusters.read_best_left[read_num] << " best right: " << snarl_clusters.read_best_right[read_num] << endl;
            bool got_read_left=false;
            bool got_read_right = false;
            bool any_clusters = false;
            for (pair<pair<size_t, size_t>, pair<size_t,size_t>> c : snarl_clusters.read_cluster_heads) {
                if (c.first.first == read_num) {
                    any_clusters = true;
                    pair<size_t, size_t> dists = c.second;
                    cerr << "\t\t" << c.first.first << ":"<<c.first.second << ": left: " << dists.first << " right : " << dists.second << ": ";
                    bool has_seeds = false;
                    for (size_t x = 0 ; x < tree_state.all_seeds->at(c.first.first)->size() ; x++) {
                        if (tree_state.read_union_find[c.first.first].find_group(x) == c.first.second) {
                            cerr << tree_state.all_seeds->at(c.first.first)->at(x).pos << " ";
                            has_seeds = true;
                        }
                    }
                    assert(dists.first == std::numeric_limits<size_t>::max() || dists.first >= snarl_clusters.read_best_left[read_num]);
                    assert(dists.second == std::numeric_limits<size_t>::max() || dists.second >= snarl_clusters.read_best_right[read_num]);
                    assert(dists.first == std::numeric_limits<size_t>::max() || dists.first >= snarl_clusters.fragment_best_left);
                    assert(dists.second == std::numeric_limits<size_t>::max() || dists.second >= snarl_clusters.fragment_best_right);
                    if (dists.first == snarl_clusters.fragment_best_left) {got_left = true;}
                    if (dists.second == snarl_clusters.fragment_best_right) {got_right = true;}
                    if (dists.first == snarl_clusters.read_best_left[read_num]) {got_read_left = true;}
                    if (dists.second == snarl_clusters.read_best_right[read_num]) {got_read_right = true;}
                    cerr << endl;
                    assert(has_seeds);
                }
            }
            assert(!any_clusters ||got_read_left ||  snarl_clusters.read_best_left[read_num] == std::numeric_limits<size_t>::max());
            assert(!any_clusters ||got_read_right ||  snarl_clusters.read_best_right[read_num] == std::numeric_limits<size_t>::max());
        }
        assert(got_left);
        assert(got_right);

        for (pair<pair<size_t, size_t>, pair<size_t, size_t>> group_id : snarl_clusters.read_cluster_heads) {
            assert (group_id.first.second == tree_state.read_union_find[group_id.first.first].find_group(group_id.first.second));
        }
#endif
        return snarl_clusters;
    };



    NewSnarlSeedClusterer::NodeClusters NewSnarlSeedClusterer::cluster_one_chain(TreeState& tree_state, net_handle_t chain_handle) const {
        if (distance_index.is_trivial_chain(chain_handle)){
            //If this is just a node pretending to be a chain, cluster the node and claim it's a chain
#ifdef DEBUG_CLUSTER
            cerr << "\tClustering a chain that is really just a node, so just cluster the node" << endl;
#endif

            
            //This will hold the clusters of the chain
            assert(tree_state.chain_to_children[chain_handle].size() == 1); 

            NodeClusters& child_cluster = tree_state.chain_to_children[chain_handle].back();
            cluster_one_node(tree_state, child_cluster);

            NodeClusters chain_clusters = std::move(child_cluster);
            chain_clusters.containing_net_handle = chain_handle;

            return chain_clusters;

        }


#ifdef DEBUG_CLUSTER
        cerr << "Cluster chain " << distance_index.net_handle_as_string(chain_handle) << endl;
        cerr << "\t chain has " << tree_state.chain_to_children[chain_handle].size() << " children" << endl;
#endif

        /*Go through the chain child by child 
        * 
        * For each child, 
        *  - check if the clusters in the child can be combined with each other by walking out then back in through the chain
        *       -Update distances to the ends of the child (taking into account the distances to loop around in the chain)
        *  - compare and combine with clusters of the chain that get build as we walk along it
        *       Each child needs to be compared everything found since the last node, since if it could reach any seed on 
        *       the chain to the left of the last node found, then the last node would also be able to reach that seed so
        *       it would be combined anyway
        *  - update the distances to the ends of the chain
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

        //This will hold the clusters of the chain
        NodeClusters chain_clusters(chain_handle, tree_state.all_seeds->size());


        // As we walk along the chain, we need to remember all clusters since the last node
        // This stores all clusters we found since the last node
        vector<NodeClusters*> last_child_clusters;

        //Get the children of this chain from the tree state. They will be ordered by their order in the chain
        //Maps the offset in the chain records to the node cluster (offset doesn't matter for anything other than putting them in order)
        vector<NodeClusters>& children_in_chain = tree_state.chain_to_children[chain_handle];

        for (NodeClusters& child_clusters: children_in_chain) {
            /*
             * Snarls and nodes are in the order that they are traversed in the chain
             */

            if (distance_index.is_node(child_clusters.containing_net_handle)){
                cluster_one_node(tree_state, child_clusters);
            }

            //See if this clusters of this child can be combined with each other
            compare_and_combine_cluster_on_child_structures(tree_state, child_clusters, 
                        child_clusters, chain_clusters);

            //Go through all child clusters since we last saw a node cluster
            for (NodeClusters* previous_child_node : last_child_clusters) {

                compare_and_combine_cluster_on_child_structures(tree_state, child_clusters, 
                        *previous_child_node, chain_clusters);
            }

            if (distance_index.is_node(child_clusters.containing_net_handle)) {
                //If the current child is a node, then forget all the previous children
                if (tree_state.fragment_distance_limit == 0) {
                    last_child_clusters.clear();
                } else {
                    //If we're also looking for fragment clusters, then only do this if the node had clusters on both reads 
                    vector<bool> has_cluster_of_read (tree_state.all_seeds->size(), false);
                    for (const pair<pair<size_t, size_t>, pair<size_t, size_t>>& cluster_indices : child_clusters.read_cluster_heads) {
                        has_cluster_of_read.at(cluster_indices.first.first) = true;
                        if (std::accumulate(has_cluster_of_read.begin(), has_cluster_of_read.end(), 0) == has_cluster_of_read.size()) {
                            break;
                        }
                    }
                    if (std::accumulate(has_cluster_of_read.begin(), has_cluster_of_read.end(), 0) == has_cluster_of_read.size()) {
                        last_child_clusters.clear();
                    }
                }

            }
            //Remember that we saw this child so we can compare against it later
            last_child_clusters.emplace_back(&child_clusters);
        }

        //If the chain loops, then we also have to compare the first thing we saw to the last things
        if (distance_index.is_looping_chain(chain_handle)) {
            NodeClusters& first_child_cluster = children_in_chain.front(); 
            for (NodeClusters* previous_child_node : last_child_clusters) {
                compare_and_combine_cluster_on_child_structures(tree_state, first_child_cluster, 
                        *previous_child_node, chain_clusters);
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
            for (pair<pair<size_t, size_t>, pair<size_t,size_t>> c : chain_clusters.read_cluster_heads) {
                if (c.first.first == read_num) {
                    any_clusters = true;
                    pair<size_t, size_t> dists = c.second;
                    cerr << "\t\t" << c.first.first << ":"<<c.first.second << ": left: " << dists.first << " right : " << dists.second << ": ";
                    bool has_seeds = false;
                    for (size_t x = 0 ; x < tree_state.all_seeds->at(c.first.first)->size() ; x++) {
                        if (tree_state.read_union_find[c.first.first].find_group(x) == c.first.second) {
                            cerr << tree_state.all_seeds->at(c.first.first)->at(x).pos << " ";
                            has_seeds = true;
                        }
                    }
                    assert(dists.first == std::numeric_limits<size_t>::max() || dists.first >= chain_clusters.read_best_left[read_num]);
                    assert(dists.second == std::numeric_limits<size_t>::max() || dists.second >= chain_clusters.read_best_right[read_num]);
                    assert(dists.first == std::numeric_limits<size_t>::max() || dists.first >= chain_clusters.fragment_best_left);
                    assert(dists.second == std::numeric_limits<size_t>::max() || dists.second >= chain_clusters.fragment_best_right);
                    if (dists.first == chain_clusters.fragment_best_left) {got_left = true;}
                    if (dists.second == chain_clusters.fragment_best_right) {got_right = true;}
                    if (dists.first == chain_clusters.read_best_left[read_num]) {got_read_left = true;}
                    if (dists.second == chain_clusters.read_best_right[read_num]) {got_read_right = true;}
                    cerr << endl;
                    assert(has_seeds);
                }
            }
            assert(!any_clusters ||got_read_left ||  chain_clusters.read_best_left[read_num] == std::numeric_limits<size_t>::max());
            assert(!any_clusters ||got_read_right ||  chain_clusters.read_best_right[read_num] == std::numeric_limits<size_t>::max());
        }
        assert(got_left);
        assert(got_right);

        for (pair<pair<size_t, size_t>, pair<size_t, size_t>> group_id : chain_clusters.read_cluster_heads) {
            assert (group_id.first.second == tree_state.read_union_find[group_id.first.first].find_group(group_id.first.second));
        }
#endif
        return chain_clusters;
    }

    //Cluster the root
    //all children of the root will be in tree_state.root_children
    //This is basically cluster_one_snarl except the snarl is the root, which has no boundary nodes
    void NewSnarlSeedClusterer::cluster_root(TreeState& tree_state) const { 
#ifdef DEBUG_CLUSTER
        cerr << "Finding clusters on the root" << endl;
#endif

        //Keep track of all clusters on the root
        NodeClusters root_clusters(distance_index.get_root(), tree_state.all_seeds->size());



        for (size_t i = 0; i < tree_state.root_children.size() ; i++) {
            //Go through each child node of the netgraph

            NodeClusters& child_clusters_i = tree_state.root_children[i];

            for (size_t j = 0 ; j <= i ; j++){
                //Go through other child net graph nodes up to and including i

                //Get the other node and its clusters
                NodeClusters& child_clusters_j = tree_state.root_children[j];

#ifdef DEBUG_CLUSTER
                cerr << "\tComparing two children of the root: " 
                     << distance_index.net_handle_as_string(child_clusters_i.containing_net_handle) << " and " 
                     << distance_index.net_handle_as_string(child_clusters_j.containing_net_handle) << endl;
#endif

                compare_and_combine_cluster_on_child_structures(tree_state, child_clusters_i,
                        child_clusters_j, root_clusters, true);

            }
        }
#ifdef DEBUG_CLUSTER
        cerr << "\tFound clusters on the root" << endl;
        for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
            cerr << "\t for read num " << read_num << endl;
            for (pair<pair<size_t, size_t>, pair<size_t,size_t>> c : root_clusters.read_cluster_heads) {
                if (c.first.first == read_num) {
                    cerr << "\t\t" << c.first.first << ":"<<c.first.second << ":  ";
                    for (size_t x = 0 ; x < tree_state.all_seeds->at(c.first.first)->size() ; x++) {
                        if (tree_state.read_union_find[c.first.first].find_group(x) == c.first.second) {
                            cerr << tree_state.all_seeds->at(c.first.first)->at(x).pos << " ";
                        }
                    }
                    cerr << endl;
                }
            }
        }

        for (pair<pair<size_t, size_t>, pair<size_t, size_t>> group_id : root_clusters.read_cluster_heads) {
            assert (group_id.first.second == tree_state.read_union_find[group_id.first.first].find_group(group_id.first.second));
        }
#endif
    }

    void NewSnarlSeedClusterer::add_child_to_vector(hash_map<net_handle_t, vector<NodeClusters>>& parent_to_child_map, const net_handle_t& parent,
            NodeClusters child_cluster) const {
        if (parent_to_child_map.count(parent) == 0) {
            //Make sure that parent shows up in the hash_map
            vector<NodeClusters> empty_vector;
            parent_to_child_map.emplace(parent, std::move(empty_vector));
        }
        if (distance_index.is_snarl(parent)){
            //If the parent is a snarl, then the order doesn't matter
            parent_to_child_map.at(parent).emplace_back(std::move(child_cluster));
        } else {
            assert(distance_index.is_chain(parent));
            //Helper function to insert a NodeClusters into a sorted vector<NodeClusters> where the vector contains only 
            //children of the same chain. The items will be sorted by their position in the chain
            auto insert_in_order = [&](vector<NodeClusters>& input_vector, NodeClusters& item) {
                //Get an iterator to where the thing we're inserting should go in the vector 
                std::vector<NodeClusters>::iterator insert_itr = std::upper_bound(input_vector.begin(), input_vector.end(),
                   0, [&](size_t val, NodeClusters vector_item) {
                       //This should return true if the vector_item goes after the item with rank val
                       return distance_index.is_ordered_in_chain(item.containing_net_handle, vector_item.containing_net_handle);
                   });
                //Add the item to the vector in sorted order
                input_vector.insert(insert_itr, std::move(item));
            };
            insert_in_order(parent_to_child_map.at(parent), child_cluster);
        }
    }
/*
    void NewSnarlSeedClusterer::cluster_only_top_level_chain_seeds(TreeState& tree_state, 
                                                                  vector<pair<size_t, size_t>>& seed_clusters) const {
#ifdef DEBUG_CLUSTER
        cerr << "Clustering top-level seeds" << endl;
#endif


        std::sort(seed_clusters.begin(), seed_clusters.end(), [&](pair<size_t, size_t> item1, pair<size_t, size_t> item2) {
                   size_t offset1 = tree_state.all_seeds->at(item1.first)->at(item1.second).offset ;
                   size_t offset2 = tree_state.all_seeds->at(item2.first)->at(item2.second).offset ;
                   return offset1 < offset2; 
               });

        tuple<size_t, size_t, size_t> last_fragment (-1, 0, 0);

        vector<tuple<size_t, size_t, size_t>> last_by_read (tree_state.all_seeds->size(), make_tuple(-1, 0, 0));

        for (pair<size_t, size_t> seed_cluster : seed_clusters) {
            size_t offset = tree_state.all_seeds->at(seed_cluster.first)->at(seed_cluster.second).offset ;

#ifdef DEBUG_CLUSTER
            cerr << tree_state.all_seeds->at(seed_cluster.first)->at(seed_cluster.second).pos << " on read " << seed_cluster.first << " at offset " <<
                    tree_state.all_seeds->at(seed_cluster.first)->at(seed_cluster.second).offset << endl; 
#endif
            if (tree_state.fragment_distance_limit != 0) {
                if (std::get<0>(last_fragment) != -1 && offset - std::get<0>(last_fragment) <= tree_state.fragment_distance_limit) {
                    //If this is close enough to the last thing we saw
                    tree_state.fragment_union_find.union_groups(
                            std::get<2>(last_fragment) + tree_state.read_index_offsets[std::get<1>(last_fragment)], 
                            seed_cluster.second + tree_state.read_index_offsets[seed_cluster.first]);
#ifdef DEBUG_CLUSTER
                    cerr << "\tCombine fragment with last seed seen, " << tree_state.all_seeds->at(std::get<1>(last_fragment))->at(std::get<2>(last_fragment)).pos <<
                            " at offset " << std::get<0>(last_fragment) << endl; 
#endif
                }
                last_fragment = make_tuple(offset, seed_cluster.first, seed_cluster.second);

            }

            if (std::get<0>(last_by_read[seed_cluster.first]) != -1 
                && offset - std::get<0>(last_by_read[seed_cluster.first]) <= tree_state.read_distance_limit) {
                //If this is close enough to the last one we saw for this read
                tree_state.read_union_find[seed_cluster.first].union_groups(std::get<2>(last_by_read[seed_cluster.first]),
                                                                             seed_cluster.second);
#ifdef DEBUG_CLUSTER
                    cerr << "\tCombine read with last seed seen, " << 
                        tree_state.all_seeds->at(std::get<1>(last_by_read[seed_cluster.first]))->at(std::get<2>(last_by_read[seed_cluster.first])).pos << 
                        " at offset " << std::get<0>(last_by_read[seed_cluster.first]) << endl; 
#endif
            }
            last_by_read[seed_cluster.first] = make_tuple(offset, seed_cluster.first, seed_cluster.second);
            
        }
    };


    hash_set<pair<size_t, size_t>> NewSnarlSeedClusterer::cluster_simple_snarl(TreeState& tree_state, vector<tuple<id_t, bool, size_t, size_t, size_t>> nodes, 
                                size_t loop_left, size_t loop_right, size_t snarl_length) const {
        //Cluster a top-level simple snarl and save the distances to the ends of the node in tree_state.read_cluster_dists 
        //Returns a vector of cluster heads (<read_num, seed num>)
#ifdef DEBUG_CLUSTER
        cerr << "cluster simple snarl " << endl;
        cerr << " loop distances: " << loop_left << " " << loop_right  << " snarl length " << snarl_length << endl;
#endif


 
        //Get the clusters of each node and find the best distances of any of them       
        vector<pair<size_t, size_t>> all_cluster_heads;
        size_t fragment_best_left = -1; size_t fragment_best_right = -1;
        vector<size_t> best_left (tree_state.all_seeds->size(), -1);
        vector<size_t> best_right (tree_state.all_seeds->size(), -1);


        for (tuple<id_t, bool, size_t, size_t, size_t> node : nodes) {
            id_t node_id = std::get<0>(node);
            size_t node_len = std::get<4>(node);

            //Cluster this node
            NewSnarlSeedClusterer::NodeClusters clusters =  cluster_one_node(tree_state, node_id, node_len);

            //Add the new cluster heads we found and update the best distances
            for (pair<size_t, size_t> cluster_head : clusters.read_cluster_heads) {
                all_cluster_heads.emplace_back(cluster_head);
            }
            fragment_best_left = std::min(fragment_best_left, clusters.fragment_best_left);
            fragment_best_right = std::min(fragment_best_right, clusters.fragment_best_right);
            for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num ++) {
                best_left[read_num] = std::min(best_left[read_num], clusters.read_best_left[read_num]);
                best_right[read_num] = std::min(best_right[read_num], clusters.read_best_right[read_num]); 
            }
        }

        //Now go through all the clusters we just found and see if any of them can be combined
        //by taking loops in the chain


        hash_set<pair<size_t, size_t>> result;
            
        //Cluster heads of clusters formed by combining clusters
        size_t combined_fragment_left = -1; size_t combined_fragment_right = -1;
        vector<size_t> combined_read_left (tree_state.all_seeds->size(), -1);
        vector<size_t> combined_read_right (tree_state.all_seeds->size(), -1);

        for (pair<size_t, size_t>& cluster_head: all_cluster_heads) {
#ifdef DEBUG_CLUSTER
            cerr << "\tfound intermediate cluster head " << tree_state.all_seeds->at(cluster_head.first)->at(cluster_head.second).pos 
                << ": dist left: " << tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].first 
                << ", dist right: " << tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].second << endl;
#endif

            result.emplace(cluster_head);

            size_t dist_left = tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].first;
            size_t dist_right = tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].second;
            //Check for loop to the left
            if (loop_left != -1) {

                //Combine fragment 
                if (tree_state.fragment_distance_limit != 0 && dist_left != -1 && fragment_best_left != -1 &&
                        dist_left + loop_left + fragment_best_left - 1 <= tree_state.fragment_distance_limit) {
                    if (combined_fragment_left == -1) {
                        combined_fragment_left = tree_state.read_index_offsets[cluster_head.first] + cluster_head.second; 
                    } else {
                        tree_state.fragment_union_find.union_groups(combined_fragment_left, 
                                    tree_state.read_index_offsets[cluster_head.first] + cluster_head.second);
                    }
                }

                //Combine read
                if (dist_left != -1 && best_left[cluster_head.first] != -1 &&
                        dist_left + loop_left + best_left[cluster_head.first] - 1 <= tree_state.read_distance_limit) {
                    if (combined_read_left[cluster_head.first] == -1) {
                        combined_read_left[cluster_head.first] = cluster_head.second;
                        tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].second = 
                                std::min(dist_right, dist_left + loop_left + snarl_length);
                    } else {
                        size_t best_left = std::min(dist_left, 
                                tree_state.read_cluster_dists[cluster_head.first][combined_read_left[cluster_head.first]].first);
                        size_t best_right =std::min(dist_right, 
                                tree_state.read_cluster_dists[cluster_head.first][combined_read_left[cluster_head.first]].second);

                        tree_state.read_union_find[cluster_head.first].union_groups(cluster_head.second, combined_read_left[cluster_head.first]);
                        if (tree_state.read_union_find[cluster_head.first].find_group(cluster_head.second) == cluster_head.second) {
                            result.erase(make_pair(cluster_head.first, combined_read_left[cluster_head.first]));
                            combined_read_left[cluster_head.first] = cluster_head.second;
                        } else {
                            result.erase(cluster_head);
                        }


                        tree_state.read_cluster_dists[cluster_head.first][combined_read_left[cluster_head.first]].first = 
                            std::min(best_left, loop_right == -1 ? -1 : best_right + loop_right + snarl_length);
                        tree_state.read_cluster_dists[cluster_head.first][combined_read_left[cluster_head.first]].second = 
                            std::min(best_right, best_left + loop_left + snarl_length);
                    }
                } else {
                    tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].second = 
                            std::min(dist_right, dist_left + loop_left + snarl_length);
                }
            }

            //Check for loop to the right
            if (loop_right != -1) {

                //Combine fragment
                if (tree_state.fragment_distance_limit != 0 && dist_right != -1 && fragment_best_right != -1 &&
                        dist_right + loop_right + fragment_best_right - 1 <= tree_state.fragment_distance_limit) {
                    if (combined_fragment_right == -1) {
                        combined_fragment_right = tree_state.read_index_offsets[cluster_head.first] + cluster_head.second;
                    } else {
                        tree_state.fragment_union_find.union_groups(combined_fragment_right, 
                                    tree_state.read_index_offsets[cluster_head.first] + cluster_head.second);
                    }
                }

                //Combine read
                if (dist_right != -1 && best_right[cluster_head.first] != -1 &&
                    dist_right  + loop_right + best_right[cluster_head.first] - 1 <= tree_state.read_distance_limit) {
                    if (combined_read_right[cluster_head.first] == -1) {
                        combined_read_right[cluster_head.first] = cluster_head.second;
                        tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].first = 
                                std::min(dist_left, dist_right + loop_right + snarl_length);
                    } else {
                        size_t best_left = std::min(dist_left, 
                                tree_state.read_cluster_dists[cluster_head.first][combined_read_right[cluster_head.first]].first);
                        size_t best_right = std::min(dist_right, 
                                tree_state.read_cluster_dists[cluster_head.first][combined_read_right[cluster_head.first]].second);
                        tree_state.read_union_find[cluster_head.first].union_groups(cluster_head.second, combined_read_right[cluster_head.first]);
                        if (tree_state.read_union_find[cluster_head.first].find_group(cluster_head.second) == cluster_head.second) {
                            result.erase(make_pair(cluster_head.first, combined_read_right[cluster_head.first]));
                            combined_read_right[cluster_head.first] = cluster_head.second;
                        } else {
                            result.erase(cluster_head);
                        }


                        tree_state.read_cluster_dists[cluster_head.first][combined_read_right[cluster_head.first]].first = 
                            std::min(best_left, best_right + loop_right + snarl_length);
                        tree_state.read_cluster_dists[cluster_head.first][combined_read_right[cluster_head.first]].second = 
                            std::min(best_right, loop_left == -1 ? -1 : best_left + loop_left + snarl_length);
                    }
                } else {
                    tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].first = 
                            std::min(dist_left, dist_right + loop_right + snarl_length);
                }
            }
        }
#ifdef DEBUG_CLUSTER
        cerr << "Found clusters on simple snarl: " << endl;
        for ( pair<size_t, size_t> cluster_head : result) {
            cerr << "\t" << tree_state.all_seeds->at(cluster_head.first)->at(cluster_head.second).pos 
                << ": dist left: " << tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].first 
                << ", dist right: " << tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].second << endl;
        }
#endif

        return result;

    };
    */
}
