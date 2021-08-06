#include "seed_clusterer.hpp"

#include <algorithm>

//#define DEBUG_CLUSTER
namespace vg {

    SnarlSeedClusterer::SnarlSeedClusterer( SnarlDistanceIndex& dist_index) :
                                            dist_index(dist_index){
    };

    vector<SnarlSeedClusterer::Cluster> SnarlSeedClusterer::cluster_seeds (const vector<Seed>& seeds, int64_t read_distance_limit) const {
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

    vector<vector<SnarlSeedClusterer::Cluster>> SnarlSeedClusterer::cluster_seeds (
                  const vector<vector<Seed>>& all_seeds, int64_t read_distance_limit,
                  int64_t fragment_distance_limit) const {

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


    tuple<vector<structures::UnionFind>, structures::UnionFind> SnarlSeedClusterer::cluster_seeds_internal (
                  const vector<const vector<Seed>*>& all_seeds, int64_t read_distance_limit,
                  int64_t fragment_distance_limit) const {
        /* Given a vector of seeds and a limit, find a clustering of seeds where
         * seeds that are closer than the limit cluster together.
         * Returns a vector of cluster assignments
         */
#ifdef DEBUG_CLUSTER
cerr << endl << endl << endl << endl << "New cluster calculation:" << endl;
#endif
        if (fragment_distance_limit != 0 &&
            fragment_distance_limit < read_distance_limit) {
            throw std::runtime_error("Fragment distance limit must be greater than read distance limit");
        }

        //For each level of the snarl tree, maps chains (net_handle_t of the snarl) at that level to children 
        //(net_handle_t of the child) belonging to the chain. Starts out with nodes from get_nodes, 
        //then adds snarls as we walk up the snarl tree
        //This is later used to populate snarl_to_node in the tree state
        vector<hash_map<net_handle_t, vector<pair<net_handle_t, NodeClusters>>>> chain_to_children_by_level;
        chain_to_children_by_level.reserve(dist_index.get_max_tree_depth()+1);



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
        tree_state.chain_to_children = std::move(chain_to_children_by_level[snarl_to_nodes_by_level.size() - 1]);

        for (int depth = snarl_to_nodes_by_level.size() - 1 ; depth >= 0 ; depth --) {
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
            cluster_chain_level(tree_state, depth);

            //And cluster all the snarls, record the parents of these snarls
            cluster_snarl_level(tree_state, depth);


            // Swap buffer over for the next level
            tree_state.chain_to_children = std::move(tree_state.parent_chain_to_children);
            tree_state.snarl_to_children.clear();
        }


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
                    pos_t rev1 = make_pos_t(get_id(pos1), !is_rev(pos1), dist_index.node_length(get_id(pos1)) - get_offset(pos1) - 1);

                    for (size_t i2 = 0 ; i2 < i1 ; i2++) {
                    
                        size_t d = group[i2];
                    
                        pos_t pos2 = tree_state.all_seeds->at(read_num)->at(d).pos;
                        pos_t rev2 = make_pos_t(get_id(pos2), !is_rev(pos2), dist_index.node_length(get_id(pos2))- get_offset(pos2) - 1);
                        int64_t d1 = dist_index.min_distance(pos1, pos2);
                        int64_t d2 = std::min(d1, dist_index.min_distance(pos1, rev2));
                        int64_t d3 = std::min(d2, dist_index.min_distance(rev1, rev2));
                        int64_t d4 = std::min(d3, dist_index.min_distance(rev1, pos2));
                        if (d4 != -1 && d4 <= tree_state.read_distance_limit) {
                        
                             uf.union_groups(i1, i2);
                        }
                    }
                    for (size_t g2 = 0 ; g2 < all_groups.size() ; g2 ++) {
                        if (g2 != g1) {
                            auto group2 = all_groups[g2];
                            for (size_t d : group2) {
                               pos_t pos2 = tree_state.all_seeds->at(read_num)->at(d).pos;
                               pos_t rev2 = make_pos_t(get_id(pos2), !is_rev(pos2), dist_index.node_length(get_id(pos2)) - get_offset(pos2) - 1);
                               
                               int64_t d1 = dist_index.min_distance(pos1, pos2);
                               int64_t d2 = std::min(d1, dist_index.min_distance(pos1, rev2));
                               int64_t d3 = std::min(d2, dist_index.min_distance(rev1, rev2));
                               int64_t d4 = std::min(d3, dist_index.min_distance(rev1, pos2));
                               
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
    //tree. snarl_to_nodes_by_level will be used to populate snarl_to_nodes
    //in the tree state as each level is processed
    void SnarlSeedClusterer::get_nodes( TreeState& tree_state,
              vector<hash_map<size_t,vector<pair<net_handle_t, NodeClusters>>>>&
                                                 snarl_to_nodes_by_level) const {
#ifdef DEBUG_CLUSTER
cerr << "Nested positions: " << endl << "\t";
#endif

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
                     net_handle_t parent = dist_index.get_parent(id);
                     size_t depth = distance_index.get_depth(parent);
                     if (depth > snarl_to_nodes_by_level.size()) {
                         snarl_to_nodes_by_level.resize(depth);
                     }
                     chain_to_children_by_level[depth][parent].insert(NodeClusters(get_node_net_handle(id), tree_state.all_seeds->size()));
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

                    size_t component = dist_index.get_connected_component(id);
                    if (tree_state.component_to_index.count(component) == 0) {
                        tree_state.component_to_index.emplace(component, tree_state.top_level_seed_clusters.size());
                        tree_state.top_level_seed_clusters.emplace_back();
                        tree_state.simple_snarl_to_nodes_by_component.emplace_back();
                    }
                    hash_map<size_t, vector<tuple<id_t, bool, int64_t, int64_t, int64_t>>>& simple_snarl_to_nodes = 
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

        if (snarl_to_nodes_by_level.empty()) {
            snarl_to_nodes_by_level.resize(1);
        }
#ifdef DEBUG_CLUSTER
        cerr << endl << "Top-level seeds:" << endl << "\t";
        for (size_t component_num = 0 ; component_num < tree_state.top_level_seed_clusters.size() ; component_num++) {
            for (pair<size_t, size_t> seed_index : tree_state.top_level_seed_clusters[component_num]) {
                cerr << seed_index.first << ":" << tree_state.all_seeds->at(seed_index.first)->at(seed_index.second).pos << ", ";
            }
        }
        cerr << endl;
#endif
    }


    //Cluster all of the snarls in tree_state at the given depth
    //Assumes that all the children of the snarls have been clustered already and are present in tree_state.snarls_to_children
    void SnarlSeedClusterer::cluster_snarl_level(TreeState& tree_state, size_t depth) const {

        for (auto& kv : tree_state.snarl_to_children){
            //Go through each of the snarls at this level, cluster them,
            //and find which chains they belong to, if any
            //key is the index of the snarl and value is a vector of pair of
            // net_handle_t, NodeClusters

            net_handle_t snarl_handle = kv.first;

#ifdef DEBUG_CLUSTER
            cerr << "At depth " << depth << distance_index.net_handle_as_string(snarl_handle)
                << " headed by " << snarl_index.id_in_parent
                << " with children " << endl;
            for (auto it2 : kv.second) {
                cerr << "\t" << typeToString(it2.first.node_type)
                     << " number " << it2.first.node_id << endl;
            }
#endif
            if (!distance_index.is_root(distance_index.parent(snarl_handle))){
                //If this snarl is in a chain, cluster and add let the
                //tree state know which chain it belongs to

                tree_state.parent_chain_to_children.at(distance_index.parent(snarl_handle)).emplace_back(
                        cluster_one_snarl(tree_state, snarl_handle));

#ifdef DEBUG_CLUSTER
                cerr << "Recording snarl number " << snarl_i << " headed by "
                      << snarl_index.id_in_parent  << " as a child of chain number "
                      << chain_assignment << " headed by " << snarl_index.parent_id << endl;
#endif

            } 
        }
    }

    void SnarlSeedClusterer::cluster_chain_level(TreeState& tree_state, size_t depth) const {
        for (auto& kv : tree_state.chain_to_children) {
            //For each chain at this level that has relevant child snarls in it,
            //find the clusters.


            // Get the chain's number
            net_handle_t chain_handle = kv.first;

#ifdef DEBUG_CLUSTER
            cerr << "At depth " << depth <<" " <<  distance_index.net_handle_as_string(chain_handle)
                 << " with children " << endl;
            for (auto it2 : kv.second) {
                cerr << "\t snarl number " << it2.second.first << endl;
            }
#endif

            // Compute the clusters for the chain
            if (distance_index.is_root(distance_index.get_parent(chain_handle))) {
                cluster_one_chain(tree_state, chain_handle, depth);
            } else {
                tree_state.snarl_to_children.at(distance_index.get_parent(chain_handle)).emplace_back(
                    cluster_one_chain(tree_state, chain_i, depth));
#ifdef DEBUG_CLUSTER
                cerr << "Recording chain number " << chain_i << " headed by "
                     << dist_index.chain_indexes[chain_i].id_in_parent
                    << " as a child of snarl number " << parent_snarl_i
                    << " headed by " << parent_id << endl;
                cerr << "Snarl number " << parent_snarl_i << " has "
                    << tree_state.parent_snarl_to_nodes[parent_snarl_i].size()
                    << " children now" << endl;
#endif
            }
        }
    }


    SnarlSeedClusterer::NodeClusters SnarlSeedClusterer::cluster_one_node(
                       TreeState& tree_state,
                       id_t node_id, int64_t node_length) const {
#ifdef DEBUG_CLUSTER
        cerr << "Finding clusters on node " << node_id << " which has length " <<
        node_length << endl;
#endif

        //Final clusters on the node that we will be returning
        NodeClusters node_clusters(distance_index.get_node_net_handle(node_id), tree_state.all_seeds->size());

        if (tree_state.read_distance_limit > node_length) {
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
                        size_t& seed_i = iter->second
                        pos_t seed = tree_state.all_seeds->at(read_num)->at(seed_i).pos;
                        //And the distances for this seed
                        int64_t dist_left = is_rev(seed) ? node_length- get_offset(seed) : get_offset(seed) + 1;
                        int64_t dist_right = is_rev(seed) ? get_offset(seed) + 1 : node_length - get_offset(seed);

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
                    tree_state.read_cluster_dists[read_num][group_id] = make_pair(node_clusters.read_best_left[read_num],
                                                        node_clusters.read_best_right[read_num]);
                    node_clusters.read_cluster_heads.emplace(read_num, group_id);

                    if (tree_state.fragment_distance_limit != 0) {
                        fragment_group_id = tree_state.fragment_union_find.find_group(fragment_group_id);
                    }
                }
            }
#ifdef DEBUG_CLUSTER
            cerr << "Found single cluster on node " << node_id << " with fragment dists " 
                    << node_clusters.fragment_best_left << " " << node_clusters.fragment_best_right << endl;

            bool got_left = false;
            bool got_right = false;
            for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
                cerr << " for read num " << read_num << " best left: " << node_clusters.read_best_left[read_num] << " best right: " << node_clusters.read_best_right[read_num] << endl;
                bool got_read_left=false;
                bool got_read_right = false;
                for (pair<size_t,size_t> c : node_clusters.read_cluster_heads) {
                    if (c.first == read_num) {
                        pair<int64_t, int64_t> dists = tree_state.read_cluster_dists[c.first][c.second];
                        cerr << "\t" << c.first << ":"<<c.second << ": left: " << dists.first << " right : " << dists.second << ": ";
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
                assert(got_read_left || node_clusters.read_best_left[read_num] == std::numeric_limits<size_t>::max());
                assert(got_read_right || node_clusters.read_best_right[read_num] == std::numeric_limits<size_t>::max());
            }
            assert(got_left);
            assert(got_right);

#endif
            return node_clusters;
        }


        //The seeds may for multiple clusters on the node
        //Sort the seeds by their offset in the node and split into clusters
            
        //<index of read, index of seed, offset of seed> for all seeds
        vector<tuple<size_t,size_t, int64_t>> seed_offsets;
        for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
            auto seed_range_start = std::lower_bound(
                tree_state.node_to_seeds[read_num].begin(),tree_state.node_to_seeds[read_num].end(),
                std::pair<id_t, size_t>(node_id, 0));
            if (seed_range_start != tree_state.node_to_seeds[read_num].end() && seed_range_start->first == node_id) {
                for (auto iter = seed_range_start; iter != tree_state.node_to_seeds[read_num].end() 
                                                    && iter->first == node_id; ++iter) {
                    //For each seed, find its offset
                    pos_t seed = tree_state.all_seeds->at(read_num)->at(iter->second).pos;
                    int64_t offset = is_rev(seed) ? node_length - get_offset(seed) : get_offset(seed) + 1;

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

        vector<int64_t> read_last_offset (tree_state.all_seeds->size(), std::numeric_limits<size_t>::max());
        size_t fragment_last_offset = std::numeric_limits<size_t>::max();
        size_t fragment_last_cluster = std::numeric_limits<size_t>::max();
        vector<size_t> read_last_cluster (tree_state.all_seeds->size(), std::numeric_limits<size_t>::max());

        for ( tuple<size_t, size_t, int64_t> s : seed_offsets) {
            //For each seed, in order of position in the node,
            //see if it belongs to a new read/fragment cluster - if it is
            //close enough to the previous seed
            size_t read_num = std::get<0>(s);

            if (read_last_offset[read_num] != std::numeric_limits<size_t>::max() &&
                std::get<2>(s) - read_last_offset[read_num] <= tree_state.read_distance_limit) {
                //If this seed is in the same read cluster as the previous one,
                //union them

                int64_t prev_dist_left = tree_state.read_cluster_dists[read_num][read_last_cluster[read_num]].first;
                tree_state.read_union_find[read_num].union_groups(std::get<1>(s), read_last_cluster[read_num]);
                read_last_cluster[read_num] = tree_state.read_union_find[read_num].find_group(std::get<1>(s));
                tree_state.read_cluster_dists[read_num][read_last_cluster[read_num]] = 
                                       make_pair(prev_dist_left,node_length-std::get<2>(s)+1);
                read_last_offset[read_num] = std::get<2>(s);

                if (tree_state.fragment_distance_limit != 0) {
                    //If we are also clustering paired end reads by fragment distance,
                    //cluster these together
                    tree_state.fragment_union_find.union_groups(std::get<1>(s)+tree_state.read_index_offsets[read_num], fragment_last_cluster);
                    fragment_last_cluster = tree_state.fragment_union_find.find_group(std::get<1>(s)+tree_state.read_index_offsets[read_num]);
                    fragment_last_offset = std::get<2>(s);
                }
            } else {
                //This becomes a new read cluster
                if (read_last_cluster[read_num] != std::numeric_limits<size_t>::max()) {
                    //Record the previous cluster
                    node_clusters.read_cluster_heads.emplace(read_num, read_last_cluster[read_num]);
                }
                read_last_cluster[read_num] = std::get<1>(s);
                read_last_offset[read_num] = std::get<2>(s);
                tree_state.read_cluster_dists[read_num][read_last_cluster[read_num]] = 
                        make_pair(read_last_offset[read_num], node_length - read_last_offset[read_num] + 1);
                if (tree_state.fragment_distance_limit != 0) {
                    if (fragment_last_offset != std::numeric_limits<size_t>::max() &&
                        std::get<2>(s) - fragment_last_offset <= tree_state.fragment_distance_limit) {
                        //If this is a new read cluster but the same fragment cluster
                        tree_state.fragment_union_find.union_groups(std::get<1>(s)+tree_state.read_index_offsets[read_num], fragment_last_cluster);
                        fragment_last_cluster = tree_state.fragment_union_find.find_group(fragment_last_cluster);
                        fragment_last_offset = std::get<2>(s);

                    } else {
                        //If this is a new fragment cluster as well
                        fragment_last_cluster = std::get<1>(s)+tree_state.read_index_offsets[read_num];
                        fragment_last_offset = std::get<2>(s);
                    }
                }
            }
        }
        for (size_t i = 0 ; i < read_last_cluster.size() ; i++) {
            if (read_last_cluster[i] != std::numeric_limits<size_t>::max()) {
                node_clusters.read_cluster_heads.emplace(i, read_last_cluster[i]);
            }
        }

#ifdef DEBUG_CLUSTER

        cerr << "Found read clusters on node " << node_id << endl;

        bool got_left = false;
        bool got_right = false;
        for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
            cerr << " for read num " << read_num << " best left: " << node_clusters.read_best_left[read_num] << " best right: " << node_clusters.read_best_right[read_num] << endl;
            bool got_read_left=false;
            bool got_read_right = false;
            for (pair<size_t,size_t> c : node_clusters.read_cluster_heads) {
                if (c.first == read_num) {
                    pair<int64_t, int64_t> dists = tree_state.read_cluster_dists[c.first][c.second];
                    cerr << "\t" << c.first << ":"<<c.second << ": left: " << dists.first << " right : " << dists.second << ": ";
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
            assert(got_read_left || node_clusters.read_best_left[read_num] == std::numeric_limits<size_t>::max());
            assert(got_read_right || node_clusters.read_best_right[read_num] == std::numeric_limits<size_t>::max());
        }
        assert(got_left);
        assert(got_right);
        for (pair<size_t, size_t> group_id : node_clusters.read_cluster_heads) {
            assert (group_id.second == tree_state.read_union_find[group_id.first].find_group(group_id.second));
        }
#endif
        
        return node_clusters;

    };



    SnarlSeedClusterer::NodeClusters SnarlSeedClusterer::cluster_one_snarl(
                    TreeState& tree_state, net_handle_t snarl_handle) const {
        //Get the clusters on this snarl, assumes that all of the snarls children have been clustered already.
        
#ifdef DEBUG_CLUSTER
        cerr << "Finding clusters on snarl " << distance_index.net_handle_as_string(snarl_handle) << endl;
#endif

        //Keep track of all clusters on this snarl
        NodeClusters snarl_clusters(snarl_handle, tree_state.all_seeds->size());


        auto combine_clusters = [&] (size_t& new_group, size_t& combined_group, size_t& fragment_combined_group,
                                    int64_t fragment_dist, int64_t read_dist, size_t read_num){
            //Helper function to compare and combine clusters in two nodes of the same snarl
            //If the distance between two clusters is small enough, then combine them
            //for the read clusters and, if applicable, for the fragment clusters
            //Updates the distances stored for the read clusters
            if (read_dist != std::numeric_limits<size_t>::max() && read_dist <= tree_state.read_distance_limit) {
                //If the clusters are close enough to combine in the read
                if (tree_state.fragment_distance_limit != 0) {
                    if (fragment_combined_group != std::numeric_limits<size_t>::max()) {
                        //Also combine fragment clusters
                        tree_state.fragment_union_find.union_groups(new_group+tree_state.read_index_offsets[read_num], 
                                                                    fragment_combined_group);
                    }
                    fragment_combined_group = tree_state.fragment_union_find.find_group(new_group+tree_state.read_index_offsets[read_num]);
                }
                pair<int64_t, int64_t>& end_dists = tree_state.read_cluster_dists[read_num][new_group];

                if (combined_group == std::numeric_limits<size_t>::max()) {
                    snarl_clusters.read_cluster_heads.emplace(read_num,new_group);
                    tree_state.read_cluster_dists[read_num][new_group] = end_dists;
                    combined_group = new_group;
                } else {
                    //Combine the clusters within the same read

                    combined_group = tree_state.read_union_find[read_num].find_group(combined_group);
                    pair<int64_t, int64_t>old_dists = tree_state.read_cluster_dists[read_num][combined_group];
                    tree_state.read_union_find[read_num].union_groups(new_group, combined_group);

                    //Update distances and cluster head of new cluster
                    size_t new_g = tree_state.read_union_find[read_num].find_group(new_group);
                    if (new_g != new_group) snarl_clusters.read_cluster_heads.erase(make_pair(read_num,new_group));
                    if (new_g != combined_group) snarl_clusters.read_cluster_heads.erase(make_pair(read_num,combined_group));
                    
                    snarl_clusters.read_cluster_heads.emplace(read_num,new_g);
                    end_dists = make_pair( std::min(end_dists.first, old_dists.first),
                                           std::min(end_dists.second, old_dists.second));
                    tree_state.read_cluster_dists[read_num][new_g] = end_dists;
                    new_group = new_g;
                    combined_group = new_g;
                }

            } else if (tree_state.fragment_distance_limit != 0
                  && fragment_dist <= tree_state.fragment_distance_limit) {

                //Same fragment
                if (fragment_combined_group != std::numeric_limits<size_t>::max()) {
                    //Also combine fragment clusters
                    tree_state.fragment_union_find.union_groups(new_group+tree_state.read_index_offsets[read_num], 
                                                                fragment_combined_group);
                }
                fragment_combined_group = tree_state.fragment_union_find.find_group(new_group+tree_state.read_index_offsets[read_num]);
            }
            return;
        };


        //Get the children of this snarl and their clusters
        vector< NodeClusters>& children = tree_state.snarl_to_children[snarl_handle];
        int64_t start_length = snarl_index.node_length(0);
        int64_t end_length = snarl_index.node_length(snarl_index.num_nodes*2 -1);


        //Maps each child to its old left and right distances
        hash_map<net_handle_t, pair<int64_t, int64_t>> old_dists;
        old_dists.reserve(children.size());

        for (size_t i = 0; i < children.size() ; i++) {
            //Go through each child node of the netgraph

            NodeClusters& curr_child_clusters = child_nodes[i];
            net_handle_t& curr_child_handle = curr_child_handle.containing_net_handle;

#ifdef DEBUG_CLUSTER
            cerr << "Finding distances to parent snarl " << snarl_index_i
                 << " ends from child " << i << "/" << child_nodes.size() << endl;
            cerr << "Child is " << typeToString(child.node_type) << " number "
                 << child.node_id << " headed by " << child_node_id << endl;
            cerr << "Node rank is " << node_rank << " fwd, " << rev_rank
                 << " rev of " << snarl_index.num_nodes * 2 << endl;
            cerr << "Clusters at this child:" << endl;
            for (pair<size_t, size_t> c : curr_child_clusters.read_cluster_heads) {
                cerr << "\tdist left: " << tree_state.read_cluster_dists[c.first][c.second].first
                << " dist right: " << tree_state.read_cluster_dists[c.first][c.second].second << endl;
                cerr << "\t\t";
                for (size_t x = 0 ; x < tree_state.all_seeds->at(c.first)->size() ; x++) {
                    if (tree_state.read_union_find[c.first].find_group(x) == c.second) {
                        cerr << tree_state.all_seeds->at(c.first)->at(x).pos << " ";
                    }
                }
                cerr << endl;
            }

#endif

            vector<pair<size_t, size_t>> children_i(
                  make_move_iterator(curr_child_clusters.read_cluster_heads.begin()),
                  make_move_iterator(curr_child_clusters.read_cluster_heads.end()));
            for (size_t c_i = 0 ; c_i < children_i.size() ; c_i ++) {
                //for each cluster of child node i, find the distances to the
                //ends of the snarl

                //Cluster head of this cluster (read_num, index of the seed)
                size_t read_num = children_i[c_i].first;
                size_t seed_num = children_i[c_i].second;
                pair<size_t, size_t> child_cluster_head = children_i[c_i];

                //Distances of the cluster on the child
                pair<int64_t, int64_t> dists_c = tree_state.read_cluster_dists[read_num][seed_num];
                old_dists[make_pair(read_num, seed_num)] = dists_c;

                //Update the distances to the ends of the current snarl 
                pair<int64_t, int64_t> new_dists = snarl_index.dist_to_ends(node_rank,
                                        dists_c.first,dists_c.second);
#ifdef DEBUG_CLUSTER
cerr << "\tcluster: " << c_i << "dists to ends in snarl" << snarl_index.id_in_parent
     << " : " << new_dists.first << " " << new_dists.second << endl;
#endif

                snarl_clusters.fragment_best_left =std::min( snarl_clusters.fragment_best_left,new_dists.first);
                snarl_clusters.fragment_best_right = std::min(snarl_clusters.fragment_best_right, new_dists.second);
                snarl_clusters.read_best_left[read_num] =std::min(
                                   snarl_clusters.read_best_left[read_num], new_dists.first);
                snarl_clusters.read_best_right[read_num] = std::min(
                                   snarl_clusters.read_best_right[read_num], new_dists.second);


                snarl_clusters.read_cluster_heads.insert(make_pair(read_num, seed_num));
                tree_state.read_cluster_dists[read_num][seed_num] = new_dists;
            }


            for (size_t j = 0 ; j <= i ; j++){
                //Go through other child net graph nodes up to and including i

                //Get the other node and its clusters
                NodeClusters& other_child_clusters = child_nodes[j];
                net_handle_t& other_child_handle = other_node_clusters.containing_net_handle;

#ifdef DEBUG_CLUSTER
                cerr << "Other net graph node is " << typeToString(other_node.node_type)
                    << " headed by node " << other_node_id;


#endif


                //Find distance from each end of current node (i) to
                //each end of other node (j)
                int64_t dist_l_l = distance_index.distance_in_parent(snarl_handle, distance_index.flip(curr_child_handle), other_child_handle);
                int64_t dist_l_r = distance_index.distance_in_parent(snarl_handle, distance_index.flip(curr_child_handle), distance_index.flip(other_child_handle));
                int64_t dist_r_l = distance_index.distance_in_parent(snarl_handle, curr_child_handle, other_child_handle);
                int64_t dist_r_r = distance_index.distance_in_parent(snarl_handle, curr_child_handle, distance_index.flip(other_child_handle));

#ifdef DEBUG_CLUSTER
cerr << "\t distances between ranks " << node_rank << " and " << other_rank
     << ": " << dist_l_l << " " << dist_l_r << " " << dist_r_l << " "
     << dist_r_r << endl;
#endif

                //group ids of clusters combined between node i left and
                //node j left, etc
                vector<size_t> group_l_l (tree_state.all_seeds->size(), std::numeric_limits<size_t>::max());
                vector<size_t> group_l_r (tree_state.all_seeds->size(), std::numeric_limits<size_t>::max());
                vector<size_t> group_r_l (tree_state.all_seeds->size(), std::numeric_limits<size_t>::max());
                vector<size_t> group_r_r (tree_state.all_seeds->size(), std::numeric_limits<size_t>::max());
                size_t fragment_group_l_l = std::numeric_limits<size_t>::max();
                size_t fragment_group_l_r = std::numeric_limits<size_t>::max();
                size_t fragment_group_r_l = std::numeric_limits<size_t>::max();
                size_t fragment_group_r_r = std::numeric_limits<size_t>::max();

                //Are we looking for fragment clusters?
                bool looking_for_fragment_clusters = tree_state.fragment_distance_limit != 0;
                //Are the two nodes reachable with each other within the read distance limit
                bool reachable_within_read_distance_limit = std::min(std::min(dist_l_l, dist_l_r), std::min(dist_r_l, dist_r_r))-2 <= tree_state.read_distance_limit;
                //And are they reachable within the fragment distance limit
                bool reachable_within_fragment_distance_limit = std::min(std::min(dist_l_l, dist_l_r), std::min(dist_r_l, dist_r_r))-2 <= tree_state.fragment_distance_limit;
                
                //TODO: I'm not sure about this one
                bool within_fragment_distance_limit = std::min(curr_child_clusters.fragment_best_left, curr_child_clusters.fragment_best_right)-2
                                                <= tree_state.read_distance_limit;

                if (max({dist_l_l, dist_l_r, dist_r_l, dist_r_r}) != std::numeric_limits<size_t>::max()                            //If they are reachable
                    && within_fragment_distance_limit                                              //and they reach the ends of the snarl  
                    && ((!looking_for_fragment_clusters && reachable_within_read_distance_limit )  //and they are within the limit we're looking for
                        ||
                        (looking_for_fragment_clusters && reachable_within_fragment_distance_limit)
                                                )) {

                    for (size_t c_i = 0 ; c_i < children_i.size() ; c_i ++) {
                        //for each cluster of child node i

                        pair<size_t, size_t> child_cluster_head = children_i[c_i];
                        size_t read_num = child_cluster_head.first;
                        size_t c_group = tree_state.read_union_find[read_num].find_group(child_cluster_head.second);
                        pair<int64_t, int64_t> dists_c = old_dists[child_cluster_head];


                        if (dist_l_l != std::numeric_limits<size_t>::max() && dists_c.first != std::numeric_limits<size_t>::max() && other_node_clusters.fragment_best_left != std::numeric_limits<size_t>::max() ) {
                            //If cluster child_cluster_head can be combined with clusters in j
                            //from the left of both of them
                            int64_t read_dist = other_node_clusters.read_best_left[read_num] == std::numeric_limits<size_t>::max() ? std::numeric_limits<size_t>::max() :
                                  dist_l_l + dists_c.first + other_node_clusters.read_best_left[read_num]-1;
                            int64_t fragment_dist = dist_l_l + dists_c.first + other_node_clusters.fragment_best_left-1;
                            combine_clusters(c_group, group_l_l[read_num], fragment_group_l_l,
                                  fragment_dist, read_dist,  read_num);
                        }

                        if (dist_l_r != std::numeric_limits<size_t>::max() && dists_c.first != std::numeric_limits<size_t>::max() && other_node_clusters.fragment_best_right != std::numeric_limits<size_t>::max() ) {
                            //If it can be combined from the left to the right of j
                            int64_t fragment_dist = dist_l_r + dists_c.first + other_node_clusters.fragment_best_right-1;
                            int64_t read_dist = other_node_clusters.read_best_right[read_num] == std::numeric_limits<size_t>::max() ? std::numeric_limits<size_t>::max() :
                                 dist_l_r + dists_c.first + other_node_clusters.read_best_right[read_num]-1;
                            combine_clusters(c_group, group_l_r[read_num], fragment_group_l_r,
                                 fragment_dist, read_dist, read_num);
                        }
                        if (dist_r_l != std::numeric_limits<size_t>::max() && dists_c.second != std::numeric_limits<size_t>::max() && other_node_clusters.fragment_best_left != std::numeric_limits<size_t>::max() ) {
                            int64_t fragment_dist = dist_r_l + dists_c.second + other_node_clusters.fragment_best_left-1;
                            int64_t read_dist = other_node_clusters.read_best_left[read_num] == std::numeric_limits<size_t>::max() ? std::numeric_limits<size_t>::max() :
                                dist_r_l + dists_c.second + other_node_clusters.read_best_left[read_num]-1;
                            combine_clusters(c_group, group_r_l[read_num], fragment_group_r_l,
                                fragment_dist, read_dist,  read_num);
                        }
                        if (dist_r_r != std::numeric_limits<size_t>::max() && dists_c.second != std::numeric_limits<size_t>::max() && other_node_clusters.fragment_best_right != std::numeric_limits<size_t>::max() ) {
                            int64_t fragment_dist = dist_r_r + dists_c.second + other_node_clusters.fragment_best_right-1;
                            int64_t read_dist = other_node_clusters.read_best_right[read_num] == std::numeric_limits<size_t>::max() ? std::numeric_limits<size_t>::max() :
                                dist_r_r + dists_c.second + other_node_clusters.read_best_right[read_num]-1;
                            combine_clusters(c_group, group_r_r[read_num], fragment_group_r_r,
                                fragment_dist, read_dist, read_num);
                        }

                    }

                    //Go through clusters of child node j
                    vector<pair<size_t, size_t>> children_j(
                             make_move_iterator(other_node_clusters.read_cluster_heads.begin()),
                             make_move_iterator(other_node_clusters.read_cluster_heads.end()));

                    for (size_t k_i = 0 ; k_i < children_j.size() ; k_i++){
                        //For each cluster of child j, find which overlaps with clusters of i
                        //child_cluster_head will already be part of a cluster in snarl_cluster_heads but 
                        //since we need to know the node that the snarl is on we can't just loop through 
                        //snarl_cluster heads
                        pair<size_t,size_t> child_cluster_head = children_j[k_i];
                        size_t read_num = child_cluster_head.first;
                        pair<int64_t, int64_t>& dists_k = old_dists[child_cluster_head];
                        size_t k_group = tree_state.read_union_find[read_num].find_group(child_cluster_head.second);


                        if (dist_l_l != std::numeric_limits<size_t>::max() && curr_child_clusters.fragment_best_left != std::numeric_limits<size_t>::max() && dists_k.first != std::numeric_limits<size_t>::max() ){

                            int64_t fragment_dist = dist_l_l + curr_child_clusters.fragment_best_left + dists_k.first-1;
                            int64_t read_dist = curr_child_clusters.read_best_left[read_num] == std::numeric_limits<size_t>::max() ? std::numeric_limits<size_t>::max() :
                                dist_l_l + curr_child_clusters.read_best_left[read_num] + dists_k.first-1;
                            combine_clusters(k_group, group_l_l[read_num], fragment_group_l_l,
                                fragment_dist,read_dist, read_num);
                        }
                        if (dist_l_r != std::numeric_limits<size_t>::max() && curr_child_clusters.fragment_best_left != std::numeric_limits<size_t>::max() && dists_k.second != std::numeric_limits<size_t>::max()  ) {

                            int64_t read_dist = curr_child_clusters.read_best_left[read_num] == std::numeric_limits<size_t>::max() ? std::numeric_limits<size_t>::max() :
                               dist_l_r + curr_child_clusters.read_best_left[read_num] + dists_k.second-1;
                            int64_t fragment_dist = dist_l_r + curr_child_clusters.fragment_best_left + dists_k.second-1;
                            combine_clusters(k_group, group_l_r[read_num], fragment_group_l_r,
                               fragment_dist, read_dist, read_num);
                        }
                        if (dist_r_l != std::numeric_limits<size_t>::max() && curr_child_clusters.fragment_best_right != std::numeric_limits<size_t>::max() && dists_k.first != std::numeric_limits<size_t>::max()  ) {

                            int64_t fragment_dist = dist_r_l + curr_child_clusters.fragment_best_right + dists_k.first-1;
                            int64_t read_dist = curr_child_clusters.read_best_right[read_num] == std::numeric_limits<size_t>::max() ? std::numeric_limits<size_t>::max() :
                                dist_r_l + curr_child_clusters.read_best_right[read_num] + dists_k.first-1;
                            combine_clusters(k_group, group_r_l[read_num], fragment_group_r_l,
                                fragment_dist, read_dist, read_num);
                        }
                        if (dist_r_r != std::numeric_limits<size_t>::max() && curr_child_clusters.fragment_best_right != std::numeric_limits<size_t>::max() && dists_k.second != std::numeric_limits<size_t>::max() ) {

                            int64_t fragment_dist = dist_r_r + curr_child_clusters.fragment_best_right + dists_k.second-1;
                            int64_t read_dist = curr_child_clusters.read_best_right[read_num] == std::numeric_limits<size_t>::max() ? std::numeric_limits<size_t>::max() :
                               dist_r_r + curr_child_clusters.read_best_right[read_num] + dists_k.second-1;
                            combine_clusters(k_group, group_r_r[read_num], fragment_group_r_r, fragment_dist, read_dist, read_num);
                        }
                    }
                }
            }
        }
#ifdef DEBUG_CLUSTER
        cerr << "Found clusters on snarl number " << snarl_index_i << " headed by"
             << snarl_index.id_in_parent << endl;
        cerr << "    with best left and right values: " << snarl_clusters.fragment_best_left << " "
             << snarl_clusters.fragment_best_right << endl;
        bool got_left = false;
        bool got_right = false;
        for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
            cerr << " for read num " << read_num << " best left: " << snarl_clusters.read_best_left[read_num] << " best right: " << snarl_clusters.read_best_right[read_num] << endl;
            bool got_read_left=false;
            bool got_read_right = false;
            bool any_clusters = false;
            for (pair<size_t,size_t> c : snarl_clusters.read_cluster_heads) {
                if (c.first == read_num) {
                    any_clusters = true;
                    pair<int64_t, int64_t> dists = tree_state.read_cluster_dists[c.first][c.second];
                    cerr << "\t" << c.first << ":"<<c.second << ": left: " << dists.first << " right : " << dists.second << ": ";
                    bool has_seeds = false;
                    for (size_t x = 0 ; x < tree_state.all_seeds->at(c.first)->size() ; x++) {
                        if (tree_state.read_union_find[c.first].find_group(x) == c.second) {
                            cerr << tree_state.all_seeds->at(c.first)->at(x).pos << " ";
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

        for (pair<size_t, size_t> group_id : snarl_clusters.read_cluster_heads) {
            assert (group_id.second == tree_state.read_union_find[group_id.first].find_group(group_id.second));
        }
#endif
        return snarl_clusters;
    };



    SnarlSeedClusterer::NodeClusters SnarlSeedClusterer::cluster_one_chain(TreeState& tree_state, net_handle_t chain_handle, size_t depth) const {


#ifdef DEBUG_CLUSTER
        cerr << "Cluster chain " << chain_index.id_in_parent << endl;
#endif

        //A node cluster for each child of the chain, in the order of the snarls in the chain
        //Maps the offset in the chain records to the node cluster (offset doesn't matter for anything other than putting them in order)
        //TODO: I think this doesn't need to be a map anymore since we're sorting later anyway, depending on the sort algorithm
        hash_map<size_t, NodeClusters>& children_in_chain =tree_state.chain_to_children[chain_handle];


        /*Start by making a list of all clusters (snarl and node) ordered by their occurrence in the chain
         * Snarl clusters are only ordered relative to seed clusters
         * When we add the snarl clusters to cluster_head_indices, check if clusters in the same snarl can combine
         * and update each one's distance to the ends of the chain
         */

        //Holds all clusters as tuple of <the snarl/node the cluster is on, read num, seed num> 
        //vectors of tuples are sorted by the offset in the chain, so in the order they occur in the chain
        vector<tuple<net_handle_t, size_t, size_t>> cluster_head_indices ;
        cluster_head_indices.reserve(children_in_chain.size());


        //Helper function to insert a cluster (offset in chain, read num, seed num) into a sorted list of clusters
        auto insert_in_order = [&](vector<tuple<net_handle_t, size_t, size_t>>& input_vector, tuple<net_handle_t, size_t, size_t> item) {
            //Get the offset of the thing we're inserting
            const Seed& curr_seed = tree_state.all_seeds->at(std::get<1>(item))->at(std::get<2>(item));
            pos_t pos =  curr_seed.pos;
            net_handle_t& handle = std::get<0>(item);

            //The offset of the record in the chain
            int64_t item_offset = distance_index.rank_in_parent(handle);

            //Get an iterator to where the thing we're inserting should go in the vector
            std::vector<tuple<int64_t, size_t, size_t>>::iterator insert_itr = std::upper_bound(input_vector.begin(), input_vector.end(),
               item_offset, [&](int64_t val, tuple<size_t, size_t, size_t> vector_item) {
                   return val < distance_index.rank_in_parent(std::get<0>(vector_item)); 
               });
            //Add the item to the vector in sorted order
            input_vector.insert(insert_itr, item);
        };


        auto combine_snarl_clusters = [&] (size_t& new_group,
                        size_t& combined_group, size_t& fragment_combined_group,
                        vector<pair<size_t,size_t>>& to_erase, int64_t fragment_dist,int64_t read_dist,
                        pair<int64_t, int64_t>& dists, size_t read_num){
            //Helper function to combine clusters of the same snarl
            //Used when two clusters in the same snarl can be combined by
            //looping in the chain

            if (read_dist != std::numeric_limits<size_t>::max() && read_dist <= tree_state.read_distance_limit) {
                if (combined_group == std::numeric_limits<size_t>::max()) {
                    combined_group = new_group;
                } else {
                    //Union the two groups
                    tree_state.read_union_find[read_num].union_groups(combined_group, new_group);
                    //Find the new distances of the combined groups
                    pair<int64_t, int64_t>& old_dists = tree_state.read_cluster_dists[read_num][combined_group];
                    size_t new_combined_group = tree_state.read_union_find[read_num].find_group(new_group);
                    //Update which groups are being kept track of
                    if (new_combined_group != new_group) {
                        to_erase.emplace_back(read_num, new_group);
                    }
                    if (new_combined_group != combined_group)  {
                        to_erase.emplace_back(read_num, combined_group);
                    }
                    combined_group = new_combined_group;

                    dists = make_pair(
                          std::min(old_dists.first, dists.first),
                          std::min(old_dists.second, dists.second));
                    tree_state.read_cluster_dists[read_num][combined_group] = dists;
#ifdef DEBUG_CLUSTER
                    cerr << " New dists for read num " << read_num << ": "
                         << tree_state.read_cluster_dists[read_num][combined_group].first << " "
                         << tree_state.read_cluster_dists[read_num][combined_group].second << endl;
#endif
                }

                if (tree_state.fragment_distance_limit != 0 && fragment_dist != std::numeric_limits<size_t>::max()) {
                    if (fragment_combined_group != std::numeric_limits<size_t>::max()) {
                    //If we're also keeping track of fragment clusters
                        tree_state.fragment_union_find.union_groups(fragment_combined_group,
                                                            new_group + tree_state.read_index_offsets[read_num]);
                    }
                    fragment_combined_group = tree_state.fragment_union_find.find_group(
                                                         new_group + tree_state.read_index_offsets[read_num]);
                }
            } else if (tree_state.fragment_distance_limit != 0 && fragment_dist != std::numeric_limits<size_t>::max() &&
                        fragment_dist <= tree_state.fragment_distance_limit) {
                //If these aren't in the same read cluster but are in
                //the same fragment cluster
                if (fragment_combined_group != std::numeric_limits<size_t>::max()) {
                    tree_state.fragment_union_find.union_groups(
                                    fragment_combined_group, new_group + tree_state.read_index_offsets[read_num]);
                }
                fragment_combined_group = tree_state.fragment_union_find.find_group(
                                                     new_group + tree_state.read_index_offsets[read_num]);
            }
            return;
        };

        //TODO size_t connected_component_num = dist_index.get_connected_component(chain_index.id_in_parent);
        //Get the length of the chain and the length of its boundary nodes
        size_t chain_length = distance_index.get_min_length(chain_handle);
        size_t chain_start_length = distance_index.minimum_length(distance_index.get_bound(chain_handle, false, false));
        size_t chain_end_length = distance_index.minimum_length(distance_index.get_bound(chain_handle, true, false));


        //Go through the chain by children. 
        //
        //For each child, 
        // - check if the clusters in the child can be combined with each other by walking out then back in through the chain
        // - update the distances to the ends of the chain
        // - compare and combine with clusters of the chain that get build as we walk along it
        //      Each child needs to be compared everything found since the last node, since if it could reach any seed on 
        //      the chain to the left of the last node found, then the last node would also be able to reach that seed so
        //      it would be combined anyway
        //
        // - after combining clusters of the current child, remove redundant cluster heads from the chain clusters
        //
        // Note that every snarl cluster must be compared to the current child, even if there was a snarl cluster that came
        // after it in the chain (but not if there was a node cluster after it). If there were clusters on snarls 1, 2, and 4
        // in the chain, and 1 got combined with 2, the minimum distance from 1 to 4 might still be smaller than the minimum 
        // distance from 2 to 4, so comparing only 2 to 4 might miss the distance from 1 to 4.
        // If multiple clusters in the same snarl get combined, then they the redundant cluster head can be removed
        //


        // As we walk along the chain, we need to remember all clusters since the last node
        // This stores all clusters we found since the last node as a pair< net_handle_t to the child it's on, pair<read_num, seed_num>>
        // It'll store all seeds, regardless of which read, and
        for (auto& kv: children_in_chain) {
            /*
             * Snarls and nodes are in the order that they are traversed in the chain
             * children_in_chain is a map from rank of a snarl/node to the NodeClusters for that snarl/node
             */

            //The clusters of the current child
            NodeClusters& child_clusters = kv.second;
            //And a net handle to the current child
            net_handle_t child_handle = child_clusters.containing_net_handle;

/*TODO: I don't think I need these but I might
            int64_t start_length = snarl_index.rev_in_parent
                         ? snarl_index.node_length(snarl_index.num_nodes * 2 - 1)
                         : snarl_index.node_length(0);
            int64_t end_length = snarl_index.rev_in_parent ? snarl_index.node_length(0)
                        : snarl_index.node_length(snarl_index.num_nodes * 2 - 1);
            int64_t snarl_length = snarl_index.snarl_length();
*/


            //Get the loop distances (distance leaving the child and coming back to the same side) so we can
            //combine clusters on the same child that can be reached by looping
            int64_t loop_dist_end = snarl_distance_index.distance_in_parent(chain_handle, child_handle, distance_index.flip(child_handle));
            int64_t loop_dist_start = snarl_distance_index.distance_in_parent(chain_handle, distance_index.flip(child_handle), child_handle);


            //Distance from the ends of the current child to the ends of the chain
            int64_t add_dist_left_left = snarl_distance_index.distance_in_parent(chain_handle, snarl_distance_index.get_bound(chain_handle, false, true), child_handle) + chain_start_length;

            int64_t add_dist_right_right = snarl_distance_index.distance_in_parent(chain_handle, snarl_distance_index.get_bound(chain_handle, true, true), distance_index.flip(child_handle)) + chain_end_length;

            int64_t add_dist_left_right = snarl_distance_index.distance_in_parent(chain_handle, snarl_distance_index.get_bound(chain_handle, false, true), distance_index.flip(child_handle)) + chain_start_length;
            int64_t add_dist_right_left = snarl_distance_index.distance_in_parent(chain_handle, snarl_distance_index.get_bound(chain_handle, true, true), child_handle) + chain_end_length;

            hash_set<pair<size_t,size_t>> to_add;//new cluster group ids from snarl clusters
            vector<pair<size_t,size_t>> to_erase; //old cluster group ids

            //Combined snarl clusters by taking chain loop left/right
            vector<size_t> snarl_cluster_left (tree_state.all_seeds->size(),std::numeric_limits<size_t>::max());
            vector<size_t> snarl_cluster_right (tree_state.all_seeds->size(), std::numeric_limits<size_t>::max());
            size_t fragment_snarl_cluster_left = std::numeric_limits<size_t>::max();
            size_t fragment_snarl_cluster_right = std::numeric_limits<size_t>::max();

            for (pair<size_t, size_t> cluster_head : child_clusters.read_cluster_heads) {
                // For each of the clusters for the current child,
                // first check if it can be combined with another cluster
                // in the same snarl by taking loops in the chain,
                // then, find if it belongs to the new combined cluster
                // that includes chain clusters
                size_t read_num = cluster_head.first;

                pair<int64_t, int64_t> snarl_dists = tree_state.read_cluster_dists[read_num][cluster_head.second];

                if (loop_dist_start != std::numeric_limits<size_t>::max()) {
                    //If there is a loop going out and back into the start of
                    //the snarl, this cluster may be combined with other snarl clusters

                    //The distance to the right side of the snarl that is found by taking the leftmost seed and
                    // looping through the chain to the left
                    int64_t new_right = snarl_dists.first == std::numeric_limits<size_t>::max() || loop_dist_start == std::numeric_limits<size_t>::max()
                                        ? std::numeric_limits<size_t>::max()
                                        : snarl_dists.first + loop_dist_start + snarl_length;
                    if (snarl_dists.second == std::numeric_limits<size_t>::max() || (new_right != std::numeric_limits<size_t>::max() && new_right < snarl_dists.second)){
                        //If the best distance to the right can be update
                        snarl_dists.second = std::min(new_right, snarl_dists.second);
                        child_clusters.fragment_best_right = std::min(child_clusters.fragment_best_right, new_right);
                        child_clusters.read_best_right[read_num] = std::min(child_clusters.read_best_right[read_num], new_right);
#ifdef DEBUG_CLUSTER
cerr << "  Updating looping distance to right of snarl cluster " << read_num <<":" << cluster_head.second << ": "
     << new_right << " -> " << snarl_dists.second <<  endl;
#endif
                    }


                    if (child_clusters.fragment_best_left!= std::numeric_limits<size_t>::max() && snarl_dists.first != std::numeric_limits<size_t>::max() ) {
                        //If this cluster can be combined with another cluster
                        //from the left

#ifdef DEBUG_CLUSTER
cerr << "  Combining this cluster from the left " << endl;
#endif
                        int64_t read_dist =  child_clusters.read_best_left[read_num] == std::numeric_limits<size_t>::max() ? std::numeric_limits<size_t>::max() :  
                                    child_clusters.read_best_left[read_num] + snarl_dists.first + loop_dist_start - 1;
                        int64_t fragment_dist = child_clusters.fragment_best_left == std::numeric_limits<size_t>::max() ? std::numeric_limits<size_t>::max() :
                                    child_clusters.fragment_best_left + snarl_dists.first + loop_dist_start - 1;

                        combine_snarl_clusters(cluster_head.second, snarl_cluster_left[read_num], 
                                fragment_snarl_cluster_left,  to_erase, fragment_dist, read_dist, snarl_dists, read_num);
                    }

                }

                if (loop_dist_end != std::numeric_limits<size_t>::max()) {
                    //If there is a loop to the right
                    int64_t new_left = snarl_dists.second == std::numeric_limits<size_t>::max() || loop_dist_end == std::numeric_limits<size_t>::max() 
                                 ? std::numeric_limits<size_t>::max()
                                 : snarl_dists.second + loop_dist_end + snarl_length;
                    if (snarl_dists.first == std::numeric_limits<size_t>::max() || (new_left != std::numeric_limits<size_t>::max() && new_left < snarl_dists.first)){
                        //If this is an improvement, update distances
                        snarl_dists.first = new_left;
                        child_clusters.read_best_left[read_num] = std::min(new_left, child_clusters.read_best_left[read_num]);
                        child_clusters.fragment_best_left = std::min(new_left, child_clusters.fragment_best_left);

#ifdef DEBUG_CLUSTER
cerr << "Updating looping distance to left of snarl cluster " << read_num << ":" << cluster_head.second << ": "
     << new_left << endl;
#endif
                    }

                    if (child_clusters.fragment_best_right != std::numeric_limits<size_t>::max() && snarl_dists.second != std::numeric_limits<size_t>::max()) {
                        //If this cluster can be combined with another cluster
                        //from the right

#ifdef DEBUG_CLUSTER
cerr << "  Maybe combining this cluster from the right" << endl;
#endif
                        int64_t read_dist = child_clusters.read_best_right[read_num] == std::numeric_limits<size_t>::max() ? std::numeric_limits<size_t>::max() :
                            child_clusters.read_best_right[read_num] + snarl_dists.second  + loop_dist_end - 1;
                        int64_t fragment_dist = child_clusters.fragment_best_right == std::numeric_limits<size_t>::max() ? std::numeric_limits<size_t>::max() : 
                            child_clusters.fragment_best_right + snarl_dists.second + loop_dist_end - 1;

                        combine_snarl_clusters(cluster_head.second, snarl_cluster_right[read_num],
                             fragment_snarl_cluster_right, to_erase,fragment_dist, read_dist, snarl_dists, read_num);
                    }
                }

                to_add.emplace(cluster_head);
            }

            for (auto& cluster_head : to_erase) {
                to_add.erase(cluster_head);
            }

            //The offset of the last nucleotide in the start node of this snarl, so we can insert this snarl into
            //cluster_head_indices immediately after top-level seed clusters on the start node
            for (auto& cluster_head : to_add) {
                //Add the clusters on this snarl to our overall list of clusters and update the distances
                //for each of the clusters
                insert_in_order(cluster_head_indices, make_tuple(distance_index.get_rank_in_parent(child_handle), cluster_head.first, cluster_head.second));
                size_t cluster_distance_left = tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].first;
                size_t cluster_distance_right = tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].second;

                //Get the distance to the start side of the chain
                int64_t dist_left_left = cluster_distance_left  == std::numeric_limits<size_t>::max() ? std::numeric_limits<size_t>::max() 
                        : cluster_distance_left +  add_dist_left_left;
                int64_t dist_right_left = cluster_distance_right == std::numeric_limits<size_t>::max() 
                        || add_dist_right_left == std::numeric_limits<size_t>::max() ? std::numeric_limits<size_t>::max() 
                        : cluster_distance_right +  add_dist_right_left;


                int64_t dist_right_right = cluster_distance_right == std::numeric_limits<size_t>::max() ? std::numeric_limits<size_t>::max() 
                        : cluster_distance_right + add_dist_right_right;
                int64_t dist_left_right = cluster_distance_left == std::numeric_limits<size_t>::max() ||
                        add_dist_left_right == std::numeric_limits<size_t>::max() ? std::numeric_limits<size_t>::max() 
                        : cluster_distance_left + add_dist_left_right;

                tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].first  = std::min(dist_left_left, dist_right_left); 
                tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].second = std::min(dist_right_right, dist_left_right); 

            }

        }

        //Now that we updated the clusters and distances within the child snarl or node, we can cluster with the 
        //clusters on the chain that we've found so far



        //Keep track the snarl clusters from the current snarl and all previous snarl clusters
        //since the last time we saw a node
        //(We don't necessarily want to combine two clusters if they aren't top level seeds, since they would have
        //already been clustered if they were reachable, whereas nodes are always reachable)
        
        //cluster head, distance left, distance right
        vector<tuple<size_t, int64_t, int64_t>> prev_snarl_cluster_fragment;

        vector<tuple<size_t, int64_t, int64_t>> snarl_cluster_fragment;

        //The rightmost seed that is a top-level seed or in a cluster with a top-level seed
        //<cluster head, distance to the left of the chain
         pair<size_t, int64_t> seed_cluster_fragment (0, std::numeric_limits<size_t>::max());


        //And the same thing for each read
        vector<vector<tuple<size_t, int64_t, int64_t>>> snarl_cluster_by_read(tree_state.all_seeds->size());

        vector<vector<tuple<size_t, int64_t, int64_t>>> prev_snarl_cluster_by_read(tree_state.all_seeds->size());

        vector<pair<size_t, int64_t>> seed_cluster_by_read (tree_state.all_seeds->size(), make_pair(0, std::numeric_limits<size_t>::max()));

        //The offset of the last snarl we got clusters from
        //I"m using a root handle to denote that we haven't seen anything yet, since the root will never be a child
        vector<net_handle_t> last_snarl_by_read (tree_state.all_seeds->size(), distance_index.get_root());
        net_handle last_snarl_fragment = distance_index.get_root();


        //The clusters of the chain that are built from the snarl and minimizer clusters
        //This will get updated as we traverse through the child clusters
        NodeClusters chain_clusters(tree_state.all_seeds->size());

        //Go through the child clusters in order and combine them with chain clusters. If the cluster was a cluster on a snarl, then compare it to the
        //last seed cluster we found. If it was a seed, compare it to all previous clusters since we last saw a seed
        for (tuple<net_handle_t, size_t, size_t>& seed_index : cluster_head_indices) {

            //Cluster heads of the chain that are going to get deleted (because they got merged with another cluster and are not longer the heads)
            vector<pair<size_t, size_t>> to_erase;

            //The new cluster head for whatever this ends up clustering with
            size_t new_cluster_head = std::get<2>(seed_index);

            size_t read_num = std::get<1>(seed_index);
            if (distance_index.is_node(std::get<0>(seed_index) )) {

                //If this is a node, then try to compare it to the most recent snarl clusters and update the
                //best node-only clusters. Also compare it to the last node cluster

                //Move buffer for previous snarls
                prev_snarl_cluster_by_read[read_num].insert(prev_snarl_cluster_by_read[read_num].end(),
                    snarl_cluster_by_read[read_num].begin(), snarl_cluster_by_read[read_num].end());

                snarl_cluster_by_read[read_num].clear();

                int64_t offset = tree_state.all_seeds->at(read_num)->at(std::get<2>(seed_index)).offset;
                int64_t right_offset = chain_length - offset + 1;
                int64_t best_left = offset;
                int64_t best_right = right_offset;
#ifdef DEBUG_CLUSTER
                cerr << "At top-level seed cluster on read " << read_num << " with pos " << tree_state.all_seeds->at(read_num)->at(new_cluster_head).pos 
                     << " with dists " << offset << ", " << right_offset << endl;
                assert( tree_state.all_seeds->at(read_num)->at(new_cluster_head).component == connected_component_num); 

#endif

                //Update the chain's bes left and right distances
                chain_clusters.read_best_left[read_num] = std::min(chain_clusters.read_best_left[read_num], offset);
                chain_clusters.read_best_right[read_num] = std::min(chain_clusters.read_best_right[read_num], right_offset);
                chain_clusters.fragment_best_left = std::min( chain_clusters.fragment_best_left, offset);
                chain_clusters.fragment_best_right = std::min(chain_clusters.fragment_best_right, right_offset);

                //Remember that the last thing we saw was a node not a snarl
                last_snarl_fragment = distance_index.get_root();
                last_snarl_by_read[read_num] = distance_index.get_root();


                //try clustering by fragment
                if (tree_state.fragment_distance_limit != 0) {

                    prev_snarl_cluster_fragment.insert(prev_snarl_cluster_fragment.end(), 
                        snarl_cluster_fragment.begin(), snarl_cluster_fragment.end());

                    snarl_cluster_fragment.clear();

                    if (!prev_snarl_cluster_fragment.empty()) {
                        //Combine this seed cluster with all clusters from the most recent snarl we've seen since
#ifdef DEBUG_CLUSTER
                        cerr << "\tCompare to last snarl clusters fragment: " << endl;
#endif
                        for (size_t i = 0 ; i < prev_snarl_cluster_fragment.size() ; i++) {
                            int64_t distance_in_chain = std::get<2>(prev_snarl_cluster_fragment[i]) - right_offset;
                            if (std::get<2>(prev_snarl_cluster_fragment[i]) != std::numeric_limits<size_t>::max() && right_offset != std::numeric_limits<size_t>::max() && distance_in_chain <= tree_state.fragment_distance_limit) {
                                tree_state.fragment_union_find.union_groups(std::get<0>(prev_snarl_cluster_fragment[i]), 
                                                        std::get<2>(seed_index)+tree_state.read_index_offsets[read_num]);
#ifdef DEBUG_CLUSTER
                            cerr << "\t\tCombining fragment: top-level seed cluster on component " << connected_component_num << ", " << 
                                    tree_state.all_seeds->at(read_num)->at(std::get<2>(seed_index)).pos <<
                                    "(distances: " << offset << ", " << right_offset << ")" << 
                                    " with snarl cluster " << " " << read_num << " " << std::get<0>(prev_snarl_cluster_fragment[i]) << " " << 
                                    "(distances: " << std::get<1>(prev_snarl_cluster_fragment[i]) << ", " << std::get<2>(prev_snarl_cluster_fragment[i]) <<
                                    ") with distance between them " << distance_in_chain << endl; 
#endif
                                }
                            }
                    }

                    //Combine this seed cluster with the last seed cluster we found

                    int64_t distance_in_chain = offset - seed_cluster_fragment.second;

#ifdef DEBUG_CLUSTER
                    cerr << "\t Compare to last seed cluster fragment: " << endl;
                    //assert(distance_in_chain >= 0);
#endif

                    if (seed_cluster_fragment.second != std::numeric_limits<size_t>::max() && offset != std::numeric_limits<size_t>::max() && distance_in_chain <= tree_state.fragment_distance_limit) {
                        tree_state.fragment_union_find.union_groups(seed_cluster_fragment.first, 
                                                            std::get<2>(seed_index)+tree_state.read_index_offsets[read_num]);
#ifdef DEBUG_CLUSTER
                        cerr << "\t\tCombining fragment: top-level seed cluster on component " << connected_component_num << ", " << 
                                    tree_state.all_seeds->at(read_num)->at(std::get<2>(seed_index)).pos << 
                                    "(distances: " << offset << ", " << right_offset << ")" << 
                                    " with another seed cluster" << " " << read_num << " " << seed_cluster_fragment.first << " " << 
                                    "(distances: " << seed_cluster_fragment.second << ", " << chain_length-seed_cluster_fragment.second+1 << ") with distance between them " << distance_in_chain << endl; 
#endif
                    }
                    //Update the most recent fragment seed cluster
                    seed_cluster_fragment = make_pair(std::get<2>(seed_index)+tree_state.read_index_offsets[read_num], offset);

                    //Forget the last fragment snarl clusters we've seen
                    prev_snarl_cluster_fragment.clear();
                    snarl_cluster_fragment.clear();
                }


                //try clustering by read
                if (!prev_snarl_cluster_by_read[read_num].empty()) {
                    //Cluster with all snarl clusters we saw since last seed cluster and the most recent seed cluster
                    //
#ifdef DEBUG_CLUSTER
                        cerr << "\tCompare to last snarl clusters read" << endl;
#endif
                    for (size_t i = 0 ; i < prev_snarl_cluster_by_read[read_num].size() ; i++) {
                        int64_t distance_in_chain = std::get<2>(prev_snarl_cluster_by_read[read_num][i]) - right_offset;

#ifdef DEBUG_CLUSTER
                    //assert(distance_in_chain >= 0);
cerr <<  "\t\t cluster" << read_num << " " << std::get<0>(prev_snarl_cluster_by_read[read_num][i]) << " " << 
                                    tree_state.all_seeds->at(read_num)->at(std::get<0>(prev_snarl_cluster_by_read[read_num][i])).pos << 
                                    "(distances: " <<  std::get<1>(prev_snarl_cluster_by_read[read_num][i]) << ", " << std::get<2>(prev_snarl_cluster_by_read[read_num][i]) << ")";
#endif
                        if (std::get<2>(prev_snarl_cluster_by_read[read_num][i]) != std::numeric_limits<size_t>::max() && right_offset != std::numeric_limits<size_t>::max() && 
                                distance_in_chain <= tree_state.read_distance_limit) {
                            tree_state.read_union_find[read_num].union_groups(std::get<2>(seed_index), 
                                                                              std::get<0>(prev_snarl_cluster_by_read[read_num][i])); 
                            size_t new_group = tree_state.read_union_find[read_num].find_group(std::get<2>(seed_index));
                            if (new_group == new_cluster_head) {
                                to_erase.emplace_back(read_num, std::get<0>(prev_snarl_cluster_by_read[read_num][i]));
                            } else {
                                to_erase.emplace_back(read_num, new_cluster_head);
                            }
                            best_left = std::min(best_left, tree_state.read_cluster_dists[read_num][std::get<0>(prev_snarl_cluster_by_read[read_num][i])].first);
                            best_right = std::min(best_right, tree_state.read_cluster_dists[read_num][std::get<0>(prev_snarl_cluster_by_read[read_num][i])].second);
                            tree_state.read_cluster_dists[read_num][new_group] = make_pair(best_left, best_right);
                            std::get<0>(prev_snarl_cluster_by_read[read_num][i]) = new_group;
                            new_cluster_head = new_group;
#ifdef DEBUG_CLUSTER
                            cerr << "\t\tCombining read: top-level seed cluster on component " << connected_component_num << ", " << 
                                    tree_state.all_seeds->at(read_num)->at(std::get<2>(seed_index)).pos << 
                                    "(distances: " << offset << ", " << right_offset << ")" <<  
                                    " with snarl cluster with distance between them " << distance_in_chain << endl;
#endif
                        }
                    }
                }
#ifdef DEBUG_CLUSTER
                    cerr << "\t Compare to last seed cluster read: " << read_num << " " << seed_cluster_by_read[read_num].first << " " << 
                                    tree_state.all_seeds->at(read_num)->at(seed_cluster_by_read[read_num].first).pos << 
                                    "(distances: " << seed_cluster_by_read[read_num].second << ", " << chain_length - seed_cluster_by_read[read_num].second + 1 << 
                                    ")" << endl;
#endif
                //If the last seed we saw was a seed cluster, then compare it to that
                int64_t distance_in_chain = offset - seed_cluster_by_read[read_num].second;
#ifdef DEBUG_CLUSTER
                    //assert(distance_in_chain >= 0);
#endif
                if (seed_cluster_by_read[read_num].second != std::numeric_limits<size_t>::max() && offset != std::numeric_limits<size_t>::max() && 
                        distance_in_chain <= tree_state.read_distance_limit) {
                    tree_state.read_union_find[read_num].union_groups(std::get<2>(seed_index), 
                                                                          seed_cluster_by_read[read_num].first); 
                    size_t new_group = tree_state.read_union_find[read_num].find_group(std::get<2>(seed_index));
                    if (new_group == new_cluster_head) {
                        to_erase.emplace_back(read_num, seed_cluster_by_read[read_num].first);
                    } else {
                        to_erase.emplace_back(read_num, new_cluster_head);
                    }
                    best_left = std::min(best_left, seed_cluster_by_read[read_num].second);
                    best_right = std::min(best_right, chain_length-seed_cluster_by_read[read_num].second+1);
                    tree_state.read_cluster_dists[read_num][new_group] = make_pair(best_left, best_right );
                    seed_cluster_by_read[read_num].first = new_group;
                    new_cluster_head = new_group;
#ifdef DEBUG_CLUSTER
                    cerr << "\t\tCombining read: top-level seed cluster on component " << connected_component_num << ", " << 
                                tree_state.all_seeds->at(read_num)->at(std::get<2>(seed_index)).pos << 
                                "(distances: " << offset << ", " << right_offset << ")" << 
                                " with another seed cluster" << " with distance between them " << distance_in_chain << endl;
#endif
                }

                //Update the most recent seed cluster we've seen
                seed_cluster_by_read[read_num] = make_pair(std::get<2>(seed_index), offset);

                //And forget the snarl clusters we've seen since the last seed cluster
                snarl_cluster_by_read[read_num].clear();
                prev_snarl_cluster_by_read[read_num].clear();


            } else {
                //If this is an snarl cluster, compare it to the most recent seed-only cluster 
                //and all snarl clusters found after the most recent seed
                //and update the snarl clusters
                int64_t offset = tree_state.read_cluster_dists[read_num][std::get<2>(seed_index)].first;
                int64_t right_offset = tree_state.read_cluster_dists[read_num][std::get<2>(seed_index)].second;

                int64_t best_left = offset;
                int64_t best_right = right_offset;

                chain_clusters.read_best_left[read_num] = std::min(chain_clusters.read_best_left[read_num], offset);
                chain_clusters.read_best_right[read_num] = std::min(chain_clusters.read_best_right[read_num], right_offset);
                chain_clusters.fragment_best_left = std::min( chain_clusters.fragment_best_left, offset);
                chain_clusters.fragment_best_right = std::min(chain_clusters.fragment_best_right, right_offset);
#ifdef DEBUG_CLUSTER
                cerr << "At a snarl cluster on read " << read_num << " with head " << tree_state.all_seeds->at(read_num)->at(std::get<2>(seed_index)).pos << 
                        " and distances " << offset << " " << right_offset << endl;
#endif

                if (std::get<0>(seed_index) != last_snarl_by_read[read_num]) {

                    //If we're moving on to a new snarl, update the most recent snarl clusters
                    prev_snarl_cluster_by_read[read_num].insert(prev_snarl_cluster_by_read[read_num].end(),
                        snarl_cluster_by_read[read_num].begin(), snarl_cluster_by_read[read_num].end());

                    last_snarl_by_read[read_num] = std::get<0>(seed_index);
                }


                //try clustering by fragment
                if (tree_state.fragment_distance_limit != 0) {
                    int64_t new_best_right = std::numeric_limits<size_t>::max();
                    if (std::get<0>(seed_index) != last_snarl_fragment) {
#ifdef DEBUG_CLUSTER
                        cerr << "\tOn a new snarl with offset in chain " << std::get<0>(seed_index) << endl;
#endif
                        //And the same for the fragment clusters
                        prev_snarl_cluster_fragment.insert(prev_snarl_cluster_fragment.end(),
                            snarl_cluster_fragment.begin(), snarl_cluster_fragment.end());
                        snarl_cluster_fragment.clear();
    
    
                        last_snarl_fragment = std::get<0>(seed_index);
                    }
                    if (!prev_snarl_cluster_fragment.empty()) {
                        //If we saw a snarl cluster in a previous snarl but haven't seen a seed cluster since then
#ifdef DEBUG_CLUSTER
                        cerr << "\tCompare to last snarl clusters fragment: " << endl;
#endif
                        for (size_t i = 0 ; i < prev_snarl_cluster_fragment.size() ; i++) {
                            int64_t distance_in_chain = offset + std::get<2>(prev_snarl_cluster_fragment[i]) - chain_length - 1;
#ifdef DEBUG_CLUSTER
                    //assert(distance_in_chain >= 0);
#endif
                            if (std::get<2>(prev_snarl_cluster_fragment[i])  != std::numeric_limits<size_t>::max() && offset != std::numeric_limits<size_t>::max() && 
                                    distance_in_chain <= tree_state.fragment_distance_limit) {
                                tree_state.fragment_union_find.union_groups(std::get<0>(prev_snarl_cluster_fragment[i]), 
                                                        std::get<2>(seed_index)+tree_state.read_index_offsets[read_num]);
                                new_best_right = std::min(new_best_right, std::get<2>(prev_snarl_cluster_fragment[i]));
#ifdef DEBUG_CLUSTER
                            cerr << "\t\tCombining fragment: snarl cluster on component " << connected_component_num << ", " << 
                                    tree_state.all_seeds->at(read_num)->at(std::get<2>(seed_index)).pos << 
                                    "(distances: " << offset << ", " << right_offset << ")" <<  
                                    " with snarl cluster" << " " << read_num << " " << std::get<0>(prev_snarl_cluster_fragment[i]) << " " << 
                                    "(distances: " << std::get<1>(prev_snarl_cluster_fragment[i]) << ", " <<  std::get<2>(prev_snarl_cluster_fragment[i]) << 
                                    ") with distance between them " << distance_in_chain << endl;
#endif
                            }
                        }
                    } 
                    //Always try to compare this to the most recent seed cluster
#ifdef DEBUG_CLUSTER
                    cerr << "\tCompare to last seed cluster fragment:" << endl;
#endif
                    int64_t distance_in_chain = offset - seed_cluster_fragment.second;
                    if (seed_cluster_fragment.second != std::numeric_limits<size_t>::max() && offset != std::numeric_limits<size_t>::max() && distance_in_chain <= tree_state.fragment_distance_limit) {
                        tree_state.fragment_union_find.union_groups(seed_cluster_fragment.first, 
                                                    std::get<2>(seed_index)+tree_state.read_index_offsets[read_num]);
#ifdef DEBUG_CLUSTER
                        cerr << "\t\tCombining fragment: snarl cluster on component " << connected_component_num << ", " << 
                                tree_state.all_seeds->at(read_num)->at(std::get<2>(seed_index)).pos << 
                                "(distances: " << offset << ", " << right_offset << ")" <<  
                                " with seed cluster" << " " << read_num << " " << seed_cluster_fragment.first << " " <<   
                                "(distances: " << seed_cluster_fragment.second << ", " <<   chain_length - seed_cluster_fragment.second + 1 << 
                                ") with distance between them " << distance_in_chain << endl;
#endif
                    }
                    snarl_cluster_fragment.emplace_back(std::get<2>(seed_index)+tree_state.read_index_offsets[read_num],
                        offset, std::min(new_best_right, right_offset));
                }


                //try clustering by read
                if (!prev_snarl_cluster_by_read[read_num].empty()) {
#ifdef DEBUG_CLUSTER
                    cerr << "\t Comparing to last snarl clusters read "  << endl;
#endif
                    //If we saw a snarl cluster in a previous snarl and haven't seen a seed cluster since then
                    for (size_t i = 0 ; i < prev_snarl_cluster_by_read[read_num].size() ; i++) {
                        int64_t distance_in_chain = offset + std::get<2>(prev_snarl_cluster_by_read[read_num][i]) - chain_length - 1;
                        if (std::get<2>(prev_snarl_cluster_by_read[read_num][i]) != std::numeric_limits<size_t>::max() && offset != std::numeric_limits<size_t>::max() &&
                                distance_in_chain <= tree_state.read_distance_limit) {
                            tree_state.read_union_find[read_num].union_groups(std::get<0>(prev_snarl_cluster_by_read[read_num][i]), 
                                                                          new_cluster_head);
#ifdef DEBUG_CLUSTER
                            cerr << "\t\tSnarl cluster" << read_num << " " << std::get<0>(prev_snarl_cluster_by_read[read_num][i]) << " " << 
                                    tree_state.all_seeds->at(read_num)->at(std::get<0>(prev_snarl_cluster_by_read[read_num][i])).pos << 
                                    "(distances: " << std::get<1>(prev_snarl_cluster_by_read[read_num][i]) << ", " <<  
                                     std::get<2>(prev_snarl_cluster_by_read[read_num][i]) <<  ")"; 
#endif
                            size_t new_group = tree_state.read_union_find[read_num].find_group(new_cluster_head);
                            if (new_group == new_cluster_head) {
                                to_erase.emplace_back(read_num, std::get<0>(prev_snarl_cluster_by_read[read_num][i]));
                            } else {
                                to_erase.emplace_back(read_num, new_cluster_head);
                            }
                            best_left = std::min(best_left, tree_state.read_cluster_dists[read_num][std::get<0>(prev_snarl_cluster_by_read[read_num][i])].first);
                            best_right = std::min(best_right, tree_state.read_cluster_dists[read_num][std::get<0>(prev_snarl_cluster_by_read[read_num][i])].second);
                            tree_state.read_cluster_dists[read_num][new_group] = make_pair(best_left, best_right);
                            std::get<0>(prev_snarl_cluster_by_read[read_num][i]) = new_group;
                            new_cluster_head = new_group;
#ifdef DEBUG_CLUSTER
                            cerr << "...Combining read: snarl cluster on component " << connected_component_num << ", " << 
                                    tree_state.all_seeds->at(read_num)->at(new_cluster_head).pos << 
                                    "(distances: " << offset << ", " << right_offset << ")" <<  
                                    " with snarl cluster with distance between them " << distance_in_chain<< endl;
#endif
                        }
                    }

                }
                //Always try to cluster with the last seed cluster
                int64_t distance_in_chain = offset - seed_cluster_by_read[read_num].second; 
#ifdef DEBUG_CLUSTER
                cerr << "\tComparing to last seed cluster read "  << read_num << " " << seed_cluster_by_read[read_num].first << " " <<  
                         tree_state.all_seeds->at(read_num)->at(seed_cluster_by_read[read_num].first).pos << 
                         "(distances: " << seed_cluster_by_read[read_num].second << ", " << chain_length - seed_cluster_by_read[read_num].second + 1 << endl;
#endif
                if (seed_cluster_by_read[read_num].second != std::numeric_limits<size_t>::max() && offset != std::numeric_limits<size_t>::max() && 
                        distance_in_chain <= tree_state.read_distance_limit) {

                    tree_state.read_union_find[read_num].union_groups(seed_cluster_by_read[read_num].first, 
                                                                      new_cluster_head);
                    size_t new_group = tree_state.read_union_find[read_num].find_group(new_cluster_head);
                    if (new_group == new_cluster_head) {
                        to_erase.emplace_back(read_num, seed_cluster_by_read[read_num].first);
                    } else {
                        to_erase.emplace_back(read_num, new_cluster_head);
                    }
                    best_left = std::min(best_left, seed_cluster_by_read[read_num].second);
                    best_right = std::min(best_right, chain_length-seed_cluster_by_read[read_num].second+1);
                    tree_state.read_cluster_dists[read_num][new_group] = make_pair(best_left, best_right);
                    seed_cluster_by_read[read_num].first = new_group;
                    new_cluster_head = new_group;
#ifdef DEBUG_CLUSTER
                        cerr << "\t\tCombining read: snarl cluster on component " << connected_component_num << ", " << 
                                tree_state.all_seeds->at(read_num)->at(new_cluster_head).pos << 
                                "(distances: " << offset << ", " << right_offset << ")" <<  
                                " with seed cluster with distance between them " << distance_in_chain << endl;
#endif

                    right_offset = std::min(right_offset, chain_length-seed_cluster_by_read[read_num].second+1);
                }

                //Update the most recent cluster we've seen
                snarl_cluster_by_read[read_num].emplace_back(new_cluster_head, offset, right_offset);
                 
            }
            for (pair<size_t, size_t> cluster_head : to_erase) {
                chain_clusters.read_cluster_heads.erase(cluster_head);

            }
            //Add the new group to cluster heads
            chain_clusters.read_cluster_heads.emplace(std::get<1>(seed_index), new_cluster_head);
        }
        if (chain_index.is_looping_chain) {
#ifdef DEBUG_CLUSTER
            cerr << "Updating distances for looping chain" << endl;
#endif
            //If the chain loops, then the clusters might be connected by
            //looping around the chain
            int64_t first_length = chain_index.prefix_sum[0]-1;

            //New cluster- there will be at most one new cluster to add
            vector<size_t> combined_cluster (tree_state.all_seeds->size(), std::numeric_limits<size_t>::max());
            size_t fragment_combined_cluster = std::numeric_limits<size_t>::max();

            for (pair<size_t, size_t> cluster_head : chain_clusters.read_cluster_heads) {
                //For each chain cluster
                size_t read_num = cluster_head.first;
                pair<int64_t, int64_t>& chain_dists = tree_state.read_cluster_dists[read_num][cluster_head.second];

                if ((chain_dists.second != std::numeric_limits<size_t>::max() && chain_clusters.read_best_left[read_num] != std::numeric_limits<size_t>::max() &&
                     chain_dists.second + chain_clusters.read_best_left[read_num] - first_length - 1 <= tree_state.read_distance_limit) ||
                   (chain_dists.first != std::numeric_limits<size_t>::max() && chain_clusters.read_best_right[read_num] != std::numeric_limits<size_t>::max() &&
                      chain_dists.first + chain_clusters.read_best_right[read_num] - first_length - 1 <= tree_state.read_distance_limit)){
                    //If this chain cluster is in the combined cluster
                    if (combined_cluster[read_num] == std::numeric_limits<size_t>::max()) {
                        combined_cluster[read_num] = cluster_head.second;
                    } else {
                        tree_state.read_union_find[read_num].union_groups(combined_cluster[read_num], cluster_head.second);
                        if (tree_state.fragment_distance_limit != 0) {
                            if (fragment_combined_cluster != std::numeric_limits<size_t>::max()) {
                                tree_state.fragment_union_find.union_groups(fragment_combined_cluster, cluster_head.second+tree_state.read_index_offsets[read_num]);
                            }
                            fragment_combined_cluster = tree_state.fragment_union_find.find_group(cluster_head.second+tree_state.read_index_offsets[read_num]);
                        }
                        size_t new_group = tree_state.read_union_find[read_num].find_group(cluster_head.second);
                        combined_cluster[read_num] = new_group;
                    }

                    if (tree_state.fragment_distance_limit != 0) {
                        fragment_combined_cluster = tree_state.fragment_union_find.find_group(cluster_head.second + tree_state.read_index_offsets[read_num]);
                    }
                } else if (tree_state.fragment_distance_limit != 0 &&
                   ((chain_dists.second != std::numeric_limits<size_t>::max() && chain_clusters.fragment_best_left != std::numeric_limits<size_t>::max() &&
                     chain_dists.second + chain_clusters.fragment_best_left - first_length - 1 <= tree_state.fragment_distance_limit) ||
                   (chain_dists.first != std::numeric_limits<size_t>::max() && chain_clusters.fragment_best_right != std::numeric_limits<size_t>::max() &&
                      chain_dists.first + chain_clusters.fragment_best_right - first_length - 1 <= tree_state.fragment_distance_limit))){
                    //If we can cluster by fragment
                    if (fragment_combined_cluster != std::numeric_limits<size_t>::max()) {
                        tree_state.fragment_union_find.union_groups(fragment_combined_cluster, cluster_head.second+tree_state.read_index_offsets[read_num]);
                    }
                    fragment_combined_cluster = tree_state.fragment_union_find.find_group(cluster_head.second+tree_state.read_index_offsets[read_num]);

                }
            }
            //Don't need to update cluster heads or best left and right distances because
            //a looping chain will be the top level chain

        }

#ifdef DEBUG_CLUSTER
        cerr << "Found clusters on chain " << chain_index.id_in_parent << endl;
        cerr << "best left : " << chain_clusters.fragment_best_left << " best right : "
             << chain_clusters.fragment_best_right << endl;
        for (pair<size_t, size_t> c : chain_clusters.read_cluster_heads) {
        }
        bool got_left = false;
        bool got_right = false;
        for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
            cerr << " for read num " << read_num << " best left: " << chain_clusters.read_best_left[read_num] << " best right: " << chain_clusters.read_best_right[read_num] << endl;
            bool got_read_left=false;
            bool got_read_right = false;
            bool any_clusters = false;
            cerr << chain_clusters.read_cluster_heads.size() << " clusters on this chain: " << endl;
            for (pair<size_t,size_t> c : chain_clusters.read_cluster_heads) {
                if (c.first == read_num) {
                    any_clusters = true;
                    pair<int64_t, int64_t> dists = tree_state.read_cluster_dists[c.first][c.second];
                    cerr << "\tcluster headed by " << tree_state.all_seeds->at(c.first)->at(c.second).pos << ": left: " << dists.first << " right : " << dists.second << ": ";
                    bool has_seeds = false;
                    for (size_t x = 0 ; x < tree_state.all_seeds->at(c.first)->size() ; x++) {
                        if (tree_state.read_union_find[c.first].find_group(x) == c.second) {
                            cerr << tree_state.all_seeds->at(c.first)->at(x).pos << " ";
                            has_seeds = true;
                        }
                    }
                    if (depth != 0) {
                        assert(dists.first == std::numeric_limits<size_t>::max() || dists.first >= chain_clusters.read_best_left[read_num]);
                        assert(dists.second == std::numeric_limits<size_t>::max() || dists.second >= chain_clusters.read_best_right[read_num]);
                        assert(dists.first == std::numeric_limits<size_t>::max() || dists.first >= chain_clusters.fragment_best_left);
                        assert(dists.second == std::numeric_limits<size_t>::max() || dists.second >= chain_clusters.fragment_best_right);
                        assert(has_seeds);
                    }
                    if (dists.first == chain_clusters.fragment_best_left) {got_left = true;}
                    if (dists.second == chain_clusters.fragment_best_right) {got_right = true;}
                    if (dists.first == chain_clusters.read_best_left[read_num]) {got_read_left = true;}
                    if (dists.second == chain_clusters.read_best_right[read_num]) {got_read_right = true;}
                    cerr << endl;
                }
            }
            if (depth != 0 && !chain_index.is_looping_chain) {
                assert(!any_clusters || got_read_left || chain_clusters.read_best_left[read_num] > tree_state.read_distance_limit || chain_clusters.read_best_left[read_num] == std::numeric_limits<size_t>::max());
                assert(!any_clusters || got_read_right || chain_clusters.read_best_right[read_num] > tree_state.read_distance_limit || chain_clusters.read_best_right[read_num] == std::numeric_limits<size_t>::max());
            }
        }

        if (depth != 0 && !chain_index.is_looping_chain) {
            assert(got_left || chain_clusters.fragment_best_left > tree_state.fragment_distance_limit);
            assert(got_right ||chain_clusters.fragment_best_right > tree_state.fragment_distance_limit );
        }
        if (depth != 0) {
            for (pair<size_t, size_t> group_id : chain_clusters.read_cluster_heads) {

                assert (group_id.second == tree_state.read_union_find[group_id.first].find_group(group_id.second));
            }
        }
#endif
        return chain_clusters;

    };

/*
    void SnarlSeedClusterer::cluster_only_top_level_chain_seeds(TreeState& tree_state, 
                                                                  vector<pair<size_t, size_t>>& seed_clusters) const {
#ifdef DEBUG_CLUSTER
        cerr << "Clustering top-level seeds" << endl;
#endif


        std::sort(seed_clusters.begin(), seed_clusters.end(), [&](pair<size_t, size_t> item1, pair<size_t, size_t> item2) {
                   int64_t offset1 = tree_state.all_seeds->at(item1.first)->at(item1.second).offset ;
                   int64_t offset2 = tree_state.all_seeds->at(item2.first)->at(item2.second).offset ;
                   return offset1 < offset2; 
               });

        tuple<int64_t, size_t, size_t> last_fragment (-1, 0, 0);

        vector<tuple<int64_t, size_t, size_t>> last_by_read (tree_state.all_seeds->size(), make_tuple(-1, 0, 0));

        for (pair<size_t, size_t> seed_cluster : seed_clusters) {
            int64_t offset = tree_state.all_seeds->at(seed_cluster.first)->at(seed_cluster.second).offset ;

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


    hash_set<pair<size_t, size_t>> SnarlSeedClusterer::cluster_simple_snarl(TreeState& tree_state, vector<tuple<id_t, bool, int64_t, int64_t, int64_t>> nodes, 
                                int64_t loop_left, int64_t loop_right, int64_t snarl_length) const {
        //Cluster a top-level simple snarl and save the distances to the ends of the node in tree_state.read_cluster_dists 
        //Returns a vector of cluster heads (<read_num, seed num>)
#ifdef DEBUG_CLUSTER
        cerr << "cluster simple snarl " << endl;
        cerr << " loop distances: " << loop_left << " " << loop_right  << " snarl length " << snarl_length << endl;
#endif


 
        //Get the clusters of each node and find the best distances of any of them       
        vector<pair<size_t, size_t>> all_cluster_heads;
        int64_t fragment_best_left = -1; int64_t fragment_best_right = -1;
        vector<int64_t> best_left (tree_state.all_seeds->size(), -1);
        vector<int64_t> best_right (tree_state.all_seeds->size(), -1);


        for (tuple<id_t, bool, int64_t, int64_t, int64_t> node : nodes) {
            id_t node_id = std::get<0>(node);
            int64_t node_len = std::get<4>(node);

            //Cluster this node
            SnarlSeedClusterer::NodeClusters clusters =  cluster_one_node(tree_state, node_id, node_len);

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

            int64_t dist_left = tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].first;
            int64_t dist_right = tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].second;
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
                        int64_t best_left = std::min(dist_left, 
                                tree_state.read_cluster_dists[cluster_head.first][combined_read_left[cluster_head.first]].first);
                        int64_t best_right =std::min(dist_right, 
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
                        int64_t best_left = std::min(dist_left, 
                                tree_state.read_cluster_dists[cluster_head.first][combined_read_right[cluster_head.first]].first);
                        int64_t best_right = std::min(dist_right, 
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
