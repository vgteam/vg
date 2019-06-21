#include "seed_clusterer.hpp"

//#define DEBUG 

namespace vg {

    SnarlSeedClusterer::SnarlSeedClusterer( MinimumDistanceIndex& dist_index) : 
                                            dist_index(dist_index){
    };
        
    vector<vector<size_t>> SnarlSeedClusterer::cluster_seeds ( 
                  vector<pos_t> seeds, int64_t distance_limit ){
        /* Given a vector of seeds and a limit, find a clustering of seeds where
         * seeds that are closer than the limit cluster together. 
         * Returns a vector of cluster assignments 
         */ 
#ifdef DEBUG
cerr << endl << "New cluster calculation:" << endl;
#endif

        //For each level of the snarl tree, maps snarls (index into
        //dist_index.snarl_indexes) at that level to 
        //nodes belonging to the snarl
        //This is later used to populate snarl_to_node in the tree state
        vector<hash_map<size_t, vector<pair<NetgraphNode, NodeClusters>>>> 
                                                  snarl_to_nodes_by_level;
        snarl_to_nodes_by_level.resize(dist_index.tree_depth+1);


        //This stores all the tree relationships and cluster information
        //for a single level of the snarl tree as it is being processed
        //It also keeps track of the parents of the current level
        TreeState tree_state (&seeds, distance_limit);

        //Populate tree_state.node_to_seeds (mapping each node to the seeds it
        //contains) and snarl_to_nodes_by_level
        get_nodes(tree_state, snarl_to_nodes_by_level);

        //Initialize the tree state to be the bottom level
        if (dist_index.tree_depth >= 0) {
            tree_state.snarl_to_nodes = 
                          move(snarl_to_nodes_by_level[dist_index.tree_depth]);
        }

        for (int depth = dist_index.tree_depth ; depth >= 0 ; depth --) {
            //Go through each level of the tree, bottom up, and cluster that 
            // level. Each level includes the snarl at that level, the nodes
            // belonging to those snarls, and the chains comprised of them

            if (depth != 0) {
                // Bring in the direct child nodes that come in at this level 
                //  in the snarl tree.
                // They only ever occur below the root.
                tree_state.parent_snarl_to_nodes = 
                                       move(snarl_to_nodes_by_level[depth - 1]);
            }

            //Cluster all the snarls at this depth
            cluster_snarls(tree_state, depth);

            //And all the chains
            cluster_chains(tree_state, depth); 

            // Swap buffer over for the next level
            tree_state.snarl_to_nodes = move(tree_state.parent_snarl_to_nodes);
            tree_state.chain_to_snarls.clear();
        }

#ifdef DEBUG 
        cerr << "Found clusters : " << endl;
        for (auto group : tree_state.union_find_clusters.all_groups()){
            for (size_t c : group) {
               cerr << tree_state.seeds->at(c) << " ";
            }
            cerr << endl;
        }
#endif
        return tree_state.union_find_clusters.all_groups();
        
    };

    //Helper function to get the minimum value that is not -1 
    int64_t min_positive (int64_t n1, int64_t n2) {
        if (n1 == -1 ) {
            return n2;
        } else if (n2 == -1) {
            return n1;
        } else {
            return min(n1, n2);
        }
    }

    void SnarlSeedClusterer::get_nodes( TreeState& tree_state,
              vector<hash_map<size_t,vector<pair<NetgraphNode, NodeClusters>>>>&
                                                               snarl_to_nodes) {

        /* Find the nodes containing seeds and
         * assign each node to a level in the snarl tree*/ 

        for (size_t i = 0; i < tree_state.seeds->size(); i++) {
            //For each seed, assign it to a node and the node to a snarl 

            pos_t pos = tree_state.seeds->at(i);
            id_t id = get_id(pos);

            size_t snarl_i = dist_index.getPrimaryAssignment(id);

            MinimumDistanceIndex::SnarlIndex* snarl_index =
                                             &dist_index.snarl_indexes[snarl_i];
            size_t depth = snarl_index->depth;

            auto s = tree_state.node_to_seeds.find(id);
            if (s != tree_state.node_to_seeds.end()) {
                s->second.push_back(i);
            } else {

                //Map node to seed
                tree_state.node_to_seeds.emplace(id, vector<size_t>({i}));
                int64_t node_length = snarl_index->nodeLength(
                                                dist_index.getPrimaryRank(id));

                snarl_to_nodes[depth][snarl_i].emplace_back(
                         NetgraphNode(id, NODE), NodeClusters());
            }
        }
    }


    void SnarlSeedClusterer::cluster_snarls(TreeState& tree_state, size_t depth) {

        for (auto& kv : tree_state.snarl_to_nodes){
            //Go through each of the snarls at this level, cluster them,
            //and find which chains they belong to, if any
            //key is the index of the snarl and value is a vector of pair of
            // NetgraphNode, NodeClusters
            size_t snarl_i = kv.first;
            MinimumDistanceIndex::SnarlIndex& snarl_index = 
                                     dist_index.snarl_indexes[snarl_i];

#ifdef DEBUG
            cerr << "At depth " << depth << " snarl number " << snarl_i
                << " headed by " << snarl_index.id_in_parent
                << " with children " << endl;
            for (auto it2 : kv.second) {
                cerr << "\t" << typeToString(it2.first.node_type) << " number " << it2.first.node_id << endl;
            }
#endif
            if (snarl_index.in_chain){
                //If this snarl is in a chain, cluster and add let the
                //tree state know which chain it belongs to

                size_t chain_assignment = dist_index.getChainAssignment(
                                                    snarl_index.parent_id);
                size_t chain_rank = dist_index.getChainRank(
                                                  snarl_index.id_in_parent);

                tree_state.chain_to_snarls[chain_assignment].emplace(
                        chain_rank, make_pair(snarl_i, 
                            cluster_one_snarl(tree_state, snarl_i, 
                                              snarl_index.rev_in_parent)));

#ifdef DEBUG
                cerr << "Recording snarl number " << snarl_i << " headed by " << snarl_index.id_in_parent
                    << " as a child of chain number " << chain_assignment << " headed by " << snarl_index.parent_id << endl;
#endif
                
            } else {
                //If this snarl is not in a chain, add it as a child of the
                //parent snarl for the next level

                if (depth != 0 && snarl_index.parent_id != 0){
                    //If this has a parent, record it
#ifdef DEBUG
                    assert(snarl_index.parent_id >= dist_index.min_node_id);
                    assert(snarl_index.parent_id <= dist_index.max_node_id);
#endif
                    size_t parent_snarl_i = 
                           dist_index.getPrimaryAssignment(
                                                snarl_index.parent_id);
                    
                    tree_state.parent_snarl_to_nodes[parent_snarl_i].emplace_back(
                            NetgraphNode (snarl_i, SNARL), 
                            cluster_one_snarl(tree_state, snarl_i, false));
                                     
#ifdef DEBUG
                    cerr << "Recording snarl number " << snarl_i << " headed by " << snarl_index.id_in_parent
                        << " as a child of snarl number " << parent_snarl_i
                        << " headed by " << snarl_index.parent_id << endl;
                    cerr << "Snarl number " << parent_snarl_i << " has "
                        << tree_state.parent_snarl_to_nodes[parent_snarl_i].size() << " children now" << endl;
#endif
                } else {
                    //If it doesn't have a parent, don't need to keep track of
                    //its contents
                    cluster_one_snarl(tree_state, snarl_i, false);
                }
            }
        }
    }

    void SnarlSeedClusterer::cluster_chains(TreeState& tree_state, size_t depth) {
        for (auto& kv : tree_state.chain_to_snarls) {
            //For each chain at this level that has relevant child snarls in it,
            //find the clusters.

            // Get the chain's number
            size_t chain_i = kv.first;
            
#ifdef DEBUG
            cerr << "At depth " << depth << " chain number " << chain_i << " with children " << endl;
            for (auto it2 : kv.second) {
                cerr << "\t snarl number " << it2.second.first << endl;
            }
#endif

            // Compute the clusters for the chain
            auto chain_clusters = cluster_one_chain(tree_state, chain_i);
            
            if (depth > 0) {
                // We actually have a parent
                
                // Find the node ID that heads the parent of that chain.
                size_t parent_id = dist_index.chain_indexes[chain_i].parent_id;
                // It must be a legitimate node ID we cover.
#ifdef DEBUG 
                assert(parent_id >= dist_index.min_node_id);
                assert(parent_id <= dist_index.max_node_id);
#endif
                // Map it to the snarl number that should be represented by it 
                // (and thus also contain the chain)
                size_t parent_snarl_i =dist_index.getPrimaryAssignment(parent_id);
                
                // Register clusters as relevant for that parent snarl.
                
                tree_state.parent_snarl_to_nodes[parent_snarl_i].emplace_back(
                      NetgraphNode (chain_i, CHAIN), std::move(chain_clusters));
#ifdef DEBUG 
                cerr << "Recording chain number " << chain_i << " headed by " << dist_index.chain_indexes[chain_i].id_in_parent
                    << " as a child of snarl number " << parent_snarl_i << " headed by " << parent_id << endl;
                cerr << "Snarl number " << parent_snarl_i << " has "
                    << tree_state.parent_snarl_to_nodes[parent_snarl_i].size() << " children now" << endl;
#endif
            } 
        }
    }
    SnarlSeedClusterer::NodeClusters SnarlSeedClusterer::cluster_one_node(
                       TreeState& tree_state,
                       id_t node_id, int64_t node_length) {
#ifdef DEBUG 
        cerr << "Finding clusters on node " << node_id << " which has length " <<
        node_length << endl;
#endif
        /*Find clusters of seeds in this node. 
         * Returns a hash_set of the union find group IDs of the new clusters,
         * and the shortest distance from any seed to the left and right sides 
         * of the node*/

         vector<size_t>& seed_indices = tree_state.node_to_seeds[node_id];

        if (tree_state.distance_limit > node_length) {
            //If the limit is greater than the node length, then all the 
            //seeds on this node must be in the same cluster
            
            size_t group_id = seed_indices[0];

            int64_t best_dist_left  = -1;
            int64_t best_dist_right = -1;
            for (size_t seed_i : seed_indices) {
                //For each seed on this node, add it to the cluster
                
                pos_t seed = tree_state.seeds->at(seed_i); 
                int64_t dist_left = is_rev(seed) ? node_length- get_offset(seed)
                                                 : get_offset(seed) + 1;
                int64_t dist_right = is_rev(seed) ? get_offset(seed) + 1 
                                               : node_length - get_offset(seed);

                best_dist_left = min_positive(dist_left, best_dist_left);
                best_dist_right = min_positive(dist_right, best_dist_right);

                tree_state.union_find_clusters.union_groups(group_id, seed_i);

            }

            group_id = tree_state.union_find_clusters.find_group(group_id);
            tree_state.cluster_dists[group_id] = make_pair(best_dist_left, 
                                                best_dist_right);
            hash_set<size_t> cluster_group_ids;
            cluster_group_ids.insert(group_id);
#ifdef DEBUG 
            assert (group_id == tree_state.union_find_clusters.find_group(group_id));
            cerr << "Found single cluster on node " << node_id << endl;
            bool got_left = false;
            bool got_right = false;
            for (size_t c : cluster_group_ids) {
                pair<int64_t, int64_t> dists = tree_state.cluster_dists[c];
                assert(dists.first == -1 || dists.first >= best_dist_left);
                assert(dists.second == -1 || dists.second >= best_dist_right);
                if (dists.first == best_dist_left) {got_left = true;}
                if (dists.second == best_dist_right) {got_right = true;}
                cerr << "\t" << c << ": left: " << dists.first << " right : " <<
                                                           dists.second << endl;
            }
            assert(got_left);
            assert(got_right);
#endif
            return NodeClusters(std::move(cluster_group_ids), 
                              best_dist_left, best_dist_right);
        }

        //indices of union find group ids of clusters in this node
        hash_set<size_t> cluster_group_ids;
        int64_t best_left = -1;
        int64_t best_right = -1;

        //Create a vector of seeds with their offsets
        vector<pair<size_t, int64_t>> seed_offsets;
        for ( size_t i : seed_indices) {
            //For each seed, find its offset 
            pos_t seed = tree_state.seeds->at(i); 
            int64_t offset = is_rev(seed) ? node_length - get_offset(seed) 
                                            : get_offset(seed) + 1;

            best_left = min_positive(offset, best_left);
            best_right = min_positive(node_length-offset+1,best_right);

            seed_offsets.emplace_back(i, offset);
                        
        }
        //Sort seeds by their position in the node 
        std::sort(seed_offsets.begin(), seed_offsets.end(), 
                     [&](const auto a, const auto b) -> bool {
                          return a.second < b.second; 
                      } );

        int64_t last_offset = 0;
        size_t last_cluster = seed_offsets[0].first;
        int64_t last_left = -1; int64_t last_right = -1;
        cluster_group_ids.insert(last_cluster);

        for ( pair<size_t, int64_t> s : seed_offsets) {
            //For each seed, see if it belongs to a new cluster
            //i is initially its own group id

            size_t i_group = s.first;
            int64_t offset = s.second;

            if (abs(offset - last_offset) <= tree_state.distance_limit) {
                //If this seed is in the same cluster as the previous one,
                //union them

                tree_state.union_find_clusters.union_groups(i_group, last_cluster);
                last_cluster = tree_state.union_find_clusters.find_group(i_group);
                last_left = min_positive(last_left, offset);
                last_right = min_positive(last_right, node_length-offset+1);
                tree_state.cluster_dists[last_cluster] = make_pair(last_left, last_right);

            } else {
                //This becomes a new cluster
                cluster_group_ids.insert(i_group);
                last_cluster = i_group;
                last_left = offset;
                last_right = node_length - offset + 1;
                tree_state.cluster_dists[i_group] = make_pair(last_left, last_right);
            }
            last_offset = offset;
                        
        }
        
#ifdef DEBUG 
        cerr << "Found clusters on node " << node_id << endl;
        bool got_left = false;
        bool got_right = false;
        for (size_t c : cluster_group_ids) {
            pair<int64_t, int64_t> dists = tree_state.cluster_dists[c];
            assert(dists.first == -1 || dists.first >= best_left);
            assert(dists.first == -1 || dists.second >= best_right);
            if (dists.first == best_left) {got_left = true;}
            if (dists.second == best_right) {got_right = true;}
            cerr << "\t" << c << ": left: " << dists.first << " right : " 
                 << dists.second << endl;
        }
        assert(got_left );
        assert(got_right);
        for (size_t group_id : cluster_group_ids) {
            assert (group_id == tree_state.union_find_clusters.find_group(group_id));
        }
#endif
        return NodeClusters(std::move(cluster_group_ids), best_left, best_right);
        
    };



    SnarlSeedClusterer::NodeClusters SnarlSeedClusterer::cluster_one_chain(
                               TreeState& tree_state, size_t chain_index_i) {
        /*
         * Find all the clusters in the given chain, given clusters inside its 
         * child snarls.
         * Processes the child snarls and then combines into the chain.
         */
 
        std::map<size_t, pair<size_t, NodeClusters>>& snarls_in_chain = 
                                      tree_state.chain_to_snarls[chain_index_i];
  
        MinimumDistanceIndex::ChainIndex& chain_index = dist_index.chain_indexes[
                                                            chain_index_i];
#ifdef DEBUG 
        cerr << "Finding clusters on chain number " << chain_index_i << " headed by node " << chain_index.id_in_parent << endl;
#endif



        
        auto combine_snarl_clusters = [&] (size_t& new_group, 
                        size_t& combined_group, vector<size_t>& to_erase, 
                        pair<int64_t, int64_t>& dists){
            //Helper function to combine clusters of the same snarl

            if (combined_group == -1) {
                combined_group = new_group;
            } else {
                combined_group = tree_state.union_find_clusters.find_group(
                                                 combined_group);
                tree_state.union_find_clusters.union_groups(combined_group, new_group);
                pair<int64_t, int64_t>& old_dists =
                                   tree_state.cluster_dists[combined_group];
                size_t new_combined_group = 
                           tree_state.union_find_clusters.find_group(new_group);
                if (new_combined_group != new_group) {
                    to_erase.push_back(new_group);
                } 
                if (new_combined_group != combined_group)  {
                    to_erase.push_back(combined_group);
                }
                combined_group = new_combined_group;

                tree_state.cluster_dists[combined_group] = make_pair(
                      min_positive(old_dists.first, dists.first),
                      min_positive(old_dists.second, dists.second));
            }
            return;
        };
        hash_set<size_t> chain_cluster_ids;

        int64_t best_left = -1;
        int64_t best_right = -1;

        //The rank of the node at which the chain clusters reach
        size_t last_rank = 0;
        int64_t last_len = 0;
        id_t start_node;
        id_t end_node;

        for (auto& kv: snarls_in_chain) {
            /* For each child snarl in the chain, find the clusters of just the
             * snarl, and progressively build up clusters spanning up to that 
             * snarl
             * Snarls are in the order that they are traversed in the chain
             * snarls_in_chain is a map from rank of a snarl to the snarl index
             *  in dist_index.snarl_indexes and NodeClusters for that snarl
             */
            

            //rank of the boundary node of the snarl that occurs first in
            //the chain
            size_t start_rank = kv.first;
            size_t curr_snarl_i = kv.second.first;

            //Pull out the clusters of the current snarl
            hash_set<size_t>& snarl_clusters = kv.second.second.cluster_heads; 
            int64_t child_dist_left = kv.second.second.best_left;
            int64_t child_dist_right = kv.second.second.best_right; 

            MinimumDistanceIndex::SnarlIndex& snarl_index = 
                                         dist_index.snarl_indexes[curr_snarl_i];
            bool rev_in_chain = snarl_index.rev_in_parent;

            //Get the lengths of the start and end nodes of the snarl, relative
            //to the order of the chain
            int64_t start_length = rev_in_chain 
                         ? snarl_index.nodeLength(snarl_index.num_nodes * 2 - 1)
                         : snarl_index.nodeLength(0);
            int64_t end_length = rev_in_chain ? snarl_index.nodeLength(0)
                        : snarl_index.nodeLength(snarl_index.num_nodes * 2 - 1);

            //Distance from right end of chain clusters to the end of current
            //snarl cluster
            int64_t snarl_length = snarl_index.snarlLength();
            int64_t dist_to_end = snarl_length - start_length;


            if (last_rank != start_rank) { 
                /* If the chain clusters don't reach this snarl,
                 * extend their dist_right to the beginning of this snarl
                 */
                int64_t offset = chain_index.chainDistance(
                         make_pair(last_rank, false), 
                         make_pair(start_rank, false), last_len, start_length);
                offset = offset - last_len + start_length;

                for (size_t i : chain_cluster_ids) {
                    tree_state.cluster_dists[i].second = tree_state.cluster_dists[i].second == -1 
                                        ? -1 : tree_state.cluster_dists[i].second + offset;
                }
                best_right = best_right == -1 ? -1 : best_right + offset;
            }

            last_rank = start_rank + 1;
            last_len = end_length;
           

            //Distance from the start of chain to the start of the current snarl
            int64_t add_dist_left = start_rank == 0 ? 0 : 
                                    chain_index.prefix_sum[start_rank] - 1;


             
            //Combine snarl clusters that can be reached by looping
            int64_t loop_dist_end = chain_index.loop_fd[start_rank + 1] - 1 ;
            int64_t loop_dist_start = chain_index.loop_rev[start_rank] - 1; 

 
#ifdef DEBUG
            cerr << "  Snarl distance limits: " << child_dist_left << " " << child_dist_right << endl;
            cerr << "  Snarl clusters to add: " << endl;
            for (size_t c : snarl_clusters) {
                pair<int64_t, int64_t> dists = tree_state.cluster_dists[c];
                cerr << "\tleft: " << dists.first << " right : " << dists.second 
                     << endl;
                cerr << "\t\t";
                for (size_t x = 0 ; x < tree_state.seeds->size() ; x++) {
                    if (tree_state.union_find_clusters.find_group(x) == c) {
                        cerr << tree_state.seeds->at(x) << " ";
                    }
                }
            }
            cerr << endl;

            cerr << "  Clusters on chain: " << endl;

            cerr << "  best left: " << best_left << " best right: " << best_right << endl;
            for (size_t c : chain_cluster_ids) {
                pair<int64_t, int64_t> dists = tree_state.cluster_dists[c];
                cerr << "\tleft: " << dists.first << " right : " << dists.second 
                     << endl;
                cerr << "\t\t";
                for (size_t x = 0 ; x < tree_state.seeds->size() ; x++) {
                    if (tree_state.union_find_clusters.find_group(x) == c) {
                        cerr << tree_state.seeds->at(x) << " ";
                    }
                }
                cerr << endl;
            }
            cerr << endl;
#endif



            //Dist from end of start node to farthest seed to the left of snarl 
            //If the right dist of a chain cluster + chain_lim_right is less 
            //than the distance limit, then they are in the same cluster
            int64_t chain_lim_right = child_dist_left - start_length;
            int64_t snarl_lim_left = best_right - start_length;
            int64_t old_right = child_dist_left;//Need this in case it becomes -1
            int64_t old_left = best_right;

            vector<size_t> to_add;//new cluster group ids 
            vector<size_t> to_erase; //old cluster group ids
            //New cluster- there will be at most one new cluster to add
            size_t combined_cluster = -1;
            int64_t combined_left = -1; int64_t combined_right = -1;

            //Combined snarl clusters by taking chain loop left/right
            size_t snarl_cluster_left = -1;
            size_t snarl_cluster_right = -1;

            best_left = -1; best_right = -1;
            for (size_t j : snarl_clusters) {
                /* For each of the clusters for the current snarl,
                 * find if it belongs to the new combined cluster*/ 

                pair<int64_t, int64_t> snarl_dists = std::move(tree_state.cluster_dists[j]);



                if (loop_dist_start != -1) {
                    //If there is a loop going out and back into the start of
                    //the snarl, might combine this cluster with other snarl
                    //clusters

                    //Find new distance to the right side of the snarl
                    int64_t new_right = 
                              snarl_dists.first == -1 || loop_dist_start == -1
                                        ? -1
                                        : snarl_dists.first + loop_dist_start 
                                               + snarl_length - start_length;
                    snarl_dists.second = min_positive(new_right, 
                                                      snarl_dists.second);
                    child_dist_right =min_positive(child_dist_right, new_right);
#ifdef DEBUG
cerr << "Updating looping distance to right of snarl cluster" << j << ": " 
     << new_right << endl;
#endif
                    
                    
                    if (child_dist_left != -1 && snarl_dists.first != -1 && 
                        child_dist_left + snarl_dists.first 
                              + loop_dist_start - start_length - 1 
                              <= tree_state.distance_limit) {  
                        //If this cluster can be combined with another cluster
                        //from the left

#ifdef DEBUG
cerr << "  Combining this cluster from the left " << endl;
#endif
                        combine_snarl_clusters(j, snarl_cluster_left, 
                                               to_erase, snarl_dists);
                    }
                }

                if (loop_dist_end != -1) {
                    //If there is a loop to the right
                    int64_t new_left = 
                        snarl_dists.second == -1 || loop_dist_end == -1
                          ? -1
                          : snarl_dists.second + loop_dist_end + snarl_length 
                                         - end_length;
                    if (snarl_dists.first == -1 || (new_left != -1 &
                                                   new_left < snarl_dists.first)){
                        //If this is an improvement, update distances
                        snarl_dists.first = new_left; 
                        child_dist_left = min_positive(new_left, child_dist_left);
                        chain_lim_right = old_right == -1 
                                         ? child_dist_left-start_length 
                                         : min(chain_lim_right,
                                               child_dist_left - start_length);
                        old_right = min_positive(old_right, child_dist_left);
                    }

#ifdef DEBUG
cerr << "Updating looping distance to left of snarl cluster" << j << ": " 
     << new_left << endl;
#endif

                    if (child_dist_right != -1 && snarl_dists.second != -1 && 
                        child_dist_right + snarl_dists.second + loop_dist_end 
                                         - end_length - 1 <= tree_state.distance_limit) {  
                        //If this cluster can be combined with another cluster
                        //from the right 

#ifdef DEBUG
cerr << "  Combining this cluster from the right" << endl;
#endif
                        combine_snarl_clusters(j, snarl_cluster_right, 
                                               to_erase, snarl_dists);
                    }
                }

                //Now check if this snarl cluster can be combined with any 
                //existing chain clusters
                if (old_left != -1 && snarl_dists.first != -1 &&
                    snarl_dists.first + snarl_lim_left-1 <= tree_state.distance_limit) {
                    //If this snarl cluster is in the combined cluster
                    if (combined_cluster == -1) {
                        combined_cluster = j;
                        combined_left = snarl_dists.first == -1 ? -1 :
                                            snarl_dists.first + add_dist_left;
                        combined_right = snarl_dists.second;
                    } else {
                        tree_state.union_find_clusters.union_groups(combined_cluster, j);
                        size_t new_group = tree_state.union_find_clusters.find_group(j);
                        combined_cluster = new_group;
                        combined_left = min_positive(combined_left,
                                            snarl_dists.first == -1 ? -1 :
                                            snarl_dists.first + add_dist_left);
                        combined_right = min_positive(combined_right,
                                                           snarl_dists.second);
                    }
                        
                } else {
                    //This snarl becomes a new chain cluster
                    to_add.push_back(j);
                    pair<int64_t, int64_t> d = make_pair(snarl_dists.first == -1
                                      ? -1 : snarl_dists.first + add_dist_left, 
                                                snarl_dists.second);
                    best_left = min_positive(best_left, d.first); 
                    best_right = min_positive(best_right, d.second);

                    tree_state.cluster_dists[j] = std::move(d); 
                }
            }
            for (size_t i : chain_cluster_ids) {
                //For each old chain cluster
                pair<int64_t, int64_t>& chain_dists = tree_state.cluster_dists[i];

                if (old_right != -1 && chain_dists.second != -1 && 
                    chain_dists.second + chain_lim_right-1 <= tree_state.distance_limit){
                    //If this chain cluster is in the combined cluster
                    if (combined_cluster == -1) {
                        combined_cluster = i;
                        combined_left = chain_dists.first;
                        combined_right = chain_dists.second + dist_to_end;
                    } else {
                        tree_state.union_find_clusters.union_groups(combined_cluster, i);
                        size_t new_group = tree_state.union_find_clusters.find_group(i);
                        if (new_group == i) {
                            to_erase.push_back(combined_cluster);
                        } else {
                            to_erase.push_back(i);
                        }
                        combined_cluster = new_group;
                        combined_left = min_positive(combined_left,
                                                            chain_dists.first);
                        combined_right = min_positive(combined_right,
                                             chain_dists.second + dist_to_end);
                    }
                } else {
                    //If this chain cluster is on its own, extend its right 
                    //distance to the end of the current snarl
                    chain_dists.second += dist_to_end;
                    if (chain_dists.first - 2 >= tree_state.distance_limit &&
                        chain_dists.second - end_length-2 >= tree_state.distance_limit) {
                        //If the distance from the seeds in this cluster to
                        //either end of the chain is greater than the distance
                        //limit, then it cannot cluster with anything else
                        //so we can stop keeping track of it
#ifdef DEBUG
                        cerr << "Removing cluster " << i << endl;
#endif
                        to_erase.push_back(i);
                    } else {
                        best_left = min_positive(best_left, chain_dists.first);
                        best_right = min_positive(best_right, 
                                                  chain_dists.second);
                    }

                }
            }
            for (size_t j : to_add) {
                chain_cluster_ids.insert(j);
            }
            for (size_t j : to_erase) {
                chain_cluster_ids.erase(j);
            }
            if (combined_cluster != -1 ) {
                chain_cluster_ids.insert(combined_cluster);
                tree_state.cluster_dists[combined_cluster] = 
                                      make_pair(combined_left, combined_right);
                best_left = min_positive(best_left, combined_left);
                best_right = min_positive(best_right, combined_right);
            }
                  
#ifdef DEBUG 
            cerr << "\t at snarl number " << curr_snarl_i << ", clusters:" <<endl;

            for (size_t c : chain_cluster_ids) {
                pair<int64_t, int64_t> dists = tree_state.cluster_dists[c];
                cerr << "\t\tleft: " << dists.first << " right : " << dists.second << endl;
                cerr << "\t\t\t";
                for (size_t x = 0 ; x < tree_state.seeds->size() ; x++) {
                    if (tree_state.union_find_clusters.find_group(x) == c) {
                        cerr << tree_state.seeds->at(x) << " ";
                    }
                }
                cerr << endl;
            }
#endif
                
        }
         

        if (last_rank != chain_index.prefix_sum.size() - 2) {
            best_right = -1;
           //Extend the right bound of each cluster to the end of the chain
           int64_t last_dist = last_rank == 0 ? 0 :
                                    chain_index.prefix_sum[last_rank] - 1;
           int64_t dist_to_end = chain_index.chainLength()
                       - last_dist - last_len;
           for (size_t i : chain_cluster_ids) {
               int64_t d = tree_state.cluster_dists[i].second;
               tree_state.cluster_dists[i].second = d == -1 ? -1 : d + dist_to_end;
               best_right = min_positive(best_right, tree_state.cluster_dists[i].second);
           }
        }


        if (chain_index.is_looping_chain) {
            //If the chain loops, then the clusters might be connected by 
            //looping around the chain
            //
            int64_t first_length = chain_index.prefix_sum[0]-1;
            vector<size_t> to_erase; //old cluster group ids
            //New cluster- there will be at most one new cluster to add
            size_t combined_cluster = -1;

            for (size_t i : chain_cluster_ids) {
                //For each chain cluster
                pair<int64_t, int64_t>& chain_dists = tree_state.cluster_dists[i];

                if ((chain_dists.second != -1 && best_left != -1 &&
                     chain_dists.second + best_left - first_length - 1 
                                                <= tree_state.distance_limit) ||
                   (chain_dists.first != -1 && best_right != -1 && 
                      chain_dists.first + best_right - first_length - 1 
                                                <= tree_state.distance_limit)){
                    //If this chain cluster is in the combined cluster
                    if (combined_cluster == -1) {
                        combined_cluster = i;
                    } else {
                        tree_state.union_find_clusters.union_groups(combined_cluster, i);
                        size_t new_group = tree_state.union_find_clusters.find_group(i);
                        if (new_group == i) {
                            to_erase.push_back(combined_cluster);
                        } else {
                            to_erase.push_back(i);
                        }
                        combined_cluster = new_group;
                    }
                }
            }
            for (size_t i : to_erase) {
                chain_cluster_ids.erase(i);
            }
            //Don't need to update best left and right distances because
            //a looping chain will be the top level chain

        }

#ifdef DEBUG 
        cerr << "Found clusters on chain " << chain_index.id_in_parent << endl;
        cerr << "best left : " << best_left << " best right : " << best_right << endl;
        for (size_t c : chain_cluster_ids) {
            cerr << "\t";
            for (size_t x = 0 ; x < tree_state.seeds->size() ; x++) {
                if (tree_state.union_find_clusters.find_group(x) == c) {
                    cerr << tree_state.seeds->at(x) << " ";
                }
            }
            cerr << endl;
        }
        bool got_left = false;
        bool got_right = false;
        for (size_t c : chain_cluster_ids) {
            pair<int64_t, int64_t> dists = tree_state.cluster_dists[c];
            if (!chain_index.is_looping_chain){
                assert(dists.first == -1 || dists.first >= best_left);
                assert(dists.second == -1 || dists.second >= best_right);
            }
            if (dists.first == best_left) {got_left = true;}
            if (dists.second == best_right) {got_right = true;}
            cerr << "\t" << c << ": left: " << dists.first << " right : " 
                 << dists.second << endl;
        }
        if (!chain_index.is_looping_chain) {
            assert(got_left);
            assert(got_right);
        }
        for (size_t group_id : chain_cluster_ids) {

            assert (group_id == tree_state.union_find_clusters.find_group(group_id));
        }
#endif

        return NodeClusters(std::move(chain_cluster_ids), best_left, best_right) ; 
    };



    SnarlSeedClusterer::NodeClusters SnarlSeedClusterer::cluster_one_snarl(
                    TreeState& tree_state, size_t snarl_index_i, bool rev) {
        /*Get the clusters on this snarl. 
         * Nodes have not yet been clustered */
        MinimumDistanceIndex::SnarlIndex& snarl_index = 
                                        dist_index.snarl_indexes[snarl_index_i];
#ifdef DEBUG 
        cerr << "Finding clusters on snarl number " << snarl_index_i << " headed by node " << snarl_index.id_in_parent << endl;
#endif
        //Get the children of this snarl and their clusters
        vector<pair<NetgraphNode, NodeClusters>>& child_nodes = 
                                       tree_state.snarl_to_nodes[snarl_index_i];
        int64_t start_length = snarl_index.nodeLength(0);
        int64_t end_length = snarl_index.nodeLength(snarl_index.num_nodes*2 -1);

        //Group ids of all clusters on this snarl
        hash_set<size_t> snarl_cluster_ids;
        int64_t best_left = -1;
        int64_t best_right = -1;
 
        auto combine_clusters = [&] (size_t& new_group, size_t& combined_group, 
                                    pair<int64_t, int64_t>& dists){
            //Helper function to combine clusters in two nodes of the same snarl
            if (combined_group == -1) {
                snarl_cluster_ids.insert(new_group);
                tree_state.cluster_dists[new_group] = dists;
                combined_group = new_group;
            } else {

                combined_group = tree_state.union_find_clusters.find_group(combined_group);
                pair<int64_t, int64_t>old_dists = tree_state.cluster_dists[combined_group];
                tree_state.union_find_clusters.union_groups(new_group, combined_group);
                size_t new_g = tree_state.union_find_clusters.find_group(new_group);
                if (new_g != new_group) {
                    snarl_cluster_ids.erase(new_group);
                } 
                if (new_g != combined_group) {
                    snarl_cluster_ids.erase(combined_group);
                }
                snarl_cluster_ids.insert(new_g);
                dists = make_pair(
                            min_positive(dists.first, old_dists.first),
                            min_positive(dists.second, old_dists.second));
                tree_state.cluster_dists[new_g] = dists;
                new_group = new_g;
                combined_group = new_g;
            }
            return;
        };

        //Maps each cluster of child nodes to its left and right distances
        hash_map<size_t, pair<int64_t, int64_t>> old_dists;

        for (size_t i = 0; i < child_nodes.size() ; i++) {
            //Go through each child node of the netgraph and find clusters

            NetgraphNode& child = child_nodes [i].first;

            // Get the node id of this netgraph node in its parent snarl
            // Ranks in parents are computed from node ID, so we have to get it.
            id_t child_node_id = child.id_in_parent(dist_index);
            
            //Rank of this node in the snarl
            //If this node is a snarl/chain, then this snarl will be the
            //secondary snarl
            size_t node_rank = child.rank_in_parent(dist_index, child_node_id);
            size_t rev_rank = node_rank % 2 == 0
                           ? node_rank + 1 : node_rank - 1;

            if (child.node_type == NODE) {
                //If this node is a node, we need to find the clusters
                int64_t node_len = snarl_index.nodeLength(node_rank);

                child_nodes[i].second = cluster_one_node(
                                     tree_state, child_node_id, node_len);

            }
            NodeClusters& curr_child_clusters = child_nodes[i].second;
            int64_t child_dist_left = curr_child_clusters.best_left; 
            int64_t child_dist_right = curr_child_clusters.best_right;

#ifdef DEBUG
            cerr << "Finding distances to parent snarl " << snarl_index_i << " ends from child " << i << "/" << child_nodes.size() << endl;
            cerr << "Child is " << typeToString(child.node_type) << " number " << child.node_id << " headed by " << child_node_id << endl;
            cerr << "Node rank is " << node_rank << " fwd, " << rev_rank << " rev of " << snarl_index.num_nodes * 2 << endl;
            cerr << "Clusters at this child:" << endl; 
            for (size_t c : child_nodes[i].second.cluster_heads) {
                cerr << "\tdist left: " << tree_state.cluster_dists[c].first 
                << " dist right: " << tree_state.cluster_dists[c].second << endl;
                cerr << "\t\t";
                for (size_t x = 0 ; x < tree_state.seeds->size() ; x++) {
                    if (tree_state.union_find_clusters.find_group(x) == c) {
                        cerr << tree_state.seeds->at(x) << " ";
                    }
                }
                cerr << endl;
            }

            // Make sure the net graph node is actually in the net graph.
            assert(node_rank != numeric_limits<size_t>::max());
#endif

            vector<size_t> children_i( 
                  make_move_iterator(curr_child_clusters.cluster_heads.begin()),
                  make_move_iterator(curr_child_clusters.cluster_heads.end()));
            for (size_t c_i = 0 ; c_i < children_i.size() ; c_i ++) {
                //for each cluster of child node i, find the distances to the
                //ends of the snarl

                size_t c = children_i[c_i];
            
                pair<int64_t, int64_t> dists_c= tree_state.cluster_dists[c];
                old_dists[c] = dists_c;

                pair<int64_t, int64_t> new_dists = snarl_index.distToEnds(node_rank,
                                        dists_c.first,dists_c.second);
#ifdef DEBUG
cerr << "\tcluster: " << c_i << "dists to ends in snarl" << snarl_index.id_in_parent << " : " << new_dists.first << " " << new_dists.second << endl;
#endif

                best_left = min_positive(best_left, new_dists.first);
                best_right = min_positive(best_right, new_dists.second);


                snarl_cluster_ids.insert(c);
                tree_state.cluster_dists[c] = new_dists;
            }
            

            for (size_t j = 0 ; j <= i ; j++){
                //Go through other child net graph nodes up to and including i
                NetgraphNode& other_node = child_nodes[j].first;
                NodeClusters& other_node_clusters = child_nodes[j].second;

                // Note that other_node.second is the type, and other_node.first is the
                // *number* of the other_node in that type, *not* the heading node ID.
                // Ranks in parents are computed from node ID, so we have to get it.
                id_t other_node_id = other_node.id_in_parent(dist_index);
                
#ifdef DEBUG
                cerr << "Other net graph node is " << typeToString(other_node.node_type)
                    << " headed by node " << other_node_id;
                    
                    
#endif

                //Rank of this node in the snarl
                size_t other_rank = other_node.rank_in_parent(dist_index,  
                                                              other_node_id);
                size_t other_rev = other_rank % 2 == 0
                                    ? other_rank + 1 : other_rank - 1;

                //Find distance from each end of current node to 
                //each end of other node
                int64_t dist_l_l = snarl_index.snarlDistance(
                                                     rev_rank, other_rank);
                int64_t dist_l_r = snarl_index.snarlDistance(
                                                      rev_rank, other_rev);
                int64_t dist_r_l = snarl_index.snarlDistance(
                                                    node_rank, other_rank);
                int64_t dist_r_r = snarl_index.snarlDistance(
                                                     node_rank, other_rev);
#ifdef DEBUG
cerr << "\t distances between ranks " << node_rank << " and " << other_rank 
     << ": " << dist_l_l << " " << dist_l_r << " " << dist_r_l << " "  
     << dist_r_r << endl;
#endif

                //group ids of clusters combined between node i left and 
                //node j left, etc
                size_t group_l_l = -1;
                size_t group_l_r = -1;
                size_t group_r_l = -1;
                size_t group_r_r = -1;

                if (max({dist_l_l, dist_l_r, dist_r_l, dist_r_r}) != -1
                   && MinimumDistanceIndex::minPos({dist_l_l, dist_l_r, 
                            dist_r_l, dist_r_r})-2 <= tree_state.distance_limit
                   && min_positive(child_dist_left, child_dist_right)-2 
                                                <= tree_state.distance_limit) {
                    //If the two nodes are reachable
                    for (size_t c_i = 0 ; c_i < children_i.size() ; c_i ++) {
                        //for each cluster of child node i

                        size_t c = children_i[c_i];
                        size_t c_group = tree_state.union_find_clusters.find_group(c);

                        pair<int64_t, int64_t> new_dists;
                        pair<int64_t, int64_t> dists_c;

                         dists_c = old_dists[c];
                         new_dists = tree_state.cluster_dists[c_group];
                        

                        if (dist_l_l != -1 && dists_c.first != -1 
                                 && other_node_clusters.best_left != -1  
                                 && dist_l_l + dists_c.first 
                                   + other_node_clusters.best_left-1 <= tree_state.distance_limit) {
                            //If cluster c can be combined with clusters in j 
                            //from the left of both of them
                            combine_clusters(c_group, group_l_l, new_dists);
                        }
                        if (dist_l_r != -1 && dists_c.first != -1 
                            && other_node_clusters.best_right != -1 
                            && dist_l_r + dists_c.first + other_node_clusters.best_right-1
                                                 <= tree_state.distance_limit) {
                            combine_clusters(c_group, group_l_r, new_dists);

                        }
                        if (dist_r_l != -1 && dists_c.second != -1 
                            && other_node_clusters.best_left != -1 
                            && dist_r_l + dists_c.second + other_node_clusters.best_left-1
                                                 <= tree_state.distance_limit) {
                            combine_clusters(c_group, group_r_l, new_dists);
                        }
                        if (dist_r_r != -1 && dists_c.second != -1 
                            && other_node_clusters.best_right != -1 
                           && dist_r_r + dists_c.second + other_node_clusters.best_right-1
                                                 <= tree_state.distance_limit) {
                            combine_clusters(c_group, group_r_r, new_dists);
                        }

                    }
                    //Go through children of j
                    vector<size_t> children_j( make_move_iterator(other_node_clusters.cluster_heads.begin()), 
                             make_move_iterator(other_node_clusters.cluster_heads.end()));

                    for (size_t k_i = 0 ; k_i < children_j.size() ; k_i++){
                        size_t k = children_j[k_i];
                        //For each cluster of child j, find which overlaps with
                        //clusters of i
                        //k will already be part of a cluster in 
                        //snarl_cluster_ids but since we need to know the node 
                        //that the snarl is on we can't just loop through 
                        //snarl_cluster_ids
                        pair<int64_t, int64_t>& dist_bounds_k = old_dists[k];
                        size_t k_group = tree_state.union_find_clusters.find_group(k);
                        pair<int64_t, int64_t> dists_k = tree_state.cluster_dists[k_group];
                    

                        if (dist_l_l != -1 && child_dist_left != -1 
                           && dist_bounds_k.first != -1 
                           && dist_l_l + child_dist_left + dist_bounds_k.first-1
                                                  <= tree_state.distance_limit){

                            combine_clusters(k_group, group_l_l, dists_k);
                        }
                        if (dist_l_r != -1 && child_dist_left != -1 
                             && dist_bounds_k.second != -1  
                          && dist_l_r + child_dist_left + dist_bounds_k.second-1
                                               <= tree_state.distance_limit ) {

                            combine_clusters(k_group, group_l_r, dists_k);
                        }
                        if (dist_r_l != -1 && child_dist_right != -1 
                            && dist_bounds_k.first != -1  
                          && dist_r_l + child_dist_right + dist_bounds_k.first-1
                                                <= tree_state.distance_limit) {

                            combine_clusters(k_group, group_r_l, dists_k);
                        }
                        if (dist_r_r != -1 && child_dist_right != -1 
                           && dist_bounds_k.second != -1
                         && dist_r_r + child_dist_right + dist_bounds_k.second-1
                                                 <= tree_state.distance_limit) {

                            combine_clusters(k_group, group_r_r, dists_k);
                        }
                    }
                }
            }
        }
#ifdef DEBUG 
        cerr << "Found clusters on snarl number " << snarl_index_i << " headed by" << snarl_index.id_in_parent << endl;
        cerr << "    with best left and right values: " << best_left << " " 
             << best_right << endl;
        bool got_left = false;
        bool got_right = false;
        for (size_t c : snarl_cluster_ids) {
            pair<int64_t, int64_t> dists = tree_state.cluster_dists[c];
            if (dists.first == best_left) {got_left = true;}
            if (dists.second == best_right) {got_right = true;}
            cerr << "\t" << c << ": left: " << dists.first << " right : " 
                 << dists.second << endl;
            cerr << "\t\t";
            for (size_t x = 0 ; x < tree_state.seeds->size() ; x++) {
                if (tree_state.union_find_clusters.find_group(x) == c) {
                    cerr << tree_state.seeds->at(x) << " ";
                }
            }
            cerr << endl;
        }
        assert(got_left);
        assert(got_right);
        for (size_t group_id : snarl_cluster_ids) {
            assert (group_id == tree_state.union_find_clusters.find_group(group_id));
        }
#endif
        return NodeClusters(std::move(snarl_cluster_ids), best_left, best_right);
    };
}
