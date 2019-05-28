#include "seed_clusterer.hpp"

//#define DEBUG 

namespace vg {

    SnarlSeedClusterer::SnarlSeedClusterer( MinimumDistanceIndex& dist_index) : 
                                            dist_index(dist_index){
    };
        
    vector<vector<size_t>> SnarlSeedClusterer::cluster_seeds ( 
                  vector<pos_t> seeds, size_t distance_limit ){
        /* Given a vector of seeds and a limit, find a clustering of seeds where
         * seeds that are closer than the limit cluster together. 
         * Returns a vector of cluster assignments 
         */ 
#ifdef DEBUG
cerr << endl << "New cluster calculation:" << endl;
#endif

        //Create a union find structure to hold the cluster assignments of
        //each seed. Each seed is initially its own cluster 
        structures::UnionFind union_find_clusters (seeds.size(), false);


        //Map each node to a vector of seed indices
        hash_map<id_t, vector<size_t>> node_to_seeds;

        //For each level of the snarl tree, maps snarls (index into
        //dist_index.snarl_indexes) at that level to 
        //nodes belonging to the snarl
        vector<hash_map<size_t, vector<pair<child_node_t, child_cluster_t>>>> 
                                                                 snarl_to_nodes;
        //TODO: reserve up to the depth of the snarl tree- also need to find the max depth of the snarl tree in the distance index 

        //Populate node_to_seed and snarl_to_nodes
        get_nodes(seeds, node_to_seeds, snarl_to_nodes);

        //Maps each cluster group ID to the left and right distances
        vector<pair<int64_t, int64_t>> cluster_dists (seeds.size(), 
                                                      make_pair(-1, -1));

        int tree_depth = snarl_to_nodes.size()-1; 

        //Maps snarls in the current level to be clustered (by snarl number) to its children.
        //Children are stored as child index in its type, and child type.
        hash_map<size_t, vector<pair<child_node_t, child_cluster_t>>> 
                                                            curr_snarl_children;
                  
        if (tree_depth >= 0) {
            curr_snarl_children = move(snarl_to_nodes[tree_depth]);
        }
        //Maps children to the results of clustering - cluster heads and dists
        hash_map<pair<id_t, bool>, child_cluster_t> child_clusters;

        for (int depth = tree_depth ; depth >= 0 ; depth --) {
            /*Go through each level of the tree, bottom up, and cluster that 
             * level. Each level includes the snarl at that level, the nodes
             * belonging to those snarls, and the chains comprised of them*/

            //For each snarl in the level above this one, map to snarls/chains 
            //in this level
            hash_map<size_t, vector<pair<child_node_t, child_cluster_t>>>
                                                          parent_snarl_children;
            if (depth != 0) {
                // Bring in the direct child nodes that come in at this level in the snarl tree.
                // They only ever occur below the root.
                parent_snarl_children = move(snarl_to_nodes[depth - 1]);
            }

            //Maps each chain to the snarls that comprise it
            //Snarls are unordered in the vector
            //snarl is a tuple of the rank of the snarl in the chain,
            //the index of the snarl in dist_index.snarl_indexes, and the
            //snarl's orientation in the chain
            hash_map<size_t, vector<tuple<size_t, size_t, bool>>> chain_to_snarls; 

            for (auto& kv : curr_snarl_children){
                //Go through each of the snarls at this level, cluster them,
                //and find which chains they belong to, if any
                //key is the index of the snarl and value is a vector of pair of
                // child_node_t, child_cluster_t
                size_t snarl_i = kv.first;
                MinimumDistanceIndex::SnarlIndex& snarl_index = 
                                         dist_index.snarl_indexes[snarl_i];

#ifdef DEBUG
                cerr << "At depth " << depth << " snarl number " << snarl_i
                    << " headed by " << snarl_index.id_in_parent
                    << " with children " << endl;
                for (auto it2 : kv.second) {
                    cerr << "\t" << typeToString(it2.first.second) << " number " << it2.first.first << endl;
                }
#endif
                if (snarl_index.in_chain){
                    //If this snarl is in a chain, add it to chains
                    //but don't cluster yet

                    size_t chain_assignment = dist_index.chain_assignments[
                           snarl_index.parent_id-dist_index.min_node_id]-1;
                    size_t rank = dist_index.chain_ranks[
                       snarl_index.id_in_parent-dist_index.min_node_id]-1;

                    chain_to_snarls[chain_assignment].push_back( 
                          make_tuple(rank, snarl_i, snarl_index.rev_in_parent));
                          
#ifdef debug
                    cerr << "Recording snarl number " << snarl_i << " headed by " << snarl_index.id_in_parent
                        << " as a child of chain number " << chain_assignment << " headed by " << snarl_index.parent_id << endl;
#endif
                    
                } else {
                    //If this snarl is not in a chain, add it as a child of the
                    //parent snarl for the next level

                    if (depth != 0 && snarl_index.parent_id != 0){
                        //If this has a parent, record it
                        assert(snarl_index.parent_id >= dist_index.min_node_id);
                        assert(snarl_index.parent_id <= dist_index.max_node_id);
                        size_t parent_snarl_i = 
                               dist_index.primary_snarl_assignments[
                                  snarl_index.parent_id-dist_index.min_node_id] - 1;
                        parent_snarl_children[parent_snarl_i].push_back(
                            make_pair(make_pair(snarl_i, SNARL), 
                                  get_clusters_snarl(seeds, union_find_clusters,
                                         cluster_dists, kv.second, 
                                         node_to_seeds, distance_limit,
                                         snarl_i, false)));
                                         
#ifdef DEBUG
                        cerr << "Recording snarl number " << snarl_i << " headed by " << snarl_index.id_in_parent
                            << " as a child of snarl number " << parent_snarl_i
                            << " headed by " << snarl_index.parent_id << endl;
                        cerr << "Snarl number " << parent_snarl_i << " has "
                            << parent_snarl_children[parent_snarl_i].size() << " children now" << endl;
#endif
                    } else {
                         get_clusters_snarl(seeds, union_find_clusters,
                                    cluster_dists, kv.second, 
                                    node_to_seeds, distance_limit,
                                    snarl_i, false);
                    }
                }
            }
            for (auto& kv : chain_to_snarls) {
                //For each chain at this level that has relevant child snarls in it, find the clusters.

                // Get the chain's number
                size_t chain_i = kv.first;
                
#ifdef DEBUG
                cerr << "At depth " << depth << " chain number " << chain_i << " with children " << endl;
                for (auto it2 : kv.second) {
                    cerr << "\t snarl number " << get<1>(it2) << endl;
                }
#endif

                // Compute the clusters for the chain
                auto chain_clusters = get_clusters_chain(seeds, union_find_clusters,
                    cluster_dists, kv.second, curr_snarl_children,
                    node_to_seeds, distance_limit, chain_i);
                
                if (depth > 0) {
                    // We actually have a parent
                    
                    // Find the node ID that heads the parent of that chain.
                    size_t parent_id = dist_index.chain_indexes[chain_i].parent_id;
                    // It must be a legitimate node ID we cover.
                    assert(parent_id >= dist_index.min_node_id);
                    assert(parent_id <= dist_index.max_node_id);
                    // Map it to the snarl number that should be represented by it (and thus also contain the chain)
                    size_t parent_snarl_i = dist_index.primary_snarl_assignments[
                                parent_id-dist_index.min_node_id]-1;
                    
                    // Register clusters as relevant for that parent snarl.
                    parent_snarl_children[parent_snarl_i].emplace_back(make_pair(chain_i, CHAIN), std::move(chain_clusters));
#ifdef DEBUG
                    cerr << "Recording chain number " << chain_i << " headed by " << dist_index.chain_indexes[chain_i].id_in_parent
                        << " as a child of snarl number " << parent_snarl_i << " headed by " << parent_id << endl;
                    cerr << "Snarl number " << parent_snarl_i << " has "
                        << parent_snarl_children[parent_snarl_i].size() << " children now" << endl;
#endif
                } 
            }

            // Swap buffer over for the next level
            curr_snarl_children = move(parent_snarl_children);
        }

#ifdef DEBUG 
        cerr << "Found clusters : " << endl;
        for (auto group : union_find_clusters.all_groups()){
            for (size_t c : group) {
               cerr << seeds[c] << " ";
            }
            cerr << endl;
        }
#endif
        return union_find_clusters.all_groups();
        
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

    void SnarlSeedClusterer::get_nodes( 
              const vector<pos_t>& seeds,
              hash_map<id_t, vector<size_t>>& node_to_seeds,
              vector<hash_map<size_t, 
                     vector<pair<child_node_t, child_cluster_t>>>>& 
                                                               snarl_to_nodes) {

        /* Given a snarl tree and seeds, find the nodes containing seeds and
         * assign each node to a level in the snarl tree*/ 

        for (size_t i = 0; i < seeds.size(); i++) {
            //For each seed, assign it to a node and the node to a snarl 

            pos_t pos = seeds[i];
            id_t id = get_id(pos);

            size_t snarl_i = dist_index.primary_snarl_assignments[
                                          id - dist_index.min_node_id]-1;

            MinimumDistanceIndex::SnarlIndex* snarl_index =
                                             &dist_index.snarl_indexes[snarl_i];
            size_t depth = snarl_index->depth;

            auto s = node_to_seeds.find(id);
            if (s != node_to_seeds.end()) {
                s->second.push_back(i);
            } else {

                //Map node to seed
                node_to_seeds.emplace(id, vector<size_t>({i}));
                int64_t node_length = snarl_index->nodeLength(
                         dist_index.primary_snarl_ranks[
                                           id - dist_index.min_node_id]-1);

                //Map snarl to node
                if (snarl_to_nodes.size() < depth+1) {
                    snarl_to_nodes.resize(depth+1); 
                }
                child_cluster_t empty;
                snarl_to_nodes[depth][snarl_i].emplace_back(
                                  make_pair(id, NODE), empty);
            }
        }
    }


    SnarlSeedClusterer::child_cluster_t SnarlSeedClusterer::get_clusters_node(
                       const vector<pos_t>& seeds,
                       structures::UnionFind& union_find_clusters, 
                       vector< pair<int64_t, int64_t>>& cluster_dists,
                       vector<size_t>& seed_indices,
                       size_t distance_limit, 
                       id_t root, int64_t node_length) {
#ifdef DEBUG 
        cerr << "Finding clusters on node " << root << " which has length " <<
        node_length << endl;
#endif
        /*Find clusters of seeds in this node, root. 
         * Returns a hash_set of the union find group IDs of the new clusters,
         * and the shortest distance from any seed to the left and right sides 
         * of the node*/

        if (distance_limit > node_length) {
            //If the limit is greater than the node length, then all the 
            //seeds on this node must be in the same cluster
            
            size_t group_id = seed_indices[0];

            int64_t best_dist_left  = -1;
            int64_t best_dist_right = -1;
            for (size_t seed_i : seed_indices) {
                //For each seed on this node, add it to the cluster
                
                pos_t seed = seeds[seed_i]; 
                int64_t dist_start = get_offset(seed) + 1; 
                int64_t dist_end = node_length - get_offset(seed);
                int64_t dist_left = is_rev(seed) ? dist_end : dist_start;
                int64_t dist_right = is_rev(seed) ? dist_start : dist_end ;

                best_dist_left = min_positive(dist_left, best_dist_left);
                best_dist_right = min_positive(dist_right, best_dist_right);

                union_find_clusters.union_groups(group_id, seed_i);

            }

            group_id = union_find_clusters.find_group(group_id);
            cluster_dists[group_id] = make_pair(best_dist_left, 
                                                best_dist_right);
            hash_set<size_t> cluster_group_ids;
            cluster_group_ids.insert(group_id);
#ifdef DEBUG 
            assert (group_id == union_find_clusters.find_group(group_id));
            cerr << "Found single cluster on node " << root << endl;
            bool got_left = false;
            bool got_right = false;
            for (size_t c : cluster_group_ids) {
                pair<int64_t, int64_t> dists = cluster_dists[c];
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
            return make_tuple(move(cluster_group_ids), 
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
            pos_t seed = seeds[i]; 
            int64_t offset = is_rev(seed) ? node_length - get_offset(seed) 
                                            : get_offset(seed) + 1;

            best_left = min_positive(offset,  best_left);
            best_right = min_positive(node_length-offset+1,best_right);

            seed_offsets.push_back(make_pair(i, offset));
                        
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

            if (abs(offset - last_offset) < distance_limit) {
                //If this seed is in the same cluster as the previous one,
                //union them

                union_find_clusters.union_groups(i_group, last_cluster);
                last_cluster = union_find_clusters.find_group(i_group);
                last_left = min_positive(last_left, offset);
                last_right = min_positive(last_right, node_length-offset+1);
                cluster_dists[last_cluster] = make_pair(last_left, last_right);

            } else {
                //This becomes a new cluster
                cluster_group_ids.insert(i_group);
                last_cluster = i_group;
                last_left = offset;
                last_right = node_length - offset + 1;
                cluster_dists[i_group] = make_pair(last_left, last_right);
            }
            last_offset = offset;
                        
        }
        
#ifdef DEBUG 
        cerr << "Found clusters on node " << root << endl;
        bool got_left = false;
        bool got_right = false;
        for (size_t c : cluster_group_ids) {
            pair<int64_t, int64_t> dists = cluster_dists[c];
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
            assert (group_id == union_find_clusters.find_group(group_id));
        }
#endif
        return make_tuple(move(cluster_group_ids), best_left, best_right);
        
    };

    SnarlSeedClusterer::child_cluster_t SnarlSeedClusterer::get_clusters_chain(
                       const vector<pos_t>& seeds,
                       structures::UnionFind& union_find_clusters,
                       vector<pair<int64_t, int64_t>>& cluster_dists,
                       vector<  tuple<size_t, size_t, bool>>& snarls_in_chain,
                       hash_map<size_t, 
                                vector<pair<child_node_t, child_cluster_t>>>& 
                                                            curr_snarl_children,
                       hash_map<id_t, vector<size_t>>& node_to_seeds,
                       size_t distance_limit,  size_t chain_index_i) {
        /*
         * Find all the clusters in the given chain, given clusters inside its 
         * child snarls.
         * Processes the child snarls and then combines into the chain.
         */
 
        //Sort vector of snarls by rank
        std::sort(snarls_in_chain.begin(), snarls_in_chain.end(), 
                     [&](const auto s, const auto t) -> bool {
                          return std::get<0>(s) < std::get<0>(t); 
                      } );
  
        MinimumDistanceIndex::ChainIndex& chain_index = dist_index.chain_indexes[
                                                            chain_index_i];
#ifdef DEBUG 
        cerr << "Finding clusters on chain number " << chain_index_i << " headed by node " << chain_index.id_in_parent << endl;
#endif
        hash_set<size_t> chain_cluster_ids;

        int64_t best_left = -1;
        int64_t best_right = -1;

        //The rank of the node at which the chain clusters reach
        size_t last_rank = 0;
        int64_t last_len = 0;
        id_t start_node;
        id_t end_node;
        size_t prev_snarl_i = std::numeric_limits<size_t>::max();

        for (auto x : snarls_in_chain) {
            /* For each child snarl in the chain, find the clusters of just the
             * snarl, and progressively build up clusters spanning up to that 
             * snarl
             * Snarls are in the order that they are traversed in the chain
             */
            bool rev_in_chain = std::get<2>(x);
            //TODO: Don't need rev here, could just get it from the snarl index
            size_t curr_snarl_i = std::get<1>(x);

            //Skip duplicated snarls
            if (curr_snarl_i == prev_snarl_i) {
                continue;
            } else {
                prev_snarl_i = curr_snarl_i;
            }

            MinimumDistanceIndex::SnarlIndex& snarl_index = 
                                         dist_index.snarl_indexes[curr_snarl_i];

            //rank of the boundary node of the snarl that occurs first in
            //the chain
            size_t start_rank = std::get<0>(x);

            //Get the lengths of the start and end nodes of the snarl, relative
            //to the order of the chain
            int64_t start_length = snarl_index.nodeLength(0);
            int64_t end_length = snarl_index.nodeLength(
                                       snarl_index.num_nodes * 2 - 1);
            if (rev_in_chain) {
                int64_t tmp = start_length;
                start_length = end_length;
                end_length = tmp;
            }

            //Distance from right end of chain clusters to the end of current
            //snarl cluster
            int64_t snarl_length = snarl_index.snarlLength();
            int64_t dist_to_end = snarl_length - start_length;


            if (last_rank != start_rank) { 
                /* If the chain clusters don't reach this snarl,
                 * extend their dist_right to the beginning of this snarl
                 */
                //TODO: I think offset will always be positive
                int64_t offset = chain_index.chainDistance(
                         make_pair(last_rank, false), 
                         make_pair(start_rank, false), last_len, start_length);
                offset = offset - last_len + start_length;

                for (size_t i : chain_cluster_ids) {
                    cluster_dists[i].second = cluster_dists[i].second == -1 
                                        ? -1 : cluster_dists[i].second + offset;
                }
                best_right = best_right == -1 ? -1 : best_right + offset;
            }

#ifdef DEBUG
            cerr << "Going to get clusters for snarl number " << curr_snarl_i << " by clustering from its "
                << curr_snarl_children[curr_snarl_i].size() << " children" << endl;
#endif

            //Find the clusters of the current snarl
            hash_set<size_t> snarl_clusters; 
            int64_t child_dist_left; int64_t child_dist_right; 
            tie (snarl_clusters, child_dist_left, child_dist_right) = 
                  get_clusters_snarl(seeds, union_find_clusters,
                               cluster_dists, curr_snarl_children[curr_snarl_i],
                               node_to_seeds, distance_limit,
                               curr_snarl_i, 
                               rev_in_chain);
            last_rank = start_rank + 1;
            last_len = end_length;
           

            //Distance from the start of chain to the start of the current snarl
            int64_t add_dist_left = start_rank == 0 ? 0 : 
                                    chain_index.prefix_sum[start_rank] - 1;


             
            //Combine snarl clusters that can be reached by looping
            //TODO: Make the secondary snarl rank indicate whether it is reversed or not
            int64_t loop_dist_end = chain_index.loop_fd[start_rank + 1] - 1 ;
            int64_t loop_dist_start = chain_index.loop_rev[start_rank] - 1; 

            if (loop_dist_start != -1 || loop_dist_end != -1) {
                //If there is a loop in the chain, then it may
                //be possible for two clusters on the same snarl to
                //be combined
#ifdef DEBUG
cerr << "      Getting distances from loops in chain " << endl;
#endif
                vector<size_t> to_remove;
                size_t combined_left = -1;
                size_t combined_right = -1;
                for (size_t c : snarl_clusters) {
                    pair<int64_t, int64_t> dists = cluster_dists[c];
                    //Find new distances to ends using loops
                    int64_t new_left = dists.second == -1 || loop_dist_end == -1
                          ? -1
                          : dists.second + loop_dist_end + snarl_length 
                                         - end_length;
                    int64_t new_right=dists.first == -1 || loop_dist_start == -1
                                        ? -1
                                        : dists.first + loop_dist_start 
                                               + snarl_length - start_length;
                    //If this cluster was not combined with the another, 
                    //Update its left and right distances with loop dists
                    cluster_dists[c] = make_pair(
                        min_positive(dists.first, new_left),
                        min_positive(dists.second, new_right));
#ifdef DEBUG
cerr << "    Distances to ends of snarl including loops: " << cluster_dists[c].first << " " << cluster_dists[c].second << endl;
#endif


                    if (child_dist_left != -1 && loop_dist_start != -1 &&
                        dists.first != -1 && 
                        child_dist_left + dists.first 
                              + loop_dist_start - start_length 
                              < distance_limit) {  
                        //If this cluster can be combined with another cluster
                        //from the left

                        if (combined_left == -1) {
                            combined_left = c;
                        } else {
                            combined_left = union_find_clusters.find_group(
                                                                 combined_left);
                            pair<int64_t, int64_t>& old_dists =
                                                   cluster_dists[combined_left];
                            union_find_clusters.union_groups(combined_left, c);
                            size_t combined_group = 
                                              union_find_clusters.find_group(c);
                            if (combined_group != combined_left) {
                                to_remove.push_back(combined_left);
                            }
                            if (combined_group != c) {
                                to_remove.push_back(c);
                            }
                            combined_left = combined_group;
                            cluster_dists[combined_left] = make_pair(
                                  min_positive(old_dists.first,dists.first),
                                  min_positive(old_dists.second, dists.second));
                        }
                    }
                    if (child_dist_right != -1 && loop_dist_end != -1 &&
                        dists.second != -1 && 
                        child_dist_right + dists.second + loop_dist_end 
                                         - end_length < distance_limit) {  
                        //If this cluster can be combined with another cluster
                        //from the right 

                        if (combined_right == -1) {
                            combined_right = c;
                        } else {
                            combined_right = union_find_clusters.find_group(
                                                                combined_right);
                            union_find_clusters.union_groups(combined_right, c);
                            pair<int64_t, int64_t>& old_dists =
                                                  cluster_dists[combined_right];
                            size_t combined_group = 
                                              union_find_clusters.find_group(c);
                            if (combined_group != c) {
                                to_remove.push_back(c);
                            } 
                            if (combined_group != combined_right)  {
                                to_remove.push_back(combined_right);
                            }
                            combined_right = combined_group;
                            cluster_dists[combined_right] = make_pair(
                                  min_positive(old_dists.first, dists.first),
                                  min_positive(old_dists.second, dists.second));
                        }
                    }

                }
                for (size_t c : to_remove) {
                    snarl_clusters.erase(c);
                }
            }
 
#ifdef DEBUG
            cerr << "  Snarl distance limits: " << child_dist_left << " " << child_dist_right << endl;
            cerr << "  Snarl clusters to add: " << endl;
            for (size_t c : snarl_clusters) {
                pair<int64_t, int64_t> dists = cluster_dists[c];
                cerr << "\tleft: " << dists.first << " right : " << dists.second 
                     << endl;
                cerr << "\t\t";
                for (size_t x = 0 ; x < seeds.size() ; x++) {
                    if (union_find_clusters.find_group(x) == c) {
                        cerr << seeds[x] << " ";
                    }
                }
            }
            cerr << endl;

            cerr << "  Clusters on chain: " << endl;

            cerr << "  best left: " << best_left << " best right: " << best_right << endl;
            for (size_t c : chain_cluster_ids) {
                pair<int64_t, int64_t> dists = cluster_dists[c];
                cerr << "\tleft: " << dists.first << " right : " << dists.second 
                     << endl;
                cerr << "\t\t";
                for (size_t x = 0 ; x < seeds.size() ; x++) {
                    if (union_find_clusters.find_group(x) == c) {
                        cerr << seeds[x] << " ";
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
            int64_t old_left = child_dist_left;
            int64_t old_right = best_right;

            vector<size_t> to_add;//new cluster group ids 
            vector<size_t> to_erase; //old cluster group ids
            //New cluster- there will be at most one new cluster to add
            size_t combined_cluster = -1;
            int64_t combined_left = -1; int64_t combined_right = -1;
            best_left = -1; best_right = -1;
            for (size_t j : snarl_clusters) {
                /* For each of the clusters for the current snarl,
                 * find if it belongs to the new combined cluster*/ 

                pair<int64_t, int64_t> snarl_dists = move(cluster_dists[j]);

                if (old_right != -1 && snarl_dists.first != -1 &&
                    snarl_dists.first + snarl_lim_left < distance_limit) {
                    //If this snarl cluster is in the combined cluster
                    if (combined_cluster == -1) {
                        combined_cluster = j;
                        combined_left = snarl_dists.first == -1 ? -1 :
                                            snarl_dists.first + add_dist_left;
                        combined_right = snarl_dists.second;
                    } else {
                        union_find_clusters.union_groups(combined_cluster, j);
                        size_t new_group = union_find_clusters.find_group(j);
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

                    cluster_dists[j] = move(d); 
                }
            }
            for (size_t i : chain_cluster_ids) {
                //For each old chain cluster
                pair<int64_t, int64_t>& chain_dists = cluster_dists[i];

                if (old_left != -1 && chain_dists.second != -1 && 
                    chain_dists.second + chain_lim_right < distance_limit){
                    //If this chain cluster is in the combined cluster
                    if (combined_cluster == -1) {
                        combined_cluster = i;
                        combined_left = chain_dists.first;
                        combined_right = chain_dists.second + dist_to_end;
                    } else {
                        union_find_clusters.union_groups(combined_cluster, i);
                        size_t new_group = union_find_clusters.find_group(i);
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
                    if (chain_dists.first > distance_limit &&
                        chain_dists.second - end_length > distance_limit) {
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
                cluster_dists[combined_cluster] = 
                                      make_pair(combined_left, combined_right);
                best_left = min_positive(best_left, combined_left);
                best_right = min_positive(best_right, combined_right);
            }
                  
#ifdef DEBUG 
            cerr << "\t at snarl number " << curr_snarl_i << ", clusters:" <<endl;

            for (size_t c : chain_cluster_ids) {
                pair<int64_t, int64_t> dists = cluster_dists[c];
                cerr << "\t\tleft: " << dists.first << " right : " << dists.second << endl;
                cerr << "\t\t\t";
                for (size_t x = 0 ; x < seeds.size() ; x++) {
                    if (union_find_clusters.find_group(x) == c) {
                        cerr << seeds[x] << " ";
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
               int64_t d = cluster_dists[i].second;
               cluster_dists[i].second = d == -1 ? -1 : d + dist_to_end;
               best_right = min_positive(best_right, cluster_dists[i].second);
           }
        }


        if (chain_index.is_looping_chain) {
            //If the chain loops, then the clusters might be connected by 
            //looping around the chain
            //
            int64_t first_length = chain_index.prefix_sum[0];
            vector<size_t> to_erase; //old cluster group ids
            //New cluster- there will be at most one new cluster to add
            size_t combined_cluster = -1;

            for (size_t i : chain_cluster_ids) {
                //For each chain cluster
                pair<int64_t, int64_t>& chain_dists = cluster_dists[i];

                if ((chain_dists.second != -1 && best_left != -1 &&
                     chain_dists.second + best_left - first_length 
                                                          < distance_limit) ||
                   (chain_dists.first != -1 && best_right != -1 && 
                      chain_dists.first + best_right - first_length 
                                                          < distance_limit)){
                    //If this chain cluster is in the combined cluster
                    if (combined_cluster == -1) {
                        combined_cluster = i;
                    } else {
                        union_find_clusters.union_groups(combined_cluster, i);
                        size_t new_group = union_find_clusters.find_group(i);
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
            for (size_t x = 0 ; x < seeds.size() ; x++) {
                if (union_find_clusters.find_group(x) == c) {
                    cerr << seeds[x] << " ";
                }
            }
            cerr << endl;
        }
        bool got_left = false;
        bool got_right = false;
        for (size_t c : chain_cluster_ids) {
            pair<int64_t, int64_t> dists = cluster_dists[c];
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

            assert (group_id == union_find_clusters.find_group(group_id));
        }
#endif

        return make_tuple(move(chain_cluster_ids), best_left, best_right) ; 
    };



    SnarlSeedClusterer::child_cluster_t SnarlSeedClusterer::get_clusters_snarl(
                       const vector<pos_t>& seeds,
                       structures::UnionFind& union_find_clusters,
                       vector<pair<int64_t, int64_t>>& cluster_dists,
                       vector<pair<child_node_t, child_cluster_t>>& child_nodes,
                       hash_map<id_t, vector<size_t>>& node_to_seeds,
                       size_t distance_limit, 
                       size_t snarl_index_i, bool rev) {
        /*Get the clusters on this snarl. child_nodes stores each of the
         * nodes of the net graph that contain seeds and their clusters
         * Nodes have not yet been clustered */
        MinimumDistanceIndex::SnarlIndex& snarl_index = dist_index.snarl_indexes[
                                                    snarl_index_i];
#ifdef DEBUG 
        cerr << "Finding clusters on snarl number " << snarl_index_i << " headed by node " << snarl_index.id_in_parent << endl;
#endif
        int64_t start_length = snarl_index.nodeLength(0);
        int64_t end_length = snarl_index.nodeLength(snarl_index.num_nodes*2 -1);

        size_t num_children = child_nodes.size();

        //clusters of children, indices assume child node, snarl, and chain
        //vectors are contiguous 
        vector<vector<size_t>> child_clusters(num_children);

        //Return value- group ids of all clusters on this snarl
        hash_set<size_t> snarl_cluster_ids;
        int64_t best_left = -1;
        int64_t best_right = -1;
 
        auto combine_clusters = [&] (size_t& new_group, size_t& combined_group, 
                                    pair<int64_t, int64_t>& dists){
            //Helper function to combine clusters in two nodes of the same snarl
            if (combined_group == -1) {
                snarl_cluster_ids.insert(new_group);
                cluster_dists[new_group] = dists;
                combined_group = new_group;
            } else {

                combined_group = union_find_clusters.find_group(combined_group);
                pair<int64_t, int64_t>old_dists = cluster_dists[combined_group];
                union_find_clusters.union_groups(new_group, combined_group);
                size_t new_g = union_find_clusters.find_group(new_group);
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
                cluster_dists[new_g] = dists;
                new_group = new_g;
                combined_group = new_g;
            }
            return;
        };

        //Maps each cluster of child nodes to its left and right distances
        hash_map<size_t, pair<int64_t, int64_t>> old_dists;
        //Maps each child node to the left and right bounds of its clusters
        vector<pair<int64_t, int64_t>> dist_bounds(num_children, 
                                            make_pair(-1, -1));

        for (size_t i = 0; i < num_children ; i++) {
            //Go through each child node of the netgraph and find clusters

            auto& children = child_nodes [i];
            child_node_t& child = children.first;
            child_cluster_t& curr_child_clusters = children.second;

            // Note that child.second is the type, and child.first is the
            // *number* of the child in that type, *not* the heading node ID.
            // Ranks in parents are computed from node ID, so we have to get it.
            id_t child_node_id;
            switch (child.second) {
            case NODE:
                child_node_id = child.first;
                break;
            case SNARL:
                child_node_id = dist_index.snarl_indexes[child.first].id_in_parent;
                break;
            case CHAIN:
                child_node_id = dist_index.chain_indexes[child.first].id_in_parent;
                break;
            }
            
            //Rank of this node in the snarl
            //If this node is a snarl/chain, then this snarl will be the
            //secondary snarl
            size_t node_rank = child.second == NODE 
                    ? dist_index.primary_snarl_ranks[
                                  child_node_id - dist_index.min_node_id]-1
                    : dist_index.secondary_snarl_ranks[
                                 child_node_id - dist_index.min_node_id]-1;
            size_t rev_rank = node_rank % 2 == 0
                           ? node_rank + 1 : node_rank - 1;
            if ((child.second == SNARL && dist_index.snarl_indexes[
                          dist_index.primary_snarl_assignments[child_node_id-dist_index.min_node_id]-1].rev_in_parent) || 
                 (child.second == CHAIN && dist_index.chain_indexes[
                     dist_index.chain_assignments[child_node_id-dist_index.min_node_id]-1].rev_in_parent)) {
                //If this node (child snarl/chain) is reversed in the snarl
                //TODO: Make the secondary snarl rank indicate whether it is reversed or not
                size_t temp = node_rank;
                node_rank = rev_rank;
                rev_rank = temp;
            }

            //The clusters furthest to the left and right for this child node
            int64_t child_dist_left; int64_t child_dist_right;
            if (child.second == NODE) {
                //If this node is a node, we need to find the clusters
                int64_t node_len = snarl_index.nodeLength(node_rank);

                hash_set<size_t> c; 
                tie (c, child_dist_left, child_dist_right) = get_clusters_node(
                         seeds, union_find_clusters, cluster_dists, 
                         node_to_seeds[child_node_id], distance_limit, 
                         child_node_id, node_len);
                child_clusters[i].insert(child_clusters[i].end(),
                           make_move_iterator(c.begin()),
                           make_move_iterator(c.end())); 

            } else  {
                //If this node is a snarl or chain, the clusters have already
                //been found, we just need to put them in a vector

                hash_set<size_t>& c = get<0>(curr_child_clusters);
                child_dist_left = get<1>(curr_child_clusters);
                child_dist_right = get<2>(curr_child_clusters);
                child_clusters[i].insert(child_clusters[i].end(),
                                            make_move_iterator(c.begin()),
                                            make_move_iterator(c.end())); 
            }
            dist_bounds[i] = make_pair(child_dist_left, child_dist_right);

#ifdef DEBUG
            cerr << "Finding distances to parent snarl " << snarl_index_i << " ends from child " << i << "/" << num_children << endl;
            cerr << "Child is " << typeToString(child.second) << " number " << child.first << " headed by " << child_node_id << endl;
            cerr << "Node rank is " << node_rank << " fwd, " << rev_rank << " rev of " << snarl_index.num_nodes * 2 << endl;
            cerr << "Clusters at this child:" << endl; 
            for (size_t c : child_clusters[i]) {
                cerr << "\tdist left: " << cluster_dists[c].first << " dist right: " << cluster_dists[c].second << endl;
                cerr << "\t\t";
                for (size_t x = 0 ; x < seeds.size() ; x++) {
                    if (union_find_clusters.find_group(x) == c) {
                        cerr << seeds[x] << " ";
                    }
                }
                cerr << endl;
            }

            // Make sure the net graph node is actually in the net graph.
            assert(node_rank != numeric_limits<size_t>::max());
#endif

            vector<size_t>& children_i = child_clusters[i];
            for (size_t c_i = 0 ; c_i < children_i.size() ; c_i ++) {
                //for each cluster of child node i, find the distances to the
                //ends of the snarl

                size_t c = children_i[c_i];
            
                pair<int64_t, int64_t> dists_c= cluster_dists[c];
                old_dists[c] = dists_c;

                pair<int64_t, int64_t> new_dists = snarl_index.distToEnds(node_rank,
                                        dists_c.first,dists_c.second);
#ifdef DEBUG
cerr << "\tcluster: " << c_i << "dists to ends in snarl" << snarl_index.id_in_parent << " : " << new_dists.first << " " << new_dists.second << endl;
#endif
                /*
                if (node_rank == 0 || node_rank == 1){
                    //If this is the first node of the chain then the dist
                    //left is just the distance to the start of the snarl
                    new_dists.first = snarl_index.rev_in_parent
                                ? dists_c.second : dists_c.first ;
                } else if (node_rank == snarl_index.num_nodes * 2-1 ||
                    node_rank == snarl_index.num_nodes * 2 - 2) {
                    new_dists.second = snarl_index.rev_in_parent
                           ? dists_c.first : dists_c.second ;
                } 
                */

                best_left = min_positive(best_left, new_dists.first);
                best_right = min_positive(best_right, new_dists.second);


                snarl_cluster_ids.insert(c);
                cluster_dists[c] = new_dists;
            }
            

            for (size_t j = 0 ; j <= i ; j++){
                //Go through other child net graph nodes up to and including i
                auto& other_node = child_nodes[j].first;

                // Note that other_node.second is the type, and other_node.first is the
                // *number* of the other_node in that type, *not* the heading node ID.
                // Ranks in parents are computed from node ID, so we have to get it.
                id_t other_node_id;
                switch (other_node.second) {
                case NODE:
                    other_node_id = other_node.first;
                    break;
                case SNARL:
                    other_node_id = dist_index.snarl_indexes[other_node.first].id_in_parent;
                    break;
                case CHAIN:
                    other_node_id = dist_index.chain_indexes[other_node.first].id_in_parent;
                    break;
                }
                
#ifdef DEBUG
                cerr << "Other net graph node is " << typeToString(other_node.second) << " number "
                    << other_node.first << " headed by node " << other_node_id
                    << " with dists: " << dist_bounds[j].first << " " << dist_bounds[j].second << endl;
                    
                    
#endif

                //Rank of this node in the snarl
                size_t other_rank = other_node.second == NODE ? 
                       dist_index.primary_snarl_ranks[
                          other_node_id - dist_index.min_node_id]-1
                    : dist_index.secondary_snarl_ranks[
                           other_node_id - dist_index.min_node_id]-1;
                size_t other_rev = other_rank % 2 == 0
                                    ? other_rank + 1 : other_rank - 1;

            if ((other_node.second == SNARL && dist_index.snarl_indexes[
                          dist_index.primary_snarl_assignments[other_node_id-dist_index.min_node_id]-1].rev_in_parent) || 
                 (other_node.second == CHAIN && dist_index.chain_indexes[
                     dist_index.chain_assignments[other_node_id-dist_index.min_node_id]-1].rev_in_parent)) {
                //If this node (child snarl/chain) is reversed in the snarl
                //TODO: Make the secondary snarl rank indicate whether it is reversed or not
                size_t temp = other_rank;
                other_rank = other_rev;
                other_rev = temp;
            }
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

                pair<int64_t, int64_t> dist_bounds_j = dist_bounds[j];

                if (max({dist_l_l, dist_l_r, dist_r_l, dist_r_r}) != -1
                   && MinimumDistanceIndex::minPos({dist_l_l, dist_l_r, 
                                           dist_r_l, dist_r_r}) < distance_limit
                   && min_positive(child_dist_left, child_dist_right) 
                                                             < distance_limit) {
                    //If the two nodes are reachable
                    for (size_t c_i = 0 ; c_i < children_i.size() ; c_i ++) {
                        //for each cluster of child node i

                        size_t c = children_i[c_i];
                        size_t c_group = union_find_clusters.find_group(c);

                        pair<int64_t, int64_t> new_dists;
                        pair<int64_t, int64_t> dists_c;

                         dists_c = old_dists[c];
                         new_dists = cluster_dists[c_group];
                        

                        if (dist_l_l != -1 && dists_c.first != -1 
                                 && dist_bounds_j.first != -1  
                                 && dist_l_l + dists_c.first 
                                      + dist_bounds_j.first < distance_limit) {
                            //If cluster c can be combined with clusters in j 
                            //from the left of both of them
                            combine_clusters(c_group, group_l_l, new_dists);
                        }
                        if (dist_l_r != -1 && dists_c.first != -1 
                            && dist_bounds_j.second != -1 
                            && dist_l_r + dists_c.first + dist_bounds_j.second 
                                                             < distance_limit) {
                            combine_clusters(c_group, group_l_r, new_dists);

                        }
                        if (dist_r_l != -1 && dists_c.second != -1 
                            && dist_bounds_j.first != -1 
                            && dist_r_l + dists_c.second + dist_bounds_j.first 
                                                             < distance_limit) {
                            combine_clusters(c_group, group_r_l, new_dists);
                        }
                        if (dist_r_r != -1 && dists_c.second != -1 
                            && dist_bounds_j.second != -1 
                            && dist_r_r + dists_c.second + dist_bounds_j.second 
                                                             < distance_limit) {
                            combine_clusters(c_group, group_r_r, new_dists);
                        }

                    }
                    //Go through children of j
                    vector<size_t>& children_j =  child_clusters[j];

                    for (size_t k_i = 0 ; k_i < children_j.size() ; k_i++){
                        size_t k = children_j[k_i];
                        //For each cluster of child j, find which overlaps with
                        //clusters of i
                        //k will already be part of a cluster in 
                        //snarl_cluster_ids but since we need to know the node 
                        //that the snarl is on we can't just loop through 
                        //snarl_cluster_ids
                        pair<int64_t, int64_t>& dist_bounds_k = old_dists[k];
                        size_t k_group = union_find_clusters.find_group(k);
                        pair<int64_t, int64_t> dists_k = cluster_dists[k_group];
                    

                        if (dist_l_l != -1 && child_dist_left != -1 
                             && dist_bounds_k.first != -1 
                             && dist_l_l + child_dist_left + dist_bounds_k.first
                                                             < distance_limit){

                            combine_clusters(k_group, group_l_l, dists_k);
                        }
                        if (dist_l_r != -1 && child_dist_left != -1 
                             && dist_bounds_k.second != -1  
                            && dist_l_r + child_dist_left + dist_bounds_k.second
                                                            < distance_limit ) {

                            combine_clusters(k_group, group_l_r, dists_k);
                        }
                        if (dist_r_l != -1 && child_dist_right != -1 
                            && dist_bounds_k.first != -1  
                            && dist_r_l + child_dist_right + dist_bounds_k.first
                                                             < distance_limit) {

                            combine_clusters(k_group, group_r_l, dists_k);
                        }
                        if (dist_r_r != -1 && child_dist_right != -1 
                           && dist_bounds_k.second != -1
                           && dist_r_r + child_dist_right + dist_bounds_k.second
                                                             < distance_limit) {

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
            pair<int64_t, int64_t> dists = cluster_dists[c];
            if (dists.first == best_left) {got_left = true;}
            if (dists.second == best_right) {got_right = true;}
            cerr << "\t" << c << ": left: " << dists.first << " right : " 
                 << dists.second << endl;
            cerr << "\t\t";
            for (size_t x = 0 ; x < seeds.size() ; x++) {
                if (union_find_clusters.find_group(x) == c) {
                    cerr << seeds[x] << " ";
                }
            }
            cerr << endl;
        }
        assert(got_left);
        assert(got_right);
        for (size_t group_id : snarl_cluster_ids) {
            assert (group_id == union_find_clusters.find_group(group_id));
        }
#endif
        return make_tuple(move(snarl_cluster_ids), best_left, best_right);
    };
}
