#include "seed_clusterer.hpp"

//#define DEBUG 

namespace vg {

    SnarlSeedClusterer::SnarlSeedClusterer() {
    };
        
    vector<vector<size_t>> SnarlSeedClusterer::cluster_seeds ( 
                  vector<pos_t> seeds, size_t distance_limit,
                  SnarlManager& snarl_manager, DistanceIndex& dist_index ){
        /* Given a vector of seeds and a limit, find a clustering of seeds where
         * seeds that are closer than the limit cluster together. 
         * Returns a vector of cluster assignments 
         */ 

        //Create a union find structure to hold the cluster assignments of
        //each seed. Each seed is initially its own cluster 
        structures::UnionFind union_find_clusters (seeds.size());


        //Map each node to a vector of seed indices
        hash_map<id_t, vector<size_t>> node_to_seeds;

        //For each level of the snarl tree, maps snarls at that level to 
        //node + node lengths belonging to the snarl
        vector<hash_map<const Snarl*, 
                        vector<pair<child_node_t, child_cluster_t>>>> 
                                                                 snarl_to_nodes;
        //TODO: reserve up to the dept of the snarl tree- also need to find the depth of the snarl tree in the SNarlManager

        //Populate node_to_seed and snarl_to_nodes
        get_nodes(seeds, snarl_manager, dist_index, node_to_seeds, 
                  snarl_to_nodes);

        //Maps each cluster group ID to the left and right distances
        vector<pair<int64_t, int64_t>> cluster_dists (seeds.size(), 
                                                      make_pair(-1, -1));

        int tree_depth = snarl_to_nodes.size()-1; 

        //Maps snarls in the current level to be clustered to its children
        hash_map<const Snarl*, vector<pair<child_node_t, child_cluster_t>>> 
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
            hash_map<const Snarl*, vector<pair<child_node_t, child_cluster_t>>>
                                                          parent_snarl_children;
            if (depth != 0) {
                parent_snarl_children = move(snarl_to_nodes[depth - 1]);
            }

            //Maps each chain to the snarls that comprise it
            //Snarls are unordered in the vector
            hash_map<const Chain*, vector<pair<size_t, 
              pair<const Snarl*, DistanceIndex::SnarlIndex*>>>> chain_to_snarls;


            for (auto it : curr_snarl_children){
                //Go through each of the snarls at this level, cluster them,
                //and find which chains they belong to, if any
                const Snarl* snarl = it.first;
                const Snarl* parent_snarl = snarl_manager.parent_of(snarl);

#ifdef DEBUG
cerr << "At depth " << depth << " snarl " << snarl->start() << "with children "
<< endl;
for (auto it2 : it.second) {
    auto type = "";
    switch(it2.first.second) {case 0 : type = "chain";
                                break;
                       case 1 : type = "snarl";
                                break;
                       case 2 : type = "node";
                                 break;};
    cerr << " " << type << ": " << it2.first.first.first << " " << it2.first.first.second << endl; 
}
#endif
                DistanceIndex::SnarlIndex& snarl_index = 
                            dist_index.snarlDistances.at(
                               make_pair(snarl->start().node_id(),
                                         snarl->start().backward()));

                if (snarl_manager.in_nontrivial_chain(snarl)){
    
                    //If this snarl is in a chain, add it to chains
                    const Chain* chain = snarl_manager.chain_of(snarl);
                    //TODO: COuld lookup chain index less
                    DistanceIndex::ChainIndex& chain_index = 
                                   dist_index.chainDistances.at(
                                                get_start_of(*chain).node_id());

                    DistanceIndex::SnarlIndex& snarl_index = 
                                    dist_index.snarlDistances.at(
                                       make_pair(snarl->start().node_id(),
                                                 snarl->start().backward()));

                    bool rev_in_chain=snarl_manager.chain_orientation_of(snarl);
                    id_t start_node = rev_in_chain ? snarl->end().node_id() :
                                                     snarl->start().node_id(); 
                    size_t rank = chain_index.snarlToIndex[start_node];

                    chain_to_snarls[chain].push_back( 
                               make_pair(rank, make_pair(snarl, &snarl_index)));
                    
                } else {
                    //If this snarl is not in a chain, add it as a child of the
                    //parent snarl for the next level

                    if (depth != 0 && parent_snarl != nullptr){
                        parent_snarl_children[parent_snarl].push_back(
                            make_pair(make_pair(
                                  make_pair(snarl->start().node_id(),
                                            snarl->start().backward()), SNARL), 
                             get_clusters_snarl(seeds, union_find_clusters,
                                    cluster_dists, it.second, 
                                    node_to_seeds, distance_limit,
                                    snarl_manager, snarl_index, snarl, false)));
                    } else {
                         get_clusters_snarl(seeds, union_find_clusters,
                                    cluster_dists, it.second, 
                                    node_to_seeds, distance_limit,
                                    snarl_manager, snarl_index, snarl, false);
                    }
                }
            }
            for (auto& kv : chain_to_snarls) {
                //For each chain in this level, find the clusters 

                const Chain* chain = kv.first;
                auto chain_start = get_start_of(*chain);
                const Snarl* parent_snarl = snarl_manager.parent_of(
                              kv.second[0].second.first);
                parent_snarl_children[parent_snarl].push_back(
                  make_pair(make_pair(
                                  make_pair(get_start_of(*chain).node_id(),
                                          get_start_of(*chain).backward()),
                                          CHAIN),
                        get_clusters_chain( seeds, union_find_clusters,
                             cluster_dists, kv.second, curr_snarl_children,
                             node_to_seeds, distance_limit, snarl_manager,
                             dist_index, chain)));
            }

            curr_snarl_children= move(parent_snarl_children);
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


    void SnarlSeedClusterer::get_nodes( 
              const vector<pos_t>& seeds,  const SnarlManager& snarl_manager, 
              DistanceIndex& dist_index, 
              hash_map<id_t, vector<size_t>>& node_to_seeds,
              vector<hash_map<const Snarl*, 
                     vector<pair<child_node_t, child_cluster_t>>>>& snarl_to_nodes) {

        /* Given a snarl tree and seeds, find the nodes containing seeds and
         * assign each node to a level in the snarl tree*/ 

        hash_map<id_t, const Snarl*> node_to_snarl;
        for (size_t i = 0; i < seeds.size(); i++) {
            //For each seed, assign it to a node and the node to a snarl 

            pos_t pos = seeds[i];
            id_t id = get_id(pos);
            const Snarl* snarl;
            auto found_snarl = node_to_snarl.find(id);
            if (found_snarl != node_to_snarl.end()) {
                snarl = found_snarl->second;
            } else {
                snarl = dist_index.snarlOf(id);

                if (snarl_manager.in_nontrivial_chain(snarl)) {
                    //If the snarl is in a chain, find which snarl the node will
                    //be assigned to
                    const Chain* chain = snarl_manager.chain_of(snarl);
    
                    bool rev_in_chain=snarl_manager.chain_orientation_of(snarl);
                    id_t chain_start = get_start_of(*chain).node_id();
    
                    //A boundary node on a chain is assigned to the snarl 
                    //preceding it in the chain
                    //First node in the chain is assigned to the first snarl 
                    if ((id == snarl->start().node_id() && !rev_in_chain) ||
                        (id == snarl->end().node_id() && rev_in_chain)){ 
                        //If seed is on first boundary node of snarl 
                        //relative to traversal in the chain
                        if (id != chain_start) {
                            //Unless this is the first node in the chain,
                            //switch to the other snarl sharing this node
                            snarl = move(snarl_manager.into_which_snarl(
                                       snarl->start().node_id(),
                                      !snarl->start().backward()));
                        }
                    }
                }
                node_to_snarl.emplace(id, snarl);
            }

            size_t depth = snarl_manager.get_depth(snarl);

            auto s = node_to_seeds.find(id);
            if (s != node_to_seeds.end()) {
                s->second.push_back(i);
            } else {

                //Map node to seed
                node_to_seeds.emplace(id, vector<size_t>({i}));
                int64_t node_length = dist_index.snarlDistances.at(
                          make_pair(snarl->start().node_id(),
                                    snarl->start().backward())).nodeLength(id);

                //Map snarl to node
                if (snarl_to_nodes.size() < depth+1) {
                    snarl_to_nodes.resize(depth+1); 
                }
                child_cluster_t empty;
                snarl_to_nodes[depth][snarl].emplace_back(
                                  make_pair(make_pair(id,false), NODE), empty);
            }
        }
    }


    SnarlSeedClusterer::child_cluster_t SnarlSeedClusterer::get_clusters_node(
                       const vector<pos_t>& seeds,
                       structures::UnionFind& union_find_clusters, 
                       vector< pair<int64_t, int64_t>>& cluster_dists,
                       vector<size_t>& seed_indices,
                       size_t distance_limit, const SnarlManager& snarl_manager,
                       id_t root, int64_t node_length) {
#ifdef DEBUG 
cerr << "Finding clusters on node " << root << " Which has length " <<
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

                best_dist_left = DistanceIndex::minPos({
                                        dist_left, best_dist_left});
                best_dist_right = DistanceIndex::minPos({
                                        dist_right, best_dist_right});

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
                                               dists.second << "\t seeds: ";
    vector<size_t> seed_is = union_find_clusters.group(c);
    for (size_t s : seed_is) {
        cerr << seeds[s] << " ";
    }
    cerr << endl;
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

            best_left = DistanceIndex::minPos({offset,  best_left});
            best_right = DistanceIndex::minPos({node_length-offset+1,best_right});

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
                last_left = DistanceIndex::minPos({last_left, offset});
                last_right = DistanceIndex::minPos({last_right, 
                                                    node_length-offset+1});
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
         << dists.second << "\t seeds: ";
    vector<size_t> seed_is = union_find_clusters.group(c);
    for (size_t s : seed_is) {
        cerr << seeds[s] << " ";
    }
    cerr << endl;
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
                       vector<  pair<size_t, 
                            pair<const Snarl*, DistanceIndex::SnarlIndex*>>>& 
                                                                snarls_in_chain,
                       hash_map<const Snarl*, 
                                 vector<pair<child_node_t, child_cluster_t>>>& 
                                                            curr_snarl_children,
                       hash_map<id_t, vector<size_t>>& node_to_seeds,
                       size_t distance_limit, const SnarlManager& snarl_manager,
                       DistanceIndex& dist_index, const Chain* root) {
        /*Find all the clusters in the given chain */
#ifdef DEBUG 
cerr << "Finding clusters on chain " << get_start_of(*root).node_id() << endl;
#endif
 
        //Sort vector of snarls by rank and remove duplicates
        std::sort(snarls_in_chain.begin(), snarls_in_chain.end(), 
                     [&](const auto s, const auto t) -> bool {
                          return s.first < t.first; 
                      } );
  
        DistanceIndex::ChainIndex& chain_index = dist_index.chainDistances.at(
                                                 get_start_of(*root).node_id());
        hash_set<size_t> chain_cluster_ids;

        int64_t best_left = -1;
        int64_t best_right = -1;
        ChainIterator chain_e = chain_end(*root);
        //The node at which the chain clusters reach
        id_t last_snarl = get_start_of(*root).node_id();
        int64_t last_len = 0;
        id_t start_node;
        id_t end_node;
        const Snarl* prev_snarl = nullptr;

        for (auto x : snarls_in_chain) {
            /* For each child snarl in the chain, find the clusters of just the
             * snarl, and progressively build up clusters spanning up to that 
             * snarl
             * Snarls are in the order that they are traversed in the chain
             */
            const Snarl* curr_snarl = x.second.first;
            if (curr_snarl == prev_snarl) {
                continue;
            } else {
                prev_snarl = curr_snarl;
            }
            DistanceIndex::SnarlIndex& snarl_index = *x.second.second;

            bool rev_in_chain = snarl_manager.chain_orientation_of( curr_snarl);


            //Start and end node of snarl relative to chain
            id_t start_node = rev_in_chain ? curr_snarl->end().node_id() :
                                                 curr_snarl->start().node_id();
            end_node = rev_in_chain ? curr_snarl->start().node_id() :
                                                 curr_snarl->end().node_id();
            int64_t start_length = snarl_index.nodeLength(start_node);
            int64_t end_length = snarl_index.nodeLength(end_node);

            //Distance from right end of chain clusters to the end of current
            //snarl cluster
            int64_t dist_to_end = snarl_index.snarlLength() - start_length;


            if (last_snarl != start_node) { 
                /* If the chain clusters don't reach this snarl,
                 * extend their dist_right to the beginning of this snarl
                 */  
                int64_t offset = chain_index.prefixSum[
                                   chain_index.snarlToIndex[start_node]] -
                                chain_index.prefixSum[
                                      chain_index.snarlToIndex[last_snarl]] +
                                start_length - last_len;

                for (size_t i : chain_cluster_ids) {
                    cluster_dists[i].second = cluster_dists[i].second == -1 
                                        ? -1 : cluster_dists[i].second + offset;
                }
                best_right = best_right == -1 ? -1 : best_right + offset;
            }

            //Find the clusters of the current snarl
            hash_set<size_t> snarl_clusters; 
            int64_t child_dist_left; int64_t child_dist_right; 
            tie (snarl_clusters, child_dist_left, child_dist_right) = 
                  get_clusters_snarl(seeds, union_find_clusters,
                                 cluster_dists, curr_snarl_children[curr_snarl],
                                 node_to_seeds, distance_limit,
                                 snarl_manager, snarl_index, curr_snarl, 
                                 rev_in_chain);
            last_snarl = end_node;
            last_len = end_length;
           

            //Distance from the start of chain to the start of the current snarl
            int64_t add_dist_left = chain_index.prefixSum[
                                      chain_index.snarlToIndex[start_node]] - 1;


             
            //Combine snarl clusters that can be reached by looping
            int64_t loop_dist_end = chain_index.loopFd[
                                      chain_index.snarlToIndex[end_node]] - 1 ;
            int64_t loop_dist_start = chain_index.loopRev[
                                     chain_index.snarlToIndex[start_node]] - 1; 

            if (loop_dist_start != -1 || loop_dist_end != -1) {
                //If there is a loop in the chain, then it may
                //be possible for two clusters on the same snarl to
                //be combined
                vector<size_t> to_remove;
                size_t combined_left = -1;
                size_t combined_right = -1;
                for (size_t c : snarl_clusters) {
                    pair<int64_t, int64_t> dists = cluster_dists[c];
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
                            if (combined_group == c) {
                                to_remove.push_back(combined_left);
                            } else {
                                to_remove.push_back(c);
                            }
                            combined_left = combined_group;
                            cluster_dists[combined_left] = make_pair(
                                      DistanceIndex::minPos({old_dists.first, 
                                                             dists.first}),
                                      DistanceIndex::minPos({old_dists.second, 
                                                             dists.second}));
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
                            if (combined_group == c) {
                                to_remove.push_back(combined_right);
                            } else {
                                to_remove.push_back(c);
                            }
                            combined_right = combined_group;
                            cluster_dists[combined_right] = make_pair(
                                      DistanceIndex::minPos({old_dists.first, 
                                                             dists.first}),
                                      DistanceIndex::minPos({old_dists.second, 
                                                             dists.second}));
                        }
                    }

                }
                for (size_t c : to_remove) {
                    snarl_clusters.erase(c);
                }
            }
 
#ifdef DEBUG
cerr << "  Clusters on chain: " << endl;

cerr << "  best left: " << best_left << " best right: " << best_right << endl;
for (size_t c : chain_cluster_ids) {
    pair<int64_t, int64_t> dists = cluster_dists[c];
    cerr << "\tleft: " << dists.first << " right : " << dists.second 
         << "\t seeds: ";
    vector<size_t> ss = union_find_clusters.group(c);
    for (size_t s : ss) {
        cerr << seeds[s] << " ";
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
                        combined_left = DistanceIndex::minPos({combined_left,
                                            snarl_dists.first == -1 ? -1 :
                                            snarl_dists.first + add_dist_left});
                        combined_right = DistanceIndex::minPos({combined_right,
                                                           snarl_dists.second});
                    }
                        
                } else {
                    //This snarl becomes a new chain cluster
                    to_add.push_back(j);
                    pair<int64_t, int64_t> d = make_pair(snarl_dists.first == -1
                                      ? -1 : snarl_dists.first + add_dist_left, 
                                                snarl_dists.second);
                    best_left = DistanceIndex::minPos({best_left, d.first}); 
                    best_right = DistanceIndex::minPos({best_right, d.second});

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
                        combined_left = DistanceIndex::minPos({combined_left,
                                                            chain_dists.first});
                        combined_right = DistanceIndex::minPos({combined_right,
                                             chain_dists.second + dist_to_end});
                    }
                } else {
                    //If this chain cluster is on its own, extend its right 
                    //distance to the end of the current snarl
                    chain_dists.second += dist_to_end;
                    if (chain_dists.first > distance_limit &&
                        chain_dists.second > distance_limit) {
                        //If the distance from the seeds in this cluster to
                        //either end of the chain is greater than the distance
                        //limit, then it cannot cluster with anything else
                        //so we can stop keeping track of it
                        to_erase.push_back(i);
                    } else {
                        best_left = DistanceIndex::minPos({best_left, 
                                                           chain_dists.first});
                        best_right = DistanceIndex::minPos({best_right, 
                                                        chain_dists.second});
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
                best_left = DistanceIndex::minPos({best_left, combined_left});
                best_right =DistanceIndex::minPos({best_right, combined_right});
            }
                  
#ifdef DEBUG 
cerr << "\t at snarl " << curr_snarl->start().node_id() << ", clusters:" <<endl;

for (size_t c : chain_cluster_ids) {
    pair<int64_t, int64_t> dists = cluster_dists[c];
    cerr << "\tleft: " << dists.first << " right : " << dists.second 
    << "\t seeds: ";
    vector<size_t> ss = union_find_clusters.group(c);
    for (size_t s : ss) {
        cerr << seeds[s] << " ";
    }
    
    cerr << endl;
}
#endif
                
        }
         

        if (last_snarl != get_end_of(*root).node_id()) {
            best_right = -1;
           //Extend the right bound of each cluster to the end of the chain
           int64_t dist_to_end = chain_index.chainLength()
                       - chain_index.prefixSum[
                                      chain_index.snarlToIndex[last_snarl]] + 1 
                                - last_len;
           for (size_t i : chain_cluster_ids) {
               int64_t& d = cluster_dists[i].second;
               d += dist_to_end;
               best_right = DistanceIndex::minPos({best_right, d});
           }
        }

#ifdef DEBUG 
cerr << "Found clusters on chain " << get_start_of(*root).node_id() << endl;
cerr << "best left : " << best_left << " best right : " << best_right << endl;
bool got_left = false;
bool got_right = false;
for (size_t c : chain_cluster_ids) {
    pair<int64_t, int64_t> dists = cluster_dists[c];
    assert(dists.first == -1 || dists.first >= best_left);
    assert(dists.second == -1 || dists.second >= best_right);
    if (dists.first == best_left) {got_left = true;}
    if (dists.second == best_right) {got_right = true;}
    cerr << "\t" << c << ": left: " << dists.first << " right : " 
         << dists.second << "\t seeds: ";
    vector<size_t> seed_is = union_find_clusters.group(c);
    for (size_t s : seed_is) {
        cerr << seeds[s] << " ";
    }
    cerr << endl;
}
assert(got_left);
assert(got_right);
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
                       size_t distance_limit, const SnarlManager& snarl_manager,
                       DistanceIndex::SnarlIndex& snarl_index, 
                       const Snarl* root, bool rev) {
        /*Get the clusters on this snarl. snarl_children stores the children of
         * this snarl, 0 for chain, 1 for snarl, 2 for node. child_clusters 
         * stores the output of clustering for snarls and chains. 
         * Nodes have not yet been clustered */
        #ifdef DEBUG 
        cerr << "Finding clusters on snarl " << root->start() << endl;
        #endif
        int64_t start_length = snarl_index.nodeLength(
                                                  snarl_index.snarlStart.first);
        int64_t end_length = snarl_index.nodeLength(snarl_index.snarlEnd.first);

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
                       DistanceIndex::minPos({dists.first, old_dists.first}),
                       DistanceIndex::minPos({dists.second, old_dists.second}));
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
            pair<pair<id_t, bool>, ChildNodeType>& child = children.first;
            child_cluster_t& curr_child_clusters = children.second;
            id_t curr_node;
            bool curr_rev;
            tie (curr_node, curr_rev) = child.first;
            //The clusters furthest to the left and right for this child node
            int64_t child_dist_left; int64_t child_dist_right;
            if (child.second == NODE) {
                //If this node is a node, we need to find the clusters
                int64_t node_len = snarl_index.nodeLength(curr_node);

                hash_set<size_t> c; 
                tie (c, child_dist_left, child_dist_right) = get_clusters_node(
                         seeds, union_find_clusters, cluster_dists, 
                         node_to_seeds[curr_node], distance_limit, 
                         snarl_manager, curr_node, node_len);
                child_clusters[i].insert(child_clusters[i].end(),
                           make_move_iterator(c.begin()),
                           make_move_iterator(c.end())); 

            } else  {
                //If this node is a snarl or chain, the clusters have already
                //been found, we just need to put them in a vector

                const Snarl* curr_snarl = snarl_manager.into_which_snarl(
                                                           curr_node, curr_rev);
                hash_set<size_t>& c = get<0>(curr_child_clusters);
                child_dist_left = get<1>(curr_child_clusters);
                child_dist_right = get<2>(curr_child_clusters);
                child_clusters[i].insert(child_clusters[i].end(),
                                            make_move_iterator(c.begin()),
                                            make_move_iterator(c.end())); 
            }
            dist_bounds[i] = make_pair(child_dist_left, child_dist_right);

            //For child node i, find distances to ends of root snarl
            int64_t dist_s_l = snarl_index.snarlDistance(
                       snarl_index.snarlStart, make_pair(curr_node, curr_rev));
            int64_t dist_s_r = snarl_index.snarlDistance(
                      snarl_index.snarlStart, make_pair(curr_node, !curr_rev));
            int64_t dist_e_l = snarl_index.snarlDistance(
                             make_pair(snarl_index.snarlEnd.first, 
                                       !snarl_index.snarlEnd.second),
                                      make_pair(curr_node, curr_rev));
             int64_t dist_e_r = snarl_index.snarlDistance(
                                    make_pair(snarl_index.snarlEnd.first, 
                                              !snarl_index.snarlEnd.second),
                                         make_pair(curr_node, !curr_rev));



            vector<size_t>& children_i = child_clusters[i];
            for (size_t c_i = 0 ; c_i < children_i.size() ; c_i ++) {
                //for each cluster of child node i, find the distances to the
                //ends of the snarl

                size_t c = children_i[c_i];
            
                pair<int64_t, int64_t> dists_c= cluster_dists[c];
                old_dists[c] = dists_c;

    
                //find distances to ends of snarl for this cluster 
                int64_t new_dist_s_l = dist_s_l == -1  ||  
                                       dists_c.first == -1 
                            ? -1 : dist_s_l + dists_c.first;
                int64_t new_dist_s_r = dist_s_r == -1 ||  
                                       dists_c.second == -1 
                            ? -1 : dist_s_r + dists_c.second;
                int64_t new_dist_e_l = dist_e_l == -1 || 
                                       dists_c.first == -1
                            ? -1 : dist_e_l + dists_c.first;
                int64_t new_dist_e_r = dist_e_r == -1 || 
                                       dists_c.second == -1
                            ? -1 : dist_e_r + dists_c.second;
      
                pair<int64_t, int64_t> new_dists = make_pair (
                         DistanceIndex::minPos({new_dist_s_l,  
                                                new_dist_s_r}), 
                         DistanceIndex::minPos({new_dist_e_l,
                                                new_dist_e_r}));
                if (curr_node == snarl_index.snarlStart.first){
                    new_dists.first = snarl_index.snarlStart.second
                                ? dists_c.second : dists_c.first ;
                }
                if (curr_node == snarl_index.snarlEnd.first) {
                    new_dists.second = snarl_index.snarlEnd.second 
                           ? dists_c.first : dists_c.second ;
                } 
                best_left = DistanceIndex::minPos({best_left, 
                                                  new_dists.first});
                best_right = DistanceIndex::minPos({best_right, 
                                                new_dists.second});


                snarl_cluster_ids.insert(c);
                cluster_dists[c] = new_dists;
            }
            

            for (size_t j = 0 ; j <= i ; j++){
                //Go through child nodes up to and including i
                id_t other_node;
                bool other_rev;
                tie (other_node, other_rev) = child_nodes[j].first.first;
                //Find distance from each end of current node to 
                //each end of other node
                int64_t dist_l_l = snarl_index.snarlDistanceShort(
                         make_pair(curr_node, !curr_rev), 
                         make_pair(other_node, other_rev));
                int64_t dist_l_r = snarl_index.snarlDistanceShort(
                         make_pair(curr_node, !curr_rev), 
                         make_pair(other_node, !other_rev));
                int64_t dist_r_l = snarl_index.snarlDistanceShort(
                          make_pair(curr_node, curr_rev), 
                          make_pair(other_node, other_rev));
                int64_t dist_r_r = snarl_index.snarlDistanceShort(
                          make_pair(curr_node, curr_rev),
                          make_pair(other_node, !other_rev));

                //group ids of clusters combined between node i left and 
                //node j left, etc
                size_t group_l_l = -1;
                size_t group_l_r = -1;
                size_t group_r_l = -1;
                size_t group_r_r = -1;

                pair<int64_t, int64_t> dist_bounds_j = dist_bounds[j];

                if (DistanceIndex::minPos({dist_bounds_j.first, 
                                           dist_bounds_j.second})
                                        < distance_limit) {
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
                }
                if (max({dist_l_l, dist_l_r, dist_r_l, dist_r_r}) != -1
                   && DistanceIndex::minPos({dist_l_l, dist_l_r, 
                                           dist_r_l, dist_r_r}) < distance_limit
                   && DistanceIndex::minPos({child_dist_left, 
                                          child_dist_right}) < distance_limit) {
                    //If the two nodes are reachable
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

#ifdef DEBUG 
cerr << "At i node " << curr_node << " j node " << other_node << endl;
cerr << "    with best left and right values: " << best_left << " " 
     << best_right << endl;
bool got_left = false;
bool got_right = false;
for (size_t c : snarl_cluster_ids) {
    pair<int64_t, int64_t> dists = cluster_dists[c];
    if (dists.first == best_left) {got_left = true;}
    if (dists.second == best_right) {got_right = true;}
    cerr << "\t" << c << ": left: " << dists.first << " right : " 
         << dists.second << "\t seeds: ";
    vector<size_t> seed_is = union_find_clusters.group(c);
    assert(dists.first == -1 || dists.first >= best_left);
    assert(dists.second == -1 || dists.second >= best_right);
    for (size_t s : seed_is) {
        cerr << seeds[s] << " ";
    }
    cerr << endl;
}
assert(got_left);
assert(got_right);
for (size_t group_id : snarl_cluster_ids) {

assert (group_id == union_find_clusters.find_group(group_id));
}
#endif
            }
        }
#ifdef DEBUG 
cerr << "Found clusters on snarl " << root->start() << endl;
cerr << "    with best left and right values: " << best_left << " " 
     << best_right << endl;
bool got_left = false;
bool got_right = false;
for (size_t c : snarl_cluster_ids) {
    pair<int64_t, int64_t> dists = cluster_dists[c];
    if (dists.first == best_left) {got_left = true;}
    if (dists.second == best_right) {got_right = true;}
    cerr << "\t" << c << ": left: " << dists.first << " right : " 
         << dists.second << "\t seeds: ";
    vector<size_t> seed_is = union_find_clusters.group(c);
    assert(dists.first == -1 || dists.first >= best_left);
    assert(dists.second == -1 || dists.second >= best_right);
    for (size_t s : seed_is) {
        cerr << seeds[s] << " ";
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
