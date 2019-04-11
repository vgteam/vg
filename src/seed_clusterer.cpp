#include "seed_clusterer.hpp"

//#define DEBUG 

namespace vg {

    SnarlSeedClusterer::SnarlSeedClusterer() {
    };


    void SnarlSeedClusterer::seed2subtree( 
               const vector<pos_t>& seeds,  const SnarlManager& snarl_manager, 
               DistanceIndex& dist_index, chains_to_snarl_t& chains_to_snarl,
               snarls_to_node_t& snarls_to_node, node_to_seed_t& node_to_seed) {

       
        /* Given a snarl tree and seeds, find the subtree of the snarl tree
         * that contains seeds */ 
        //Maps chains to the snarls containing seeds
        hash_map< const Chain*, map<size_t,const Snarl*>> chains_to_snarl_rank;
        for (size_t i = 0; i < seeds.size(); i++) {
            //For each seed, add its containing snarls and chains to the
            //subtree of the snarl tree
            pos_t pos = seeds[i];
            id_t id = get_id(pos);
            const Snarl* snarl = dist_index.snarlOf(id);
            auto s = node_to_seed.find(id);
            if (s != node_to_seed.end()) {
                s->second.push_back(i);
            } else {
                node_to_seed.emplace(id, vector<size_t>({i}));
                
            }
            pair<id_t, bool> prev (id, false);
            while (snarl != nullptr) { 
                bool seen = false;
                //Add ancestors of snarl to the snarl subtree
                pair<id_t, bool> prev_snarl = move(prev);
                if (snarl_manager.in_nontrivial_chain(snarl)) {
                    //Add chain to snarl subtree
                    const Chain* chain = snarl_manager.chain_of(snarl);
                    prev = make_pair(get_start_of(*chain).node_id(), 
                                     get_start_of(*chain).backward());
                    bool rev_in_chain = snarl_manager.chain_orientation_of(snarl);

                    //A boundary node on a chain is assigned to the snarl 
                    //preceding it in the chain
                    //First node in the chain is assigned to the first snarl 
                    if ((id == snarl->start().node_id() && !rev_in_chain) ||
                        (id == snarl->end().node_id() && rev_in_chain)){ 
                        //If seed is on first boundary node of snarl in chain
                        if (id != get_start_of(*chain).node_id()) {
                            //Unless this is the first node in the chain,
                            //The snarl is switched 
                            snarl = move(snarl_manager.into_which_snarl(
                                       snarl->start().node_id(),
                                      !snarl->start().backward()));
                        }
                    }

                     
                    id_t start_node = rev_in_chain ? snarl->end().node_id() :
                                                     snarl->start().node_id(); 
  
                    DistanceIndex::ChainIndex& chain_index = 
                                dist_index.chainDistances.at(
                                                get_start_of(*chain).node_id());
                    //rank of first node in snarl relative to orientation in chain 
                    size_t rank = chain_index.snarlToIndex[start_node];

                    if (chains_to_snarl_rank.count(chain) == 0) {
                        //If this chain has no seeds yet
                        map<size_t, const Snarl*> chain_snarls;
                        chain_snarls.emplace(rank, snarl);
                        chains_to_snarl_rank.emplace(chain, move(chain_snarls));
                    } else {
                        if (chains_to_snarl_rank[chain].count(rank) == 0) {
                            //Add snarl to chain
                            chains_to_snarl_rank[chain].emplace(rank, snarl); 
                        } 
                    }
                } else {
                    prev = make_pair( snarl->start().node_id(), 
                                      snarl->start().backward());

                }
                if (snarls_to_node.count(snarl) == 0) {
                    hash_set<pair<id_t, bool>> nodes;
                    nodes.insert(prev_snarl);
                    snarls_to_node.emplace(snarl, move(nodes));
                } else {
                    snarls_to_node.at(snarl).insert(prev_snarl);
                    seen = true;
                }
                if (seen) {
                    snarl = nullptr;
                } else {
                    snarl = snarl_manager.parent_of(snarl);
                }

            }
        }

        //Get chains_to_snarl from chains_to_snarl_rank by ordering snarls 
        //by rank
        for (pair<const Chain*, map<size_t, const Snarl*>> ranked_snarls : 
                                                           chains_to_snarl_rank){
            vector<const Snarl*> snarl_list;
            for (pair<size_t, const Snarl*> snarl : ranked_snarls.second){

                //Rank (order) of snarl in the chain 
                snarl_list.push_back( snarl.second );
            } 
            chains_to_snarl.emplace(ranked_snarls.first, move(snarl_list));
        }
#ifdef DEBUG

cerr << "CHAINS: " << endl;
for (auto& c : chains_to_snarl) {
    cerr << get_start_of(*c.first).node_id() << ": ";
    for (const Snarl* s : c.second) {
        cerr << s->start() << " ";
    }
    cerr << endl;
}

cerr << "SNARLS: " << endl; 
for (auto& s : snarls_to_node) {
    cerr << s.first->start() << "  : " ;
    for (auto& n : s.second) {
        cerr << n.first << " " << n.second << ", "; 
    }
    cerr << endl;
}

cerr << "NODES: " << endl;
for (auto& n : node_to_seed){
    cerr << n.first << ": " << endl;
    for (size_t s : n.second ) {
        cerr << seeds[s] << ", ";
    }
    cerr << endl;
} 
#endif
    }

    vector<vector<size_t>> SnarlSeedClusterer::cluster_seeds ( 
                  vector<pos_t> seeds, size_t distance_limit,
                  SnarlManager& snarl_manager, DistanceIndex& dist_index ){
        /* Given a vector of seeds and a limit, find a clustering of seeds where
         * seeds that are closer than the limit cluster together. 
         * Returns a vector of cluster assignments 
         */ 


        /*Store the subtree of the snarl tree that contains seeds*/
        
        //Maps each chain to an ordered vector of child snarls with seeds
        chains_to_snarl_t chains_to_snarl;

        //Maps each snarl to the nodes in its netgraph that contain seeds
        //To prevent the same node from being assigned to multiple snarls,
        //  if a snarl is in a chain, then only the second boundary node of the
        //  snarl (relative to the order of the chain) will be assigned to
        //  that snarl. The first node in a chain is assigned to the first snarl
        snarls_to_node_t snarls_to_node; 

        //Maps nodes to the indexes of the seeds it contains
        node_to_seed_t node_to_seed;

        //Create a union find structure to hold the cluster assignments of
        //each seed. Each seed is initially its own cluster 
        structures::UnionFind union_find_clusters (seeds.size());

        //Maps each cluster group ID to the left and right distances
        hash_map<size_t, pair<int64_t, int64_t>> cluster_dists;

        SnarlSeedClusterer::seed2subtree(seeds, snarl_manager, dist_index, 
                       chains_to_snarl, snarls_to_node, node_to_seed); 
        
        // We track seen-ness by chain, so we don't have to do O(N) work to tag
        // all snarls in chromosome-spanning chains as seen.
        unordered_set<const Chain*> seen_chains;
        for (auto& snarl : snarls_to_node) {
            //Find clusters for each disconnected snarl/chain
            const Snarl* root_snarl = snarl.first;
            if (snarl_manager.parent_of(root_snarl) != nullptr) {
                //If this is not a top level snarl
                continue;
            }
            
            // Look up the (possibly trivial) chain for the snarl
            const Chain* root_chain = snarl_manager.chain_of(root_snarl);

            if (seen_chains.count(root_chain) == 0) {
                // This is the first snarl we are visiting in its chain
                if (snarl_manager.in_nontrivial_chain(root_snarl)){
                    get_clusters_chain( seeds, union_find_clusters,
                                cluster_dists, chains_to_snarl, snarls_to_node, 
                                node_to_seed, distance_limit,
                                snarl_manager, dist_index, root_chain);
                } else {
                    get_clusters_snarl( seeds, union_find_clusters,
                                  cluster_dists, chains_to_snarl, snarls_to_node, 
                                  node_to_seed, distance_limit, 
                              snarl_manager, dist_index, root_snarl, false);
                }
                seen_chains.insert(root_chain);
            }

        }
        return union_find_clusters.all_groups();
        
    };


    tuple<hash_set<size_t>, int64_t, int64_t> SnarlSeedClusterer::get_clusters_node(
                       const vector<pos_t>& seeds,
                       structures::UnionFind& union_find_clusters, 
                       hash_map<size_t, pair<int64_t, int64_t>>& cluster_dists,
                       const node_to_seed_t& node_to_seed,
                       size_t distance_limit, const SnarlManager& snarl_manager,
                       DistanceIndex& dist_index, id_t root,
                       int64_t node_length) {
#ifdef DEBUG 
cerr << "Finding clusters on node " << root << " Which has length " << node_length << endl;
#endif
        /*Find clusters of seeds in this node, root. 
         * rev is true if the left and right distances should be reversed
         * Returns a hash_set of the union find group IDs of the new clusters,
         * the group id of the cluster with seeds furthest to the left, and
         * the group id of the cluster furthest to the right*/

        const vector<size_t>& seed_indices = node_to_seed.at(root);
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
                group_id = union_find_clusters.find_group(seed_i);

                                
            }
            cluster_dists[group_id] = make_pair(best_dist_left, 
                                                best_dist_right);
            group_id = union_find_clusters.find_group(group_id);
            hash_set<size_t> result;
            result.insert(group_id);
#ifdef DEBUG 
assert (group_id == union_find_clusters.find_group(group_id));
cerr << "Found clusters on node " << root << endl;
bool got_left = false;
bool got_right = false;
for (size_t c : result) {
    pair<int64_t, int64_t> dists = cluster_dists[c];
    assert(dists.first >= best_dist_left);
    assert(dists.second >= best_dist_right);
    if (dists.first == best_dist_left) {got_left = true;}
    if (dists.second == best_dist_right) {got_right = true;}
    cerr << "\t" << c << ": left: " << dists.first << " right : " << dists.second << "\t seeds: ";
    vector<size_t> seed_is = union_find_clusters.group(c);
    for (size_t s : seed_is) {
        cerr << seeds[s] << " ";
    }
    cerr << endl;
}
assert(got_left);
assert(got_right);
#endif
            return make_tuple(result, best_dist_left, best_dist_right);
        }

        //indices of union find group ids of clusters in this node
        hash_set<size_t> cluster_group_ids;
        int64_t best_left = -1;
        int64_t best_right = -1;

        for ( size_t i : seed_indices) {
            //For each seed, see if it belongs to a new cluster
            //i is also its own group id
            pos_t seed = seeds[i]; 
            size_t i_group = i;

            int64_t dist_start = get_offset(seed) + 1; 
            int64_t dist_end = node_length - get_offset(seed);
            int64_t dist_left = is_rev(seed) ? dist_end : dist_start;
            int64_t dist_right = is_rev(seed) ? dist_start : dist_end ;
            best_left = DistanceIndex::minPos({dist_left,  best_left});
            best_right = DistanceIndex::minPos({dist_right , best_right});
                        
            bool combined = false;//True if i got incorporated into another cluster

            //group_ids in cluster_group_ids that are no longer the heads of groups
            vector<size_t> to_remove;
            vector<size_t> to_add;
            for (size_t j : cluster_group_ids) {
                //Check which new clusters this seed belongs in 


                pair<int64_t, int64_t> j_dists = cluster_dists[j];
                if (abs(j_dists.first - dist_left) < distance_limit || 
                    abs(j_dists.second - dist_right) < distance_limit ||
                    abs(j_dists.first + dist_right - node_length) < distance_limit || 
                    abs(j_dists.second + dist_left - node_length) < distance_limit ||
                    (j_dists.first < dist_left && (node_length - j_dists.second) > dist_left )  ||
                    (dist_left < j_dists.first && (node_length > dist_right) > j_dists.first)){
                    //If i belongs in this cluster
                    combined = true;
                    union_find_clusters.union_groups(i, j);
                    size_t new_group_id = union_find_clusters.find_group(j);

                    if (new_group_id == i_group) {
                        to_add.push_back(i_group);
                        to_remove.push_back(j);
                        cluster_dists.erase(j);
                    } else {
                        to_remove.push_back(i_group);
                        cluster_dists.erase(i_group);
                    }

                    dist_left = DistanceIndex::minPos({dist_left, j_dists.first});
                    dist_right = DistanceIndex::minPos({dist_right, j_dists.second});

                    cluster_dists[new_group_id] =  make_pair(dist_left, dist_right);
                    i_group = new_group_id;
                }
            }
            for (size_t j : to_add) {
                cluster_group_ids.insert(j);
            }
            for (size_t j : to_remove) {
                //Remove old cluster group ids from cluster_group_ids
                cluster_group_ids.erase(j);
            }
            if (!combined) {
                //If i was not added to any clusters
                cluster_group_ids.insert(i);
                cluster_dists[i] = make_pair(dist_left, dist_right);
            }
        }
#ifdef DEBUG 
cerr << "Found clusters on node " << root << endl;
bool got_left = false;
bool got_right = false;
for (size_t c : cluster_group_ids) {
    pair<int64_t, int64_t> dists = cluster_dists[c];
    assert(dists.first >= best_left);
    assert(dists.second >= best_right);
    if (dists.first == best_left) {got_left = true;}
    if (dists.second == best_right) {got_right = true;}
    cerr << "\t" << c << ": left: " << dists.first << " right : " << dists.second << "\t seeds: ";
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
        return make_tuple(cluster_group_ids, best_left, best_right);
        
    };

    tuple<hash_set<size_t>, int64_t, int64_t>  SnarlSeedClusterer::get_clusters_chain(
                         const vector<pos_t>& seeds,
                         structures::UnionFind& union_find_clusters,
                         hash_map<size_t, pair<int64_t, int64_t>>& cluster_dists,
                         const chains_to_snarl_t& chains_to_snarl,
                         const snarls_to_node_t& snarls_to_node,
                         const node_to_seed_t& node_to_seed,
                         size_t distance_limit, const SnarlManager& snarl_manager,
                         DistanceIndex& dist_index, const Chain* root) {
        #ifdef DEBUG 
        cerr << "Finding clusters on chain " << get_start_of(*root).node_id() << endl;
        #endif
    
        hash_set<size_t> chain_cluster_ids;

        int64_t best_left = -1;
        int64_t best_right = -1;
        DistanceIndex::ChainIndex& chain_index = dist_index.chainDistances.at(
                                                 get_start_of(*root).node_id());
        ChainIterator chain_e = chain_end(*root);
        //The node at which the chain clusters reach
        id_t last_snarl = get_start_of(*root).node_id();
        int64_t last_len = 0;
        id_t start_node;
        id_t end_node;
        const vector<const Snarl*>& snarls_in_chain = chains_to_snarl.at(root);

        for (const Snarl* curr_snarl : snarls_in_chain) {
            /* For each child snarl in the chain, find the clusters of just the
             * snarl, and progressively build up clusters spanning up to that 
             * snarl
             * Snarls are in the order that they are traversed in the chain
             */

            //Clusters of the current snarl/boundry node to be combined with
            //existing chain clusters

            bool rev_in_chain = snarl_manager.chain_orientation_of( curr_snarl);

            DistanceIndex::SnarlIndex& snarl_index = 
                            dist_index.snarlDistances.at(
                               make_pair(curr_snarl->start().node_id(),
                                         curr_snarl->start().backward()));

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
                    cluster_dists[i].second = cluster_dists[i].second == -1 ? -1 :
                                              cluster_dists[i].second + offset;
                }
                best_right = best_right == -1 ? -1 : best_right + offset;
            }

            hash_set<size_t> snarl_clusters; 
            int64_t child_dist_left; int64_t child_dist_right; 
            tie (snarl_clusters, child_dist_left, child_dist_right) = 
                  get_clusters_snarl( seeds, union_find_clusters, cluster_dists,
                               chains_to_snarl, snarls_to_node, node_to_seed, 
                               distance_limit, snarl_manager, dist_index, curr_snarl, 
                               rev_in_chain);
            last_snarl = end_node;
            last_len = end_length;
           

            //Combine the clusters (snarl_clusters) of the current snarl with
            //existing chain clusters. Chain clusters dist_right extends up to
            //the same node as dist_left of snarl_clusters, both include the 
            //length of the shared node.

            //Distance from the start of the chain to the start of the current snarl
            int64_t add_dist_left = chain_index.prefixSum[
                                      chain_index.snarlToIndex[start_node]] - 1;


             
            //Combine snarl clusters that can be reached by looping
            int64_t loop_dist_end = chain_index.loopFd[
                                      chain_index.snarlToIndex[end_node]] - 1 ;
            int64_t loop_dist_start = chain_index.loopRev[
                                     chain_index.snarlToIndex[start_node]] - 1; 

            if (loop_dist_start != -1 || loop_dist_end != -1) {
                vector<size_t> to_remove;
                hash_set<size_t> seen;
                seen.reserve(snarl_clusters.size());
                for (size_t i : snarl_clusters) {
                    for (size_t j : seen ) {
                        //If there is a loop in the chain, then it may
                        //be possible for two clusters on the same snarl to
                        //be combined
                        pair<int64_t, int64_t>& dists_i = cluster_dists[i];
                        pair<int64_t, int64_t>& dists_j = cluster_dists[j];

                        int64_t dist_l = dists_i.first == -1 || dists_j.first == -1
                                          || loop_dist_start == -1 ? -1 
                                   : dists_i.first + dists_j.first 
                                     + loop_dist_start - start_length;
                        int64_t dist_r = dists_i.second == -1 
                                          || dists_j.second == -1 
                                          || loop_dist_end == -1 ? -1 :
                                        dists_i.second + dists_j.second 
                                             + loop_dist_end - end_length;

                        if ((dist_r !=-1 && dist_r < distance_limit)
                                ||(dist_l !=-1 && dist_l < distance_limit)){

                            union_find_clusters.union_groups(i,j);
                            size_t group = union_find_clusters.find_group(i);
                            cluster_dists[group] = make_pair(
                                   DistanceIndex::minPos({dists_i.first, dists_j.first}),
                                   DistanceIndex::minPos({dists_i.second, dists_j.second}));
                            if (group == i) {
                                to_remove.push_back(j);
                            } else {
                                to_remove.push_back(i);
                            }
                        }
                    }
                    seen.insert(i); 
                }
                for (size_t i : to_remove) {
                    cluster_dists.erase(i);
                    snarl_clusters.erase(i);
                }
            }
 
#ifdef DEBUG
cerr << "  Clusters on chain: " << endl;

cerr << "  best left: " << best_left << " best right: " << best_right << endl;
for (size_t c : chain_cluster_ids) {
    pair<int64_t, int64_t> dists = cluster_dists[c];
    cerr << "\tleft: " << dists.first << " right : " << dists.second << "\t seeds: ";
    vector<size_t> ss = union_find_clusters.group(c);
    for (size_t s : ss) {
        cerr << seeds[s] << " ";
    }
    
    cerr << endl;
}
cerr << endl;
#endif



            //Dist from end of start node to farthest seed to the left of snarl 
            //If the right dist of a chain cluster + chain_lim_right is less than
            //the distance limit, then they are in the same cluster
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
                    cluster_dists.erase(j);
                        
                } else {
                    //This snarl becomes a new chain cluster
                    to_add.push_back(j);
                    pair<int64_t, int64_t> d = make_pair(snarl_dists.first == -1 ? -1 :
                                                    snarl_dists.first + add_dist_left, 
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
                            cluster_dists.erase(combined_cluster);
                        } else {
                            to_erase.push_back(i);
                        }
                        combined_cluster = new_group;
                        combined_left = DistanceIndex::minPos({combined_left,
                                                               chain_dists.first});
                        combined_right = DistanceIndex::minPos({combined_right,
                                               chain_dists.second + dist_to_end}); 
                    }
                    cluster_dists.erase(i);
                } else {
                    //If this chain cluster is on its own
                    chain_dists.second += dist_to_end;
                    best_left = DistanceIndex::minPos({best_left, chain_dists.first});
                    best_right = DistanceIndex::minPos({best_right, chain_dists.second});

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
                cluster_dists[combined_cluster] = make_pair(combined_left, combined_right);
                best_left = DistanceIndex::minPos({best_left, combined_left});
                best_right = DistanceIndex::minPos({best_right, combined_right});
            }
                  
#ifdef DEBUG 
cerr << "\t at snarl " << curr_snarl->start().node_id() << ", clusters:" << endl;

for (size_t c : chain_cluster_ids) {
    pair<int64_t, int64_t> dists = cluster_dists[c];
    cerr << "\tleft: " << dists.first << " right : " << dists.second << "\t seeds: ";
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
    assert(dists.first >= best_left);
    assert(dists.second >= best_right);
    if (dists.first == best_left) {got_left = true;}
    if (dists.second == best_right) {got_right = true;}
    cerr << "\t" << c << ": left: " << dists.first << " right : " << dists.second << "\t seeds: ";
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
        return make_tuple(chain_cluster_ids, best_left, best_right) ; 
    };



    tuple<hash_set<size_t>, int64_t, int64_t> SnarlSeedClusterer::get_clusters_snarl(
                         const vector<pos_t>& seeds,
                         structures::UnionFind& union_find_clusters,
                         hash_map<size_t, pair<int64_t, int64_t>>& cluster_dists,
                         const chains_to_snarl_t& chains_to_snarl,
                         const snarls_to_node_t& snarls_to_node,
                         const node_to_seed_t& node_to_seed,
                         size_t distance_limit, const SnarlManager& snarl_manager,
                         DistanceIndex& dist_index, const Snarl* root, 
                         bool rev) {
        #ifdef DEBUG 
        cerr << "Finding clusters on snarl " << root->start() << endl;
        #endif
    
        DistanceIndex::SnarlIndex& snarl_index = dist_index.snarlDistances.at(
                           make_pair(root->start().node_id(),
                                     root->start().backward()));
        int64_t start_length = snarl_index.nodeLength(snarl_index.snarlStart.first);
        int64_t end_length = snarl_index.nodeLength(snarl_index.snarlEnd.first);
//TODO: Find a maybe better way to classify all nodes in the snarl
        vector<id_t> child_nodes;
        vector<const Snarl*> child_snarls; 
        vector<const Chain*> child_chains;
        for (const pair<id_t, bool> child : snarls_to_node.at(root)) {
            //Get and classify all the children of the snarl that have seeds on them
            const Snarl* snarl = snarl_manager.into_which_snarl(child.first, child.second);
            if (snarl == nullptr) {
                snarl = snarl_manager.into_which_snarl(child.first, !child.second);
            }
            if (snarl == nullptr && node_to_seed.find(child.first) != node_to_seed.end()) {
                //If this is a node
                child_nodes.push_back(child.first);
                
            } else if ( child.first == snarl_index.snarlStart.first ||
                        child.first == snarl_index.snarlEnd.first ) {

                if (node_to_seed.find(child.first) != node_to_seed.end()) { 
                    //If this is a boundary node
                    child_nodes.push_back(child.first);
                }
            } else {
                //This is a snarl or a chain
                if (snarl_manager.in_nontrivial_chain(snarl)) {
                    //This is a chain
                    const Chain* chain = snarl_manager.chain_of(snarl); 
                    if (chains_to_snarl.count(chain) != 0) {
                        child_chains.push_back(chain);
                    }
                } else {
                    //This is a snarl
                    if (snarls_to_node.count(snarl) != 0) { 
                        child_snarls.push_back(snarl);
                    }
                }
            }
        }
        size_t num_children = child_nodes.size() + child_snarls.size() + 
                              child_chains.size();

        //clusters of children, indices assume child node, snarl, and chain
        //vectors are contiguous 
        vector<vector<size_t>> child_clusters(num_children);

        //Return value- group ids of all clusters on this snarl
        hash_set<size_t> snarl_cluster_ids;
        int64_t best_left = -1;
        int64_t best_right = -1;
 
        auto combine_clusters = [&] (size_t& new_group, size_t& combined_group, 
                                    pair<int64_t, int64_t> dists){
            //Helper function to combine clusters in two nodes of the same snarl 
            if (combined_group == -1) {
                snarl_cluster_ids.insert(new_group);
                cluster_dists[new_group] = dists;
                combined_group = new_group;
            } else {

                combined_group = union_find_clusters.find_group(combined_group);
                pair<int64_t, int64_t> old_dists = cluster_dists[combined_group];
                union_find_clusters.union_groups(new_group, combined_group);
                size_t new_g = union_find_clusters.find_group(new_group);
                if (new_g != new_group) {
                    snarl_cluster_ids.erase(new_group);
                    cluster_dists.erase(new_group); 
                } 
                if (new_g != combined_group) {
                    snarl_cluster_ids.erase(combined_group);
                    cluster_dists.erase(combined_group);
                }
                snarl_cluster_ids.insert(new_g);
                cluster_dists[new_g] = make_pair(
                              DistanceIndex::minPos({dists.first, old_dists.first}),
                              DistanceIndex::minPos({dists.second, old_dists.second}));
                new_group = new_g;
                combined_group = new_g;
            }
            return;
        };
        //Maps each cluster of child nodes to its left and right distances
        hash_map<size_t, pair<int64_t, int64_t>> old_dists;
        //Maps each child node to the left and right bounds of its clusters
        hash_map<size_t, pair<int64_t, int64_t>> dist_bounds;
        dist_bounds.resize(num_children);

        for (size_t i = 0; i < num_children ; i++) {
            //Go through each child node of the netgraph and find clusters

            id_t curr_node;
            bool curr_rev;
            //The clusters furthest to the left and right for this child node
            int64_t child_dist_left; int64_t child_dist_right;
            if (i < child_nodes.size()) {
                curr_node = child_nodes[i];
                curr_rev = false;
                int64_t node_len = snarl_index.nodeLength(curr_node);
                //If this node has seeds, cluster them
                hash_set<size_t> c; 
                tie (c, child_dist_left, child_dist_right) = get_clusters_node(
                         seeds, union_find_clusters, cluster_dists, 
                         node_to_seed, distance_limit, snarl_manager, 
                         dist_index, curr_node, node_len);
                old_dists.resize(old_dists.size() + c.size());
                child_clusters[i].insert(child_clusters[i].end(),
                           make_move_iterator(c.begin()),
                           make_move_iterator(c.end())); 

            } else if (i < child_snarls.size() + child_nodes.size()) {

                const Snarl* curr_snarl = child_snarls[i - child_nodes.size()];
                curr_node = curr_snarl->start().node_id();
                curr_rev = curr_snarl->start().backward();
                hash_set<size_t> c; 
                tie (c, child_dist_left, child_dist_right) = 
                      get_clusters_snarl(
                         seeds, union_find_clusters, cluster_dists, chains_to_snarl,
                         snarls_to_node, node_to_seed,distance_limit, snarl_manager,
                         dist_index, curr_snarl, false);
                old_dists.resize(old_dists.size() + c.size());
                child_clusters[i].insert(child_clusters[i].end(),
                                            make_move_iterator(c.begin()),
                                            make_move_iterator(c.end())); 
            } else {
                const Chain* curr_chain =  child_chains[i - child_nodes.size() -
                                             child_snarls.size()];
                Visit start = get_start_of(*curr_chain);
                curr_node = start.node_id();
                curr_rev = start.backward();
                hash_set<size_t> c; 
                tie (c, child_dist_left, child_dist_right) =
                       get_clusters_chain(seeds, 
                             union_find_clusters, cluster_dists, chains_to_snarl,
                             snarls_to_node, node_to_seed, distance_limit,
                             snarl_manager, dist_index, curr_chain);

                old_dists.resize(old_dists.size() + c.size());
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
            bool first_pass = true;

            for (size_t j = 0 ; j <= i ; j++){
                //Go through child nodes up to and including i
                id_t other_node;
                bool other_rev;
                if (j < child_nodes.size()) {
                    other_node = child_nodes[j];
                    other_rev = false;
                } else if (j < child_snarls.size() + child_nodes.size()) {
                    const Snarl* other_snarl =  
                               child_snarls[j - child_nodes.size()];
                    other_node = other_snarl->start().node_id();
                    other_rev = other_snarl->start().backward();
                } else {
                    const Chain* other_chain = child_chains [
                      j - child_nodes.size() - child_snarls.size()];
                    Visit start = get_start_of(*other_chain);
                    other_node = start.node_id();
                    other_rev = start.backward();
                }
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

                //group ids of clusters combined between node i left and node j left, etc
                size_t group_l_l = -1;
                size_t group_l_r = -1;
                size_t group_r_l = -1;
                size_t group_r_r = -1;

                pair<int64_t, int64_t> dist_bounds_j = dist_bounds[j];

                for (size_t c_i = 0 ; c_i < children_i.size() ; c_i ++) {
                    //TODO: Make children clusters a set so only go through each combined cluster once
                    //for each cluster of child node i
                    size_t c = children_i[c_i];
                    size_t c_group = union_find_clusters.find_group(c);

                    //True if this cluster got combined with an existing cluster
                    bool combined = false;


                    pair<int64_t, int64_t> new_dists;
                    pair<int64_t, int64_t> dists_c;

                    if (first_pass) {
                        dists_c = cluster_dists[c];
                        old_dists[c] = dists_c;

    
                        //find distances to ends of snarl for this cluster 
                        int64_t new_dist_s_l = dist_s_l == -1  ||  dists_c.first == -1 
                                    ? -1 : dist_s_l + dists_c.first;
                        int64_t new_dist_s_r = dist_s_r == -1 ||  dists_c.second == -1 
                                    ? -1 : dist_s_r + dists_c.second;
                        int64_t new_dist_e_l = dist_e_l == -1 || dists_c.first == -1
                                    ? -1 : dist_e_l + dists_c.first;
                        int64_t new_dist_e_r = dist_e_r == -1 || dists_c.second == -1
                                    ? -1 : dist_e_r + dists_c.second;
      
                        new_dists = make_pair (
                                 DistanceIndex::minPos({new_dist_s_l, new_dist_s_r}), 
                                 DistanceIndex::minPos({new_dist_e_l,new_dist_e_r}));
                        if (curr_node == snarl_index.snarlStart.first){
                            new_dists.first = snarl_index.snarlStart.second ? 
                                    dists_c.second : dists_c.first ;
                        }
                        if (curr_node == snarl_index.snarlEnd.first) {
                            new_dists.second = snarl_index.snarlEnd.second ? 
                                     dists_c.first : dists_c.second ;
                        } 
                        best_left = DistanceIndex::minPos({best_left, new_dists.first});
                        best_right = DistanceIndex::minPos({best_right, new_dists.second});


                    } else {
                        dists_c = old_dists[c];
                        new_dists = cluster_dists[c_group];
                        combined = true;
                    }
                    

                    if (dist_l_l != -1 && dists_c.first != -1 && dist_bounds_j.first != -1  
                             && dist_l_l + dists_c.first +  dist_bounds_j.first < distance_limit) {
                        //If cluster c can be combined with clusters in j from the left
                        //of both of them
                        combined = true;
                        combine_clusters(c_group, group_l_l, new_dists);
                    }
                    if (dist_l_r != -1 && dists_c.first != -1 && dist_bounds_j.second != -1  
                                && dist_l_r + dists_c.first + dist_bounds_j.second < distance_limit) {
                        combined = true;
                        combine_clusters(c_group, group_l_r, new_dists);

                    }
                    if (dist_r_l != -1 && dists_c.second != -1 && dist_bounds_j.first != -1 
                                && dist_r_l + dists_c.second + dist_bounds_j.first < distance_limit) {
                        combined = true;
                        combine_clusters(c_group, group_r_l, new_dists);
                    }
                    if (dist_r_r != -1 && dists_c.second != -1 && dist_bounds_j.second != -1 
                               && dist_r_r + dists_c.second + dist_bounds_j.second < distance_limit) {
                        combined = true;
                        combine_clusters(c_group, group_r_r, new_dists);
                    }
                    if (!combined) {
                        snarl_cluster_ids.insert(c_group);
                        cluster_dists[c_group] = new_dists;
                    }

                }
                first_pass = false;
                if (max({dist_l_l, dist_l_r, dist_r_l, dist_r_r}) != -1) {
                    //If the two nodes are reachable
                    vector<size_t>& children_j =  child_clusters[j];

                    for (size_t k_i = 0 ; k_i < children_j.size() ; k_i++){
                        size_t k = children_j[k_i];
                        //For each other cluster, see which this child cluster belongs to
                        //k will already be part of a cluster in snarl_cluster_ids
                        //but since we need to know the node that the snarl is on we
                        //can't just loop through snarl_cluster_ids
                        pair<int64_t, int64_t>& dist_bounds_k = old_dists[k];
                        size_t k_group = union_find_clusters.find_group(k);
                        pair<int64_t, int64_t> dists_k = cluster_dists[k_group];
                    

                        if (dist_l_l != -1 && child_dist_left != -1 && dist_bounds_k.first != -1 
                              && dist_l_l + child_dist_left +  dist_bounds_k.first < distance_limit){

                            combine_clusters(k_group, group_l_l, dists_k);
                        }
                        if (dist_l_r != -1 && child_dist_left != -1 && dist_bounds_k.second != -1  
                            && dist_l_r + child_dist_left + dist_bounds_k.second < distance_limit ) {

                            combine_clusters(k_group, group_l_r, dists_k);
                        }
                        if (dist_r_l != -1 && child_dist_right != -1 && dist_bounds_k.first != -1  
                              && dist_r_l + child_dist_right + dist_bounds_k.first < distance_limit) {

                            combine_clusters(k_group, group_r_l, dists_k);
                        }
                        if (dist_r_r != -1 && child_dist_right != -1 && dist_bounds_k.second != -1
                             && dist_r_r + child_dist_right + dist_bounds_k.second < distance_limit) {

                            combine_clusters(k_group, group_r_r, dists_k);
                        }
                    }
                }

#ifdef DEBUG 
cerr << "At i node " << curr_node << " j node " << other_node << endl;
cerr << "    with best left and right values: " << best_left << " " << best_right << endl;
bool got_left = false;
bool got_right = false;
for (size_t c : snarl_cluster_ids) {
    pair<int64_t, int64_t> dists = cluster_dists[c];
    if (dists.first == best_left) {got_left = true;}
    if (dists.second == best_right) {got_right = true;}
    cerr << "\t" << c << ": left: " << dists.first << " right : " << dists.second << "\t seeds: ";
    vector<size_t> seed_is = union_find_clusters.group(c);
    assert(dists.first >= best_left);
    assert(dists.second >= best_right);
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
cerr << "    with best left and right values: " << best_left << " " << best_right << endl;
bool got_left = false;
bool got_right = false;
for (size_t c : snarl_cluster_ids) {
    pair<int64_t, int64_t> dists = cluster_dists[c];
    if (dists.first == best_left) {got_left = true;}
    if (dists.second == best_right) {got_right = true;}
    cerr << "\t" << c << ": left: " << dists.first << " right : " << dists.second << "\t seeds: ";
    vector<size_t> seed_is = union_find_clusters.group(c);
    assert(dists.first >= best_left);
    assert(dists.second >= best_right);
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
        return make_tuple(snarl_cluster_ids, best_left, best_right);
    };
}
