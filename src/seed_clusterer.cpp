#include "seed_clusterer.hpp"

//#define DEBUG 

namespace vg {

    SnarlSeedClusterer::SnarlSeedClusterer() {
    };


    void SnarlSeedClusterer::seed2subtree( 
               const vector<pos_t>& seeds,  const SnarlManager& snarl_manager, 
               DistanceIndex& dist_index, chains_to_snarl_t& chains_to_snarl,
               snarls_to_node_t& snarls_to_node, node_to_seed_t& node_to_seed,
               vector<SnarlSeedClusterer::seed_index_list>& seed_list ) {

       
        /* Given a snarl tree and seeds, find the subtree of the snarl tree
         * that contains seeds */ 
        //Maps chains to the snarls containing seeds
        hash_map< const Chain*, map<size_t,const Snarl*>> chains_to_snarl_rank;
        for (size_t i = 0; i < seeds.size(); i++) {
            //For each seed, add its containing snarls and chains to the
            //subtree of the snarl tree
            pos_t pos = seeds[i];
            id_t id = get_id(pos);
            seed_list[i].index = i;
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
;
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
for (pair<const Chain*, vector<const Snarl*>>& c : chains_to_snarl) {
    cerr << get_start_of(*c.first).node_id() << ": ";
    for (const Snarl* s : c.second) {
        cerr << s->start() << " ";
    }
    cerr << endl;
}

cerr << "SNARLS: " << endl; 
for (pair<const Snarl*, hash_set<pair<id_t, bool>>>& s : snarls_to_node) {
    cerr << s.first->start() << "  : " ;
    for (pair<id_t, bool>& n : s.second) {
        cerr << n.first << " " << n.second << ", "; 
    }
    cerr << endl;
}

cerr << "NODES: " << endl;
for (pair<id_t, vector<size_t>>& n : node_to_seed){
    cerr << n.first << ": " << endl;
    for (size_t s : n.second ) {
        cerr << seeds[s] << ", ";
    }
    cerr << endl;
} 
#endif
    }

    vector<hash_set<size_t>> SnarlSeedClusterer::cluster_seeds ( 
                  vector<pos_t> seeds, size_t distance_limit,
                  SnarlManager& snarl_manager, DistanceIndex& dist_index ){
        /* Given a vector of seeds and a limit, find a clustering of seeds where
         * seeds that are closer than the limit cluster together. 
         * Returns a vector of cluster assignments 
         */ 


        /*Store the subtree of the snarl tree that contains seeds*/
        

        //This will be filled in based on chains_to_snarl_rank
        //snarls will be ordered by order in chain
        chains_to_snarl_t chains_to_snarl;

        //Maps each snarl to the nodes in its netgraph that contain seeds
        //To prevent the same node from being assigned to multiple snarls,
        //  if a snarl is in a chain, then only the second boundary node of the
        //  snarl (relative to the order of the chain) will be assigned to
        //  that snarl. The first node in a chain is assigned to the first snarl
        snarls_to_node_t snarls_to_node; 

        //Maps nodes to the indexes of the seeds it contains
        node_to_seed_t node_to_seed;

        vector<SnarlSeedClusterer::seed_index_list> seed_list(seeds.size());

        SnarlSeedClusterer::seed2subtree(  seeds, snarl_manager, dist_index, 
                       chains_to_snarl, snarls_to_node, node_to_seed,seed_list); 

  
        // This will hold all top-level snarls that have seeds under them
        //vector<const Snarl*> top_snarls;
        //for (auto& snarl : snarls_to_node) {
        //    if (snarl_manager.parent_of(snarl.first) == nullptr) {
        //        // This is top level
        //        top_snarls.push_back(snarl.first);
        //    }
        //}
        
        // We track seen-ness by chain, so we don't have to do O(N) work to tag
        // all snarls in chromosome-spanning chains as seen.
        unordered_set<const Chain*> seen_chains;
        vector<hash_set<size_t>> cluster_assignments; 
        //for (const Snarl* root_snarl : top_snarls) {
        for (auto& snarl : snarls_to_node) {
	    const Snarl* root_snarl = snarl.first;
            if (snarl_manager.parent_of(root_snarl) != nullptr) {
		continue;
	    }
            //Find clusters for each disconnected snarl/chain
            
            // Look up the (possibly trivial) chain for the snarl
            const Chain* root_chain = snarl_manager.chain_of(root_snarl);

            vector<SnarlSeedClusterer::cluster_t> clusters;
            if (seen_chains.count(root_chain) == 0) {
                // This is the first snarl we are visiting in its chain
                if (snarl_manager.in_nontrivial_chain(root_snarl)){
                    clusters =  get_clusters_chain( seeds, seed_list,
                                 chains_to_snarl, snarls_to_node, 
                                node_to_seed, distance_limit,
                                snarl_manager, dist_index, root_chain);
                } else {
                    clusters =  get_clusters_snarl( seeds, seed_list,
                                  chains_to_snarl, snarls_to_node, 
                                  node_to_seed, distance_limit, 
                              snarl_manager, dist_index, root_snarl, false);
                }
                seen_chains.insert(root_chain);
            }

            for (size_t i = 0 ; i < clusters.size() ; i++) {
                SnarlSeedClusterer::cluster_t& cluster = clusters[i]; 
                cluster_assignments.emplace_back(); 
                SnarlSeedClusterer::seed_index_list* j = cluster.seeds_start;
                while (j != nullptr) {
                    cluster_assignments.back().insert(j->index);
                    j = j->next;
                }
            }
        }
        return cluster_assignments;
        
    };


    vector<SnarlSeedClusterer::cluster_t> SnarlSeedClusterer::get_clusters_node(
                        const vector<pos_t>& seeds,
                        vector<SnarlSeedClusterer::seed_index_list>& seed_list, 
                        const node_to_seed_t& node_to_seed,
                        size_t distance_limit, const SnarlManager& snarl_manager,
                        DistanceIndex& dist_index, id_t root,
                        int64_t node_length) {
#ifdef DEBUG 
cerr << "Finding clusters on node " << root << " Which has length " << node_length << endl;
#endif
        /*Find clusters of seeds in this node, root. rev is true if the left and right
         * distances should be reversed */ 

        const vector<size_t>& seed_indices = node_to_seed.at(root);
        if (distance_limit > node_length) {
            //If the limit is greater than the node length, then all the seeds
            //on this node must be in the same cluster
            vector<SnarlSeedClusterer::cluster_t> clusters(1);
            for (size_t seed_i : seed_indices) {
                //Distance from seed to start of node 
                
                pos_t seed = seeds[seed_i]; 
                int64_t dist_start = get_offset(seed) + 1; 
                int64_t dist_end = node_length - get_offset(seed);
                int64_t dist_left = is_rev(seed) ? dist_end : dist_start;
                int64_t dist_right = is_rev(seed) ? dist_start : dist_end ;

                clusters[0].dist_left = DistanceIndex::minPos({
                                        dist_left, clusters[0].dist_left});
                clusters[0].dist_right = DistanceIndex::minPos({
                                        dist_right, clusters[0].dist_right});

                SnarlSeedClusterer::seed_index_list* curr_seed_i = &seed_list[seed_i];
                if (clusters[0].seeds_start == nullptr) {
                    clusters[0].seeds_start = curr_seed_i; 
                    clusters[0].seeds_end = curr_seed_i; 
                } else {
                    clusters[0].seeds_end->next = curr_seed_i;
                    clusters[0].seeds_end = curr_seed_i;
                }
                
            }
            return clusters;
        }

        //vector of clusters where each cluster contains indices of seeds
        vector<hash_set<size_t>> cluster_indices;

        for ( size_t i : seed_indices) {
            //For each seed, see if it belongs to a new cluster
            pos_t seed = seeds[i]; 

            //Which clusters in cluster_indices this seed belongs to
            hash_set<size_t> curr_cluster_assignments;

            //Distance from seed to start of node 
            int64_t offset = is_rev(seed) ? node_length - get_offset(seed) : 
                                            get_offset(seed)+1;
                        
            for (size_t j = 0 ; j < cluster_indices.size(); j++) {
                //Check which new clusters this seed belongs in 
                hash_set<size_t>& curr_cluster = cluster_indices[j];
                for (size_t k : curr_cluster) {  
                    //Check current seed against seeds in curr cluster
                    pos_t other_seed = seeds[k];
                    int64_t other_offset = is_rev(other_seed) ?  
                                       node_length - get_offset(other_seed) : 
                                       get_offset(other_seed) + 1;
                    if (abs(offset - other_offset) < distance_limit) {
                        curr_cluster_assignments.insert(j);
                        break;
                    } 
                }
            }
            //Combine all clusters that the seed belongs to and make new cluster
            if (curr_cluster_assignments.size() == 0) {

                cluster_indices.emplace_back(); 
                cluster_indices.back().insert(i);

            } else {

                vector<hash_set<size_t>> combined_cluster;
                combined_cluster.emplace_back();
                combined_cluster.front().insert(i);
                for (size_t j = 0 ; j < cluster_indices.size(); j++) {

                    if (curr_cluster_assignments.count(j) != 0) {
                        //If this seed belongs to this cluster

                        hash_set<size_t>& old_cluster = cluster_indices[j];
                        combined_cluster.front().insert( old_cluster.begin(), 
                                            old_cluster.end());
                    } else {

                        combined_cluster.push_back(cluster_indices[j]);
                    }
                }
                cluster_indices = move(combined_cluster);
            }
        }

        //Create actual clusters from indices
        vector<SnarlSeedClusterer::cluster_t> clusters(cluster_indices.size());
        for (int j = 0; j < cluster_indices.size(); ++j) {
	    hash_set<size_t>& indices = cluster_indices[j];

            for (size_t i : indices) {
                SnarlSeedClusterer::seed_index_list* curr_seed_i = &seed_list[i];
                if (clusters[j].seeds_start == nullptr) {
                    clusters[j].seeds_start = curr_seed_i; 
                    clusters[j].seeds_end = curr_seed_i; 
                } else {
                    clusters[j].seeds_end->next = curr_seed_i;
                    clusters[j].seeds_end = curr_seed_i;
                }

                pos_t seed = seeds[i];
                int64_t dist_start = get_offset(seed) + 1; 
                int64_t dist_end = node_length - get_offset(seed);
                int64_t dist_left = is_rev(seed) ? dist_end : dist_start;
                int64_t dist_right = is_rev(seed) ? dist_start : dist_end ;
                clusters[j].dist_left = DistanceIndex::minPos({
                                          dist_left, clusters[j].dist_left});
                clusters[j].dist_right = DistanceIndex::minPos(
                                         {dist_right, clusters[j].dist_right});
            }
        }
#ifdef DEBUG 
cerr << "Found clusters on node " << root << endl;
for (cluster_t& c : clusters) {
    cerr << "\tleft: " << c.dist_left << " right : " << c.dist_right << "\t seeds: ";
    SnarlSeedClusterer::seed_index_list* s = c.seeds_start;
    while (s != nullptr){
        cerr << seeds[s->index] << "\t";
        s = s->next;
    }
    cerr << endl;
}
#endif
        return clusters;
        
    };

    vector<SnarlSeedClusterer::cluster_t> 
                  SnarlSeedClusterer::get_clusters_chain(
                         const vector<pos_t>& seeds,
                         vector<SnarlSeedClusterer::seed_index_list>& seed_list,
                         const chains_to_snarl_t& chains_to_snarl,
                         const snarls_to_node_t& snarls_to_node,
                         const node_to_seed_t& node_to_seed,
                         size_t distance_limit, const SnarlManager& snarl_manager,
                         DistanceIndex& dist_index, const Chain* root) {
        #ifdef DEBUG 
        cerr << "Finding clusters on chain " << get_start_of(*root).node_id() << endl;
        #endif
    
        vector<SnarlSeedClusterer::cluster_t> chain_clusters;

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
                vector<cluster_t> new_chain_clusters(chain_clusters.size());

                for (int i = 0; i < chain_clusters.size(); ++i) {
		    cluster_t& curr = chain_clusters[i];
                    //Update the clusters to reach the current snarl
                    new_chain_clusters[i].dist_right = 
                                                      curr.dist_right + offset;
                    new_chain_clusters[i].dist_left = curr.dist_left;

                    if (new_chain_clusters[i].seeds_start == nullptr) {

                        new_chain_clusters[i].seeds_start =curr.seeds_start;
                        new_chain_clusters[i].seeds_end = curr.seeds_end; 
                    } else {

                        new_chain_clusters[i].seeds_end->next = 
                                                               curr.seeds_start;
                        new_chain_clusters[i].seeds_end = curr.seeds_end;
                    }
                }
                chain_clusters = move(new_chain_clusters);
            }

            vector<SnarlSeedClusterer::cluster_t> curr_clusters = 
                  get_clusters_snarl( seeds, seed_list, chains_to_snarl,
                               snarls_to_node, node_to_seed, distance_limit,
                               snarl_manager, dist_index, curr_snarl, 
                               rev_in_chain);
            last_snarl = end_node;
            last_len = end_length;
            
#ifdef DEBUG
cerr << "  Cluster to add to chain: " << endl;

for (SnarlSeedClusterer::cluster_t& c : curr_clusters) {
    cerr << "        cluster left " << c.dist_left << " right: " << c.dist_right << "seeds : " ;
    SnarlSeedClusterer::seed_index_list *s = c.seeds_start;
    while (s != nullptr){
        cerr << seeds[s->index] << "\t";
        s = s->next;
    }
    
    cerr << endl;
}
cerr << endl;
#endif

            //Combine the clusters (curr_clusters) of the current snarl with
            //existing chain clusters. Chain clusters dist_right extends up to
            //the same node as dist_left of curr_clusters, both include the 
            //length of the shared node.
         
            vector<pair<hash_set<size_t>, hash_set<size_t>>> 
                                                    new_chain_cluster_indices(chain_clusters.size());
            for (size_t i = 0 ; i < chain_clusters.size() ; i++) {
                //Pair of old cluster indices, new cluster indices
                //Old clusters are all initially separate
                new_chain_cluster_indices[i].first.insert(i);
            }  

            for (size_t j = 0; j < curr_clusters.size() ; j++) {
                /* For each of the clusters for the current snarl,
                 * find which chain clusters it belongs to and
                 * combine shared chain clusters */ 
                SnarlSeedClusterer::cluster_t& snarl_cluster = curr_clusters[j];
                hash_set<size_t> curr_cluster_assignments;
                //Which chain cluster the snarl cluster is near

                for (size_t i = 0; i < new_chain_cluster_indices.size(); i++ ) {
                    //For each old chain cluster
                  
                    int64_t loop_dist_end = chain_index.loopFd[
                                       chain_index.snarlToIndex[end_node]] - 1 ;
                    int64_t loop_dist_start = chain_index.loopRev[
                                     chain_index.snarlToIndex[start_node]] - 1; 
                    pair<hash_set<size_t>, hash_set<size_t>>& 
                              old_chain_clusters = new_chain_cluster_indices[i];

                    for (size_t c : old_chain_clusters.first) {
                        //Check snarl cluster against existing chain clusters
                        SnarlSeedClusterer::cluster_t& chain_cluster = 
                                                          chain_clusters[c];
                        if (chain_cluster.dist_right + 
                                    snarl_cluster.dist_left - start_length <
                                   distance_limit ) {
                           //If two clusters are close enough to combine
                            curr_cluster_assignments.insert(i);
                            break;
			}
                    }
                    if (loop_dist_start != -1 || loop_dist_end != -1) {
                        //If there is a loop in the chain, then it may
                        //be possible for two clusters on the same snarl to
                        //be combined
                        for (size_t c : old_chain_clusters.second) {
                             //Check snarl cluster against other snarl clusters
                            SnarlSeedClusterer::cluster_t& other_cluster=
                                                curr_clusters[c];       
                            int64_t dist_l = snarl_cluster.dist_left == -1 
                                          || other_cluster.dist_left == -1 
                                          || loop_dist_start == -1 ? -1 
                                   : snarl_cluster.dist_left 
                                     + other_cluster.dist_left
                                     + loop_dist_start - start_length;
                            int64_t dist_r = snarl_cluster.dist_right == -1 
                                          || other_cluster.dist_right == -1 
                                          || loop_dist_end == -1 ? -1 :
                                        snarl_cluster.dist_right + 
                                   other_cluster.dist_right + loop_dist_end
                                     - end_length;
                            if ((dist_r !=-1 && dist_r < distance_limit)
                                ||(dist_l !=-1 && dist_l < distance_limit)){

                                 curr_cluster_assignments.insert(i);
                                 break;
			    }
                        }
                    }
                }

                //Combine chain clusters that are in the same cluster as this snarl cluster
                vector<pair<hash_set<size_t>, hash_set<size_t>>> new_chain_indices;
                pair<hash_set<size_t>, hash_set<size_t>> new_chain_cluster;
                new_chain_cluster.second.insert(j);
                for (size_t i = 0; i < new_chain_cluster_indices.size(); i++ ) {
                    /* Combine chain clusters that share this snarl cluster*/
                    if (curr_cluster_assignments.count(i) != 0) {
                        pair<hash_set<size_t>, hash_set<size_t>>& old_cluster
                                             = new_chain_cluster_indices[i];
                        new_chain_cluster.first.insert( 
                                                      make_move_iterator(old_cluster.first.begin()),
                                                      make_move_iterator(old_cluster.first.end()));
                        new_chain_cluster.second.insert(
                                                  make_move_iterator(old_cluster.second.begin()),
                                                  make_move_iterator(old_cluster.second.end()));
                    } else {
                        new_chain_indices.push_back( 
                                                 move(new_chain_cluster_indices[i]) );
                    }
                }

                new_chain_indices.push_back(new_chain_cluster);
                new_chain_cluster_indices = move(new_chain_indices);
            }
 
            /* Each element in new_chain_cluster_indices contains 
             * the indices of each chain cluster and snarl cluster 
             * that belong together. Combine these into actual clusters
             */
 
            vector<SnarlSeedClusterer::cluster_t> new_chain_clusters(new_chain_cluster_indices.size());
            for (size_t k = 0; k < new_chain_cluster_indices.size(); ++k) {
                //Create a new cluster (back of new_chain_clusters)
                //from this combined cluster
                pair<hash_set<size_t>, hash_set<size_t>>& 
                          new_chain_cluster = new_chain_cluster_indices[k];

                for (size_t i : new_chain_cluster.first) {
                    //Combine the chain clusters
                    SnarlSeedClusterer::cluster_t& old_chain_cluster = 
                                                              chain_clusters[i];

                    new_chain_clusters[k].dist_left = DistanceIndex::minPos(
                                  {old_chain_cluster.dist_left, 
                                   new_chain_clusters[k].dist_left});
                    int64_t old_dist_right =  old_chain_cluster.dist_right == -1 ? 
                                -1 : old_chain_cluster.dist_right + dist_to_end;

                    new_chain_clusters[k].dist_right = 
                               DistanceIndex::minPos({ old_dist_right, 
                                         new_chain_clusters[k].dist_right});

                    if (new_chain_clusters[k].seeds_start == nullptr) {
                        new_chain_clusters[k].seeds_start = 
                                                old_chain_cluster.seeds_start; 
                        new_chain_clusters[k].seeds_end = 
                                                old_chain_cluster.seeds_end; 
                    } else {
                        new_chain_clusters[k].seeds_end->next =
                                                  old_chain_cluster.seeds_start;
                        new_chain_clusters[k].seeds_end = 
                                                  old_chain_cluster.seeds_end;
                    }
                }


                for (size_t j : new_chain_cluster.second) {
                    //combine snarl clusters
                    SnarlSeedClusterer::cluster_t& old_snarl_cluster = 
                                                        curr_clusters[j];
                    new_chain_clusters[k].dist_right =DistanceIndex::minPos(
                                   {old_snarl_cluster.dist_right, 
                                      new_chain_clusters[k].dist_right});

                    int64_t old_dist_left = old_snarl_cluster.dist_left == -1 ? 
                         -1 :  old_snarl_cluster.dist_left 
                               + chain_index.prefixSum[
                                      chain_index.snarlToIndex[start_node]] - 1;
                    new_chain_clusters[k].dist_left = DistanceIndex::minPos(
                                      {old_dist_left, 
                                       new_chain_clusters[k].dist_left});

                    if (new_chain_clusters[k].seeds_start == nullptr) {
                        new_chain_clusters[k].seeds_start =
                                                  old_snarl_cluster.seeds_start;
                        new_chain_clusters[k].seeds_end =
                                                   old_snarl_cluster.seeds_end; 
                    } else {
                        new_chain_clusters[k].seeds_end->next = 
                                                  old_snarl_cluster.seeds_start;
                        new_chain_clusters[k].seeds_end = 
                                                    old_snarl_cluster.seeds_end;
                    }
                }
            }
            chain_clusters = move(new_chain_clusters);
#ifdef DEBUG 
cerr << "\t at snarl " << curr_snarl->start().node_id() << ", clusters:" << endl;

for (SnarlSeedClusterer::cluster_t& curr_cluster :chain_clusters) {
    cerr << "\t \t left: " << curr_cluster.dist_left << " right: " << curr_cluster.dist_right << " seeds; ";
    SnarlSeedClusterer::seed_index_list *s = curr_cluster.seeds_start;
    while (s != nullptr){
        cerr << seeds[s->index] << "\t";
        s = s->next;
    }
    
    cerr << endl;
}
#endif
                

            
        }
         
#ifdef DEBUG 
cerr << "Found clusters on chain " << get_start_of(*root).node_id() << endl;
for (cluster_t& c : chain_clusters) {
    cerr << "\tleft: " << c.dist_left << " right : " << c.dist_right << "\tseeds: ";
    SnarlSeedClusterer::seed_index_list *s = c.seeds_start;
    while (s != nullptr){
        cerr << seeds[s->index] << "\t";
        s = s->next;
    }
    cerr << endl;
}
#endif

        if (last_snarl != get_end_of(*root).node_id()) {
           //Extend the right bound of each cluster to the end of the chain
           int64_t dist_to_end = chain_index.chainLength()
                       - chain_index.prefixSum[
                                      chain_index.snarlToIndex[last_snarl]] +
                                - last_len;
           for (size_t i = 0; i < chain_clusters.size() ; i++) {
               chain_clusters[i].dist_right += dist_to_end;
           }
        }

        return chain_clusters; 
    };





    vector<SnarlSeedClusterer::cluster_t> SnarlSeedClusterer::get_clusters_snarl(
                         const vector<pos_t>& seeds,
                         vector<SnarlSeedClusterer::seed_index_list>& seed_list,
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
        vector<vector<SnarlSeedClusterer::cluster_t>> child_clusters(num_children);
        //Vector of clusters where a cluster is represented by indices into
        //child_nodes/snarls/chains
        vector<pair<vector<hash_set<size_t>>, pair<int64_t, int64_t>>>
                                                     combined_clusters_indices;
         
        for (size_t i = 0; i < num_children ; i++) {
            //Go through each child node of the netgraph and find clusters

            id_t curr_node;
            bool curr_rev;
            if (i < child_nodes.size()) {
                curr_node = child_nodes[i];
                curr_rev = false;
                int64_t node_len = snarl_index.nodeLength(curr_node);
                    //If this node has seeds, cluster them
		    vector<SnarlSeedClusterer::cluster_t> c = get_clusters_node(
                             seeds, seed_list, node_to_seed,
                             distance_limit, snarl_manager, dist_index, 
                             curr_node, node_len);
                    child_clusters[i].insert(child_clusters[i].end(),
					     make_move_iterator(c.begin()),
					     make_move_iterator(c.end())); 

            } else if (i < child_snarls.size() + child_nodes.size()) {

                const Snarl* curr_snarl = child_snarls[i - child_nodes.size()];
                curr_node = curr_snarl->start().node_id();
                curr_rev = curr_snarl->start().backward();
		    vector<SnarlSeedClusterer::cluster_t> c = get_clusters_snarl(
                         seeds, seed_list, chains_to_snarl, snarls_to_node, 
                         node_to_seed,distance_limit, snarl_manager, dist_index,
                         curr_snarl, false);
                    child_clusters[i].insert(child_clusters[i].end(),
					     make_move_iterator(c.begin()),
					     make_move_iterator(c.end())); 
            } else {
                const Chain* curr_chain =  child_chains[i - child_nodes.size() -
                                             child_snarls.size()];
                Visit start = get_start_of(*curr_chain);
                curr_node = start.node_id();
                curr_rev = start.backward();
		    vector<SnarlSeedClusterer::cluster_t> c = get_clusters_chain(seeds, 
                             seed_list, chains_to_snarl,
                             snarls_to_node, node_to_seed, distance_limit,
                             snarl_manager, dist_index, curr_chain);
                    child_clusters[i].insert(child_clusters[i].end(),
					     make_move_iterator(c.begin()),
					     make_move_iterator(c.end())); 
            }
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



            for (size_t c = 0 ; c < child_clusters[i].size() ; c++) {
                //for each cluster of child node i
                SnarlSeedClusterer::cluster_t& child_cluster = 
                                                    child_clusters[i][c];
 
                hash_set<size_t> assignments;
                for (size_t k = 0; k < combined_clusters_indices.size() ; k++){
                    //For each existing cluster, see which this child cluster belongs to
                    vector<hash_set<size_t>>& old_cluster_indices =
                                            combined_clusters_indices[k].first;

                    for (size_t j = 0 ; j <= i ; j++){
                        //Go through child nodes up to i
                            
                        //Clusters of child j in this combined cluster
                        hash_set<size_t>& other_clusters = 
                                                     old_cluster_indices[j];
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
                        if (max({dist_l_l, dist_l_r, dist_r_l, dist_r_r}) != -1) {
                            for (size_t l : other_clusters) {
                                //For each cluster l of child j in this combined cluster k 
                                SnarlSeedClusterer::cluster_t& other_cluster
                                                     = child_clusters[j][l];
                                int64_t new_dist_l_l = dist_l_l == -1 || 
                                     child_cluster.dist_left == -1 ||
                                      other_cluster.dist_left == -1 ? -1 : 
                                          dist_l_l + child_cluster.dist_left + 
                                          other_cluster.dist_left;
                                int64_t new_dist_l_r = dist_l_r == -1 ||
                                         child_cluster.dist_left == -1 ||
                                         other_cluster.dist_right == -1 ? -1 : 
                                         dist_l_r + child_cluster.dist_left + 
                                          other_cluster.dist_right;
                                int64_t new_dist_r_l = dist_r_l == -1 ||
                                          child_cluster.dist_right == -1 ||
                                          other_cluster.dist_left == -1  ? -1 : 
                                          dist_r_l + child_cluster.dist_right + 
                                          other_cluster.dist_left;
                                int64_t new_dist_r_r = dist_r_r == -1 ||
                                         child_cluster.dist_right == -1 ||
                                         other_cluster.dist_right == -1 ? -1 : 
                                         dist_r_r + child_cluster.dist_right + 
                                          other_cluster.dist_right;
                                int64_t min_dist = 
                                     DistanceIndex::minPos({new_dist_l_l, 
                                     new_dist_l_r, new_dist_r_l, new_dist_r_r});
                                if (min_dist != -1 && min_dist < distance_limit) {
                                    //If these are in the same cluster
                                    assignments.insert(k);
                                    break;
                                } 
                            }
                        } 
                    
                    }
                }
                //find distances to ends of snarl for this cluster 
                int64_t new_dist_s_l = dist_s_l == -1  || 
                                       child_cluster.dist_left == -1 
                       ? -1 : dist_s_l + child_cluster.dist_left;
                int64_t new_dist_s_r = dist_s_r == -1 || 
                                        child_cluster.dist_right == -1 
                      ? -1 : dist_s_r + child_cluster.dist_right;
                int64_t new_dist_e_l = dist_e_l == -1 || 
                                       child_cluster.dist_left == -1
                       ? -1 : dist_e_l + child_cluster.dist_left;
                int64_t new_dist_e_r = dist_e_r == -1 || 
                                       child_cluster.dist_right == -1
                      ? -1 : dist_e_r + child_cluster.dist_right;
  
                pair<int64_t, int64_t> new_dists (
                         DistanceIndex::minPos({new_dist_s_l, new_dist_s_r}), 
                         DistanceIndex::minPos({new_dist_e_l,new_dist_e_r}));
                if (curr_node == snarl_index.snarlStart.first){
                    new_dists.first = snarl_index.snarlStart.second ? 
                            child_cluster.dist_right : child_cluster.dist_left ;
                }
                if (curr_node == snarl_index.snarlEnd.first) {
                    new_dists.second = snarl_index.snarlEnd.second ? 
                             child_cluster.dist_left : child_cluster.dist_right ;
                } 

                //Combine clusters that child cluster belongs to
                vector<pair<vector<hash_set<size_t>>, pair<int64_t, int64_t>>>
                                                     new_cluster_indices(1);
                new_cluster_indices.front().first.resize(i+1);

                //Combine clusters
                for (size_t k = 0; k < combined_clusters_indices.size() ; k++){
               
                    if (assignments.count(k) > 0) {
                        //If this cluster belongs with the current cluster
                        vector<hash_set<size_t>>& old_cluster_indices =
                                      combined_clusters_indices[k].first;
                        pair<int64_t, int64_t> old_cluster_dists = 
                                      combined_clusters_indices[k].second;
                        new_dists.first = DistanceIndex::minPos(
                                   {new_dists.first, old_cluster_dists.first});
                        new_dists.second = DistanceIndex::minPos(
                                   {new_dists.second,old_cluster_dists.second});
                        for (size_t j = 0 ; j <= i ; j++) {
                            //Combine clusters
                            if (old_cluster_indices[j].size() > 0) {
                                new_cluster_indices.front().first[j].insert(
                                      make_move_iterator(old_cluster_indices[j].begin()),
                                      make_move_iterator(old_cluster_indices[j].end()));
                            }
                        }
                    } else {
                        combined_clusters_indices[k].first.emplace_back();
                        new_cluster_indices.push_back(
                                                  move(combined_clusters_indices[k]));
                    } 
                    
                }
                new_cluster_indices.front().first[i].insert(c);
                new_cluster_indices.front().first.emplace_back();
                new_cluster_indices.front().second = new_dists;
                combined_clusters_indices = move(new_cluster_indices);
                
            }
        }
        
        //Create actual clusters from indices
        vector<SnarlSeedClusterer::cluster_t> snarl_clusters(combined_clusters_indices.size());
        for (size_t j = 0; j < combined_clusters_indices.size(); ++j) {
            pair<vector<hash_set<size_t>>, pair<int64_t, int64_t>>&
                           curr_cluster_indices = combined_clusters_indices[j];

            snarl_clusters[j].dist_left = rev ? 
                                      curr_cluster_indices.second.second : 
                                      curr_cluster_indices.second.first;
            snarl_clusters[j].dist_right = rev ?  
                                        curr_cluster_indices.second.first :
                                        curr_cluster_indices.second.second;
            for (size_t i = 0; i < curr_cluster_indices.first.size() ; i++) {
                for (size_t c : curr_cluster_indices.first[i]) {
                    if (snarl_clusters[j].seeds_start == nullptr) {
                            snarl_clusters[j].seeds_start = child_clusters[i][c].seeds_start; 
                            snarl_clusters[j].seeds_end = child_clusters[i][c].seeds_end; 
                        } else {
                            snarl_clusters[j].seeds_end->next = child_clusters[i][c].seeds_start;
                            snarl_clusters[j].seeds_end = child_clusters[i][c].seeds_end;
                        }
                }
            } 

        }
#ifdef DEBUG 
cerr << "Found clusters on snarl " << root->start() << endl;
for (cluster_t& c : snarl_clusters) {
    cerr << "\tleft: " << c.dist_left << " right : " << c.dist_right << "\tseeds: ";
    SnarlSeedClusterer::seed_index_list *s = c.seeds_start;
    while (s != nullptr){
        cerr << seeds[s->index] << "\t";
        s = s->next;
    }
    cerr << endl;
}
#endif
        return snarl_clusters;
    };
}
