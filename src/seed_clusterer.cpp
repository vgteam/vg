#include "seed_clusterer.hpp"

namespace vg {
    //TODO: Put this somewhere else???
    SnarlSeedClusterer::SnarlSeedClusterer() {
    };

    vector<size_t> SnarlSeedClusterer::cluster_seeds ( 
                  vector<pos_t> seeds, size_t distance_limit,
                  SnarlManager& snarl_manager, DistanceIndex& dist_index ){
        /* Given a vector of seeds and a limit, find a clustering of seeds where
         * seeds that are closer than the limit cluster together. 
         * Returns a vector of cluster assignments 
         */ 

        vector<size_t> cluster_assignments (seeds.size(), 0); 

        //Map snarl tree node to positions it contains, empty set if its child has seeds
        hash_set< const Chain* > chains_with_seeds;
        hash_set< const Snarl*> snarls_with_seeds; 
        hash_map< id_t, vector<size_t> > node_to_seed;
        for (size_t i = 0; i < seeds.size(); i++) {
            pos_t pos = seeds[i];
            id_t id = get_id(pos);
            const Snarl* snarl = dist_index.snarlOf(id);
            auto s = node_to_seed.find(id);
            if (s != node_to_seed.end()) {
                vector<size_t> new_seeds = s->second;
                new_seeds.push_back(i);
                node_to_seed.emplace(id, new_seeds);
            }
            while (snarl != NULL && snarls_with_seeds.count(snarl) == 0) {
                snarls_with_seeds.insert(snarl);
                if (snarl_manager.in_nontrivial_chain(snarl)) {
                    const Chain* chain = snarl_manager.chain_of(snarl);
                    chains_with_seeds.insert(chain);
                }
            }
        }
  
        const vector<const Snarl*>& top_snarls = 
                                               snarl_manager.top_level_snarls();
        unordered_set<const Snarl*> seen_snarls;
        for (const Snarl* root_snarl : top_snarls) {
            //Find clusters for each disconnected snarl/chain

            vector<SnarlSeedClusterer::cluster_t> clusters;
            if (seen_snarls.count(root_snarl) == 0) {
                if (snarl_manager.in_nontrivial_chain(root_snarl)){
                    const Chain* root_chain =snarl_manager.chain_of(root_snarl);
                    clusters =  get_clusters_chain( seeds, chains_with_seeds, 
                               snarls_with_seeds, node_to_seed, distance_limit,
                                snarl_manager, dist_index, root_chain);
                    for (auto s : *root_chain) {
                        seen_snarls.insert(s.first);
                    }
                } else {
                    clusters =  get_clusters_snarl( seeds, chains_with_seeds, 
                              snarls_with_seeds, node_to_seed, distance_limit, 
                              snarl_manager, dist_index, root_snarl, true);
                    seen_snarls.insert(root_snarl);
                }
            }

            for (size_t i = 0 ; i < clusters.size() ; i++) {
                SnarlSeedClusterer::cluster_t cluster = clusters[i]; 
                for (size_t seed_i : cluster.seeds) {
                    assert(cluster_assignments[seed_i] == 0 ||
                            cluster_assignments[seed_i] == i);//TODO: remove this, also the second part shouldn't be true but ok if it is
                    cluster_assignments[seed_i] = i;
                }
            }
        }
        return cluster_assignments;
        
    };


    vector<SnarlSeedClusterer::cluster_t> SnarlSeedClusterer::get_clusters_node(
                             vector<pos_t>& seeds,
                             hash_map< id_t, vector<size_t> >& node_to_seed,
                             size_t distance_limit, SnarlManager& snarl_manager,
                             DistanceIndex& dist_index, id_t root, bool rev,
                             int64_t node_length) {
        /*Find clusters of seeds in this node. rev is true if the left and right
         * distances should be reversed */ 
        vector<size_t> seed_indices = node_to_seed[root];                

        //vector of clusters where each cluster contains indices of seeds
        vector<hash_set<size_t>> cluster_indices;

        for ( size_t i : seed_indices) {
            //For each seed, see if it belongs to a new cluster
            pos_t seed = seeds[i]; 

            //Which clusters in cluster_indices this seed belongs to
            hash_set<size_t> curr_cluster_assignments;

            //Distance from seed to start of node 
            int64_t offset = is_rev(seed) ? node_length - get_offset(seed) : 
                                            get_offset(seed);
                        
            for (size_t j = 0 ; j < cluster_indices.size(); j++) {
                //Check which new clusters this seed belongs in 
                bool in_this_cluster = false;
                hash_set<size_t> curr_cluster = cluster_indices[j];
                for (size_t k : curr_cluster) {  
                    //Check current seed against seeds in curr cluster
                    if (!in_this_cluster) { 
                        pos_t other_seed = seeds[k];
                        int64_t other_offset = is_rev(other_seed) ?  
                                       node_length - get_offset(other_seed) : 
                                       get_offset(seed);
                        if (abs(offset - other_offset) < distance_limit) {
                            in_this_cluster = true;
                            curr_cluster_assignments.insert(j);
                        } 
                    }
                }
            }
            //Combine all clusters that the seed belongs to and make new cluster
            hash_set<size_t>new_cluster;
            new_cluster.insert(i);
            if (curr_cluster_assignments.size() == 0) {

                cluster_indices.push_back(new_cluster); 

            } else {

                vector<hash_set<size_t>> combined_cluster;
                for (size_t j = 0 ; j < cluster_indices.size(); j++) {

                    hash_set<size_t> old_cluster = cluster_indices[j];
                    if (curr_cluster_assignments.count(j) != 0) {
                        //If this seed belongs to this cluster

                        new_cluster.insert( old_cluster.begin(), 
                                            old_cluster.end());
                    } else {

                        combined_cluster.push_back(old_cluster);
                    }
                }
                cluster_indices = combined_cluster;
            }
        }

        //Create actual clusters from indices
        vector<SnarlSeedClusterer::cluster_t> clusters;
        for (hash_set<size_t> indices : cluster_indices) {
            SnarlSeedClusterer::cluster_t curr_cluster;
            for (size_t i : indices) {
                curr_cluster.seeds.insert(i);

                pos_t seed = seeds[i];
                int64_t dist_start = get_offset(seed) ; 
                int64_t dist_end = node_length - get_offset(seed);
                int64_t dist_left = rev == is_rev(seed) ? dist_start : dist_end;
                int64_t dist_right = rev == is_rev(seed) ? dist_end : 
                                                           dist_start;

                curr_cluster.dist_left = curr_cluster.dist_left == -1 ?
                             dist_left : min(dist_left, curr_cluster.dist_left);
                curr_cluster.dist_right = curr_cluster.dist_right == -1 ?
                          dist_right : min(dist_right, curr_cluster.dist_right);
            }
            clusters.push_back(curr_cluster);
        }
        return clusters;
        
    };

    vector<SnarlSeedClusterer::cluster_t> 
                  SnarlSeedClusterer::get_clusters_chain(
                             vector<pos_t>& seeds,
                             hash_set< const Chain* >& chains_with_seeds,
                             hash_set< const Snarl* >& snarls_with_seeds,
                             hash_map< id_t, vector<size_t> >& node_to_seed,
                             size_t distance_limit, SnarlManager& snarl_manager,
                             DistanceIndex& dist_index, const Chain* root) {
    
        vector<SnarlSeedClusterer::cluster_t> chain_clusters;

        DistanceIndex::ChainIndex& chain_index = dist_index.chainDistances.at(
                                                 get_start_of(*root).node_id());
        ChainIterator chain_e = chain_end(*root);
        for (ChainIterator c = chain_begin(*root); c != chain_e; ++c) {
            /* For each child snarl in the chain, find the clusters of just the
             * snarl, add the seeds that mapped to boundary nodes on the chain,
             * and progressively build up clusters spanning up to that snarl
             */
            const Snarl* curr_snarl = c->first;
            bool rev_in_chain = c->second;
            DistanceIndex::SnarlIndex& snarl_index = 
                        dist_index.snarlDistances.at(
                           make_pair(curr_snarl->start().node_id(),
                                     curr_snarl->start().backward()));


             //Start and end node of snarl relative to chain
            id_t start_node = rev_in_chain ? curr_snarl->end().node_id() :
                                                 curr_snarl->start().node_id();
            id_t end_node = rev_in_chain ? curr_snarl->start().node_id() :
                                                 curr_snarl->end().node_id();
            //TODO: This should be true if the node in the chain is traversed backward
            //TODO: Actually might not need this- always traverse a node forward in a chain????????
            id_t start_rev = rev_in_chain ? curr_snarl->end().backward() :
                                           ! curr_snarl->start().backward();
            id_t end_rev = rev_in_chain ? ! curr_snarl->start().backward() :
                                                curr_snarl->end().backward();
            int64_t start_length = snarl_index.nodeLength(start_node);
            int64_t end_length = snarl_index.nodeLength(end_node);
            if (c == chain_begin(*root)) {
                /* If this is the first snarl in the chain, then the 
                 * clusters in the first node of the chain can be considered
                 * as the clusters up to this snarl */ 

                chain_clusters =  get_clusters_node(seeds, node_to_seed,
                             distance_limit, snarl_manager, dist_index, 
                             start_node, false, start_length); 
            }
            if (snarls_with_seeds.count(curr_snarl) != 0) {
                //If this snarl has seeds on its children
                

                //Get clusters on current snarl and 2nd boundary node 
                vector<SnarlSeedClusterer::cluster_t> snarl_clusters =
                    get_clusters_snarl( seeds, chains_with_seeds,
                               snarls_with_seeds, node_to_seed, distance_limit,
                               snarl_manager, dist_index, curr_snarl, false);
   
                vector<SnarlSeedClusterer::cluster_t> boundary_clusters = 
                           get_clusters_node( seeds, node_to_seed,
                             distance_limit, snarl_manager, dist_index, 
                             end_node, false, end_length) ;
           

                //Vector of clusters, each cluster is pair of 
                //<cluster indices of node, cluster indices of snarl>
                vector<pair<hash_set<size_t>, hash_set<size_t>>> 
                                                     new_cluster_indices;
                for (size_t i = 0; i < snarl_clusters.size(); i++) {
                    /* Loop through each of the snarl's clusters and find 
                     * which boundary node clusters it belongs to, combining 
                     * clusters at each loop */

                    SnarlSeedClusterer::cluster_t cluster = snarl_clusters[i]; 
                    /* The shortest distance to the start node could actually be
                     * the shortest distance to the end node plus the distance 
                     * to loop in the chain back to the start node */ 
                    int64_t loop_dist_start = chain_index.loopFd[
                                            chain_index.snarlToIndex[end_node]] 
                               + snarl_index.snarlLength() - end_length;
                    int64_t loop_dist_end = chain_index.loopRev[
                                           chain_index.snarlToIndex[start_node]]
                               + snarl_index.snarlLength() - start_length; 
                    cluster.dist_left = min(cluster.dist_left,
                                          cluster.dist_right + loop_dist_start);
                    cluster.dist_right = min(cluster.dist_right,
                                          cluster.dist_left + loop_dist_end);
                    if (rev_in_chain) {
                        int64_t temp = cluster.dist_left;
                        cluster.dist_left = cluster.dist_right;
                        cluster.dist_right = temp;
                    }

                    //Which new clusters in new_cluster_indices 
                    //this cluster belongs to
                    hash_set<size_t> curr_cluster_assignments;

                    for (size_t j = 0 ; j < new_cluster_indices.size(); j++) {
                        //Check which new clusters this cluster belongs in 
                        bool in_this_cluster = false;
                        pair<hash_set<size_t>, hash_set<size_t>> 
                                     new_cluster = new_cluster_indices[j];
                        for (size_t k = 0; k < new_cluster.first.size(); k++) { 
                            SnarlSeedClusterer::cluster_t boundary_cluster = 
                                                           boundary_clusters[k];
                            //Check current seed against seeds in cluster
                            if (!in_this_cluster) {

                                if (abs((end_length - cluster.dist_right) -
                                        boundary_cluster.dist_left) 
                                      < distance_limit) {
                                        in_this_cluster = true;
                                        curr_cluster_assignments.insert(j);
                                } 
                            }
                        }
                        pair<hash_set<size_t>, hash_set<size_t>> comb_cluster;
                        comb_cluster.second.insert(i);
                        if (curr_cluster_assignments.size() == 0) {
                              new_cluster_indices.push_back(comb_cluster); 
                        } else {
                            vector<pair<hash_set<size_t>, hash_set<size_t>>>
                                                      combined_cluster_indices;
                            for (size_t j = 0 ; j < new_cluster_indices.size(); j++) {

                                pair<hash_set<size_t>, hash_set<size_t>> 
                                     old_cluster = new_cluster_indices[j];
                                if (curr_cluster_assignments.count(j) != 0) {
                                    //If this seed belongs to this cluster
                                    comb_cluster.first.insert(
                                                    old_cluster.first.begin(), 
                                                    old_cluster.first.end());
                                    comb_cluster.second.insert(
                                                    old_cluster.second.begin(), 
                                                    old_cluster.second.end());
                                 } else {
                                    combined_cluster_indices.push_back(old_cluster);
                                }
                            }
                            new_cluster_indices = combined_cluster_indices;
                        }
                    }
                }  
                //Build new clusters for this snarl that include boundary seeds
                vector<SnarlSeedClusterer::cluster_t> new_snarl_clusters;
                for (pair<hash_set<size_t>, hash_set<size_t>> cluster_indices :
                                                          new_cluster_indices) {
                    SnarlSeedClusterer::cluster_t curr_cluster;
                    int64_t boundary_offset = snarl_index.snarlLength() - 
                                              end_length;
                    for (size_t i : cluster_indices.first) {
                        //Go through boundary clusters in combined cluster and 
                        //combine them
                        SnarlSeedClusterer::cluster_t boundary_cluster = 
                                                           boundary_clusters[i];
                        curr_cluster.seeds.insert(boundary_cluster.seeds.begin(),
                                                  boundary_cluster.seeds.end());
                        int64_t old_dist_left = 
                                  boundary_cluster.dist_left + boundary_offset;
                        curr_cluster.dist_left = curr_cluster.dist_left == -1 ?
                                 old_dist_left : 
                                 min(old_dist_left, curr_cluster.dist_left);
                        curr_cluster.dist_right= curr_cluster.dist_right == -1 ?
                                 boundary_cluster.dist_right : 
                                 min(boundary_cluster.dist_right, 
                                                      curr_cluster.dist_right);
                    }
                    for (size_t i : cluster_indices.second) {
                        //Go through snarl clusters and combine
                        SnarlSeedClusterer::cluster_t old_cluster = 
                                                              snarl_clusters[i];
                        curr_cluster.seeds.insert(old_cluster.seeds.begin(),
                                                  old_cluster.seeds.end());

                        curr_cluster.dist_left = curr_cluster.dist_left == -1 ?
                             old_cluster.dist_left : 
                             min(old_cluster.dist_left, curr_cluster.dist_left);
                        curr_cluster.dist_right= curr_cluster.dist_right == -1 ?
                           old_cluster.dist_right : 
                           min(old_cluster.dist_right, curr_cluster.dist_right);
                    }
                    new_snarl_clusters.push_back(curr_cluster);
                }

//TODO: This is pretty repetative
//TODO: Could find distance to boundaries earlier
                //Combine clusters in chain up to this snarl with clusters in 
                //this snarl 
                //Each cluster is a pair of indices of old clusters and snarl clusters 
                vector<pair<hash_set<size_t>, hash_set<size_t>>> 
                                                      new_chain_cluster_indices;
                for (size_t i = 0 ; i < chain_clusters.size() ; i++) {
                    hash_set<size_t> old_clusters;
                    old_clusters.insert(i);
                    hash_set<size_t> new_clusters;
                    new_chain_cluster_indices.push_back(make_pair(
                                      old_clusters, new_clusters));
                }
                for (size_t j = 0; j < new_snarl_clusters.size() ; j++) {
                    /* For each of the clusters for the current snarl,
                     * find which chain clusters it belongs to and
                     * combine shared chain clusters */ 
                    SnarlSeedClusterer::cluster_t snarl_cluster = 
                                                          new_snarl_clusters[j];
                    hash_set<size_t> curr_cluster_assignments;
                                //Which chain cluster the snarl cluster is near 
                    bool in_this_cluster = false;
                    for (size_t i = 0; i < new_chain_cluster_indices.size();
                                        i ++ ) {
                        if (!in_this_cluster) {
                            pair<hash_set<size_t>, hash_set<size_t>> 
                                  curr_clusters = new_chain_cluster_indices[i];
                            for (size_t c : curr_clusters.first) {
                                SnarlSeedClusterer::cluster_t chain_cluster = 
                                           chain_clusters[c];         
                                if (chain_cluster.dist_right + 
                                        snarl_cluster.dist_left - start_length <
                                       distance_limit ) {
                                   //If two clusters are close enough to combine
                                    curr_cluster_assignments.insert(i);
                                    in_this_cluster = true;
                                }
                            }
                        }
                    }
                    vector<pair<hash_set<size_t>, hash_set<size_t>>> 
                                                              new_chain_indices;
                    pair<hash_set<size_t>, hash_set<size_t>> new_chain_cluster;
                    new_chain_cluster.second.insert(j);
                    for (size_t i = 0; i < new_chain_cluster_indices.size();
                                 i++ ) {
                        /* Combine chain clusters that share this snarl
                         * cluster*/
                        pair<hash_set<size_t>, hash_set<size_t>> old_cluster
                                             = new_chain_cluster_indices[i];
                        if (curr_cluster_assignments.count(i) != 0) {
                            new_chain_cluster.first.insert(
                                                 old_cluster.first.begin(),
                                                      old_cluster.first.end());
                            new_chain_cluster.second.insert(
                                                 old_cluster.second.begin(),
                                                  old_cluster.second.end());
                        } else {
                            new_chain_indices.push_back( old_cluster );
                        }
                    }
                    new_chain_indices.push_back(new_chain_cluster);
                    new_chain_cluster_indices = new_chain_indices;
                }
 
                /* Each element in new_chain_cluster_indices contains 
                 * the indices of each chain cluster and snarl cluster 
                 * that belong together. Combine these into actual clusters
                 */
 
                vector<SnarlSeedClusterer::cluster_t> new_chain_clusters;
                for (pair<hash_set<size_t>, hash_set<size_t>> 
                              new_chain_cluster : new_chain_cluster_indices) {

                    SnarlSeedClusterer::cluster_t curr_cluster;
                    for (size_t i : new_chain_cluster.first) {
                        SnarlSeedClusterer::cluster_t old_chain_cluster = 
                                                              chain_clusters[i];

                        curr_cluster.dist_left =  curr_cluster.dist_left == -1 ?
                                 old_chain_cluster.dist_left : 
                                 min(old_chain_cluster.dist_left, 
                                                 curr_cluster.dist_left);
                        int64_t old_dist_right =  old_chain_cluster.dist_right +
                                     snarl_index.snarlLength() - start_length;
                        curr_cluster.dist_right = 
                               curr_cluster.dist_right == -1 ?  old_dist_right :
                                min(old_dist_right, curr_cluster.dist_right);
                        curr_cluster.seeds.insert(
                                              old_chain_cluster.seeds.begin(),
                                              old_chain_cluster.seeds.end());
                    }
                    for (size_t j : new_chain_cluster.second) {
                        SnarlSeedClusterer::cluster_t old_snarl_cluster = 
                                                          new_snarl_clusters[j];
                        curr_cluster.dist_right =curr_cluster.dist_right == -1 ?
                                 old_snarl_cluster.dist_right : 
                                 min(old_snarl_cluster.dist_right, 
                                                 curr_cluster.dist_right);
                        int64_t old_dist_left = old_snarl_cluster.dist_left
                                  + chain_index.prefixSum[
                                      chain_index.snarlToIndex[start_node]];
                        curr_cluster.dist_left = 
                                 curr_cluster.dist_left == -1 ? old_dist_left : 
                                    min(old_dist_left, curr_cluster.dist_left);

                        curr_cluster.seeds.insert(
                                              old_snarl_cluster.seeds.begin(),
                                              old_snarl_cluster.seeds.end());
                    }
                    new_chain_clusters.push_back(curr_cluster);
                }
                chain_clusters = new_chain_clusters;
                

            } else { 
                //Extend the boundaries of the clusters to end of current snarl

                int64_t dist_to_next = snarl_index.snarlLength() - start_length;
                for (SnarlSeedClusterer::cluster_t curr_cluster :
                                                               chain_clusters) {

                    curr_cluster.dist_right = curr_cluster.dist_right + 
                                              dist_to_next;
                }
            }
            
        }
         
        return chain_clusters; 
    };

    vector<SnarlSeedClusterer::cluster_t> SnarlSeedClusterer::get_clusters_snarl(
                             vector<pos_t>& seeds,
                             hash_set< const Chain* >& chains_with_seeds,
                             hash_set< const Snarl* >& snarls_with_seeds,
                             hash_map< id_t, vector<size_t> >& node_to_seed,
                             size_t distance_limit, SnarlManager& snarl_manager,
                             DistanceIndex& dist_index, const Snarl* root,
                             bool include_boundaries) {
    
        DistanceIndex::SnarlIndex& snarl_index = dist_index.snarlDistances.at(
                           make_pair(root->start().node_id(),
                                     root->start().backward()));
//TODO: Find a maybe better way to classify all nodes in the snarl
        vector<id_t> child_nodes;
        vector<const Snarl*> child_snarls; 
        vector<const Chain*> child_chains;
        for (auto v : snarl_index.visitToIndex) {

            pair<id_t, bool> visit = v.first;
            if (visit.second) {
                const Snarl* snarl = snarl_manager.into_which_snarl(visit.first,
                                                                    false);
                if (snarl == NULL) {
                    snarl = snarl_manager.into_which_snarl(visit.first, true);
                }
                if (snarl == NULL) {
                    //If this is a node and 
                    //if not including boundary nodes, not a boundary node
                    if ( include_boundaries || 
                               (visit.first != snarl_index.snarlStart.first &&
                                   visit.first != snarl_index.snarlEnd.first)) {
                        child_nodes.push_back(visit.first);
                    }
                } else {
                    //This is a snarl or a chain
                    if (snarl_manager.in_nontrivial_chain(snarl)) {
                        //This is a chain
                        const Chain* chain = snarl_manager.chain_of(snarl); 
                        child_chains.push_back(chain);
                    } else {
                        //This is a snarl
                        child_snarls.push_back(snarl);
                    }
                }
            }
        }
        size_t num_children = child_nodes.size() + child_snarls.size() + 
                              child_chains.size();

        //clusters of children, indices assume child node, snarl, and chain
        //vectors are contiguous 
        vector<vector<SnarlSeedClusterer::cluster_t>> child_clusters;
        //Vector of clusters where a cluster is represented by indices into
        //child_nodes/snarls/chains
        vector<pair<vector<hash_set<size_t>>, pair<int64_t, int64_t>>>
                                                     combined_clusters_indices;
         
        for (size_t i = 0; i < num_children ; i++) {
            vector<SnarlSeedClusterer::cluster_t> curr_clusters; 
            id_t curr_node;
            bool curr_rev;
            if (i < child_nodes.size()) {
                curr_node = child_nodes[i];
                curr_rev = false;
                int64_t node_len = snarl_index.nodeLength(curr_node);
                if (node_to_seed.find(curr_node) != node_to_seed.end()) {
                    //If this node has seeds, cluster them
                    curr_clusters = get_clusters_node(
                             seeds, node_to_seed,
                             distance_limit, snarl_manager, dist_index, 
                             curr_node, false, node_len); 
                }

            } else if (i < child_snarls.size()) {

                const Snarl* curr_snarl = child_snarls[i - child_nodes.size()];
                curr_node = curr_snarl->start().node_id();
                curr_rev = curr_snarl->start().backward();
                if (snarls_with_seeds.count(curr_snarl) != 0) {
                    curr_clusters = get_clusters_snarl(
                         seeds, chains_with_seeds, snarls_with_seeds, 
                         node_to_seed,distance_limit, snarl_manager, dist_index,
                         curr_snarl, true);
                }
            } else {
                const Chain* curr_chain =  child_chains[i - child_nodes.size() -
                                             child_snarls.size()];
                Visit start = get_start_of(*curr_chain);
                curr_node = start.node_id();
                curr_rev = start.backward();
                if (chains_with_seeds.count(curr_chain) != 0) {
                    curr_clusters = get_clusters_chain(seeds, chains_with_seeds,
                             snarls_with_seeds, node_to_seed, distance_limit,
                             snarl_manager, dist_index, curr_chain);
                }
            }
            child_clusters.push_back(curr_clusters);

            for (size_t c = 0 ; c < curr_clusters.size() ; c++) {
                //for each cluster of child i
                SnarlSeedClusterer::cluster_t child_cluster = curr_clusters[c];
 
                hash_set<size_t> assignments;
                for (size_t k = 0; k < combined_clusters_indices.size() ; k++){
                    //For each existing cluster, see which this child cluster belongs to
                    bool in_this_cluster = false;
                    vector<hash_set<size_t>> old_cluster_indices =
                                            combined_clusters_indices[k].first;
                    if (old_cluster_indices.size() < i-1) {
                        //If this cluster doesn't have a set for this child node
                        //TODO: Maybe a better way of doing this??

                        hash_set<size_t> new_cluster ;
                        old_cluster_indices.push_back(new_cluster);
                    }
                    for (size_t j = 0 ; j < i-1; j++){
                        if (!in_this_cluster) {
                            
                            //Clusters of child j in this combined cluster
                            hash_set<size_t> other_clusters = 
                                                         old_cluster_indices[j];
                            id_t other_node;
                            bool other_rev;
                            if (j < child_nodes.size()) {
                                other_node = child_nodes[j];
                                other_rev = false;
                            } else if (j < child_snarls.size()) {
                                const Snarl* other_snarl =  
                                           child_snarls[j - child_nodes.size()];
                                other_node = other_snarl->start().node_id();
                            } else {
                                const Chain* other_chain = child_chains [
                                  j - child_nodes.size() - child_snarls.size()];
                                Visit start = get_start_of(*other_chain);
                                other_node = start.node_id();
                                other_rev = start.backward();
                            }
                            //Find distance from each end of current node to 
                            //each end of other node
                            int64_t dist_l_l = snarl_index.distances[
                              snarl_index.index(
                                     make_pair(curr_node, !curr_rev), 
                                     make_pair(other_node, other_rev))];
                            int64_t dist_l_r = snarl_index.distances[
                              snarl_index.index(
                                     make_pair(curr_node, !curr_rev), 
                                     make_pair(other_node, !other_rev))];
                            int64_t dist_r_l = snarl_index.distances[   
                                  snarl_index.index(
                                      make_pair(curr_node, curr_rev), 
                                      make_pair(other_node, other_rev))];
                            int64_t dist_r_r = snarl_index.distances[
                                    snarl_index.index(
                                      make_pair(curr_node, curr_rev), 
                                      make_pair(other_node, !other_rev))];
                            for (size_t l : other_clusters) {
                                SnarlSeedClusterer::cluster_t other_cluster = 
                                                           child_clusters[j][l];
                                dist_l_l = dist_l_l == -1 ? -1 : dist_l_l + 
                                              child_cluster.dist_left + 
                                              other_cluster.dist_left;
                                dist_l_r = dist_l_r == -1 ? -1 : dist_l_r + 
                                              child_cluster.dist_left + 
                                              other_cluster.dist_right;
                                dist_r_l = dist_r_l == -1 ? -1 : dist_r_l + 
                                              child_cluster.dist_right + 
                                              other_cluster.dist_left;
                                dist_r_r = dist_r_r == -1 ? -1 : dist_r_r + 
                                              child_cluster.dist_right + 
                                              other_cluster.dist_right;
                                int64_t min_dist = 
                                     DistanceIndex::minPos({dist_l_l, dist_l_r,
                                                 dist_r_l, dist_r_r});
                                if (min_dist != -1 && min_dist < distance_limit) {
                                    //If these are in the same cluster
                                    in_this_cluster = true;
                                    assignments.insert(k);
                                } 
                            } 
                        }
                    }
                }
                //Combine clusters that child cluster belongs to
                vector<pair<vector<hash_set<size_t>>, pair<int64_t, int64_t>>>
                                                     new_cluster_indices;
                vector<hash_set<size_t>> combined_cluster;
                for (size_t j = 0; j < i-1 ; j++) {
                    hash_set<size_t> set;
                    combined_cluster.push_back(set);
                }
                //find distances to ends of snarl for this cluster 
                int64_t dist_s_l = snarl_index.distances[
                                      snarl_index.index(snarl_index.snarlStart,
                                              make_pair(curr_node, !curr_rev))];
                int64_t dist_s_r = snarl_index.distances[
                                   snarl_index.index(snarl_index.snarlStart,
                                              make_pair(curr_node, curr_rev))];
                int64_t dist_e_l = snarl_index.distances[
                                      snarl_index.index(snarl_index.snarlEnd,
                                              make_pair(curr_node, !curr_rev))];
                int64_t dist_e_r = snarl_index.distances[
                                     snarl_index.index(snarl_index.snarlEnd,
                                             make_pair(curr_node, curr_rev))];
                dist_s_l = dist_s_l == -1 ? -1 : dist_s_l +
                                                child_cluster.dist_left;
                dist_s_r = dist_s_r == -1 ? -1 : dist_s_r +
                                                child_cluster.dist_right;
                dist_e_l = dist_e_l == -1 ? -1 : dist_e_l +
                                                child_cluster.dist_left;
                dist_e_r = dist_e_r == -1 ? -1 : dist_e_r +
                                                child_cluster.dist_right;
                pair<int64_t, int64_t> new_dists (
                         DistanceIndex::minPos({dist_s_l, dist_s_r}), 
                         DistanceIndex::minPos({dist_e_l,dist_e_r}));

                //Combine clusters
                for (size_t k = 0; k < combined_clusters_indices.size() ; k++){
               
                    if (assignments.count(k) > 0) {
                        vector<hash_set<size_t>> old_cluster_indices =
                                      combined_clusters_indices[k].first;
                        pair<int64_t, int64_t> old_cluster_dists = 
                                      combined_clusters_indices[k].second;
                        for (size_t j = 0 ; j < i-1 ; j++) {
                            //Combine clusters
                            combined_cluster[j].insert(
                                old_cluster_indices[j].begin(),
                                old_cluster_indices[j].end());
                            new_dists.first = DistanceIndex::minPos(
                                   {new_dists.first, old_cluster_dists.first});
                            new_dists.second = DistanceIndex::minPos(
                                   {new_dists.second,old_cluster_dists.second});
                        }
                    } else {

                        new_cluster_indices.push_back(
                                                  combined_clusters_indices[k]);
                    } 
                    
                }
                combined_cluster[i].insert(c);
                new_cluster_indices.push_back(make_pair(combined_cluster, 
                                                         new_dists));
                combined_clusters_indices = new_cluster_indices;
                
            }
        }
        
        //Create actual clusters from indices
        vector<SnarlSeedClusterer::cluster_t> snarl_clusters;
        for (pair<vector<hash_set<size_t>>, pair<int64_t, int64_t>>
                           curr_cluster_indices : combined_clusters_indices) {

            SnarlSeedClusterer::cluster_t cluster;
            cluster.dist_left = curr_cluster_indices.second.first;
            cluster.dist_right = curr_cluster_indices.second.second;
            for (size_t i = 0; i < curr_cluster_indices.first.size() ; i++) {
                for (size_t c : curr_cluster_indices.first[i]) {
                    cluster.seeds.insert(child_clusters[i][c].seeds.begin(),
                                          child_clusters[i][c].seeds.end() );
                }
            } 

            snarl_clusters.push_back(cluster); 
        }
        return snarl_clusters;
    };
}
