#include "seed_clusterer.hpp"

#define DEBUG 

namespace vg {
    //TODO: Put this somewhere else???
    SnarlSeedClusterer::SnarlSeedClusterer() {
    };

    vector<hash_set<size_t>> SnarlSeedClusterer::cluster_seeds ( 
                  vector<pos_t> seeds, size_t distance_limit,
                  SnarlManager& snarl_manager, DistanceIndex& dist_index ){
        /* Given a vector of seeds and a limit, find a clustering of seeds where
         * seeds that are closer than the limit cluster together. 
         * Returns a vector of cluster assignments 
         */ 


        /*Store the subtree of the snarl tree that contains seeds*/
        
        //Maps chains to the snarls and boundary nodes containing seeds
        //Snarls are paired with its second boundary node in order of chain,
        //First boundary node in the chain has NULL snarl
        //Snarl/node pairs are mapped to by the rank of node in the chain
        hash_map< const Chain*, map<size_t, pair<const Snarl*, int64_t>>>
                                                chains_to_snarl_rank;
        //This will be filled in based on chains_to_snarl_rank
        //snarl/boundary nodes will be ordered by order in chain
        //boundary node is node id of node, -1 for placeholder value
        hash_map< const Chain*, vector<pair<const Snarl*, int64_t>>>
                                                            chains_to_snarl;
        //Maps each snarl to the nodes in its netgraph that contain seeds
        hash_map< const Snarl*, hash_set<pair<id_t, bool>>> snarls_to_node; 

        //Maps nodes to the indexes of the seeds it contains
        hash_map< id_t, vector<size_t> > node_to_seed;

        vector<SnarlSeedClusterer::seed_index_list> seed_list;
        for (size_t i = 0; i < seeds.size(); i++) {
            //For each seed, add its containing snarls and chains to the
            //subtree of the snarl tree
            pos_t pos = seeds[i];
            id_t id = get_id(pos);
            SnarlSeedClusterer::seed_index_list seed_l;
            seed_l.index = i;
            seed_list.push_back(seed_l);
            const Snarl* snarl = dist_index.snarlOf(id);
            auto s = node_to_seed.find(id);
            if (s != node_to_seed.end()) {
                s->second.push_back(i);
            } else {
                vector<size_t> new_seeds({i});
                node_to_seed.emplace(id, new_seeds);
            }
            pair<id_t, bool> prev (id, false);
cerr << endl << "Ancestors of " << id << endl;
            while (snarl != nullptr) { 
                //Add ancestors of snarl to the snarl subtree
cerr << snarl->start() << ", ";
                if (snarls_to_node.count(snarl) == 0) {
                    hash_set<pair<id_t, bool>> nodes;
                    nodes.insert(prev);
                    snarls_to_node.emplace(snarl, nodes);
                } else {
                    snarls_to_node.at(snarl).insert(prev);
                } 
                prev = make_pair(snarl->start().node_id(), 
                                 snarl->start().backward());
                if (snarl_manager.in_nontrivial_chain(snarl)) {
                    //Add chain to snarl subtree
                    const Chain* chain = snarl_manager.chain_of(snarl);
                    bool rev_in_chain = snarl_manager.chain_orientation_of(snarl);
                    //Rank in chain of the second boundary node of the snarl 
                    //relative to the order of the chain
                    size_t rank = snarl_manager.chain_rank_of(snarl);//TODO: Change this to get the rank from the distance index instead - rank of first node in snarl relative to orientation in chain 

//TODO: Could not add snarl if seed is on boundary node in a chain?
                    bool on_bound = false;//seed on boundary node, not in snarl
                    if ((id == snarl->start().node_id() && !rev_in_chain) ||
                        (id == snarl->end().node_id() && rev_in_chain)){ 
                        //If seed is on first boundary node of snarl in chain
                        on_bound = true;
                    } else {
                        rank += 1;//rank is order of second end of snarl
                         if ((id == snarl->start().node_id() && rev_in_chain)
                        || (id == snarl->end().node_id() && !rev_in_chain)){ 
                            //If seed is on second boundary node of snarl in chain
                            on_bound = true;
                        }
                    }
                    pair<const Snarl*, id_t> node_with_seed (nullptr, -1);
                    if (on_bound) {
                        node_with_seed.second = id;
                    } else {
                        node_with_seed.first = snarl;
                    }

                    auto found_chain = chains_to_snarl_rank.find(chain);
                    if (found_chain == chains_to_snarl_rank.end()) {
                        //If this chain has no seeds yet
                        map<size_t, pair<const Snarl*, int64_t>> chain_snarls;
                        chain_snarls.emplace(rank, node_with_seed);
                        chains_to_snarl_rank.emplace(chain, chain_snarls);
                    } else {

                        map<size_t, pair<const Snarl*, int64_t>> 
                                   old_chain_snarls = move(found_chain->second);

                        pair<const Snarl*, id_t> old_node_with_seed = 
                                 old_chain_snarls.count(rank) == 0 ?
                                make_pair(nullptr, -1) : old_chain_snarls[rank];


                        pair<const Snarl*, id_t> new_node_with_seed (
                               node_with_seed.first == nullptr ? 
                                         old_node_with_seed.first : 
                                         node_with_seed.first, 
                                node_with_seed.second == -1 ? 
                                          old_node_with_seed.second :
                                          node_with_seed.second);

                        old_chain_snarls.erase(rank);
                        old_chain_snarls.emplace(rank, new_node_with_seed); 
                        chains_to_snarl_rank.erase(chain);
                        chains_to_snarl_rank.emplace(chain, old_chain_snarls);
                        //TODO:Not sure why i need to replace items in the map 
                    }
                }
                snarl = snarl_manager.parent_of(snarl);
            }
        }

        //Get chains_to_snarl from chains_to_snarl_rank
        for (pair<const Chain*, map<size_t, pair<const Snarl*, int64_t>>> 
                 ranked_snarls : chains_to_snarl_rank){
            vector<pair<const Snarl*, int64_t>> snarl_list;
            for (pair<size_t, pair<const Snarl*, int64_t>> ranked : ranked_snarls.second){
                //TODO: Traverse map of rank -> (snarl, boundary node) and add
                //snarl, boundary node to list
                snarl_list.push_back( ranked.second );
            } 
            chains_to_snarl.emplace(ranked_snarls.first, snarl_list);
        }
#ifdef DEBUG

cerr << "CHAINS: " << endl;
for (pair<const Chain*, vector<pair<const Snarl*, int64_t>>> c : chains_to_snarl) {
    cerr << get_start_of(*c.first).node_id() << ": ";
    for (pair<const Snarl*, int64_t> s : c.second) {
        id_t s_id = s.first == nullptr ? 0 : s.first->start().node_id();
        cerr << "snarl: " << s_id << ", node: " << s.second << " ";
    }
    cerr << endl;
}

cerr << "SNARLS: " << endl; 
for (pair<const Snarl*, hash_set<pair<id_t, bool>>> s : snarls_to_node) {
    cerr << s.first->start() << "  : " ;
    for (pair<id_t, bool> n : s.second) {
        cerr << n.first << " " << n.second << ", "; 
    }
    cerr << endl;
}

cerr << "NODES: " << endl;
for (pair<id_t, vector<size_t>> n : node_to_seed){
    cerr << n.first << ": " << endl;
    for (size_t s : n.second ) {
        cerr << seeds[s] << ", ";
    }
    cerr << endl;
} 
#endif
  
        // This will hold all top-level snarls that have seeds under them
        vector<const Snarl*> top_snarls;
        for (auto& snarl : snarls_to_node) {
            if (snarl_manager.parent_of(snarl.first) == nullptr) {
                // This is top level
                top_snarls.push_back(snarl.first);
            }
        }
        
        // We track seen-ness by chain, so we don't have to do O(N) work to tag
        // all snarls in chromosome-spanning chains as seen.
        unordered_set<const Chain*> seen_chains;
        vector<hash_set<size_t>> cluster_assignments; 
        for (const Snarl* root_snarl : top_snarls) {
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
                              snarl_manager, dist_index, root_snarl, true);
                }
                seen_chains.insert(root_chain);
            }

            for (size_t i = 0 ; i < clusters.size() ; i++) {
                SnarlSeedClusterer::cluster_t cluster = clusters[i]; 
                hash_set<size_t> seeds;
                SnarlSeedClusterer::seed_index_list* j = cluster.seeds_start;
                while (j != nullptr) {
                    seeds.insert(j->index);
                    j = j->next;
                }
                cluster_assignments.push_back(seeds); 
            }
        }
        return cluster_assignments;
        
    };


    vector<SnarlSeedClusterer::cluster_t> SnarlSeedClusterer::get_clusters_node(
                        vector<pos_t>& seeds,
                        vector<SnarlSeedClusterer::seed_index_list>& seed_list, 
                        hash_map< id_t, vector<size_t> >& node_to_seed,
                        size_t distance_limit, SnarlManager& snarl_manager,
                        DistanceIndex& dist_index, id_t root, bool rev,
                        int64_t node_length) {
#ifdef DEBUG 
cerr << "Finding clusters on node " << root << endl;
#endif
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
                                            get_offset(seed)+1;
                        
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
                                       get_offset(other_seed) + 1;
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
                cluster_indices = move(combined_cluster);
            }
        }

        //Create actual clusters from indices
        vector<SnarlSeedClusterer::cluster_t> clusters;
        for (hash_set<size_t> indices : cluster_indices) {
            SnarlSeedClusterer::cluster_t curr_cluster;
            for (size_t i : indices) {
                SnarlSeedClusterer::seed_index_list* curr_seed_i = &seed_list[i];
                if (curr_cluster.seeds_start == nullptr) {
                    curr_cluster.seeds_start = curr_seed_i; 
                    curr_cluster.seeds_end = curr_seed_i; 
                } else {
                    curr_cluster.seeds_end->next = curr_seed_i;
                    curr_cluster.seeds_end = curr_seed_i;
                }

                pos_t seed = seeds[i];
                int64_t dist_start = get_offset(seed) + 1; 
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
#ifdef DEBUG 
cerr << "Found clusters on node " << root << endl;
for (cluster_t c : clusters) {
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
                         vector<pos_t>& seeds,
                         vector<SnarlSeedClusterer::seed_index_list>& seed_list,
                         hash_map< const Chain*, 
                                  vector<pair<const Snarl*, int64_t>> >& 
                                              chains_to_snarl,
                         hash_map< const Snarl*, hash_set<pair<id_t, bool>> >& 
                                              snarls_to_node,
                         hash_map< id_t, vector<size_t> >& node_to_seed,
                         size_t distance_limit, SnarlManager& snarl_manager,
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
        vector<pair<const Snarl*, int64_t>> snarls_in_chain = chains_to_snarl.at(root);

        for (pair<const Snarl*, int64_t> snarl_node : snarls_in_chain) {
            /* For each child snarl in the chain, find the clusters of just the
             * snarl, add the seeds that mapped to boundary nodes on the chain,
             * and progressively build up clusters spanning up to that snarl
             */
            const Snarl* curr_snarl = snarl_node.first;
            int64_t boundary_node_id = snarl_node.second;

            //Clusters of the current snarl/boundry node to be combined with
            //existing chain clusters
            vector<SnarlSeedClusterer::cluster_t> curr_clusters;
            int64_t dist_to_end;//Dist to update chain clsters to end of new cluster 
            int64_t start_length; // length of first node in current cluster
            if (curr_snarl != nullptr){
                /* If the snarl has seeds on it, 
                 * extend the existing clusters
                 *    up to the end of the first node in the snarl (relative to 
                 *    traversal of the chain),
                 * find clusters of seeds on this snarl 
                 * if the boundary node of snarl has seeds, combine snarl
                 *     clusters with clusters of seeds on the node -> 
                 *     curr_clusters 
                 */
                bool rev_in_chain = snarl_manager.chain_orientation_of(
                                                                    curr_snarl);

                //TODO: This should be true if the node in the chain is traversed backward
                //TODO: Actually might not need this- always traverse a node forward in a chain????????
                id_t start_rev = rev_in_chain ? curr_snarl->end().backward() :
                                               ! curr_snarl->start().backward();
                id_t end_rev = rev_in_chain ? ! curr_snarl->start().backward() :
                                                 curr_snarl->end().backward();
             //Start and end node of snarl relative to chain
                id_t start_node = rev_in_chain ? curr_snarl->end().node_id() :
                                                 curr_snarl->start().node_id();
                id_t end_node = rev_in_chain ? curr_snarl->start().node_id() :
                                                 curr_snarl->end().node_id();
                DistanceIndex::SnarlIndex& snarl_index = 
                            dist_index.snarlDistances.at(
                               make_pair(curr_snarl->start().node_id(),
                                         curr_snarl->start().backward()));
                start_length = snarl_index.nodeLength(start_node);
                int64_t end_length = snarl_index.nodeLength(end_node);


                if (last_snarl != start_node) { 
                    /* If the chain clusters don't reach this snarl,
                     * extend their dist_right to the beginning of this snarl
                     */  
                    int64_t offset = chain_index.prefixSum[
                                      chain_index.snarlToIndex[start_node]] -
                                chain_index.prefixSum[
                                      chain_index.snarlToIndex[last_snarl]] +
                                start_length - last_len;
                    vector<cluster_t> new_chain_clusters;
                    for (cluster_t curr : chain_clusters) {
                        //Update the clusters to reach the current snarl
                        cluster_t new_c;
                        new_c.dist_right = curr.dist_right + offset;
                        new_c.dist_left = curr.dist_left;
                        if (new_c.seeds_start == nullptr) {
                            new_c.seeds_start = curr.seeds_start; 
                            new_c.seeds_end = curr.seeds_end; 
                        } else {
                            new_c.seeds_end->next = curr.seeds_start;
                            new_c.seeds_end = curr.seeds_end;
                        }
                        new_chain_clusters.push_back(new_c);
                    }
                    chain_clusters = move(new_chain_clusters);
                }
                last_snarl = end_node;
                last_len = end_length;

                if ( boundary_node_id == -1) {
                    //If there are seeds on the snarl but not the boundary node
                    //Get clusters on just the snarl 
                    curr_clusters = get_clusters_snarl( 
                               seeds, seed_list, chains_to_snarl,
                               snarls_to_node, node_to_seed, distance_limit,
                               snarl_manager, dist_index, curr_snarl, false);
                } else {
                    //If there are seeds on both the snarl and boundary node,
                    //combine these into curr_clusters
                    vector<SnarlSeedClusterer::cluster_t> snarl_clusters =
                         get_clusters_snarl( seeds, seed_list, chains_to_snarl,
                               snarls_to_node, node_to_seed, distance_limit,
                               snarl_manager, dist_index, curr_snarl, false);
                     vector<SnarlSeedClusterer::cluster_t> boundary_clusters = 
                          get_clusters_node( seeds, seed_list, node_to_seed,
                             distance_limit, snarl_manager, dist_index, 
                             end_node, false, end_length) ;

                    //Vector of clusters, each cluster is pair of 
                    //<index into boundary_clusters, index into snarl_clusters>
                    vector<pair<hash_set<size_t>, hash_set<size_t>>> 
                                                     new_cluster_indices;
                    for (size_t b = 0; b < boundary_clusters.size() ; b++) {
                        //These are all from the same node so none of the 
                        //clusters should be combined
                        hash_set<size_t> s ;
                        s.insert(b);
                        new_cluster_indices.push_back(make_pair(
                                                    s, hash_set<size_t>()));
                    }
                    for (size_t i = 0; i < snarl_clusters.size(); i++) {
                        /* Loop through each of the snarl's clusters and find 
                         * which boundary node clusters it belongs with, 
                         * combining clusters at each loop */

                        SnarlSeedClusterer::cluster_t cluster = snarl_clusters[i]; 
                        if (rev_in_chain) {
                            //Orient this cluster to be the same as the chain
                            int64_t temp = cluster.dist_left;
                            cluster.dist_left = cluster.dist_right;
                            cluster.dist_right = temp;
                        }
                        /* The shortest distance to the start node could 
                         * actually be the shortest distance to the end node 
                         * plus the distance to loop in the chain back to the
                         * start node */ 
                        int64_t loop_dist_start = chain_index.loopFd[
                                       chain_index.snarlToIndex[end_node]] - 1 ;
                        int64_t loop_dist_end = chain_index.loopRev[
                                     chain_index.snarlToIndex[start_node]] - 1; 
                        cluster.dist_left = loop_dist_start == -1 
                                         ? cluster.dist_left 
                                         : min(cluster.dist_left,
                                          cluster.dist_right + loop_dist_start 
                               + snarl_index.snarlLength() - end_length);
                        cluster.dist_right = loop_dist_end == -1 
                                      ? cluster.dist_right 
                                      : min(cluster.dist_right,
                                          cluster.dist_left + loop_dist_end
                               + snarl_index.snarlLength() - start_length);

                        //Which new clusters in new_cluster_indices 
                        //this cluster belongs to
                        hash_set<size_t> curr_cluster_assignments;

                        for (size_t j = 0 ; j < new_cluster_indices.size(); j++) {
                            //Check which new clusters this cluster belongs in 
                            bool in_this_cluster = false;
                            pair<hash_set<size_t>, hash_set<size_t>> 
                                         new_cluster = new_cluster_indices[j];
                            for (size_t k = 0; k < new_cluster.first.size(); k++) { 
                                //Go through each boundary cluster in this combined cluster
                                SnarlSeedClusterer::cluster_t boundary_cluster = 
                                                           boundary_clusters[k];
                                //Check current seed against seeds in cluster
                                if (!in_this_cluster) {

                                    if (abs(cluster.dist_right -
                                            boundary_cluster.dist_right) 
                                          < distance_limit) {
                                        in_this_cluster = true;
                                        curr_cluster_assignments.insert(j);
                                    } 
                                }
                            }
                        }

                        //For each of the boundary clusters that this snarl 
                        //cluster belongs to, combine the clusters
                        pair<hash_set<size_t>, hash_set<size_t>> comb_cluster;
                        comb_cluster.second.insert(i);
                        if (curr_cluster_assignments.size() == 0) {
                              //If this snarl clusters is on its own,
                              //create new cluster
                              new_cluster_indices.push_back(comb_cluster); 
                        } else {
                            vector<pair<hash_set<size_t>, hash_set<size_t>>>
                                                      combined_cluster_indices;
                            for (size_t j = 0 ; j < new_cluster_indices.size() ; j++){

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
                            combined_cluster_indices.push_back(comb_cluster);
                            new_cluster_indices = move(combined_cluster_indices);

                        }
                    }  
                    //Build new clusters for this snarl that include boundary seeds
                    for (pair<hash_set<size_t>, hash_set<size_t>> 
                                       cluster_indices : new_cluster_indices) {
                        SnarlSeedClusterer::cluster_t curr_cluster;
                        int64_t boundary_offset = snarl_index.snarlLength() - 
                                                  end_length;
                        for (size_t i : cluster_indices.first) {
                            //Go through boundary clusters in combined cluster and 
                            //combine them
                            SnarlSeedClusterer::cluster_t boundary_cluster = 
                                                           boundary_clusters[i];
                            if (curr_cluster.seeds_start == nullptr) {
                                curr_cluster.seeds_start =
                                                   boundary_cluster.seeds_start;
                                curr_cluster.seeds_end =
                                                    boundary_cluster.seeds_end; 
                            } else {
                                curr_cluster.seeds_end->next = 
                                                   boundary_cluster.seeds_start;
                                curr_cluster.seeds_end = boundary_cluster.seeds_end;
                            }
                            int64_t old_dist_left = 
                                  boundary_cluster.dist_left + boundary_offset;
                            curr_cluster.dist_left = 
                                    curr_cluster.dist_left == -1 ?
                                     old_dist_left : 
                                     min(old_dist_left, curr_cluster.dist_left);
                            curr_cluster.dist_right = 
                                     curr_cluster.dist_right == -1 ?
                                     boundary_cluster.dist_right : 
                                     min(boundary_cluster.dist_right, 
                                                      curr_cluster.dist_right);
                        }
                        for (size_t i : cluster_indices.second) {
                            //Go through snarl clusters and combine
                            SnarlSeedClusterer::cluster_t old_cluster = 
                                                              snarl_clusters[i];

                            if (curr_cluster.seeds_start == nullptr) {
                                curr_cluster.seeds_start = old_cluster.seeds_start; 
                                curr_cluster.seeds_end = old_cluster.seeds_end; 
                            } else {
                                curr_cluster.seeds_end->next = 
                                                   old_cluster.seeds_start;
                                curr_cluster.seeds_end = old_cluster.seeds_end;
                            }

                            curr_cluster.dist_left = curr_cluster.dist_left == -1 ?
                                 old_cluster.dist_left : 
                                 min(old_cluster.dist_left, curr_cluster.dist_left);
                            curr_cluster.dist_right= curr_cluster.dist_right == -1 ?
                               old_cluster.dist_right : 
                               min(old_cluster.dist_right, curr_cluster.dist_right);
                        }

                        curr_clusters.push_back(curr_cluster);
                    }

                }
                dist_to_end = snarl_index.snarlLength() - start_length;

            } else {
                /* If there is no snarl- only boundary node
                 * Find clusters on boundary node -> curr_clusters
                 * Extend existing chain clusters to the end of boundary node
                 */ 
                const Snarl* sn = dist_index.snarlOf(boundary_node_id);
                DistanceIndex::SnarlIndex& snarl_index = 
                      dist_index.snarlDistances.at(make_pair(
                               sn->start().node_id(), sn->start().backward()));

                int64_t node_length = snarl_index.nodeLength(boundary_node_id);
                curr_clusters = 
                     get_clusters_node( seeds, seed_list, node_to_seed,
                             distance_limit, snarl_manager, dist_index, 
                             boundary_node_id, false, node_length) ;
                if (last_snarl != boundary_node_id) { 
                    /* If the chain clusters don't reach this node,
                     * extend their dist_right to the end of this node
                     */  
                    int64_t offset = chain_index.prefixSum[
                                chain_index.snarlToIndex[boundary_node_id]] -
                                chain_index.prefixSum[
                                      chain_index.snarlToIndex[last_snarl]] +
                                node_length - last_len;
                    vector<cluster_t> new_chain_clusters;
                    for (cluster_t curr : chain_clusters) {
                        //Update the clusters to reach the current snarl
                        cluster_t new_c;
                        new_c.dist_right = curr.dist_right + offset;
                        new_c.dist_left = curr.dist_left;
                        if (new_c.seeds_start == nullptr) {
                            new_c.seeds_start = curr.seeds_start; 
                            new_c.seeds_end = curr.seeds_end; 
                        } else {
                            new_c.seeds_end->next = curr.seeds_start;
                            new_c.seeds_end = curr.seeds_end;
                        }
                        new_chain_clusters.push_back(new_c);
                    }
                    chain_clusters = move(new_chain_clusters);
                }
                last_snarl = boundary_node_id;
                last_len = node_length;
                dist_to_end = 0;
                start_length = node_length;
            }
#ifdef DEBUG
cerr << "  Cluster to add to chain: " << endl;

for (SnarlSeedClusterer::cluster_t c : curr_clusters) {
    cerr << "        cluster: " ;
    SnarlSeedClusterer::seed_index_list *s = c.seeds_start;
    while (s != nullptr){
        cerr << seeds[s->index] << "\t";
        s = s->next;
    }
    
    cerr << endl;
}
cerr << endl;
#endif

            //Combine the curr_clusters of both snarl and boundary node with
            //existing chain clusters. Chain clusters dist_right extends up to
            //the same node as dist_left of curr_clusters - both include the 
            //length of the shared node.
         
            vector<pair<hash_set<size_t>, hash_set<size_t>>> 
                                                    new_chain_cluster_indices;
            for (size_t i = 0 ; i < chain_clusters.size() ; i++) {
                hash_set<size_t> old_clusters;
                old_clusters.insert(i);
                hash_set<size_t> new_clusters;
                new_chain_cluster_indices.push_back(make_pair(
                                                  old_clusters, new_clusters));
            }  

            for (size_t j = 0; j < curr_clusters.size() ; j++) {
                /* For each of the clusters for the current snarl,
                 * find which chain clusters it belongs to and
                 * combine shared chain clusters */ 
                SnarlSeedClusterer::cluster_t snarl_cluster = curr_clusters[j];
                hash_set<size_t> curr_cluster_assignments;
                //Which chain cluster the snarl cluster is near

                for (size_t i = 0; i < new_chain_cluster_indices.size(); i++ ) {
                    //For each old chain cluster
                    bool in_this_cluster = false;
                    pair<hash_set<size_t>, hash_set<size_t>> 
                              old_chain_clusters = new_chain_cluster_indices[i];
                    for (size_t c : old_chain_clusters.first) {
                        if (!in_this_cluster) {
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

                //Combine chain clusters that are in the same cluster as this snarl cluster
                vector<pair<hash_set<size_t>, hash_set<size_t>>> new_chain_indices;
                pair<hash_set<size_t>, hash_set<size_t>> new_chain_cluster;
                new_chain_cluster.second.insert(j);
                for (size_t i = 0; i < new_chain_cluster_indices.size(); i++ ) {
                    /* Combine chain clusters that share this snarl cluster*/
                    pair<hash_set<size_t>, hash_set<size_t>> old_cluster
                                             = new_chain_cluster_indices[i];
                    if (curr_cluster_assignments.count(i) != 0) {
                        new_chain_cluster.first.insert( old_cluster.first.begin(),
                                                      old_cluster.first.end());
                        new_chain_cluster.second.insert( old_cluster.second.begin(),
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
                //Create a new cluster, curr_cluster, from this combined cluster

                SnarlSeedClusterer::cluster_t curr_cluster;

                for (size_t i : new_chain_cluster.first) {
                    //Combine the chain clusters
                    SnarlSeedClusterer::cluster_t old_chain_cluster = 
                                                              chain_clusters[i];

                    curr_cluster.dist_left =  curr_cluster.dist_left == -1 ?
                                         old_chain_cluster.dist_left : 
                                        min(old_chain_cluster.dist_left, 
                                                     curr_cluster.dist_left);
                    int64_t old_dist_right =  old_chain_cluster.dist_right +
                                     dist_to_end;
                    curr_cluster.dist_right = 
                               curr_cluster.dist_right == -1 ?  old_dist_right :
                                min(old_dist_right, curr_cluster.dist_right);

                    if (curr_cluster.seeds_start == nullptr) {
                        curr_cluster.seeds_start = old_chain_cluster.seeds_start; 
                        curr_cluster.seeds_end = old_chain_cluster.seeds_end; 
                    } else {
                        curr_cluster.seeds_end->next =old_chain_cluster.seeds_start;
                        curr_cluster.seeds_end = old_chain_cluster.seeds_end;
                    }
                          
                }


                for (size_t j : new_chain_cluster.second) {
                    //combine snarl clusters
                    SnarlSeedClusterer::cluster_t old_snarl_cluster = 
                                                        curr_clusters[j];
                    curr_cluster.dist_right =curr_cluster.dist_right == -1 ?
                                 old_snarl_cluster.dist_right : 
                              min(old_snarl_cluster.dist_right, 
                                                 curr_cluster.dist_right);
                    int64_t old_dist_left = old_snarl_cluster.dist_left
                               + chain_index.prefixSum[
                                      chain_index.snarlToIndex[last_snarl]] - 1;
                    curr_cluster.dist_left = 
                                curr_cluster.dist_left == -1 ? old_dist_left : 
                                    min(old_dist_left, curr_cluster.dist_left);

                    if (curr_cluster.seeds_start == nullptr) {
                        curr_cluster.seeds_start = old_snarl_cluster.seeds_start; 
                        curr_cluster.seeds_end = old_snarl_cluster.seeds_end; 
                    } else {
                        curr_cluster.seeds_end->next = 
                                                  old_snarl_cluster.seeds_start;
                        curr_cluster.seeds_end = old_snarl_cluster.seeds_end;
                    }
                }
                new_chain_clusters.push_back(curr_cluster);
            }
            chain_clusters = new_chain_clusters;
#ifdef DEBUG 
int64_t curr_id = curr_snarl == nullptr ? boundary_node_id : 
              curr_snarl->start().node_id();
cerr << "\t at snarl " << curr_id << ", clusters:" << endl;

for (SnarlSeedClusterer::cluster_t curr_cluster :chain_clusters) {
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
for (cluster_t c : chain_clusters) {
    cerr << "\tleft: " << c.dist_left << " right : " << c.dist_right << "\tseeds: ";
    SnarlSeedClusterer::seed_index_list *s = c.seeds_start;
    while (s != nullptr){
        cerr << seeds[s->index] << "\t";
        s = s->next;
    }
    cerr << endl;
}
#endif

        return chain_clusters; 
    };





    vector<SnarlSeedClusterer::cluster_t> SnarlSeedClusterer::get_clusters_snarl(
                         vector<pos_t>& seeds,
                         vector<SnarlSeedClusterer::seed_index_list>& seed_list,
                         hash_map< const Chain*, 
                                   vector<pair<const Snarl*, int64_t>> >& 
                                                        chains_to_snarl,
                         hash_map< const Snarl*, hash_set<pair<id_t, bool>> >& 
                                                           snarls_to_node,
                         hash_map< id_t, vector<size_t> >& node_to_seed,
                         size_t distance_limit, SnarlManager& snarl_manager,
                         DistanceIndex& dist_index, const Snarl* root,
                         bool include_boundaries) {
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
        for (auto v : snarl_index.visitToIndex) {

            id_t visit = v.first;
            const Snarl* snarl = snarl_manager.into_which_snarl(visit, false);
            if (snarl == nullptr) {
                snarl = snarl_manager.into_which_snarl(visit, true);
            }
            if (snarl == nullptr) {
                //If this is a node
                child_nodes.push_back(visit);
                
            } else if ( visit == snarl_index.snarlStart.first ||
                        visit == snarl_index.snarlEnd.first ) {
                //If this is a boundary node
                if ( include_boundaries ) {
                    child_nodes.push_back(visit);
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
            //Go through each child node of the netgraph and find clusters

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
                             seeds, seed_list, node_to_seed,
                             distance_limit, snarl_manager, dist_index, 
                             curr_node, false, node_len); 
                }

            } else if (i < child_snarls.size() + child_nodes.size()) {

                const Snarl* curr_snarl = child_snarls[i - child_nodes.size()];
                curr_node = curr_snarl->start().node_id();
                curr_rev = curr_snarl->start().backward();
                if (snarls_to_node.count(curr_snarl) != 0) {
                    curr_clusters = get_clusters_snarl(
                         seeds, seed_list, chains_to_snarl, snarls_to_node, 
                         node_to_seed,distance_limit, snarl_manager, dist_index,
                         curr_snarl, true);
                }
            } else {
                const Chain* curr_chain =  child_chains[i - child_nodes.size() -
                                             child_snarls.size()];
                Visit start = get_start_of(*curr_chain);
                curr_node = start.node_id();
                curr_rev = start.backward();
                if (chains_to_snarl.count(curr_chain) != 0) {
                    curr_clusters = get_clusters_chain(seeds, seed_list,
                             chains_to_snarl,
                             snarls_to_node, node_to_seed, distance_limit,
                             snarl_manager, dist_index, curr_chain);
                }
            }
            child_clusters.push_back(curr_clusters);

            for (size_t c = 0 ; c < curr_clusters.size() ; c++) {
                //for each cluster of child node i
                SnarlSeedClusterer::cluster_t child_cluster = curr_clusters[c];
 
                hash_set<size_t> assignments;
                for (size_t k = 0; k < combined_clusters_indices.size() ; k++){
                    //For each existing cluster, see which this child cluster belongs to
                    bool in_this_cluster = false;
                    vector<hash_set<size_t>> old_cluster_indices =
                                            combined_clusters_indices[k].first;
                    while (old_cluster_indices.size() <= i) {
                        //If this cluster doesn't have a set for this child node
                        //TODO: Maybe a better way of doing this??

                        hash_set<size_t> new_cluster ;
                        old_cluster_indices.push_back(new_cluster);
                        combined_clusters_indices[k].first = old_cluster_indices;
                    }
                    for (size_t j = 0 ; j < i ; j++){
                        //Go through child nodes up to but not including i
                        if (!in_this_cluster) {
                            
                            //Clusters of child j in this combined cluster
                            hash_set<size_t> other_clusters = 
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
                                     make_pair(other_node, other_rev))] - 1;
                            int64_t dist_l_r = snarl_index.distances[
                              snarl_index.index(
                                     make_pair(curr_node, !curr_rev), 
                                     make_pair(other_node, !other_rev))] - 1;
                            int64_t dist_r_l = snarl_index.distances[   
                                  snarl_index.index(
                                      make_pair(curr_node, curr_rev), 
                                      make_pair(other_node, other_rev))] - 1;
                            int64_t dist_r_r = snarl_index.distances[
                                    snarl_index.index(
                                      make_pair(curr_node, curr_rev), 
                                      make_pair(other_node, !other_rev))] - 1;
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
                //find distances to ends of snarl for this cluster 
                int64_t dist_s_l = snarl_index.distances[
                                      snarl_index.index(snarl_index.snarlStart,
                                          make_pair(curr_node, curr_rev))] - 1;
                int64_t dist_s_r = snarl_index.distances[
                                   snarl_index.index(snarl_index.snarlStart,
                                          make_pair(curr_node, !curr_rev))] - 1;
                int64_t dist_e_l = snarl_index.distances[
                                      snarl_index.index(
                                    make_pair(snarl_index.snarlEnd.first, 
                                                  !snarl_index.snarlEnd.second),
                                          make_pair(curr_node, curr_rev))] - 1;
                int64_t dist_e_r = snarl_index.distances[
                                     snarl_index.index(
                                    make_pair(snarl_index.snarlEnd.first, 
                                                  !snarl_index.snarlEnd.second),
                                         make_pair(curr_node, !curr_rev))] - 1;
                dist_s_l = dist_s_l == -1  || child_cluster.dist_left == -1 
                       ? -1 : start_length + dist_s_l + child_cluster.dist_left;
                dist_s_r = dist_s_r == -1 || child_cluster.dist_right == -1 
                      ? -1 : start_length + dist_s_r + child_cluster.dist_right;
                dist_e_l = dist_e_l == -1 || child_cluster.dist_left == -1
                       ? -1 : end_length + dist_e_l + child_cluster.dist_left;
                dist_e_r = dist_e_r == -1 || child_cluster.dist_right == -1
                      ? -1 : end_length + dist_e_r + child_cluster.dist_right;
  
                pair<int64_t, int64_t> new_dists (
                         DistanceIndex::minPos({dist_s_l, dist_s_r}), 
                         DistanceIndex::minPos({dist_e_l,dist_e_r}));
                if (snarl_index.snarlEnd.first == curr_node) {
                    new_dists.second = child_cluster.dist_right ;
                } 
                if (curr_node == snarl_index.snarlStart.first){
                    new_dists.second = child_cluster.dist_left ;
                }

                //Combine clusters that child cluster belongs to
                vector<pair<vector<hash_set<size_t>>, pair<int64_t, int64_t>>>
                                                     new_cluster_indices;
                vector<hash_set<size_t>> combined_cluster;
                for (size_t j = 0; j <= i ; j++) {
                    //Add a set of clusters for each node up to i
                    hash_set<size_t> set;
                    combined_cluster.push_back(set);
                }
                //Combine clusters
                for (size_t k = 0; k < combined_clusters_indices.size() ; k++){
               
                    if (assignments.count(k) > 0) {
                        //If this cluster belongs with the current cluster
                        vector<hash_set<size_t>> old_cluster_indices =
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
                            combined_cluster[j].insert(
                                old_cluster_indices[j].begin(),
                                old_cluster_indices[j].end());
                            }
                        }
                    } else {

                        new_cluster_indices.push_back(
                                                  combined_clusters_indices[k]);
                    } 
                    
                }
                combined_cluster[i].insert(c);
                new_cluster_indices.push_back(make_pair(combined_cluster, 
                                                         new_dists));
                combined_clusters_indices = move(new_cluster_indices);
                
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
                    if (cluster.seeds_start == nullptr) {
                            cluster.seeds_start = child_clusters[i][c].seeds_start; 
                            cluster.seeds_end = child_clusters[i][c].seeds_end; 
                        } else {
                            cluster.seeds_end->next = child_clusters[i][c].seeds_start;
                            cluster.seeds_end = child_clusters[i][c].seeds_end;
                        }
                }
            } 

            snarl_clusters.push_back(cluster); 
        }
#ifdef DEBUG 
cerr << "Found clusters on snarl " << root->start() << endl;
for (cluster_t c : snarl_clusters) {
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
