#include "seed_clusterer.hpp"

#include <algorithm>

//#define DEBUG_CLUSTER

namespace vg {

    SnarlSeedClusterer::SnarlSeedClusterer( MinimumDistanceIndex& dist_index) :
                                            dist_index(dist_index){
    };

    SnarlSeedClusterer::cluster_group_t SnarlSeedClusterer::cluster_seeds (vector<pos_t> seeds, int64_t read_distance_limit) const {
        vector<vector<pos_t>> all_seeds;
        all_seeds.push_back(seeds);
        tuple<vector<SnarlSeedClusterer::cluster_group_t>,SnarlSeedClusterer::cluster_group_t> all_clusters =
            cluster_seeds(all_seeds, read_distance_limit, 0);
        return std::get<0>(all_clusters)[0];
    };

    tuple<vector<SnarlSeedClusterer::cluster_group_t>,SnarlSeedClusterer::cluster_group_t> SnarlSeedClusterer::cluster_seeds (
                  vector<vector<pos_t>>& all_seeds, int64_t read_distance_limit,
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
        size_t seed_count = 0;
        for (auto& v : all_seeds) seed_count+= v.size();
        TreeState tree_state (&all_seeds, read_distance_limit, fragment_distance_limit, seed_count);


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
            //
            // tree_state knows all children of the snarls at this level

            if (depth != 0) {
                // Bring in the direct child nodes that come in at this level
                //  in the snarl tree.
                // They only ever occur below the root.
                tree_state.parent_snarl_to_nodes =
                                       move(snarl_to_nodes_by_level[depth - 1]);
            }

#ifdef DEBUG_CLUSTER
assert(tree_state.read_index_offsets[0] == 0);
for (size_t i = 1 ; i < tree_state.all_seeds->size() ; i++) {
    assert (tree_state.read_index_offsets[i] + tree_state.all_seeds->at(i).size() == tree_state.read_index_offsets[i+1]);
}
#endif
            //Cluster all the snarls at this depth
            //Also records which snarls are in chains and the parents of these
            //snarls in tree_state.parent_snarl_to_node
            cluster_snarls(tree_state, depth);

            //And cluster all the chains, record the parents of these chains
            cluster_chains(tree_state, depth);

            // Swap buffer over for the next level
            tree_state.snarl_to_nodes = move(tree_state.parent_snarl_to_nodes);
            tree_state.chain_to_snarls.clear();
        }

#ifdef DEBUG_CLUSTER
        cerr << "Found read clusters : " << endl;
        for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
            cerr << "\t read num " << read_num << ": " ;
            for (auto group : tree_state.read_union_find[read_num].all_groups()){
                cerr << "\t\t";
                for (size_t c : group) {
                   cerr << tree_state.all_seeds->at(read_num)[c] << " ";
                }
                cerr << endl;
            }
            cerr << endl;
        }
        vector<pos_t> ordered_seeds;
        for (size_t i = 0 ; i < tree_state.all_seeds->size() ; i++) {
            for (auto x : tree_state.all_seeds->at(i)) {
                ordered_seeds.push_back(x);
            }
        }
        cerr << "Found fragment clusters : " << endl;
        for (auto group : tree_state.fragment_union_find.all_groups()){
            cerr << "\t";
            for (size_t c : group) {
               cerr << ordered_seeds[c] << " ";
            }
            cerr << endl;
        }
#endif
        vector<SnarlSeedClusterer::cluster_group_t> read_clusters;
        for (auto& uf : tree_state.read_union_find) {
            read_clusters.emplace_back(uf.all_groups());
        }
        return make_tuple(std::move(read_clusters),
                          tree_state.fragment_union_find.all_groups());

    };

    //Helper function to get the minimum value that is not -1
    int64_t min_not_minus_one(int64_t n1, int64_t n2) {
        return static_cast<int64_t>(std::min(static_cast<uint64_t>(n1), static_cast<uint64_t>(n2)));
    }

    void SnarlSeedClusterer::get_nodes( TreeState& tree_state,
              vector<hash_map<size_t,vector<pair<NetgraphNode, NodeClusters>>>>&
                                                               snarl_to_nodes_by_level) const {

        // Assign each seed to a node.
        for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++){ 
            vector<pos_t>& seeds = tree_state.all_seeds->at(read_num);
            for (size_t i = 0; i < seeds.size(); i++) {
                id_t id = get_id(seeds.at(i));
                tree_state.node_to_seeds[read_num].emplace_back(id, i);
                //For each seed, assign it to a node and the node to a snarl
            }
            std::sort(tree_state.node_to_seeds[read_num].begin(), tree_state.node_to_seeds[read_num].end());
        }

        // Assign each node to a snarl.
        hash_set<id_t> seen_nodes;
        for (auto& read_node :tree_state.node_to_seeds) {
            for (auto& mapping : read_node) {
                if (seen_nodes.count(mapping.first) < 1) {
                    seen_nodes.insert( mapping.first);
                    size_t snarl_i = dist_index.getPrimaryAssignment(mapping.first);
                    size_t depth = dist_index.snarl_indexes[snarl_i].depth;
                    snarl_to_nodes_by_level[depth][snarl_i].emplace_back(
                             NetgraphNode(mapping.first, NODE), NodeClusters(tree_state.all_seeds->size()));
                }
            }
        }
    }


    void SnarlSeedClusterer::cluster_snarls(TreeState& tree_state, size_t depth) const {

        for (auto& kv : tree_state.snarl_to_nodes){
            //Go through each of the snarls at this level, cluster them,
            //and find which chains they belong to, if any
            //key is the index of the snarl and value is a vector of pair of
            // NetgraphNode, NodeClusters

            size_t snarl_i = kv.first;
            MinimumDistanceIndex::SnarlIndex& snarl_index =
                                     dist_index.snarl_indexes[snarl_i];

#ifdef DEBUG_CLUSTER
            cerr << "At depth " << depth << " snarl number " << snarl_i
                << " headed by " << snarl_index.id_in_parent
                << " with children " << endl;
            for (auto it2 : kv.second) {
                cerr << "\t" << typeToString(it2.first.node_type)
                     << " number " << it2.first.node_id << endl;
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
                            cluster_one_snarl(tree_state, snarl_i)));

#ifdef DEBUG_CLUSTER
                cerr << "Recording snarl number " << snarl_i << " headed by "
                      << snarl_index.id_in_parent  << " as a child of chain number "
                      << chain_assignment << " headed by " << snarl_index.parent_id << endl;
#endif

            } else {
                //If this snarl is not in a chain, add it as a child of the
                //parent snarl for the next level

                if (depth != 0 && snarl_index.parent_id != 0){
                    //If this has a parent, record it
#ifdef DEBUG_CLUSTER
                    assert(snarl_index.parent_id >= dist_index.min_node_id);
                    assert(snarl_index.parent_id <= dist_index.max_node_id);
#endif
                    size_t parent_snarl_i =
                           dist_index.getPrimaryAssignment(
                                                snarl_index.parent_id);

                    tree_state.parent_snarl_to_nodes[parent_snarl_i].emplace_back(
                            NetgraphNode (snarl_i, SNARL),
                            cluster_one_snarl(tree_state, snarl_i));

#ifdef DEBUG_CLUSTER
                    cerr << "Recording snarl number " << snarl_i
                        << " headed by " << snarl_index.id_in_parent
                        << " as a child of snarl number " << parent_snarl_i
                        << " headed by " << snarl_index.parent_id << endl;
                    cerr << "Snarl number " << parent_snarl_i << " has "
                        << tree_state.parent_snarl_to_nodes[parent_snarl_i].size()
                        << " children now" << endl;
#endif
                } else {
                    //If it doesn't have a parent, don't need to keep track of
                    //its contents
                    cluster_one_snarl(tree_state, snarl_i);
                }
            }
        }
    }

    void SnarlSeedClusterer::cluster_chains(TreeState& tree_state, size_t depth) const {
        for (auto& kv : tree_state.chain_to_snarls) {
            //For each chain at this level that has relevant child snarls in it,
            //find the clusters.

            // Get the chain's number
            size_t chain_i = kv.first;

#ifdef DEBUG_CLUSTER
            cerr << "At depth " << depth << " chain number " << chain_i
                 << " with children " << endl;
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
#ifdef DEBUG_CLUSTER
                assert(parent_id >= dist_index.min_node_id);
                assert(parent_id <= dist_index.max_node_id);
#endif
                // Map it to the snarl number that should be represented by it
                // (and thus also contain the chain)
                size_t parent_snarl_i =dist_index.getPrimaryAssignment(parent_id);

                // Register clusters as relevant for that parent snarl.

                tree_state.parent_snarl_to_nodes[parent_snarl_i].emplace_back(
                      NetgraphNode (chain_i, CHAIN), std::move(chain_clusters));
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
        /*Find clusters of seeds in this node.
         * Result contains hash_set of the union find group IDs of the new clusters,
         * and the shortest distance from any seed to the left and right sides
         * of the node*/

        //indices of union find group ids of clusters in this node
        NodeClusters node_clusters(tree_state.all_seeds->size());

        if (tree_state.read_distance_limit > node_length) {
            //If the limit is greater than the node length, then all the
            //seeds on this node must be in the same cluster

            size_t fragment_group_id = -1;
            for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
                auto seed_range_start = std::lower_bound(
                    tree_state.node_to_seeds[read_num].begin(),
                    tree_state.node_to_seeds[read_num].end(),
                    std::pair<id_t, size_t>(node_id, 0));
                if (seed_range_start != tree_state.node_to_seeds[read_num].end() && seed_range_start->first == node_id) {

                    size_t group_id = seed_range_start->second;

                    for (auto iter = seed_range_start; iter != tree_state.node_to_seeds[read_num].end() && iter->first == node_id; ++iter) {
                        //For each seed on this node, add it to the cluster
                        //And find the shortest distance from any seed to both
                        //ends of the node

                        pos_t seed = tree_state.all_seeds->at(read_num)[iter->second];
                        int64_t dist_left = is_rev(seed) ? node_length- get_offset(seed)
                                                         : get_offset(seed) + 1;
                        int64_t dist_right = is_rev(seed) ? get_offset(seed) + 1
                                                       : node_length - get_offset(seed);

                        node_clusters.read_best_left[read_num] = min_not_minus_one(dist_left,
                                                              node_clusters.read_best_left[read_num]);
                        node_clusters.read_best_right[read_num] = min_not_minus_one(dist_right,
                                                              node_clusters.read_best_right[read_num]);
                        node_clusters.fragment_best_left = min_not_minus_one(dist_left,
                                                              node_clusters.fragment_best_left);
                        node_clusters.fragment_best_right = min_not_minus_one(dist_right,
                                                              node_clusters.fragment_best_right);

                        tree_state.read_union_find[read_num].union_groups(group_id, iter->second);
                        if (tree_state.fragment_distance_limit != 0 ) {
                            if (fragment_group_id == -1 ) fragment_group_id = seed_range_start->second + tree_state.read_index_offsets[read_num];
                            tree_state.fragment_union_find.union_groups(fragment_group_id, iter->second + tree_state.read_index_offsets[read_num]);
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
            cerr << "Found single cluster on node " << node_id << " with fragment dists " << node_clusters.fragment_best_left << " " << node_clusters.fragment_best_right << endl;

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
                        for (size_t x = 0 ; x < tree_state.all_seeds->at(c.first).size() ; x++) {
                            if (tree_state.read_union_find[c.first].find_group(x) == c.second) {
                                cerr << tree_state.all_seeds->at(c.first)[x] << " ";
                                has_seeds = true;
                            }
                        }
                        assert(dists.first == -1 || dists.first >= node_clusters.read_best_left[read_num]);
                        assert(dists.second == -1 || dists.second >= node_clusters.read_best_right[read_num]);
                        assert(dists.first == -1 || dists.first >= node_clusters.fragment_best_left);
                        assert(dists.second == -1 || dists.second >= node_clusters.fragment_best_right);
                        if (dists.first == node_clusters.fragment_best_left) {got_left = true;}
                        if (dists.second == node_clusters.fragment_best_right) {got_right = true;}
                        if (dists.first == node_clusters.read_best_left[read_num]) {got_read_left = true;}
                        if (dists.second == node_clusters.read_best_right[read_num]) {got_read_right = true;}
                        cerr << endl;
                        assert(has_seeds);
                    }
                }
                assert(got_read_left || node_clusters.read_best_left[read_num] == -1);
                assert(got_read_right || node_clusters.read_best_right[read_num] == -1);
            }
            assert(got_left);
            assert(got_right);

#endif
            return node_clusters;
        }

        vector<tuple<size_t,size_t, int64_t>> seed_offsets;
        for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
            //<index of read, index of seed, offset of seed> for all seeds
                auto seed_range_start = std::lower_bound(
                    tree_state.node_to_seeds[read_num].begin(),
                    tree_state.node_to_seeds[read_num].end(),
                    std::pair<id_t, size_t>(node_id, 0));
            if (seed_range_start != tree_state.node_to_seeds[read_num].end()) {
                for (auto iter = seed_range_start; iter != tree_state.node_to_seeds[read_num].end() && iter->first == node_id; ++iter) {
                    //For each seed, find its offset
                    pos_t seed = tree_state.all_seeds->at(read_num)[iter->second];
                    int64_t offset = is_rev(seed) ? node_length - get_offset(seed)
                                                    : get_offset(seed) + 1;

                    node_clusters.fragment_best_left = min_not_minus_one(offset, node_clusters.fragment_best_left);
                    node_clusters.fragment_best_right = min_not_minus_one(node_length-offset+1, node_clusters.fragment_best_right);
                    node_clusters.read_best_left[read_num] = min_not_minus_one(offset, node_clusters.read_best_left[read_num]);
                    node_clusters.read_best_right[read_num] = min_not_minus_one(node_length-offset+1, node_clusters.read_best_right[read_num]);

                    seed_offsets.emplace_back(read_num, iter->second, offset);

                }
            }
        }
        //Sort seeds by their position in the node
        std::sort(seed_offsets.begin(), seed_offsets.end(),
                     [&](const auto a, const auto b) -> bool {
                          return  std::get<2>(a) < std::get<2>(b);
                      } );

        vector<int64_t> read_last_offset (tree_state.all_seeds->size(), -1);
        int64_t fragment_last_offset = -1;
        size_t fragment_last_cluster = -1;
        vector<size_t> read_last_cluster (tree_state.all_seeds->size(), -1);

        for ( tuple<size_t, size_t, int64_t> s : seed_offsets) {
            //For each seed, in order of position in the node,
            //see if it belongs to a new read/fragment cluster - if it is
            //close enough to the previous seed
            size_t read_num = std::get<0>(s);

            if (read_last_offset[read_num] != -1 &&
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
                if (read_last_cluster[read_num] != -1) {
                    node_clusters.read_cluster_heads.emplace(read_num, read_last_cluster[read_num]);
                }
                read_last_cluster[read_num] = std::get<1>(s);
                read_last_offset[read_num] = std::get<2>(s);
                tree_state.read_cluster_dists[read_num][read_last_cluster[read_num]] = 
                        make_pair(read_last_offset[read_num], node_length - read_last_offset[read_num] + 1);
                if (tree_state.fragment_distance_limit != 0) {
                    if (fragment_last_offset != -1 &&
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
            if (read_last_cluster[i] != -1) {
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
                    for (size_t x = 0 ; x < tree_state.all_seeds->at(c.first).size() ; x++) {
                        if (tree_state.read_union_find[c.first].find_group(x) == c.second) {
                            cerr << tree_state.all_seeds->at(c.first)[x] << " ";
                            has_seeds = true;
                        }
                    }
                    assert(dists.first == -1 || dists.first >= node_clusters.read_best_left[read_num]);
                    assert(dists.second == -1 || dists.second >= node_clusters.read_best_right[read_num]);
                    assert(dists.first == -1 || dists.first >= node_clusters.fragment_best_left);
                    assert(dists.second == -1 || dists.second >= node_clusters.fragment_best_right);
                    if (dists.first == node_clusters.fragment_best_left) {got_left = true;}
                    if (dists.second == node_clusters.fragment_best_right) {got_right = true;}
                    if (dists.first == node_clusters.read_best_left[read_num]) {got_read_left = true;}
                    if (dists.second == node_clusters.read_best_right[read_num]) {got_read_right = true;}
                    cerr << endl;
                    assert(has_seeds);
                }
            }
            assert(got_read_left || node_clusters.read_best_left[read_num] == -1);
            assert(got_read_right || node_clusters.read_best_right[read_num] == -1);
        }
        assert(got_left);
        assert(got_right);
        for (pair<size_t, size_t> group_id : node_clusters.read_cluster_heads) {
            assert (group_id.second == tree_state.read_union_find[group_id.first].find_group(group_id.second));
        }
#endif
        
        return node_clusters;

    };



    SnarlSeedClusterer::NodeClusters SnarlSeedClusterer::cluster_one_chain(
                               TreeState& tree_state, size_t chain_index_i) const {
        /*
         * Find all the clusters in the given chain
         */

        std::map<size_t, pair<size_t, NodeClusters>>& snarls_in_chain =
                                      tree_state.chain_to_snarls[chain_index_i];

        MinimumDistanceIndex::ChainIndex& chain_index = dist_index.chain_indexes[
                                                            chain_index_i];
#ifdef DEBUG_CLUSTER
        cerr << "Finding clusters on chain number " << chain_index_i
             << " headed by node " << chain_index.id_in_parent << endl;
#endif

        auto combine_snarl_clusters = [&] (size_t& new_group,
                        size_t& combined_group, size_t& fragment_combined_group,
                        vector<pair<size_t,size_t>>& to_erase, int64_t fragment_dist,int64_t read_dist,
                        pair<int64_t, int64_t>& dists, size_t read_num){
            //Helper function to combine clusters of the same snarl
            //Used when two clusters in the same snarl can be combined by
            //looping in the chain

            if (read_dist != -1 && read_dist <= tree_state.read_distance_limit) {
                if (combined_group == -1) {
                    combined_group = new_group;
                } else {
                    //Union the two groups
                    combined_group = tree_state.read_union_find[read_num].find_group(combined_group);
                    tree_state.read_union_find[read_num].union_groups(combined_group, new_group);
                    //Find the new distances of the combined groups
                    pair<int64_t, int64_t>& old_dists =
                                           tree_state.read_cluster_dists[read_num][combined_group];
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
                          min_not_minus_one(old_dists.first, dists.first),
                          min_not_minus_one(old_dists.second, dists.second));
                    tree_state.read_cluster_dists[read_num][new_group] = dists;
                    tree_state.read_cluster_dists[read_num][combined_group] = dists;
#ifdef DEBUG_CLUSTER
                    cerr << " New dists for read num " << read_num << ": "
                         << tree_state.read_cluster_dists[read_num][combined_group].first << " "
                         << tree_state.read_cluster_dists[read_num][combined_group].second << endl;
#endif
                }

                if (tree_state.fragment_distance_limit != 0) {
                    if (fragment_combined_group != -1) {
                    //If we're keeping track of fragment clusters, union this
                        tree_state.fragment_union_find.union_groups(fragment_combined_group,
                                                            new_group + tree_state.read_index_offsets[read_num]);
                    }
                    fragment_combined_group = tree_state.fragment_union_find.find_group(
                                                         new_group + tree_state.read_index_offsets[read_num]);
                }
            } else if (tree_state.fragment_distance_limit != 0 &&
                        fragment_dist <= tree_state.fragment_distance_limit) {
                //If these aren't in the same read cluster but are in
                //the same fragment cluster
                if (fragment_combined_group != -1) {
                    tree_state.fragment_union_find.union_groups(fragment_combined_group,
                                                        new_group + tree_state.read_index_offsets[read_num]);
                }
                fragment_combined_group = tree_state.fragment_union_find.find_group(
                                                     new_group + tree_state.read_index_offsets[read_num]);
            }
            return;
        };
        //The clusters of the chain that are built from the snarl clusters
        //This will get updated as we traverse through the snarls
        NodeClusters chain_clusters(tree_state.all_seeds->size());

        //The rank of the node at which the chain clusters reach
        // (the last snarl that was traversed)
        size_t last_rank = 0;
        int64_t last_len = 0;
        id_t start_node;
        id_t end_node;

        for (auto& kv: snarls_in_chain) {
            /* For each child snarl in the chain, get the clusters from the
             * tree_state and progressively build up clusters spanning up
             * to that snarl
             * Snarls are in the order that they are traversed in the chain
             * snarls_in_chain is a map from rank of a snarl to the snarl index
             *  in dist_index.snarl_indexes and NodeClusters for that snarl
             */


            //rank of the boundary node of the snarl that occurs first in
            //the chain
            size_t start_rank = kv.first;
            size_t curr_snarl_i = kv.second.first;

            //The clusters of the current snarl
            NodeClusters& snarl_clusters = kv.second.second;

            MinimumDistanceIndex::SnarlIndex& snarl_index =
                                         dist_index.snarl_indexes[curr_snarl_i];

            //Get the lengths of the start and end nodes of the snarl, relative
            //to the order of the chain
            int64_t start_length = snarl_index.rev_in_parent
                         ? snarl_index.nodeLength(snarl_index.num_nodes * 2 - 1)
                         : snarl_index.nodeLength(0);
            int64_t end_length = snarl_index.rev_in_parent ? snarl_index.nodeLength(0)
                        : snarl_index.nodeLength(snarl_index.num_nodes * 2 - 1);

            //Distance from right end of chain clusters that have been made so far
            //to the end of current snarl cluster
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

                for (pair<size_t, size_t> c : chain_clusters.read_cluster_heads) {
                    tree_state.read_cluster_dists[c.first][c.second].second =
                            tree_state.read_cluster_dists[c.first][c.second].second == -1
                            ? -1 : tree_state.read_cluster_dists[c.first][c.second].second + offset;
                }
                chain_clusters.fragment_best_right = chain_clusters.fragment_best_right == -1 ? -1
                                            : chain_clusters.fragment_best_right + offset;
                for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
                    chain_clusters.read_best_right[read_num] = chain_clusters.read_best_right[read_num] == -1 ? -1
                                            : chain_clusters.read_best_right[read_num] + offset;
                }
            }

            last_rank = start_rank + 1;
            last_len = end_length;


            //Distance from the start of chain to the start of the current snarl
            int64_t add_dist_left = start_rank == 0 ? 0 :
                                    chain_index.prefix_sum[start_rank] - 1;



            //Combine snarl clusters that can be reached by looping
            int64_t loop_dist_end = chain_index.loop_fd[start_rank + 1] - 1 ;
            int64_t loop_dist_start = chain_index.loop_rev[start_rank] - 1;


#ifdef DEBUG_CLUSTER
            cerr << "Looking at snarl rank " << start_rank << " representing " << snarl_index.id_in_parent << endl;
            cerr << "  Snarl fragment distance limits: " << snarl_clusters.fragment_best_left
                 << " " << snarl_clusters.fragment_best_right << endl;
            cerr << "  Snarl clusters to add: " << endl;
            for (pair<size_t, size_t> c : snarl_clusters.read_cluster_heads) {
                pair<int64_t, int64_t> dists = tree_state.read_cluster_dists[c.first][c.second];
                cerr << "\tread " << c.first << ",cluster " << c.second << " left: " << dists.first << " right : " << dists.second
                     << endl;
                cerr << "\t\t";
                bool has_seeds = false;
                for (size_t x = 0 ; x < tree_state.all_seeds->at(c.first).size() ; x++) {
                    if (tree_state.read_union_find[c.first].find_group(x) == c.second) {
                        cerr << tree_state.all_seeds->at(c.first)[x] << " ";
                        has_seeds = true;
                    }
                }
                assert(has_seeds);
                cerr << endl;
            }
            cerr << endl;

            cerr << "  Clusters on chain: " << endl;

            cerr << "  best left: " << chain_clusters.fragment_best_left << " best right: "
                  << chain_clusters.fragment_best_right << endl;
            for (pair<size_t, size_t> c : chain_clusters.read_cluster_heads) {
                pair<int64_t, int64_t> dists = tree_state.read_cluster_dists[c.first][c.second];
                cerr << "\tleft: " << dists.first << " right : " << dists.second
                     << endl;
                cerr << "\t\t";
                for (size_t x = 0 ; x < tree_state.all_seeds->at(c.first).size() ; x++) {
                    if (tree_state.read_union_find[c.first].find_group(x) == c.second) {
                        cerr << tree_state.all_seeds->at(c.first)[x] << " ";
                    }
                }
                cerr << endl;
            }
            cerr << endl;
#endif


            //Need to remember this to check if snarl clusters overlap the old
            //best distance
            int64_t fragment_chain_right = chain_clusters.fragment_best_right;
            vector<int64_t> read_chain_right = std::move(chain_clusters.read_best_right);

            vector<pair<size_t,size_t>> to_add;//new cluster group ids from snarl clusters
            vector<pair<size_t,size_t>> to_erase; //old cluster group ids
            //New cluster- there will be at most one new cluster to add
            vector< size_t> combined_cluster (tree_state.all_seeds->size(), -1);
            size_t fragment_combined_cluster = -1;
            vector<int64_t> combined_left (tree_state.all_seeds->size(), -1); 
            vector<int64_t> combined_right (tree_state.all_seeds->size(), -1); 

            //Combined snarl clusters by taking chain loop left/right
            vector<size_t> snarl_cluster_left (tree_state.all_seeds->size(),-1);
            vector<size_t> snarl_cluster_right (tree_state.all_seeds->size(), -1);
            size_t fragment_snarl_cluster_left = -1;
            size_t fragment_snarl_cluster_right = -1;

            chain_clusters.fragment_best_right = -1;
            chain_clusters.read_best_right.assign(tree_state.all_seeds->size(), -1);
            for (pair<size_t, size_t> cluster_head : snarl_clusters.read_cluster_heads) {
                // For each of the clusters for the current snarl,
                // first check if it can be combined with any other
                // snarl clusters by taking loops in the chain,
                // then, find if it belongs to the new combined cluster
                // that includes chain clusters
                size_t read_num = cluster_head.first;

                pair<int64_t, int64_t> snarl_dists =
                                        std::move(tree_state.read_cluster_dists[read_num][cluster_head.second]);

                if (loop_dist_start != -1) {
                    //If there is a loop going out and back into the start of
                    //the snarl, might combine this cluster with other snarl
                    //clusters

                    //The distance to the right side of the snarl
                    // that is found by taking the leftmost seed and
                    // looping through the chain to the left
                    int64_t new_right = snarl_dists.first == -1 || loop_dist_start == -1
                                        ? -1
                                        : snarl_dists.first + loop_dist_start + snarl_length - start_length;
                    snarl_dists.second = min_not_minus_one(new_right, snarl_dists.second);
                    snarl_clusters.fragment_best_right =
                               min_not_minus_one(snarl_clusters.fragment_best_right, new_right);
                    snarl_clusters.read_best_right[read_num] =
                               min_not_minus_one(snarl_clusters.read_best_right[read_num], new_right);
#ifdef DEBUG_CLUSTER
cerr << "  (Possibly) updating looping distance to right of snarl cluster " << read_num <<":" << cluster_head.second << ": "
     << new_right << " -> " << snarl_dists.second <<  endl;
#endif


                    if (snarl_clusters.fragment_best_left!= -1 && snarl_dists.first != -1 ) {
                        //If this cluster can be combined with another cluster
                        //from the left

#ifdef DEBUG_CLUSTER
cerr << "  Combining this cluster from the left " ;
#endif
                        int64_t read_dist =  snarl_clusters.read_best_left[read_num] == -1 ? -1 :  
                                     snarl_clusters.read_best_left[read_num] + snarl_dists.first + loop_dist_start - start_length - 1;
                        combine_snarl_clusters(cluster_head.second, snarl_cluster_left[read_num], fragment_snarl_cluster_left,
                                     to_erase, snarl_clusters.fragment_best_left + snarl_dists.first + loop_dist_start - start_length - 1,
                                     read_dist, snarl_dists, read_num);
                    }

                }

                if (loop_dist_end != -1) {
                    //If there is a loop to the right
                    int64_t new_left = snarl_dists.second == -1 || loop_dist_end == -1
                        ? -1
                          : snarl_dists.second + loop_dist_end + snarl_length - end_length;
                    if (snarl_dists.first == -1 || (new_left != -1 & new_left < snarl_dists.first)){
                        //If this is an improvement, update distances
                        snarl_dists.first = new_left;
                        snarl_clusters.read_best_left[read_num] = 
                                min_not_minus_one(new_left, snarl_clusters.read_best_left[read_num]);
                        snarl_clusters.fragment_best_left = min_not_minus_one(new_left, snarl_clusters.fragment_best_left);

#ifdef DEBUG_CLUSTER
cerr << "Updating looping distance to left of snarl cluster " << read_num << ":" << cluster_head.second << ": "
     << new_left << endl;
#endif
                    }

                    if (snarl_clusters.fragment_best_right != -1 && snarl_dists.second != -1 ) {
                        //If this cluster can be combined with another cluster
                        //from the right

#ifdef DEBUG_CLUSTER
cerr << "  Combining this cluster from the right" << endl;
#endif
                        int64_t read_dist = snarl_clusters.read_best_right[read_num] == -1 ? -1 :
                            snarl_clusters.read_best_right[read_num] + snarl_dists.second  + loop_dist_end - end_length - 1;
                        combine_snarl_clusters(cluster_head.second, snarl_cluster_right[read_num],
                             fragment_snarl_cluster_right, to_erase,
                            snarl_clusters.fragment_best_right + snarl_dists.second + loop_dist_end - end_length - 1,
                            read_dist, snarl_dists, read_num);
                    }
                }

                //Now check if this snarl cluster can be combined with any
                //existing chain clusters
                if (read_chain_right[read_num] != -1 && snarl_dists.first != -1 &&
                    snarl_dists.first + read_chain_right[read_num] - start_length-1
                                                <= tree_state.read_distance_limit) {
                    //If this snarl cluster's leftmost seed is close enough to
                    //the rightmost seed in the chain (up to this point), then
                    //this snarl cluster is in the combined cluster

                    if (combined_cluster[read_num] == -1) {
                        combined_cluster[read_num] = cluster_head.second;
                        combined_left[read_num] = snarl_dists.first == -1 ? -1 :
                                            snarl_dists.first + add_dist_left;
                        combined_right[read_num] = snarl_dists.second;
                    } else {
                        //Cluster
                        tree_state.read_union_find[read_num].union_groups(combined_cluster[read_num], cluster_head.second);
                        size_t new_group  = tree_state.read_union_find[read_num].find_group(cluster_head.second);

                        if (new_group == cluster_head.second) {
                            to_erase.emplace_back(read_num,combined_cluster[read_num]);
                        } else {
                            to_erase.push_back(cluster_head);
                        }

                        combined_cluster[read_num] = new_group;
                        combined_left[read_num] = min_not_minus_one(combined_left[read_num],
                                    snarl_dists.first == -1 ? -1 : snarl_dists.first + add_dist_left);
                        combined_right[read_num] = min_not_minus_one(combined_right[read_num],snarl_dists.second);
                    }
                    if (tree_state.fragment_distance_limit != 0) {
                        if (fragment_combined_cluster != -1) {
                            //Also cluster by fragment
                            tree_state.fragment_union_find.union_groups(fragment_combined_cluster, 
                                                    cluster_head.second+tree_state.read_index_offsets[read_num]);
                        }
                        fragment_combined_cluster = tree_state.fragment_union_find.find_group(cluster_head.second+tree_state.read_index_offsets[read_num]);
                    }
                } else {
                    //If the snarl cluster does not get combined with any of
                    //the existing chain clusters, then it becomes a new chain cluster
                    if (tree_state.fragment_distance_limit != 0 && fragment_chain_right != -1 && snarl_dists.first != -1 &&
                           snarl_dists.first+fragment_chain_right-start_length-1 <= tree_state.fragment_distance_limit) {
                        //Cluster in the same fragment but not the same read
                        if (fragment_combined_cluster != -1) {
                            //Also cluster by fragment
                            tree_state.fragment_union_find.union_groups(fragment_combined_cluster, 
                                                    cluster_head.second+tree_state.read_index_offsets[read_num]);
                        }
                        fragment_combined_cluster = tree_state.fragment_union_find.find_group(cluster_head.second+tree_state.read_index_offsets[read_num]);
                    }
                    to_add.push_back(cluster_head);
                    //Update its distances to the correct nodes in the chain
                    pair<int64_t, int64_t> d = make_pair(snarl_dists.first == -1 ? -1 : snarl_dists.first + add_dist_left,
                                                snarl_dists.second);
                    chain_clusters.fragment_best_left = min_not_minus_one(chain_clusters.fragment_best_left,d.first);
                    chain_clusters.fragment_best_right = min_not_minus_one(chain_clusters.fragment_best_right,d.second);
                    chain_clusters.read_best_left[read_num] = min_not_minus_one(chain_clusters.read_best_left[read_num],
                                                            d.first);
                    chain_clusters.read_best_right[read_num] = min_not_minus_one(chain_clusters.read_best_right[read_num],
                                                             d.second);

                    tree_state.read_cluster_dists[read_num][cluster_head.second] = std::move(d);
                }
            }

            //Next, go through each of the clusters of the chain and decide
            //if they get combined with snarl clusters
            for (pair<size_t, size_t> cluster_head : chain_clusters.read_cluster_heads) {
                //For each old chain cluster
                size_t read_num = cluster_head.first;
                pair<int64_t, int64_t>& chain_dists = tree_state.read_cluster_dists[read_num][cluster_head.second];

                if (snarl_clusters.read_best_left[read_num] != -1 && chain_dists.second != -1
                     && chain_dists.second + snarl_clusters.read_best_left[read_num]
                                - start_length-1 <= tree_state.read_distance_limit){
                    //If this chain cluster's rightmost seed is close enough
                    //to the leftmost seed of any cluster in the snarl, then
                    //this chain cluster is in the combined cluster

                    if (combined_cluster[read_num] == -1) {
                        //New chain cluster
                        combined_cluster[read_num] = cluster_head.second;
                        combined_left[read_num] = chain_dists.first;
                        combined_right[read_num] = chain_dists.second + dist_to_end;
                    } else {
                        //Combine
                        tree_state.read_union_find[read_num].union_groups(combined_cluster[read_num], cluster_head.second);
                        size_t new_group = tree_state.read_union_find[read_num].find_group(cluster_head.second);
                        if (new_group == cluster_head.second) {
                            to_erase.emplace_back(read_num,combined_cluster[read_num]);
                        } else {
                            to_erase.push_back(cluster_head);
                        }
                        combined_cluster[read_num] = new_group;
                        combined_left[read_num] = min_not_minus_one(combined_left[read_num], chain_dists.first);
                        combined_right[read_num] = min_not_minus_one(combined_right[read_num], chain_dists.second + dist_to_end);
                    }
                    if (tree_state.fragment_distance_limit != 0) {
                        if (fragment_combined_cluster != -1) {
                            tree_state.fragment_union_find.union_groups(fragment_combined_cluster, cluster_head.second+tree_state.read_index_offsets[read_num]);
                        }
                        fragment_combined_cluster = tree_state.fragment_union_find.find_group(cluster_head.second+tree_state.read_index_offsets[read_num]);
                    }
                } else {
                    //If this chain cluster is on its own, extend its right
                    //distance to the end of the current snarl

                    if (tree_state.fragment_distance_limit != 0 &&
                        snarl_clusters.fragment_best_left != -1 && chain_dists.second != -1
                        && chain_dists.second + snarl_clusters.fragment_best_left
                                - start_length-1 <= tree_state.fragment_distance_limit) {
                        //If this is a new read cluster but the same fragment cluster
                        if (fragment_combined_cluster != -1) {
                            tree_state.fragment_union_find.union_groups(fragment_combined_cluster, cluster_head.second+tree_state.read_index_offsets[read_num]);
                        }
                        fragment_combined_cluster = tree_state.fragment_union_find.find_group(cluster_head.second+tree_state.read_index_offsets[read_num]);
                    }
                    chain_dists.second += dist_to_end;
                    if ((tree_state.fragment_distance_limit == 0 &&
                         chain_dists.first - 2 >= tree_state.read_distance_limit &&
                         chain_dists.second - end_length-2 >= tree_state.read_distance_limit) ||
                        (tree_state.fragment_distance_limit != 0 &&
                         chain_dists.first - 2 >= tree_state.fragment_distance_limit &&
                         chain_dists.second - end_length-2 >= tree_state.fragment_distance_limit)) {
                        //If the distance from the seeds in this cluster to
                        //either end of the chain is greater than the distance
                        //limit, then it cannot cluster with anything else
                        //so we can stop keeping track of it
#ifdef DEBUG_CLUSTER
                        cerr << "Removing cluster " << cluster_head.first << ":" << cluster_head.second << endl;
#endif
                        to_erase.push_back(cluster_head);
                    } else {
                        chain_clusters.fragment_best_left =  min_not_minus_one(chain_clusters.fragment_best_left, chain_dists.first);
                        chain_clusters.fragment_best_right = min_not_minus_one(chain_clusters.fragment_best_right, chain_dists.second);
                        chain_clusters.read_best_left[read_num] = min_not_minus_one(chain_clusters.read_best_left[read_num], chain_dists.first);
                        chain_clusters.read_best_right[read_num] = min_not_minus_one(chain_clusters.read_best_right[read_num], chain_dists.second);
                    }
                }
            }
            //Update the chain cluster heads
            for (auto c : to_add) {
                chain_clusters.read_cluster_heads.insert(c);
            }
            for (auto c : to_erase) {
                chain_clusters.read_cluster_heads.erase(c);
            }
            for (size_t read_num = 0 ; read_num < tree_state.all_seeds->size() ; read_num++) {
                if (combined_cluster[read_num] != -1 ) {
                    chain_clusters.read_cluster_heads.emplace(read_num, combined_cluster[read_num]);
                    tree_state.read_cluster_dists[read_num][combined_cluster[read_num]] =
                                          make_pair(combined_left[read_num], combined_right[read_num]);
                    chain_clusters.fragment_best_left = min_not_minus_one(chain_clusters.fragment_best_left, combined_left[read_num]);
                    chain_clusters.fragment_best_right = min_not_minus_one(chain_clusters.fragment_best_right, combined_right[read_num]);
                    chain_clusters.read_best_left[read_num] = min_not_minus_one(chain_clusters.read_best_left[read_num], combined_left[read_num]);
                    chain_clusters.read_best_right[read_num] = min_not_minus_one(chain_clusters.read_best_right[read_num], combined_right[read_num]);

                }
            }

#ifdef DEBUG_CLUSTER
            cerr << "\t finished with snarl " << snarl_index.id_in_parent
                 << "with best distances " << chain_clusters.fragment_best_left
                 << " " << chain_clusters.fragment_best_right
                 << ", clusters:" <<endl;

            for (pair<size_t,size_t> c : chain_clusters.read_cluster_heads) {
                pair<int64_t, int64_t> dists = tree_state.read_cluster_dists[c.first][c.second];
                cerr << "\t\tleft: " << dists.first << " right : " << dists.second << endl;
                cerr << "\t\t\t";
                bool has_seeds = false;
                for (size_t x = 0 ; x < tree_state.all_seeds->at(c.first).size() ; x++) {
                    if (tree_state.read_union_find[c.first].find_group(x) == c.second) {
                        cerr << tree_state.all_seeds->at(c.first)[x] << " ";
                        has_seeds = true;
                    }
                }
                assert (has_seeds);
                cerr << endl;
            }
#endif

        }

        //Finished looping through all the snarls in the chain

        if (last_rank != chain_index.prefix_sum.size() - 2) {
            //If the last snarl we traversed was not the end of the chain,
            //Extend the right bound of each cluster to the end of the chain
            chain_clusters.fragment_best_right = -1;
            chain_clusters.read_best_right.assign(tree_state.all_seeds->size(), -1);
            int64_t last_dist = last_rank == 0 ? 0 : chain_index.prefix_sum[last_rank] - 1;
            int64_t dist_to_end = chain_index.chainLength() - last_dist - last_len;
            for (pair<size_t, size_t> cluster_head : chain_clusters.read_cluster_heads) {
                int64_t d = tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].second;
                tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].second = d == -1 ? -1: d + dist_to_end;
                chain_clusters.fragment_best_right = min_not_minus_one(chain_clusters.fragment_best_right,
                                       tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].second);
                chain_clusters.read_best_right[cluster_head.first] = min_not_minus_one(chain_clusters.read_best_right[cluster_head.first],
                                       tree_state.read_cluster_dists[cluster_head.first][cluster_head.second].second);
            }
        }


        if (chain_index.is_looping_chain) {
            //If the chain loops, then the clusters might be connected by
            //looping around the chain
            //
            int64_t first_length = chain_index.prefix_sum[0]-1;
            vector<pair<size_t, size_t>> to_erase; //old cluster group ids
            //New cluster- there will be at most one new cluster to add
            vector<size_t> combined_cluster (tree_state.all_seeds->size(), -1);
            size_t fragment_combined_cluster = -1;

            for (pair<size_t, size_t> cluster_head : chain_clusters.read_cluster_heads) {
                //For each chain cluster
                size_t read_num = cluster_head.first;
                pair<int64_t, int64_t>& chain_dists = tree_state.read_cluster_dists[read_num][cluster_head.second];

                if ((chain_dists.second != -1 && chain_clusters.read_best_left[read_num] != -1 &&
                     chain_dists.second + chain_clusters.read_best_left[read_num] - first_length - 1
                                                <= tree_state.read_distance_limit) ||
                   (chain_dists.first != -1 && chain_clusters.read_best_right[read_num] != -1 &&
                      chain_dists.first + chain_clusters.read_best_right[read_num] - first_length - 1
                                                <= tree_state.read_distance_limit)){
                    //If this chain cluster is in the combined cluster
                    if (combined_cluster[read_num] == -1) {
                        combined_cluster[read_num] = cluster_head.second;
                    } else {
                        tree_state.read_union_find[read_num].union_groups(combined_cluster[read_num], cluster_head.second);
                        if (tree_state.fragment_distance_limit != 0) {
                            if (fragment_combined_cluster != -1) {
                                tree_state.fragment_union_find.union_groups(fragment_combined_cluster, cluster_head.second+tree_state.read_index_offsets[read_num]);
                            }
                            fragment_combined_cluster = tree_state.fragment_union_find.find_group(cluster_head.second+tree_state.read_index_offsets[read_num]);
                        }
                        size_t new_group = tree_state.read_union_find[read_num].find_group(cluster_head.second);
                        if (new_group == cluster_head.second) {
                            to_erase.emplace_back(read_num, combined_cluster[read_num]);
                        } else {
                            to_erase.emplace_back(read_num, cluster_head.second);
                        }
                        combined_cluster[read_num] = new_group;
                    }

                    if (tree_state.fragment_distance_limit != 0) {
                        fragment_combined_cluster = tree_state.fragment_union_find.find_group(cluster_head.second + tree_state.all_seeds->size());
                    }
                } else if (tree_state.fragment_distance_limit != 0 &&
                   ((chain_dists.second != -1 && chain_clusters.fragment_best_left != -1 &&
                     chain_dists.second + chain_clusters.fragment_best_left - first_length - 1
                                                <= tree_state.fragment_distance_limit) ||
                   (chain_dists.first != -1 && chain_clusters.fragment_best_right != -1 &&
                      chain_dists.first + chain_clusters.fragment_best_right - first_length - 1
                                                <= tree_state.fragment_distance_limit))){
                    //If we can cluster by fragment
                    if (fragment_combined_cluster != -1) {
                        tree_state.fragment_union_find.union_groups(fragment_combined_cluster, cluster_head.second+tree_state.read_index_offsets[read_num]);
                    }
                    fragment_combined_cluster = tree_state.fragment_union_find.find_group(cluster_head.second+tree_state.read_index_offsets[read_num]);

                }
            }
            for (auto c : to_erase) {
                chain_clusters.read_cluster_heads.erase(c);
            }
            //Don't need to update best left and right distances because
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
            for (pair<size_t,size_t> c : chain_clusters.read_cluster_heads) {
                if (c.first == read_num) {
                    any_clusters = true;
                    pair<int64_t, int64_t> dists = tree_state.read_cluster_dists[c.first][c.second];
                    cerr << "\t" << c.first << ":"<<c.second << ": left: " << dists.first << " right : " << dists.second << ": ";
                    bool has_seeds = false;
                    for (size_t x = 0 ; x < tree_state.all_seeds->at(c.first).size() ; x++) {
                        if (tree_state.read_union_find[c.first].find_group(x) == c.second) {
                            cerr << tree_state.all_seeds->at(c.first)[x] << " ";
                            has_seeds = true;
                        }
                    }
                    assert(dists.first == -1 || dists.first >= chain_clusters.read_best_left[read_num]);
                    assert(dists.second == -1 || dists.second >= chain_clusters.read_best_right[read_num]);
                    assert(dists.first == -1 || dists.first >= chain_clusters.fragment_best_left);
                    assert(dists.second == -1 || dists.second >= chain_clusters.fragment_best_right);
                    if (dists.first == chain_clusters.fragment_best_left) {got_left = true;}
                    if (dists.second == chain_clusters.fragment_best_right) {got_right = true;}
                    if (dists.first == chain_clusters.read_best_left[read_num]) {got_read_left = true;}
                    if (dists.second == chain_clusters.read_best_right[read_num]) {got_read_right = true;}
                    cerr << endl;
                    assert(has_seeds);
                }
            }
            if (!chain_index.is_looping_chain) {
                assert(!any_clusters || got_read_left || chain_clusters.read_best_left[read_num] > tree_state.read_distance_limit || chain_clusters.read_best_left[read_num] == -1);
                assert(!any_clusters || got_read_right || chain_clusters.read_best_right[read_num] > tree_state.read_distance_limit || chain_clusters.read_best_right[read_num] == -1);
            }
        }

        if (!chain_index.is_looping_chain) {
            assert(got_left || chain_clusters.fragment_best_left > tree_state.fragment_distance_limit);
            assert(got_right ||chain_clusters.fragment_best_right > tree_state.fragment_distance_limit );
        }
        for (pair<size_t, size_t> group_id : chain_clusters.read_cluster_heads) {

            assert (group_id.second == tree_state.read_union_find[group_id.first].find_group(group_id.second));
        }
#endif

        return chain_clusters ;
    };



    SnarlSeedClusterer::NodeClusters SnarlSeedClusterer::cluster_one_snarl(
                    TreeState& tree_state, size_t snarl_index_i) const {
        /*Get the clusters on this snarl.
         * Nodes have not yet been clustered */
        MinimumDistanceIndex::SnarlIndex& snarl_index =
                                        dist_index.snarl_indexes[snarl_index_i];
#ifdef DEBUG_CLUSTER
        cerr << "Finding clusters on snarl number " << snarl_index_i
             << " headed by node " << snarl_index.id_in_parent << endl;
#endif

        //Keep track of all clusters on this snarl
        NodeClusters snarl_clusters(tree_state.all_seeds->size());

        auto combine_clusters = [&] (size_t& new_group, size_t& combined_group,
                                    size_t& fragment_combined_group,
                                    int64_t fragment_dist, int64_t read_dist, size_t read_num){
            //Helper function to compare and combine clusters in two nodes of the same snarl
            //If the distance between two clusters is small enough, then combine them
            //for the read clusters and, if applicable, for the fragment clusters
            //Updates the distances stored for the read clusters
            if (read_dist != -1 && read_dist <= tree_state.read_distance_limit) {
                //If the clusters are close enough to combine in the read
                if (tree_state.fragment_distance_limit != 0) {
                    if (fragment_combined_group != -1) {
                        //Also combine fragment clusters
                        tree_state.fragment_union_find.union_groups(new_group+tree_state.read_index_offsets[read_num], 
                                                                    fragment_combined_group);
                    }
                    fragment_combined_group = tree_state.fragment_union_find.find_group(new_group+tree_state.read_index_offsets[read_num]);
                }
                pair<int64_t, int64_t>& end_dists = tree_state.read_cluster_dists[read_num][new_group];

                if (combined_group == -1) {
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
                    if (new_g != new_group) {
                        snarl_clusters.read_cluster_heads.erase(make_pair(read_num,new_group));
                    }
                    if (new_g != combined_group) {
                        snarl_clusters.read_cluster_heads.erase(make_pair(read_num,combined_group));
                    }
                    snarl_clusters.read_cluster_heads.emplace(read_num,new_g);
                    end_dists = make_pair(
                                min_not_minus_one(end_dists.first, old_dists.first),
                                min_not_minus_one(end_dists.second, old_dists.second));
                    tree_state.read_cluster_dists[read_num][new_g] = end_dists;
                    new_group = new_g;
                    combined_group = new_g;
                }

            } else if (tree_state.fragment_distance_limit != 0
                  && fragment_dist <= tree_state.fragment_distance_limit) {

                //Same fragment
                if (fragment_combined_group != -1) {
                    //Also combine fragment clusters
                    tree_state.fragment_union_find.union_groups(new_group+tree_state.read_index_offsets[read_num], 
                                                                fragment_combined_group);
                }
                fragment_combined_group = tree_state.fragment_union_find.find_group(new_group+tree_state.read_index_offsets[read_num]);
            }
            return;
        };


        //Get the children of this snarl and their clusters
        vector<pair<NetgraphNode, NodeClusters>>& child_nodes =
                                       tree_state.snarl_to_nodes[snarl_index_i];
        int64_t start_length = snarl_index.nodeLength(0);
        int64_t end_length = snarl_index.nodeLength(snarl_index.num_nodes*2 -1);


        //Maps each cluster of child nodes to its left and right distances
        //of the node its on
        hash_map<pair<size_t,size_t>, pair<int64_t, int64_t>> old_dists;

        for (size_t i = 0; i < child_nodes.size() ; i++) {
            //Go through each child node of the netgraph and get clusters

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
            //Represents all the clusters on this child node
            NodeClusters& curr_child_clusters = child_nodes[i].second;

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
                for (size_t x = 0 ; x < tree_state.all_seeds->at(c.first).size() ; x++) {
                    if (tree_state.read_union_find[c.first].find_group(x) == c.second) {
                        cerr << tree_state.all_seeds->at(c.first)[x] << " ";
                    }
                }
                cerr << endl;
            }

            // Make sure the net graph node is actually in the net graph.
            assert(node_rank != numeric_limits<size_t>::max());
#endif

            vector<pair<size_t, size_t>> children_i(
                  make_move_iterator(curr_child_clusters.read_cluster_heads.begin()),
                  make_move_iterator(curr_child_clusters.read_cluster_heads.end()));
            for (size_t c_i = 0 ; c_i < children_i.size() ; c_i ++) {
                //for each cluster of child node i, find the distances to the
                //ends of the snarl

                pair<size_t, size_t> child_cluster_head = children_i[c_i];

                pair<int64_t, int64_t> dists_c = tree_state.read_cluster_dists[child_cluster_head.first][child_cluster_head.second];
                old_dists[child_cluster_head] = dists_c;

                pair<int64_t, int64_t> new_dists = snarl_index.distToEnds(node_rank,
                                        dists_c.first,dists_c.second);
#ifdef DEBUG_CLUSTER
cerr << "\tcluster: " << c_i << "dists to ends in snarl" << snarl_index.id_in_parent
     << " : " << new_dists.first << " " << new_dists.second << endl;
#endif

                snarl_clusters.fragment_best_left =min_not_minus_one(
                                   snarl_clusters.fragment_best_left,new_dists.first);
                snarl_clusters.fragment_best_right = min_not_minus_one(
                                   snarl_clusters.fragment_best_right, new_dists.second);
                snarl_clusters.read_best_left[child_cluster_head.first] =min_not_minus_one(
                                   snarl_clusters.read_best_left[child_cluster_head.first], new_dists.first);
                snarl_clusters.read_best_right[child_cluster_head.first] = min_not_minus_one(
                                   snarl_clusters.read_best_right[child_cluster_head.first], new_dists.second);


                snarl_clusters.read_cluster_heads.insert(child_cluster_head);
                tree_state.read_cluster_dists[child_cluster_head.first][child_cluster_head.second] = new_dists;
            }


            for (size_t j = 0 ; j <= i ; j++){
                //Go through other child net graph nodes up to and including i

                //Get the other node and its clusters
                NetgraphNode& other_node = child_nodes[j].first;
                NodeClusters& other_node_clusters = child_nodes[j].second;

                id_t other_node_id = other_node.id_in_parent(dist_index);
                //Rank of this node in the snarl
                size_t other_rank = other_node.rank_in_parent(dist_index,
                                                              other_node_id);
                size_t other_rev = other_rank % 2 == 0
                                    ? other_rank + 1 : other_rank - 1;

#ifdef DEBUG_CLUSTER
                cerr << "Other net graph node is " << typeToString(other_node.node_type)
                    << " headed by node " << other_node_id;


#endif


                //Find distance from each end of current node (i) to
                //each end of other node (j)
                int64_t dist_l_l = snarl_index.snarlDistance(
                                                     rev_rank, other_rank);
                int64_t dist_l_r = snarl_index.snarlDistance(
                                                      rev_rank, other_rev);
                int64_t dist_r_l = snarl_index.snarlDistance(
                                                    node_rank, other_rank);
                int64_t dist_r_r = snarl_index.snarlDistance(
                                                     node_rank, other_rev);
#ifdef DEBUG_CLUSTER
cerr << "\t distances between ranks " << node_rank << " and " << other_rank
     << ": " << dist_l_l << " " << dist_l_r << " " << dist_r_l << " "
     << dist_r_r << endl;
#endif

                //group ids of clusters combined between node i left and
                //node j left, etc
                vector<size_t> group_l_l (tree_state.all_seeds->size(), -1);
                vector<size_t> group_l_r (tree_state.all_seeds->size(), -1);
                vector<size_t> group_r_l (tree_state.all_seeds->size(), -1);
                vector<size_t> group_r_r (tree_state.all_seeds->size(), -1);
                size_t fragment_group_l_l = -1;
                size_t fragment_group_l_r = -1;
                size_t fragment_group_r_l = -1;
                size_t fragment_group_r_r = -1;

                if (max({dist_l_l, dist_l_r, dist_r_l, dist_r_r}) != -1
                   && ((tree_state.fragment_distance_limit == 0 &&
                         MinimumDistanceIndex::minPos({dist_l_l, dist_l_r,
                            dist_r_l, dist_r_r})-2 <= tree_state.read_distance_limit
                   && min_not_minus_one(curr_child_clusters.fragment_best_left, curr_child_clusters.fragment_best_right)-2
                                                <= tree_state.read_distance_limit) ||
                       (tree_state.fragment_distance_limit != 0 &&
                            MinimumDistanceIndex::minPos({dist_l_l, dist_l_r,
                            dist_r_l, dist_r_r})-2 <= tree_state.fragment_distance_limit
                   && min_not_minus_one(curr_child_clusters.fragment_best_left, curr_child_clusters.fragment_best_right)-2
                                                <= tree_state.fragment_distance_limit)
                                                )) {
                    //If the two nodes are reachable
                    for (size_t c_i = 0 ; c_i < children_i.size() ; c_i ++) {
                        //for each cluster of child node i

                        pair<size_t, size_t> child_cluster_head = children_i[c_i];
                        size_t read_num = child_cluster_head.first;
                        size_t c_group = tree_state.read_union_find[read_num].find_group(child_cluster_head.second);
                        pair<int64_t, int64_t> dists_c = old_dists[child_cluster_head];


                        if (dist_l_l != -1 && dists_c.first != -1
                                 && other_node_clusters.fragment_best_left != -1 ) {
                            //If cluster child_cluster_head can be combined with clusters in j
                            //from the left of both of them
                            int64_t read_dist = other_node_clusters.read_best_left[read_num] == -1 ? -1 :
                                  dist_l_l + dists_c.first + other_node_clusters.read_best_left[read_num]-1;
                            int64_t fragment_dist = dist_l_l + dists_c.first + other_node_clusters.fragment_best_left-1;
                            combine_clusters(c_group, group_l_l[read_num], fragment_group_l_l,
                                  fragment_dist, read_dist,  read_num);
                        }

                        if (dist_l_r != -1 && dists_c.first != -1
                            && other_node_clusters.fragment_best_right != -1 ) {
                            //If it can be combined from the left to the right of j
                            int64_t fragment_dist = dist_l_r + dists_c.first + other_node_clusters.fragment_best_right-1;
                            int64_t read_dist = other_node_clusters.read_best_right[read_num] == -1 ? -1 :
                                 dist_l_r + dists_c.first + other_node_clusters.read_best_right[read_num]-1;
                            combine_clusters(c_group, group_l_r[read_num], fragment_group_l_r,
                                 fragment_dist, read_dist, read_num);
                        }
                        if (dist_r_l != -1 && dists_c.second != -1
                            && other_node_clusters.fragment_best_left != -1 ) {
                            int64_t fragment_dist = dist_r_l + dists_c.second + other_node_clusters.fragment_best_left-1;
                            int64_t read_dist = other_node_clusters.read_best_left[read_num] == -1 ? -1 :
                                dist_r_l + dists_c.second + other_node_clusters.read_best_left[read_num]-1;
                            combine_clusters(c_group, group_r_l[read_num], fragment_group_r_l,
                                fragment_dist, read_dist,  read_num);
                        }
                        if (dist_r_r != -1 && dists_c.second != -1
                            && other_node_clusters.fragment_best_right != -1 ) {
                            int64_t fragment_dist = dist_r_r + dists_c.second + other_node_clusters.fragment_best_right-1;
                            int64_t read_dist = other_node_clusters.read_best_right[read_num] == -1 ? -1 :
                                dist_r_r + dists_c.second + other_node_clusters.read_best_right[read_num]-1;
                            combine_clusters(c_group, group_r_r[read_num], fragment_group_r_r,
                                fragment_dist, read_dist, read_num);
                        }

                    }
                    //Go through children of j
                    vector<pair<size_t, size_t>> children_j(
                             make_move_iterator(other_node_clusters.read_cluster_heads.begin()),
                             make_move_iterator(other_node_clusters.read_cluster_heads.end()));

                    for (size_t k_i = 0 ; k_i < children_j.size() ; k_i++){
                        //For each cluster of child j, find which overlaps with
                        //clusters of i
                        //child_cluster_head will already be part of a cluster in
                        //snarlcluster heads but since we need to know the node
                        //that the snarl is on we can't just loop through
                        //snarl_cluster heads
                        pair<size_t,size_t> child_cluster_head = children_j[k_i];
                        size_t read_num = child_cluster_head.first;
                        pair<int64_t, int64_t>& dists_k = old_dists[child_cluster_head];
                        size_t k_group = tree_state.read_union_find[read_num].find_group(child_cluster_head.second);


                        if (dist_l_l != -1 && curr_child_clusters.fragment_best_left != -1
                           && dists_k.first != -1 ){

                            int64_t fragment_dist = dist_l_l + curr_child_clusters.fragment_best_left + dists_k.first-1;
                            int64_t read_dist = curr_child_clusters.read_best_left[read_num] == -1 ? -1 :
                                dist_l_l + curr_child_clusters.read_best_left[read_num] + dists_k.first-1;
                            combine_clusters(k_group, group_l_l[read_num], fragment_group_l_l,
                                fragment_dist,read_dist, read_num);
                        }
                        if (dist_l_r != -1 && curr_child_clusters.fragment_best_left != -1
                             && dists_k.second != -1  ) {

                            int64_t read_dist = curr_child_clusters.read_best_left[read_num] == -1 ? -1 :
                               dist_l_r + curr_child_clusters.read_best_left[read_num] + dists_k.second-1;
                            int64_t fragment_dist = dist_l_r + curr_child_clusters.fragment_best_left + dists_k.second-1;
                            combine_clusters(k_group, group_l_r[read_num], fragment_group_l_r,
                               fragment_dist, read_dist, read_num);
                        }
                        if (dist_r_l != -1 && curr_child_clusters.fragment_best_right != -1
                            && dists_k.first != -1  ) {

                            int64_t fragment_dist = dist_r_l + curr_child_clusters.fragment_best_right + dists_k.first-1;
                            int64_t read_dist = curr_child_clusters.read_best_right[read_num] == -1 ? -1 :
                                dist_r_l + curr_child_clusters.read_best_right[read_num] + dists_k.first-1;
                            combine_clusters(k_group, group_r_l[read_num], fragment_group_r_l,
                                fragment_dist, read_dist, read_num);
                        }
                        if (dist_r_r != -1 && curr_child_clusters.fragment_best_right != -1
                           && dists_k.second != -1 ) {

                            int64_t fragment_dist = dist_r_r + curr_child_clusters.fragment_best_right + dists_k.second-1;
                            int64_t read_dist = curr_child_clusters.read_best_right[read_num] == -1 ? -1 :
                               dist_r_r + curr_child_clusters.read_best_right[read_num] + dists_k.second-1;
                            combine_clusters(k_group, group_r_r[read_num], fragment_group_r_r,
                               fragment_dist, read_dist, read_num);
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
                    for (size_t x = 0 ; x < tree_state.all_seeds->at(c.first).size() ; x++) {
                        if (tree_state.read_union_find[c.first].find_group(x) == c.second) {
                            cerr << tree_state.all_seeds->at(c.first)[x] << " ";
                            has_seeds = true;
                        }
                    }
                    assert(dists.first == -1 || dists.first >= snarl_clusters.read_best_left[read_num]);
                    assert(dists.second == -1 || dists.second >= snarl_clusters.read_best_right[read_num]);
                    assert(dists.first == -1 || dists.first >= snarl_clusters.fragment_best_left);
                    assert(dists.second == -1 || dists.second >= snarl_clusters.fragment_best_right);
                    if (dists.first == snarl_clusters.fragment_best_left) {got_left = true;}
                    if (dists.second == snarl_clusters.fragment_best_right) {got_right = true;}
                    if (dists.first == snarl_clusters.read_best_left[read_num]) {got_read_left = true;}
                    if (dists.second == snarl_clusters.read_best_right[read_num]) {got_read_right = true;}
                    cerr << endl;
                    assert(has_seeds);
                }
            }
            assert(!any_clusters ||got_read_left ||  snarl_clusters.read_best_left[read_num] == -1);
            assert(!any_clusters ||got_read_right ||  snarl_clusters.read_best_right[read_num] == -1);
        }
        assert(got_left);
        assert(got_right);

        for (pair<size_t, size_t> group_id : snarl_clusters.read_cluster_heads) {
            assert (group_id.second == tree_state.read_union_find[group_id.first].find_group(group_id.second));
        }
#endif
        return snarl_clusters;
    };
}
